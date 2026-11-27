// group_simulate_parallel.cpp
// One draw per chain for a group, run across chains with RcppParallel. Same
// accept-reject steps as rNormalGLM_std (rNormalGLM.cpp), on raw views of
// the standardized model and envelope, with the thread-local RNG.

#include "RcppArmadillo.h"
#include "RcppParallel.h"
#include "famfuncs.h"
#include "rng_utils.h"
#include "simfuncs_groups_split.h"

#include <Rmath.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <string>
#include <vector>

namespace glmbayes {
namespace sim {
namespace group {

using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

// rnorm_ct (rnorm_ct.cpp) without its qnorm lock: Rf_qnorm5 keeps no state.
double rnorm_ct_unlocked(double lgrt, double lglt, double mu, double sigma) {
  const double U = glmbayes::rng::runif_safe();
  if (lgrt >= lglt) {
    const double u1 = -std::expm1(lgrt);
    const double lgu1 = std::log(u1);
    const double lgU2 = std::log(U) + lglt + std::log(-std::expm1(lgu1 - lglt));
    const double lgU3 = lgU2 + std::log1p(std::exp(lgu1 - lgU2));
    return ::Rf_qnorm5(lgU3, mu, sigma, 1, 1);
  }
  const double e1mu2 = -std::expm1(lglt);
  const double lg1mu2 = std::log(e1mu2);
  const double lgU2 = std::log(U) + lgrt + std::log(-std::expm1(lg1mu2 - lgrt));
  const double lgU3 = lgU2 + std::log1p(std::exp(lg1mu2 - lgU2));
  return ::Rf_qnorm5(lgU3, mu, sigma, 0, 1);
}

// out = L2Inv (L3Inv b) + mu; t holds l1 doubles.
void back_transform(const GroupDrawIn& in, const double* b, double* t,
                    double* out) {
  const int l1 = in.l1;
  for (int a = 0; a < l1; ++a) {
    double s = 0.0;
    for (int c = 0; c < l1; ++c) s += in.L3Inv[a + l1 * c] * b[c];
    t[a] = s;
  }
  for (int a = 0; a < l1; ++a) {
    double s = 0.0;
    for (int c = 0; c < l1; ++c) s += in.L2Inv[a + l1 * c] * t[c];
    out[a] = s + in.mu[a];
  }
}

// The NumericMatrix f2 rNormalGLM_std uses for family/link.
double f2_reference(const std::string& family, const std::string& link,
                    const NumericMatrix& b, const GroupStd& g) {
  using namespace glmbayes::fam;
  NumericVector v;
  if (family == "binomial" || family == "quasibinomial") {
    if (link == "logit") {
      v = f2_binomial_logit(b, g.y, g.x2, g.mu2, g.P2, g.alpha, g.wt2, 0);
    } else if (link == "probit") {
      v = f2_binomial_probit(b, g.y, g.x2, g.mu2, g.P2, g.alpha, g.wt2, 0);
    } else {
      v = f2_binomial_cloglog(b, g.y, g.x2, g.mu2, g.P2, g.alpha, g.wt2, 0);
    }
  } else if (family == "poisson" || family == "quasipoisson") {
    v = f2_poisson(b, g.y, g.x2, g.mu2, g.P2, g.alpha, g.wt2, 0);
  } else if (family == "Gamma") {
    v = f2_gamma(b, g.y, g.x2, g.mu2, g.P2, g.alpha, g.wt2, 0);
  } else {
    v = f2_gaussian(b, g.y, g.x2, g.mu2, g.P2, g.alpha, g.wt2);
  }
  return v[0];
}

struct GroupSimWorker : public RcppParallel::Worker {
  const std::vector<GroupDrawIn>& in;
  glmbayes::fam::F2ValueFn f2v;
  double* coef;
  double* iters;
  std::size_t nc;

  GroupSimWorker(const std::vector<GroupDrawIn>& in_,
                 glmbayes::fam::F2ValueFn f2v_, double* coef_, double* iters_)
      : in(in_), f2v(f2v_), coef(coef_), iters(iters_), nc(in_.size()) {}

  void operator()(std::size_t begin, std::size_t end) {
    if (begin >= end) return;
    const int l1 = in[begin].l1;
    std::vector<double> work(group_draw_work_size(l1, in[begin].nobs));
    std::vector<double> out(l1);
    for (std::size_t c = begin; c < end; ++c) {
      const std::size_t need = group_draw_work_size(in[c].l1, in[c].nobs);
      if (work.size() < need) work.resize(need);
      iters[c] = group_draw_one(in[c], f2v, work.data(), out.data());
      for (int j = 0; j < l1; ++j) coef[c + nc * j] = out[j];
    }
  }
};

}  // namespace

int group_draw_work_size(int l1, int nobs) {
  return glmbayes::fam::f2_value_work_size(nobs, l1) + 2 * l1;
}

glmbayes::fam::F2ValueFn group_f2_value(const std::string& family,
                                        const std::string& link) {
  return glmbayes::fam::f2_value_select(family, link);
}

int group_sim_method() {
  SEXP opt = Rf_GetOption1(Rf_install("glmbayesCore.group_sim"));
  if (Rf_isNull(opt)) return SIM_PARALLEL;
  if (!Rf_isString(opt) || Rf_length(opt) != 1) {
    Rcpp::stop("options(glmbayesCore.group_sim) must be \"serial\", "
               "\"parallel\" or \"compare\".");
  }
  const std::string s = CHAR(STRING_ELT(opt, 0));
  if (s == "serial") return SIM_SERIAL;
  if (s == "parallel") return SIM_PARALLEL;
  if (s == "compare") return SIM_COMPARE;
  Rcpp::stop("options(glmbayesCore.group_sim) must be \"serial\", "
             "\"parallel\" or \"compare\".");
  return SIM_SERIAL;
}

GroupDrawIn group_draw_in(const GroupStd& g, const List& envelope) {
  GroupDrawIn in;
  in.l1 = g.x2.ncol();
  in.nobs = g.x2.nrow();
  in.y = g.y.begin();
  in.x2 = g.x2.begin();
  in.alpha = g.alpha.begin();
  in.wt2 = g.wt2.begin();
  in.mu2 = g.mu2.begin();
  in.P2 = g.P2.begin();
  in.L2Inv = g.L2Inv.begin();
  in.L3Inv = g.L3Inv.begin();
  in.mu = g.mu.begin();
  const NumericVector PLSD = envelope["PLSD"];
  const NumericVector LLconst = envelope["LLconst"];
  const NumericMatrix loglt = envelope["loglt"];
  const NumericMatrix logrt = envelope["logrt"];
  const NumericMatrix cbars = envelope["cbars"];
  in.n_cells = PLSD.size();
  if (LLconst.size() != in.n_cells || cbars.nrow() != in.n_cells ||
      cbars.ncol() != in.l1 || loglt.nrow() != in.n_cells ||
      logrt.nrow() != in.n_cells) {
    Rcpp::stop("group_draw_in: envelope dimensions do not match the model");
  }
  in.PLSD = PLSD.begin();
  in.LLconst = LLconst.begin();
  in.loglt = loglt.begin();
  in.logrt = logrt.begin();
  in.cbars = cbars.begin();
  return in;
}

int group_draw_one(const GroupDrawIn& in, glmbayes::fam::F2ValueFn f2v,
                   double* work, double* coef_out) {
  const int l1 = in.l1;
  const int nc = in.n_cells;
  double* b = work;
  double* t = work + l1;
  double* fwork = work + 2 * l1;
  int candidates = 0;
  for (;;) {
    ++candidates;
    double U = glmbayes::rng::runif_safe();
    int J = 0;
    while (J < nc - 1 && U > in.PLSD[J]) {
      U -= in.PLSD[J];
      ++J;
    }
    double bc = 0.0;
    for (int j = 0; j < l1; ++j) {
      const std::size_t k = static_cast<std::size_t>(J) +
                            static_cast<std::size_t>(nc) * j;
      b[j] = rnorm_ct_unlocked(in.logrt[k], in.loglt[k], -in.cbars[k], 1.0);
      bc += b[j] * in.cbars[k];
    }
    const double U2 = glmbayes::rng::runif_safe();
    const double f = f2v(in.nobs, l1, in.x2, in.y, in.alpha, in.wt2, in.mu2,
                         in.P2, b, fwork);
    if (in.LLconst[J] + bc - std::log(U2) - f >= 0.0) break;
  }
  back_transform(in, b, t, coef_out);
  return candidates;
}

void group_simulate_chains(const std::vector<GroupDrawIn>& in,
                           glmbayes::fam::F2ValueFn f2v, double* coef,
                           double* iters) {
  if (in.empty()) return;
  GroupSimWorker worker(in, f2v, coef, iters);
  RcppParallel::parallelFor(0, in.size(), worker, 8);
}

void sim_compare_update(SimCompare& sc, const GroupStd& g,
                        const List& envelope, const GroupDrawIn& in,
                        glmbayes::fam::F2ValueFn f2v,
                        const std::string& family, const std::string& link) {
  const int l1 = in.l1;
  std::vector<double> work(group_draw_work_size(l1, in.nobs));
  std::vector<double> pt(l1);
  NumericMatrix bm(l1, 1);
  double max_kf = 0.0;

  auto check_point = [&](const double* p) {
    for (int j = 0; j < l1; ++j) {
      pt[j] = p[j];
      bm(j, 0) = p[j];
    }
    const double f_ref = f2_reference(family, link, bm, g);
    const double f_new = f2v(in.nobs, l1, in.x2, in.y, in.alpha, in.wt2,
                             in.mu2, in.P2, pt.data(), work.data());
    double d = std::fabs(f_new - f_ref) / (1.0 + std::fabs(f_ref));
    if (!std::isfinite(d)) d = R_PosInf;
    max_kf = std::max(max_kf, d);
  };

  check_point(g.bstar2.begin());
  const NumericMatrix thetabars = envelope["thetabars"];
  const int n_cells = thetabars.nrow();
  if (n_cells > 0 && thetabars.ncol() == l1) {
    std::vector<double> row(l1);
    for (int r : {0, n_cells - 1}) {
      for (int j = 0; j < l1; ++j) row[j] = thetabars(r, j);
      check_point(row.data());
    }
  }

  std::vector<double> t(l1), out(l1);
  back_transform(in, g.bstar2.begin(), t.data(), out.data());
  double max_tr = 0.0;
  for (int j = 0; j < l1; ++j) {
    const double mode = g.b2a[j] + g.mu[j];
    double d = std::fabs(out[j] - mode) / (1.0 + std::fabs(mode));
    if (!std::isfinite(d)) d = R_PosInf;
    max_tr = std::max(max_tr, d);
  }

  sc.calls += 1.0;
  if (max_kf > 1e-10 || max_tr > 1e-10) sc.flagged += 1.0;
  sc.max_kf = std::max(sc.max_kf, max_kf);
  sc.max_transform = std::max(sc.max_transform, max_tr);
}

}  // namespace group
}  // namespace sim
}  // namespace glmbayes
