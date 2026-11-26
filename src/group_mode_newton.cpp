// group_mode_newton.cpp
// Mode finding for one group by damped Newton on the shifted negative log
// posterior. The optimizer holds no family math: it calls an evaluator,
// either env's f2_f3_non_opencl with exact_hessian (NEWTON_ENV) or
// fam::f2_f3_f4_* (NEWTON_FAM). Every supported family/link pair is
// log-concave, so H = X' diag(w) X + P is positive definite and the final H
// is the A1 that glmb_Standardize_Model needs.

#include "RcppArmadillo.h"
#include "famfuncs.h"
#include "simfuncs_groups_split.h"

#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <exception>
#include <string>
#include <vector>

// Defined at global scope in EnvelopeEval.cpp. b is l1 x points; returns qf
// (points) and grad (points x l1).
Rcpp::List f2_f3_non_opencl(std::string family, std::string link,
                            Rcpp::NumericMatrix b, Rcpp::NumericVector y,
                            Rcpp::NumericMatrix x, Rcpp::NumericMatrix mu,
                            Rcpp::NumericMatrix P, Rcpp::NumericVector alpha,
                            Rcpp::NumericVector wt, int progbar);

namespace glmbayes {
namespace sim {
namespace group {

using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

// Same thresholds as MODE_COMPARE in group_mode_vmmin.cpp.
const double kTolD = 0.01;
const double kTolDf = 1e-6;
const double kTolHess = 1e-3;
const double kTolHessExact = 1e-4;

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
      .count();
}

// In-place Cholesky A = L L' (lower triangle of column-major A).
bool chol_lower(int n, double* A) {
  for (int j = 0; j < n; ++j) {
    double s = A[j + n * j];
    for (int k = 0; k < j; ++k) s -= A[j + n * k] * A[j + n * k];
    if (!(s > 0.0)) return false;
    const double d = std::sqrt(s);
    A[j + n * j] = d;
    for (int i = j + 1; i < n; ++i) {
      double t = A[i + n * j];
      for (int k = 0; k < j; ++k) t -= A[i + n * k] * A[j + n * k];
      A[i + n * j] = t / d;
    }
  }
  return true;
}

// Solve L L' x = v in place.
void chol_solve(int n, const double* L, double* v) {
  for (int i = 0; i < n; ++i) {
    double s = v[i];
    for (int k = 0; k < i; ++k) s -= L[i + n * k] * v[k];
    v[i] = s / L[i + n * i];
  }
  for (int i = n - 1; i >= 0; --i) {
    double s = v[i];
    for (int k = i + 1; k < n; ++k) s -= L[k + n * i] * v[k];
    v[i] = s / L[i + n * i];
  }
}

// NEWTON_ENV: the f2_f3_non_opencl call that vmmin_mode makes, cached at the
// last point, and exact_hessian.
struct EnvCtx {
  const GroupStd* g;
  std::string family;
  std::string link;
  int l1;
  NumericMatrix b;      // l1 x 1
  NumericMatrix mu1m;   // zero prior mean, l1 x 1
  std::vector<double> par;
  std::vector<double> cache_par;
  std::vector<double> cache_grad;
  double cache_val = 0.0;
  bool cache_ok = false;
  arma::mat He;
  double t_kernel = 0.0;
  double n_kernel = 0.0;
};

bool env_eval(void* ctx, const double* p, double* f, double* grad,
              double* H) {
  EnvCtx* d = static_cast<EnvCtx*>(ctx);
  const auto t0 = std::chrono::steady_clock::now();
  try {
    bool same = d->cache_ok;
    for (int i = 0; same && i < d->l1; ++i) same = d->cache_par[i] == p[i];
    if (!same) {
      for (int i = 0; i < d->l1; ++i) d->b[i] = p[i];
      const List out = f2_f3_non_opencl(d->family, d->link, d->b, d->g->y,
                                        d->g->x, d->mu1m, d->g->P,
                                        d->g->alpha, d->g->wt2, 0);
      const NumericVector qf = out["qf"];
      const NumericMatrix gr = out["grad"];
      d->cache_par.assign(p, p + d->l1);
      d->cache_val = qf[0];
      for (int i = 0; i < d->l1; ++i) d->cache_grad[i] = gr[i];
      d->cache_ok = true;
      d->n_kernel += 1.0;
    }
    *f = d->cache_val;
    if (grad) std::copy(d->cache_grad.begin(), d->cache_grad.end(), grad);
    if (H) {
      d->par.assign(p, p + d->l1);
      if (!exact_hessian(*d->g, d->par, d->family, d->link, d->He)) {
        d->t_kernel += seconds_since(t0);
        return false;
      }
      std::copy(d->He.begin(), d->He.end(), H);
    }
  } catch (std::exception&) {
    d->t_kernel += seconds_since(t0);
    return false;
  } catch (...) {
    d->t_kernel += seconds_since(t0);
    return false;
  }
  d->t_kernel += seconds_since(t0);
  return true;
}

// NEWTON_FAM: fam::f2_f3_f4_* on the group's data (l1 observations and l2
// coefficients in fam's naming).
struct FamCtx {
  glmbayes::fam::F2F3F4Fn fn;
  int nobs;
  int ncoef;
  const double* x;
  const double* y;
  const double* alpha;
  const double* wt;
  const double* P;
  const double* mu;     // zeros
  double* work;
  double t_kernel = 0.0;
  double n_kernel = 0.0;
};

bool fam_eval(void* ctx, const double* p, double* f, double* grad,
              double* H) {
  FamCtx* d = static_cast<FamCtx*>(ctx);
  const auto t0 = std::chrono::steady_clock::now();
  *f = d->fn(d->nobs, d->ncoef, d->x, d->y, d->alpha, d->wt, d->mu, d->P, p,
             grad, H, d->work);
  d->t_kernel += seconds_since(t0);
  d->n_kernel += 1.0;
  return true;
}

}  // namespace

int newton_minimize(int n, const NewtonEval& ev, double* b, double* H,
                    double* f_out, double* work, const NewtonOpts& opts,
                    NewtonStats& st) {
  double* g = work;
  double* d = work + n;
  double* bt = work + 2 * n;
  double* L = work + 3 * n;
  double f = 0.0;
  ++st.fevals;
  ++st.gevals;
  if (!ev.fn(ev.ctx, b, &f, g, H) || !R_FINITE(f)) return NEWTON_EVAL;
  *f_out = f;
  for (int it = 0;; ++it) {
    for (int k = 0; k < n * n; ++k) L[k] = H[k];
    if (!chol_lower(n, L)) return NEWTON_CHOL;
    for (int i = 0; i < n; ++i) d[i] = -g[i];
    chol_solve(n, L, d);
    double dec = 0.0;
    for (int i = 0; i < n; ++i) dec -= g[i] * d[i];
    if (!R_FINITE(dec)) return NEWTON_EVAL;
    if (0.5 * dec <= opts.tol) return NEWTON_OK;
    if (it >= opts.maxit) return NEWTON_MAXIT;
    ++st.iters;

    // Rounding in f near the mode must not reject a descent step.
    const double slack = 8.0 * DBL_EPSILON * (1.0 + std::fabs(f));
    double t = 1.0;
    bool accepted = false;
    for (int k = 0; k <= opts.max_halvings; ++k) {
      for (int i = 0; i < n; ++i) bt[i] = b[i] + t * d[i];
      double ft = 0.0;
      ++st.fevals;
      if (ev.fn(ev.ctx, bt, &ft, nullptr, nullptr) && R_FINITE(ft) &&
          ft <= f - opts.armijo * t * dec + slack) {
        accepted = true;
        break;
      }
      t *= 0.5;
      ++st.backtracks;
    }
    if (!accepted) return NEWTON_LINESEARCH;
    for (int i = 0; i < n; ++i) b[i] = bt[i];
    ++st.fevals;
    ++st.gevals;
    if (!ev.fn(ev.ctx, b, &f, g, H) || !R_FINITE(f)) return NEWTON_EVAL;
    *f_out = f;
  }
}

NewtonResult newton_mode(const GroupStd& g, const std::string& family,
                         const std::string& link, const double* par0,
                         int backend, std::vector<double>& work,
                         ModeTiming* timing) {
  auto t0 = std::chrono::steady_clock::now();
  const int l1 = g.x.ncol();
  const int l2 = g.x.nrow();
  const int nw = newton_work_size(l1);
  const int fw =
      backend == NEWTON_FAM ? glmbayes::fam::f2_f3_f4_work_size(l2, l1) : 0;
  const std::size_t need = static_cast<std::size_t>(nw) + fw + l1;
  if (work.size() < need) work.resize(need);
  double* mu0 = work.data() + nw + fw;
  std::fill(mu0, mu0 + l1, 0.0);

  NewtonResult res;
  res.b2a = NumericMatrix(l1, 1);
  res.A1 = NumericMatrix(l1, l1);
  double* b = res.b2a.begin();
  if (par0) {
    std::copy(par0, par0 + l1, b);
  } else {
    std::fill(b, b + l1, 0.0);
  }

  EnvCtx ec;
  FamCtx fc;
  NewtonEval ev;
  if (backend == NEWTON_FAM) {
    fc.fn = glmbayes::fam::f2_f3_f4_select(family, link);
    if (fc.fn == nullptr) return res;
    fc.nobs = l2;
    fc.ncoef = l1;
    fc.x = g.x.begin();
    fc.y = g.y.begin();
    fc.alpha = g.alpha.begin();
    fc.wt = g.wt2.begin();
    fc.P = g.P.begin();
    fc.mu = mu0;
    fc.work = work.data() + nw;
    ev.fn = fam_eval;
    ev.ctx = &fc;
  } else {
    ec.g = &g;
    ec.family = family;
    ec.link = link;
    ec.l1 = l1;
    ec.b = NumericMatrix(l1, 1);
    ec.mu1m = NumericMatrix(l1, 1);
    ec.cache_par.assign(l1, 0.0);
    ec.cache_grad.assign(l1, 0.0);
    ev.fn = env_eval;
    ev.ctx = &ec;
  }

  const NewtonOpts opts = NewtonOpts();
  res.status = newton_minimize(l1, ev, b, res.A1.begin(), &res.value,
                               work.data(), opts, res.stats);
  if (timing) {
    timing->newton += seconds_since(t0);
    timing->newton_kernel += backend == NEWTON_FAM ? fc.t_kernel : ec.t_kernel;
    timing->kernel_calls += backend == NEWTON_FAM ? fc.n_kernel : ec.n_kernel;
    timing->newton_iters += res.stats.iters;
    timing->newton_backtracks += res.stats.backtracks;
  }
  if (res.status != NEWTON_OK) return res;

  t0 = std::chrono::steady_clock::now();
  Rcpp::IntegerVector counts = Rcpp::IntegerVector::create(
      Rcpp::Named("function") = res.stats.fevals,
      Rcpp::Named("gradient") = res.stats.gevals);
  res.fit = List::create(
      Rcpp::Named("par") = NumericVector(b, b + l1),
      Rcpp::Named("value") = res.value, Rcpp::Named("counts") = counts,
      Rcpp::Named("convergence") = res.status,
      Rcpp::Named("message") = R_NilValue,
      Rcpp::Named("hessian") = res.A1);
  if (timing) timing->fit_list += seconds_since(t0);
  return res;
}

void newton_compare_update(NewtonCompare& nc, const GroupStd& g,
                           const std::string& family, const std::string& link,
                           const VmminResult& vm, const NewtonResult& nt,
                           double t_vm, double t_nt, bool kernel_checks) {
  const int l1 = g.x.ncol();
  const int l2 = g.x.nrow();
  nc.calls += 1.0;
  nc.iters += nt.stats.iters;
  nc.fevals += nt.stats.fevals;
  nc.gevals += nt.stats.gevals;
  nc.backtracks += nt.stats.backtracks;
  nc.t_vmmin += t_vm;
  nc.t_newton += t_nt;
  if (nt.status != NEWTON_OK) {
    nc.flagged += 1.0;
    return;
  }

  NumericMatrix bm(l1, 1);
  std::copy(nt.b2a.begin(), nt.b2a.end(), bm.begin());
  const List out = f2_f3_non_opencl(family, link, bm, g.y, g.x,
                                    NumericMatrix(l1, 1), g.P, g.alpha, g.wt2,
                                    0);
  const double f_env = Rcpp::as<NumericVector>(out["qf"])[0];
  const NumericMatrix gr_env = out["grad"];

  arma::vec delta(l1);
  for (int i = 0; i < l1; ++i) delta[i] = nt.b2a[i] - vm.b2a[i];
  const arma::mat Av(const_cast<double*>(vm.A1.begin()), l1, l1, false);
  const double q = arma::as_scalar(delta.t() * Av * delta);
  const double dist = std::sqrt(q > 0.0 ? q : 0.0);
  const double df = std::fabs(f_env - vm.value) / (1.0 + std::fabs(vm.value));
  if (f_env < vm.value) nc.newton_lower += 1.0;
  const arma::mat An(const_cast<double*>(nt.A1.begin()), l1, l1, false);
  const double hn = arma::norm(Av, "fro");
  const double hd = arma::norm(An - Av, "fro") / (hn > 0.0 ? hn : 1.0);

  bool bad = !(dist <= kTolD) || !(df <= kTolDf) || !(hd <= kTolHess);
  if (!(dist <= nc.max_d)) nc.max_d = dist;
  if (!(df <= nc.max_df)) nc.max_df = df;
  if (!(hd <= nc.max_hess)) nc.max_hess = hd;

  if (kernel_checks) {
    const glmbayes::fam::F2F3F4Fn fn =
        glmbayes::fam::f2_f3_f4_select(family, link);
    double kf = R_PosInf, kg = R_PosInf, kH = R_PosInf;
    if (fn != nullptr) {
      std::vector<double> w(glmbayes::fam::f2_f3_f4_work_size(l2, l1));
      std::vector<double> mu0(l1, 0.0), gk(l1);
      const double fk = fn(l2, l1, g.x.begin(), g.y.begin(), g.alpha.begin(),
                           g.wt2.begin(), mu0.data(), g.P.begin(),
                           nt.b2a.begin(), gk.data(), nullptr, w.data());
      kf = std::fabs(fk - f_env) / (1.0 + std::fabs(f_env));
      arma::vec dg(l1);
      for (int i = 0; i < l1; ++i) dg[i] = gk[i] - gr_env[i];
      arma::vec s;
      if (arma::solve(s, An, dg, arma::solve_opts::no_approx)) {
        const double qq = arma::dot(dg, s);
        kg = std::sqrt(qq > 0.0 ? qq : 0.0);
      }
      arma::mat He;
      std::vector<double> par(nt.b2a.begin(), nt.b2a.end());
      if (exact_hessian(g, par, family, link, He)) {
        const double en = arma::norm(He, "fro");
        kH = arma::norm(An - He, "fro") / (en > 0.0 ? en : 1.0);
      }
    }
    bad = bad || !(kf <= kTolDf) || !(kg <= kTolD) || !(kH <= kTolHessExact);
    if (!(kf <= nc.max_kf)) nc.max_kf = kf;
    if (!(kg <= nc.max_kgrad)) nc.max_kgrad = kg;
    if (!(kH <= nc.max_kH)) nc.max_kH = kH;
  }
  if (bad) nc.flagged += 1.0;
}

}  // namespace group
}  // namespace sim
}  // namespace glmbayes
