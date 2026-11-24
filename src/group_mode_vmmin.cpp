// group_mode_vmmin.cpp
// Mode finding for one group by calling vmmin (the BFGS routine behind
// optim(method = "BFGS")) directly, with either the R f2/f3 closures or the
// C++ family kernels, followed by a replica of optimhess (or the exact
// Hessian of group_mode_hessian.cpp). With the R closures and the replica
// the result is bit-identical to optim(); the C++ kernels avoid the R
// callbacks.

#include "RcppArmadillo.h"
#include "simfuncs_groups_split.h"

#include <R_ext/Applic.h>

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

using Rcpp::Function;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

// optim() defaults for method = "BFGS".
const int kMaxit = 100;
const int kTrace = 0;
const int kReport = 10;
const double kNdeps = 1e-3;

// MODE_COMPARE flag thresholds for vmmin with the C++ f2/f3.
const double kTolD = 0.01;      // sqrt(delta' A1 delta)
const double kTolDf = 1e-6;     // |f_vmmin - f_optim| / (1 + |f_optim|)
const double kTolHess = 1e-3;   // ||H_vmmin - H_optim||_F / ||H_optim||_F

// HESS_COMPARE flag threshold, ||H_exact - H_numeric||_F / ||H_numeric||_F.
// optimhess's central differences (step 1e-3) are accurate to about 1e-6.
const double kTolHessExact = 1e-4;

struct ModeData {
  const GroupStd* g;
  Function* f2;
  Function* f3;
  std::string family;
  std::string link;
  bool r_callbacks;
  int l1;
  NumericVector mu1;    // zero prior mean, as passed to optim
  NumericMatrix mu1m;   // same, l1 x 1 for the C++ kernels
  NumericMatrix b;      // l1 x 1 scratch point for the C++ kernels
  std::vector<double> cache_par;
  std::vector<double> cache_grad;
  double cache_val = 0.0;
  bool cache_ok = false;
  bool failed = false;
  std::string err;
  double t_kernel = 0.0;
  double n_kernel = 0.0;
};

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
      .count();
}

void eval_cpp(ModeData* d, const double* p) {
  if (d->cache_ok) {
    bool same = true;
    for (int i = 0; i < d->l1; ++i) {
      if (d->cache_par[i] != p[i]) { same = false; break; }
    }
    if (same) return;
  }
  for (int i = 0; i < d->l1; ++i) d->b[i] = p[i];
  const auto t0 = std::chrono::steady_clock::now();
  const List out = f2_f3_non_opencl(d->family, d->link, d->b, d->g->y,
                                    d->g->x, d->mu1m, d->g->P, d->g->alpha,
                                    d->g->wt2, 0);
  d->t_kernel += seconds_since(t0);
  d->n_kernel += 1.0;
  const NumericVector qf = out["qf"];
  const NumericMatrix gr = out["grad"];
  d->cache_par.assign(p, p + d->l1);
  d->cache_val = qf[0];
  for (int i = 0; i < d->l1; ++i) d->cache_grad[i] = gr[i];
  d->cache_ok = true;
}

double eval_r_fn(ModeData* d, const double* p) {
  NumericVector par(p, p + d->l1);
  const NumericVector s = (*d->f2)(
      par, Rcpp::_["y"] = d->g->y, Rcpp::_["x"] = d->g->x,
      Rcpp::_["mu"] = d->mu1, Rcpp::_["P"] = d->g->P,
      Rcpp::_["alpha"] = d->g->alpha, Rcpp::_["wt"] = d->g->wt2);
  return s[0];
}

double vm_fn(int n, double* p, void* ex) {
  ModeData* d = static_cast<ModeData*>(ex);
  if (d->failed) return R_PosInf;
  try {
    if (d->r_callbacks) return eval_r_fn(d, p);
    eval_cpp(d, p);
    return d->cache_val;
  } catch (std::exception& e) {
    d->failed = true;
    d->err = e.what();
  } catch (...) {
    d->failed = true;
    d->err = "unknown error";
  }
  return R_PosInf;
}

void vm_gr(int n, double* p, double* df, void* ex) {
  ModeData* d = static_cast<ModeData*>(ex);
  for (int i = 0; i < n; ++i) df[i] = 0.0;
  if (d->failed) return;
  try {
    if (d->r_callbacks) {
      NumericVector par(p, p + n);
      const NumericVector s = (*d->f3)(
          par, Rcpp::_["y"] = d->g->y, Rcpp::_["x"] = d->g->x,
          Rcpp::_["mu"] = d->mu1, Rcpp::_["P"] = d->g->P,
          Rcpp::_["alpha"] = d->g->alpha, Rcpp::_["wt"] = d->g->wt2);
      for (int i = 0; i < n; ++i) df[i] = s[i];
      return;
    }
    eval_cpp(d, p);
    for (int i = 0; i < n; ++i) df[i] = d->cache_grad[i];
    return;
  } catch (std::exception& e) {
    d->failed = true;
    d->err = e.what();
  } catch (...) {
    d->failed = true;
    d->err = "unknown error";
  }
}

// Same steps and arithmetic as optimhess() in R's optim.c (parscale and
// fnscale are 1), including the final symmetrization.
NumericMatrix hessian_replica(ModeData* d, const std::vector<double>& par) {
  const int n = d->l1;
  NumericMatrix H(n, n);
  std::vector<double> dpar(par);
  std::vector<double> df1(n), df2(n);
  for (int i = 0; i < n; ++i) {
    const double eps = kNdeps / 1.0;
    dpar[i] = dpar[i] + eps;
    vm_gr(n, dpar.data(), df1.data(), d);
    dpar[i] = dpar[i] - 2 * eps;
    vm_gr(n, dpar.data(), df2.data(), d);
    for (int j = 0; j < n; ++j) {
      H[i * n + j] = 1.0 * (df1[j] - df2[j]) / (2 * eps * 1.0 * 1.0);
    }
    dpar[i] = dpar[i] + eps;
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      const double tmp = 0.5 * (H[i * n + j] + H[j * n + i]);
      H[i * n + j] = H[j * n + i] = tmp;
    }
  }
  return H;
}

bool bits_equal(const double* a, const double* b, R_xlen_t len) {
  for (R_xlen_t i = 0; i < len; ++i) {
    if (std::isnan(a[i]) && std::isnan(b[i])) continue;
    if (a[i] != b[i]) return false;
  }
  return true;
}

}  // namespace

bool vmmin_supported(const std::string& family, const std::string& link) {
  if (family == "binomial" || family == "quasibinomial") {
    return link == "logit" || link == "probit" || link == "cloglog";
  }
  return family == "poisson" || family == "quasipoisson" ||
         family == "Gamma" || family == "gaussian";
}

VmminResult vmmin_mode(const GroupStd& g, Function f2, Function f3,
                       const std::string& family, const std::string& link,
                       bool r_callbacks, ModeTiming* timing,
                       const double* par0, int hess_method,
                       HessCompare* hc) {
  auto t0 = std::chrono::steady_clock::now();
  const int l1 = g.x.ncol();

  ModeData d;
  d.g = &g;
  d.f2 = &f2;
  d.f3 = &f3;
  d.family = family;
  d.link = link;
  d.r_callbacks = r_callbacks;
  d.l1 = l1;
  d.mu1 = NumericVector(l1);
  d.mu1m = NumericMatrix(l1, 1);
  d.b = NumericMatrix(l1, 1);
  d.cache_par.assign(l1, 0.0);
  d.cache_grad.assign(l1, 0.0);

  std::vector<double> par(l1, 0.0);
  if (par0) par.assign(par0, par0 + l1);
  std::vector<int> mask(l1, 1);
  double fmin = 0.0;
  int fncount = 0, grcount = 0, fail = 0;

  // vmmin raises an R error on a non-finite initial value; stop cleanly
  // instead.
  const double f0 = vm_fn(l1, par.data(), &d);
  if (d.failed) Rcpp::stop("vmmin_mode: f2/f3 evaluation failed: %s", d.err);
  if (!R_FINITE(f0)) Rcpp::stop("vmmin_mode: initial value is not finite");

  vmmin(l1, par.data(), &fmin, vm_fn, vm_gr, kMaxit, kTrace, mask.data(),
        R_NegInf, std::sqrt(DBL_EPSILON), kReport, &d, &fncount, &grcount,
        &fail);
  if (d.failed) Rcpp::stop("vmmin_mode: f2/f3 evaluation failed: %s", d.err);
  const double k_vmmin = d.t_kernel;
  if (timing) {
    timing->vmmin += seconds_since(t0);
    timing->vmmin_kernel += k_vmmin;
  }
  t0 = std::chrono::steady_clock::now();

  NumericMatrix H;
  if (hess_method == HESS_NUMERIC) {
    H = hessian_replica(&d, par);
  } else {
    auto th = std::chrono::steady_clock::now();
    arma::mat He;
    const bool ok = exact_hessian(g, par, family, link, He);
    const double t_exact = seconds_since(th);
    if (ok && hess_method == HESS_EXACT) {
      H = NumericMatrix(l1, l1);
      std::copy(He.begin(), He.end(), H.begin());
    } else {
      th = std::chrono::steady_clock::now();
      H = hessian_replica(&d, par);
      const double t_numeric = seconds_since(th);
      if (hc) {
        if (!ok) hc->fallback += 1.0;
        if (hess_method == HESS_COMPARE) {
          hc->calls += 1.0;
          hc->t_exact += t_exact;
          hc->t_numeric += t_numeric;
          double hd = R_PosInf;
          if (ok) {
            const arma::mat Hn(H.begin(), l1, l1, false);
            const double hn = arma::norm(Hn, "fro");
            hd = arma::norm(He - Hn, "fro") / (hn > 0.0 ? hn : 1.0);
          }
          if (!(hd <= kTolHessExact)) hc->flagged += 1.0;
          if (ok && !(hd <= hc->max_hess)) hc->max_hess = hd;
        }
      }
    }
  }
  if (d.failed) Rcpp::stop("vmmin_mode: f3 evaluation failed: %s", d.err);
  if (timing) {
    timing->hessian += seconds_since(t0);
    timing->hessian_kernel += d.t_kernel - k_vmmin;
    timing->kernel_calls += d.n_kernel;
  }
  t0 = std::chrono::steady_clock::now();

  VmminResult res;
  res.b2a = NumericMatrix(l1, 1);
  for (int i = 0; i < l1; ++i) res.b2a[i] = par[i];
  res.A1 = H;
  res.value = fmin;
  res.convergence = fail;
  res.fncount = fncount;
  res.grcount = grcount;

  Rcpp::IntegerVector counts = Rcpp::IntegerVector::create(
      Rcpp::Named("function") = fncount, Rcpp::Named("gradient") = grcount);
  res.fit = List::create(
      Rcpp::Named("par") = NumericVector(par.begin(), par.end()),
      Rcpp::Named("value") = fmin, Rcpp::Named("counts") = counts,
      Rcpp::Named("convergence") = fail,
      Rcpp::Named("message") = R_NilValue, Rcpp::Named("hessian") = H);
  if (timing) timing->fit_list += seconds_since(t0);
  return res;
}

void warm_compare_update(WarmCompare& wc, const VmminResult& cold,
                         const VmminResult& warm, double t_cold,
                         double t_warm) {
  const int l1 = cold.b2a.nrow();
  wc.calls += 1.0;
  wc.fn_cold += cold.fncount;
  wc.gr_cold += cold.grcount;
  wc.fn_warm += warm.fncount;
  wc.gr_warm += warm.grcount;
  wc.t_cold += t_cold;
  wc.t_warm += t_warm;

  arma::vec delta(l1);
  for (int i = 0; i < l1; ++i) delta[i] = warm.b2a[i] - cold.b2a[i];
  const arma::mat A1(const_cast<double*>(cold.A1.begin()), l1, l1, false);
  const double q = arma::as_scalar(delta.t() * A1 * delta);
  const double dist = std::sqrt(q > 0.0 ? q : 0.0);
  const double df =
      std::fabs(warm.value - cold.value) / (1.0 + std::fabs(cold.value));
  const arma::mat Hw(const_cast<double*>(warm.A1.begin()), l1, l1, false);
  const double hn = arma::norm(A1, "fro");
  const double hd = arma::norm(Hw - A1, "fro") / (hn > 0.0 ? hn : 1.0);

  if (!(dist <= kTolD) || !(df <= kTolDf) || !(hd <= kTolHess) ||
      warm.convergence != 0) {
    wc.flagged += 1.0;
  }
  if (!(dist <= wc.max_d)) wc.max_d = dist;
  if (!(df <= wc.max_df)) wc.max_df = df;
  if (!(hd <= wc.max_hess)) wc.max_hess = hd;
}

void mode_compare_update(ModeCompare& cmp, const GroupStd& g,
                         const List& opt, const VmminResult& vr,
                         const VmminResult& vc, Function f2) {
  const int l1 = g.x.ncol();
  const NumericVector opar = Rcpp::as<NumericVector>(opt["par"]);
  const double oval = Rcpp::as<double>(opt["value"]);
  const NumericMatrix oH = Rcpp::as<NumericMatrix>(opt["hessian"]);
  const Rcpp::IntegerVector ocounts =
      Rcpp::as<Rcpp::IntegerVector>(opt["counts"]);

  cmp.calls += 1.0;
  cmp.fn_optim += ocounts[0];
  cmp.gr_optim += ocounts[1];
  cmp.fn_vmmin += vc.fncount;
  cmp.gr_vmmin += vc.grcount;

  const bool exact =
      opar.size() == l1 && oH.size() == vr.A1.size() &&
      bits_equal(opar.begin(), vr.b2a.begin(), l1) &&
      bits_equal(&oval, &vr.value, 1) &&
      bits_equal(oH.begin(), vr.A1.begin(), oH.size());
  if (!exact) cmp.exact_mismatch += 1.0;

  arma::vec delta(l1);
  for (int i = 0; i < l1; ++i) delta[i] = vc.b2a[i] - opar[i];
  const arma::mat A1(const_cast<double*>(oH.begin()), l1, l1, false);
  const double q = arma::as_scalar(delta.t() * A1 * delta);
  const double dist = std::sqrt(q > 0.0 ? q : 0.0);

  NumericVector cpar(vc.b2a.begin(), vc.b2a.end());
  const NumericVector fc = f2(
      cpar, Rcpp::_["y"] = g.y, Rcpp::_["x"] = g.x,
      Rcpp::_["mu"] = NumericVector(l1), Rcpp::_["P"] = g.P,
      Rcpp::_["alpha"] = g.alpha, Rcpp::_["wt"] = g.wt2);
  const double df = std::fabs(fc[0] - oval) / (1.0 + std::fabs(oval));

  const arma::mat Hc(const_cast<double*>(vc.A1.begin()), l1, l1, false);
  const double hn = arma::norm(A1, "fro");
  const double hd = arma::norm(Hc - A1, "fro") / (hn > 0.0 ? hn : 1.0);

  if (!(dist <= kTolD) || !(df <= kTolDf) || !(hd <= kTolHess) ||
      vc.convergence != 0) {
    cmp.flagged += 1.0;
  }
  if (!(dist <= cmp.max_d)) cmp.max_d = dist;
  if (!(df <= cmp.max_df)) cmp.max_df = df;
  if (!(hd <= cmp.max_hess)) cmp.max_hess = hd;
}

}  // namespace group
}  // namespace sim
}  // namespace glmbayes
