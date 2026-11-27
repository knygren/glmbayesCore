// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file famfuncs.h
 * @brief Family-specific likelihood and envelope components for glmbayes.
 *
 * @namespace glmbayes::fam
 * @brief Core log-likelihood, log-posterior, and corresponding log gradient components for each GLM family.
 *
 * @section ImplementedIn
 *   These declarations are implemented in:
 *     - famfuncs_binomial.cpp
 *     - famfuncs_Gamma.cpp
 *     - famfuncs_gaussian.cpp
 *     - famfuncs_poisson.cpp
 *
 * @section UsedBy
 *   These functions are consumed by:
 *     - Envelopefuncs.cpp
 *     - EnvelopeEval.cpp
 *     - EnvelopeBuild.cpp
 *     - rNormalGLM.cpp
 *     - rIndepNormalGammaReg.cpp
 *
 * @section Responsibilities
 *   Provides f1, f2, f3 components for all GLM families and link functions.
 *   Supplies both standard and parallel-safe (RMatrix/RVector) variants.
 *   
 *   f1 functions --> Negative log-likelihoods
 *   f2 functions --> Negative log-likelihood+ small prior component (used for models in standard form)
 *   f3 functions --> Gradients for negative log-likelihoods + small prior components (used for models in standard form)
 *   f4 functions --> Hessians of f2, X' diag(w) X + P (observed information for probit and cloglog, so w depends on y)
 *   inv_f3       --> Functions used to compute inverse gradients (with respect to dispersion) for gaussian families
 *
 *   f2_f3_f4 functions are the allocation-free single-point variant: raw
 *   pointers and a caller-supplied workspace, no Rcpp objects, so the same
 *   steps can be ported to OpenCL. They return the f2 value and fill the f3
 *   gradient and f4 Hessian when those outputs are non-null. f2_value
 *   functions are their value-only form (grad and H null).
 *
 *   Functions are used during both envelope construction and simulation
 *   
 */


#ifndef GLMBAYES_FAM_H
#define GLMBAYES_FAM_H

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppParallel.h>

#include <cstddef>
#include <string>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;




namespace glmbayes{

namespace fam {

NumericVector dbinom_glmb( NumericVector x, NumericVector N, NumericVector means, int lg);

enum class BinomLink { logit, probit, cloglog };

// Binomial log density (lg = 1) or density at linear predictor eta, with p and
// q = 1 - p computed directly from eta for the link (accurate when p is near 0 or 1).
NumericVector dbinom_glmb_eta(const NumericVector& x, const NumericVector& N,
                              const NumericVector& eta, BinomLink link, int lg);
NumericVector dpois_glmb( NumericVector x, NumericVector means, int lg);
NumericVector dgamma_glmb( NumericVector x, NumericVector shape, NumericVector scale, int lg);
NumericVector dnorm_glmb( NumericVector x, NumericVector means, NumericVector sds,int lg);

void neg_dpois_glmb_rmat(const RVector<double>& x,     // observed counts
                         const std::vector<double>& means, // Poisson rates
                         std::vector<double>& res,         // output buffer (preallocated)
                         const int lg);                     // log=TRUE?


//----------------- binomial_Logit -------------------------------------------------


NumericVector  f1_binomial_logit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
arma::vec f2_binomial_logit_rmat(
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    int progbar);
arma::mat  f3_binomial_logit(NumericMatrix b,
                             NumericVector y, 
                             NumericMatrix x,
                             NumericMatrix mu,
                             NumericMatrix P,
                             NumericVector alpha,
                             NumericVector wt,
                             int progbar);

// Combined f2/f3 evaluator (binomial–logit)
// Computes negative log‑posterior (f2) and gradient (f3) in one pass.
// Returns a list with:
//   $qf   : NumericVector length = m1   (same as f2_binomial_logit)
//   $grad : NumericMatrix m1 × l2       (same as f3_binomial_logit)
Rcpp::List f2_f3_binomial_logit(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);

//----------------- binomial_Probit -------------------------------------------------
  
NumericVector  f1_binomial_probit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
arma::vec f2_binomial_probit_rmat(
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    int progbar);
arma::mat  f3_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

Rcpp::List f2_f3_binomial_probit(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);

//----------------- binomial_cloglog -------------------------------------------------
  

NumericVector  f1_binomial_cloglog(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);

    arma::vec f2_binomial_cloglog_rmat(
      const RMatrix<double>& b,
      const RVector<double>& y,
      const RMatrix<double>& x,
      const RMatrix<double>& mu,
      const RMatrix<double>& P,
      const RVector<double>& alpha,
      const RVector<double>& wt,
      int progbar);
arma::mat  f3_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);


Rcpp::List f2_f3_binomial_cloglog(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);
  
//----------------- Poisson -------------------------------------------------

NumericVector  f1_poisson(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
arma::vec f2_poisson_rmat(const RMatrix<double>& b,       // candidate coefficients
                             const RVector<double>& y,       // observed counts
                             const RMatrix<double>& x,       // design matrix
                             const RMatrix<double>& mu,      // mode vector
                             const RMatrix<double>& P,       // precision matrix
                             const RVector<double>& alpha,   // predictor offset
                             const RVector<double>& wt,      // observation weights
                             int progbar );          // progress toggle

arma::mat  f3_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

Rcpp::List f2_f3_poisson(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);

//----------------- Gamma -------------------------------------------------
  
NumericVector  f1_gamma(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
arma::vec f2_gamma_rmat(const RMatrix<double>& b,       // candidate coefficients
                          const RVector<double>& y,       // observed counts
                          const RMatrix<double>& x,       // design matrix
                          const RMatrix<double>& mu,      // mode vector
                          const RMatrix<double>& P,       // precision matrix
                          const RVector<double>& alpha,   // predictor offset
                          const RVector<double>& wt,      // observation weights
                          int progbar );          // progress toggle

arma::mat  f3_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

Rcpp::List f2_f3_gamma(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);

//----------------- Gaussian -------------------------------------------------
  
NumericVector  f1_gaussian(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_gaussian(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
arma::vec f2_gaussian_rmat(const RcppParallel::RMatrix<double>& b,       // candidate coefficients
                           const RcppParallel::RVector<double>& y,       // observed counts
                           const RcppParallel::RMatrix<double>& x,       // design matrix
                           const RcppParallel::RMatrix<double>& mu,      // mode vector
                           const RcppParallel::RMatrix<double>& P,       // precision matrix
                           const RcppParallel::RVector<double>& alpha,   // predictor offset
                           const RcppParallel::RVector<double>& wt,      // observation weights
                           int progbar);                                 // progress toggle

// New parallel-safe variant (wt as RMatrix)
arma::vec f2_gaussian_rmat_mat(const RcppParallel::RMatrix<double>& b,   // candidate coefficients
                               const RcppParallel::RVector<double>& y,   // observed counts
                               const RcppParallel::RMatrix<double>& x,   // design matrix
                               const RcppParallel::RMatrix<double>& mu,  // mode vector
                               const RcppParallel::RMatrix<double>& P,   // precision matrix
                               const RcppParallel::RVector<double>& alpha, // predictor offset
                               const RcppParallel::RMatrix<double>& wt,  // observation weights (matrix view)
                               int progbar);                             // progress toggle

arma::mat  f3_gaussian(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);


Rcpp::List f2_f3_gaussian(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);

//----------------- Inv_f3 --Gaussian -------------------------------------------------
  
Rcpp::List Inv_f3_precompute_disp(NumericMatrix cbars,
                                  NumericVector y,
                                  NumericMatrix x,
                                  NumericMatrix mu,
                                  NumericMatrix P,
                                  NumericVector alpha,
                                  NumericVector wt);

// Dispersion-aware envelope solver
arma::mat Inv_f3_with_disp(Rcpp::List cache,
                           double dispersion,
                           Rcpp::NumericMatrix cbars_small);


arma::mat Inv_f3_with_disp_rmat(
    const RcppParallel::RMatrix<double>& Pmat_r,
    const RcppParallel::RMatrix<double>& Pmu_r,
    const RcppParallel::RVector<double>& base_B0_r,
    const RcppParallel::RMatrix<double>& base_A_r,
    double dispersion,
    const RcppParallel::RMatrix<double>& cbars_r // p × m
);

//----------------- f2_f3_f4 (allocation-free, single point) --------------------
//
// Dimensions follow the f2/f3 functions: l1 observations, l2 coefficients,
// x is l1 x l2 column-major, b, mu are l2 values, P is l2 x l2. grad (l2) and
// H (l2 x l2, column-major) may be nullptr. work holds
// f2_f3_f4_work_size(l1, l2) doubles and must not alias any other argument.
// The value, gradient and residuals use the same scalar formulas as the
// matching f2_f3_* function; the Hessian weights are the second derivatives
// in eta (clamped at 0, since every family/link here is log-concave).

typedef double (*F2F3F4Fn)(int l1, int l2, const double* x, const double* y,
                           const double* alpha, const double* wt,
                           const double* mu, const double* P, const double* b,
                           double* grad, double* H, double* work);

inline int f2_f3_f4_work_size(int l1, int l2) { return 3 * l1 + l2; }

double f2_f3_f4_binomial_logit(int l1, int l2, const double* x, const double* y,
                               const double* alpha, const double* wt,
                               const double* mu, const double* P, const double* b,
                               double* grad, double* H, double* work);
double f2_f3_f4_binomial_probit(int l1, int l2, const double* x, const double* y,
                                const double* alpha, const double* wt,
                                const double* mu, const double* P, const double* b,
                                double* grad, double* H, double* work);
double f2_f3_f4_binomial_cloglog(int l1, int l2, const double* x, const double* y,
                                 const double* alpha, const double* wt,
                                 const double* mu, const double* P, const double* b,
                                 double* grad, double* H, double* work);
double f2_f3_f4_poisson(int l1, int l2, const double* x, const double* y,
                        const double* alpha, const double* wt,
                        const double* mu, const double* P, const double* b,
                        double* grad, double* H, double* work);
double f2_f3_f4_gamma(int l1, int l2, const double* x, const double* y,
                      const double* alpha, const double* wt,
                      const double* mu, const double* P, const double* b,
                      double* grad, double* H, double* work);
double f2_f3_f4_gaussian(int l1, int l2, const double* x, const double* y,
                         const double* alpha, const double* wt,
                         const double* mu, const double* P, const double* b,
                         double* grad, double* H, double* work);

// The f2_f3_f4 function for a family/link pair (same pairs as
// f2_f3_non_opencl), or nullptr.
inline F2F3F4Fn f2_f3_f4_select(const std::string& family,
                                const std::string& link) {
  if (family == "binomial" || family == "quasibinomial") {
    if (link == "logit") return f2_f3_f4_binomial_logit;
    if (link == "probit") return f2_f3_f4_binomial_probit;
    if (link == "cloglog") return f2_f3_f4_binomial_cloglog;
    return nullptr;
  }
  if (family == "poisson" || family == "quasipoisson") return f2_f3_f4_poisson;
  if (family == "Gamma") return f2_f3_f4_gamma;
  if (family == "gaussian") return f2_f3_f4_gaussian;
  return nullptr;
}

//----------------- f2_value (allocation-free, value only) ----------------------
//
// The f2 value alone: the f2_f3_f4 function with grad and H null, so the
// same arguments, workspace (f2_value_work_size doubles) and scalar
// formulas. No Rcpp objects and no locks, so safe inside RcppParallel
// workers.

typedef double (*F2ValueFn)(int l1, int l2, const double* x, const double* y,
                            const double* alpha, const double* wt,
                            const double* mu, const double* P, const double* b,
                            double* work);

inline int f2_value_work_size(int l1, int l2) {
  return f2_f3_f4_work_size(l1, l2);
}

inline double f2_value_binomial_logit(int l1, int l2, const double* x,
                                      const double* y, const double* alpha,
                                      const double* wt, const double* mu,
                                      const double* P, const double* b,
                                      double* work) {
  return f2_f3_f4_binomial_logit(l1, l2, x, y, alpha, wt, mu, P, b, nullptr,
                                 nullptr, work);
}
inline double f2_value_binomial_probit(int l1, int l2, const double* x,
                                       const double* y, const double* alpha,
                                       const double* wt, const double* mu,
                                       const double* P, const double* b,
                                       double* work) {
  return f2_f3_f4_binomial_probit(l1, l2, x, y, alpha, wt, mu, P, b, nullptr,
                                  nullptr, work);
}
inline double f2_value_binomial_cloglog(int l1, int l2, const double* x,
                                        const double* y, const double* alpha,
                                        const double* wt, const double* mu,
                                        const double* P, const double* b,
                                        double* work) {
  return f2_f3_f4_binomial_cloglog(l1, l2, x, y, alpha, wt, mu, P, b, nullptr,
                                   nullptr, work);
}
inline double f2_value_poisson(int l1, int l2, const double* x,
                               const double* y, const double* alpha,
                               const double* wt, const double* mu,
                               const double* P, const double* b,
                               double* work) {
  return f2_f3_f4_poisson(l1, l2, x, y, alpha, wt, mu, P, b, nullptr, nullptr,
                          work);
}
inline double f2_value_gamma(int l1, int l2, const double* x, const double* y,
                             const double* alpha, const double* wt,
                             const double* mu, const double* P,
                             const double* b, double* work) {
  return f2_f3_f4_gamma(l1, l2, x, y, alpha, wt, mu, P, b, nullptr, nullptr,
                        work);
}
inline double f2_value_gaussian(int l1, int l2, const double* x,
                                const double* y, const double* alpha,
                                const double* wt, const double* mu,
                                const double* P, const double* b,
                                double* work) {
  return f2_f3_f4_gaussian(l1, l2, x, y, alpha, wt, mu, P, b, nullptr, nullptr,
                           work);
}

// The f2_value function for a family/link pair (the f2_f3_f4_select pairs),
// or nullptr.
inline F2ValueFn f2_value_select(const std::string& family,
                                 const std::string& link) {
  if (family == "binomial" || family == "quasibinomial") {
    if (link == "logit") return f2_value_binomial_logit;
    if (link == "probit") return f2_value_binomial_probit;
    if (link == "cloglog") return f2_value_binomial_cloglog;
    return nullptr;
  }
  if (family == "poisson" || family == "quasipoisson") return f2_value_poisson;
  if (family == "Gamma") return f2_value_gamma;
  if (family == "gaussian") return f2_value_gaussian;
  return nullptr;
}

// Steps shared by the f2_f3_f4 functions.

// bmu = b - mu; returns 0.5 bmu' P bmu and sets grad = P bmu when non-null.
inline double f2_f3_f4_prior(int l2, const double* P, const double* mu,
                             const double* b, double* bmu, double* grad) {
  for (int a = 0; a < l2; ++a) bmu[a] = b[a] - mu[a];
  double q = 0.0;
  for (int a = 0; a < l2; ++a) {
    double s = 0.0;
    for (int c = 0; c < l2; ++c) s += P[a + l2 * c] * bmu[c];
    q += bmu[a] * s;
    if (grad) grad[a] = s;
  }
  return 0.5 * q;
}

// eta = alpha + x b.
inline void f2_f3_f4_eta(int l1, int l2, const double* x, const double* alpha,
                         const double* b, double* eta) {
  for (int j = 0; j < l1; ++j) eta[j] = alpha[j];
  for (int k = 0; k < l2; ++k) {
    const double bk = b[k];
    const double* xk = x + static_cast<std::size_t>(l1) * k;
    for (int j = 0; j < l1; ++j) eta[j] += xk[j] * bk;
  }
}

// grad += x' r (grad already holds P bmu); H = x' diag(w) x + P.
inline void f2_f3_f4_finish(int l1, int l2, const double* x, const double* P,
                            const double* r, const double* w, double* grad,
                            double* H) {
  if (grad) {
    for (int k = 0; k < l2; ++k) {
      const double* xk = x + static_cast<std::size_t>(l1) * k;
      double s = 0.0;
      for (int j = 0; j < l1; ++j) s += xk[j] * r[j];
      grad[k] += s;
    }
  }
  if (H) {
    for (int c = 0; c < l2; ++c) {
      const double* xc = x + static_cast<std::size_t>(l1) * c;
      for (int a = c; a < l2; ++a) {
        const double* xa = x + static_cast<std::size_t>(l1) * a;
        double s = 0.0;
        for (int j = 0; j < l1; ++j) s += xa[j] * w[j] * xc[j];
        H[a + l2 * c] = s + 0.5 * (P[a + l2 * c] + P[c + l2 * a]);
        H[c + l2 * a] = H[a + l2 * c];
      }
    }
  }
}

// Multi-point f4: the l2 x l2 Hessian at each column of b (l2 x m1).
inline arma::cube f4_multi(F2F3F4Fn fn, const NumericMatrix& b,
                           const NumericVector& y, const NumericMatrix& x,
                           const NumericMatrix& mu, const NumericMatrix& P,
                           const NumericVector& alpha, const NumericVector& wt) {
  const int l1 = x.nrow();
  const int l2 = x.ncol();
  const int m1 = b.ncol();
  arma::cube out(l2, l2, m1);
  std::vector<double> work(f2_f3_f4_work_size(l1, l2));
  for (int i = 0; i < m1; ++i) {
    fn(l1, l2, x.begin(), y.begin(), alpha.begin(), wt.begin(), mu.begin(),
       P.begin(), b.begin() + static_cast<std::size_t>(l2) * i, nullptr,
       out.slice(i).memptr(), work.data());
  }
  return out;
}

arma::cube f4_binomial_logit(NumericMatrix b, NumericVector y, NumericMatrix x,
                             NumericMatrix mu, NumericMatrix P,
                             NumericVector alpha, NumericVector wt);
arma::cube f4_binomial_probit(NumericMatrix b, NumericVector y, NumericMatrix x,
                              NumericMatrix mu, NumericMatrix P,
                              NumericVector alpha, NumericVector wt);
arma::cube f4_binomial_cloglog(NumericMatrix b, NumericVector y, NumericMatrix x,
                               NumericMatrix mu, NumericMatrix P,
                               NumericVector alpha, NumericVector wt);
arma::cube f4_poisson(NumericMatrix b, NumericVector y, NumericMatrix x,
                      NumericMatrix mu, NumericMatrix P,
                      NumericVector alpha, NumericVector wt);
arma::cube f4_gamma(NumericMatrix b, NumericVector y, NumericMatrix x,
                    NumericMatrix mu, NumericMatrix P,
                    NumericVector alpha, NumericVector wt);
arma::cube f4_gaussian(NumericMatrix b, NumericVector y, NumericMatrix x,
                       NumericMatrix mu, NumericMatrix P,
                       NumericVector alpha, NumericVector wt);


} // famfuncs

} //glmbayes

#endif