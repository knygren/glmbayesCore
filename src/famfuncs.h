// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;


namespace famfuncs {

NumericVector dbinom_glmb( NumericVector x, NumericVector N, NumericVector means, int lg);
NumericVector  f1_binomial_logit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
// arma::vec  f2_binomial_logit_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
arma::vec f2_binomial_logit_rmat(
    // NumericMatrix b, NumericVector y,
    //                              NumericMatrix x, NumericMatrix mu,
    //                              NumericMatrix P, NumericVector alpha,
    //                              NumericVector wt, int progbar = 0
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    int progbar);

arma::mat  f3_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);



List f2_binomial_logit_prep_v2(
    NumericMatrix b,       // l2 × m1
    NumericVector y,       // length l1
    NumericMatrix x,       // l1 × l2
    NumericMatrix mu,      // l2 × 1
    NumericMatrix P,       // l2 × l2
    NumericVector alpha,   // length l1
    NumericVector wt,      // length l1
    int progbar = 0
);


NumericVector f2_binomial_logit_accum(
    NumericMatrix xb,        // n × m matrix of π = P(y=1)
    NumericVector qf,        // length m: 0.5*(b-μ)'P(b-μ)
    NumericVector y,         // length n observed {0,1}
    NumericVector wt,        // length n weights
    int progbar          // 0 = no bar, 1 = show bar
);



NumericVector f2_accum(
    std::string family,
    std::string link,
    NumericMatrix xb,        // n × m matrix of π = P(y=1)
    NumericVector qf,        // length m: 0.5*(b-μ)'P(b-μ)
    NumericVector y,         // length n observed {0,1}
    NumericVector wt,        // length n weights
    int progbar          // 0 = no bar, 1 = show bar
);


NumericVector  f1_binomial_probit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
// arma::vec  f2_binomial_probit_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);
arma::vec f2_binomial_probit_rmat(
    // NumericMatrix b, NumericVector y,
    //                              NumericMatrix x, NumericMatrix mu,
    //                              NumericMatrix P, NumericVector alpha,
    //                              NumericVector wt, int progbar = 0
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    int progbar);

arma::mat  f3_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

NumericVector  f1_binomial_cloglog(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
// arma::vec  f2_binomial_cloglog_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
arma::vec f2_binomial_cloglog_rmat(
    // NumericMatrix b, NumericVector y,
    //                              NumericMatrix x, NumericMatrix mu,
    //                              NumericMatrix P, NumericVector alpha,
    //                              NumericVector wt, int progbar = 0
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    int progbar);


arma::mat  f3_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

void neg_dpois_glmb_rmat(const RVector<double>& x,     // observed counts
                         const std::vector<double>& means, // Poisson rates
                         std::vector<double>& res,         // output buffer (preallocated)
                         const int lg);                     // log=TRUE?
  
NumericVector dpois_glmb( NumericVector x, NumericVector means, int lg);
NumericVector  f1_poisson(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
// arma::vec   f2_poisson_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
// arma::vec f2_poisson_rmat_old(const NumericMatrix& b,
//                           const NumericVector& y,
//                           const NumericMatrix& x,
//                           const NumericMatrix& mu,
//                           const NumericMatrix& P,
//                           const NumericVector& alpha,
//                           const NumericVector& wt,
//                           int progbar);

arma::vec f2_poisson_rmat(const RMatrix<double>& b,       // candidate coefficients
                             const RVector<double>& y,       // observed counts
                             const RMatrix<double>& x,       // design matrix
                             const RMatrix<double>& mu,      // mode vector
                             const RMatrix<double>& P,       // precision matrix
                             const RVector<double>& alpha,   // predictor offset
                             const RVector<double>& wt,      // observation weights
                             int progbar );          // progress toggle

arma::mat  f3_poisson(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);


NumericVector dgamma_glmb( NumericVector x, NumericVector shape, NumericVector scale, int lg);
NumericVector  f1_gamma(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
// arma::vec  f2_gamma_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar);
arma::vec f2_gamma_rmat(const RMatrix<double>& b,       // candidate coefficients
                          const RVector<double>& y,       // observed counts
                          const RMatrix<double>& x,       // design matrix
                          const RMatrix<double>& mu,      // mode vector
                          const RMatrix<double>& P,       // precision matrix
                          const RVector<double>& alpha,   // predictor offset
                          const RVector<double>& wt,      // observation weights
                          int progbar );          // progress toggle

arma::mat  f3_gamma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt,int progbar);

NumericVector dnorm_glmb( NumericVector x, NumericVector means, NumericVector sds,int lg);
NumericVector  f1_gaussian(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt);
NumericVector  f2_gaussian(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
// arma::vec  f2_gaussian_arma(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt);
// Original version (unchanged)
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


arma::mat Inv_f3_with_disp_rmat_v2(
    const RcppParallel::RMatrix<double>& Pmat_r,
    const RcppParallel::RMatrix<double>& Pmu_r,
    const RcppParallel::RVector<double>& base_B0_r,
    const RcppParallel::RMatrix<double>& base_A_r,
    double dispersion,
    const RcppParallel::RMatrix<double>& cbars_r // p × m
);


}