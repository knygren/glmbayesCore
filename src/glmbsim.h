// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;



Rcpp::List  rnnorm_reg_std_cpp(int n,
                               NumericVector y,
                               NumericMatrix x,
                               NumericMatrix mu,
                               NumericMatrix P,
                               NumericVector alpha,
                               NumericVector wt,
                               Function f2,
                               Rcpp::List  Envelope,
                               Rcpp::CharacterVector   family,
                               Rcpp::CharacterVector   link, 
                               int progbar=1);

Rcpp::List rnnorm_reg_cpp(int n,
                          NumericVector y,
                          NumericMatrix x,
                          NumericVector mu,
                          NumericMatrix P,
                          NumericVector offset2,
                          NumericVector wt,
                          double dispersion,
                          Rcpp::List famfunc,
                          Function f1,
                          Function f2,
                          Function f3,
                          NumericVector start,
                          std::string family="binomial",
                          std::string link="logit",
                          int Gridtype=2);


Rcpp::List rnorm_reg_cpp(int n,NumericVector y,NumericMatrix x, 
                         NumericVector mu,NumericMatrix P,
                         NumericVector offset,NumericVector wt,
                         double dispersion,
                         Function f2,Function f3,
                         NumericVector start,
                         std::string family="gaussian",
                         std::string link="identity",
                         int Gridtype=2      
);



Rcpp::List  rindep_norm_gamma_reg_std_cpp(int n,NumericVector y,NumericMatrix x,
                                          NumericMatrix mu, /// This is typically standardized to be a zero vector
                                          NumericMatrix P, /// Part of prior precision shifted to the likelihood
                                          NumericVector alpha,NumericVector wt,
                                          Function f2,Rcpp::List  Envelope,
                                          Rcpp::List  gamma_list,
                                          Rcpp::List  UB_list,
                                          Rcpp::CharacterVector   family,Rcpp::CharacterVector   link,
                                          bool progbar=true,
                                          bool verbose=false);

Rcpp::List rindep_norm_gamma_reg_std_parallel_cpp(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,   // typically standardized to be a zero vector
    Rcpp::NumericMatrix P,    // part of prior precision shifted to the likelihood
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    Rcpp::Function f2,        // kept for signature parity
    Rcpp::List Envelope,
    Rcpp::List gamma_list,
    Rcpp::List UB_list,
    Rcpp::CharacterVector family,
    Rcpp::CharacterVector link,
    bool progbar = true,
    bool verbose = false
);


Rcpp::List glmb_Standardize_Model(
    NumericVector y, 
    NumericMatrix x,   // Original design matrix (to be adjusted)
    NumericMatrix P,   // Prior Precision Matrix (to be adjusted)
    NumericMatrix bstar, // Posterior Mode from optimization (to be adjusted)
    NumericMatrix A1  // Precision for Log-Posterior at posterior mode (to be adjusted)
);