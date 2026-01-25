// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;

List EnvelopeBuild_cpp(NumericVector bStar,
                     NumericMatrix A,
                     NumericVector y,
                     NumericMatrix x,
                     NumericMatrix mu,
                     NumericMatrix P,
                     NumericVector alpha,
                     NumericVector wt,
                     std::string family = "binomial",
                     std::string link   = "logit",
                     int Gridtype       = 2,
                     int n              = 1,
                     int n_envopt       = -1,   // NEW: effective sample size for EnvelopeOpt (defaults to n if -1)
                     bool sortgrid      = false,
                     bool use_opencl    = false, // Enables OpenCL acceleration during envelope construction
                     bool verbose       = false  // Enables diagnostic output
);



double rss_face_at_disp(double dispersion,
                        Rcpp::List cache,
                        Rcpp::NumericVector cbars_j,
                        Rcpp::NumericVector y,
                        Rcpp::NumericMatrix x,
                        Rcpp::NumericVector alpha,
                        Rcpp::NumericVector wt) ;

Rcpp::List EnvelopeSize(const arma::vec& a,
                        const Rcpp::NumericMatrix& G1,
                        int Gridtype   = 2,
                        int n          = 1000,
                        int n_envopt   = -1,
                        bool use_opencl = false,
                        bool verbose    = false);

Rcpp::List EnvelopeEval(const Rcpp::NumericMatrix& G4,   // grid (parameters × grid points)
                        const Rcpp::NumericVector& y,
                        const Rcpp::NumericMatrix& x,
                        const Rcpp::NumericMatrix& mu,
                        const Rcpp::NumericMatrix& P,
                        const Rcpp::NumericVector& alpha,
                        const Rcpp::NumericVector& wt,
                        const std::string& family,
                        const std::string& link,
                        bool use_opencl = false,
                        bool verbose = false);

List EnvelopeBuild_Ind_Normal_Gamma(NumericVector bStar,
                                    NumericMatrix A,
                                    NumericVector y, 
                                    NumericMatrix x,
                                    NumericMatrix mu,
                                    NumericMatrix P,
                                    NumericVector alpha,
                                    NumericVector wt,
                                    std::string family="binomial",
                                    std::string link="logit",
                                    int Gridtype=2, 
                                    int n=1,
                                    int n_envopt=-1,
                                    bool sortgrid=false,
                                    bool use_opencl    = false,
                                    bool verbose       = false);


NumericVector RSS(NumericVector y, NumericMatrix x,NumericMatrix b,NumericVector alpha,NumericVector wt);

double UB2(double dispersion,
           Rcpp::List cache,
           Rcpp::NumericVector cbars_j,
           Rcpp::NumericVector y,
           Rcpp::NumericMatrix x,
           Rcpp::NumericVector alpha,
           Rcpp::NumericVector wt,
           double rss_min_global);

List EnvelopeDispersionBuild_cpp(
    List Env,
    double Shape,
    double Rate,
    NumericMatrix P,
    NumericVector y,
    NumericMatrix x,
    NumericVector alpha,
    int n_obs,
    double RSS_post,
    double RSS_ML,
    NumericMatrix mu,         // ← new
    NumericVector wt,         // ← new
    double max_disp_perc = 0.99,
    Nullable<double> disp_lower = R_NilValue,
    Nullable<double> disp_upper = R_NilValue,
    bool verbose = false,
    bool use_parallel = true   // ← add flag here
  );


List Set_Grid(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint);
Rcpp::List Set_Grid_C(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint,Rcpp::NumericMatrix Down,Rcpp::NumericMatrix Up,Rcpp::NumericMatrix lglt,Rcpp::NumericMatrix lgrt,Rcpp::NumericMatrix lgct,Rcpp::NumericMatrix logU,Rcpp::NumericMatrix logP);
void Set_Grid_C2(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint,Rcpp::NumericMatrix Down,Rcpp::NumericMatrix Up,Rcpp::NumericMatrix lglt,Rcpp::NumericMatrix lgrt,Rcpp::NumericMatrix lgct,Rcpp::NumericMatrix logU,Rcpp::NumericMatrix logP);
void Set_Grid_C2_pointwise(Rcpp::NumericMatrix GIndex,  Rcpp::NumericMatrix cbars, Rcpp::NumericMatrix Lint,Rcpp::NumericMatrix Down,Rcpp::NumericMatrix Up,Rcpp::NumericMatrix lglt,Rcpp::NumericMatrix lgrt,Rcpp::NumericMatrix lgct,Rcpp::NumericMatrix logU,Rcpp::NumericMatrix logP);


List   setlogP(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3);
// Rcpp::List   setlogP_C(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst);
void setlogP_C2(NumericMatrix logP,NumericVector NegLL,NumericMatrix cbars,NumericMatrix G3,NumericMatrix LLconst);
