#include "RcppArmadillo.h"
#include "Envelopefuncs.h"
#include "openclPort.h"
#include "simfuncs.h"

using namespace openclPort;
using namespace glmbayes::env;
using namespace glmbayes::sim;



////////////// openclPort namespace ///////////////////


// [[Rcpp::export]]
std::string load_kernel_source_wrapper_export(
    const std::string& relative_path,
    const std::string& package = "glmbayes"
) {
  return load_kernel_source_wrapper(relative_path, package);
}


// [[Rcpp::export]]
std::string load_kernel_library_wrapper_export(
    const std::string& subdir,
    const std::string& package = "glmbayes",
    bool verbose = false
) {
  return load_kernel_library_wrapper(subdir, package, verbose);
}


// [[Rcpp::export]]
bool has_opencl_export() {
  return has_opencl();
}

// [[Rcpp::export]]
int get_opencl_core_count_export() {
  return get_opencl_core_count();
}


// [[Rcpp::export]]
Rcpp::CharacterVector gpu_names_export() {
  return gpu_names();
}

////////////// glmbayes::env namespace ////////////////

// [[Rcpp::export]]
Rcpp::List EnvelopeSize_export(
    const arma::vec& a,
    const Rcpp::NumericMatrix& G1,
    int Gridtype,
    int n,
    int n_envopt,
    bool use_opencl,
    bool verbose
) {
  return glmbayes::env::EnvelopeSize(
    a, G1, Gridtype, n, n_envopt, use_opencl, verbose
  );
}


// [[Rcpp::export]]
Rcpp::List EnvelopeBuild_cpp_export(
    Rcpp::NumericVector bStar,
    Rcpp::NumericMatrix A,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    std::string family,
    std::string link,
    int Gridtype,
    int n,
    int n_envopt,
    bool sortgrid,
    bool use_opencl,
    bool verbose
) {
  return glmbayes::env::EnvelopeBuild(
    bStar, A, y, x, mu, P, alpha, wt,
    family, link, Gridtype, n, n_envopt,
    sortgrid, use_opencl, verbose
  );
}


// [[Rcpp::export]]
Rcpp::List EnvelopeEval_export(
    const Rcpp::NumericMatrix& G4,      // grid (parameters × grid points)
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const std::string& family,
    const std::string& link,
    bool use_opencl = false,
    bool verbose = false
) {
  return EnvelopeEval(
    G4, y, x, mu, P, alpha, wt,
    family, link,
    use_opencl, verbose
  );
}


// [[Rcpp::export]]
Rcpp::List EnvelopeBuild_Ind_Normal_Gamma_export(
    const Rcpp::NumericVector& bStar,
    const Rcpp::NumericMatrix& A,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const std::string& family = "binomial",
    const std::string& link   = "logit",
    int Gridtype              = 2,
    int n                     = 1,
    int n_envopt              = -1,
    bool sortgrid             = false,
    bool use_opencl           = false,
    bool verbose              = false
) {
  return EnvelopeBuild_Ind_Normal_Gamma(
    bStar, A, y, x, mu, P, alpha, wt,
    family, link,
    Gridtype, n, n_envopt,
    sortgrid, use_opencl, verbose
  );
}


// [[Rcpp::export]]
Rcpp::List EnvelopeDispersionBuild_cpp_export(
    const Rcpp::List& Env,
    double Shape,
    double Rate,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    int n_obs,
    double RSS_post,
    double RSS_ML,
    const Rcpp::NumericMatrix& mu,        // new
    const Rcpp::NumericVector& wt,        // new
    double max_disp_perc = 0.99,
    Rcpp::Nullable<double> disp_lower = R_NilValue,
    Rcpp::Nullable<double> disp_upper = R_NilValue,
    bool verbose = false,
    bool use_parallel = true
) {
  return EnvelopeDispersionBuild(
    Env,
    Shape,
    Rate,
    P,
    y,
    x,
    alpha,
    n_obs,
    RSS_post,
    RSS_ML,
    mu,
    wt,
    max_disp_perc,
    disp_lower,
    disp_upper,
    verbose,
    use_parallel
  );
}


// [[Rcpp::export]]
Rcpp::List EnvelopeOrchestrator_cpp_export(
    const Rcpp::NumericVector& bstar2,
    const Rcpp::NumericMatrix& A,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x2,
    const Rcpp::NumericMatrix& mu2,
    const Rcpp::NumericMatrix& P2,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    
    int n,
    int Gridtype,
    Rcpp::Nullable<int> n_envopt,
    
    double shape,
    double rate,
    double RSS_Post2,
    double RSS_ML,
    
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    
    bool use_parallel,
    bool use_opencl,
    bool verbose
) {
  return EnvelopeOrchestrator(
    bstar2,
    A,
    y,
    x2,
    mu2,
    P2,
    alpha,
    wt,
    
    n,
    Gridtype,
    n_envopt,
    
    shape,
    rate,
    RSS_Post2,
    RSS_ML,
    
    max_disp_perc,
    disp_lower,
    disp_upper,
    
    use_parallel,
    use_opencl,
    verbose
  );
}


// [[Rcpp::export]]
Rcpp::List setlogP_export(
    const Rcpp::NumericMatrix& logP,
    const Rcpp::NumericVector& NegLL,
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericMatrix& G3
) {
  return EnvelopeSet_LogP(
    logP,
    NegLL,
    cbars,
    G3
  );
}

// [[Rcpp::export]]
Rcpp::List Set_Grid_export(
    const Rcpp::NumericMatrix& GIndex,
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericMatrix& Lint
) {
  return EnvelopeSet_Grid(
    GIndex,
    cbars,
    Lint
  );
}


// [[Rcpp::export]]
double rss_face_at_disp_export(
    double dispersion,
    const Rcpp::List& cache,
    const Rcpp::NumericVector& cbars_j,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt
) {
  return rss_face_at_disp(
    dispersion,
    cache,
    cbars_j,
    y,
    x,
    alpha,
    wt
  );
}

// [[Rcpp::export]]
double UB2_export(
    double dispersion,
    const Rcpp::List& cache,
    const Rcpp::NumericVector& cbars_j,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    double rss_min_global
) {
  return UB2(
    dispersion,
    cache,
    cbars_j,
    y,
    x,
    alpha,
    wt,
    rss_min_global
  );
}



////////////openclPort

// =====================================================================
//  rnnorm_reg_std_cpp_export
// =====================================================================

// [[Rcpp::export]]
Rcpp::List rnnorm_reg_std_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2,
    const Rcpp::List& Envelope,
    const Rcpp::CharacterVector& family,
    const Rcpp::CharacterVector& link,
    int progbar = 1,
    bool verbose = false
) {
  return rNormalGLM_std(
    n, y, x, mu, P, alpha, wt,
    f2, Envelope, family, link,
    progbar, verbose
  );
}

// =====================================================================
//  rnnorm_reg_cpp_export
// =====================================================================


// [[Rcpp::export]]
Rcpp::List rnnorm_reg_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double dispersion,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::NumericVector& start,
    const std::string& family = "binomial",
    const std::string& link   = "logit",
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
) {
  return rNormalGLM(
    n, y, x, mu, P, offset, wt,
    dispersion,
    f2, f3, start,
    family, link, Gridtype,
    n_envopt, use_parallel, use_opencl, verbose
  );
}

// =====================================================================
//  rnorm_reg_cpp_export
// =====================================================================

// [[Rcpp::export]]
Rcpp::List rnorm_reg_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double dispersion,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::NumericVector& start,
    const std::string& family = "gaussian",
    const std::string& link   = "identity",
    int Gridtype = 2
) {
  return rNormalReg(
    n, y, x, mu, P, offset, wt,
    dispersion, f2, f3, start,
    family, link, Gridtype
  );
}


// =====================================================================
//  rindep_norm_gamma_reg_std_cpp_export
// =====================================================================

// [[Rcpp::export]]
Rcpp::List rindep_norm_gamma_reg_std_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2,
    const Rcpp::List& Envelope,
    const Rcpp::List& gamma_list,
    const Rcpp::List& UB_list,
    const Rcpp::CharacterVector& family,
    const Rcpp::CharacterVector& link,
    bool progbar = true,
    bool verbose = false
) {
  return rIndepNormalGammaReg_std(
    n, y, x, mu, P, alpha, wt,
    f2, Envelope, gamma_list, UB_list,
    family, link, progbar, verbose
  );
}


// =====================================================================
//  rindep_norm_gamma_reg_std_parallel_cpp_export
// =====================================================================

// [[Rcpp::export]]
Rcpp::List rindep_norm_gamma_reg_std_parallel_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2,
    const Rcpp::List& Envelope,
    const Rcpp::List& gamma_list,
    const Rcpp::List& UB_list,
    const Rcpp::CharacterVector& family,
    const Rcpp::CharacterVector& link,
    bool progbar = true,
    bool verbose = false
) {
  return rIndepNormalGammaReg_std_parallel(
    n, y, x, mu, P, alpha, wt,
    f2, Envelope, gamma_list, UB_list,
    family, link, progbar, verbose
  );
}


// [[Rcpp::export]]
Rcpp::List rindep_norm_gamma_reg_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<Rcpp::NumericVector> disp_lower,
    Rcpp::Nullable<Rcpp::NumericVector> disp_upper,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    bool progbar
) {
  return rIndepNormalGammaReg(
    n,
    y,
    x,
    mu,
    P,
    offset,
    wt,
    shape,
    rate,
    max_disp_perc,
    disp_lower,
    disp_upper,
    Gridtype,
    n_envopt,
    use_parallel,
    use_opencl,
    verbose,
    progbar
  );
}


// =====================================================================
//  rGammaGaussian_cpp_export
// =====================================================================

// [[Rcpp::export]]
Rcpp::List rGammaGaussian_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& wt,
    const Rcpp::NumericVector& alpha,
    double shape,
    double rate,
    Rcpp::Nullable<double> disp_lower = R_NilValue,
    Rcpp::Nullable<double> disp_upper = R_NilValue,
    bool verbose = false
) {
  return glmbayes::sim::rGammaGaussian(
    n, y, x, beta, wt, alpha,
    shape, rate,
    disp_lower, disp_upper,
    verbose
  );
}


// =====================================================================
//  rGammaGamma_cpp_export
// =====================================================================

// [[Rcpp::export]]
Rcpp::List rGammaGamma_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& wt,
    const Rcpp::NumericVector& alpha,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower = R_NilValue,
    Rcpp::Nullable<double> disp_upper = R_NilValue,
    bool verbose = false
) {
  return glmbayes::sim::rGammaGamma(
    n, y, x, beta, wt, alpha,
    shape, rate, max_disp_perc,
    disp_lower, disp_upper,
    verbose
  );
}

// =====================================================================
//  glmb_Standardize_Model_export
// =====================================================================

// [[Rcpp::export]]
Rcpp::List glmb_Standardize_Model_export(
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericMatrix& bstar,
    const Rcpp::NumericMatrix& A1
) {
  return glmb_Standardize_Model(
    y, x, P, bstar, A1
  );
}




// extern "C" SEXP _glmbayes_EnvelopeSize(SEXP aSEXP,
//                                       SEXP G1SEXP,
//                                       SEXP GridtypeSEXP,
//                                       SEXP nSEXP,
//                                       SEXP n_envoptSEXP,
//                                       SEXP use_openclSEXP,
//                                       SEXP verboseSEXP)
// {
//   try {
//     // Convert inputs
//     arma::vec a = Rcpp::as<arma::vec>(aSEXP);
//     Rcpp::NumericMatrix G1 = Rcpp::as<Rcpp::NumericMatrix>(G1SEXP);
//     
//     int Gridtype   = Rcpp::as<int>(GridtypeSEXP);
//     int n          = Rcpp::as<int>(nSEXP);
//     int n_envopt   = Rcpp::as<int>(n_envoptSEXP);
//     bool use_opencl = Rcpp::as<bool>(use_openclSEXP);
//     bool verbose    = Rcpp::as<bool>(verboseSEXP);
//     
//     // Call the implementation (namespaced or not)
//     Rcpp::List out = glmbayes::envelopefuncs::EnvelopeSize(
//       a, G1, Gridtype, n, n_envopt, use_opencl, verbose
//     );
// 
//     return out;
//   }
//   catch (std::exception &ex) {
//     forward_exception_to_r(ex);
//   }
//   catch (...) {
//     Rcpp::stop("Unknown C++ exception in _glmbayes_EnvelopeSize");
//   }
//   
//   return R_NilValue; // never reached
// }
// 


// extern "C" SEXP _glmbayes_EnvelopeBuild_cpp(
//     SEXP bStarSEXP, SEXP ASEXP, SEXP ySEXP, SEXP xSEXP,
//     SEXP muSEXP, SEXP PSEXP, SEXP alphaSEXP, SEXP wtSEXP,
//     SEXP familySEXP, SEXP linkSEXP,
//     SEXP GridtypeSEXP, SEXP nSEXP, SEXP n_envoptSEXP,
//     SEXP sortgridSEXP, SEXP use_openclSEXP, SEXP verboseSEXP)
// {
//   try {
//     Rcpp::NumericVector bStar = Rcpp::as<Rcpp::NumericVector>(bStarSEXP);
//     Rcpp::NumericMatrix A     = Rcpp::as<Rcpp::NumericMatrix>(ASEXP);
//     Rcpp::NumericVector y     = Rcpp::as<Rcpp::NumericVector>(ySEXP);
//     Rcpp::NumericMatrix x     = Rcpp::as<Rcpp::NumericMatrix>(xSEXP);
//     Rcpp::NumericMatrix mu    = Rcpp::as<Rcpp::NumericMatrix>(muSEXP);
//     Rcpp::NumericMatrix P     = Rcpp::as<Rcpp::NumericMatrix>(PSEXP);
//     Rcpp::NumericVector alpha = Rcpp::as<Rcpp::NumericVector>(alphaSEXP);
//     Rcpp::NumericVector wt    = Rcpp::as<Rcpp::NumericVector>(wtSEXP);
//     
//     std::string family = Rcpp::as<std::string>(familySEXP);
//     std::string link   = Rcpp::as<std::string>(linkSEXP);
//     
//     int Gridtype = Rcpp::as<int>(GridtypeSEXP);
//     int n        = Rcpp::as<int>(nSEXP);
//     int n_envopt = Rcpp::as<int>(n_envoptSEXP);
//     
//     bool sortgrid   = Rcpp::as<bool>(sortgridSEXP);
//     bool use_opencl = Rcpp::as<bool>(use_openclSEXP);
//     bool verbose    = Rcpp::as<bool>(verboseSEXP);
//     
//     Rcpp::List out = glmbayes::envelopefuncs::EnvelopeBuild_cpp(
//       bStar, A, y, x, mu, P, alpha, wt,
//       family, link, Gridtype, n, n_envopt,
//       sortgrid, use_opencl, verbose
//     );
//     
//     return out;
//   }
//   catch (std::exception &ex) { forward_exception_to_r(ex); }
//   catch (...) { Rcpp::stop("Unknown C++ exception in _glmbayes_EnvelopeBuild_cpp"); }
//   
//   return R_NilValue;
// }
// 
// 



