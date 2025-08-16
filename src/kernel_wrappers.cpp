
//#include <Rcpp.h>
#include <vector>
#include <string>
#include "kernel_loader.h"
#include "OpenCL_helper.h"
#include "kernel_wrappers.h"
#include "kernel_runners.h"
#include <RcppArmadillo.h>
#include "famfuncs.h"
using namespace Rcpp;

using namespace OpenCLHelper;


// // forward declare the runner
// void f2_binomial_logit_prep_kernel_runner(
//     const std::string& kernel_source,
//     const char*        kernel_name,
//     int                l1,
//     int                l2,
//     int                m1,
//     const std::vector<double>& X_flat,
//     const std::vector<double>& B_flat,
//     const std::vector<double>& mu_flat,
//     const std::vector<double>& P_flat,
//     const std::vector<double>& alpha_flat,
//     std::vector<double>&       qf_flat,
//     std::vector<double>&       xb_flat,
//     int progbar
// );

// [[Rcpp::export]]
Rcpp::List f2_binomial_logit_prep_opencl(
    Rcpp::NumericMatrix b,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    int progbar 
) {
  int l1 = x.nrow(), l2 = x.ncol(), m1 = b.ncol();
  
  // prepare flat inputs
  auto X_flat     = flattenMatrix(x);
  auto B_flat     = flattenMatrix(b);
  auto mu_flat    = flattenMatrix(mu);
  auto P_flat     = flattenMatrix(P);
  auto alpha_flat = copyVector(alpha);
  

  
  // prepare outputs
  std::vector<double> qf_flat(m1);
  std::vector<double> xb_flat((size_t)l1 * m1);
  
  Rcpp::NumericVector qf(m1);
  Rcpp::NumericMatrix xb(l1, m1);
  
   
#ifdef USE_OPENCL

  
  // load kernels
  //std::string core_src = load_kernel_library("f2_binomial_logit_prep"); 
  std::string ksrc     = load_kernel_source("src/f2_binomial_logit_prep.cl");
//  std::string all_src  = core_src + "\n" + ksrc;
  std::string all_src  = ksrc;
  
  // call runner
  f2_binomial_logit_prep_kernel_runner(
    all_src,
    "f2_binomial_logit_prep",
    l1, l2, m1,
    X_flat, B_flat, mu_flat, P_flat, alpha_flat,
    qf_flat, xb_flat,
    progbar
  );
  

  
  // reconstruct R outputs
  for (int j = 0; j < m1; ++j) qf[j] = qf_flat[j];
  
  for (int j = 0; j < m1; ++j)
    for (int i = 0; i < l1; ++i)
      xb(i, j) = xb_flat[(size_t)j * l1 + i];
#else
  Rcpp::Rcout << "[INFO] OpenCL not available — returning zero vector/matrices.\n";
  
#endif 
  
  return Rcpp::List::create(
    Rcpp::Named("xb") = xb,
    Rcpp::Named("qf") = qf
  );
}



// [[Rcpp::export]]
Rcpp::List f2_binomial_logit_prep_grad_opencl(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar
) {
  int l1 = x.nrow(), l2 = x.ncol(), m1 = b.ncol();
  
  // flatten inputs
  auto X_flat     = flattenMatrix(x);
  auto B_flat     = flattenMatrix(b);
  auto mu_flat    = flattenMatrix(mu);
  auto P_flat     = flattenMatrix(P);
  auto alpha_flat = copyVector(alpha);
  auto y_flat     = copyVector(y);
  auto wt_flat    = copyVector(wt);
  
  // allocate outputs
  std::vector<double> qf_flat(m1);
  std::vector<double> xb_flat(static_cast<size_t>(l1) * m1);
  std::vector<double> grad_flat(static_cast<size_t>(m1) * l2);
  
  Rcpp::NumericMatrix xb(l1, m1);
  Rcpp::NumericVector qf(m1);
  
#ifdef USE_OPENCL
  
  // load & call kernel runner
  std::string ksrc    = load_kernel_source("src/f2_binomial_logit_prep.cl");
  std::string all_src = ksrc;
  
  f2_binomial_logit_prep_grad_kernel_runner(
    all_src,
    "f2_binomial_logit_prep_grad",
    l1, l2, m1,
    X_flat, B_flat, mu_flat, P_flat, alpha_flat,
    y_flat, wt_flat,
    qf_flat, xb_flat, grad_flat,
    progbar
  );
  
  // rebuild xb, qf exactly as before
  for (int j = 0; j < m1; ++j) {
    qf[j] = qf_flat[j];
    for (int i = 0; i < l1; ++i) {
      xb(i, j) = xb_flat[static_cast<size_t>(j) * l1 + i];
    }
  }
  
#else
  
  Rcpp::Rcout << "[INFO] OpenCL not available — returning zero vector/matrices.\n";
  
#endif
  
  // wrap gradient directly as arma::mat (m1 rows × l2 cols),
  // no copies, column-major data matches Armadillo & R
  arma::mat grad_arma(
      grad_flat.data(),  // pointer to your flat array
      m1,                // n_rows - dimensions
      l2,                // n_cols - grid points
      false,             // don't copy memory
      false              // strict = false
  );
  
  return Rcpp::List::create(
    Rcpp::Named("xb")   = xb,
    Rcpp::Named("qf")   = qf,
    Rcpp::Named("grad") = grad_arma
  );
}

//////////////////////////  GENERALIZATION /////////////////////





// [[Rcpp::export]]
Rcpp::List f2_prep_grad_opencl(
    std::string family,
    std::string link,
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar
) {
  int l1 = x.nrow(), l2 = x.ncol(), m1 = b.ncol();
  
  // flatten inputs
  auto X_flat     = flattenMatrix(x);
  auto B_flat     = flattenMatrix(b);
  auto mu_flat    = flattenMatrix(mu);
  auto P_flat     = flattenMatrix(P);
  auto alpha_flat = copyVector(alpha);
  auto y_flat     = copyVector(y);
  auto wt_flat    = copyVector(wt);
  
  // allocate outputs
  std::vector<double> qf_flat(m1);
  std::vector<double> xb_flat(static_cast<size_t>(l1) * m1);
  std::vector<double> grad_flat(static_cast<size_t>(m1) * l2);
  
  Rcpp::NumericMatrix xb(l1, m1);
  Rcpp::NumericVector qf(m1);
  
  // Dispatch kernel name and source
  std::string kernel_name;
  std::string kernel_file;
  std::string all_src;

  

  
#ifdef USE_OPENCL
  
   std::string OPENCL_source     = load_kernel_source("OPENCL.CL");
//   std::string rmath_source     = load_kernel_library("rmath");
//   std::string nmath_source     = load_kernel_library("nmath");
//   std::string dpq_source     = load_kernel_library("dpq");
   
   std::string rmath_source2     = load_kernel_source("rmath/Rmath.cl");
   std::string dpq_prelude_source     = load_kernel_source("dpq/dpq_prelude.cl");
   std::string dpq_source2     = load_kernel_source("dpq/dpq.cl");
//   - dpq_prelude
//   - dpq
   
   std::string nmath_source2     = load_kernel_source("nmath/nmath.cl");
   
   
   std::string chebyshev_source     = load_kernel_source("nmath/chebyshev.cl");
   std::string d1mach_source     = load_kernel_source("nmath/d1mach.cl");
   std::string dnorm_source     = load_kernel_source("nmath/dnorm.cl");
   std::string fmax2_source     = load_kernel_source("nmath/fmax2.cl");
   std::string gammalims_source     = load_kernel_source("nmath/gammalims.cl");
   std::string lgammacor_source     = load_kernel_source("nmath/lgammacor.cl");
   std::string log1p_source     = load_kernel_source("nmath/log1p.cl");
   std::string pnorm_source     = load_kernel_source("nmath/pnorm.cl");
   std::string stirlerr_large_source     = load_kernel_source("nmath/stirlerr_large.cl");
   std::string expm1_source     = load_kernel_source("nmath/expm1.cl");
   std::string gamma_source     = load_kernel_source("nmath/gamma.cl");
   std::string lgamma_source     = load_kernel_source("nmath/lgamma.cl");
   std::string lgamma1p_source     = load_kernel_source("nmath/lgamma1p.cl");
   std::string stirlerr_small_source     = load_kernel_source("nmath/stirlerr_small.cl");
   std::string dgamma_source     = load_kernel_source("nmath/dgamma.cl");
   std::string stirlerr_source     = load_kernel_source("nmath/stirlerr.cl");
   std::string bd0_source     = load_kernel_source("nmath/bd0.cl");
   std::string dbinom_source     = load_kernel_source("nmath/dbinom.cl");
   std::string dpois_source     = load_kernel_source("nmath/dpois.cl");
   std::string pgamma_source     = load_kernel_source("nmath/pgamma.cl");
   
   
   

  if (family == "binomial") {
    if (link == "logit") {
      kernel_name = "f2_binomial_logit_prep_grad";
      kernel_file = "src/f2_binomial_logit_prep.cl";
    } 
    else if (link == "probit") {
      kernel_name = "f2_binomial_probit_prep_grad";
      kernel_file = "src/f2_binomial_probit_prep.cl";
    }
    else if (link == "cloglog") {
      kernel_name = "f2_binomial_cloglog_prep_grad";
      kernel_file = "src/f2_binomial_cloglog_prep.cl";
    }
    else {
      Rcpp::stop("Unsupported link function for binomial family: " + link);
    }
  }
  
  else if (family =="poisson"){
    kernel_name = "f2_poisson_prep_grad";
    kernel_file  = "src/f2_poisson_prep.cl";
  }
  
  else if (family=="Gamma"){
    kernel_name = "f2_gamma_prep_grad";
    kernel_file  = "src/f2_gamma_prep.cl";
    
  }
  
  
  else {
    Rcpp::stop("Unsupported family: " + family);
  }  
    
    

  // load & call kernel runner
  std::string ksrc    = load_kernel_source(kernel_file);
  
  // For the probit model, we include some basic OpenCL enablement and inline functions 
  if (family == "binomial"&&link == "probit") {

     all_src = OPENCL_source +
//       "\n" +   rmath_source + 
       "\n" +   rmath_source2 + 
       //       "\n" + dpq_source +
       "\n" +dpq_prelude_source+
       "\n" +dpq_source2+
      "\n" +nmath_source2   
      + "\n" + chebyshev_source
      + "\n" + d1mach_source
      + "\n" + dnorm_source
      + "\n" + fmax2_source
      + "\n" + gammalims_source
      + "\n" + lgammacor_source
      + "\n" + log1p_source
      + "\n" + pnorm_source
//      + "\n" + stirlerr_large_source
//      + "\n" + expm1_source
//      + "\n" + gamma_source
//      + "\n" + lgamma_source
//      + "\n" + lgamma1p_source
//      + "\n" + stirlerr_small_source
//      + "\n" + stirlerr_source
//      + "\n" + bd0_source
//      + "\n" + dbinom_source
//      + "\n" + dpois_source
//      + "\n" + dgamma_source
      + "\n" +   ksrc;

 //    all_src= ksrc;
      }    
  else{
  all_src = ksrc;
  }
  
  
  f2_binomial_logit_prep_grad_kernel_runner(
    all_src,
    kernel_name.c_str(),  // ✅ convert to const char*
    l1, l2, m1,
    X_flat, B_flat, mu_flat, P_flat, alpha_flat,
    y_flat, wt_flat,
    qf_flat, xb_flat, grad_flat,
    progbar
  );
  
  // rebuild xb, qf exactly as before
  for (int j = 0; j < m1; ++j) {
    qf[j] = qf_flat[j];
    for (int i = 0; i < l1; ++i) {
      xb(i, j) = xb_flat[static_cast<size_t>(j) * l1 + i];
    }
  }
  
#else
  
  Rcpp::Rcout << "[INFO] OpenCL not available — returning zero vector/matrices.\n";
  
#endif
  
  // wrap gradient directly as arma::mat (m1 rows × l2 cols),
  // no copies, column-major data matches Armadillo & R
  arma::mat grad_arma(
      grad_flat.data(),  // pointer to your flat array
      m1,                // n_rows - dimensions
      l2,                // n_cols - grid points
      false,             // don't copy memory
      false              // strict = false
  );
  
  return Rcpp::List::create(
    Rcpp::Named("xb")   = xb,
    Rcpp::Named("qf")   = qf,
    Rcpp::Named("grad") = grad_arma
  );
}




// [[Rcpp::export]]
Rcpp::List f2_f3_opencl(
    std::string family,
    std::string link,
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar
) {
  int l1 = x.nrow(), l2 = x.ncol(), m1 = b.ncol();
  
  // flatten inputs
  auto X_flat     = flattenMatrix(x);
  auto B_flat     = flattenMatrix(b);
  auto mu_flat    = flattenMatrix(mu);
  auto P_flat     = flattenMatrix(P);
  auto alpha_flat = copyVector(alpha);
  auto y_flat     = copyVector(y);
  auto wt_flat    = copyVector(wt);
  
  

    
    Rcpp::Rcout << "Allocating Required Memory: "
                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                << "\n";

  
  
  
  // allocate outputs
  std::vector<double> qf_flat(m1);
  std::vector<double> grad_flat(static_cast<size_t>(m1) * l2);

  
  Rcpp::NumericVector qf(m1);
  

    Rcpp::Rcout << "Allocating xb Memory: "
                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                << "\n";

  
//    std::vector<double> xb_flat(static_cast<size_t>(l1) * m1);
  // Rcpp::NumericMatrix xb(l1, m1);
  

    Rcpp::Rcout << "Dispatch kernel: "
                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                << "\n";

  
    
  
  
  // Dispatch kernel name and source
  std::string kernel_name;
  std::string kernel_file;
  std::string all_src;
  

#ifdef USE_OPENCL
  
  std::string OPENCL_source     = load_kernel_source("OPENCL.CL");
     std::string rmath_source     = load_kernel_library("rmath","glmbayes", false);
     std::string nmath_source     = load_kernel_library("nmath","glmbayes", false);
     std::string dpq_source     = load_kernel_library("dpq","glmbayes", false);
  

  if (family == "binomial"||family == "quasibinomial") {
    if (link == "logit") {
      kernel_name = "f2_f3_binomial_logit";
      kernel_file = "src/f2_f3_binomial_logit.cl";
    } 
    else if (link == "probit") {
      kernel_name = "f2_f3_binomial_probit";
      kernel_file = "src/f2_f3_binomial_probit.cl";
      
          }
    else if (link == "cloglog") {

      kernel_name = "f2_f3_binomial_cloglog";
      kernel_file = "src/f2_f3_binomial_cloglog.cl";
      
                }
    else {
      Rcpp::stop("Unsupported link function for binomial family: " + link);
    }
  }
  
  else if (family =="poisson"||family =="quasipoisson"){
    kernel_name = "f2_f3_poisson";
    kernel_file  = "src/f2_f3_poisson.cl";
    
      }
  
  else if (family=="Gamma"){

    kernel_name = "f2_f3_gamma";
    kernel_file  = "src/f2_f3_gamma.cl";
    
        
  }
  
  
  else {
    Rcpp::stop("Unsupported family: " + family);
  }  
  
  
  
  // load & call kernel runner
  std::string ksrc    = load_kernel_source(kernel_file);
  
  
  /// Updated to use same "Program" logic for all models
  
  all_src = OPENCL_source +
           "\n" +   rmath_source + 
           "\n" + dpq_source +
    "\n" +nmath_source   
  + "\n" +   ksrc;
  
  
  
  f2_f3_kernel_runner(
    all_src,
    kernel_name.c_str(),  // ✅ convert to const char*
    l1, l2, m1,
    X_flat, B_flat, mu_flat, P_flat, alpha_flat,
    y_flat, wt_flat,
    qf_flat, 
    // xb_flat,
    grad_flat,
    progbar
  );
  
  // rebuild xb, qf exactly as before
  for (int j = 0; j < m1; ++j) {
    qf[j] = qf_flat[j];
    // for (int i = 0; i < l1; ++i) {
    //   xb(i, j) = xb_flat[static_cast<size_t>(j) * l1 + i];
    // }
  }
  
#else
  
  Rcpp::Rcout << "[INFO] OpenCL not available — returning zero vector/matrices.\n";
  
#endif
  
  // wrap gradient directly as arma::mat (m1 rows × l2 cols),
  // no copies, column-major data matches Armadillo & R
  arma::mat grad_arma(
      grad_flat.data(),  // pointer to your flat array
      m1,                // n_rows - dimensions
      l2,                // n_cols - grid points
      false,             // don't copy memory
      false              // strict = false
  );
  
  return Rcpp::List::create(
    // Rcpp::Named("xb")   = xb,
    Rcpp::Named("qf")   = qf,
    Rcpp::Named("grad") = grad_arma
  );
}





// [[Rcpp::export]]
Rcpp::NumericVector f2_accum(
    std::string family,
    std::string link,
    Rcpp::NumericMatrix xb,        // n × m matrix of π = P(y=1)
    Rcpp::NumericVector qf,        // length m: 0.5*(b-μ)'P(b-μ)
    Rcpp::NumericVector y,         // length n observed {0,1}
    Rcpp::NumericVector wt,        // length n weights
    int progbar                    // 0 = no bar, 1 = show bar
) {
  int n = xb.nrow();
  int m = xb.ncol();
  Rcpp::NumericVector res(m);
  
  for (int i = 0; i < m; ++i) {
    Rcpp::checkUserInterrupt();
    
    // extract column i of xb
    Rcpp::NumericVector xbi = xb(_, i);
    
    // dispatch to appropriate likelihood function
    Rcpp::NumericVector ll;
    
    if (family == "binomial") {
      ll = dbinom_glmb(y, wt, xbi, true);
    }
    
    else if (family=="poisson"){
    
      ll=dpois_glmb(y,xbi,true);
      
      for(int j=0;j<n;j++){ ll[j]=ll[j]*wt[j];    }
      
      
    }
    else if (family== "Gamma"){
      
      ll=dgamma_glmb(y,wt,xbi,true);
      
      

    }
    
    
    
    else {
      Rcpp::stop("Unsupported family: " + family);
    }
    
    // sum of log-likelihoods
    double sumll = std::accumulate(ll.begin(), ll.end(), 0.0);
    
    // total negative log-lik = quadratic form + (− sum log-lik)
    res[i] = qf[i] - sumll;
    // 
    // if (i == 1) {
    //   Rcpp::Rcout << "  sumxb = " << std::setprecision(17) << sumxb << "\n";
    //   Rcpp::Rcout << "  sum11_raw = " << std::setprecision(17) << sumll2 << "\n";
    //   Rcpp::Rcout << "  sum11 = " << std::setprecision(17) << sumll << "\n";
    //   Rcpp::Rcout << "  res[1] = " << std::setprecision(17) << res(i) << "\n";
    // }
    
  }
  
  return res;
}