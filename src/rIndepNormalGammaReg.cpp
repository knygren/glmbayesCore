// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include "simfuncs.h"
#include "progress_utils.h"

#include <cmath>         // for std::log or std::exp if used
#include <math.h>
#include "rng_utils.h"  // for safe_runif()

#include "nmath_local.h"
#include "dpq_local.h"

// Required headers
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
#include <tbb/mutex.h>
static tbb::mutex f2_mutex;
#endif
#include <string>
#include <atomic>
#include <memory>



using namespace Rcpp;
using namespace glmbayes::fam;
using namespace glmbayes::env;
using namespace glmbayes::sim;
using namespace glmbayes::rng;
using namespace glmbayes::progress;




// ------------------------------------------------------------
// g1_face_at_disp
// ------------------------------------------------------------
// Computes the quadratic–linear face energy:
//
//   g1_j(d) = -0.5 * θ_j(d)^T P θ_j(d) + c_j^T θ_j(d)
//
// where θ_j(d) is obtained exactly the same way as in the main
// code, i.e. via:
//
//   arma::mat theta2 = Inv_f3_with_disp(cache, dispersion, transpose(cbars_small));
//
// with cbars_small being a 1 × p NumericMatrix slice of cbars.
// ------------------------------------------------------------

double g1_face_at_disp(
    double dispersion,
    int j,
    const Rcpp::List& cache,
    const arma::mat& P2,
    const Rcpp::NumericMatrix& cbars
) {
  // ---- 1) Build cbars_small exactly as in the main code: 1 × p ----
  int p = cbars.ncol();
  Rcpp::NumericMatrix cbars_small(1, p);
  for (int k = 0; k < p; ++k) {
    cbars_small(0, k) = cbars(j, k);
  }
  
  // ---- 2) Call Inv_f3_with_disp with transpose(cbars_small), as in main ----
  Rcpp::NumericMatrix cbars_small_t = Rcpp::transpose(cbars_small);
  arma::mat theta_mat = Inv_f3_with_disp(cache, dispersion, cbars_small_t);
  // theta_mat is 1 × p (Out.t() in Inv_f3_with_disp)
  
  // Convert 1 × p → p × 1, matching main logic
  arma::rowvec theta_row = theta_mat.row(0);
  arma::vec theta = theta_row.t();
  
  // ---- 3) Build c_j as arma::vec (p × 1) ----
  arma::vec c_j(p);
  for (int k = 0; k < p; ++k) {
    c_j(k) = cbars(j, k);
  }
  
  // ---- 4) Quadratic–linear face energy g1_j(d) ----
  double g1 = arma::as_scalar(
    -0.5 * theta.t() * P2 * theta
  +       c_j.t()  * theta
  );
  
  return g1;
}


double g2_face_at_disp(
    double dispersion,
    int j,
    const Rcpp::List& cache,
    const arma::mat& P2,
    const Rcpp::NumericMatrix& cbars,
    const double& d1_star,
    const Rcpp::NumericVector& New_LL_Slope
) {

  // Linear approximation
  double g2 = g1_face_at_disp(d1_star,j,cache,P2,cbars)+New_LL_Slope[j]*(dispersion-d1_star);
  
  return g2;
}


//-----------------------------------------------------------------------------
// rIndepNormalGammaReg_worker: parallel Normal–Gamma simulation with envelope
//-----------------------------------------------------------------------------
struct rIndepNormalGammaReg_worker : public RcppParallel::Worker {
  // --- Inputs ---
  int n;
  
  // Likelihood inputs (thread-safe views)
  RcppParallel::RVector<double>       y_r;
  RcppParallel::RMatrix<double>       x_r;
  RcppParallel::RMatrix<double>       mu_r;
  RcppParallel::RMatrix<double>       P_r;
  RcppParallel::RVector<double>       alpha_r;
  RcppParallel::RVector<double>       wt_r;
  
  // Envelope components
  RcppParallel::RMatrix<double>       cbars_r;
  RcppParallel::RVector<double>       PLSD_r;
  RcppParallel::RMatrix<double>       loglt_r;
  RcppParallel::RMatrix<double>       logrt_r;
  
  // UB vectors
  RcppParallel::RVector<double>       lg_prob_factor_r;
  RcppParallel::RVector<double>       UB2min_r;
  
  // Scalars
  double shape3, rate2, disp_upper, disp_lower, RSS_Min;
  double max_New_LL_UB, max_LL_log_disp, lm_log1, lm_log2, lmc1, lmc2;
  
  // Cache (precomputed upstream)
  RcppParallel::RMatrix<double>       Pmat_r;
  RcppParallel::RMatrix<double>       Pmu_r;
  RcppParallel::RVector<double>       base_B0_r;
  RcppParallel::RMatrix<double>       base_A_r;
  
  // --- Outputs ---
  RcppParallel::RMatrix<double>       beta_out_r;   // n × l1
  RcppParallel::RVector<double>       disp_out_r;   // length n
  RcppParallel::RVector<double>       iters_out_r;  // length n
  RcppParallel::RVector<double>       weight_out_r; // length n
  
  // --- Constructor ---
  rIndepNormalGammaReg_worker(
    int n_,
    const RcppParallel::RVector<double>& y_r_,
    const RcppParallel::RMatrix<double>& x_r_,
    const RcppParallel::RMatrix<double>& mu_r_,
    const RcppParallel::RMatrix<double>& P_r_,
    const RcppParallel::RVector<double>& alpha_r_,
    const RcppParallel::RVector<double>& wt_r_,
    const RcppParallel::RMatrix<double>& cbars_r_,
    const RcppParallel::RVector<double>& PLSD_r_,
    const RcppParallel::RMatrix<double>& loglt_r_,
    const RcppParallel::RMatrix<double>& logrt_r_,
    const RcppParallel::RVector<double>& lg_prob_factor_r_,
    const RcppParallel::RVector<double>& UB2min_r_,
    double shape3_, double rate2_,
    double disp_upper_, double disp_lower_,
    double RSS_Min_,
    double max_New_LL_UB_, double max_LL_log_disp_,
    double lm_log1_, double lm_log2_,
    double lmc1_, double lmc2_,
    const RcppParallel::RMatrix<double>& Pmat_r_,
    const RcppParallel::RMatrix<double>& Pmu_r_,
    const RcppParallel::RVector<double>& base_B0_r_,
    const RcppParallel::RMatrix<double>& base_A_r_,
    RcppParallel::RMatrix<double>& beta_out_r_,
    RcppParallel::RVector<double>& disp_out_r_,
    RcppParallel::RVector<double>& iters_out_r_,
    RcppParallel::RVector<double>& weight_out_r_)
    : n(n_),
      y_r(y_r_), x_r(x_r_), mu_r(mu_r_), P_r(P_r_), alpha_r(alpha_r_), wt_r(wt_r_),
      cbars_r(cbars_r_), PLSD_r(PLSD_r_), loglt_r(loglt_r_), logrt_r(logrt_r_),
      lg_prob_factor_r(lg_prob_factor_r_), UB2min_r(UB2min_r_),
      shape3(shape3_), rate2(rate2_), disp_upper(disp_upper_), disp_lower(disp_lower_),
      RSS_Min(RSS_Min_), max_New_LL_UB(max_New_LL_UB_), max_LL_log_disp(max_LL_log_disp_),
      lm_log1(lm_log1_), lm_log2(lm_log2_), lmc1(lmc1_), lmc2(lmc2_),
      Pmat_r(Pmat_r_), Pmu_r(Pmu_r_), base_B0_r(base_B0_r_), base_A_r(base_A_r_),
      beta_out_r(beta_out_r_), disp_out_r(disp_out_r_), iters_out_r(iters_out_r_), weight_out_r(weight_out_r_) {}
  
  // --- Parallel Loop ---
  void operator()(std::size_t begin, std::size_t end);
};
  
  
// --- rIndepNormalGammaReg_worker implementation ---
void rIndepNormalGammaReg_worker::operator()(std::size_t begin, std::size_t end) {
  const int l2 = x_r.nrow();
  const int l1 = x_r.ncol();

  for (std::size_t i = begin; i < end; ++i) {
    // Thread-local buffers and views (no shared state)
    std::vector<double> out_buf(static_cast<std::size_t>(l1), 0.0);
    RcppParallel::RMatrix<double> out_row(out_buf.data(), 1,  l1);  // 1×l1
    RcppParallel::RMatrix<double> out_col(out_buf.data(), l1, 1);   // l1×1

    std::vector<double> theta_buf(static_cast<std::size_t>(l1), 0.0);
    RcppParallel::RMatrix<double> theta_row(theta_buf.data(), 1,  l1); // 1×l1
    RcppParallel::RMatrix<double> theta_col(theta_buf.data(), l1, 1);  // l1×1

    std::vector<double> cbars_col_buf(static_cast<std::size_t>(l1), 0.0);
    RcppParallel::RMatrix<double> cbars_small_col(cbars_col_buf.data(), l1, 1); // l1×1


    // Scaled weights: classical logic requires wt2 = wt / dispersion before likelihood
    //   Rcpp::NumericVector wt2_nv(l2);                  // thread-local
    //  RcppParallel::RVector<double> wt2_r_old(wt2_nv);     // matches f2_gaussian_rmat signature


    std::vector<double> wt2_buf(static_cast<std::size_t>(l2), 0.0);
    // Wrap as a 1‑column matrix (l1 × 1)
    RcppParallel::RMatrix<double> wt2_r(wt2_buf.data(), l2, 1);


    iters_out_r[i]  = 1.0;
    weight_out_r[i] = 1.0;

    int accept = 0;



    while (accept == 0) {
      // 1) Slice/component selection via PLSD
      double U = runif_safe();
      int J_idx = 0;
      double U_left = U;

      while (true) {
        if (U_left <= PLSD_r[J_idx]) break;
        U_left -= PLSD_r[J_idx];
        ++J_idx;
      }

      // 2) Draw truncated-normal beta row
      for (int j = 0; j < l1; ++j) {
        out_row(0, j) = rnorm_ct(
          logrt_r(J_idx, j),
          loglt_r(J_idx, j),
          -cbars_r(J_idx, j),
          1.0
        );
      }

      // 3) Draw dispersion
      

      // Now call the sampler
      double dispersion = rinvgamma_ct_safe(shape3, rate2, disp_upper, disp_lower);
      
      // DEBUG: print the result

      // Guard
      if (!std::isfinite(dispersion)) {
        Rcpp::stop("Error: rinvgamma_ct_safe returned non-finite dispersion.");
      }      
      // 4) Solve theta (strict row-only)
      for (int j = 0; j < l1; ++j) cbars_small_col(j, 0) = cbars_r(J_idx, j);


      arma::mat theta_sol = Inv_f3_with_disp_rmat(Pmat_r, Pmu_r, base_B0_r, base_A_r,
                                                     dispersion, cbars_small_col);



      for (int j = 0; j < l1; ++j) {
        theta_row(0, j) = theta_sol(0, j);
      }

      // 5) Scale weights before likelihood
      for (int r = 0; r < l2; ++r) {
        wt2_r(r, 0) = wt_r[r] / dispersion;
      }



#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
      tbb::mutex::scoped_lock lock(f2_mutex);
#endif

      // Rcout << "Entering f2_gaussian_rmat_mat 1"  << std::endl;

      // 6) Likelihood calls (column views, pre-scaled weights)
      double LL_New2_scalar =-f2_gaussian_rmat_mat(theta_col, y_r, x_r, mu_r, P_r, alpha_r, wt2_r, 0)[0];


        double LL_Test_scalar =-f2_gaussian_rmat_mat(out_col,   y_r, x_r, mu_r, P_r, alpha_r, wt2_r, 0)[0];



        // 7) Upper bounds
        double U2     = runif_safe();
        double log_U2 = std::log(U2);

        double UB1 = LL_New2_scalar;
        for (int j = 0; j < l1; ++j)
          UB1 -= cbars_r(J_idx, j) * (out_row(0, j) - theta_row(0, j));

        double quad_sum = 0.0;
        for (int r = 0; r < l2; ++r) {
          double x_theta = 0.0;
          for (int c = 0; c < l1; ++c) x_theta += x_r(r, c) * theta_row(0, c);
          double resid  = (y_r[r] - alpha_r[r] - x_theta);
          double scaled = resid * std::sqrt(wt_r[r]);
          quad_sum += scaled * scaled;
        }
        double UB2 = 0.5 * (1.0 / dispersion) * (quad_sum - RSS_Min);
        UB2 -= UB2min_r[J_idx];

        double theta_P_theta = 0.0;
        for (int r = 0; r < l1; ++r) {
          double acc = 0.0;
          for (int c = 0; c < l1; ++c) acc += P_r(r, c) * theta_row(0, c);
          theta_P_theta += theta_row(0, r) * acc;
        }
        double c_theta = 0.0;
        for (int j = 0; j < l1; ++j) c_theta += cbars_r(J_idx, j) * theta_row(0, j);
        double New_LL_J = -0.5 * theta_P_theta + c_theta;

        double UB3A = lg_prob_factor_r[J_idx] + lmc1 + lmc2 * dispersion - New_LL_J;
        double New_LL_log_disp = lm_log1 + lm_log2 * std::log(dispersion);

        // Old UB3B definition
        
       double UB3B = (max_New_LL_UB - max_LL_log_disp + New_LL_log_disp)
         - (lmc1 + lmc2 * dispersion);

        
        // 3B: dispersion tilt bound (new definition)
        // double UB3B =
        // (lm_log1 + lm_log2 * std::log(dispersion))
        //   - (lmc1   + lmc2   * dispersion);
        
        double test1 = (LL_Test_scalar - UB1);
        double test  = test1 - (UB2 + UB3A + UB3B);
        test = test - log_U2;
        
        // Sanity checks: all must satisfy their sign constraints
        bool bad = false;
        std::ostringstream msg;
        
        if (test1 > 0.0) {
          bad = true;
          msg << "Sign violation: test1 = " << test1 << " > 0\n";
        }
        if (UB2 < 0.0) {
          bad = true;
          msg << "Sign violation: UB2 = " << UB2 << " < 0\n";
        }
        if (UB3A < 0.0) {
          bad = true;
          msg << "Sign violation: UB3A = " << UB3A << " < 0\n";
        }
        if (UB3B < 0.0) {
          bad = true;
          msg << "Sign violation: UB3B = " << UB3B << " < 0\n";
        }
        
        if (bad) {
          // Provide context for debugging
          msg << "Dispersion=" << dispersion
              << " LL_Test=" << LL_Test_scalar
              << " UB1=" << UB1
              << " UB2=" << UB2
              << " UB3A=" << UB3A
              << " UB3B=" << UB3B
              << " test=" << test;
          // Stop execution with informative error
          throw std::runtime_error(msg.str());
          
        }
        
        
        

        // Rcout << "Entering Output assignment"  << std::endl;

        // 8) Record outputs and accept/reject
        disp_out_r[i] = dispersion;
        for (int j = 0; j < l1; ++j) beta_out_r(i, j) = out_row(0, j);

        if (test >= 0.0) {
          accept = 1;
        } else {
          iters_out_r[i] = iters_out_r[i] + 1;
        }
        
    } // while (accept == 0)
  }   // for i
}



namespace glmbayes {

namespace sim {

// 414--> 839
Rcpp::List  rIndepNormalGammaReg_std(int n,NumericVector y,NumericMatrix x,
                                             NumericMatrix mu, /// This is typically standardized to be a zero vector
                                             NumericMatrix P, /// Part of prior precision shifted to the likelihood
                                             NumericVector alpha,NumericVector wt,
                                             Function f2,Rcpp::List  Envelope,
                                             Rcpp::List  gamma_list,
                                             Rcpp::List  UB_list,
                                             Rcpp::CharacterVector   family,Rcpp::CharacterVector   link,
                                             bool progbar,
                                            bool verbose
)
{

  
  int l1 = mu.nrow();
  int l2 = x.nrow();
  
  
  // Get various inputs frm the provided lists
  
  double shape3 =gamma_list["shape3"];
  double rate2 =gamma_list["rate2"];
  double disp_upper =gamma_list["disp_upper"];
  double disp_lower =gamma_list["disp_lower"];
  // double RSS_ML =UB_list["RSS_ML"];
  double max_New_LL_UB =UB_list["max_New_LL_UB"];
  double max_LL_log_disp =UB_list["max_LL_log_disp"];
  double lm_log1 =UB_list["lm_log1"];
  double lm_log2 =UB_list["lm_log2"];
  double lmc1 =UB_list["lmc1"];
  double lmc2 =UB_list["lmc2"];
  NumericVector lg_prob_factor =UB_list["lg_prob_factor"];
  NumericMatrix cbars=Envelope["cbars"];
  

  NumericVector iters_out(n);
  NumericVector disp_out(n);
  NumericVector weight_out(n);
  NumericMatrix beta_out(n,l1);
  double dispersion;
  NumericVector wt2(l1);
  
  
  arma::vec wt1b(wt.begin(), x.nrow());
  
  NumericMatrix thetabars_new(1,cbars.ncol());

  
  arma::mat thetabarsb_new(thetabars_new.begin(), thetabars_new.nrow(), thetabars_new.ncol(), false);

  arma::vec y2(y.begin(),l2);
  arma::vec alpha2(alpha.begin(),l2);
  arma::mat x2(x.begin(),l2,l1);
  arma::mat P2(P.begin(),l1,l1);
  
  double UB1;
  double UB2;
  double UB3A;
  double UB3B;
  double New_LL_log_disp;
    
  int a1=0;
  double test1=0;
  double test=0;
  NumericVector J(n);
  // NumericVector draws(n);
  NumericMatrix out(1,l1);
  double a2=0;
  double U=0;
  double U2=0;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
  
  double RSS_Min=UB_list["RSS_Min"];
  NumericVector UB2min=UB_list["UB2min"];
  
  
  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg_std:Inv_f3_precompute_disp] Entering: "
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  

  
  // Build cache once outside the loop
  Rcpp::List cache = Inv_f3_precompute_disp(cbars, y, x, mu, P, alpha, wt);
  
  
  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg_std:Inv_f3_precompute_disp] Exiting: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  

  
  for(int i=0;i<n;i++){

    Rcpp::checkUserInterrupt();
    
   if(progbar==1){
     // progress_bar3(i, n-1);
     progress_bar(i, n-1);
     
     if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
   }
    

    // Rcpp::Rcout << "[rIndepNormalGammaReg_std] Entering accept/reject: \n";
   
    
    a1=0;
    iters_out[i]=1;  
    while(a1==0){

          
      
      // Simulate from discrete distribution
      
      U=R::runif(0.0, 1.0);
      a2=0;
      J(0)=0;    
      while(a2==0){
        if(U<=PLSD(J(0))) a2=1;
        if(U>PLSD(J(0))){ 
          U=U-PLSD(J(0));
          J(0)=J(0)+1;
          
        }
      }
      

            
      // Simulate for beta
      
      for(int j=0;j<l1;j++){  out(0,j)=rnorm_ct(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);          }
      
      dispersion=rinvgamma_ct_safe(shape3,rate2,disp_upper,disp_lower);
      
      
      wt2=wt/dispersion;
      NumericMatrix cbars_small = cbars( Range(J(0),J(0)) , Range(0,cbars.ncol()-1) );
      
      // Compute Adjusted theta (accounting for changed dispersion) - New tangency points
    
      arma::mat theta2 = Inv_f3_with_disp(cache, dispersion, transpose(cbars_small));
      thetabarsb_new = theta2;
      

      // Recompoute LL at the new gradient point
      NumericVector LL_New2=-f2_gaussian(transpose(thetabars_new),  y, x, mu, P, alpha, wt2);  
      
    
      
      U2=R::runif(0.0, 1.0);
      
      double log_U2=log(U2);
      NumericVector J_out=J;
      NumericVector b_out=out(0,_);
      arma::rowvec b_out2(b_out.begin(),l1,false);
      NumericVector thetabars_temp=thetabars_new(0,_); // Changed
      
      arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
      NumericVector cbars_temp=cbars(J_out(0),_);
      arma::vec  cbars_temp2(cbars_temp.begin(), l1);

      NumericVector LL_Test=-f2_gaussian(transpose(out),  y, x, mu, P, alpha, wt2);

      // Block 1: UB1 
      //   Same form as in fixed dispersion case but thetabar is a function of the dispersion
      //   So all components that include thetabar must now be bounded as well
      
      arma::colvec betadiff=trans(b_out2)-thetabars_temp2;
      UB1=LL_New2(0) -arma::as_scalar(trans(cbars_temp2)*betadiff);
      
      //Block 2: UB2 [RSS Term bounded by shifting it to the gamma candidate]
      
      
      arma::colvec yxbeta=(y2-alpha2-x2*thetabars_temp2)%sqrt(wt1b); 

      
      // Extract the current cbars row as a NumericVector
      // NumericVector cbars_j = cbars(J_out(0), _);
      

    
      // Continue with existing UB2 calculation without modification
      UB2=0.5*(1.0/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_Min);
      
      // Subtract UB2min --> Should improve acceptance
      
      UB2=UB2-UB2min[J_out(0)];
      

      // Compute g1_j(d) for the chosen face j = J_out(0)

      
      
      int j = J_out(0);
      double g1j = g1_face_at_disp(
        dispersion,   // current dispersion draw
        j,            // face index
        cache,        // cached matrices for Inv_f3_with_disp
        P2,           // precision matrix
        cbars         // matrix of c_j rows
      );
      

      // Modified UB3A 
      
      UB3A= lg_prob_factor(J_out(0))+lmc1+lmc2*dispersion-g1j;
      
    
      New_LL_log_disp=lm_log1+lm_log2*log(dispersion);
      UB3B=(max_New_LL_UB-max_LL_log_disp+New_LL_log_disp)-(lmc1+lmc2*dispersion);

      test1=LL_Test[0]-UB1;
        
      test= test1-(UB2+UB3A+UB3B);  // Should be all negative 
      test = test - log_U2;

    
      // Sanity checks: all must satisfy their sign constraints
      bool bad = false;
      std::ostringstream msg;
      
      if (test1 > 0.0) {
        bad = true;
        msg << "Sign violation: test1 = " << test1 << " > 0\n";
      }
      if (UB2 < 0.0) {
        bad = true;
        msg << "Sign violation: UB2 = " << UB2 << " < 0\n";
      }
      if (UB3A < 0.0) {
        bad = true;
        msg << "Sign violation: UB3A = " << UB3A << " < 0\n";
      }

            
      if (UB3B < 0.0) {
        bad = true;
        msg << "Sign violation: UB3B = " << UB3B << " < 0\n";
      }
    

      
      if (bad) {
        // Provide context for debugging
        msg << "Dispersion=" << dispersion
            << " LL_Test=" << LL_Test[0]
            << " UB1=" << UB1
            << " UB2=" << UB2
            << " UB3A=" << UB3A
            << " UB3B=" << UB3B
            << " lg_prob_factor[]=" << lg_prob_factor[J_out(0)]
            << " g1j=" << g1j
            << " test=" << test;
        
        // Stop execution with informative error
        throw std::runtime_error(msg.str());
        
      }
      
      
      
      disp_out[i] = dispersion;
      beta_out(i, _) = out(0, _);
      

      if(test>=0){

        
        a1=1;
        
      }
      else{
        iters_out[i]=iters_out[i]+1;
        }    
      

    }  
    
    
  }
  

  return Rcpp::List::create(Rcpp::Named("beta_out")=beta_out,Rcpp::Named("disp_out")=disp_out,
                            Rcpp::Named("iters_out")=iters_out,Rcpp::Named("weight_out")=weight_out);  
  
  
}



Rcpp::List  rIndepNormalGammaReg_std_face(int n,NumericVector y,NumericMatrix x,
                                     NumericMatrix mu, /// This is typically standardized to be a zero vector
                                     NumericMatrix P, /// Part of prior precision shifted to the likelihood
                                     NumericVector alpha,NumericVector wt,
                                     Function f2,Rcpp::List  Envelope,
                                     Rcpp::List  gamma_list,
                                     Rcpp::List  UB_list,
                                     Rcpp::CharacterVector   family,Rcpp::CharacterVector   link,
                                     bool progbar,
                                     bool verbose
)
{

  Rcpp::Rcout << "[rIndepNormalGammaReg_std] Face version: \n";
  
  // 1. Grab the base environment
  Rcpp::Environment base = Rcpp::Environment::base_env();
  
  // 2. Pull out the 'interactive' function
  Rcpp::Function interactive = base["interactive"];
  
  
  int l1 = mu.nrow();
  int l2 = x.nrow();
  
  
  // Get various inputs frm the provided lists
  
  double shape3 =gamma_list["shape3"];
  double rate2 =gamma_list["rate2"];
  double disp_upper =gamma_list["disp_upper"];
  double disp_lower =gamma_list["disp_lower"];
  // double RSS_ML =UB_list["RSS_ML"];
  double max_New_LL_UB =UB_list["max_New_LL_UB"];
  double max_LL_log_disp =UB_list["max_LL_log_disp"];
  double lm_log1 =UB_list["lm_log1"];
  double lm_log2 =UB_list["lm_log2"];
  double lmc1 =UB_list["lmc1"];
  double lmc2 =UB_list["lmc2"];
  NumericVector lg_prob_factor =UB_list["lg_prob_factor"];
  NumericMatrix cbars=Envelope["cbars"];
  
 
  
  // --- New: face-specific UB3A/UB3B and Gamma proposal geometry ---
  
  
  // NEW: extract anchor geometry
  double d1_star = UB_list["d1_star"];        // UB3A anchor
  //  double d2_star = UB_list["d2_star"];        // UB3B anchor
  // NumericVector thetabar_const_base = UB_list["thetabar_const_base"];  // should be g1(d1*)
  // NumericVector New_LL_Slope       = UB_list["New_LL_Slope"];          // should g1'(d1*)
  
  
  NumericVector shape3_face  = gamma_list["shape3_face"];   // per-face Gamma shapes
  
  
  NumericVector iters_out(n);
  NumericVector disp_out(n);
  NumericVector weight_out(n);
  NumericMatrix beta_out(n,l1);
  double dispersion;
  NumericVector wt2(l1);
  
  
  arma::vec wt1b(wt.begin(), x.nrow());
  
  
  NumericMatrix cbarst(cbars.ncol(),cbars.nrow());
  NumericMatrix thetabars(cbars.nrow(),cbars.ncol());
  NumericMatrix thetabars_new(1,cbars.ncol());
  
  NumericVector New_LL(cbars.nrow());
  
  
  
  
  arma::mat cbarsb(cbars.begin(), cbars.nrow(), cbars.ncol(), false);
  arma::mat cbarstb(cbarst.begin(), cbarst.nrow(), cbarst.ncol(), false);
  
  arma::mat thetabarsb(thetabars.begin(), thetabars.nrow(), thetabars.ncol(), false);
  arma::mat thetabarsb_new(thetabars_new.begin(), thetabars_new.nrow(), thetabars_new.ncol(), false);
  cbarstb=trans(cbarsb);
  
  arma::vec y2(y.begin(),l2);
  arma::vec alpha2(alpha.begin(),l2);
  arma::mat x2(x.begin(),l2,l1);
  arma::mat P2(P.begin(),l1,l1);
  
  double UB1;
  double UB2;
  double UB3A;
  double UB3B;
  double New_LL_log_disp;
  
  double UB3A_face;
  
  
  int a1=0;
  double test1=0;
  double test=0;
  NumericVector J(n);
  // NumericVector draws(n);
  NumericMatrix out(1,l1);
  double a2=0;
  double U=0;
  double U2=0;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
  
  double RSS_Min=UB_list["RSS_Min"];
  NumericVector UB2min=UB_list["UB2min"];
  
  // Adding to enable face specific dispersion bounds
  
  NumericVector thetabar_const_base=UB_list["thetabar_const_base"];
  NumericVector New_LL_Slope=UB_list["New_LL_Slope"];
  
  
  //  NumericVector ub2_min=;
  
  
  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg_std:Inv_f3_precompute_disp] Entering: "
                << glmbayes::progress::timestamp_cpp()
                << "\n";
  }
  
  
  
  
  // Build cache once outside the loop
  Rcpp::List cache = Inv_f3_precompute_disp(cbars, y, x, mu, P, alpha, wt);
  
  
  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg_std:Inv_f3_precompute_disp] Exiting: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  
  
  for(int i=0;i<n;i++){
    
    Rcpp::checkUserInterrupt();
    
    if(progbar==1){
      // progress_bar3(i, n-1);
      progress_bar(i, n-1);
      
      if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
    }
    
    // 3. Test progbar *and* interactive()
    
    
    // Rcpp::Rcout << "[rIndepNormalGammaReg_std] Entering accept/reject: \n";
    
    
    a1=0;
    iters_out[i]=1;  
    while(a1==0){
      
      // Simulate from discrete distribution
      
      U=R::runif(0.0, 1.0);
      a2=0;
      J(0)=0;    
      while(a2==0){
        if(U<=PLSD(J(0))) a2=1;
        if(U>PLSD(J(0))){ 
          U=U-PLSD(J(0));
          J(0)=J(0)+1;
          
        }
      }
      
      
      
      // Simulate for beta
      
      for(int j=0;j<l1;j++){  out(0,j)=rnorm_ct(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);          }
      
      
      
      // Update this to make distribution contingent on component of the grid
      
      //      dispersion=rinvgamma_ct(shape3,rate2,disp_upper,disp_lower);
//      dispersion=rinvgamma_ct_safe(shape3,rate2,disp_upper,disp_lower);
      dispersion=rinvgamma_ct_safe(shape3_face[J(0)],rate2,disp_upper,disp_lower);
      
      
      wt2=wt/dispersion;
      NumericMatrix cbars_small = cbars( Range(J(0),J(0)) , Range(0,cbars.ncol()-1) );
      
      // Compute Adjusted theta (accounting for changed dispersion) - New tangency points
      
      arma::mat theta2 = Inv_f3_with_disp(cache, dispersion, transpose(cbars_small));
      thetabarsb_new = theta2;
      
      
      // Recompoute LL at the new gradient point
      NumericVector LL_New2=-f2_gaussian(transpose(thetabars_new),  y, x, mu, P, alpha, wt2);  
      
      
      
      U2=R::runif(0.0, 1.0);
      
      double log_U2=log(U2);
      NumericVector J_out=J;
      NumericVector b_out=out(0,_);
      arma::rowvec b_out2(b_out.begin(),l1,false);
      NumericVector thetabars_temp=thetabars_new(0,_); // Changed
      
      arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
      NumericVector cbars_temp=cbars(J_out(0),_);
      arma::vec  cbars_temp2(cbars_temp.begin(), l1);
      
      
      
      NumericVector LL_Test=-f2_gaussian(transpose(out),  y, x, mu, P, alpha, wt2);
      
      
      
      // Block 1: UB1 
      //   Same form as in fixed dispersion case but thetabar is a function of the dispersion
      //   So all components that include thetabar must now be bounded as well
      
      arma::colvec betadiff=trans(b_out2)-thetabars_temp2;
      UB1=LL_New2(0) -arma::as_scalar(trans(cbars_temp2)*betadiff);
      
      //Block 2: UB2 [RSS Term bounded by shifting it to the gamma candidate]
      
      
      arma::colvec yxbeta=(y2-alpha2-x2*thetabars_temp2)%sqrt(wt1b); 
      
      
      // Extract the current cbars row as a NumericVector
      // NumericVector cbars_j = cbars(J_out(0), _);
      
      
      
      // Continue with existing UB2 calculation without modification
      UB2=0.5*(1.0/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_Min);
      
      // Subtract UB2min --> Should improve acceptance
      
      UB2=UB2-UB2min[J_out(0)];
      
      
      // Block 3: UB3A (adjusts because probabilities of components in grid are different from original grid)
      // Investigate whether changing probabilities of grid components for proposal
      // allows us to do away with this term and to thereby improve the acceptance rate
      
      // This is likely time consuming part
      
      
      // Compute g1_j(d) for the chosen face j = J_out(0)
      
      
      
      int j = J_out(0);
      double g1j = g1_face_at_disp(
        dispersion,   // current dispersion draw
        j,            // face index
        cache,        // cached matrices for Inv_f3_with_disp
        P2,           // precision matrix
        cbars         // matrix of c_j rows
      );
      
      
      double g2j=g2_face_at_disp(
        dispersion,   // current dispersion draw
        j,            // face index
        cache,        // cached matrices for Inv_f3_with_disp
        P2,           // precision matrix
        cbars,         // matrix of c_j rows
        d1_star,
        New_LL_Slope
      );
      
      // Store it exactly where New_LL(j) was previously used
      New_LL(j) = g1j;
      
      // Modified UB3A 
      
//      UB3A= lg_prob_factor(J_out(0))+lmc1+lmc2*dispersion-New_LL(J_out(0));
      UB3A= lg_prob_factor(J_out(0))+lmc1+lmc2*dispersion-New_LL(J_out(0));
      
      
      // UB3A_face =thetabar_const_upp_apprx[J_out(0)]+ lmc1_face[J_out(0)] + lmc2_face[J_out(0)] * dispersion - New_LL(J_out(0));
      //double UB3A = lg_prob_factor_r[J_idx] + lmc1 + lmc2 * dispersion - New_LL_J;
      
      UB3A_face=g2j-g1j;
      
      // Block 4: UB3B  
      
      //--------------------------------------------------------
      // 4. UB3B slope (log-tilt slope)
      //    A7-compatible slope relation:
      //        lm_log2_j = slope_j * d2
      
      // where d2 = (upp - low) / log(upp/low)
      
      double lmc1_j = g1j - New_LL_Slope[j] * d1_star;
      double lmc2_j=New_LL_Slope[j];
      
      double lm_log2_j = New_LL_Slope[j] *(disp_upper-disp_lower)/log(disp_upper/disp_lower);
      
      double lm_log1_j =lmc1_j + New_LL_Slope[j] *disp_lower-lm_log2_j*(log(disp_lower)) ;
      
      //    ⇒ lm_log1_j = lmc1_j + slope_j * low
      //                   - lm_log2_j * log(low)
      
      
      double UB3B_face=lm_log1_j+lm_log2_j*log(dispersion)-(lmc1_j+lmc2_j*dispersion);
      
      
      
      // Old UB3B definition
      
      
      
      New_LL_log_disp=lm_log1+lm_log2*log(dispersion);
      
      UB3B=(max_New_LL_UB-max_LL_log_disp+New_LL_log_disp)-(lmc1+lmc2*dispersion);
      
      
      
      test1=LL_Test[0]-UB1;
      
      test= test1-(UB2+UB3A+UB3B);  // Should be all negative 
      
      
      test = test - log_U2;
      
      
      double test_face=test1-(UB2+UB3A_face+UB3B_face)- log_U2;
      
      
      // Sanity checks: all must satisfy their sign constraints
      bool bad = false;
      std::ostringstream msg;
      
      if (test1 > 0.0) {
        bad = true;
        msg << "Sign violation: test1 = " << test1 << " > 0\n";
      }
      if (UB2 < 0.0) {
        bad = true;
        msg << "Sign violation: UB2 = " << UB2 << " < 0\n";
      }
      if (UB3A < 0.0) {
        bad = true;
        msg << "Sign violation: UB3A = " << UB3A << " < 0\n";
      }
      if (UB3A_face < 0.0) {
        bad = true;
        msg << "Sign violation: UB3A_face = " << UB3A_face << " < 0\n";
      }
      
      
      if (UB3B < 0.0) {
        bad = true;
        msg << "Sign violation: UB3B = " << UB3B << " < 0\n";
      }
      
      if (UB3B_face < 0.0) {
        bad = true;
        msg << "Sign violation: UB3B_face = " << UB3B_face << " < 0\n";
      }
      
      
      if (bad) {
        // Provide context for debugging
        msg << "Dispersion=" << dispersion
            << " LL_Test=" << LL_Test[0]
            << " UB1=" << UB1
            << " UB2=" << UB2
            << " UB3A=" << UB3A
            << " UB3A_face=" << UB3A_face
            << " UB3B=" << UB3B
            << " UB3B_face=" << UB3B_face
        // << " thetabar_const_upp_apprx=" << thetabar_const_upp_apprx[J_out(0)]
        // << " thetabar_const_low_apprx=" << thetabar_const_low_apprx[J_out(0)]
           << " lg_prob_factor[]=" << lg_prob_factor[J_out(0)]
           << " lg_prob_factor()=" << lg_prob_factor(J_out(0))
           << " New_LL(j)=" << New_LL(j)
           << " g1j=" << g1j
           << " test=" << test;
        
        // Stop execution with informative error
        throw std::runtime_error(msg.str());
        
      }
      
      
      
      disp_out[i] = dispersion;
      beta_out(i, _) = out(0, _);
      
      
      //      if(test>=0){
      if(test_face>0){  
        
        
        a1=1;
        
      }
      else{
        iters_out[i]=iters_out[i]+1;
      }    
      
      
    }  
    
    
  }
  
  // Temporarily just return non-sense constants equal to all 1
  
  return Rcpp::List::create(Rcpp::Named("beta_out")=beta_out,Rcpp::Named("disp_out")=disp_out,
                            Rcpp::Named("iters_out")=iters_out,Rcpp::Named("weight_out")=weight_out);  
  
  
  
}




Rcpp::List rIndepNormalGammaReg_std_parallel(
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
    bool progbar ,
    bool verbose 
) {


    // Base env (kept as-is)
  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function interactive = base["interactive"];

  const int l1 = mu.nrow();
  // const int l2 = x.nrow();

  // Scalars from lists
  double shape3          = gamma_list["shape3"];
  double rate2           = gamma_list["rate2"];
  double disp_upper      = gamma_list["disp_upper"];
  double disp_lower      = gamma_list["disp_lower"];
  // double RSS_ML          = UB_list["RSS_ML"];
  double max_New_LL_UB   = UB_list["max_New_LL_UB"];
  double max_LL_log_disp = UB_list["max_LL_log_disp"];
  double lm_log1         = UB_list["lm_log1"];
  double lm_log2         = UB_list["lm_log2"];
  double lmc1            = UB_list["lmc1"];
  double lmc2            = UB_list["lmc2"];

  // --- New: face-specific UB3A/UB3B and Gamma proposal geometry ---
  // NumericVector lmc1_face    = gamma_list["lmc1_face"];     // UB3A intercepts
  // NumericVector lmc2_face    = gamma_list["lmc2_face"];     // UB3A slopes
  // NumericVector lm_log1_face = gamma_list["lm_log1_face"];  // UB3B intercepts
  // NumericVector lm_log2_face = gamma_list["lm_log2_face"];  // UB3B slopes
  // NumericVector shape3_face  = gamma_list["shape3_face"];   // per-face Gamma shapes
  
  Rcpp::NumericVector lg_prob_factor = UB_list["lg_prob_factor"];
  Rcpp::NumericMatrix cbars          = Envelope["cbars"];
  Rcpp::NumericVector PLSD           = Envelope["PLSD"];
  Rcpp::NumericMatrix loglt          = Envelope["loglt"];
  Rcpp::NumericMatrix logrt          = Envelope["logrt"];
  double RSS_Min                     = UB_list["RSS_Min"];
  Rcpp::NumericVector UB2min         = UB_list["UB2min"];

  // Outputs
  Rcpp::NumericVector iters_out(n);
  Rcpp::NumericVector disp_out(n);
  Rcpp::NumericVector weight_out(n);
  Rcpp::NumericMatrix beta_out(n, l1);

  // Build cache once outside the loop


  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:Inv_f3_precompute_disp] Entering: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  
  Rcpp::List cache = Inv_f3_precompute_disp(cbars, y, x, mu, P, alpha, wt);


  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:Inv_f3_precompute_disp] Exiting: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  
    Rcpp::NumericMatrix Pmat_nm    = cache["Pmat"];
  Rcpp::NumericMatrix Pmu_nm     = cache["Pmu"];
  Rcpp::NumericVector base_B0_nv = cache["base_B0"];
  Rcpp::NumericMatrix base_A_nm  = cache["base_A"];

  // Wrap outputs with RcppParallel views
  RcppParallel::RMatrix<double> beta_out_r(beta_out);
  RcppParallel::RVector<double> disp_out_r(disp_out);
  RcppParallel::RVector<double> iters_out_r(iters_out);
  RcppParallel::RVector<double> weight_out_r(weight_out);

  // Wrap inputs with RcppParallel views
  RcppParallel::RVector<double> y_r(y);
  RcppParallel::RMatrix<double> x_r(x);
  RcppParallel::RMatrix<double> mu_r(mu);
  RcppParallel::RMatrix<double> P_r(P);
  RcppParallel::RVector<double> alpha_r(alpha);
  RcppParallel::RVector<double> wt_r(wt);

  RcppParallel::RMatrix<double> cbars_r(cbars);
  RcppParallel::RVector<double> PLSD_r(PLSD);
  RcppParallel::RMatrix<double> loglt_r(loglt);
  RcppParallel::RMatrix<double> logrt_r(logrt);

  RcppParallel::RVector<double> lg_prob_factor_r(lg_prob_factor);
  RcppParallel::RVector<double> UB2min_r(UB2min);

  RcppParallel::RMatrix<double> Pmat_r(Pmat_nm);
  RcppParallel::RMatrix<double> Pmu_r(Pmu_nm);
  RcppParallel::RVector<double> base_B0_r(base_B0_nv);
  RcppParallel::RMatrix<double> base_A_r(base_A_nm);



  // Construct worker
  rIndepNormalGammaReg_worker worker(
      n,
      y_r, x_r, mu_r, P_r, alpha_r, wt_r,
      cbars_r, PLSD_r, loglt_r, logrt_r,
      lg_prob_factor_r, UB2min_r,
      shape3, rate2, disp_upper, disp_lower,
      RSS_Min, max_New_LL_UB, max_LL_log_disp,
      lm_log1, lm_log2, lmc1, lmc2,
      Pmat_r, Pmu_r, base_B0_r, base_A_r,
      beta_out_r, disp_out_r, iters_out_r, weight_out_r
  );

  
  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:Simulation:Pilot] Entering: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  

  // --- Single-draw test (serial) ---
  int m_test = 1;
  auto t0 = std::chrono::steady_clock::now();
  worker(0, m_test);  // run worker serially for 1 observation
  auto t1 = std::chrono::steady_clock::now();
  double elapsed_test_sec = std::chrono::duration<double>(t1 - t0).count();

  if (verbose) Rcpp::Rcout << "[rIndepNormalGammaReg:Simulation:Pilot] Completed in " << elapsed_test_sec << "s.\n";

  // --- Conservative calibration sizing (time-bounded) ---
  // Use single test to bound worst-case per-observation time in ms
  double per_obs_ms_serial = elapsed_test_sec * 1000.0 / std::max(1, m_test);
  
  // Aim for ~1% of n but cap at ~5 minutes worst-case based on serial bound
  int m1 = std::max(1, (int)std::ceil(0.01 * (double)n));
  int m2 = std::max(1, (int)std::floor(300000.0 / std::max(1.0, per_obs_ms_serial))); // 300k ms ≈ 5 min
  int m_stage = std::min(m1, m2);   // <-- defined here, before use
  
  if (verbose) {Rcpp::Rcout << "[rIndepNormalGammaReg:Simulation:Calibration] Using " << m_stage
              << " observations at "
              << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")()))
              << "\n";}
  
  // --- Calibration run (parallel) ---
  auto t_cal0 = std::chrono::steady_clock::now();
  RcppParallel::parallelFor(0, m_stage, worker);   // <-- now m_stage is defined
  auto t_cal1 = std::chrono::steady_clock::now();
  double cal_elapsed_sec = std::chrono::duration<double>(t_cal1 - t_cal0).count();
  
  // Per‑observation cost from calibration
  double per_obs_sec   = cal_elapsed_sec / std::max(1.0, (double)m_stage);
  double est_total_sec = per_obs_sec * (double)n;
  
  // Diagnostics
  if (verbose){  Rcpp::Rcout << "[rIndepNormalGammaReg:Simulation:Calibration] Completed in " << cal_elapsed_sec
              << " s for " << m_stage << " observations.\n";
  Rcpp::Rcout << "[rIndepNormalGammaReg:Simulation:Calibration] per_obs_sec = " << per_obs_sec
              << " s; estimated total = " << est_total_sec << " s\n";}
  
  auto fmt_hms = [](double seconds) {
    long long s = static_cast<long long>(std::round(seconds));
    long long h = s / 3600; s %= 3600;
    long long m = s / 60;   s %= 60;
    std::ostringstream oss;
    if (h > 0) oss << h << " h  ";
    if (m > 0 || h > 0) oss << m << " m  ";
    oss << s << " s";
    return oss.str();
  };
  
  if (verbose) {Rcpp::Rcout << "[rIndepNormalGammaReg:Simulation]  " << n << " observations: "
                            << fmt_hms(est_total_sec) << " (" << est_total_sec << " seconds).\n";}
  
  
  // --- Interactive safeguard if estimate exceeds 5 minutes ---
  if (est_total_sec > 300.0) {
    std::string prompt = "Estimated simulation exceeds 5 minutes. Continue? [y/N]: ";
    Rcpp::Function r_interactive("interactive");
    bool is_interactive = Rcpp::as<bool>(r_interactive());
    
    if (is_interactive) {
      Rcpp::Function readline("readline");
      while (true) {
        std::string ans = Rcpp::as<std::string>(readline(Rcpp::wrap(prompt)));
        // trim whitespace
        auto ltrim = [](std::string &s) {
          s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                          [](unsigned char ch){ return !std::isspace(ch); }));
        };
        auto rtrim = [](std::string &s) {
          s.erase(std::find_if(s.rbegin(), s.rend(),
                               [](unsigned char ch){ return !std::isspace(ch); }).base(), s.end());
        };
        ltrim(ans); rtrim(ans);
        
        if (ans == "y" || ans == "yes" || ans == "1" || ans == "continue") {
          Rcpp::Rcout << "[INFO] User chose to continue full run.\n";
          break; // proceed
        } else if (ans == "n" || ans == "no" || ans == "2" || ans.empty()) {
          Rcpp::Rcout << "[INFO] User declined. Stopping simulation.\n";
          Rcpp::stop("Simulation stopped by user after time estimate.");
        } else {
          Rcpp::Rcout << "Invalid input. Please enter y (continue) or N (stop).\n";
        }
      }
    } else {
      Rcpp::Rcout << "[NOTE] Non-interactive session: proceeding automatically.\n";
      Rcpp::Rcout << "[INFO] Proceeding with full run.\n";
    }
  }  

  Rcpp::Function fmt("format");
  Rcpp::Function systime("Sys.time");
  Rcpp::CharacterVector now = fmt(systime(), Rcpp::Named("format") = "%H:%M:%S");
  
   if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:Simulation]  Entering "
                << Rcpp::as<std::string>(now[0]) << " <<<\n";
   }
  
  // --- Capture start time ---
  double sim_start = Rcpp::as<double>(
    Rcpp::Function("as.numeric")(Rcpp::Function("Sys.time")())
  );
  
  // Parallel loop
  RcppParallel::parallelFor(0, n, worker);
  // worker(0,n);
  
  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:Simulation] Exiting "
                << Rcpp::as<std::string>(now[0]) << " <<<\n";
  }
  
  // --- Capture end time ---
  double sim_end = Rcpp::as<double>(
    Rcpp::Function("as.numeric")(Rcpp::Function("Sys.time")())
  );
  
  double sim_elapsed = sim_end - sim_start;
  int h_elapsed = static_cast<int>(sim_elapsed / 3600);
  int m_elapsed = static_cast<int>((sim_elapsed - h_elapsed*3600) / 60);
  int s_elapsed = static_cast<int>(sim_elapsed - h_elapsed*3600 - m_elapsed*60);
  
   if (verbose) {
    now = fmt(systime(), Rcpp::Named("format") = "%H:%M:%S");
    Rcpp::Rcout << "[rIndepNormalGammaReg:Simulation] Completed in: "
                << h_elapsed << " h  " << m_elapsed << " m  " << s_elapsed << " s.\n";
   }  



  return Rcpp::List::create(
    Rcpp::Named("beta_out")   = beta_out,
    Rcpp::Named("disp_out")   = disp_out,
    Rcpp::Named("iters_out")  = iters_out,
    Rcpp::Named("weight_out") = weight_out
  );
}





Rcpp::List rIndepNormalGammaReg(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector offset,
    Rcpp::NumericVector wt,
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
){
  
  
  int disp_grid_type=2;
  
  if(use_parallel) disp_grid_type=2;
  
  
  // Base R functions
  Rcpp::Function lm_wfit("lm.wfit");
  Rcpp::Function optim("optim");
  Rcpp::Function gaussian("gaussian");
  
  
  // glmbayes namespace function
  Rcpp::Environment glmbayes_ns = Rcpp::Environment::namespace_env("glmbayes");
  Rcpp::Function glmbfamfunc = glmbayes_ns["glmbfamfunc"];
  
  
  // Build y* = y - offset
  int n_obs = y.size();
  Rcpp::NumericVector ystar(n_obs);
  for (int i = 0; i < n_obs; i++) {
    ystar[i] = y[i] - offset[i];
  }
  
  double n_w = 0.0;
  for (int i = 0; i < wt.size(); ++i)   n_w += wt[i];
  
  
  // Call lm.wfit(X, y*, w)
  Rcpp::List fit = lm_wfit(
    Rcpp::_["x"] = x,
    Rcpp::_["y"] = ystar,
    Rcpp::_["w"] = wt
  );
  
  // Extract residuals
  Rcpp::NumericVector res = fit["residuals"];
  
  // Compute RSS
  double RSS = 0.0;
  for (int i = 0; i < res.size(); i++) {
    RSS += res[i] * res[i];
  }
  
  // Extract rank
  int p = Rcpp::as<int>(fit["rank"]);
  
  // Compute dispersion
  double dispersion2 = RSS / (n_obs - p);
  
  
  // Call glmbfamfunc(gaussian())
  Rcpp::List famfunc = glmbfamfunc( gaussian() );
  
  // Extract f2 and f3
  Rcpp::Function f2 = famfunc["f2"];
  Rcpp::Function f3 = famfunc["f3"];
  
  
  // Armadillo views
  arma::mat X   = Rcpp::as<arma::mat>(x);          // n_obs × p
  arma::vec Y   = Rcpp::as<arma::vec>(y);          // n_obs
  arma::rowvec y_row = Y.t();                      // 1 × n_obs
  arma::rowvec off_row = Rcpp::as<arma::rowvec>(offset); // 1 × n_obs
  arma::rowvec wt_row  = Rcpp::as<arma::rowvec>(wt);     // 1 × n_obs
  
  Rcpp::List cpp_out;   // declare here so it survives the loop
  double RSS_Post2 = NA_REAL;   // declare before the loop
  
  
  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:rNormalReg] Entering Loop: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  for (int j = 0; j < 10; ++j) {
    
    // --- Call rNormalReg (C++ version of .rnorm_reg) ---
    cpp_out = rNormalReg(
      10000,          // n
      y,              // y
      x,              // x
      mu,             // mu
      P,              // P
      offset,         // offset
      wt,             // wt
      dispersion2,    // dispersion
      f2,             // f2
      f3,             // f3
      mu,             // start
      "gaussian",     // family
      "identity",     // link
      Gridtype        // Gridtype
    );


    // Posterior draws: matrix (n_draws × p)
    arma::mat beta_draws = Rcpp::as<arma::mat>(cpp_out["coefficients"]);

    // lp_mat: n_draws × n_obs = beta_draws %*% t(x)
    arma::mat lp_mat = beta_draws * X.t();
    
    // eta_mat = lp_mat + offset (broadcasted by row)
    arma::mat eta_mat = lp_mat.each_row() + off_row;
    
    // For Gaussian identity link, mu_mat = eta_mat
    arma::mat mu_mat = eta_mat;
    
    // diff = mu_mat - y (broadcast y by row)
    arma::mat diff = mu_mat.each_row() - y_row;
    
    // elementwise square
    arma::mat res_sq = diff % diff;
    
    // weight each column by wt
    arma::mat res_sq_weighted = res_sq;
    res_sq_weighted.each_row() %= wt_row;
    
    // RSS_k = sum_i w_i (y_i - mu_ik)^2  (per draw)
    arma::vec RSS_temp = arma::sum(res_sq_weighted, 1);
    
    // Posterior mean of RSS
    RSS_Post2 = arma::mean(RSS_temp);
    
    // Mode from cpp_out
    arma::vec b_old = Rcpp::as<arma::vec>(cpp_out["coef.mode"]);
    
    // x %*% b_old
    arma::vec xbetastar = X * b_old;
    
    // RSS2_post = (y - X b_old)' (y - X b_old)
    arma::vec resid_ml = Y - xbetastar;
    //double RSS2_post = arma::dot(resid_ml, resid_ml);
    

    
        // Update shape, rate, dispersion
    double shape2 = shape + n_w / 2.0;
    double rate2  = rate  + RSS_Post2 / 2.0;
    
    dispersion2 = rate2 / (shape2 - 1.0);
  }
  
  
  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:rNormalReg] Exiting Loop: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  
  
  // Rcpp::Rcout << "dispersion2 = " << dispersion2 << "\n";
  
  
  // -------------------------------
  // Posterior Mode + Hessian Block
  // -------------------------------
  
  // Extract posterior mode from sampler
  arma::vec betastar = Rcpp::as<arma::vec>(cpp_out["coef.mode"]);
  double dispstar = dispersion2;
  
  // wt2 = wt / dispstar
  Rcpp::NumericVector wt2(n_obs);
  for (int i = 0; i < n_obs; ++i)    wt2[i] = wt[i] / dispstar;
  

  
  // alpha = X %*% mu + offset
  arma::vec alpha_vec = X * Rcpp::as<arma::vec>(mu) + Rcpp::as<arma::vec>(offset);
  Rcpp::NumericVector alpha = Rcpp::wrap(alpha_vec);
  
  // mu2 = 0 * mu
  Rcpp::NumericVector mu2(mu.size());
  for (int i = 0; i < mu.size(); ++i)
    mu2[i] = 0.0;
  
  // parin = mu - mu  (zero vector)
  Rcpp::NumericVector parin(mu.size());
  for (int i = 0; i < mu.size(); ++i)
    parin[i] = 0.0;
  
  // ---- Posterior Mode Optimization ----
  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:optim] Entering: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  // Call R's optim() directly
  Rcpp::List opt_out = optim(
    Rcpp::_["par"] = parin,
    Rcpp::_["fn"]  = f2,
    Rcpp::_["gr"]  = f3,
    Rcpp::_["y"]   = Rcpp::as<Rcpp::NumericVector>(y),
    Rcpp::_["x"]   = Rcpp::as<Rcpp::NumericMatrix>(x),
    Rcpp::_["mu"]  = mu2,
    Rcpp::_["P"]   = Rcpp::as<Rcpp::NumericMatrix>(P),
    Rcpp::_["alpha"] = alpha,
    Rcpp::_["wt"]    = wt2,
    Rcpp::_["method"] = "BFGS",
    Rcpp::_["hessian"] = true
  );


  

  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:optim] Exiting: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
    
  // Extract posterior mode and Hessian
  Rcpp::NumericVector bstar  = opt_out["par"];
  Rcpp::NumericMatrix A1     = opt_out["hessian"];
  
  
  // -------------------------------
  // Step 4: Standardization (glmb_Standardize_Model)
  // -------------------------------
  
  // bstar is a NumericVector from optim; turn it into a p×1 matrix
  int p_dim = bstar.size();
  Rcpp::NumericMatrix bstar_mat(p_dim, 1);
  for (int i = 0; i < p_dim; ++i) {
    bstar_mat(i, 0) = bstar[i];
  }
  
  // A1 is already a p×p NumericMatrix from optim
  Rcpp::NumericMatrix A1_mat = A1;
  
  // x2 <- x; P2 <- P; mu2 <- 0
  Rcpp::NumericMatrix x2_mat = x;
  Rcpp::NumericMatrix P2_mat = P;
  Rcpp::NumericMatrix mu2_mat(p_dim, 1);
  for (int i = 0; i < p_dim; ++i) {
    mu2_mat(i, 0) = 0.0;
  }

  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:glmb_Standardize_Model] Entering: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  
    
  // Call C++ standardization
  Rcpp::List Standard_Mod = glmb_Standardize_Model(
    Rcpp::as<Rcpp::NumericVector>(y),   // y
    x2_mat,                             // x
    P2_mat,                             // P
    bstar_mat,                          // bstar
    A1_mat                              // A1
  );

  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:glmb_Standardize_Model] Exiting: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
    
  // Extract standardized components
  Rcpp::NumericVector bstar2 = Standard_Mod["bstar2"];
  Rcpp::NumericMatrix A      = Standard_Mod["A"];
  Rcpp::NumericMatrix x2_std = Standard_Mod["x2"];
  Rcpp::NumericMatrix mu2_std= Standard_Mod["mu2"];
  Rcpp::NumericMatrix P2_std = Standard_Mod["P2"];
  Rcpp::NumericMatrix L2Inv  = Standard_Mod["L2Inv"];
  Rcpp::NumericMatrix L3Inv  = Standard_Mod["L3Inv"];
  
  // -------------------------------
  // Step 5: Envelope (EnvelopeOrchestrator)
  // -------------------------------
  
  // RSS_Post2 from your last dispersion loop iteration
  double RSS_ML = NA_REAL;  // matches R: RSS_ML = NA

  
  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:EnvelopeOrchestrator] Entering: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  // Pass wt-Compute wt2 internally  
    
  // Call C++ envelope orchestrator
  Rcpp::List env_out = EnvelopeOrchestrator(
    bstar2,
    A,
    Rcpp::as<Rcpp::NumericVector>(y),
    x2_std,
    mu2_std,
    P2_std,
    alpha,
    wt,
    // wt2,
    n,
    Gridtype,
    
    // n_envopt: treat negative as NULL
    (n_envopt < 0 ? R_NilValue : Rcpp::wrap(n_envopt)),
                
                shape,
                rate,
                RSS_Post2,
                RSS_ML,
                max_disp_perc,
                
                // disp_lower: Nullable<NumericVector> -> Nullable<double>
                (disp_lower.isNull()
                   ? R_NilValue
                   : Rcpp::wrap(Rcpp::as<Rcpp::NumericVector>(disp_lower)[0])),
                     
                     // disp_upper: same logic
                     (disp_upper.isNull()
                        ? R_NilValue
                        : Rcpp::wrap(Rcpp::as<Rcpp::NumericVector>(disp_upper)[0])),
                          
                          use_parallel,
                          use_opencl,
                          verbose
  );
  
  
  if (verbose) {
    Rcpp::Rcout << "[rIndepNormalGammaReg:EnvelopeOrchestrator] Exiting: "
    //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << glmbayes::progress::timestamp_cpp()
                  << "\n";
  }
  
  
    
  // Extract outputs (matching your R code)
  Rcpp::List Env3          = env_out["Env"];
  Rcpp::List gamma_list_new= env_out["gamma_list"];
  Rcpp::List UB_list_new   = env_out["UB_list"];
  double low               = env_out["low"];
  double upp               = env_out["upp"];
  Rcpp::List diagnostics   = env_out["diagnostics"];
  
  
  // -------------------------------
  // Step 6: Simulation (standardized space)
  // -------------------------------
  
  // family / link as CharacterVector, matching R
  Rcpp::CharacterVector family = Rcpp::CharacterVector::create("gaussian");
  Rcpp::CharacterVector link   = Rcpp::CharacterVector::create("identity");
  
  // Choose serial vs parallel simulator
  Rcpp::List sim_temp;

  
  
  
  if (!use_parallel || n == 1) {

    if (verbose) {
      Rcpp::Rcout << "[rIndepNormalGammaReg:rIndepNormalGammaReg_std] Entering: "
      //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << glmbayes::progress::timestamp_cpp()
                    << "\n";
    }
    
    
    if(disp_grid_type==2){
        
        // serial version (assumes same signature)
    sim_temp = rIndepNormalGammaReg_std(
      n,
      Rcpp::as<Rcpp::NumericVector>(y),
      x2_std,
      mu2_std,
      P2_std,
      alpha,
      wt,
      f2,
      Env3,
      gamma_list_new,
      UB_list_new,
      family,
      link,
      progbar,
      verbose
    );
    
    }
    
    if(disp_grid_type==1){
      
    // serial version (assumes same signature)
    sim_temp = rIndepNormalGammaReg_std_face(
      n,
      Rcpp::as<Rcpp::NumericVector>(y),
      x2_std,
      mu2_std,
      P2_std,
      alpha,
      wt,
      f2,
      Env3,
      gamma_list_new,
      UB_list_new,
      family,
      link,
      progbar,
      verbose
    );
    
    }
    
    
    if (verbose) {
      Rcpp::Rcout << "[rIndepNormalGammaReg:rIndepNormalGammaReg_std] Exiting: "
      //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << glmbayes::progress::timestamp_cpp()
                    << "\n";
    }
    
  } else {
    
    if (verbose) {
      Rcpp::Rcout << "[rIndepNormalGammaReg:rIndepNormalGammaReg_std_parallel] Entering: "
      //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << glmbayes::progress::timestamp_cpp()
                    << "\n";
    }
    
    // parallel version (the one you pasted)
    sim_temp = rIndepNormalGammaReg_std_parallel(
      n,
      Rcpp::as<Rcpp::NumericVector>(y),
      x2_std,
      mu2_std,
      P2_std,
      alpha,
      wt,
      f2,
      Env3,
      gamma_list_new,
      UB_list_new,
      family,
      link,
      progbar,
      verbose
    );
    if (verbose) {
      Rcpp::Rcout << "[rIndepNormalGammaReg:rIndepNormalGammaReg_std_parallel] Exiting: "
      //            << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << glmbayes::progress::timestamp_cpp()
                    << "\n";
    }
    
    
  }

    
  // -------------------------------
  // Step 7: Back-transform
  // -------------------------------
  
  Rcpp::NumericMatrix beta_out   = sim_temp["beta_out"];   // n × p
  Rcpp::NumericVector disp_out   = sim_temp["disp_out"];
  Rcpp::NumericVector iters_out  = sim_temp["iters_out"];
  Rcpp::NumericVector weight_out = sim_temp["weight_out"];
  
  int n_draws = beta_out.nrow();

  // Armadillo views
  arma::mat L2Inv_arma(L2Inv.begin(), L2Inv.nrow(), L2Inv.ncol(), false);
  arma::mat L3Inv_arma(L3Inv.begin(), L3Inv.nrow(), L3Inv.ncol(), false);
  arma::mat beta_std(beta_out.begin(), n_draws, p_dim, false); // n × p
  
  // out = L2Inv %*% L3Inv %*% t(beta_out)  (p × n)
  arma::mat out_arma = L2Inv_arma * L3Inv_arma * beta_std.t();
  
  // Add mu to each column: for (i in 1:n) out[, i] <- out[, i] + mu
  arma::vec mu_vec(mu.begin(), mu.size(), false);
  for (int i = 0; i < n_draws; ++i) {
    out_arma.col(i) += mu_vec;
  }
  
  // Convert back to NumericMatrix (p × n)
  Rcpp::NumericMatrix out(p_dim, n_draws);
  std::copy(out_arma.begin(), out_arma.end(), out.begin());
  
  // -------------------------------
  // Final return (mirror R core)
  // -------------------------------
  
  return Rcpp::List::create(
    Rcpp::Named("out")        = out,
    Rcpp::Named("betastar")   = bstar+mu,       // posterior mode from optim() - Add back in prior mean
    Rcpp::Named("disp_out")   = disp_out,
    Rcpp::Named("iters_out")  = iters_out,
    Rcpp::Named("weight_out") = weight_out,
    Rcpp::Named("low")        = low,
    Rcpp::Named("upp")        = upp
  );
}  
  
 
} //sim
} //glmbayes