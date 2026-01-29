// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include "glmbsim.h"
#include "miscfuncs.h"

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
using namespace famfuncs;


void progress_bar3(double x, double N)
{
  // how wide you want the progress meter to be
  int totaldotz=40;
  double fraction = x / N;
  // part of the progressmeter that's already "full"
  int dotz = round(fraction * totaldotz);
  
  Rcpp::Rcout.precision(3);
  Rcout << "\r                                                                 " << std::flush ;
  Rcout << "\r" << std::flush ;
  Rcout << std::fixed << fraction*100 << std::flush ;
  Rcout << "% [" << std::flush ;
  int ii=0;
  for ( ; ii < dotz;ii++) {
    Rcout << "=" << std::flush ;
  }
  // remaining part (spaces)
  for ( ; ii < totaldotz;ii++) {
    Rcout << " " << std::flush ;
  }
  // and back to line begin 
  
  Rcout << "]" << std::flush ;
  
  // and back to line begin 
  
  Rcout << "\r" << std::flush ;
  
}


double p_inv_gamma(double dispersion,double shape,double rate){
  
  return(1- R::pgamma(1/dispersion,shape,1/rate,TRUE,FALSE));
}



double  q_inv_gamma(double p,double shape,double rate,double disp_upper,double disp_lower){
  double p_upp=p_inv_gamma(disp_upper,shape,rate);
  double p_low=p_inv_gamma(disp_lower,shape,rate);
  double p1=p_low+p*(p_upp-p_low);
  double p2=1-p1;
  return(1/ R::qgamma(p2,shape,1/rate,TRUE,FALSE));
}

double r_invgamma(double shape,double rate,double disp_upper,double disp_lower){
  double p= R::runif(0,1);
  return(q_inv_gamma(p,shape,rate,disp_upper,disp_lower));
}



namespace {


//-----------------------------------------------------------------------------
// rindep_norm_gamma_worker: parallel Normal–Gamma simulation with envelope
//-----------------------------------------------------------------------------
struct rindep_norm_gamma_worker : public RcppParallel::Worker {
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
  rindep_norm_gamma_worker(
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
  
  
// --- rindep_norm_gamma_worker implementation ---
void rindep_norm_gamma_worker::operator()(std::size_t begin, std::size_t end) {
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
      double U = safe_runif();
      int J_idx = 0;
      double U_left = U;
      while (true) {
        if (U_left <= PLSD_r[J_idx]) break;
        U_left -= PLSD_r[J_idx];
        ++J_idx;
      }

      // 2) Draw truncated-normal beta row
      for (int j = 0; j < l1; ++j) {
        out_row(0, j) = ctrnorm_cpp(
          logrt_r(J_idx, j),
          loglt_r(J_idx, j),
          -cbars_r(J_idx, j),
          1.0
        );
      }

      // 3) Draw dispersion
      double dispersion = r_invgamma_safe(shape3, rate2, disp_upper, disp_lower);

      // 4) Solve theta (strict row-only)
      for (int j = 0; j < l1; ++j) cbars_small_col(j, 0) = cbars_r(J_idx, j);

      // RcppParallel::RMatrix<double> theta_sol_r =
      //   Inv_f3_with_disp_rmat(Pmat_r, Pmu_r, base_B0_r, base_A_r,
      //                         dispersion, cbars_small_col);

      arma::mat theta_sol = Inv_f3_with_disp_rmat_v2(Pmat_r, Pmu_r, base_B0_r, base_A_r,
                                                     dispersion, cbars_small_col);


      //for (int j = 0; j < l1; ++j) theta_row(0, j) = theta_sol_r(0, j);

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
      double LL_New2_scalar =
        -f2_gaussian_rmat_mat(theta_col, y_r, x_r, mu_r, P_r, alpha_r, wt2_r, 0)[0];

        // Rcout << "Entering f2_gaussian_rmat_mat 2"  << std::endl;

        double LL_Test_scalar =
        -f2_gaussian_rmat_mat(out_col,   y_r, x_r, mu_r, P_r, alpha_r, wt2_r, 0)[0];



        // 7) Upper bounds
        double U2     = safe_runif();
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
        double UB3B = (max_New_LL_UB - max_LL_log_disp + New_LL_log_disp)
          - (lmc1 + lmc2 * dispersion);

        double test1 = (LL_Test_scalar - UB1);
        double test  = test1 - (UB2 + UB3A + UB3B);
        test = test - log_U2;


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


}


// [[Rcpp::export(".rindep_norm_gamma_reg_std_cpp")]]

Rcpp::List  rindep_norm_gamma_reg_std_cpp(int n,NumericVector y,NumericMatrix x,
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
  
  int a1=0;
  double test1=0;
  double test=0;
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(1,l1);
  double a2=0;
  double U=0;
  double U2=0;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
  
  double RSS_Min=UB_list["RSS_Min"];
  NumericVector UB2min=UB_list["UB2min"];
  
//  NumericVector ub2_min=;
  
  
  
  
  // Build cache once outside the loop
  Rcpp::List cache = Inv_f3_precompute_disp(cbars, y, x, mu, P, alpha, wt);
  
  
  for(int i=0;i<n;i++){

    Rcpp::checkUserInterrupt();
    
   if(progbar==1){
     progress_bar3(i, n-1);
     if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
   }
    
    // 3. Test progbar *and* interactive()



    
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
      
      for(int j=0;j<l1;j++){  out(0,j)=ctrnorm_cpp(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);          }
      
      

      // Update this to make distribution contingent on component of the grid
      
      dispersion=r_invgamma(shape3,rate2,disp_upper,disp_lower);
      
      
      
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
      NumericVector cbars_j = cbars(J_out(0), _);
      
      // Call rss_face_at_disp with current dispersion and other parameters
//      double rss_val = rss_face_at_disp(dispersion, cache, cbars_j, y, x, alpha, wt);
      
      // Print or log both RSS values for comparison
//      Rcpp::Rcout << "rss_face_at_disp: " << rss_val << ", inline RSS: " << arma::as_scalar(trans(yxbeta)*yxbeta) << std::endl;
      

      

      // Call the UB2 function with current dispersion and other parameters
  
      
  //    Rcpp::Rcout << "Weights (wt): "; for (int i = 0; i < wt.size(); ++i) { Rcpp::Rcout << wt[i] << " "; } Rcpp::Rcout << std::endl;
      
      
      
      
      
      // Continue with existing UB2 calculation without modification
      UB2 = 0.5 * (1.0 / dispersion) * (arma::as_scalar(trans(yxbeta)*yxbeta) - RSS_Min);
      
      // Print or log both UB2 values for comparison
//      Rcpp::Rcout << "UB2 function: " << ub2_val << ", inline UB2: " << UB2 << std::endl;
      
      // Print or log both UB2 values for comparison
//      Rcpp::Rcout << "UB2_alt function: " << ub2_val_alt << ", inline UB2: " << UB2 << std::endl;
      
//      Rcpp::Rcout << "Index: " << J_out(0) << ", UB2min value: " << UB2min[J_out(0)] << std::endl;
      
//      UB2=0.5*(1/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_ML);
      UB2=0.5*(1/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_Min);
      
      // Subtract UB2min --> Should improve acceptance
      
      UB2=UB2-UB2min[J_out(0)];
      
      
      // Block 3: UB3A (adjusts because probabilities of components in grid are different from original grid)
      // Investigate whether changing probabilities of grid components for proposal
      // allows us to do away with this term and to thereby improve the acceptance rate
      
      // This is likely time consuming part
      

      
      for(int j=J_out(0);j<(J_out(0)+1);j++){
        thetabars_temp=thetabars_new(0,_); // Changed
        
        
        cbars_temp=cbars(j,_);
        arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
        arma::vec  cbars_temp2(cbars_temp.begin(), l1);
        
        New_LL(j)=arma::as_scalar(-0.5*trans(thetabars_temp2)*P2*thetabars_temp2
                                    +trans(cbars_temp2)*thetabars_temp2);
        
      }
      

      // Modified UB3A 
      
      UB3A= lg_prob_factor(J_out(0))+lmc1+lmc2*dispersion-New_LL(J_out(0));
      
      // Block 4: UB3B  
      
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
  
  // Temporarily just return non-sense constants equal to all 1
  
  return Rcpp::List::create(Rcpp::Named("beta_out")=beta_out,Rcpp::Named("disp_out")=disp_out,
                            Rcpp::Named("iters_out")=iters_out,Rcpp::Named("weight_out")=weight_out);  
  
  
  
}



// [[Rcpp::export(".rindep_norm_gamma_reg_std_parallel_cpp")]]

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



  Rcpp::List cache = Inv_f3_precompute_disp(cbars, y, x, mu, P, alpha, wt);


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
  rindep_norm_gamma_worker worker(
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

  
  // --- Single-draw test (serial) ---
  int m_test = 1;
  auto t0 = std::chrono::steady_clock::now();
  worker(0, m_test);  // run worker serially for 1 observation
  auto t1 = std::chrono::steady_clock::now();
  double elapsed_test_sec = std::chrono::duration<double>(t1 - t0).count();
  
  if (verbose) Rcpp::Rcout << "[Pilot] Single test run took " << elapsed_test_sec << "s.\n";
  
  // --- Conservative calibration sizing (time-bounded) ---
  // Use single test to bound worst-case per-observation time in ms
  double per_obs_ms_serial = elapsed_test_sec * 1000.0 / std::max(1, m_test);
  
  // Aim for ~1% of n but cap at ~5 minutes worst-case based on serial bound
  int m1 = std::max(1, (int)std::ceil(0.01 * (double)n));
  int m2 = std::max(1, (int)std::floor(300000.0 / std::max(1.0, per_obs_ms_serial))); // 300k ms ≈ 5 min
  int m_stage = std::min(m1, m2);   // <-- defined here, before use
  
  if (verbose) {Rcpp::Rcout << "Calibrating simulation time estimate using " << m_stage
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
  if (verbose){  Rcpp::Rcout << "[CALIB] Calibration elapsed = " << cal_elapsed_sec
              << " s for " << m_stage << " observations.\n";
  Rcpp::Rcout << "[CALIB] per_obs_sec = " << per_obs_sec
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
  
  if (verbose) {Rcpp::Rcout << "[Estimate] Simulation time for " << n << " observations: "
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
    Rcpp::Rcout << "[Simulation] >>> Starting full run at "
                << Rcpp::as<std::string>(now[0]) << " <<<\n";
   }
  
  // --- Capture start time ---
  double sim_start = Rcpp::as<double>(
    Rcpp::Function("as.numeric")(Rcpp::Function("Sys.time")())
  );
  
  // Parallel loop
  RcppParallel::parallelFor(0, n, worker);
  // worker(0,n);
  
  if (verbose) Rcpp::Rcout << "Exiting Parallel Worker" << std::endl;
  
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
    Rcpp::Rcout << "[Simulation] >>> Exiting full run at "
                << Rcpp::as<std::string>(now[0]) << " <<<\n";
    Rcpp::Rcout << "[Simulation] Simulation completed in: "
                << h_elapsed << " h  " << m_elapsed << " m  " << s_elapsed << " s.\n";
   }  



  return Rcpp::List::create(
    Rcpp::Named("beta_out")   = beta_out,
    Rcpp::Named("disp_out")   = disp_out,
    Rcpp::Named("iters_out")  = iters_out,
    Rcpp::Named("weight_out") = weight_out
  );
}




// [[Rcpp::export]]

Rcpp::List rindep_norm_gamma_reg_cpp(
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
  
  for (int j = 0; j < 10; ++j) {
    
    // --- Call rnorm_reg_cpp (C++ version of .rnorm_reg_cpp) ---
    cpp_out = rnorm_reg_cpp(
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
    double RSS2_post = arma::dot(resid_ml, resid_ml);
    
    // Update shape, rate, dispersion
    double shape2 = shape + static_cast<double>(n_obs) / 2.0;
    double rate2  = rate  + RSS_Post2 / 2.0;
    
    dispersion2 = rate2 / (shape2 - 1.0);
  }
  
  
  // -------------------------------
  // Posterior Mode + Hessian Block
  // -------------------------------
  
  // Extract posterior mode from sampler
  arma::vec betastar = Rcpp::as<arma::vec>(cpp_out["coef.mode"]);
  double dispstar = dispersion2;
  
  // wt2 = wt / dispstar
  Rcpp::NumericVector wt2(n_obs);
  for (int i = 0; i < n_obs; ++i)
    wt2[i] = wt[i] / dispstar;
  
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
    Rcpp::Rcout << "[PosteriorMode] >>> Entering optim() call <<<\n";
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
  
  // Call C++ standardization
  Rcpp::List Standard_Mod = glmb_Standardize_Model(
    Rcpp::as<Rcpp::NumericVector>(y),   // y
    x2_mat,                             // x
    P2_mat,                             // P
    bstar_mat,                          // bstar
    A1_mat                              // A1
  );
  
  // Extract standardized components
  Rcpp::NumericVector bstar2 = Standard_Mod["bstar2"];
  Rcpp::NumericMatrix A      = Standard_Mod["A"];
  Rcpp::NumericMatrix x2_std = Standard_Mod["x2"];
  Rcpp::NumericMatrix mu2_std= Standard_Mod["mu2"];
  Rcpp::NumericMatrix P2_std = Standard_Mod["P2"];
  Rcpp::NumericMatrix L2Inv  = Standard_Mod["L2Inv"];
  Rcpp::NumericMatrix L3Inv  = Standard_Mod["L3Inv"];
  
  // -------------------------------
  // Step 5: Envelope (EnvelopeOrchestrator_cpp)
  // -------------------------------
  
  // RSS_Post2 from your last dispersion loop iteration
  double RSS_ML = NA_REAL;  // matches R: RSS_ML = NA
  

  
  
  
  // Call C++ envelope orchestrator
  Rcpp::List env_out = EnvelopeOrchestrator_cpp(
    bstar2,
    A,
    Rcpp::as<Rcpp::NumericVector>(y),
    x2_std,
    mu2_std,
    P2_std,
    alpha,
    wt,
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
    // serial version (assumes same signature)
    sim_temp = rindep_norm_gamma_reg_std_cpp(
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
  } else {
    // parallel version (the one you pasted)
    sim_temp = rindep_norm_gamma_reg_std_parallel_cpp(
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
    Rcpp::Named("betastar")   = bstar,       // posterior mode from optim()
    Rcpp::Named("disp_out")   = disp_out,
    Rcpp::Named("iters_out")  = iters_out,
    Rcpp::Named("weight_out") = weight_out,
    Rcpp::Named("low")        = low,
    Rcpp::Named("upp")        = upp
  );
}  
  
 
