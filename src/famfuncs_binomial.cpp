// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <cmath>
#include "RcppArmadillo.h"
#include <limits>
#include <RcppParallel.h>
#define MATHLIB_STANDALONE
#include "nmath_local.h"
#include "dpq_local.h"
#include "openclPort.h"
#include "famfuncs.h"
#include "progress_utils.h"

// Namespaces

using namespace Rcpp;
using namespace RcppParallel;
using namespace openclPort;
using namespace glmbayes::fam;
using namespace glmbayes::progress;



/////////////////////////////////////////////////////////////


inline double dbinom_raw_local(double x, double n, double p, double q, int give_log) {
  // Pass-through version for now
  return dbinom_raw(x, n, p, q, give_log);
}


///////////////////////////////////////////////////////////

namespace glmbayes{

namespace fam {
// Neg log binomial likelihood, vectorized
void neg_dbinom_glmb_rmat(const RVector<double>& x,          // success proportion
                          const RVector<double>& N,          // trial count
                          const std::vector<double>& p_vec,  // success probabilities
                          const std::vector<double>& q_vec,  // failure probabilities
                          std::vector<double>& res,          // output buffer
                          const int lg)                      // log=TRUE?
{
  std::size_t n = x.size();
  if (res.size() != n)
    res.resize(n);
  
  for (std::size_t i = 0; i < n; ++i) {
    double trials  = std::round(N[i]);
    double success = std::round(x[i] * N[i]);
    double p       = p_vec[i];
    double q       = q_vec[i];
    
    // Thread-safe backend log-density
    res[i] = -dbinom_raw_local(success, trials, p, q, lg);
    
    // Optional diagnostics
    // if (!std::isfinite(res[i])) {
    //     Rcpp::Rcout << "[WARN] neg_dbinom_glmb_rmat: non-finite res[" << i << "]\n"
    //                 << "  trials  = " << trials  << "\n"
    //                 << "  success = " << success << "\n"
    //                 << "  p       = " << p       << "\n"
    //                 << "  q       = " << q       << "\n"
    //                 << "  raw x   = " << x[i]    << "\n"
    //                 << "  raw N   = " << N[i]    << "\n";
    // }
  }
}




NumericVector dbinom_glmb(NumericVector x, NumericVector N, NumericVector means, int lg) {
  int n = x.size();
  NumericVector res(n);
  
  for (int i = 0; i < n; i++) {
    // Round to nearest integer for trial count and success count
    int trials  = static_cast<int>(std::round(N[i]));
    int success = static_cast<int>(std::round(x[i] * N[i]));
    
    // Clamp probabilities to avoid extreme values
    double p = std::min(1.0, std::max(0.0, means[i]));
    
    // Evaluate binomial log-likelihood
    res[i] = R::dbinom(success, trials, p, lg);

    
      }
  
  return res;
}


///////////////////////// Logit Functions ///////////////////////////////////////

NumericVector  f1_binomial_logit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);


    for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
  
    
    xb2=exp(-alpha2- x2 * b2);

    for(int j=0;j<l1;j++){
    xb(j)=1/(1+xb(j));
    }

    yy=-dbinom_glmb(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);

    }
    
    return res;      
}




NumericVector  f2_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
     //   Note: Does not seem to be used-- Editing out         
  //  NumericVector invwt=1/sqrt(wt);

    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);

    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 


    double res1=0;


    for(int i=0;i<m1;i++){
      
      
      Rcpp::checkUserInterrupt();
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 

    bmu2=b2-mu2;
        
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
    xb2=exp(-alpha2- x2 * b2);

    for(int j=0;j<l1;j++){
    xb(j)=1/(1+xb(j));
    }


    
    
    yy=-dbinom_glmb(y,wt,xb,true);
    //yy=-dbinom_glmb_raw(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), res1);

    }
    
    return res;      
}



arma::vec f2_binomial_logit_rmat(
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    const int progbar /*=0*/
) {
  std::size_t l1 = x.nrow();
  std::size_t l2 = x.ncol();
  std::size_t m1 = b.ncol();
  
  // Armadillo views over RMatrix memory
  arma::mat b2full(const_cast<double*>(&*b.begin()), l2, m1, false);
  arma::mat x2(const_cast<double*>(&*x.begin()), l1, l2, false);
  arma::mat mu2(const_cast<double*>(&*mu.begin()), l2, 1, false);
  arma::mat P2(const_cast<double*>(&*P.begin()), l2, l2, false);
  arma::mat alpha2(const_cast<double*>(&*alpha.begin()), l1, 1, false);
  
  arma::vec res(m1, arma::fill::none);
  arma::mat bmu(l2, 1, arma::fill::none);
  
  // Buffers (declare once, reuse per i)
  std::vector<double> p_vec(l1), q_vec(l1), yy_temp(l1);
  
  for (std::size_t i = 0; i < m1; ++i) {
    arma::mat b_i(b2full.colptr(i), l2, 1, false);
    
    // Mahalanobis penalty
    bmu = b_i - mu2;
    double mahal = 0.5 * arma::as_scalar(bmu.t() * P2 * bmu);
    
    // Compute p and q stably for each observation
    for (std::size_t j = 0; j < l1; ++j) {
      double eta = alpha2(j,0) + arma::dot(x2.row(j), b_i);
      double p, q;
      if (eta >= 0) {
        double e = std::exp(-eta);
        p = 1.0 / (1.0 + e);
        q = e / (1.0 + e);
      } else {
        double e = std::exp(eta);
        p = e / (1.0 + e);
        q = 1.0 / (1.0 + e);
      }
      p_vec[j] = p;
      q_vec[j] = q;
    }
    
    // Call refactored backend with both p and q
    neg_dbinom_glmb_rmat(y, wt, p_vec, q_vec, yy_temp, /*lg=*/1);
    
    // Accumulate with penalty
    res(i) = std::accumulate(yy_temp.begin(), yy_temp.end(), mahal);
  }
  
  return res;
}




arma::mat  f3_binomial_logit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat y2(y.begin(), l1, 1, false);
    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
       
   
       
    NumericMatrix Ptemp(l1,l1);  
      
    for(int i=0;i<l1;i++){
     Ptemp(i,i)=wt(i); 
    }  
    
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);
    NumericMatrix out(l2,m1);


    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 
    arma::mat Ptemp2(Ptemp.begin(), l1, l1, false);
    arma::mat out2(out.begin(), l2, m1, false);
    
    NumericMatrix::Column outtemp=out(_,0);
    //NumericMatrix::Row outtempb=out(0,_);

    arma::mat outtemp2(outtemp.begin(),1,l2,false);
    //arma::mat outtempb2(outtempb.begin(),1,l2,false);
    
    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
    if(progbar==1){ 
    progress_bar2(i, m1-1);
    if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
    };  
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);

    bmu2=b2-mu2;

//  		p<-1/(1+t(exp(-alpha-x%*%b)))
//		t(x)%*%((t(p)-y)*wt)+P%*%(b-mu)


    xb2=exp(-alpha2- x2 * b2);
    
    

    for(int j=0;j<l1;j++){
    xb(j)=1/(1+xb(j));  
    xb(j)=(xb(j)-y(j))*wt(j);
    }



    outtemp2= P2 * bmu2+x2.t() * xb2;
    }
    
   // return  b;
  
    return trans(out2);      
}




///////////////////////// Probit Functions ///////////////////////////////////////

NumericVector  f1_binomial_probit(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);


    for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
 
    xb2=alpha2+ x2 * b2;   
    xb=pnorm(xb,0.0,1.0);
	

    yy=-dbinom_glmb(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);

    }
    
    return res;      
}




NumericVector  f2_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
     //   Note: Does not seem to be used-- Editing out        
//    NumericVector invwt=1/sqrt(wt);

    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);

    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 

    double res1=0;


    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 

    bmu2=b2-mu2;
        
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
  
    xb2=alpha2+ x2 * b2;   
    xb=pnorm(xb,0.0,1.0);

    yy=-dbinom_glmb(y,wt,xb,true);


    res(i) =std::accumulate(yy.begin(), yy.end(), res1);

    }
    
    return res;      
}


arma::vec f2_binomial_probit_rmat(
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    const int progbar
) {
  std::size_t l1 = x.nrow();
  std::size_t l2 = x.ncol();
  std::size_t m1 = b.ncol();
  
  // Armadillo views over RMatrix memory
  arma::mat b2full(const_cast<double*>(&*b.begin()), l2, m1, false);
  arma::mat x2(const_cast<double*>(&*x.begin()), l1, l2, false);
  arma::mat mu2(const_cast<double*>(&*mu.begin()), l2, 1, false);
  arma::mat P2(const_cast<double*>(&*P.begin()), l2, l2, false);
  arma::mat alpha2(const_cast<double*>(&*alpha.begin()), l1, 1, false);
  
  arma::vec res(m1, arma::fill::none);
  arma::mat bmu(l2, 1, arma::fill::none);
  
  // Buffers
  std::vector<double> p_vec(l1), q_vec(l1), yy_temp(l1);
  
  for (std::size_t i = 0; i < m1; ++i) {
    arma::mat b_i(b2full.colptr(i), l2, 1, false);
    
    bmu = b_i - mu2;
    double mahal = 0.5 * arma::as_scalar(bmu.t() * P2 * bmu);
    
    // Compute linear predictor
    arma::colvec eta = alpha2 + x2 * b_i;
    
    // Compute p and q stably via probit link
    for (std::size_t j = 0; j < l1; j++) {
      double z = eta(j);
      double p = pnorm5_local(z, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/0);
      double q = pnorm5_local(z, 0.0, 1.0, /*lower_tail=*/0, /*log_p=*/0);
      p_vec[j] = p;
      q_vec[j] = q;
    }
    
    // Call refactored backend
    neg_dbinom_glmb_rmat(y, wt, p_vec, q_vec, yy_temp, 1);
    
    res(i) = std::accumulate(yy_temp.begin(), yy_temp.end(), mahal);
  }
  
  return res;
}



arma::mat  f3_binomial_probit(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat y2(y.begin(), l1, 1, false);
    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
       
   
       
    NumericMatrix Ptemp(l1,l1);  
      
    for(int i=0;i<l1;i++){
     Ptemp(i,i)=wt(i); 
    }  
    
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);
    NumericMatrix out(l2,m1);
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector d1(l1);


    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 
    arma::mat Ptemp2(Ptemp.begin(), l1, l1, false);
    arma::mat out2(out.begin(), l2, m1, false);
    
    NumericMatrix::Column outtemp=out(_,0);
    //NumericMatrix::Row outtempb=out(0,_);

    arma::mat outtemp2(outtemp.begin(),1,l2,false);
    //arma::mat outtempb2(outtempb.begin(),1,l2,false);
    
    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);

    bmu2=b2-mu2;
    
    xb2=alpha2+ x2 * b2;
    p1=pnorm(xb,0.0,1.0);
    p2=pnorm(-xb,0.0,1.0);
    d1=dnorm(xb,0.0,1.0);
    
//    -t(x)%*%as.matrix(((y*dnorm(alpha+x%*%b)/p1)-(1-y)*dnorm(alpha+x%*%b)/p2)*wt)+P%*%(b-mu)


    for(int j=0;j<l1;j++){
    xb(j)=(y(j)*d1(j)/p1(j)-(1-y(j))*d1(j)/p2(j))*wt(j);    
    }


    outtemp2= P2 * bmu2-x2.t() * xb2;
    }
    
   // return  b;
  
    return trans(out2);      
}



///////////////////////// cLOGLOG FUNCTION ///////////////////////////////////////

NumericVector  f1_binomial_cloglog(NumericMatrix b,NumericVector y,NumericMatrix x,NumericVector alpha,NumericVector wt)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);


    for(int i=0;i<m1;i++){
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
 
    xb2=alpha2+ x2 * b2;   
    xb=exp(-exp(xb));
  
    for(int j=0;j<l1;j++){
    xb(j)=1-xb(j);
    }
  

    yy=-dbinom_glmb(y,wt,xb,true);
    

    res(i) =std::accumulate(yy.begin(), yy.end(), 0.0);

    }
    
    return res;      
}




NumericVector  f2_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
     
//   Note: Does not seem to be used-- Editing out      
//    NumericVector invwt=1/sqrt(wt);

    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);

    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 

    double res1=0;


    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 

    bmu2=b2-mu2;
        
    res1=0.5*arma::as_scalar(bmu2.t() * P2 *  bmu2);
    
    xb2=alpha2+ x2 * b2;   
    
    for(int j=0;j<l1;j++){
    xb(j)=1-exp(-exp(xb(j)));
    }

    yy=-dbinom_glmb(y,wt,xb,true);


    res(i) =std::accumulate(yy.begin(), yy.end(), res1);

    }
    
    return res;      
}

arma::vec f2_binomial_cloglog_rmat(
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    const int progbar /*=0*/
) {
  std::size_t l1 = x.nrow();
  std::size_t l2 = x.ncol();
  std::size_t m1 = b.ncol();
  
  // Armadillo views over RMatrix memory
  arma::mat b2full(const_cast<double*>(&*b.begin()), l2, m1, false);
  arma::mat x2(const_cast<double*>(&*x.begin()), l1, l2, false);
  arma::mat mu2(const_cast<double*>(&*mu.begin()), l2, 1, false);
  arma::mat P2(const_cast<double*>(&*P.begin()), l2, l2, false);
  arma::mat alpha2(const_cast<double*>(&*alpha.begin()), l1, 1, false);
  
  arma::vec res(m1, arma::fill::none);
  arma::mat bmu(l2, 1, arma::fill::none);
  
  // Buffers reused per coefficient index i
  std::vector<double> p_vec(l1), q_vec(l1), yy_temp(l1);
  
  for (std::size_t i = 0; i < m1; ++i) {
    arma::mat b_i(b2full.colptr(i), l2, 1, false);
    
    // Mahalanobis penalty
    bmu = b_i - mu2;
    double mahal = 0.5 * arma::as_scalar(bmu.t() * P2 * bmu);
    
    // Linear predictor: cloglog uses η = α + x b
    arma::colvec eta = alpha2 + x2 * b_i;
    
    // p = 1 - exp(-exp(η)), q = exp(-exp(η))
    for (std::size_t j = 0; j < l1; ++j) {
      double t = std::exp(eta(j));              // t = exp(η), may be Inf for large η
      double q = std::exp(-t);                  // q = exp(-t), safe even if t = Inf -> q = 0
      double p;
      // For tiny t, use expm1 to avoid cancellation: 1 - exp(-t) ≈ t
      if (t < 1e-6) {
        p = -std::expm1(-t);                  // p ≈ t for small t
      } else {
        p = 1.0 - q;
      }
      p_vec[j] = p;
      q_vec[j] = q;
    }
    
    // Tail-stable backend with both p and q
    neg_dbinom_glmb_rmat(y, wt, p_vec, q_vec, yy_temp, /*lg=*/1);
    
    // Accumulate with penalty
    res(i) = std::accumulate(yy_temp.begin(), yy_temp.end(), mahal);
  }
  
  return res;
}



arma::mat  f3_binomial_cloglog(NumericMatrix b,NumericVector y, NumericMatrix x,NumericMatrix mu,NumericMatrix P,NumericVector alpha,NumericVector wt, int progbar=0)
{
 
    // Get dimensions of x - Note: should match dimensions of
    //  y, b, alpha, and wt (may add error checking)
    
    // May want to add method for dealing with alpha and wt when 
    // constants instead of vectors
    
    int l1 = x.nrow(), l2 = x.ncol();
    int m1 = b.ncol();
    
//    int lalpha=alpha.nrow();
//    int lwt=wt.nrow();

    Rcpp::NumericMatrix b2temp(l2,1);

    arma::mat y2(y.begin(), l1, 1, false);
    arma::mat x2(x.begin(), l1, l2, false); 
    arma::mat alpha2(alpha.begin(), l1, 1, false); 

    Rcpp::NumericVector xb(l1);
    arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
       
   
       
    NumericMatrix Ptemp(l1,l1);  
      
    for(int i=0;i<l1;i++){
     Ptemp(i,i)=wt(i); 
    }  
    
    // Moving Loop inside the function is key for speed

    NumericVector yy(l1);
    NumericVector res(m1);
    NumericMatrix bmu(l2,1);
    NumericMatrix out(l2,m1);
    NumericVector p1(l1);
    NumericVector p2(l1);
    NumericVector atemp(l1);


    arma::mat mu2(mu.begin(), l2, 1, false); 
    arma::mat bmu2(bmu.begin(), l2, 1, false); 
    arma::mat P2(P.begin(), l2, l2, false); 
    arma::mat Ptemp2(Ptemp.begin(), l1, l1, false);
    arma::mat out2(out.begin(), l2, m1, false);
    
    NumericMatrix::Column outtemp=out(_,0);
    //NumericMatrix::Row outtempb=out(0,_);

    arma::mat outtemp2(outtemp.begin(),1,l2,false);
    //arma::mat outtempb2(outtempb.begin(),1,l2,false);
    
    for(int i=0;i<m1;i++){
      Rcpp::checkUserInterrupt();
      
      if(progbar==1){ 
        progress_bar2(i, m1-1);
        if(i==m1-1) {Rcpp::Rcout << "" << std::endl;}
      };  
      
      
    b2temp=b(Range(0,l2-1),Range(i,i));
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    
    NumericMatrix::Column outtemp=out(_,i);
    arma::mat outtemp2(outtemp.begin(),1,l2,false);

    bmu2=b2-mu2;
    
    xb2=alpha2+ x2 * b2;

    for(int j=0;j<l1;j++){
    p1(j)=1-exp(-exp(xb(j)));
    p2(j)=exp(-exp(xb(j)));
    atemp(j)=exp(xb(j)-exp(xb(j)));
    xb(j)=((y(j)*atemp(j)/p1(j))-((1-y(j))*atemp(j)/p2(j)))*wt(j);    
    }


    outtemp2= P2 * bmu2-x2.t() * xb2;
    }
    
  
    return trans(out2);      
}

}  //famfuncs
}  //glmbayes
