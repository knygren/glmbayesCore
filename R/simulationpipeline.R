#' Low-Level Simulation Pipeline for Bayesian GLMs
#'
#' @description
#' A detailed overview of the low-level simulation pipeline used by
#' \code{rglmb()} and related functions. These routines implement the
#' optimization → standardization → envelope sizing → envelope construction →
#' sampling → back-transformation workflow described in Nygren & Nygren (2006).
#'
#' @details
#' (summaries of each step)
#'
#' @example inst/examples/Ex_rnnorm_reg_std.R
#'
#' @name SimulationPipeline
NULL




#' Return family functions used during simulation and post processing
#'
#' This function takes as input a \code{\link{family}} object and returns a 
#' set of functions that are used during simulation and summarization of models 
#' using the \code{\link{glmb}}, and \code{\link{rglmb}} functions.
#' 
#' @name glmbfamfunc
#' @aliases
#' glmbfamfunc
#' print.glmbfamfunc
#' @param family an object of class \code{\link{family}}
#' @param x an object of class \code{"glmbfamfunc"} for which a printed output is desired.
#' @param \ldots additional optional arguments
#' @return A list with the following components
#' \item{f1}{Negative log-likelihood function}
#' \item{f2}{Negative log-posterior function}
#' \item{f3}{Gradient function for negative log-posterior function}
#' \item{f4}{Deviance function}
#' \item{f7}{Another function}
#' @details  This function takes as input a family and returns a set of functions related to the family.
#' @example inst/examples/Ex_glmbfamfunc.R
#' @export
#' @rdname glmbfamfunc
#' @order 1

glmbfamfunc<-function(family){
  
  # need to add handling for offsets 
  
  # f1-(negative)Log-Likelihood function
  # f2-(negative)Log-Posterior density
  # f3-(negative) Gradient for Log-Posterior density 
  # f4-Deviance function (note multiplies the difference in two times log-likelihood by dispersion)
  
  if(family$family=="gaussian")
  {
    f1<-function(b,y,x,alpha=0,wt=1){
      Xb<-alpha+x%*%b
      -sum(dnorm(y, mean=Xb,sd=sqrt(1/wt),log=TRUE))
    }
    f2<-function(b,y,x,mu,P,alpha=0,wt=1){
      Xb<-alpha+x%*%b
      -sum(dnorm(y, mean=Xb,sd=sqrt(1/wt),log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
    }
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      l2<-length(y)
      ltemp<-length(wt)
      yxb2<-NULL
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      xb<-alpha+x%*%b-y
      yXb2<-P%*%(b-mu)+t(x)%*%(Ptemp%*%xb)
      yXb2
    }
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      
      ## adds 2*log-likelihood at mean=y which measn this is 2* Negative log-likelihood difference from saturated
      ## Multiplication by dispersion likely is/was incorrect
      ## As long as single value for the parameters and dispersion, can do the below
      
      #dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dnorm(y, mean=y,sd=sqrt(1/wt),log=TRUE)))
      2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dnorm(y, mean=y,sd=sqrt(dispersion/wt),log=TRUE))
      
    }
    
    # Edit This
    #	f5<-f2_gaussian
    #	f6<-f3_gaussian
    
    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
      l2<-length(y)
      ltemp<-length(wt)
      yxb2<-NULL
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      Pout<-t(x)%*%(Ptemp)%*%x
      Pout
    }
    
    
    
    
    
    
  }
  
  # Check if Poisson weight should be outside or inside function
  
  
  if(family$family=="poisson"||family$family=="quasipoisson")
  {
    f1<-function(b,y,x,alpha=0,wt=1){
      lambda<-t(exp(alpha+x%*%b))
      #-sum(dpois(y, lambda,log=TRUE)*wt)
      -sum(dpois2(y, lambda,log=TRUE)*wt)
      
      
      
    }
    f2<-function(b,y,x,mu,P,alpha=0,wt=1){
      lambda<-t(exp(alpha+x%*%b))
      #-sum(dpois(y, lambda,log=TRUE)*wt)+0.5*t((b-mu))%*%P%*%(b-mu)
      -sum(dpois2(y, lambda,log=TRUE)*wt)+0.5*t((b-mu))%*%P%*%(b-mu)
    }
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      -t(x)%*%((y-exp(alpha+x%*%b))*wt)+P%*%(b-mu)
      
    }
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      #dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dpois(y, y,log=TRUE)*wt))
      (2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dpois2(y, y,log=TRUE)*(wt/dispersion)))
    }
    #  f5<-f2_poisson
    #  f6<-f3_poisson
    
    #    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
    #      l2<-length(y)
    #      ltemp<-length(wt)
    #      yxb2<-NULL
    #      if(ltemp==1){
    #        Ptemp<-wt*diag(l2)
    #      }
    #      else {
    #        Ptemp<-diag(wt)      
    #      }
    #      Pout<-t(x)%*%(Ptemp)%*%x*exp(alpha+x%*%b)
    #      Pout
    #    }
    
    f7<-function(b, y, x, mu, P, alpha = 0, wt = 1) {
      l2 <- length(y)
      l1 <- length(b)
      ltemp <- length(wt)
      
      # Construct diagonal weight matrix
      if (ltemp == 1) {
        Ptemp <- wt * diag(l2)
      } else {
        Ptemp <- diag(wt)
      }
      
      # Compute expected counts
      mu_i <- exp(alpha + x %*% b)
      
      # Initialize output matrix
      Pout <- matrix(0, nrow = l1, ncol = l1)
      
      # Loop over observations
      for (i in 1:l2) {
        xi <- x[i, , drop = FALSE]  # row as matrix
        weight_i <- Ptemp[i, i] * mu_i[i]
        Pout <- Pout + weight_i * (t(xi) %*% xi)
      }
      
      Pout
    }
    
    
    
    
  }
  
  if(family$family %in%  c("binomial","quasibinomial") && family$link=="logit")
  {
    f1<-function(b,y,x,alpha=0,wt=1){
      lambda<-t(exp(alpha+x%*%b))
      p<-lambda/(1+lambda)
      -sum(dbinom(round(wt*y),round(wt), p,log=TRUE))
    }
    f2<-function(b,y,x,mu,P,alpha=0,wt=1){
      lambda<-t(exp(alpha+x%*%b))
      p<-lambda/(1+lambda)
      -sum(dbinom(round(wt*y),round(wt),p,log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
    }
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      p<-1/(1+t(exp(-alpha-x%*%b)))
      t(x)%*%((t(p)-y)*wt)+P%*%(b-mu)
    }
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      #dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dbinom(round(wt*y),round(wt),y,log=TRUE)))
      (2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dbinom(round((wt/dispersion)*y),round(wt/dispersion),y,log=TRUE)))
    }
    #  f5<-f2_binomial_logit
    #  f6<-f3_binomial_logit
    
    
    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
      
      l2<-length(y)
      l1<-length(b)
      ltemp<-length(wt)
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      
      p<-1/(1+t(exp(-alpha-x%*%b)))
      
      
      Pout<-matrix(0,nrow=l1,ncol=l1)
      
      for(i in 1:l2){
        Pout<-Pout+Ptemp[i,i]*p[i]*(1-p[i])*x[i,]%*%t(x[i,])
      }
      
      Pout
    }
    
    
  }
  
  if(family$family %in%  c("binomial","quasibinomial") && family$link=="probit")
  {
    f1<-function(b,y,x,alpha=0,wt=1){
      p<-pnorm(alpha+x%*%b)
      -sum(dbinom(round(wt*y),round(wt), p,log=TRUE))
    }
    f2<-function(b,y,x,mu,P,alpha=0,wt=1){
      p<-pnorm(alpha+x%*%b)
      -sum(dbinom(round(wt*y),round(wt),p,log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
    }
    
    
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      p1<-pnorm(alpha+x%*%b)
      p2<-pnorm(-alpha-x%*%b)
      -t(x)%*%as.matrix(((y*dnorm(alpha+x%*%b)/p1)-(1-y)*dnorm(alpha+x%*%b)/p2)*wt)+P%*%(b-mu)
    }
    
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      #dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dbinom(round(wt*y),round(wt),y,log=TRUE)))
      (2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dbinom(round((wt/dispersion)*y),round(wt/dispersion),y,log=TRUE)))
    }
    #  f5<-f2_binomial_probit
    #  f6<-f3_binomial_probit
    
    
    ##########  This should be replaced!
    
    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
      
      l2<-length(y)
      l1<-length(b)
      ltemp<-length(wt)
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      
      p<-1/(1+t(exp(-alpha-x%*%b)))
      
      
      Pout<-matrix(0,nrow=l1,ncol=l1)
      
      for(i in 1:l2){
        Pout<-Pout+Ptemp[i,i]*p[i]*(1-p[i])*x[i,]%*%t(x[i,])
      }
      
      Pout
    }
    
    
    
    
    
  }
  if(family$family %in%  c("binomial","quasibinomial") && family$link=="cloglog")
  {
    f1<-function(b,y,x,alpha=0,wt=1){
      Xb<-alpha+x%*%b
      p<-1-exp(-exp(Xb))
      -sum(dbinom(round(wt*y),round(wt), p,log=TRUE))
    }
    f2<-function(b,y,x,mu,P,alpha=0,wt=1){
      Xb<-alpha+x%*%b
      p<-1-exp(-exp(Xb))
      -sum(dbinom(round(wt*y),round(wt),p,log=TRUE))+0.5*t((b-mu))%*%P%*%(b-mu)
    }
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      Xb<-alpha+x%*%b
      p1<-1-exp(-exp(Xb))
      p2<-exp(-exp(Xb))
      atemp<-exp(Xb-exp(Xb))
      -t(x)%*%as.matrix(((y*atemp/p1)-(1-y)*atemp/p2)*wt)+P%*%(b-mu)
    }
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      dispersion*(2*f1(b,y,x,alpha,wt)+2*sum(dbinom(round(wt*y),round(wt),y,log=TRUE)))
      (2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dbinom(round((wt/dispersion)*y),round(wt/dispersion),y,log=TRUE)))
    }
    #  f5<-f2_binomial_cloglog
    #  f6<-f3_binomial_cloglog
    
    ##########  This should be replaced!
    
    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
      
      l2<-length(y)
      l1<-length(b)
      ltemp<-length(wt)
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      
      p<-1/(1+t(exp(-alpha-x%*%b)))
      
      
      Pout<-matrix(0,nrow=l1,ncol=l1)
      
      for(i in 1:l2){
        Pout<-Pout+Ptemp[i,i]*p[i]*(1-p[i])*x[i,]%*%t(x[i,])
      }
      
      Pout
    }
    
    
    
  }
  if(family$family=="Gamma" && family$link=="log")
  {
    
    f1<-function(b,y,x,alpha=0,wt=1){
      mu<-t(exp(alpha+x%*%b))
      disp2<-1/wt
      -sum(dgamma(y,shape=1/disp2,scale=mu*disp2,log=TRUE))
    }
    
    #  f1<-f1_gamma
    f2 <- function(b, y, x, mu, P, alpha = 0, wt = 1) {
      
      eta    <- alpha + x %*% b
      scale2 <- t(exp(eta - log(wt)))
      disp2  <- 1 / wt
      
      ## TEMPORARY DIAGNOSTIC (optional)
      #      if (any(!is.finite(scale2)) || any(scale2 == 0)) {
      #        cat("\n*** scale2 is zero or non-finite in f2 ***\n")
      #        cat("alpha:\n")
      #        print(alpha)
      #        cat("range(eta):", range(eta), "\n")
      #        cat("range(scale2):", range(scale2), "\n")
      #      }
      
      -sum(dgamma(y, shape = 1/disp2, scale = scale2, log = TRUE)) +
        0.5 * t((b - mu)) %*% P %*% (b - mu)
    }
    
    
    f3<-function(b,y,x,mu,P,alpha=0,wt=1){
      mu2<-t(exp(alpha+x%*%b))
      t(x)%*%(t(1-y/mu2)*wt)+P%*%(b-mu)
    }
    
    f4<-function(b,y,x,alpha=0,wt=1,dispersion=1){
      disp2<-1/(wt/dispersion)
      
      (2*f1(b,y,x,alpha,wt/dispersion)+2*sum(dgamma(y,shape=1/disp2,scale=y*disp2,log=TRUE)))
      
    }
    
    #  f5<-f2_gamma
    #  f6<-f3_gamma
    
    ##########  This should be replaced!
    
    f7<-function(b,y,x,mu,P,alpha=0,wt=1){
      
      l2<-length(y)
      l1<-length(b)
      ltemp<-length(wt)
      if(ltemp==1){
        Ptemp<-wt*diag(l2)
      }
      else {
        Ptemp<-diag(wt)      
      }
      
      p<-1/(1+t(exp(-alpha-x%*%b)))
      
      
      Pout<-matrix(0,nrow=l1,ncol=l1)
      
      for(i in 1:l2){
        Pout<-Pout+Ptemp[i,i]*p[i]*(1-p[i])*x[i,]%*%t(x[i,])
      }
      
      Pout
    }
    
    
    
  }	
  
  out=list(f1=f1,f2=f2,f3=f3,f4=f4,
           #f5=f5,
           #f6=f6,
           f7=f7)
  
  out$call<-match.call()
  
  class(out)<-"glmbfamfunc"
  out
  
  
}


#' @rdname glmbfamfunc
#' @order 2
#' @method print glmbfamfunc
#' @export


print.glmbfamfunc<-function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nNegative Log-Likelihood Function:\nf1<-")
  print(x$f1)
  cat("\n(Negative Log-Posterior Function:\nf2<-")
  print(x$f2)
  cat("\nNegative Log-Posterior Gradient Function:\nf3<-")
  print(x$f3)
  cat("\nDeviance Function:\nf4<-")
  print(x$f4)
  
  
}


#' Standardize A Non-Gaussian Model
#'
#' Standardizes a Non-Gaussian Model prior to Envelope Creation
#' @param y a vector of observations of length m
#' @param x a design matrix of dimension m*p
#' @param P a positive-definite symmetric matrix specifying the prior precision matrix of the variables.
#' @param bstar a matrix containing the posterior mode from an optimization step
#' @param A1 a matrix containing the posterior precision matrix at the posterior mode
#' @return A list with the following components
#' \item{bstar2}{Standardized Posterior Mode}
#' \item{A}{Standardized Data Precision Matrix}
#' \item{x2}{Standardized Design Matrix}
#' \item{mu2}{Standardized Prior Mean vector}
#' \item{P2}{Standardized Precision Matrix Added to log-likelihood}
#' \item{L2Inv}{A matrix used when undoing the first step in standardization described below}
#' \item{L3Inv}{A matrix used when undoing the second step in standardization described below}
#' @details This functions starts with basic information about the model in the argument list and then
#' uses the following steps to further standardize the model (the model is already assumed to have a 0 prior mean vector
#' when this step is applied).
#' 
#' 1) An eigenvalue composition is applied to the posterior precision matrix, and the model is (as an interim step)
#' standardized to have a posterior precision matrix equal to the identity matrix. Please note that this means
#' that the prior precision matrix after this step is \code{"smaller"} than the identity matrix.
#' 
#' 
#' 2) A diagonal matrix epsilon is pulled out from the standardized prior precision matrix so that the remaining
#' part of the prior precision matrix still is positive definite. That part is then treated as part of the likelihood
#' for the rest of the standardization and simulation and only the part connected to epsilon is treated as part of the prior. 
#' Note that the exact epsilon chosen seems not to matter. Hence there are many possible ways of doing this 
#' standardization and future versions of this package may tweak the current approach 
#' if it helps improve numerical accuracy or acceptance rates.
#' 
#' 
#' 3) The model is next standardized (using a second eigenvalue decomposition) so that the prior (i.e., the portion connected to epsilon) is the identity 
#' matrix. The standardized model then simutaneously has the feature that the prior precision matrix is the 
#' identity matrix and that the data precision A (at the posterior mode) is a diagonal matrix. Hence the variables
#' in the standardized model are approximately independent at the posterior mode.
#' 
#' The steps here are based on the procedure descrived in \insertCite{Nygren2006}{glmbayes}.
#' 
#'  @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_glmb_Standardize_Model.R
#' @export


glmb_Standardize_Model<-function(y, x, P, bstar, A1){
  
  return(.glmb_Standardize_Model_cpp(y, x, P, bstar, A1))
  
}


#' Envelope Sizing and Optimization
#'
#' \code{EnvelopeSize()} is the high-level entry point that constructs
#' per-dimension grids and expected draw counts, while \code{EnvelopeOpt()}
#' performs the adaptive optimization used when \code{Gridtype = 2}.
#'
#' These functions implement the grid sizing logic used in envelope construction
#' for rejection sampling. They make use of the theory described in
#' \insertCite{Nygren2006}{glmbayes} and the general implementation outlined in
#' \insertCite{glmbayesSimmethods}{glmbayes}. 
#' 
#' @param a Numeric vector of diagonal precisions for the log-likelihood
#'   (posterior precision is \eqn{1 + a_i}).
#' @param G1 Numeric matrix of candidate grid points (3 * l1).
#' @param Gridtype Integer code controlling grid sizing logic:
#'   \itemize{
#'     \item 1 = static threshold test
#'     \item 2 = adaptive optimization via \code{EnvelopeOpt()}
#'     \item 3 = always three-point grid
#'     \item 4 = always single-point grid
#'   }
#' @param n Integer; number of posterior draws to generate (used for grid sizing).
#' @param n_envopt Integer; effective sample size passed to \code{EnvelopeOpt}.
#'   Defaults to -1, which means "use n".
#' @param use_opencl Logical; if \code{TRUE}, attempt GPU acceleration.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#'
#' @param a1 Numeric vector of diagonal elements of the data precision matrix
#'   (used by \code{EnvelopeOpt}).
#' @param core_cnt Integer; number of OpenCL cores or parallel workers available
#'   (default 1). When >1, envelope build cost is scaled down to reflect
#'   parallel construction.
#'
#' @section Gridtype Logic and Candidates per Draw:
#'
#' The envelope sizing logic follows the analysis of \insertCite{Nygren2006}{glmbayes}.
#'
#' \describe{
#'
#'   \item{Gridtype 1: Static Threshold}{
#'     For each dimension \eqn{i}, if
#'     \eqn{\sqrt{1 + a_i} \leq 2/\sqrt{\pi} \approx 1.128379},
#'     then a single tangent at the posterior mode suffices.
#'     Expected candidates per draw in that dimension:
#'     \eqn{\sqrt{1 + a_i}}.  
#'
#'     Otherwise, a symmetric three-point envelope is used at
#'     \eqn{\{\theta^\star_i - \omega_i, \theta^\star_i, \theta^\star_i + \omega_i\}},
#'     with expected candidates per draw bounded above by
#'     \eqn{2/\sqrt{\pi}}.
#'   }
#'
#'   \item{Gridtype 2: Adaptive Optimization}{
#'     Each dimension is assigned either a single-point or three-point envelope
#'     by minimizing
#'     \deqn{T_\mathrm{total}(g_i) = T_\mathrm{build}(g_i) + T_\mathrm{sample}(n, acc_i(g_i)).}
#'
#'     The optimizer balances build cost (grows with number of tangents) against
#'     sampling cost (decreases as acceptance improves).  
#'
#'     Expected candidates per draw:
#'     \eqn{\prod_j \text{scaleest}[i,j]}, where each factor is either
#'     \eqn{\sqrt{1+a_j}} (single-point) or \eqn{2/\sqrt{\pi}} (three-point),
#'     depending on the optimization outcome.
#'   }
#'
#'   \item{Gridtype 3: Always Three-Point}{
#'     Every dimension uses a symmetric three-point envelope.  
#'     Expected candidates per draw:
#'     \deqn{\left(\tfrac{2}{\sqrt{\pi}}\right)^k}
#'     for \eqn{k} dimensions, as shown in Theorem 3 of
#'     \insertCite{Nygren2006}{glmbayes}.
#'   }
#'
#'   \item{Gridtype 4: Always Single-Point}{
#'     Every dimension uses a single tangent at the posterior mode.  
#'     Expected candidates per draw:
#'     \deqn{\prod_{i=1}^k \sqrt{1 + a_i}}
#'     (Example 1 in \insertCite{Nygren2006}{glmbayes}).
#'   }
#'
#' }
#'
#' @details
#' \code{EnvelopeSize()} returns the constructed grid (\code{G2}),
#' index vectors (\code{GIndex1}), expected draw count (\code{E_draws}),
#' and the per-dimension grid index.  
#'
#' \code{EnvelopeOpt()} implements the adaptive optimization used in
#' \code{Gridtype = 2}, ranking dimensions by posterior variance and
#' promoting them to three-point tangents when the tradeoff is favorable.
#'
#' @return
#' \describe{
#'   \item{\code{EnvelopeSize()}}{A list with components \code{G2}, \code{GIndex1},
#'   \code{E_draws}, and \code{gridindex}.}
#'   \item{\code{EnvelopeOpt()}}{An integer vector of length \eqn{l1} with entries
#'   1 (single-point) or 3 (three-point).}
#' }
#'
#' @seealso \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeSort}}
#'
#' @references
#' \insertAllCited{}
#' @example inst/examples/Ex_EnvelopeOpt.R
#' @export
#' @usage EnvelopeSize(a, G1, Gridtype = 2L, n = 1000L, n_envopt = -1,
#'                     use_opencl = FALSE, verbose = FALSE)
#' @rdname EnvelopeSize
#' @export
EnvelopeSize <- function(a,
                         G1,
                         Gridtype   = 2L,
                         n          = 1000L,   # <-- updated default
                         n_envopt   = -1,
                         use_opencl = FALSE,
                         verbose    = FALSE) {
  .EnvelopeSize_cpp(a, G1, Gridtype, n, n_envopt, use_opencl, verbose)
}


#' @usage EnvelopeOpt(a1,n,core_cnt=1L)
#' @rdname EnvelopeSize
#' @export


EnvelopeOpt<-function(a1,n,core_cnt=1L){
  
  core_cnt <- as.integer(core_cnt)
  if (is.na(core_cnt) || core_cnt < 1L) core_cnt <- 1L
  
  
  a1rank<-rank(1/(1+a1))
  l1<-length(a1)
  
  dimcount<-matrix(0,(l1+1),l1)
  scaleest<-matrix(0,(l1+1),l1)
  intest<-c(1:(l1+1))
  slopeest<-c(1:(l1+1))
  
  dimcount[1,]<-diag(diag(l1))
  scaleest[1,]<-sqrt(1+a1)
  slopeest[1]<-prod(scaleest[1,])
  
  for(i in 2:(l1+1)){
    dimcount[i,]<-dimcount[i-1,]
    scaleest[i,]<-scaleest[i-1,]
    for(j in 1:l1){
      if(a1rank[j]==i-1){ 
        dimcount[i,j]<-3
        scaleest[i,j]<-2/sqrt(pi) 
      }
    }
    ##    intest[i]<-3^(i-1)
    intest[i]<-(3^(i-1))
    slopeest[i]<-prod(scaleest[i,])
  }
  evalest<-(intest/core_cnt)+n*slopeest
  minindex<-0
  for(j in 1:(l1+1)){if(evalest[j]==min(evalest)){minindex<-j}}
  
  ##  message("Estimated draws per Acceptance: ", slopeest[minindex])
  
  dimcount[minindex,]
  
}



#' GPU-Accelerated Envelope Construction for Posterior Simulation
#'
#' @name EnvelopeBuild
#'
#' @details
#' Constructs an enveloping function for posterior simulation using a grid of
#' tangency points. The envelope is used in accept-reject sampling to guarantee
#' iid draws from the posterior distribution. The implementation follows
#' \insertCite{Nygren2006}{glmbayes}, with extensions for GPU acceleration
#' (via OpenCL), dynamic grid optimization, and parallelized evaluation.
#'
#' The envelope is typically built around the posterior mode \eqn{\theta^\star} for a model in standard 
#' form (which in this context means a model with a diagonal posterior precision matrix
#' and prior identity precision matrix - \code{glmb_Standardize_Model}). It uses dimension-specific width parameters 
#' \eqn{\omega_i} derived from the precision matrix. Tangency points are selected per dimension, and the full grid is
#' formed via Cartesian expansion. Negative log-likelihood and gradient values
#' are computed at each grid point, either on CPU or GPU depending on the
#' \code{use_opencl} flag. These values are used to construct a piecewise
#' envelope function that dominates the posterior density.
#'
#' @section Models in standard form:
#'
#' Following Nygren & Nygren (2006, Section 3.3), the envelope construction
#' assumes the model has been reparameterized into *standard form*. In this form:
#' \itemize{
#'   \item The prior precision matrix is the identity.
#'   \item The posterior precision matrix, evaluated at the posterior mode
#'         \eqn{\theta^\star}, is diagonal.
#' }
#'
#' Example 2 in the paper illustrates the special case of a zero-mean normal
#' prior with identity covariance. In this setting:
#' \itemize{
#'   \item The generalized likelihood-subgradient density at a tangency point
#'         \eqn{\bar\theta} has mean vector \eqn{-c(\bar\theta)} and covariance \eqn{I}.
#'   \item The normalizing integrals over restricted sets factorize across
#'         dimensions. Specifically, for a rectangular set
#'         \eqn{A = \{\theta : l_L \leq \theta \leq l_U\}}, we have
#'         \deqn{
#'           \int_{\theta \in A} q^{\bar{\theta}}(\theta)\, d\theta
#'           \;=\;
#'           \prod_{r=1}^{p}
#'           \Big[
#'             \Phi\!\big(l_{U,r} + c_{r}(\bar{\theta})\big)
#'             - \Phi\!\big(l_{L,r} + c_{r}(\bar{\theta})\big)
#'           \Big],
#'         }
#'         and the truncated expectation in coordinate \eqn{r} is
#'         \deqn{
#'           \mathbb{E}_{\tilde{q}^{\bar{\theta}}}[\theta_r \mid \theta \in A]
#'           =
#'           -\,c_{r}(\bar{\theta})
#'           +
#'           \frac{
#'             \phi\!\big(l_{L,r} + c_{r}(\bar{\theta})\big)
#'             - \phi\!\big(l_{U,r} + c_{r}(\bar{\theta})\big)
#'           }{
#'             \Phi\!\big(l_{U,r} + c_{r}(\bar{\theta})\big)
#'             - \Phi\!\big(l_{L,r} + c_{r}(\bar{\theta})\big)
#'           }.
#'         }
#' }
#'
#' These closed-form expressions explain why \code{EnvelopeBuild} evaluates
#' \code{logU}, \code{loglt}, and \code{logrt} using univariate normal CDFs
#' and densities, rather than numerical integration. The gradients
#' (\code{cbars}) directly determine the shifted means of the restricted
#' densities, and the separability across dimensions makes the grid-based
#' construction computationally tractable.
#'
#' Models with Zellner's \eqn{g}-priors are essentially in standard form,
#' since in the whitened design space both the prior and likelihood precisions
#' are diagonal. Each dimension still needs to be scaled so that the prior
#' precision is exactly the identity matrix. For other models, standard form
#' can be achieved by reparameterization (e.g. via Cholesky of the posterior
#' precision) or by shifting part of the prior quadratic form into the
#' likelihood.

#' @section Construction of restricted subgradient densities:
#'
#' Following Remark 5 in Nygren & Nygren (2006), each unrestricted
#' likelihood-subgradient density \eqn{q_{\bar{\theta}}(\cdot)} can be
#' restricted to a subset \eqn{A \subset \Theta}. The restricted density is
#' defined as
#' \deqn{
#'   \tilde{q}_{\bar{\theta}}(\theta)
#'   =
#'   \frac{ q_{\bar{\theta}}(\theta)\,\mathbf{1}_{\{\theta \in A\}} }
#'        { \int_{\theta' \in A} q_{\bar{\theta}}(\theta')\, d\theta' },
#' }
#' and the corresponding constant is
#' \deqn{
#'   \tilde{a}(\theta_bar)
#'   =
#'   a(\theta_bar)
#'   \int_{\theta \in A} q_{\theta_bar}(\theta)\, d\theta,
#' }
#' where \eqn{a(\theta_bar)} is the global normalizing constant from
#' Theorem 1. For every \eqn{\theta \in A}, the identity
#' \deqn{
#'   \tilde{a}(\theta_bar) \cdot h_{\theta_bar}(\theta) \cdot
#'   \tilde{q}_{\theta_bar}(\theta)
#'   = \pi(\theta \mid y)
#' }
#' holds, ensuring that the restricted densities reproduce the posterior
#' when combined.
#'
#' The constant \eqn{a(\bar{\theta})} is defined in Theorem 1 as
#' \deqn{
#'   a(\bar{\theta}) =
#'   \frac{g(\bar{\theta})\,\mathrm{MGF}\!\big(-c(\bar{\theta})\big)}
#'        {f(y)\,\exp\!\big(-c(\bar{\theta})^{T}\bar{\theta}\big)},
#' }
#' where \eqn{g(\bar{\theta})} is the reference density at the tangency point,
#' \eqn{c(\bar{\theta})} is the subgradient of the log-likelihood, and
#' \eqn{\mathrm{MGF}(-c(\bar{\theta}))} is the moment-generating function of
#' the prior evaluated at \eqn{-c(\bar{\theta})}.
#'
#' The envelope function \eqn{h_{\bar{\theta}}(\theta)} is defined as
#' \deqn{
#'   h_{\bar{\theta}}(\theta) =
#'   \frac{\exp\!\big(-c(\bar{\theta})^{T}\bar{\theta}\big)\,f(y\mid\theta)}
#'        {\exp\!\big(-c(\bar{\theta})^{T}\theta\big)\,g(\bar{\theta})},
#' }
#' and satisfies
#' \deqn{
#'   0 \le h_{\bar{\theta}}(\theta) \le 1
#'   \quad \forall\,\theta \in \Theta,
#'   \qquad
#'   h_{\bar{\theta}}(\bar{\theta}) = 1
#'   \quad \text{if } f(y\mid\bar{\theta}) = g(\bar{\theta}).
#' }
#'
#' In the standardized model (zero-mean normal prior with identity covariance
#' and diagonal posterior precision at the mode), the restricted integral
#' \eqn{\int_{A} q_{\bar{\theta}}(\theta)\, d\theta} factorizes across
#' dimensions and can be evaluated in closed form using normal CDFs:
#' \deqn{
#'   \int_{\theta \in A} q_{\bar{\theta}}(\theta)\, d\theta
#'   =
#'   \prod_{r=1}^{p}
#'   \Big[
#'     \Phi\!\big(l_{U,r} + c_{r}(\bar{\theta})\big)
#'     - \Phi\!\big(l_{L,r} + c_{r}(\bar{\theta})\big)
#'   \Big].
#' }
#'
#' In standardized models, where the prior is \eqn{\mathcal{N}(0, I)}, the
#' moment-generating function simplifies to
#' \deqn{
#'   \mathrm{MGF}\big(-c(\bar{\theta})\big)
#'   = \exp\left(\tfrac{1}{2} c(\bar{\theta})^{T} c(\bar{\theta})\right),
#' }
#' so the constant \eqn{a(\bar{\theta})} from Theorem 1 becomes
#' \deqn{
#'   a(\bar{\theta}) =
#'   \frac{g(\bar{\theta})\,\exp\left(\tfrac{1}{2} c(\bar{\theta})^{T} c(\bar{\theta})\right)}
#'        {f(y)\,\exp\!\big(-c(\bar{\theta})^{T} \bar{\theta}\big)}.
#' }
#'
#' This is why \code{EnvelopeBuild} computes and stores \code{logU},
#' \code{loglt}, and \code{logrt} using univariate normal CDF evaluations.
#' The constants \eqn{\tilde{a}(\bar{\theta})} are then obtained by scaling
#' the global constant with these integrals, and the mixture weights
#' (\code{PLSD}) are normalized accordingly. In practice:
#' \itemize{
#'   \item \code{EnvelopeSet_Grid_C2_pointwise} evaluates restricted densities
#'         at each grid point.
#'   \item \code{LLconst} stores the log of the restricted integrals.
#'   \item \code{EnvelopeSet_LogP_C2} computes \eqn{\tilde{a}(\bar{\theta})} and
#'         normalizes mixture weights.
#' }
#' 
#' @section Mixture construction and tractable probabilities:
#'
#' Claim 2 in Nygren & Nygren (2006) shows that the posterior density
#' \eqn{\pi(\theta \mid y)} can be expressed as a mixture of restricted
#' likelihood-subgradient densities. Let \eqn{A_1, \dots, A_m} be a partition
#' of the parameter space \eqn{\Theta}, and define
#' \deqn{
#'   \tilde{q}_{\bar{\theta}}(\theta)
#'   =
#'   \sum_{i=1}^{m} \tilde{p}_i \, q^{\bar{\theta}}_{A_i}(\theta),
#'   \qquad
#'   \tilde{p}_i=  \frac{ \tilde{a}_i }{ \sum_{j=1}^{k} \tilde{a}_j } 
#'   }
#' 
#' 
#'
#' In Remark 6 of Nygren & Nygren (2006), the mixture weights
#' \eqn{\tilde{p}_i} for each restricted likelihood-subgradient density
#' \eqn{q^{\bar{\theta}_i}_{A_i}} are defined as
#' \deqn{
#'   \tilde{p}_i
#'   =
#'   \frac{
#'     g(\bar{\theta}_i)\,\mathrm{MGF}\big(-c(\bar{\theta}_i)\big)
#'     \int_{\theta \in A_i} q^{\bar{\theta}_i}(\theta)\, d\theta
#'     \,/\, \exp\big(-c(\bar{\theta}_i)^{T} \bar{\theta}_i\big)
#'   }{
#'     \sum_{j=1}^{k}
#'     g(\bar{\theta}_j)\,\mathrm{MGF}\big(-c(\bar{\theta}_j)\big)
#'     \int_{\theta \in A_j} q^{\bar{\theta}_j}(\theta)\, d\theta
#'     \,/\, \exp\big(-c(\bar{\theta}_j)^{T} \bar{\theta}_j\big)
#'   }.
#' }
#' This expression reflects the full normalization of the mixture, where each
#' \eqn{\tilde{p}_i} is proportional to the restricted constant
#' \eqn{\tilde{a}_i(\bar{\theta}_i)} from Remark 5, and the denominator sums
#' over all such constants across the partition. The resulting mixture
#' \eqn{\tilde{q}^{\bar{\theta}}(\theta)} is a valid approximation to the
#' posterior density \eqn{\pi(\theta \mid y)}.
#'
#' Remark 6 emphasizes that these mixture weights are tractable to compute
#' in standardized models. When the prior is \eqn{\mathcal{N}(0, I)} and the
#' posterior precision is diagonal at the mode, each integral
#' \eqn{\int_{A_i} q^{\bar{\theta}}(\theta)\, d\theta} factorizes across
#' dimensions and can be evaluated using normal CDFs:
#' \deqn{
#'   \int_{\theta \in A_i} q^{\bar{\theta}}(\theta)\, d\theta
#'   =
#'   \prod_{r=1}^{p}
#'   \Big[
#'     \Phi\!\big(l_{U,r}^{(i)} + c_{r}(\bar{\theta})\big)
#'     - \Phi\!\big(l_{L,r}^{(i)} + c_{r}(\bar{\theta})\big)
#'   \Big],
#' }
#' where \eqn{l_{L}^{(i)}} and \eqn{l_{U}^{(i)}} are the bounds defining
#' region \eqn{A_i}.
#'
#' This tractability is central to the envelope construction. It allows
#' \code{EnvelopeBuild} to compute:
#' \itemize{
#'   \item \code{LLconst}, which stores the log of each restricted integral
#'         \eqn{\log \int_{A_i} q^{\bar{\theta}}(\theta)\, d\theta}.
#'   \item \code{PLSD}, which stores the normalized mixture weights
#'         \eqn{\tilde{p}_i}.
#' }
#' These quantities are used to construct the combined density
#' \eqn{\tilde{q}_{\bar{\theta}}(\theta)} and to evaluate the envelope
#' approximation to the posterior. Because all components are normalized and
#' tractable, the mixture is both valid and computationally efficient.
#' @section Log-scale properties of the envelope function:
#'
#' The envelope function \eqn{h_{\bar{\theta}}(\theta)} is defined in
#' Theorem 1 as
#' \deqn{
#'   h_{\bar{\theta}}(\theta) =
#'   \frac{\exp\!\big(-c(\bar{\theta})^{T} \bar{\theta}\big)\,f(y \mid \theta)}
#'        {\exp\!\big(-c(\bar{\theta})^{T} \theta\big)\,g(\bar{\theta})}.
#' }
#' When \eqn{g(\bar{\theta}) = f(y \mid \bar{\theta})}, this simplifies to
#' \deqn{
#'   h_{\bar{\theta}}(\theta) =
#'   \exp\!\big( c(\bar{\theta})^{T}(\theta - \bar{\theta}) \big)
#'   \cdot \frac{f(y \mid \theta)}{f(y \mid \bar{\theta})}.
#' }
#' Taking logarithms yields
#' \deqn{
#'   \log h_{\bar{\theta}}(\theta)
#'   =
#'   c(\bar{\theta})^{T}(\theta - \bar{\theta})
#'   + \log f(y \mid \theta) - \log f(y \mid \bar{\theta}),
#' }
#' which is tractable as long as the log-likelihood \eqn{\log f(y \mid \theta)}
#' is. In particular, if the log-likelihood is concave or piecewise affine,
#' then \eqn{\log h_{\bar{\theta}}(\theta)} inherits that structure and can be
#' efficiently evaluated across grid regions.
#'
#' This tractability is central to the envelope construction: it allows
#' \code{EnvelopeBuild} to evaluate the envelope function pointwise using
#' \code{EnvelopeSet_Grid_C2_pointwise}, and ensures that the resulting approximation
#' remains bounded between 0 and 1. At the tangency point, we recover
#' \deqn{
#'   h_{\bar{\theta}}(\bar{\theta}) = 1,
#' }
#' confirming that the envelope touches the posterior density exactly.
#' The key inequality that ensures envelope dominance follows from the
#' subgradient inequality for concave functions. If \eqn{\log f(y \mid \theta)}
#' is concave and \eqn{c(\bar{\theta})} is a subgradient at \eqn{\bar{\theta}},
#' then
#' \deqn{
#'   \log f(y \mid \theta)
#'   \le \log f(y \mid \bar{\theta}) + c(\bar{\theta})^{T}(\theta - \bar{\theta}),
#' }
#' which implies
#' \deqn{
#'   \log h_{\bar{\theta}}(\theta)
#'   \le 0,
#'   \qquad
#'   h_{\bar{\theta}}(\theta) \le 1.
#' }
#' This inequality guarantees that the envelope function dominates the posterior
#' density pointwise, as required by Theorem 1. Equality holds at the tangency
#' point \eqn{\theta = \bar{\theta}}, where
#' \deqn{
#'   h_{\bar{\theta}}(\bar{\theta}) = 1.
#' }
#' @section Use of the envelope during sampling:
#'
#' The functions \code{rnnorm_reg_std_cpp()} and \code{rnnorm_reg_std_cpp_parallel()}
#' use the envelope to generate posterior samples via rejection sampling. Although not exported,
#' these functions are called internally by \code{rnnorm_reg_cpp()}, which in turn is invoked by
#' the user-facing function \code{rNormal_reg()}. Together, these routines implement
#' envelope-based sampling for generalized linear models with log-concave likelihood functions
#' and multivariate normal priors.
#'
#' The envelope provides a mixture of restricted likelihood-subgradient densities,
#' each defined over a region \eqn{A_i}, with associated mixture weights
#' \eqn{\tilde{p}_i} stored in \code{PLSD}. The sampling proceeds as follows:
#'
#' \enumerate{
#'   \item A region index \eqn{J(i)} is drawn from the discrete distribution
#'         defined by \code{PLSD}.
#'   \item A candidate \eqn{\theta_i} is drawn from the restricted density
#'         \eqn{q^{\bar{\theta}_{J(i)}}_{A_{J(i)}}}, using the normal CDF bounds
#'         \code{loglt} and \code{logrt}, and subgradient vector \code{cbars}.
#'         Simulation for each dimension uses the internal C++ function
#'         \code{ctrnorm_cpp()}, which explicitly uses these inputs.
#'   \item The log-likelihood \eqn{\log f(y \mid \theta_i)} is computed and
#'         stored in \code{testll[0]} using the appropriate likelihood function
#'         \code{f2}.
#' }
#'
#' The acceptance test is performed using the inequality
#' \deqn{
#'   \log(U_2) \le \mathrm{LLconst}[J(i)] + \mathrm{cbars}[J(i), ]^{T} \theta_i
#'   + \log f(y \mid \theta_i),
#' }
#' which is equivalent to
#' \deqn{
#'   \log(U_2) \le \log f(y \mid \theta_i) - \left( \log f(y \mid \bar{\theta}_{J(i)}) - c(\bar{\theta}_{J(i)})^{T}(\theta_i - \bar{\theta}_{J(i)}) \right),
#' }
#' 
#' where:
#' \itemize{
#'   \item \code{LLconst[J(i)]} stores the precomputed quantity
#'         \eqn{-\log f(y \mid \bar{\theta}_{J(i)}) - c(\bar{\theta}_{J(i)})^{T} \bar{\theta}_{J(i)}},
#'         computed during envelope construction via \code{EnvelopeSet_LogP_C2()}.
#'   \item \code{cbars[J(i), ]} is the precomputed subgradient vector
#'         \eqn{c(\bar{\theta}_{J(i)})}, extracted via \code{cbars(J(i), _)}.
#'         It defines the exponential tilt direction used to evaluate the envelope.
#'   \item \code{testll[0]} is the log-likelihood at the candidate draw
#'         \eqn{\theta_i}, evaluated using the model specified by
#'         \code{family} and \code{link}.
#'   \item \eqn{-\log(U_2)} is the threshold from a uniform draw
#'         \eqn{U_2 \sim \mathrm{Unif}(0,1)}.
#' }
#'
#' The right-hand side of this inequality is always non-positive, and equals zero
#' when \eqn{\theta_i = \bar{\theta}_{J(i)}}. This reflects the fact that the envelope
#' is tangent to the log-likelihood at each \eqn{\bar{\theta}_j}, and lies above it elsewhere.
#'
#' This procedure guarantees that accepted samples are drawn from the posterior
#' \eqn{\pi(\theta \mid y)}. The envelope ensures bounded rejection probability,
#' and the mixture structure allows efficient sampling across regions. The output
#' \code{out} contains accepted draws, and \code{draws} records the number of
#' attempts per sample.
#'
#' The components returned by \code{EnvelopeBuild()} are used in specific steps of the
#' sampling procedure as follows:
#' \itemize{
#'   \item \code{PLSD} is used to randomly select a region index \eqn{J(i)} from the envelope mixture.
#'   \item \code{loglt} and \code{logrt} define the truncated normal bounds for each dimension,
#'         used together with \code{cbars} to generate candidate values \eqn{\theta_i}.
#'   \item \code{cbars} provides the subgradient vectors \eqn{c(\bar{\theta}_j)} used both for
#'         candidate generation and for computing the acceptance test.
#'   \item \code{LLconst} stores precomputed constants used in the acceptance inequality,
#'         avoiding recomputation of posterior terms at tangency points.
#'   \item \code{logU} stores the per-dimension log-density contributions for each region,
#'         computed during envelope setup. These values are summed to produce \code{logP},
#'         which determines the mixture weights \code{PLSD}.
#'   \item \code{logP} contains the total log-probabilities for each grid component,
#'         which are normalized to form the mixture weights \code{PLSD}.
#'   \item \code{thetabars} stores the tangency points \eqn{\bar{\theta}_j} used to define
#'         subgradients and region-specific densities.
#'   \item \code{GridIndex} encodes the sampling type (tail, center, line) used for each
#'         dimension and region, guiding how each coordinate is simulated.
#' }
#' @section Algorithmic steps (linked to theory):
#'
#' The implementation of \code{EnvelopeBuild} follows the envelope construction
#' in Nygren & Nygren (2006) for models in standard form (See Section 3-3.3). 
#' Each computational step corresponds to a theoretical guarantee:
#'
#' 1. **Compute width parameters \eqn{\omega_i} from the diagonal precision matrix.**  In particular,
#' let \eqn{\theta^{\ast}} denote the unique posterior mode. For each dimension \eqn{i},
#'  define
#' \deqn{
#'   \omega_{i} :=
#'   \frac{\sqrt{2} - \exp\!\big(-1.20491 - 0.7321\,\sqrt{0.5 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2}}\big)}
#'        {\sqrt{1 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2}}}.
#' }
#' 
#'    As seen from the above, the widths \eqn{\omega_i} are derived from the local curvature of the 
#'    log-likelihood at the posterior mode. This ensures that the three-interval construction per
#'    dimension below yields an envelope whose efficiency does not deteriorate with sample size.
#'
#' 2. **Use the width parameters to construct intervals around the posterior mode \eqn{\theta^\star}.**  Specifically,
#' we set
#' \deqn{
#'   \ell_{i,1} = \theta^{\ast}_{i} - 0.5\,\omega_{i}, \quad
#'   \ell_{i,2} = \theta^{\ast}_{i} + 0.5\,\omega_{i},
#' }
#' and construct three intervals per dimension:
#' \deqn{
#'   A_{i,1} = (-\infty,\ell_{i,1}), \quad
#'   A_{i,2} = [\ell_{i,1},\ell_{i,2}], \quad
#'   A_{i,3} = (\ell_{i,2},\infty).
#' }
#'
#' For each dimension \eqn{i}, let \eqn{J_{i} = \{1,2,3\}} and define
#' \eqn{J = \prod_{i=1}^{p} J_{i}}, which has \eqn{3^{p}} elements. Each
#' \eqn{j \in J} is a vector \eqn{(j_{1},\ldots,j_{p})}, and we define
#' \deqn{
#'   A^{\ast}_{j} = \prod_{i=1}^{p} A_{i,j_{i}}.
#' }
#' The collection \eqn{A^{\ast} = \{A^{\ast}_{j} : j \in J\}} forms a partition of \eqn{\Theta}.
#'
#'    
#' 3. **For each member of the partition, select tangency points \eqn{\theta^\star \pm \omega_i}.**
#' 
#'  For each \eqn{j \in J}, define index sets
#' \deqn{
#'   C_{j1} = \{i : j_{i} = 1\}, \quad
#'   C_{j2} = \{i : j_{i} = 2\}, \quad
#'   C_{j3} = \{i : j_{i} = 3\}.
#' }
#' 
#' The tangency points \eqn{\bar{\theta}_{j}} are then defined componentwise by
#' \deqn{
#'   \bar{\theta}_{j,i} =
#'   \begin{cases}
#'     \theta^{\ast}_{i} - \omega_{i}, & i \in C_{j1}, \\
#'     \theta^{\ast}_{i},              & i \in C_{j2}, \\
#'     \theta^{\ast}_{i} + \omega_{i}, & i \in C_{j3}.
#'   \end{cases}
#' }
#'   
#'    The tangency points are hence chosen so that the envelope touches the log-likelihood at
#'    representative points in each interval, guaranteeing dominance and tightness.
#'
#' 4. **Build the full grid of tangency points (Cartesian product across dimensions).**  
#'    
#'    The Cartesian product of per-dimension partitions yields the \eqn{3^p} restricted
#'    densities described in the paper, ensuring coverage of the full parameter space.
#'
#' 5. **Evaluate negative log-likelihood and gradients at each grid point to construct the 
#'    likelihood subgradient densities and to facilitate accept rejection sampling**  
#'      
#'    The subgradients \eqn{c(\bar{\theta})} are part of the definitions of Likelihood subgradient densities 
#'    (see definition 2 below) while both the subgradients and the negative log-likelihoods (through the function 
#'    \eqn{h_{\bar{\theta}}(.)}) are used during the accept-reject procedure. 
#'    CPU and GPU routines compute these values efficiently.
#'    - On CPU: via \code{f2_*} and \code{f3_*} routines.
#'    - On GPU: via \code{f2_f3_opencl}, which computes both together
#'    
#' 6. **Call \code{EnvelopeSet_Grid_C2_pointwise} to evaluate restricted multivariate normal log-densities.**  
#'    (Claim 2; Remark 5).  
#'    Each restricted density corresponds to a subset of the partition, normalized
#'    as in Remark 5.
#'
#' 7. **Call \code{EnvelopeSet_LogP_C2} to compute component log-probabilities and constants.**  
#'    (Remark 6).  
#'    The constants \eqn{\tilde{a}} and mixture weights \eqn{\tilde{p}_i} are computed
#'    explicitly as in Remark 6, ensuring that the mixture envelope is properly normalized.
#'
#' 8. **Normalize probabilities (\code{PLSD}) and optionally sort grid components.**  
#'    (Claim 2).  
#'    Normalization ensures that the mixture of restricted densities forms a valid
#'    dominating density for the posterior. Sorting is an implementation detail to
#'    improve sampling efficiency.
#' 
#' 
#' 
#' 
#' @section Formal definition and key claim:
#' \strong{Definition 2.} A probability density function \eqn{q(\cdot)} is a
#' generalized likelihood-subgradient probability density for a posterior density
#' \eqn{\pi(\cdot\mid y)} with prior density \eqn{\pi(\cdot)} and likelihood
#' function \eqn{f(y\mid\cdot)} at a point \eqn{\bar{\theta} \in \Theta} if there
#' exists a subgradient \eqn{c(\bar{\theta})} for the negative of the log of a
#' function \eqn{g} at \eqn{\bar{\theta}} such that:
#' \itemize{
#'   \item (a) \eqn{g(\cdot)} bounds \eqn{f(y\mid\cdot)} from above,
#'   \item (b) \eqn{\mathrm{MGF}(-c(\bar{\theta})) = \int_{\Theta} \exp\!\big(-c(\bar{\theta})^{T}\theta\big)\,\pi(\theta)\,d\theta} is finite,
#'   \item (c) \eqn{\forall\,\theta \in \Theta:~ q(\theta) = \exp\!\big(-c(\bar{\theta})^{T}\theta\big)\,\pi(\theta)\,/\,\mathrm{MGF}\!\big(-c(\bar{\theta})\big)}.
#' }
#'
#' \emph{Special cases:}
#' \itemize{
#'   \item If \eqn{g(\theta)} is the likelihood function, then we call \eqn{q(\theta)} a
#'         likelihood-subgradient density. Log-concave likelihood functions are ubiquitous
#'         in statistical modeling. Models with log-concave likelihood functions include
#'         the Poisson and logit regression models, as well as some survival models.
#'   \item If \eqn{g(\theta) = \bar{f} \ge f(y\mid\theta)} for all \eqn{\theta}, then we are in the
#'         Bayesian context discussed earlier and the prior is a generalized likelihood-subgradient density.
#' }
#'
#' The Appendix provides a more detailed discussion on the existence of likelihood-subgradient densities.
#' These are guaranteed to exist at every point if the prior is a finite mixture of multivariate normals
#' and the likelihood function is log-concave. The generalized likelihood-subgradient density is particularly
#' simple in the case of finite mixtures of multivariate normal priors.
#'
#' \strong{Claim 1.} Suppose that the prior \eqn{\pi(\cdot)} in Definition 2 is a finite mixture
#' of multivariate normals \eqn{\sum_{i=1}^{k} p_{i}\,\pi_{i}(\cdot\mid \mu_{i}, \Sigma_{i})}. Let
#' \eqn{c(\bar{\theta})} be a subgradient for \eqn{-\ln(g(\cdot))} at \eqn{\bar{\theta}}. Then
#' \eqn{\mathrm{MGF}(-c(\bar{\theta}))} is finite and is in the form:
#' \deqn{
#'   \mathrm{MGF}\!\big(-c(\bar{\theta})\big)
#'   = \sum_{i=1}^{k} p_{i}\,\exp\!\Big(-c(\bar{\theta})^{T}\mu_{i} + \tfrac{1}{2}\,c(\bar{\theta})^{T}\Sigma_{i}\,c(\bar{\theta})\Big).
#' }
#'
#' The corresponding generalized likelihood-subgradient density is also a mixture of multivariate normals in the form:
#' \deqn{
#'   q(\theta) = \sum_{i=1}^{k} \tilde{p}_{i}\,\pi_{i}\!\big(\theta \mid \tilde{\mu}_{i}, \Sigma_{i}\big),
#' }
#' where
#' \deqn{
#'   \tilde{p}_{i}
#'   = \frac{p_{i}\,\exp\!\Big(-c(\bar{\theta})^{T}\mu_{i} + \tfrac{1}{2}\,c(\bar{\theta})^{T}\Sigma_{i}\,c(\bar{\theta})\Big)}
#'          {\mathrm{MGF}\!\big(-c(\bar{\theta})\big)}
#'   \quad\text{and}\quad
#'   \tilde{\mu}_{i} = \mu_{i} - \Sigma_{i}\,c(\bar{\theta}).
#' }
#'
#' @section Theorem 1 (envelope dominance and equality at tangency):
#' Let \eqn{q_{\bar{\theta}}(\cdot)} be a generalized likelihood-subgradient
#' density as in Definition 2. Define:
#' \deqn{
#'   a(\bar{\theta}) =
#'   \frac{g(\bar{\theta})\,\mathrm{MGF}\!\big(-c(\bar{\theta})\big)}
#'        {f(y)\,\exp\!\big(-c(\bar{\theta})^{T}\bar{\theta}\big)}
#' }
#' and
#' \deqn{
#'   h_{\bar{\theta}}(\theta) =
#'   \frac{\exp\!\big(-c(\bar{\theta})^{T}\bar{\theta}\big)\,f(y\mid\theta)}
#'        {\exp\!\big(-c(\bar{\theta})^{T}\theta\big)\,g(\bar{\theta})}.
#' }
#'
#' Then
#' \deqn{
#'   a(\bar{\theta})\,q_{\bar{\theta}}(\theta)
#'   \;\ge\;
#'   a(\bar{\theta})\,h_{\bar{\theta}}(\theta)\,q_{\bar{\theta}}(\theta)
#'   \;=\;
#'   \pi(\theta\mid y),
#' }
#' and
#' \deqn{
#'   0 \le h_{\bar{\theta}}(\theta) \le 1
#'   \quad \forall\,\theta \in \Theta.
#' }
#'
#' Finally, if \eqn{f(y\mid\bar{\theta}) = g(\bar{\theta})}, then
#' \deqn{
#'   h_{\bar{\theta}}(\bar{\theta}) = 1.
#' }
#'
#' @section Claim 2 (mixtures over a partition) and Remark 5 (restricted densities):
#' \strong{Claim 2.} Let \eqn{A_{1},A_{2},\ldots,A_{k}} be a finite partition of \eqn{\Theta} and
#' let \eqn{\tilde{q}_{1}(\cdot),\tilde{q}_{2}(\cdot),\ldots,\tilde{q}_{k}(\cdot)} be associated restricted densities
#' such that there exist associated constants \eqn{\tilde{a}_{1},\tilde{a}_{2},\ldots,\tilde{a}_{k}} and functions
#' \eqn{\tilde{h}_{1}(\cdot),\tilde{h}_{2}(\cdot),\ldots,\tilde{h}_{k}(\cdot)} satisfying the following:
#' \itemize{
#'   \item (a) \eqn{\theta \in A_{i} \Rightarrow \tilde{a}_{i}\,\tilde{h}_{i}(\theta)\,\tilde{q}_{i}(\theta) = \pi(\theta \mid y)}, \eqn{i=1,2,\ldots,k}.
#'   \item (b) \eqn{\theta \in A_{i} \Rightarrow 0 \le \tilde{h}_{i}(\theta) \le 1}, \eqn{i=1,2,\ldots,k}.
#' }
#' Define a new density by
#' \deqn{
#'   \tilde{q}(\theta) = \sum_{i=1}^{k} \tilde{p}_{i}\,\tilde{q}_{i}(\theta),
#' }
#' where \eqn{\tilde{p}_{i} = \tilde{a}_{i} \big/ \big(\sum_{j=1}^{k} \tilde{a}_{j}\big)}.
#' Let \eqn{\tilde{a} = \sum_{j=1}^{k} \tilde{a}_{j}}, and let \eqn{\tilde{h}(\theta)} be a function satisfying
#' \eqn{\theta \in A_{i} \Rightarrow \tilde{h}(\theta) = \tilde{h}_{i}(\theta)}, \eqn{i=1,2,\ldots,k}.
#' Then
#' \deqn{
#'   \tilde{a}\,\tilde{q}(\theta) \;\ge\; \tilde{a}\,\tilde{h}(\theta)\,\tilde{q}(\theta) \;=\; \pi(\theta \mid y)
#' }
#' and
#' \deqn{
#'   0 \le \tilde{h}(\theta) \le 1 \quad \forall\,\theta \in \Theta.
#' }
#' In other words, a bounding function for the full space can be constructed by combining bounding
#' functions for the individual elements of the partition. Our main interest is in mixtures of restricted
#' generalized likelihood-subgradient densities. \emph{Remark 5} shows how to use generalized
#' likelihood-subgradient densities to construct restricted densities of the form required for Claim 2.
#'
#' \strong{Remark 5.} Let \eqn{q_{\bar{\theta}}(\cdot)} be a generalized likelihood-subgradient density as in Theorem 1.
#' Define a restricted density \eqn{\tilde{q}_{\bar{\theta}}(\cdot)} and a corresponding constant \eqn{\tilde{a}(\bar{\theta})}
#' on a set \eqn{A \subset \Theta} by
#' \deqn{
#'   \tilde{q}_{\bar{\theta}}(\theta) =
#'   \frac{q_{\bar{\theta}}(\theta)}{\int_{\theta \in A} q_{\bar{\theta}}(\theta)\,d\theta}
#' }
#' and
#' \deqn{
#'   \tilde{a}(\bar{\theta}) =
#'   a(\bar{\theta}) \int_{\theta \in A} q_{\bar{\theta}}(\theta)\,d\theta.
#' }
#' Then for every \eqn{\theta \in A}, we have
#' \deqn{
#'   \tilde{a}(\bar{\theta})\,h_{\bar{\theta}}(\theta)\,\tilde{q}_{\bar{\theta}}(\theta) = \pi(\theta \mid y).
#' }
#'
#'
#' @section Remark 6 (explicit mixture weights and constants):
#' Suppose that the restricted densities in Claim 2 are restricted generalized
#' likelihood-subgradient densities as in Remark 5. Then both the overall
#' constant \eqn{\tilde{a}} and the mixture probabilities \eqn{\tilde{p}_{i}}
#' admit explicit formulas:
#'
#' \deqn{
#'   \tilde{a}
#'   = \frac{1}{f(y)}
#'     \sum_{j=1}^{k}
#'       g(\bar{\theta}_{j})\,
#'       \mathrm{MGF}\!\big(-c(\bar{\theta}_{j})\big)\,
#'       \int_{\theta \in A_{j}} q_{\bar{\theta}_{j}}(\theta)\,d\theta\,
#'       \exp\!\big(-c(\bar{\theta}_{j})^{T}\bar{\theta}_{j}\big)
#' }
#'
#' and
#'
#' \deqn{
#'   \tilde{p}_{i}
#'   =
#'   \frac{
#'     g(\bar{\theta}_{i})\,
#'     \mathrm{MGF}\!\big(-c(\bar{\theta}_{i})\big)\,
#'     \int_{\theta \in A_{i}} q_{\bar{\theta}_{i}}(\theta)\,d\theta\,
#'     \exp\!\big(-c(\bar{\theta}_{i})^{T}\bar{\theta}_{i}\big)
#'   }{
#'     \sum_{j=1}^{k}
#'       g(\bar{\theta}_{j})\,
#'       \mathrm{MGF}\!\big(-c(\bar{\theta}_{j})\big)\,
#'       \int_{\theta \in A_{j}} q_{\bar{\theta}_{j}}(\theta)\,d\theta\,
#'       \exp\!\big(-c(\bar{\theta}_{j})^{T}\bar{\theta}_{j}\big)
#'   }.
#' }
#'
#' The expression for \eqn{\tilde{a}} provides insight into how the partition
#' \eqn{\{A_{j}\}} and the positioning of the tangency points
#' \eqn{\bar{\theta}_{j}, j=1,2,\ldots,k}, affect its value.
#' In fact, the optimal placement of the tangencies for each element of the
#' partition should satisfy the property that the tangency points are the
#' expectations of the resulting restricted likelihood-subgradient densities.

#' @section Example 2 (standard normal prior, restricted set):
#' Suppose that the prior \eqn{\pi(\cdot)} in Definition 2 is a
#' \eqn{p}-dimensional multivariate normal density with mean vector 0
#' and variance-covariance matrix \eqn{I}, the identity matrix.
#' Let \eqn{q_{\bar{\theta}}(\cdot)} be a generalized likelihood-subgradient
#' density at \eqn{\bar{\theta}}.
#'
#' It is straightforward to verify that \eqn{q_{\bar{\theta}}(\cdot)} has
#' mean vector \eqn{-c(\bar{\theta})} and variance-covariance matrix \eqn{I}.
#'
#' Define a restricted set
#' \deqn{
#'   A = \{\theta \in \Theta : \ell^{L} \le \theta \le \ell^{U}\},
#' }
#' for some vectors \eqn{\ell^{L}} and \eqn{\ell^{U}}.
#'
#' Then
#' \deqn{
#'   \int_{\theta \in A} q_{\bar{\theta}}(\theta)\,d\theta
#'   = \prod_{r=1}^{p} \Big[ \Phi(\ell^{U}_{r} + c_{r}(\bar{\theta}))
#'                          - \Phi(\ell^{L}_{r} + c_{r}(\bar{\theta})) \Big],
#' }
#' where \eqn{\Phi(\cdot)} denotes the standard normal cumulative distribution function.
#'
#' The expectation under the restricted density is
#' \deqn{
#'   \mathbb{E}_{\tilde{q}_{\bar{\theta}}}[\theta]
#'   = -c(\bar{\theta})
#'     - \lambda\!\big(\ell^{L} + c(\bar{\theta}),\,\ell^{U} + c(\bar{\theta})\big),
#' }
#' where \eqn{\lambda(\cdot)} is a vector-valued function with
#' \eqn{r}th component given by Mills' ratio:
#' \deqn{
#'   \lambda_{r}
#'   = \frac{\varphi(\ell^{L}_{r} + c_{r}(\bar{\theta}))
#'          - \varphi(\ell^{U}_{r} + c_{r}(\bar{\theta}))}
#'          {\Phi(\ell^{U}_{r} + c_{r}(\bar{\theta}))
#'          - \Phi(\ell^{L}_{r} + c_{r}(\bar{\theta}))},
#' }
#' with \eqn{\varphi(\cdot)} the standard normal density.
#'
#' This example illustrates that in the standard normal prior case,
#' the generalized likelihood-subgradient density remains normal with
#' shifted mean \eqn{-c(\bar{\theta})}, and that restriction to a box
#' set \eqn{A} yields closed-form expressions for both the normalizing
#' constant and the expectation.
#' @section Remarks on sampling from restricted normals:
#' \strong{Remark 7.} Sampling from the restricted normal densities in
#' Example 2 can be implemented using the inverse-transform method
#' (see, e.g., Fishman 1999).
#'
#' \strong{Remark 8.} In many applications, the inverse-transform method
#' of sampling from the restricted density in Example 2 will require
#' evaluating the cumulative normal distribution function (or its logarithm)
#' in the extreme tail of a normal distribution. Accurate computation in
#' this regime requires numerical procedures with uniformly small relative
#' errors. Authors presenting such procedures include Hart (1957, 1966)
#' and Bryc (2002).
#' 
#' 
#' @section Theorem 2 (log-concave univariate models with normal priors):
#' In the univariate case with a normal prior (variance = 1) and normal data,
#' the partitioning approach requires careful positioning of the intervals and
#' corresponding restricted likelihood-subgradient densities.
#'
#' Empirical investigation shows:
#' \itemize{
#'   \item A single optimally positioned likelihood-subgradient density deteriorates
#'         in performance as the number of data points increases.
#'   \item The same deterioration occurs for optimally positioned two-interval partitions.
#'   \item Remarkably, the optimal three-interval partition does not suffer this deterioration:
#'         the enveloping function remains a close approximation even as the sample size grows.
#' }
#'
#' Define the posterior mode \eqn{\theta^{\ast}} and set
#' \deqn{
#'   \omega :=
#'   \frac{\sqrt{2} - \exp\!\big(-1.20491 - 0.7321\,(0.5 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta^{2})\big)}
#'        {1 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta^{2}}.
#' }
#'
#' Then define the partition points
#' \deqn{
#'   \ell_{1} := \theta^{\ast} - 0.5\,\omega, \quad
#'   \ell_{2} := \theta^{\ast} + 0.5\,\omega,
#' }
#' and the three regions
#' \deqn{
#'   A_{1} := (-\infty,\ell_{1}), \quad
#'   A_{2} := [\ell_{1},\ell_{2}], \quad
#'   A_{3} := (\ell_{2},\infty).
#' }
#' The tangency points are chosen as
#' \deqn{
#'   \bar{\theta}_{1} = \theta^{\ast} - \omega, \quad
#'   \bar{\theta}_{2} = \theta^{\ast}, \quad
#'   \bar{\theta}_{3} = \theta^{\ast} + \omega.
#' }
#'
#' Note: if the data represent \eqn{N} observations from a normal density with unit variance,
#' then \eqn{-\partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta^{2} = N}.
#'
#' \strong{Theorem 2.} Consider this normal data model and let \eqn{\tilde{a}^{\ast}(N)}
#' denote the value of \eqn{\tilde{a}} at sample size \eqn{N}. Then
#' \deqn{
#'   \lim_{N \to \infty} \tilde{a}^{\ast}(N) = \frac{2}{\sqrt{\pi}}.
#' }
#'
#' \emph{Proof sketch.} By symmetry, \eqn{\tilde{a}_{1}(\bar{\theta}_{1}) = \tilde{a}_{3}(\bar{\theta}_{3})}.
#' It follows that
#' \deqn{
#'   \lim_{N \to \infty} \tilde{a}^{\ast}(N)
#'   = \lim_{N \to \infty} \big(\tilde{a}_{2}(\bar{\theta}_{2}) + 2\,\tilde{a}_{3}(\bar{\theta}_{3})\big)
#'   = \frac{1}{\sqrt{\pi}} + \frac{2}{2\sqrt{\pi}}
#'   = \frac{2}{\sqrt{\pi}},
#' }
#' where the second-to-last equality follows from Claims A.1 and A.4 in the Appendix.
#' @section Log-concave models with multivariate normal priors:
#' To ensure that the quality of the enveloping function remains good in the
#' multivariate case, models are first reparameterized into a standard form.
#'
#' \strong{Definition 3.} A probability model with a multivariate normal prior
#' and log-concave likelihood function is in \emph{standard form} if:
#' \itemize{
#'   \item (a) the prior variance-covariance matrix is the identity matrix, and
#'   \item (b) the Hessian of the log-posterior density evaluated at the posterior mode
#'         is a diagonal matrix.
#' }
#'
#' \emph{Remarks on reparameterization:}
#' \itemize{
#'   \item \strong{Remark 11.} If the log-likelihood function is concave and twice
#'         continuously differentiable, then a Cholesky decomposition of the posterior
#'         precision at the unique posterior mode can be used to reparameterize the model
#'         so that the posterior precision at the mode is diagonal.
#'   \item \strong{Remark 12.} For any positive definite matrix \eqn{P}, there exists a
#'         positive definite diagonal matrix \eqn{D} such that \eqn{P - D} is also positive definite.
#'   \item \strong{Remark 13.} Let \eqn{P} and \eqn{D} be as above. Then the following two
#'         models have the same posterior density:
#'         \enumerate{
#'           \item A model with prior mean vector \eqn{\mu}, prior precision matrix \eqn{P},
#'                 and log-likelihood function \eqn{LL(\theta)}.
#'           \item A model with prior mean vector \eqn{\mu}, prior precision matrix \eqn{D},
#'                 and log-likelihood function
#'                 \deqn{
#'                   LL^{\ast}(\theta) = -\tfrac{1}{2}(\theta - \mu)^{T}(P - D)(\theta - \mu) + LL(\theta).
#'                 }
#'         }
#'   \item \strong{Remark 14.} If a model has a multivariate normal prior with a diagonal
#'         variance-covariance matrix and the posterior precision at the posterior mode is
#'         also diagonal, then the model can be reparameterized into standard form.
#'   \item \strong{Remark 15.} If a probability model with a multivariate normal prior and
#'         log-concave likelihood function has a twice continuously differentiable log-likelihood,
#'         then it can be reparameterized into standard form.
#' }
#'
#' \emph{Multivariate partition construction:}
#'
#' Let \eqn{\theta^{\ast}} denote the unique posterior mode. For each dimension \eqn{i},
#' define
#' \deqn{
#'   \omega_{i} :=
#'   \frac{\sqrt{2} - \exp\!\big(-1.20491 - 0.7321\,(0.5 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2})\big)}
#'        {1 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2}}.
#' }
#'
#' Then set
#' \deqn{
#'   \ell_{i,1} = \theta^{\ast}_{i} - 0.5\,\omega_{i}, \quad
#'   \ell_{i,2} = \theta^{\ast}_{i} + 0.5\,\omega_{i},
#' }
#' and construct three intervals per dimension:
#' \deqn{
#'   A_{i,1} = (-\infty,\ell_{i,1}), \quad
#'   A_{i,2} = [\ell_{i,1},\ell_{i,2}], \quad
#'   A_{i,3} = (\ell_{i,2},\infty).
#' }
#'
#' For each dimension \eqn{i}, let \eqn{J_{i} = \{1,2,3\}} and define
#' \eqn{J = \prod_{i=1}^{p} J_{i}}, which has \eqn{3^{p}} elements. Each
#' \eqn{j \in J} is a vector \eqn{(j_{1},\ldots,j_{p})}, and we define
#' \deqn{
#'   A^{\ast}_{j} = \prod_{i=1}^{p} A_{i,j_{i}}.
#' }
#' The collection \eqn{A^{\ast} = \{A^{\ast}_{j} : j \in J\}} forms a partition of \eqn{\Theta}.
#'
#' For each \eqn{j \in J}, define index sets
#' \deqn{
#'   C_{j1} = \{i : j_{i} = 1\}, \quad
#'   C_{j2} = \{i : j_{i} = 2\}, \quad
#'   C_{j3} = \{i : j_{i} = 3\}.
#' }
#'
#' The tangency points \eqn{\bar{\theta}_{j}} are then defined componentwise by
#' \deqn{
#'   \bar{\theta}_{j,i} =
#'   \begin{cases}
#'     \theta^{\ast}_{i} - \omega_{i}, & i \in C_{j1}, \\
#'     \theta^{\ast}_{i},              & i \in C_{j2}, \\
#'     \theta^{\ast}_{i} + \omega_{i}, & i \in C_{j3}.
#'   \end{cases}
#' }
#'
#' \strong{Remark 16.} The mixture-generalized likelihood-subgradient density
#' that results from this construction is a mixture of restricted multivariate
#' normal densities, for which a straightforward sampling procedure exists.
#' 
#' 
#' @section Subgradient density formulation:
#' Each grid component corresponds to a tilted multivariate normal density,
#' normalized using the moment-generating function (MGF). In the single-point
#' case, centered at the posterior mode \eqn{\theta^\star}, the density is:
#' \deqn{
#' f(\theta) = \frac{1}{(2\pi)^{p/2} |A|^{-1/2} \cdot \text{MGF}_A(c)} \exp\left( -\frac{1}{2} (\theta - \mu)^T A (\theta - \mu) + c^T (\theta - \theta^\star) \right)
#' }
#' where:
#' - \eqn{A} is the precision matrix,
#' - \eqn{\mu} is the prior mean vector,
#' - \eqn{c} is the gradient of the log-likelihood at \eqn{\theta^\star},
#' - \eqn{\text{MGF}_A(c)} is the moment-generating function:
#' \deqn{
#' \text{MGF}_A(c) = \exp\left( \frac{1}{2} c^T A^{-1} c \right)
#' }
#'
#' This closed-form density dominates the posterior locally and is used when
#' \code{Gridtype = 1}. For richer envelopes, multiple such components are
#' constructed at tangency points \eqn{\theta_j}, each with its own gradient
#' \eqn{c_j}, and combined into a mixture:
#' \deqn{
#' f_{\text{env}}(\theta) = \sum_{j=1}^{K} p_j f_j(\theta)
#' }
#' where the weights \eqn{p_j} are computed using log-CDF differences and constants:
#' \deqn{
#' \log p_j = \log \Phi(U_j) - \log \Phi(L_j) - \text{NegLL}_j + \text{LLconst}_j
#' }
#' @section Algorithmic steps:
#' 1. Compute width parameters \eqn{\omega_i} from the diagonal of the precision matrix.
#' 2. Construct intervals around the posterior mode \eqn{\theta^\star}.
#' 3. Select tangency points at the mode and at \eqn{\theta^\star \pm \omega_i}.
#' 4. Build the full grid of tangency points (Cartesian product across dimensions).
#' 5. Evaluate negative log-likelihood and gradients at each grid point:
#'    - On CPU: via \code{f2_*} and \code{f3_*} routines.
#'    - On GPU: via \code{f2_f3_opencl}, which computes both in parallel.
#' 6. Call \code{EnvelopeSet_Grid_C2_pointwise} to evaluate restricted multivariate normal
#'    log-densities in parallel.
#' 7. Call \code{EnvelopeSetlogP_C2} to compute component log-probabilities and constants.
#' 8. Normalize probabilities (\code{PLSD}) and optionally sort grid components
#'    by probability if \code{sortgrid = TRUE}.
#' @section Gridtype logic:
#' The \code{Gridtype} argument controls how many tangency points are used per dimension:
#' - 1: Threshold rule. If \eqn{1 + a_i \le 2/\sqrt{\pi}}, use a single-point envelope at the mode;
#'      otherwise use three points.
#' - 2: Dynamic optimization via \code{EnvelopeOpt}, which balances grid build cost and
#'      expected acceptance rate. Grid size is scaled by \code{n} and the number of
#'      OpenCL cores when GPU is enabled.
#' - 3: Always use three points per dimension.
#' - 4: Always use a single point (mode only).
#' @section Supported families and links:
#' The following families and link functions are supported:
#' - Binomial: logit, probit, cloglog
#' - Quasibinomial: logit, probit
#' - Poisson: log
#' - Quasipoisson: log
#' - Gamma: log
#' - Gaussian: identity
#'
#' GPU acceleration (\code{use_opencl = TRUE}) is available for all of the above
#' except Gaussian, which is always evaluated on CPU.
#' @section GPU acceleration:
#' When \code{use_opencl = TRUE}, likelihood and gradient evaluations are
#' offloaded to the GPU using OpenCL. This can substantially reduce runtime for
#' high-dimensional models or large grids. Results are mathematically equivalent
#' to the CPU version, but small numerical differences may occur due to
#' floating-point arithmetic. If reproducibility across hardware is critical,
#' prefer the CPU path.
#'
#' If OpenCL support was not detected at compile time, the flag is ignored and
#' the CPU implementation is used. Diagnostic messages are printed when
#' \code{verbose = TRUE}.
#' @section Verbose output:
#' When \code{verbose = TRUE}, the function prints:
#' - Grid type, number of draws, OpenCL usage, and detected core count.
#' - Grid size after expansion.
#' - Time-stamped messages when entering the grid loop, starting likelihood
#'   evaluations, starting gradient evaluations, and invoking GPU kernels.
#' - Messages when setting grid values, computing log-probabilities, and sorting.
#' 
#' Any constants needed by the sampling are added to a list and returned.
#'
#' @param bStar     Point at which envelope should be centered (typically posterior mode).
#' @param A         Diagonal precision matrix for the log-likelihood in standard form.
#' @param y         A vector of observations of length \code{m}.
#' @param x A design matrix of dimension \code{m * p}.
#' @param mu        A vector giving the prior means of the variables.
#' @param P         Prior precision matrix of the variables (positive-definite).
#' @param alpha     Offset vector.
#' @param wt        A vector of weights.
#' @param family    Family for the envelope: \code{binomial}, \code{quasibinomial}, \code{poisson}, \code{quasipoisson}, or \code{Gamma}.
#' @param link      Link function ("logit", "probit", "cloglog" for binomial; "log" for Poisson/Gamma).
#' @param Gridtype  Method to determine the number of subgradient densities in the grid.
#' @param n         Number of draws from the posterior (used for grid sizing).
#' @param n_envopt Effective sample size passed to EnvelopeOpt for grid construction.
#'   Defaults to match `n`. Larger values encourage tighter envelopes.
#' @param sortgrid  Logical; if \code{TRUE}, sort the envelope descending by component probability.
#' @param use_opencl Logical; if \code{TRUE}, use OpenCL for gradient evaluations.
#' @param verbose   Logical; if \code{TRUE}, print progress messages.
#'
#' @param GridIndex A matrix indicating, for each grid component, whether the component
#'   lies in the left tail, center, or right tail of the density. Rows correspond to
#'   grid components; columns correspond to standardized variables.
#' @param cbars     A matrix containing the subgradient of the (adjusted) negative log-likelihood
#'   at each grid component.
#' @param Lint      A matrix storing the lower and upper bounds for each grid component,
#'   depending on whether sampling is from the left, center, or right.
#'
#' @param logP      A matrix (typically two columns) with information for each grid component.
#'   The first column usually holds the output from \code{EnvelopeSetGrid()}, corresponding to
#'   the restricted normal density.
#' @param NegLL     A vector of negative log-likelihood evaluations at each grid component.
#' @param G3        A matrix of tangency points used in the grid.
#'    
#' @return
#' \describe{
#'
#'   \item{\code{EnvelopeBuild()}}{A list of envelope components used for accept-reject sampling:
#'     \describe{
#'       \item{\code{GridIndex}}{Integer matrix encoding sampling type (tail, center, line) per dimension and region.}
#'       \item{\code{thetabars}}{Matrix of tangency points \eqn{\bar{\theta}_j} for each grid region.}
#'       \item{\code{cbars}}{Matrix of subgradients \eqn{c(\bar{\theta}_j)} of the negative log-likelihood at tangency.}
#'       \item{\code{loglt}}{Matrix of log left-tail probabilities per dimension and region.}
#'       \item{\code{logrt}}{Matrix of log right-tail probabilities per dimension and region.}
#'       \item{\code{logU}}{Matrix of selected per-dimension log-density contributions (tail/center) for each region.}
#'       \item{\code{logP}}{Matrix of total log-probabilities per region (first column); used to derive mixture weights.}
#'       \item{\code{PLSD}}{Vector of normalized mixture weights over grid regions used to draw region indices.}
#'       \item{\code{LLconst}}{Vector of acceptance-test constants per region used in the inequality for rejection sampling.}
#'     }
#'   }
#'
#'   \item{\code{EnvelopeSetGrid()}}{A list of matrices computed for grid-based log-density evaluation:
#'     \describe{
#'       \item{\code{Down}}{Lower bounds for truncated-normal evaluation per dimension and region.}
#'       \item{\code{Up}}{Upper bounds for truncated-normal evaluation per dimension and region.}
#'       \item{\code{lglt}}{Log left-tail probabilities (from \eqn{(-\infty, \mathrm{Up}]}) per dimension and region.}
#'       \item{\code{lgrt}}{Log right-tail probabilities (from \eqn{[\mathrm{Down}, \infty)}) per dimension and region.}
#'       \item{\code{lgct}}{Log central-interval probabilities (from \eqn{[\mathrm{Down}, \mathrm{Up}]}) per dimension and region.}
#'       \item{\code{logU}}{Selected log-probability per grid cell based on \code{GridIndex} (tail or center).}
#'       \item{\code{logP}}{Matrix with row-wise sums of \code{logU} (first column) used to form mixture weights.}
#'     }
#'   }
#'
#'   \item{\code{EnvelopeSetLogP()}}{A list with updated mixture-weight and acceptance constants:
#'     \describe{
#'       \item{\code{logP}}{Input \code{logP} with its second column populated by the log of unnormalized visit probabilities per region (mixture denominators).}
#'       \item{\code{LLconst}}{Vector of acceptance constants \eqn{-\log f(y \mid \bar{\theta}_j) - c(\bar{\theta}_j)^{T}\bar{\theta}_j} used in the accept-reject test.}
#'     }
#'   }
#' }
#' 
#'  @references
#' \insertAllCited{}
#' @example inst/examples/Ex_EnvelopeBuild.R
#' @importFrom Rdpack reprompt


#' @usage EnvelopeBuild(bStar,A,y,x,mu,P,alpha,wt,family = "binomial",link = "logit", 
#' Gridtype = 2L,n = 1L,n_envopt=NULL,sortgrid = FALSE,use_opencl = FALSE,verbose = FALSE)
#' @rdname EnvelopeBuild
#' @export
EnvelopeBuild <- function(
    bStar, A, y, x, mu, P, alpha, wt,
    family     = "binomial",
    link       = "logit",
    Gridtype   = 2L,
    n          = 1L,
    n_envopt   = NULL,       # effective sample size for EnvelopeOpt
    sortgrid   = FALSE,
    use_opencl = FALSE,
    verbose    = FALSE
) {
  # normalize n_envopt: if not supplied, fall back to n
  if (is.null(n_envopt)) {
    n_envopt <- n
  }
  
  if (length(n_envopt) != 1L || is.na(n_envopt) || n_envopt < 1) {
    stop("`n_envopt` must be a positive integer scalar.")
  }
  # coerce safely to integer
  n_envopt <- as.integer(n_envopt)
  
  if (family == "gaussian") {
    return(.EnvelopeBuild_Ind_Normal_Gamma_cpp(
      bStar, A, y, x, mu, P, alpha, wt,
      family = family, link = link,
      Gridtype = Gridtype, n = n,    
      n_envopt  = n_envopt, 
      sortgrid = sortgrid,
      use_opencl = use_opencl,
      verbose    = verbose
    ))
  }
  
  .EnvelopeBuild_cpp(
    bStar, A, y, x, mu, P, alpha, wt,
    family    = family,
    link      = link,
    Gridtype  = Gridtype,
    n         = n,
    n_envopt  = n_envopt,
    sortgrid  = sortgrid,
    use_opencl = use_opencl,
    verbose    = verbose
  )
}



#' @usage EnvelopeSetGrid(GridIndex, cbars, Lint)
#' @rdname EnvelopeBuild
#' @export
EnvelopeSetGrid <- function(GridIndex, cbars, Lint) {
  .EnvelopeSet_Grid_cpp(GridIndex, cbars, Lint)
}



#' @usage EnvelopeSetLogP(logP, NegLL, cbars, G3)
#' @rdname EnvelopeBuild
#' @export
EnvelopeSetLogP <- function(logP, NegLL, cbars, G3) {
  .EnvelopeSet_LogP_cpp(logP, NegLL, cbars, G3)
}



#' Evaluate Negative Log-Likelihood and Gradients
#'
#' \code{EnvelopeEval()} evaluates the negative log-likelihood and gradients
#' at a grid of parameter values, optionally using OpenCL acceleration.
#'
#' The lower-level helpers `f2_f3_non_opencl` and `f2_f3_opencl`
#' are internal C++ kernels used by the CPU and OpenCL backends.
#' The internal routine `run_opencl_pilot` benchmarks OpenCL performance
#' on a pilot subset of the grid to estimate runtime before full evaluation.
#'
#' These functions implement the grid evaluation logic used in envelope
#' construction for rejection sampling. They make use of the theory described
#' in \insertCite{Nygren2006}{glmbayes} and the general implementation outlined
#' in \insertCite{glmbayesSimmethods}{glmbayes}.
#'
#' @param G4 Numeric matrix of parameter values (parameters * grid points).
#' @param y Numeric response vector.
#' @param x Numeric design matrix.
#' @param mu Numeric matrix of offsets or prior means.
#' @param P Numeric matrix representing the portion of the prior precision
#'   shifted into the likelihood.
#' @param alpha Numeric vector of prior shape parameters.
#' @param wt Numeric vector of weights.
#' @param family Character string; model family (e.g. \code{"gaussian"}).
#' @param link Character string; link function (e.g. \code{"identity"}).
#' @param use_opencl Logical; if \code{TRUE}, attempt OpenCL acceleration.
#' @param verbose Logical; if \code{TRUE}, print diagnostic output.
#' @details
#' The evaluation workflow has several layers:
#' **1. High-level dispatch (`EnvelopeEval`)**
#'
#' * `EnvelopeEval()` is the user-facing entry point. It accepts a grid of
#'   parameter values (`G4`) and the data (`y`, `x`, `mu`, `P`, `alpha`, `wt`).
#' * If the grid is large (>= 14 columns), it first calls
#'   `run_opencl_pilot` to benchmark OpenCL performance and optionally
#'      report estimated runtime.
#' * It then dispatches to either the CPU or GPU backend:
#'   - If `use_opencl = TRUE` and the family is not `"gaussian"`, it calls
#'     `f2_f3_opencl` (an internal C++ kernel).
#'   - Otherwise, it calls `f2_f3_non_opencl` (the CPU kernel).
#'   
#' **2. CPU backend (`f2_f3_non_opencl`)**
#'
#' * This function evaluates the negative log-likelihood and gradients using
#'   standard CPU routines.
#' * It inspects the `family` and `link` arguments and routes to the correct
#'   pair of kernels (`f2_*` for the likelihood, `f3_*` for the gradient).
#' * For example:
#'   - `"binomial"` with `"logit"` calls `f2_binomial_logit()` and
#'     `f3_binomial_logit()`.
#'   - `"poisson"` calls `f2_poisson()` and `f3_poisson()`.
#'   - `"gaussian"` calls `f2_gaussian()` and `f3_gaussian()`.
#' * These kernels ultimately rely on the same C math routines that R itself
#'   uses (from the `nmath`/`rmath` libraries), ensuring numerical consistency
#'   with base R functions like `dnorm`, `dpois`, etc.
#'
#' **3. GPU backend (`f2_f3_opencl`)**
#'
#' * This function mirrors the CPU backend but executes the likelihood and
#'   gradient calculations on an OpenCL device (GPU or CPU).
#' * It flattens the input matrices/vectors and allocates output buffers.
#' * It then constructs a full OpenCL program by concatenating:
#'   - a generic OpenCL support header (`OPENCL.CL`),
#'   - OpenCL ports of R's `rmath`, `nmath`, and `dpq` libraries,
#'   - and the family/link-specific kernel source (e.g.
#'     `f2_f3_binomial_logit.cl`).
#' * The resulting program is compiled and passed to a kernel runner
#'   (`f2_f3_kernel_runner`) which executes the likelihood and gradient
#'   calculations in parallel on the device.
#' * This ensures that the GPU backend produces results consistent with the
#'   CPU backend, but can scale to much larger grids efficiently.
#'
#' **4. Pilot timing (`run_opencl_pilot`)**
#'
#' * This helper runs a small subset of the grid through the OpenCL backend
#'   to estimate runtime.
#' * It is used by `EnvelopeEval()` to inform users (when `verbose = TRUE`)
#'   whether OpenCL acceleration is likely to be beneficial.
#'
#' **5. Returned values**
#'
#' * All backends return a list with:
#'   - `NegLL`: numeric vector of negative log-likelihood values.
#'   - `cbars`: numeric matrix of gradients (parameters * grid points).
#'   
#' **6. Role of likelihood and gradients in sampling**
#'
#' * The outputs of `EnvelopeEval()` - the negative log-likelihood values
#'   (`NegLL`) and the gradient matrix (`cbars`) - are not endpoints in
#'   themselves. They form the *envelope* used in the rejection sampler
#'   implemented by internal functions such as
#'   `rnnorm_reg_std_cpp()` and `rnnorm_reg_std_cpp_parallel()`.
#' * These routines are called by `rnnorm_reg_cpp()`, which underlies the
#'   user-facing function `rNormal_reg()`. Together they implement
#'   envelope-based posterior sampling for GLMs with log-concave likelihoods
#'   and multivariate normal priors.
#'
#' **7. Simulation execution (accept/reject procedure)**
#'
#' The acceptance test is performed using
#' \deqn{
#'   \log(U_2) \leq
#'   \log f(y \mid \theta_i) -
#'   \Big(\log f(y \mid \bar{\theta}_{J(i)}) -
#'        c(\bar{\theta}_{J(i)})^T(\theta_i - \bar{\theta}_{J(i)})\Big) \leq 0
#' }
#'
#' **Connections between code and notation:**
#' * The arguments `G4` (in `EnvelopeEval`) and `b` (in `f2_f3_*`) both
#'   represent the grid of tangency points \eqn{\bar{\theta}_j}.
#' * The output `NegLL` corresponds to
#'   \eqn{-\log f(y \mid \bar{\theta}_{J(i)})}, i.e. the negative
#'   log-likelihood evaluated at each tangency point.
#' * The output `cbars` corresponds to the subgradient vectors
#'   \eqn{c(\bar{\theta}_{J(i)})}, which define the tangent hyperplanes
#'   used in the envelope construction.
#'
#' **Precomputation for efficiency:**
#' * Both `NegLL` and `cbars` are computed once during envelope construction,
#'   prior to the simulation stage.
#' * This means the sampler does not need to recompute likelihoods or
#'   gradients at every candidate draw - it simply reuses the stored values
#'   (`NegLL`, `cbars`, and `LLconst`) in the acceptance inequality.
#'
#' This design ensures that the envelope is tangent to the log-likelihood at
#' each \eqn{\bar{\theta}_j}, lies above it elsewhere, and that the
#' accept-reject procedure can run efficiently while still producing samples
#' from the true posterior \eqn{\pi(\theta \mid y)}.#'
#'
#' @references
#' \insertAllCited{}
#' @return
#' \describe{
#'   \item{EnvelopeEval}{List with components \code{NegLL} (numeric vector of
#'   negative log-likelihood values) and \code{cbars} (numeric matrix of gradients).}
#'   \item{f2_f3_non_opencl}{List with components \code{qf} (negative log-likelihood)
#'   and \code{grad} (gradients) from the CPU kernel.}
#'   \item{f2_f3_opencl}{List with components \code{qf} and \code{grad} from the
#'   OpenCL kernel.}
#'   \item{run_opencl_pilot}{Numeric scalar giving estimated runtime (seconds)
#'   for OpenCL evaluation on a pilot subset of the grid.}
#' }
#'
#' @rdname EnvelopeEval
#' @export
#' @usage EnvelopeEval(G4, y, x, mu, P, alpha, wt,
#'                     family, link,
#'                     use_opencl = FALSE, verbose = FALSE)
EnvelopeEval <- function(G4, y, x, mu, P, alpha, wt,
                         family, link,
                         use_opencl = FALSE,
                         verbose = FALSE) {
  if (!is.matrix(G4)) stop("G4 must be a numeric matrix")
  if (!is.numeric(y)) stop("y must be numeric")
  if (!is.matrix(x)) stop("x must be a numeric matrix")
  if (!is.matrix(mu)) stop("mu must be a numeric matrix")
  if (!is.matrix(P)) stop("P must be a numeric matrix")
  if (!is.numeric(alpha)) stop("alpha must be numeric")
  if (!is.numeric(wt)) stop("wt must be numeric")
  if (!is.character(family) || length(family) != 1L) stop("family must be a string")
  if (!is.character(link) || length(link) != 1L) stop("link must be a string")
  
  
  .EnvelopeEval_cpp(G4, y, x, mu, P, alpha, wt,
               family, link,
               use_opencl, verbose)
}


#' Builds Dispersion-Aware Envelope for Simulation
#'
#' @description
#' Constructs a dispersion-aware envelope for simulation in Gaussian models with uncertain variance.
#' This function extrapolates the coefficient envelope across a high-probability interval for the
#' dispersion parameter \code{sigma^2}, and builds a global upper bound for the log-posterior remainder.
#' It also computes mixture weights for envelope faces and adjusts the Gamma proposal for precision.
#'
#' The envelope is constructed using the slopes of the face constants with respect to dispersion,
#' evaluated at an anchor point. The resulting structure supports exact i.i.d. sampling via accept-reject correction.
#'
#' The procedure follows these steps:
#' \enumerate{
#'
#'   \item \strong{Posterior precision (Gamma) parameters.}
#'   Using the prior and posterior-predictive RSS, set
#'
#'   \deqn{ \mathrm{shape2} = \mathrm{Shape} + n_{\mathrm{obs}}/2, \quad
#'          \mathrm{rate3}  = \mathrm{Rate} + \mathrm{RSS}_{\mathrm{post}}/2 }
#'
#'   These parameterize the posterior precision
#'   \eqn{v = 1/\sigma^2 \sim \mathrm{Gamma}(\mathrm{shape2}, \mathrm{rate3})}.
#'
#'   \item \strong{Central credible interval for dispersion (low, upp).}
#'   Choose a central mass level \code{max_disp_perc} (e.g., 0.99) for precision, then invert the corresponding
#'   Gamma quantiles to dispersion:
#'
#'   \deqn{ \mathrm{low} = 1 / Q_{\Gamma}(\mathrm{max\_disp\_perc}; \mathrm{shape2}, \mathrm{rate3}), \quad
#'          \mathrm{upp} = 1 / Q_{\Gamma}(1 - \mathrm{max\_disp\_perc}; \mathrm{shape2}, \mathrm{rate3}) }
#'
#'   The interval \eqn{[\mathrm{low}, \mathrm{upp}]} is the domain over which all envelopes must dominate.
#'
#'   \item \strong{Face slopes at an anchor (dispstar).}
#'
#'   \deqn{ \mathrm{dispstar} = \mathrm{rate3} / (\mathrm{shape2} - 1) }
#'
#'   (posterior mean of \eqn{\sigma^2}). Compute \eqn{\mathrm{New\_LL\_Slope}_j} for each face \eqn{j}.
#'
#'   \item \strong{Linear extrapolation of face constants.}
#'
#'   \deqn{ \theta^{\mathrm{low}}_{\bar{j}} = \theta^{\mathrm{base}}_{\bar{j}} + (\mathrm{low} - \mathrm{dispstar}) \cdot \mathrm{New\_LL\_Slope}_j }
#'
#'   \deqn{ \theta^{\mathrm{upp}}_{\bar{j}} = \theta^{\mathrm{base}}_{\bar{j}} + (\mathrm{upp} - \mathrm{dispstar}) \cdot \mathrm{New\_LL\_Slope}_j }
#'
#'   \item \strong{Global upper line and endpoint maxima.}
#'
#'   \deqn{ \mathrm{max\_low} = \max_j \theta^{\mathrm{low}}_{\bar{j}}, \quad
#'          \mathrm{max\_upp} = \max_j \theta^{\mathrm{upp}}_{\bar{j}} }
#'
#'   \deqn{ \mathrm{new\_slope} = (\mathrm{max\_upp} - \mathrm{max\_low}) / (\mathrm{upp} - \mathrm{low}), \quad
#'          \mathrm{new\_int}   = \mathrm{max\_low} - \mathrm{new\_slope} \cdot \mathrm{low} }
#'
#'   \item \strong{Face slack and mixture weights.}
#'
#'   \deqn{ \mathrm{lg\_prob\_factor}_j =
#'          \max\big(\theta^{\mathrm{upp}}_{\bar{j}} - \mathrm{max\_upp},\;
#'                    \theta^{\mathrm{low}}_{\bar{j}} - \mathrm{max\_low}\big) }
#'
#'   Combine with \eqn{\mathrm{New\_logP2}_j = \mathrm{logP}_j + \tfrac{1}{2}\|\bar{c}_j\|^2}
#'   to form mixture weights
#'   \eqn{\mathrm{PLSD}_j \propto \exp\!\bigl(\mathrm{New\_logP2}_j + \mathrm{lg\_prob\_factor}_j\bigr)}.
#'
#'   \item \strong{Gamma tilt and dispersion-axis envelope.}
#'
#'   \deqn{ \mathrm{dispstar} = (\mathrm{upp} - \mathrm{low}) / \log(\mathrm{upp}/\mathrm{low}) }
#'
#'   \deqn{ \mathrm{lm\_log2} = \mathrm{new\_slope} \cdot \mathrm{dispstar}, \quad
#'          \mathrm{lm\_log1} = \mathrm{new\_int} + \mathrm{new\_slope} \cdot \mathrm{dispstar}
#'                             - \mathrm{new\_slope} \cdot \log(\mathrm{dispstar}) }
#'
#'   Tilt the Gamma proposal via \eqn{\mathrm{shape3} = \mathrm{shape2} - \mathrm{lm\_log2}}.
#' }
#'
#' @param Env        Envelope object from \code{\link{EnvelopeBuild}}, containing tangency points and gradients
#' @param Shape      Prior shape parameter for precision \code{v = 1 / sigma^2}
#' @param Rate       Prior rate parameter for precision
#' @param mu         Prior mean parameter
#' @param P          Prior precision matrix for coefficients
#' @param y          Design matrix
#' @param x          Design matrix
#' @param alpha          Design matrix
#' @param n_obs      Number of observations
#' @param RSS_post   Posterior-predictive residual sum of squares
#' @param RSS_ML   Residual sum of squares associated with MLE estimate
#' @param max_disp_perc Truncation level for dispersion (default 0.99)
#' @param disp_lower lower bound truncation for dispersion 
#' @param disp_upper upper bound truncation for dispersion
#' @param verbose Option to have verbose output
#' @param wt          weight vector
#' @param use_parallel Logical. Whether to use parallel processing.
#'
#' @return
#' \describe{
#' \item{\code{EnvelopeDispersionBuild()}}{A list containing:
#' \describe{
#' \item{\code{Env_out}}{Envelope object with updated mixture weights (\code{PLSD})}
#' \item{\code{gamma_list}}{Posterior Gamma tilt parameters
#' \describe{
#' \item{\code{shape3}}{Adjusted shape parameter after slope correction}
#' \item{\code{rate2}}{Posterior rate parameter, defined as \code{Rate + rss_min_global/2}}
#' \item{\code{disp_upper}}{Upper bound of the dispersion interval \eqn{\sigma^2}}
#' \item{\code{disp_lower}}{Lower bound of the dispersion interval \eqn{\sigma^2}}
#' }}
#' \item{\code{UB_list}}{Upper-bound diagnostics
#' \describe{
#' \item{\code{RSS_ML}}{Residual sum of squares at the maximum-likelihood estimate}
#' \item{\code{RSS_Min}}{Minimum residual sum of squares across envelope faces}
#' \item{\code{max_New_LL_UB}}{Maximum extrapolated face constant at the upper dispersion bound}
#' \item{\code{max_LL_log_disp}}{Log-posterior upper bound evaluated at \code{disp_upper}}
#' \item{\code{lm_log1}}{Intercept term of the global upper line approximation}
#' \item{\code{lm_log2}}{Slope term of the global upper line approximation}
#' \item{\code{lg_prob_factor}}{Per-face slack factors used in mixture weighting}
#' \item{\code{lmc1}}{Linear extrapolation constant (intercept)}
#' \item{\code{lmc2}}{Linear extrapolation constant (slope)}
#' \item{\code{UB2min}}{Minimum UB2 value across faces, used for diagnostics}
#' }}
#' \item{\code{diagnostics}}{Internal diagnostic values
#' \describe{
#' \item{\code{dispstar}}{Anchor dispersion value (posterior mean or geometric mean)}
#' \item{\code{New_LL_Slope}}{Vector of slopes of face constants at \code{dispstar}}
#' \item{\code{shape2}}{Posterior shape parameter before tilt correction}
#' \item{\code{rate3}}{Posterior rate parameter before tilt correction}
#' \item{\code{shape3}}{Adjusted shape parameter (same as in \code{gamma_list})}
#' \item{\code{max_low}}{Maximum extrapolated face constant at the lower dispersion bound}
#' \item{\code{max_upp}}{Maximum extrapolated face constant at the upper dispersion bound}
#' \item{\code{new_slope}}{Slope of the global upper line across dispersion bounds}
#' \item{\code{new_int}}{Intercept of the global upper line across dispersion bounds}
#' \item{\code{prob_factor}}{Normalized mixture weights across faces}
#' \item{\code{UB2min}}{Minimum UB2 diagnostic value (duplicated for consistency)}
#' }}
#' }}
#' \item{\code{EnvBuildLinBound()}}{Numeric vector of slopes of face constants with respect to dispersion, evaluated at the anchor \code{dispstar}}
#' \item{\code{thetabar_const()}}{Numeric vector of base face constants computed from tangency points and gradient vectors under prior precision \code{P}}
#' \item{\code{Inv_f3_with_disp()}}{Numeric matrix of inverse function evaluations at a given dispersion and face subset, returned by the C++ routine \code{_glmbayes_Inv_f3_with_disp}}
#' \item{\code{UB2()}}{Numeric scalar representing the UB2 upper-bound criterion for a given dispersion and face, defined as \eqn{(1/dispersion) * (RSS - rss\_min\_global)}}
#' \item{\code{rss_face_at_disp()}}{Numeric scalar giving the residual sum of squares for a specified face at a given dispersion, computed from cached matrices and the inverse function evaluation}
#' }
#'    
#'
#' @details
#' This function is designed to complement \code{\link{EnvelopeBuild}} for Gaussian models
#' with Normal-Gamma priors. It enables exact sampling of both coefficients and dispersion
#' by constructing a joint envelope that respects posterior curvature in both dimensions.
#'
#' The dispersion anchor point is chosen as the log-scale center of the credible interval,
#' and the Gamma proposal is tilted to match the envelope slope at this point.
#' @section Use in accept/reject procedure:
#'
#' The accept/reject sampler relies on a decomposition of the log-posterior into
#' a test statistic and several bounding terms. Each component is constructed so
#' that its sign is controlled, ensuring the validity of the accept/reject step.
#'
#' \describe{
#'
#'   \item{test1 (log-likelihood bound)}{
#'     Placeholder: explain how test1 is formed and why it is non-positive.
#'   }
#'
#'   \item{UB1 (if applicable)}{
#'     Placeholder: describe UB1's role and why it is non-negative.
#'   }
#'
#'   \item{UB2 (residual sum of squares bound)}{
#'     Placeholder: explain how UB2 is constructed from RSS differences and why it is non-negative.
#'   }
#'
#'   \item{UB3A (face-wise quadratic/linear envelope surplus)}{
#'     Placeholder: explain how lg_prob_factor, lmc1, and lmc2 are derived and why UB3A >= 0.
#'   }
#'
#'   \item{UB3B (dispersion-axis envelope surplus)}{
#'     Placeholder: explain how lm_log1, lm_log2, and max_New_LL_UB are used and why UB3B >= 0.
#'   }
#'
#' }
#'
#' Together, these components define
#' \deqn{ test = test1 - UB2 - UB3A - UB3B, }
#' with \eqn{test1 \le 0} and each UB term \eqn{\ge 0}, ensuring the accept/reject
#' procedure is valid and unbiased.
#' @seealso \code{\link{EnvelopeBuild}}, \code{\link{glmb}}, \code{\link{glmbfamfunc}}
#' @usage EnvelopeDispersionBuild(
#'   Env, Shape, Rate, P, y, x, alpha, n_obs, RSS_post, RSS_ML,
#'   mu, wt, max_disp_perc = 0.99,
#'   disp_lower = NULL, disp_upper = NULL,
#'   verbose = FALSE, use_parallel = TRUE
#' )
#' @export 
#' @order 1


EnvelopeDispersionBuild <- function(Env, Shape, Rate, P, y, x, alpha, n_obs, RSS_post, RSS_ML,
                                    mu, wt, max_disp_perc = 0.99, disp_lower = NULL, disp_upper = NULL,
                                    verbose = FALSE, use_parallel = TRUE) {
  .EnvelopeDispersionBuild_cpp(
        Env, Shape, Rate, P, y, x, alpha, n_obs,
        RSS_post, RSS_ML, mu, wt,
        max_disp_perc, disp_lower, disp_upper,
        verbose, use_parallel)
}


#' Sorts Envelope function for simulation
#'
#' Sorts Enveloping function for simulation.
#' how frequently each component of the resulting grid should be sampled
#' during simulation.
#' @param l1 dimension for model (number of independent variables in X matrix)
#' @param l2 dimension for Envelope (number of components)
#' @param GIndex matrix containing information on how each dimension should be sampled (1 means left tail of a restricted normal, 2 center, 3 right tail, and 4 the entire line)
#' @param G3 A matrix containing the points of tangencies associated with each component of the grid
#' @param cbars A matrix containing the gradients for the negative log-likelihood at each tangency
#' @param logU A matrix containing the log of the cummulative probability associated with each dimension
#' @param logrt A matrix containing the log of the probability associated with the right tail (i.e. that to the right of the lower bound)
#' @param loglt A matrix containing the log of the probability associated with the left tail (i.e., that to the left of the upper bound)
#' @param logP A matrix containing log-probabilities related to the components of the grid
#' @param LLconst A vector containing constant for each component of the grid used during the accept-reject procedure
#' @param PLSD A vector containing the probability of each component in the Grid
#' @param a1 A vector containing the diagonal of the standardized precision matrix
#' @param E_draws Bound on Expected number of candidates per accepted draw
#' @param lg_prob_factor vector of lg_prob_factors used for the Envelope connected to the independent normal gamma prior
#' @param UB2min Vector containing min for UB2 for each component (relevant for EnvelopeDispersionBuild)
#' @return The function returns a list consisting of the following components (the first six of which are matrics with number of rows equal to the number of components in the Grid and columns equal to the number of parameters):
#' \item{GridIndex}{A matrix containing information on how each dimension should be sampled (1 means left tail of a restricted normal, 2 center, 3 right tail, and 4 the entire line)}
#' \item{thetabars}{A matrix containing the points of tangencies associated with each component of the grid}
#' \item{cbars}{A matrix containing the gradients for the negative log-likelihood at each tangency}
#' \item{logU}{A matrix containing the log of the cummulative probability associated with each dimension}
#' \item{logrt}{A matrix containing the log of the probability associated with the right tail (i.e. that to the right of the lower bound)}
#' \item{loglt}{A matrix containing the log of the probability associated with the left tail (i.e., that to the left of the upper bound)}
#' \item{LLconst}{A vector containing constant for each component of the grid used during the accept-reject procedure}
#' \item{logP}{A matrix containing log-probabilities related to the components of the grid}
#' \item{PLSD}{A vector containing the probability of each component in the Grid}
#' \item{E_draws}{A containing a computed theoretical bound on the expected number of draws}
#' @details This function sorts the envelope in descending order based on the 
#' probability associated with each component in the Grid. Sorting helps 
#' speed up simulation once the envelope is constructed.
#' @example inst/examples/Ex_EnvelopeSort.R
#' @export

EnvelopeSort <- function(l1, l2,
                         GIndex, G3, cbars, logU, logrt, loglt,
                         logP, LLconst, PLSD, a1, E_draws,
                         lg_prob_factor = NULL,
                         UB2min=NULL,
                         thetabar_const_base = NULL,
                         New_LL_Slope=NULL,
                         shape3_face=NULL
                         ) {   # <-- new optional arg
  # Order indices by decreasing PLSD
  ord <- order(PLSD, decreasing = TRUE)
  sel <- ord[seq_len(l2)]  # top l2 rows
  
  
  # Reorder inputs once
  GIndex <- GIndex[sel, , drop = FALSE]
  G3     <- G3[sel, , drop = FALSE]
  cbars  <- cbars[sel, , drop = FALSE]
  logU   <- logU[sel]
  logrt  <- logrt[sel, , drop = FALSE]
  loglt  <- loglt[sel, , drop = FALSE]
  logP   <- logP[sel, 1, drop = TRUE]
  LLconst<- LLconst[sel]
  PLSD   <- PLSD[sel]
  
  if (!is.null(lg_prob_factor)) {
    stopifnot(length(lg_prob_factor) == l2)
    lg_prob_factor <- lg_prob_factor[sel]
  }
  
  if (!is.null(UB2min)) {
    stopifnot(length(UB2min) == l2)
    UB2min <- UB2min[sel]
  }
  

  if (!is.null(thetabar_const_base)) {
    stopifnot(length(thetabar_const_base) == l2)
    thetabar_const_base <- thetabar_const_base[sel]
  }
  
  if (!is.null(New_LL_Slope)) {
    stopifnot(length(New_LL_Slope) == l2)
    New_LL_Slope <- New_LL_Slope[sel]
  }
  
  if (!is.null(shape3_face)) {
    stopifnot(length(shape3_face) == l2)
    shape3_face <- shape3_face[sel]
  }
  
  
  
    
  # Build output list
  outlist <- list(
    GridIndex = GIndex,
    thetabars = G3,
    cbars     = if (l1 == 1) matrix(cbars, nrow = l2, ncol = l1) else cbars,
    logU      = logU,
    logrt     = if (l1 == 1) matrix(logrt, nrow = l2, ncol = l1) else logrt,
    loglt     = if (l1 == 1) matrix(loglt, nrow = l2, ncol = l1) else loglt,
    LLconst   = LLconst,
    logP      = logP,
    PLSD      = PLSD,
    a1        = a1,
    E_draws   = E_draws
  )
  if (!is.null(lg_prob_factor)) {
    outlist$lg_prob_factor <- lg_prob_factor
  }
  if (!is.null(UB2min)) {
    outlist$UB2min <- UB2min
  }
  
  ## Adding to enable face specific dispersion bounds
  if (!is.null(thetabar_const_base)) {
    outlist$thetabar_const_base <- thetabar_const_base
  }
  
  if (!is.null(New_LL_Slope)) {
    outlist$New_LL_Slope <- New_LL_Slope
  }
    
  if (!is.null(shape3_face)) {
    outlist$shape3_face <- shape3_face
  }
  
  outlist
}


#' The Bayesian Generalized Linear Model Distribution in Standard Form
#'
#' \code{rnnorm_reg_std} is used to generate iid samplers from Non-Gaussian Generalized
#' Linear Models in standard form. The function should only be called after standardization 
#' of a Generalized Linear Model.
#' @param n number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y a vector of observations of length \code{m}.
#' @param x a design matrix of dimension \code{m * p}.
#' @param mu a vector of length \code{p} giving the prior means of the variables in the design matrix.
#' @param P a positive-definite symmetric matrix of dimension \code{p * p} specifying the prior precision
#'  matrix of the variable.
#' @param alpha this can be used to specify an \emph{a priori} known component 
#' to be included in the linear predictor during fitting. This should be 
#' \code{NULL} or a numeric vector of length equal to the number of cases. 
#' One or more offset terms can be included in the formula instead or as well,
#' and if more than one is specified their sum is used. See 
#' \code{\link{model.offset}}.
#' @param wt an optional vector of \sQuote{prior weights} to be used in the fitting process. 
#' Should be NULL or a numeric vector.
#' @param f2 function used to calculate the negative of the log-posterior function
#' @param Envelope an object of type \code{glmbenvelope}.
#' @param family family used for simulation. Used that this is different from the 
#' family used in other functions.
#' @param link link function used for simulation.
#' @param progbar dummy for flagging if a progressbar should be produced during the call
#' @return A list consisting of the following:
#' \item{out}{A matrix with simulated draws from a model in standard form. Each row represents 
#' one draw from the density}
#' \item{draws}{A vector with the number of candidates required before 
#' acceptance for each draw}
#' @details This function uses the information contained in the constructed envelope list in order to sample 
#' from a model in standard form. The simulation proceeds as follows in order to generate each draw in the required
#' number of samples.
#' 
#' 1)  A random number between 0 and 1 is generated and is used together with the information in the PLSD vector 
#' (from the envelope) in order to identify the part of the grid from which a candidate is to be generated.
#'
#' 2) For the part of the grid selected, the dimensions are looped through and a candidate component for each dimension
#' is generated from a restricted normal using information from the Envelope (in particular, the values for logrt, loglt, 
#' and cbars corresponding to that the part of the grid selected and the dimension sampled)
#' 
#' 3) The log-likelihood for the standardized model is evaluated for the generated candidate (note that the 
#' log-likelihood here includes the portion of the prior that was shifted to the log-likelihood)
#' 
#'  4) An additional random number is generated and the log of this random number is compared to a 
#'  log-acceptance rate that is calculated based on the candidate and the LLconst component from the Envelope 
#'  component selected in order to determine if the candidate should be accepted or rejected 
#' 
#' 5) If the candidate was not accepted, the process above is repeated from step 1 until a candidate is accepted
#' 
#' 
#' 
#' @example inst/examples/Ex_rnnorm_reg_std.R
#' @export



rnnorm_reg_std<-function(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link, progbar = 1L){
  
  return(.rNormalGLM_std_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link, progbar = progbar))
  
}




# Helpers --------------------------------------------------------------------


#' @noRd

dpois2<-function(x,lambda,log=TRUE){
  
  test=max(abs(round(x)-x))
  
  if(test>0){
    warning("Non-Integer Values to Poisson Density - Switching to Gamma Function to Evaluate Factorial")
    return(-lambda+x*log(lambda)-log(gamma(x+1)))
    
  } 
  
  return(dpois(x,lambda,log=TRUE))
}


#' @noRd

EnvelopeDispersionBuild_parallel_internal <- function(par0, low, upp,
                                                      cache, cbars, y, x, alpha, wt,
                                                      cores = parallel::detectCores(logical = FALSE),
                                                      use_parallel = TRUE) {
  gs <- nrow(cbars)
  avail_cores <- parallel::detectCores(logical = FALSE)
  
  # Respect CRAN check limit
  if (Sys.getenv("_R_CHECK_LIMIT_CORES_", "") == "TRUE") {
    cores <- min(cores, 2L)
  } else {
    cores <- min(cores, avail_cores, gs)
  }
  
  message("[EnvelopeDispersionBuild_parallel_internal] Available cores: ", avail_cores,
          " | Using cores: ", if (use_parallel) cores else 1L)
  
  worker_fun <- function(j) {
    cbars_j <- cbars[j, ]
    optim(par0,
          fn     = .rss_face_at_disp_cpp,
          method = "L-BFGS-B",
          lower  = low,
          upper  = upp,
          cache  = cache,
          cbars_j = cbars_j,
          y      = y,
          x      = x,
          alpha  = alpha,
          wt     = wt)
  }
  
  if (use_parallel) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      results <- parallel::parLapply(cl, seq_len(gs), worker_fun)
    } else {
      results <- parallel::mclapply(seq_len(gs), worker_fun, mc.cores = cores)
    }
  } else {
    # Serial fallback: just loop over indices
    results <- lapply(seq_len(gs), worker_fun)
  }
  
  disp_min <- vapply(results, function(res) res$par, numeric(1))
  rss_min  <- vapply(results, function(res) res$value, numeric(1))
  
  list(disp_min = disp_min, rss_min = rss_min)
}

#' @noRd

EnvelopeDispersionBuild_parallel <- function(par0, low, upp,
                                             cache, cbars, y, x, alpha, wt,
                                             cores = parallel::detectCores(logical = FALSE),
                                             use_parallel = TRUE) {
  tryCatch(
    EnvelopeDispersionBuild_parallel_internal(par0, low, upp, cache, cbars, y, x, alpha, wt, cores,use_parallel),
    error = function(e) {
      message("[EnvelopeDispersionBuild_parallel] Failed gracefully: ", e$message)
      list(disp_min = NULL, rss_min = NULL, error = TRUE)
    }
  )
}



#' @noRd

EnvelopeUB2_parallel_internal <- function(par0, low, upp,
                                          cache, cbars, y, x, alpha, wt,
                                          rss_min_global,
                                          cores = parallel::detectCores(logical = FALSE),
                                          use_parallel = TRUE) {
  gs <- nrow(cbars)
  avail_cores <- parallel::detectCores(logical = FALSE)
  
  # Respect CRAN check limit
  if (Sys.getenv("_R_CHECK_LIMIT_CORES_", "") == "TRUE") {
    cores <- min(cores, 2L)
  } else {
    cores <- min(cores, avail_cores, gs)
  }
  
  message("[EnvelopeUB2_parallel_internal] Available cores: ", avail_cores,
          " | Using cores: ", if (use_parallel) cores else 1L)
  
  worker_fun <- function(j) {
    cbars_j <- cbars[j, ]
    optim(par0,
          fn     = .UB2_cpp,                 # UB2 objective function
          method = "L-BFGS-B",
          lower  = low,
          upper  = upp,
          cache  = cache,
          cbars_j = cbars_j,
          y      = y,
          x      = x,
          alpha  = alpha,
          wt     = wt,
          rss_min_global = rss_min_global)
  }
  
  if (use_parallel) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      results <- parallel::parLapply(cl, seq_len(gs), worker_fun)
    } else {
      results <- parallel::mclapply(seq_len(gs), worker_fun, mc.cores = cores)
    }
  } else {
    # Serial fallback: run optim sequentially
    results <- lapply(seq_len(gs), worker_fun)
  }
  
  disp_min_ub2 <- vapply(results, function(res) res$par, numeric(1))
  ub2_min      <- vapply(results, function(res) res$value, numeric(1))
  
  list(disp_min = disp_min_ub2, ub2_min = ub2_min)
}

#' @noRd

EnvelopeUB2_parallel <- function(par0, low, upp,
                                 cache, cbars, y, x, alpha, wt,
                                 rss_min_global,
                                 cores = parallel::detectCores(logical = FALSE),
                                 use_parallel = TRUE) {
  tryCatch(
    EnvelopeUB2_parallel_internal(par0, low, upp, cache, cbars, y, x, alpha, wt,
                                  rss_min_global, cores),
    error = function(e) {
      message("[EnvelopeUB2_parallel] Graceful failure: ", e$message)
      list(disp_min = NULL, ub2_min = NULL, error = TRUE)
    }
  )
}

