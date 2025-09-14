#' Prior Family Objects for Bayesian Models
#'
#' Prior family objects provide a convenient way to specify the details of the priors 
#' used by functions such as \code{\link{glmb}}. See the documentations for \code{\link{lmb}},
#' \code{\link{glmb}}, \code{\link{glmb}}, and \code{\link{rglmb}} for the details of how such model fitting 
#' takes place.
#' @param object the function \code{pfamily} accesses the \code{pfamily} objects which
#' are stored within objects created by modelling functions (e.g., \code{glmb}).
#' @param mu a prior mean vector for the the modeling coefficients used in several pfamilies
#' @param Sigma a prior Variance-Covariance matrix for the model coefficients in several pfamilies
#' @param dispersion the dispersion to be assumed when it is not given a prior. Should be provided
#' when the Normal prior is for the \code{gaussian()}, \code{Gamma()}, \code{quasibinomial},
#' or \code{quasipoisson} families. The \code{binomial()} and \code{poisson()} families
#' do not have dispersion coefficients. 
#' @param shape the prior shape parameter used by the gamma component of the prior. 
#' The gamma distribution is used as a prior for the inverse dispersion coefficients.
#' @param rate the rate parameter used by the gamma component of the prior.
#' @param beta the regression coefficients to be assumed when it is not given a prior. 
#' Needs to be provided when the Gamma prior is used for the dispersion. This
#' specification is typically only used as part of Gibbs sampling where the beta and 
#' dispersion parameters are updated separately. 
#' @param x an object, a pfamily function that is to be printed
#' @param \ldots additional argument(s) for methods.
#' @details \code{pfamily} is a generic function with methods for classe \code{glmb} and 
#' \code{lmb}. Many \code{glmb} models currently only have implementations for the \code{dNormal()} 
#' prior family. The \code{Gamma()} family also works with the \code{dGamma()} prior 
#' family while the \code{gaussian()} family works with the \code{dGamma()} and 
#' \code{dNormal_Gamma()} pfamilies.  
#' @return An object of class \code{"pfamily"} (which has a concise print method). This is a
#' list with elements.
#' \item{pfamily}{character: the pfamily name}
#' \item{prior_list}{a list with the prior parameters associated with the prior specification}
#' \item{okfamilies}{currently implemented families for which the prior family can be used.}
#' \item{plinks}{a function that assigns a set of oklinks for the combination of a family and 
#' and pfamily.}
#' \item{simfun}{function: the function used to generate samples from the posterior density. 
#' All currently implemented pfamilies have simulation functions that generate iid samples
#' for the associated posterior distribution.}
#' 
#' @details
#' A `pfamily` object represents a structured prior specification for use in Bayesian generalized linear modeling. 
#' Each constructor function (e.g., `dNormal()`, `dGamma()`, `dNormal_Gamma()`) returns an object of class `"pfamily"` 
#' containing the prior parameters, supported likelihood families, compatible link functions, and a simulation function 
#' for posterior sampling.
#'
#' These priors are designed to integrate seamlessly with modeling functions such as `glmb()` and `rlmb()` in the 
#' \pkg{glmbayes} package, which consume the `pfamily` object to define the prior distribution over model parameters. 
#' The `pfamily()` generic retrieves the embedded prior from a fitted model object, while `print.pfamily()` displays its structure.
#'
#' ## Prior Families
#'
#' - **`dNormal()`**: Specifies a multivariate normal prior over regression coefficients. It is conjugate for 
#'   Gaussian likelihoods with an identity link function, and serves as the primary implemented prior for all 
#'   other supported likelihood families in the current framework. This structure facilitates efficient posterior 
#'   sampling and analytical tractability.
#'
#'   For models with log-concave likelihood functions—such as Poisson, Binomial, and Gamma families—
#'   posterior sampling under a `dNormal` prior is performed using a \insertCite{Nygren2006}{glmbayes} 
#'   likelihood subgradient approach. This method constructs tight enveloping functions around the posterior 
#'   using subgradients of the log-likelihood, enabling efficient accept-reject sampling even in high dimensions.
#'
#'   When the posterior distribution is approximately normal (typically the case for large sample sizes), the 
#'   area under the enveloping function is bounded above by a constant factor—approximately \eqn{2 / \sqrt{\pi} \approx 1.128} 
#'   in the univariate case, and \eqn{(2 / \sqrt{\pi})^k} in \eqn{k}-dimensional models. These bounds ensure that 
#'   the rejection rate remains manageable and that the sampler remains computationally efficient.
#'
#'   The concept of conjugate priors was first formalized by \insertCite{Raiffa1961}{glmbayes}, and further 
#'   developed for regression models using g-prior structures by \insertCite{zellner1986gprior}{glmbayes}.
#'
#' - **`dGamma()`**: Defines a gamma prior over a scalar precision parameter, often used in hierarchical models 
#'   or variance components. This prior is particularly relevant for Gamma likelihoods and dispersion modeling 
#'   in exponential families \insertCite{Gelman2013,Dobson1990,McCullagh1989}{glmbayes}.
#'
#' - **`dNormal_Gamma()`**: Combines a multivariate normal prior on coefficients with a gamma prior on precision, 
#'   forming a conjugate structure for Gaussian models with unknown variance. This formulation parallels classical 
#'   Normal-Gamma models and is compatible with hierarchical extensions \insertCite{Gelman2013,Raiffa1961}{glmbayes}.
#'
#' - **`dIndependent_Normal_Gamma()`**: Similar to `dNormal_Gamma()`, but assumes independence between the 
#'   coefficient and precision priors. This structure is useful for models where prior independence is desired 
#'   or analytically convenient.
#'
#' Each `pfamily` object includes:
#' - `prior_list`: A named list of prior parameters
#' - `okfamilies`: A character vector of compatible likelihood families
#' - `plinks`: A function returning valid link functions for a given family
#' - `simfun`: A simulation function used for posterior sampling
#' 
#' @references
#' \insertAllCited{}
#'   
#' @author The design of the \code{pfamily} set of functions was developed by Kjell Nygren and was 
#' inspired by the family used by the \code{\link{glmb}} function to specify the likelihood 
#' function. That design in turn was inspired by S functions of the same names described in
#' Hastie and Pregibon (1992).
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' 
#' @seealso
#' \code{\link{glmb}}, \code{\link{rlmb}}, \code{\link{lmb}}, \code{\link{rglmb}} for modeling functions that consume \code{pfamily} objects.
#'
#' \code{\link{rNormal_reg}}, \code{\link{rNormal_Gamma_reg}}, \code{\link{rGamma_reg}} for lower-level sampling functions used by \code{pfamily} constructors.
#'
#' \code{\link{Prior_Setup}}, \code{\link{Prior_Check}} for initializing and validating prior specifications.
#'
#' \code{\link{EnvelopeBuild}} for envelope construction methods used in likelihood subgradient sampling \insertCite{Nygren2006}{glmbayes}.
#'
#' See also \insertCite{Hastie1992}{glmbayes} for the original S modeling framework that inspired the design of \code{pfamily}.
#'  
#' @example inst/examples/Ex_pfamily.R
#' @export 
# #' @exportClass pfamily # Temporarily disabled - No Current exportclass
#' @rdname pfamily
#' @order 1

pfamily <- function(object, ...) UseMethod("pfamily")

#' @export 
#' @method pfamily default

pfamily.default <- function(object, ...){

  if(is.null(object$pfamily)) stop("no pfamily object found")
#  if(!class(object$pfamily)=="pfamily") stop("Object named pfamily is not of class pfamily")
  if (!inherits(object$pfamily, "pfamily"))  stop("Object named pfamily is not of class pfamily")
  
  return(object$pfamily)
}


#' @export 
#' @method print pfamily

#' @rdname pfamily
#' @order 6

print.pfamily <- function(x, ...)
{
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Prior Family:", x$pfamily, "\n\n")
  cat("Prior List:\n\n")
  print(x$prior_list)
  
  invisible(x)
}

#' @export 
#' @rdname pfamily
#' @order 2

dNormal<-function(mu,Sigma,dispersion=NULL){
  
  ## Check that the inputs are numeric
  
  if(is.numeric(mu)==FALSE||is.numeric(Sigma)==FALSE) stop("non-numeric argument to numeric function")

  mu=as.matrix(mu,ncol=1) ## Force mu to matrix
  Sigma=as.matrix(Sigma)  ## Force Sigma to matrix 
  
  nvar=length(mu)
  nvar1=nrow(Sigma)
  nvar2=ncol(Sigma)
  
  if(!nvar==nvar1||!nvar==nvar2) stop("dimensions of mu and Sigma are not consistent")

  ## Check for symmetry and positive definiteness
  if(!isSymmetric(Sigma))stop("matrix Sigma must be symmetric")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(Sigma, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  thr <- -tol * abs(ev[1])   # = -1e-06 * 12.56941 ≈ -1.256941e-05
  
  
  
  
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  
  if(is.null(dispersion)) dispersion=1
  if(!is.null(dispersion)){
    if(!is.numeric(dispersion)) stop("non-numeric argument to numeric function")
    if(!length(dispersion)==1) stop("dispersion has length>1")
    if(!length(dispersion)>0) stop("dispersion must be >0")
  }
    
  okfamilies <- c("gaussian","poisson","binomial","quasipoisson","quasibinomial","Gamma")

  plinks<-function(family){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-c("log")		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-c("logit","probit","cloglog")		
    if(family$family=="Gamma") oklinks<-c("log")	
    return(oklinks)
  }
  
  prior_list=list(mu=mu,Sigma=Sigma,dispersion=dispersion)
  attr(prior_list,"Prior Type")="dNormal"  

  outlist=list(pfamily="dNormal",prior_list=prior_list,okfamilies=okfamilies,
  plinks=plinks,             
  simfun=rNormal_reg)
  attr(outlist,"Prior Type")="dNormal"             
  class(outlist)="pfamily"
  outlist$call<-match.call()
  return(outlist)
  }

#' @export 
#' @rdname pfamily
#' @order 3

dGamma<-function(shape,rate,beta){

  if(is.numeric(shape)==FALSE||is.numeric(rate)==FALSE||is.numeric(beta)==FALSE) stop("non-numeric argument to numeric function")
  
  if(length(shape)>1) stop("shape is not of length 1")
  if(length(shape)>1) stop("rate is not of length 1")
  if(shape<=0) stop("shape must be>0")
  if(rate<=0) stop("rate must be>0")
  
  beta=as.matrix(beta,ncol=1)
  
  okfamilies <- c("gaussian","Gamma")
  
  plinks<-function(family){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-NULL		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-NULL		
    if(family$family=="Gamma") oklinks<-c("log")	
    return(oklinks)
  }
  
  prior_list=list(shape=shape,rate=rate,beta=beta)
  attr(prior_list,"Prior Type")="dGamma"  
  outlist=list(pfamily="dGamma",prior_list=prior_list,okfamilies=okfamilies,
               plinks=plinks,             
               simfun=rGamma_reg)
               
  attr(outlist,"Prior Type")="dGamma"
  class(outlist)="pfamily"
  outlist$call<-match.call()
  
  return(outlist)

}

#' @export 
#' @rdname pfamily
#' @order 4

dNormal_Gamma<-function(mu, Sigma,shape, rate){

  ############################################################  
  
  if(is.numeric(mu)==FALSE||is.numeric(Sigma)==FALSE) stop("non-numeric argument to numeric function")
  if(is.numeric(shape)==FALSE||is.numeric(rate)==FALSE) stop("non-numeric argument to numeric function")
  
  if(length(shape)>1) stop("shape is not of length 1")
  if(length(rate)>1) stop("rate is not of length 1")
  if(shape<=0) stop("shape must be>0")
  if(rate<=0) stop("rate must be>0")
  
  mu=as.matrix(mu,ncol=1) ## Force mu to matrix
  Sigma=as.matrix(Sigma)  ## Force Sigma to matrix 
    
  nvar=length(mu)
  nvar1=nrow(Sigma)
  nvar2=ncol(Sigma)
  
  if(!nvar==nvar1||!nvar==nvar2) stop("dimensions of mu and Sigma are not consistent")
  
  ## Check for symmetry and positive definiteness
  if(!isSymmetric(Sigma))stop("matrix Sigma must be symmetric")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(Sigma, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  
  
  ############################################################
  
  okfamilies <- c("gaussian") # Unclear if this could be used for Gamma  or quasi-families

  plinks<-function(family){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-NULL		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-NULL		
    if(family$family=="Gamma") oklinks<-NULL	
    return(oklinks)
  }
  
  prior_list=list(mu=mu,Sigma=Sigma,shape=shape,rate=rate)
  attr(prior_list,"Prior Type")="dNormal_Gamma"  
  outlist=list(pfamily="dNormal_Gamma",call=call,prior_list=prior_list,
    okfamilies=okfamilies,plinks=plinks,simfun=rNormal_Gamma_reg)
  
  attr(outlist,"Prior Type")="dNormal_Gamma"             
  class(outlist)="pfamily"
  outlist$call<-match.call()
  
  return(outlist)
  }



#' @export 
#' @rdname pfamily
#' @order 5

dIndependent_Normal_Gamma<-function(mu, Sigma,shape, rate){
  
  ##############################################################
  
  if(is.numeric(mu)==FALSE||is.numeric(Sigma)==FALSE) stop("non-numeric argument to numeric function")
  if(is.numeric(shape)==FALSE||is.numeric(rate)==FALSE) stop("non-numeric argument to numeric function")
  
  if(length(shape)>1) stop("shape is not of length 1")
  if(length(rate)>1) stop("rate is not of length 1")
  if(shape<=0) stop("shape must be>0")
  if(rate<=0) stop("rate must be>0")
  
  mu=as.matrix(mu,ncol=1) ## Force mu to matrix
  Sigma=as.matrix(Sigma)  ## Force Sigma to matrix 
  

  nvar=length(mu)
  nvar1=nrow(Sigma)
  nvar2=ncol(Sigma)
  
  if(!nvar==nvar1||!nvar==nvar2) stop("dimensions of mu and Sigma are not consistent")
  
  ## Check for symmetry and positive definiteness
  if(!isSymmetric(Sigma))stop("matrix Sigma must be symmetric")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(Sigma, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  
  
  ##############################################################
  
  okfamilies <- c("gaussian") # Unclear if this could be used for Gamma or quasi-families
  
  plinks<-function(family){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-NULL		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-NULL
    if(family$family=="Gamma") oklinks<-NULL	
    return(oklinks)
  }
  
  
  prior_list=list(mu=mu,Sigma=Sigma,shape=shape,rate=rate)
  attr(prior_list,"Prior Type")="dIndependent_Normal_Gamma"  
  outlist=list(pfamily="dIndependent_Normal_Gamma",prior_list=prior_list,
               okfamilies=okfamilies,plinks=plinks,simfun=rindependent_norm_gamma_reg)
  
  attr(outlist,"Prior Type")="dIndependent_Normal_Gamma"             
  class(outlist)="pfamily"
  outlist$call<-match.call()
  
  return(outlist)
  
}


