#' Prior Family Objects for Bayesian Models
#'
#' Prior family objects provide a convenient way to specify the details of the priors 
#' used by functions such as \code{\link{glmb}}. See the documentations for \code{\link{lmb}},
#' \code{\link{glmb}}, \code{\link{glmb}}, and \code{\link{rglmb}} for the details of how such model fitting 
#' takes place.
#' @name pfamily
#' @param object the function \code{pfamily} accesses the \code{pfamily} objects which
#' are stored within objects created by modelling functions (e.g., \code{glmb}).
#' @param mu a prior mean vector for the the modeling coefficients used in several pfamilies
#' @param Sigma a prior variance-covariance matrix for \code{dNormal()} and
#'   \code{dIndependent_Normal_Gamma()}.
#' @param Sigma_0 prior variance-covariance on the precision-weighted coefficient scale for
#'   \code{dNormal_Gamma()} only (Gaussian). Stored in \code{prior_list$Sigma} for compatibility
#'   with downstream samplers.
#' @param dispersion the dispersion to be assumed when it is not given a prior. Should be provided
#' when the Normal prior is for the \code{gaussian()}, \code{Gamma()}, \code{quasibinomial},
#' or \code{quasipoisson} families. The \code{binomial()} and \code{poisson()} families
#' do not have dispersion coefficients. Omitted or \code{NULL} uses the internal default
#' \code{1} and sets \code{ddef} in \code{prior_list} (see Details).
#' @param shape The prior shape parameter for the gamma piece (inverse dispersion / precision).
#'   When taking defaults from \code{\link{Prior_Setup}}, use \code{ps$shape} with
#'   \code{\link{dNormal_Gamma}()} and \code{\link{dGamma}()}, and \code{ps$shape_ING} with
#'   \code{\link{dIndependent_Normal_Gamma}()} on the Gaussian calibrated path (see Details).
#' @param rate The prior rate parameter paired with \code{shape}. With Gaussian
#'   \code{\link{Prior_Setup}}, \code{\link{dNormal_Gamma}()} and \code{\link{dIndependent_Normal_Gamma}()}
#'   use \code{ps$rate}; for \code{\link{dGamma}()} with fixed \code{beta}, prefer \code{ps$rate_gamma}
#'   when that field is non-\code{NULL} (see Details).
#' @param beta the regression coefficients to be assumed when it is not given a prior. 
#' Needs to be provided when the Gamma prior is used for the dispersion. This
#' specification is typically only used as part of Gibbs sampling where the beta and 
#' dispersion parameters are updated separately. 
#' @param max_disp_perc Specifies the percentile used to truncate the posterior dispersion 
#' distribution when constructing the envelope for accept-reject sampling. This determines 
#' the lower and upper bounds for the dispersion (\eqn{\sigma^2}) used in the simulation. A value of 0.99
#' corresponds to using the central 98 percent of the posterior dispersion mass (i.e., excluding 
#' the outer 1 percent in each tail). Smaller values yield tighter bounds and may improve acceptance 
#' rates, while larger values allow broader dispersion support but may increase envelope complexity.
#' @param disp_lower lower bound truncation for dispersion 
#' @param disp_upper upper bound truncation for dispersion
#' @param x an object, a pfamily function that is to be printed
#' @param \ldots additional argument(s) for methods.
#' @details
#' \code{pfamily} is a generic with methods for fitted objects such as \code{\link{glmb}} and
#' \code{\link{lmb}}. Many \code{glmb} models currently only implement the \code{dNormal()}
#' prior family. The \code{Gamma()} response family works with \code{dGamma()}; the
#' \code{gaussian()} family works with \code{dGamma()} and \code{dNormal_Gamma()}.
#'
#' A `pfamily` object represents a structured prior specification for use in Bayesian generalized linear modeling. 
#' Each constructor function (e.g., `dNormal()`, `dGamma()`, `dNormal_Gamma()`) returns an object of class `"pfamily"` 
#' containing the prior parameters, supported likelihood families, compatible link functions, and a simulation function 
#' for posterior sampling.
#'
#' These priors are designed to integrate seamlessly with modeling functions such as `glmb()` and `rlmb()` in the 
#' \pkg{glmbayes} package, which consume the `pfamily` object to define the prior distribution over model parameters. 
#' The `pfamily()` generic retrieves the embedded prior from a fitted model object, while `print.pfamily()` displays its structure.
#'
#' **\code{prior_list} and \code{simfun}.** The named list \code{prior_list} holds the hyperparameters for the chosen
#' prior family. When a model function draws from the posterior, it passes \code{prior_list} into the element
#' \code{simfun} (e.g., \code{\link{rNormal_reg}}, \code{\link{rGamma_reg}}) so the low-level sampler receives
#' one consistent list structure regardless of which constructor built the \code{pfamily}.
#'
#' **\code{\link{Prior_Setup}} and default hyperparameters.** \code{Prior_Setup()} fits an auxiliary GLM and returns
#' default \code{mu}, \code{Sigma} / \code{Sigma_0}, \code{dispersion}, Gamma \code{shape} and \code{rate}, and
#' related fields aligned with the data and prior-weight (\code{pwt}) choices. Those values can be supplied as
#' arguments to the \code{pfamily} constructors when you want package-default priors on the same scale as the model
#' matrix. Recommended use of \code{shape} and \code{rate} is not identical across constructors: for
#' \code{\link{dIndependent_Normal_Gamma}()}, pass \code{shape = ps$shape_ING} from \code{Prior_Setup} (not the
#' scalar \code{ps$shape} used by \code{\link{dNormal_Gamma}()}). For \code{\link{dGamma}()} with fixed coefficients
#' (\code{beta}), pass \code{rate = ps$rate_gamma} when that field is present (otherwise \code{ps$rate}); see
#' \code{\link{Prior_Setup}} and \code{\link{compute_gaussian_prior}}.
#'
#' ## Prior Families
#'
#' - **`dNormal()`**: Specifies a multivariate normal prior over regression coefficients. It is conjugate for 
#'   Gaussian likelihoods with an identity link function, and serves as the primary implemented prior for all 
#'   other supported likelihood families in the current framework. This structure facilitates efficient posterior 
#'   sampling and analytical tractability. The returned \code{prior_list} includes \code{ddef}: \code{TRUE} when
#'   \code{dispersion} was omitted or \code{NULL} (so the default \code{1} was used), \code{FALSE} when
#'   \code{dispersion} was supplied explicitly (including \code{1}).
#'
#'   For models with log-concave likelihood functions-such as Poisson, Binomial, and Gamma families-
#'   posterior sampling under a `dNormal` prior is performed using a \insertCite{Nygren2006}{glmbayes} 
#'   likelihood subgradient approach. This method constructs tight enveloping functions around the posterior 
#'   using subgradients of the log-likelihood, enabling efficient accept-reject sampling even in high dimensions.
#'
#'   When the posterior distribution is approximately normal (typically the case for large sample sizes), the 
#'   area under the enveloping function is bounded above by a constant factor-approximately \eqn{2 / \sqrt{\pi} \approx 1.128} 
#'   in the univariate case, and \eqn{(2 / \sqrt{\pi})^k} in \eqn{k}-dimensional models. These bounds ensure that 
#'   the rejection rate remains manageable and that the sampler remains computationally efficient.
#'
#'   The concept of conjugate priors was first formalized by \insertCite{Raiffa1961}{glmbayes}, and further 
#'   developed for regression models using g-prior structures by \insertCite{zellner1986gprior}{glmbayes}.
#'
#' - **`dGamma()`**: Defines a gamma prior over a scalar precision parameter, often used in hierarchical models 
#'   or variance components. This prior is particularly relevant for Gamma likelihoods and dispersion modeling 
#'   in exponential families \insertCite{Gelman2013,Dobson1990,McCullagh1989}{glmbayes}.
#'   With Gaussian \code{\link{Prior_Setup}} output, prefer \code{rate_gamma} for \code{rate} when updating
#'   dispersion with fixed \code{beta} (see Details above).
#'
#' - **`dNormal_Gamma()`**: Combines a multivariate normal prior on coefficients with a gamma prior on precision, 
#'   forming a conjugate structure for Gaussian models with unknown variance. The second argument is \code{Sigma_0}
#'   (precision-weighted scale); it is aliased internally to \code{Sigma} in \code{prior_list}.
#'   This formulation parallels classical Normal-Gamma models and is compatible with hierarchical extensions \insertCite{Gelman2013,Raiffa1961}{glmbayes}.
#'   
#'
#' - **`dIndependent_Normal_Gamma()`**: Similar to `dNormal_Gamma()`, but assumes independence between the 
#'   coefficient and precision priors. This structure is useful for models where prior independence is desired 
#'   or analytically convenient. With \code{\link{Prior_Setup}} on a Gaussian model, pass \code{shape_ING} as the
#'   \code{shape} argument (see Details above).
#'
#' Each `pfamily` object includes:
#' - `pfamily`, `prior_list`, `okfamilies`, `plinks`, and `simfun` (see Value).
#'
#' @return An object of class \code{"pfamily"} (with a concise \code{print} method). A list with elements:
#' \item{pfamily}{Character string: the constructor name (\code{"dNormal"}, \code{"dGamma"},
#'   \code{"dNormal_Gamma"}, or \code{"dIndependent_Normal_Gamma"}).}
#' \item{prior_list}{Named list of prior hyperparameters. It is passed into \code{simfun} when sampling so the
#'   relevant low-level routine receives the prior in a fixed list form. Contents depend on the constructor:
#'   \describe{
#'     \item{\code{dNormal}:}{\code{mu}, \code{Sigma}, \code{dispersion}, and logical \code{ddef}
#'       (\code{TRUE} if \code{dispersion} was omitted or \code{NULL}, so the default \code{1} was used;
#'       \code{FALSE} if set explicitly).}
#'     \item{\code{dGamma}:}{\code{shape}, \code{rate}, \code{beta}, \code{max_disp_perc},
#'       \code{disp_lower}, \code{disp_upper}.}
#'     \item{\code{dNormal_Gamma}:}{\code{mu}, \code{Sigma} (the \code{Sigma_0} precision-weighted input),
#'       \code{shape}, \code{rate}.}
#'     \item{\code{dIndependent_Normal_Gamma}:}{\code{mu}, \code{Sigma} (coefficient-scale covariance),
#'       \code{shape}, \code{rate}, \code{max_disp_perc}, \code{disp_lower}, \code{disp_upper}.}
#'   }
#' }
#' \item{okfamilies}{Character vector of implemented \code{\link[stats]{family}} names for which this
#'   \code{pfamily} may be used.}
#' \item{plinks}{Function of one \code{family} argument returning allowed link names for that family.}
#' \item{simfun}{Function used to generate posterior draws (e.g., \code{\link{rNormal_reg}},
#'   \code{\link{rGamma_reg}}, \code{\link{rGamma_Conjugate_reg}}, \code{\link{rNormalGamma_reg}}, \code{\link{rindepNormalGamma_reg}});
#'   for standard use these produce i.i.d.\ posterior samples for the implemented settings.}
#' 
#' @author The design of the \code{pfamily} set of functions was developed by Kjell Nygren and was 
#' inspired by the family used by the \code{\link{glmb}} function to specify the likelihood 
#' function. That design in turn was inspired by S functions of the same names from
#' the statistical modeling literature.
#'
#' @seealso
#' \code{\link{glmb}}, \code{\link{rlmb}}, \code{\link{lmb}}, \code{\link{rglmb}} for modeling functions that consume \code{pfamily} objects.
#'
#' \code{\link{rNormal_reg}}, \code{\link{rNormalGamma_reg}}, \code{\link{rGamma_reg}}, \code{\link{rGamma_Conjugate_reg}}, \code{\link{rindepNormalGamma_reg}} for lower-level sampling functions used by \code{pfamily} constructors.
#'
#' \code{\link{Prior_Setup}}, \code{\link{Prior_Check}} for initializing and validating prior specifications.
#'
#' \code{\link{EnvelopeBuild}} for envelope construction methods used in likelihood subgradient sampling \insertCite{Nygren2006}{glmbayes}.
#'
#' See also \insertCite{Hastie1992}{glmbayes} for the original S modeling framework that inspired the design of \code{pfamily}.
#'
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#'
#' @example inst/examples/Ex_pfamily.R
#' @export 
#' @rdname pfamily
#' @order 1

pfamily <- function(object, ...) UseMethod("pfamily")

#' @export 
#' @method pfamily default

pfamily.default <- function(object, ...){

  if(is.null(object$pfamily)) stop("no pfamily object found")
  if (!inherits(object$pfamily, "pfamily"))  stop("Object named pfamily is not of class pfamily")
  
  return(object$pfamily)
}


#' @export
#' @method print pfamily
#' @rdname pfamily
#' @order 7

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
  thr <- -tol * abs(ev[1])   # = -1e-06 * 12.56941 ~= -1.256941e-05  
  
  
  
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  
  ddef <- missing(dispersion) || is.null(dispersion)
  if (ddef) dispersion <- 1
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
  
  prior_list=list(mu=mu,Sigma=Sigma,dispersion=dispersion,ddef=ddef)
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

dGamma<-function(shape,rate,beta,max_disp_perc = 0.99,disp_lower=NULL,disp_upper=NULL){

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
  
  prior_list=list(shape=shape,rate=rate,beta=beta,max_disp_perc = max_disp_perc,disp_lower=disp_lower,disp_upper=disp_upper)
  attr(prior_list,"Prior Type")="dGamma"  
  outlist=list(pfamily="dGamma",prior_list=prior_list,okfamilies=okfamilies,
               plinks=plinks,             
               simfun=rGamma_reg)
               
  attr(outlist,"Prior Type")="dGamma"
  class(outlist)="pfamily"
  outlist$call<-match.call()
  
  return(outlist)

}

#' Conjugate Gamma prior family (stub; cloned from \code{\link{dGamma}} pending specialization).
#'
#' @export
#' @rdname pfamily
#' @order 4

dGamma_Conjugate<-function(shape,rate,beta,max_disp_perc = 0.99,disp_lower=NULL,disp_upper=NULL){

  if(is.numeric(shape)==FALSE||is.numeric(rate)==FALSE||is.numeric(beta)==FALSE) stop("non-numeric argument to numeric function")

  if(length(shape)>1) stop("shape is not of length 1")
  if(length(shape)>1) stop("rate is not of length 1")
  if(shape<=0) stop("shape must be>0")
  if(rate<=0) stop("rate must be>0")

  beta=as.matrix(beta,ncol=1)

  okfamilies <- c("poisson", "Gamma")

  ## rglmb() calls pfamily$plinks(family); return allowed links by family (conjugate / identity scale).
  plinks <- function(family) {
    oklinks <- NULL
    if (family$family %in% c("poisson", "quasipoisson")) oklinks <- c("identity")
    if (family$family == "Gamma") oklinks <- c("identity")
    if (family$family %in% c("binomial", "quasibinomial")) oklinks <- NULL
    ## gaussian / other families: oklinks stays NULL until supported here
    oklinks
  }

  prior_list=list(shape=shape,rate=rate,beta=beta,max_disp_perc = max_disp_perc,disp_lower=disp_lower,disp_upper=disp_upper)
  attr(prior_list,"Prior Type")="dGamma_Conjugate"  
  outlist=list(pfamily="dGamma_Conjugate",prior_list=prior_list,okfamilies=okfamilies,
               plinks=plinks,             
               simfun=rGamma_Conjugate_reg)

  attr(outlist,"Prior Type")="dGamma_Conjugate"
  class(outlist)="pfamily"
  outlist$call<-match.call()

  return(outlist)

}



#' @export 
#' @rdname pfamily
#' @order 5

dNormal_Gamma <- function(mu, Sigma_0, shape, rate) {
  Sigma <- Sigma_0

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
    okfamilies=okfamilies,plinks=plinks,simfun=rNormalGamma_reg)
  
  attr(outlist,"Prior Type")="dNormal_Gamma"             
  class(outlist)="pfamily"
  outlist$call<-match.call()
  
  return(outlist)
  }



#' @export 
#' @rdname pfamily
#' @order 6

dIndependent_Normal_Gamma <- function(mu, Sigma, shape, rate, max_disp_perc = 0.99,disp_lower=NULL,disp_upper=NULL) {

  ##############################################################
  
  if(is.numeric(mu)==FALSE||is.numeric(Sigma)==FALSE) stop("non-numeric argument to numeric function")
  if(is.numeric(shape)==FALSE||is.numeric(rate)==FALSE) stop("non-numeric argument to numeric function")
  
  if(length(shape)>1) stop("shape is not of length 1")
  if(length(rate)>1) stop("rate is not of length 1")
  if(shape<=0) stop("shape must be>0")
  if(rate<=0) stop("rate must be>0")
  if (!is.numeric(max_disp_perc) || length(max_disp_perc) != 1 || max_disp_perc <= 0.5 || max_disp_perc >= 1) {
    stop("max_disp_perc must be a single number between 0.5 and 1")
  }
  
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
  
  
  prior_list <- list(mu = mu, Sigma = Sigma, shape = shape, rate = rate, max_disp_perc = max_disp_perc,
                     disp_lower=disp_lower,disp_upper=disp_upper)
  attr(prior_list,"Prior Type")="dIndependent_Normal_Gamma"  
  outlist=list(pfamily="dIndependent_Normal_Gamma",prior_list=prior_list,
               okfamilies=okfamilies,plinks=plinks,simfun=rindepNormalGamma_reg)
  
  attr(outlist,"Prior Type")="dIndependent_Normal_Gamma"             
  class(outlist)="pfamily"
  outlist$call<-match.call()
  
  return(outlist)
  
}


