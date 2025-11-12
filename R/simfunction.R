#' @name simfuncs 
#' @title 
#' Simulation Functions for Bayesian Generalized Linear Models
#'
#' @description
#' Simulation functions provide a unified interface for generating posterior samples from Bayesian GLMs.
#'  These functions are typically used within model fitting routines such as \code{\link{rglmb}} and \code{\link{rlmb}}, and 
#'  are also suitable for use im Block Gibbs sampling and other simulation-based inference techniques.
#'
#'
#' @param object A fitted model object containing a \code{pfamily} component. The generic function \code{simfunction()} accesses the simulation metadata stored within such objects.
#' @param x An object of class \code{"simfunction"} or \code{"rGamma_reg"} to be printed.
#' @param n Number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y A vector of observations of length \code{m}.
#' @param x A design matrix of dimension \code{m × p}.
#' @param prior_list A list with prior parameters (e.g., shape, rate, beta) used in the simulation.
#' @param offset Optional numeric vector of length \code{m} specifying known components of the linear predictor.
#' @param weights Optional numeric vector of prior weights.
#' @param family A description of the error distribution and link function (see \code{\link{family}}).
#' @param Gridtype Optional integer specifying the method used to construct the envelope function.
#' @param n_envopt Effective sample size passed to EnvelopeOpt for grid
#'   construction. Defaults to match `n`. Larger values encourage tighter
#'   envelopes.
#' @param use_parallel Logical. Whether to use parallel processing.
#' @param use_opencl Logical. Whether to use OpenCL acceleration.
#' @param verbose Logical. Whether to print progress messages.
#' @param digits Number of significant digits to use for printed output.
#' @param progbar Logical. Whether to display a progress base during simulation.
#' @param \ldots Additional arguments passed to or from other methods.
#'
#' @return
#' \describe{
#'
#'   \item{\code{simfunction()}}{An object of class \code{"simfunction"} containing:
#'     \describe{
#'       \item{\code{name}}{Character string with the name of the simulation function.}
#'       \item{\code{call}}{The matched call used to generate the simulation.}
#'       \item{\code{args}}{A named list of arguments passed to the simulation function.}
#'     }
#'   }
#'
#'   \item{\code{rNormal_reg()}}{A list object with classes \code{"rglmb"}, \code{"glmb"}, \code{"glm"}, and \code{"lm"}.
#'     Elements include:
#'     \describe{
#'       \item{\code{coefficients}}{Matrix (\code{n × p}) of simulated regression coefficients, with column names from \code{x}.}
#'       \item{\code{coef.mode}}{Posterior mode of the coefficients. Gaussian: from \code{lm.fit}; non-Gaussian: BFGS mode shifted by prior mean.}
#'       \item{\code{dispersion}}{Scalar dispersion used. Poisson/Binomial: \code{1}; otherwise the supplied value. Quasi families: mean residual-based dispersion computed in the wrapper.}
#'       \item{\code{Prior}}{List with \code{mean} (prior mean vector) and \code{Precision} (prior precision matrix \code{P}).}
#'       \item{\code{prior.weights}}{Vector of prior weights used in the simulation (unscaled).}
#'       \item{\code{offset}}{Offset vector passed to the C++ sampler.}
#'       \item{\code{offset2}}{Offset used internally by the wrapper (copy of input or a zero vector).}
#'       \item{\code{y}}{Response vector.}
#'       \item{\code{x}}{Design matrix.}
#'       \item{\code{fit}}{Fitted/diagnostic object. Gaussian: result of \code{lm.fit} (class \code{"lm"}). Non-Gaussian: result of \code{glmb.wfit(...)}.}
#'       \item{\code{iters}}{Vector of iteration counts per sample. Gaussian: vector of ones; non-Gaussian: counts from the sampler.}
#'       \item{\code{Envelope}}{Envelope list used for accept-reject sampling (non-Gaussian); \code{NULL} for Gaussian.}
#'       \item{\code{family}}{Family object describing distribution and link.}
#'       \item{\code{famfunc}}{Processed family functions used internally (e.g., \code{f2}, \code{f3}).}
#'       \item{\code{call}}{Matched call to \code{rNormal_reg()}.}
#'       \item{\code{formula}}{Formula reconstructed from \code{y} and \code{x}.}
#'       \item{\code{model}}{Model frame corresponding to \code{formula}.}
#'       \item{\code{data}}{Data frame combining \code{y} and \code{x}.}
#'     }
#'   }
#'
#'   \item{\code{rNormal_Gamma_reg()}}{A list with class \code{"rglmb"} containing:
#'     \describe{
#'       \item{\code{coefficients}}{Matrix (\code{n × p}) of simulated regression coefficients; row \code{i} equals \code{Btilde + IR \%*\% rnorm(p) * sqrt(dispersion[i])}. Column names are set to \code{colnames(x)}.}
#'       \item{\code{coef.mode}}{Posterior mean/mode vector \code{Btilde} from \code{rNormal_reg.wfit()}.}
#'       \item{\code{dispersion}}{Numeric vector of length \code{n} with draws from the inverse-gamma posterior \code{1/rgamma(shape = shape + nobs/2, rate = rate + 0.5*S)}.}
#'       \item{\code{Prior}}{List with \code{mean} (as numeric vector \code{mu}) and \code{Precision} (matrix \code{P}).}
#'       \item{\code{offset}}{Offset vector as supplied.}
#'       \item{\code{prior.weights}}{Vector of prior weights \code{wt}.}
#'       \item{\code{y}}{Response vector.}
#'       \item{\code{x}}{Design matrix.}
#'       \item{\code{fit}}{Result from \code{rNormal_reg.wfit()}, including fields such as \code{Btilde}, \code{IR}, \code{S}, and \code{k}.}
#'       \item{\code{famfunc}}{Processed family functions for Gaussian models (from \code{glmbfamfunc(gaussian())}).}
#'       \item{\code{iters}}{Numeric vector (length \code{n}) of ones indicating per-draw iteration counts.}
#'       \item{\code{Envelope}}{\code{NULL}; no envelope is constructed in this conjugate setup.}
#'       \item{\code{call}}{Matched call to \code{rNormal_Gamma_reg()}.}
#'     }
#'   }
#'
#'   \item{\code{rindependent_norm_gamma_reg()}}{A list with class \code{"rglmb"} containing:
#'     \describe{
#'       \item{\code{coefficients}}{Matrix (\code{n × p}) of simulated regression coefficients, back-transformed to the original scale; column names set to \code{colnames(x)}.}
#'       \item{\code{coef.mode}}{Vector with the conditional posterior mode used for envelope anchoring (from the Gaussian fit).}
#'       \item{\code{dispersion}}{Numeric vector of length \code{n} with simulated dispersion draws.}
#'       \item{\code{Prior}}{List with prior components: \code{mean} (prior mean \code{mu}), \code{Sigma} (prior covariance), \code{shape} and \code{rate} (Gamma prior for dispersion), \code{Precision} (\code{solve(Sigma)}).}
#'       \item{\code{family}}{The \code{gaussian()} family object.}
#'       \item{\code{prior.weights}}{Vector of prior weights used in the simulation.}
#'       \item{\code{y}}{Response vector.}
#'       \item{\code{x}}{Design matrix.}
#'       \item{\code{call}}{Matched call to \code{rindependent_norm_gamma_reg()}.}
#'       \item{\code{famfunc}}{Processed family functions for Gaussian models (from \code{glmbfamfunc}).}
#'       \item{\code{iters}}{Vector with per-draw iteration counts returned by the joint sampler.}
#'       \item{\code{Envelope}}{\code{NULL}; envelope diagnostics are not returned by this function.}
#'       \item{\code{loglike}}{\code{NULL}; placeholder for log-likelihood values.}
#'       \item{\code{weight_out}}{Numeric vector of per-draw weights returned by the C++ routine.}
#'       \item{\code{sim_bounds}}{List with \code{low} and \code{upp}, the dispersion bounds used by the shared envelope.}
#'       \item{\code{offset2}}{Offset vector used internally (copy of input or a zero vector).}
#'     }
#'   }
#'
#'   \item{\code{rGamma_reg()}}{An object of class \code{"rGamma_reg"} containing:
#'     \describe{
#'       \item{\code{coefficients}}{A 1 × p matrix of assumed regression coefficients.}
#'       \item{\code{coef.mode}}{Currently \code{NULL}; reserved for future use.}
#'       \item{\code{dispersion}}{A vector of simulated dispersion values.}
#'       \item{\code{Prior}}{A list with prior parameters: \code{shape} and \code{rate}.}
#'       \item{\code{prior.weights}}{Vector of prior weights used in the simulation.}
#'       \item{\code{y}}{The response vector.}
#'            }
#'            }
#'            }
#'      
#'         
#' @details The low-level simulation functions **\code{rNormal_reg()}**, **\code{rNormal_Gamma_reg()}**, 
#' **\code{rindependent_norm_gamma_reg()}**, and **\code{rGamma_reg()}** generate iid samples from posterior 
#' distributions for specific model components. These model functions are used internally by the functions
#' **\code{rglmb()}** and **\code{rlmb()}** to generate samples.  
#'  
#' The \code{simfunction()} generic extracts metadata from simulation objects, including the function name, call, and arguments used. This is useful for introspection, reproducibility, and diagnostics.
#'
#' The lower-level simulation functions generate iid samples from posterior distributions for specific model components. 
#' These functions are used internally by \code{pfamily} constructors and model fitting routines.
#'
#' ## Simulation Functions
#'
#' - **\code{rNormal_reg()}**: Produces iid draws for regression coefficients in models with
#'   multivariate normal priors and log-concave likelihood functions. For Gaussian likelihoods,
#'   these are conjugate priors and standard simulation procedures for multivariate normal
#'   distributions are utilized \insertCite{LindleySmith1972,DiaconisYlvisaker1979}{glmbayes}.
#'   For all other families/link functions, the likelihood subgradient approach of
#'   \insertCite{Nygren2006}{glmbayes} is used to generate iid samples.
#'
#' - **\code{rNormal_Gamma_reg()}**: Produces iid draws for regression coefficients and the
#'   dispersion parameter in models with Normal-Gamma priors and Gaussian likelihoods, where
#'   this is a conjugate prior distribution. Standard simulation procedures for gamma and
#'   multivariate normal distributions are utilized
#'   \insertCite{Raiffa1961,LindleySmith1972}{glmbayes}.
#'
#' - **\code{rindependent_norm_gamma_reg()}**: Produces iid draws for regression coefficients
#'   and the dispersion parameter in models with independent Normal and truncated Gamma priors.
#'   This is a non-conjugate specification but can still be sampled using accept-reject procedures
#'   based on an enveloping approach (see vignette
#'   \insertCite{glmbayesIndNormGammaVignette}{glmbayes}).
#'
#' - **\code{rGamma_reg()}**: Simulates dispersion parameters for Gaussian and Gamma families
#'   using either standard gamma sampling or accept-reject methods based on likelihood
#'   subgradients \insertCite{Chen1979,glmbayesGammaVignette}{glmbayes}.
#'
#' @references
#' \insertAllCited{}
#' 
#' @author
#' The simulation framework was developed by Kjell Nygren as part of the \pkg{glmbayes} package. It builds on the likelihood subgradient approach described in \insertCite{Nygren2006}{glmbayes}, and extends classical Bayesian GLM sampling techniques.
#'
#' @seealso
#' \code{\link{pfamily}}, \code{\link{glmb}}, \code{\link{lmb}}, \code{\link{rglmb}} for modeling functions that consume simulation functions.
#'
#' \code{\link{rNormal_reg}}, \code{\link{rNormal_Gamma_reg}}, \code{\link{rGamma_reg}} for individual simulation functions.
#'
#' \code{\link{EnvelopeBuild}} for envelope construction methods used in likelihood subgradient sampling.
#'





#' @usage simfunction(object, ...)
#' @export simfunction
#' @rdname simfuncs
#' @order 1

simfunction <- function(object, ...) {
  UseMethod("simfunction")
}



#' @method simfunction default
#' @noRd
#' @export

simfunction.default <- function(object, ...) {
  if (is.null(object$pfamily)) stop("no pfamily object found")
  if (!inherits(object$pfamily, "pfamily")) stop("Object named pfamily is not of class pfamily")
  
  pf <- object$pfamily
  simfun <- pf$simfun
  
  simfun_name <- "anonymous or not found"
  fun_env <- environment(simfun)
  fun_names <- ls(fun_env)
  for (name in fun_names) {
    if (identical(simfun, get(name, envir = fun_env))) {
      simfun_name <- name
      break
    }
  }
  
  simfun_call <- if (!is.null(object$simfun_call)) object$simfun_call else NULL
  simfun_args <- if (!is.null(object$simfun_args)) object$simfun_args else list()
  
  structure(
    list(
      name = simfun_name,
      call = simfun_call,
      args = simfun_args
    ),
    class = "simfunction"
  )
}


#' @export
#' @method print simfunction
#' @rdname simfuncs
#' @order 9
print.simfunction <- function(x, ...) {
  cat("\nCall to Simulation Function:\n")
  if (!is.null(x$call)) {
    print(x$call)
  } else {
    cat("  [call not recorded]\n")
  }
  
  cat("\nSimulation Function Name:", x$name, "\n")
  
  if (!is.null(x$args) && length(x$args) > 0) {
    cat("\nArguments Passed:\n\n")
    for (argname in names(x$args)) {
      val <- x$args[[argname]]
      
      if (is.null(val)) {
        cat("  ", argname, ": [NULL]\n", sep = "")
      } else if (argname == "family") {
        cat("  ", argname, ":\n", sep = "")
        print(val)
      } else if (argname == "prior_list" && is.list(val)) {
        cat("  prior_list:\n")
        for (pname in names(val)) {
          pval <- val[[pname]]
          cat("    ", pname, ":\n", sep = "")
          if (is.null(pval)) {
            cat("      [NULL]\n")
          } else if (is.atomic(pval) || is.matrix(pval)) {
            print(pval)
          } else {
            cat("      [", class(pval), " with length ", length(pval), "]\n", sep = "")
          }
        }
      } else {
        cat("  ", argname, ":\n", sep = "")
        if (is.atomic(val) || is.matrix(val)) {
          print(val)
        } else {
          cat("    [", class(val), " with length ", length(val), "]\n", sep = "")
        }
      }
    }
  } else {
    cat("\nArguments Passed: [none recorded]\n")
  }
  
  invisible(x)
}


#' @family simfuncs 
#' @references A reference
#' @example inst/examples/Ex_rglmb_dispersion.R
#' @usage rGamma_reg(n, y, x, prior_list, offset = NULL, weights = 1, family = gaussian(),
#'            Gridtype = 2,n_envopt = NULL,
#'             use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE)
#' @export 
#' @rdname simfuncs
#' @order 5
#' @export



rGamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                     Gridtype=2,n_envopt = NULL,
                     use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE
){
  
  call <- match.call()
  

  ## Argument renaming and prior
  
  wt=weights
  alpha=offset
  
  b=prior_list$beta
  shape=prior_list$shape
  rate=prior_list$rate
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  okfamilies <- c("gaussian","Gamma")
  if(family$family %in% okfamilies){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="Gamma") oklinks<-c("log")		
    if(family$link %in% oklinks)  {}
    else{stop(gettextf("link \"%s\" not available for selected family; available links are %s", 
                       family$link , paste(sQuote(oklinks), collapse = ", ")), 
              domain = NA)
    }
  }
  
  else{
    stop(gettextf("family \"%s\" not available in glmbdisp; available families are %s", 
                  family$family , paste(sQuote(okfamilies), collapse = ", ")), 
         domain = NA)
    
  }
  
  n1<-length(y)
  
  if(family$family=="gaussian"){
    y1<-as.matrix(y)-alpha
    xb<-x%*%b
    res<-y1-xb
    SS<-res*res
    
    a1<-shape+n1/2
    b1<-rate+sum(SS)/2
    
    out<-1/rgamma(n,shape=a1,rate=b1) 
  }
  
  if(family$family=="Gamma")
  {
    
    ## Compute mu1 using the fixed coefficients
    
    mu1<-t(exp(alpha+x%*%b))
    
    ## testfunc is part of the log-likelihood (excluding a constant and the part included in the 
    ## ) 
    
    testfunc<-function(v,wt){  
      -sum(lgamma(wt*v)+0.5*log(wt*v)+wt*v-wt*v*log(wt*v))
    }
    
    ## Update the shape and rate using the fitted values from the likelihood
    
    shape2=shape + 0.5 *n1
    rate1=rate +sum(wt*((y/mu1)-log(y/mu1)-1))
    
    ## Initialize vstar1 to the ratio (i.e., posterior mode)
    
    vstar1<-shape2/rate1
    
    ## Use Newton method to update vstar1 and to solve for posterior mode
    
    vout<-function(v){
      vstar1-(v/rate1)*sum((wt*digamma(wt*v) -wt*log(wt*v) + 0.5/v) )  
    }
    
    # Initialize vstar2
    vstar<-vstar1
    
    ## Optimize vstar?
    for(j in 1:20){
      vstar<-vout(vstar)
    }
    
    ## Find value of testfunc at the posterior mode
    ## and the negative of the gradient of the testfunc at vstar
        
    testbar<-testfunc(vstar,wt)
    cbar<--sum((wt*digamma(wt*vstar) -wt*log(wt*vstar) + 0.5/vstar))
    
    ## Set the rate to correspond to the posterior mode
    
    rate2=  rate +sum(wt*((y/mu1)-log(y/mu1)-1))-sum((wt*digamma(wt*vstar) -wt*log(wt*vstar) + 0.5/vstar) )
    
    out<-matrix(0,n)
    test<-matrix(0,n)
    a<-matrix(0,n)
    
    ## Implements rejection sampling for dispersion (likelihood subgradient approach)
    ## Likely should have a short paper with this derivation
    ## Not sure if approach extends to other densities besides gamma
    
    for(i in 1:n)
    {
      while(a[i]==0){
        out[i]<-rgamma(1,shape=shape2,rate=rate2)
        
        test[i]<-testfunc(out[i],wt)-(testbar+cbar*(out[i]-vstar))-log(runif(1,0,1))
        if(test[i]>0) a[i]<-1
      }
    }
    
    ## Convert "Precision" to dispersion
    out<-1/out
    
  }
  
  outlist=list(
    coefficients=matrix(b,nrow=1,ncol=length(b)),
    coef.mode=NULL,
    dispersion=out,
    Prior=list(shape=shape,rate=rate),
    prior.weights=weights,
    y=y,
    x=x,
    famfunc=glmbfamfunc(family),
    iters=rep(1,n),
    Envelope=NULL
  )
  
  
  outlist$call<-match.call()
  
  class(outlist)<-c(outlist$class,"rGamma_reg")
  
  return(outlist)
  
}


#' @export
#' @rdname simfuncs
#' @order 6
#' @method print rGamma_reg

print.rGamma_reg<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Simulated Dispersion")
    cat(":\n")
    print.default(format(x$dispersion, digits = digits), 
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
}


#' @export
#' @rdname simfuncs
#' @order 7
#' @method summary rGamma_reg


summary.rGamma_reg<-function(object,...){
  
  n<-length(object$dispersion)  
  percentiles<-matrix(0,nrow=1,ncol=7)
  me=mean(object$dispersion)
  se<-sqrt(var(object$dispersion))
  mc<-se/n
  Priorwt<-(se/(sqrt(object$Prior$shape)/object$Prior$rate))^2
  percentiles[1,]<-quantile(object$dispersion,probs=c(0.01,0.025,0.05,0.5,0.95,0.975,0.99))
  test<-append(object$dispersion,object$Prior$shape/object$Prior$rate)
  test2<-rank(test)
  priorrank<-test2[n+1]
  pval1<-priorrank/(n+1)
  pval2<-min(pval1,1-pval1)
  
  
  Tab1<-cbind("Prior.Mean"=object$Prior$shape/object$Prior$rate,"Prior.Sd"=sqrt(object$Prior$shape)/object$Prior$rate
              ,"Approx.Prior.wt"=Priorwt
  )
  TAB<-cbind(
    #"Post.Mode"=as.numeric(object$PostMode),
    "Post.Mean"=me,
    "Post.Sd"=se,
    "MC Error"=as.numeric(mc)
    ,"Pr(tail)"=as.numeric(pval2)
  )
  TAB2<-cbind("1.0%"=percentiles[,1],"2.5%"=percentiles[,2],"5.0%"=percentiles[,3],Median=as.numeric(percentiles[,4]),"95.0%"=percentiles[,5],"97.5%"=as.numeric(percentiles[,6]),"99.0%"=as.numeric(percentiles[,7]))
  
  rownames(TAB)=c("dispersion")
  rownames(Tab1)=c("dispersion")
  rownames(TAB2)=c("dispersion")
  
  res<-list(call=object$call,
            n=n,
            coefficients1=Tab1,
            coefficients=TAB,
            Percentiles=TAB2
  )
  
  # Reuse summary.rglmb class
  
  class(res)<-"summary.rglmb"
  
  res
  
}




#' @family simfuncs 
#' @example inst/examples/Ex_rindep_norm_gamma_reg.R
#' @usage rindependent_norm_gamma_reg(n, y, x, prior_list, offset = NULL, weights = 1,
#'                              family = gaussian(), Gridtype = 2,n_envopt = NULL,
#'                               use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, 
#'                              progbar = TRUE)
#' @export 
#' @rdname simfuncs
#' @order 4



rindependent_norm_gamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                                      Gridtype=2,n_envopt = NULL,
                                      use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE,
                                      progbar=TRUE){
  
  call<-match.call()
  
  offset2=offset
  wt=weights
  
  if(length(wt)==1) wt=rep(wt,length(y))
  
  ### Initial implementation of Likelihood subgradient Sampling 
  ### Currently uses as single point for conditional tangencis
  ### (at conditional posterior modes)
  ### Verify this yields correct results and then try to implement grid approach
  
  ## Use the prior list to set the prior elements if it is not missing
  ## Error checking to verify that the correct elements are present
  ## Shold be implemented
  
  
  ## Step 1: Validate Prior Specification
  
  if(missing(prior_list)) stop("Prior Specification Missing")
  if(!missing(prior_list)){
    if(!is.null(prior_list$mu)) mu=prior_list$mu
    if(!is.null(prior_list$Sigma)) Sigma=prior_list$Sigma
    if(!is.null(prior_list$dispersion)) dispersion=prior_list$dispersion
    else dispersion=NULL
    if(!is.null(prior_list$shape)) shape=prior_list$shape
    else shape=NULL
    if(!is.null(prior_list$rate)) rate=prior_list$rate
    else rate=NULL
    if (!is.null(prior_list$max_disp_perc)) {
      max_disp_perc <- prior_list$max_disp_perc
    } else {
      max_disp_perc <- 0.99
    }
    
    ## New: extract optional low/upp from prior_list
    if (!is.null(prior_list$disp_lower))  disp_lower <- prior_list$disp_lower  else disp_lower <- NULL
    if (!is.null(prior_list$disp_upper))  disp_upper <- prior_list$disp_upper  else disp_upper <- NULL
    
    ## Validation if both are provided
    if (!is.null(disp_lower) && !is.null(disp_upper)) {
      if (!is.numeric(disp_lower) || !is.numeric(disp_upper)) {
        stop("prior_list$disp_lower and prior_list$disp_upper must be numeric.")
      }
      if (disp_lower <= 0 || disp_upper <= 0) {
        stop("prior_list$disp_lower and prior_list$disp_upper must be positive.")
      }
      if (disp_upper <= disp_lower) {
        stop("prior_list$disp_upper must be strictly greater than prior_list$disp_lower.")
      }
    }
    
  }
  

  ## Step 2: Fit Classical Model and Extract Estimates and RSS

  lm_out=lm(y ~ x-1,weights=wt,offset=offset2) # run classical regression to get maximum likelhood estimate
  
  ### 2) Initi1alize Residuals and dispersion

  RSS=sum(residuals(lm_out)^2)
  
  RSS_ML=sum(residuals(lm_out)^2)
  n_obs=length(y)

  dispersion2=dispersion
  RSS_temp<-matrix(0,nrow=1000)
  
  ## Step 3: Iterative Dispersion Anchoring (Finds good value for the dispersion)
  
  for(j in 1:10){
    glmb_out1=glmb(y~x-1,family=gaussian(),
                   dNormal(mu=mu,Sigma=Sigma,dispersion=dispersion2),weights=wt,offset=offset2)
    
    res_temp=residuals(glmb_out1)
    
    for(k in 1:1000){
      RSS_temp[k]=sum(res_temp[k,1:length(y)]*res_temp[k,1:length(y)])
    }
    
    RSS_Post2=mean(RSS_temp)
    b_old=glmb_out1$coef.mode
    xbetastar=x%*%b_old
    RSS2_post=t(y-xbetastar)%*%(y-xbetastar)  
    shape2= shape + n_obs/2
    #  rate2 =rate + RSS2_post/2  # 38 candidates per acceptance
    rate2=rate + RSS_Post2/2    # 35.7 candidates per acceptance [though some with positive acceptance]
    
    dispersion2=rate2/(shape2-1)
    
  }
  
  
  
  ## Step 4: Standardized Model 
  ## 4) Steps to standardized the model (i.e. to reorient dimensions)
  
  betastar=glmb_out1$coef.mode
  dispstar=rate2/(shape2-1)
  
  
##  dispstar <- rate3 / (shape2 - 1)
  
  
  ## Get family functions for gaussian()  
  
  famfunc<-glmbfamfunc(gaussian())
  
  f1<-famfunc$f1
  f2<-famfunc$f2
  f3<-famfunc$f3
  f5<-famfunc$f5
  f6<-famfunc$f6
  
  start <- mu
  
  if(is.null(offset2))  offset2=rep(as.numeric(0.0),length(y))
  P=solve(Sigma)
  
  ###### Adjust weight for dispersion
  
  dispersion2=dispstar
  
  if(is.null(wt))  wt=rep(1,length(y))
  if(length(wt)==1)  wt=rep(wt,length(y))
  
  wt2=wt/rep(dispersion2,length(y))
  
  ######################### Shift mean vector to offset so that adjusted model has 0 mean
  
  alpha=x%*%as.vector(mu)+offset2
  mu2=0*as.vector(mu)
  P2=P
  x2=x
  
  parin=start-mu
  #parin=start
  
  #parin=glmb_out1$glm$coefficients
  
  # This step only used to get posterior precision it seems
  # Since normal case, can likely be computed instead
  
  opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x2),mu=as.vector(mu2),
                P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
                method="BFGS",hessian=TRUE
  )
  
  bstar=opt_out$par  ## Posterior mode for adjusted model
  
  ## Temporarily use bstar as posterior mode
  
  
  #  bstar
  #  bstar+as.vector(mu)  # mode for actual model
  A1=opt_out$hessian # Approximate Precision at mode
  
  
  Standard_Mod=glmb_Standardize_Model(y=as.vector(y), x=as.matrix(x2),
                                      P=as.matrix(P2),bstar=as.matrix(bstar,ncol=1), A1=as.matrix(A1))
  
  bstar2=Standard_Mod$bstar2  
  A=Standard_Mod$A
  x2=Standard_Mod$x2
  mu2=Standard_Mod$mu2  ## Typically 0 vector 
  P2=Standard_Mod$P2    ## Part of Prior that is moved to the log-likelihood
  L2Inv=Standard_Mod$L2Inv
  L3Inv=Standard_Mod$L3Inv
  
  
  
  ## Step 5: Build Model
  
  ## Build initial Envelope based on the optimized values
  
  ## Note, use Gridtype =4 here temporarily (Single Likelihood subgradient)
  
  #Gridtype=as.integer(Gridtype)
  
  ## Pull the initial Envelope based on optimized values above
  
  Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
                     as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),
                     as.vector(wt2),
                     family="gaussian",link="identity",
                     Gridtype=Gridtype, n=as.integer(n),
                     n_envopt=n_envopt,
                     sortgrid=TRUE,use_opencl = use_opencl,verbose = verbose)
  
  cat("[DEBUG] EnvelopeBuild returned\n")
  
  
  ###  Call new function to build shared envelope
  
  cat("[DEBUG] Entering EnvelopeDispersionBuild \n")
  
  
  # disp_env_out_old <- EnvelopeDispersionBuild(
  #   Env        = Env2,
  #   Shape      = shape,
  #   Rate       = rate,
  #   P          = P2,
  #   y=y,
  #   x          = x2,
  #   alpha      =as.vector(alpha),
  #   n_obs      = n_obs,
  #   RSS_post   = RSS_Post2,
  #   RSS_ML     =RSS_ML,
  #   max_disp_perc = max_disp_perc,
  #   disp_lower=disp_lower,
  #   disp_upper=disp_upper
  # )
  
  disp_env_out <- EnvelopeDispersionBuild_cpp(
    Env        = Env2,
    Shape      = shape,
    Rate       = rate,
    P          = P2,
    y          = y,
    x          = x2,
    alpha      = as.vector(alpha),
    n_obs      = n_obs,
    RSS_post   = RSS_Post2,
    RSS_ML     = RSS_ML,
    mu=as.matrix(mu2,ncol=1),
#    wt=as.vector(wt2),
    wt=as.vector(wt),
    max_disp_perc = max_disp_perc,
    disp_lower = disp_lower,
    disp_upper = disp_upper,
    verbose    = TRUE   # optional, matches the C++ signature
  )
  
  cat("[DEBUG] Exiting EnvelopeDispersionBuild \n")
  

  Env3_raw       <- disp_env_out$Env_out
  gamma_list_new <- disp_env_out$gamma_list
  UB_list_new    <- disp_env_out$UB_list
  low            <- gamma_list_new$disp_lower
  upp            <- gamma_list_new$disp_upper
  diagnostics    <- disp_env_out$diagnostics
  
  # Decide how many rows to keep
  l1 <- ncol(x2)
  l2 <- nrow(Env3_raw$cbars)
  
  logP_mat <- if (is.null(dim(Env3_raw$logP))) {
    # Gaussian path: coerce vector to one-column matrix
    as.matrix(Env3_raw$logP)
  } else {
    Env3_raw$logP
  }
  
  Env3 <- EnvelopeSort(
    l1      = ncol(Env3_raw$cbars),
    l2      = nrow(Env3_raw$cbars),
    GIndex  = Env3_raw$GridIndex,
    G3      = Env3_raw$thetabars,
    cbars   = Env3_raw$cbars,
    logU    = Env3_raw$logU,
    logrt   = Env3_raw$logrt,
    loglt   = Env3_raw$loglt,
    logP    = logP_mat,
    LLconst = Env3_raw$LLconst,
    PLSD    = Env3_raw$PLSD,
    a1      = Env3_raw$a1,
    E_draws = Env3_raw$E_draws,
    lg_prob_factor = disp_env_out$UB_list$lg_prob_factor
  )

  
  
  
  # cat("[DEBUG] disp_lower =", low,
  #     " disp_upper =", upp, "\n")
  # 
  # cat("[DEBUG] Calling .rindep_norm_gamma_reg_std_V4_cpp \n")
  
  # sim_temp <- .rindep_norm_gamma_reg_std_V4_cpp(
  #   n=n, y=y, x=x2, mu=mu2, P=P2, alpha=alpha, wt,
  #   f2=f2, Envelope=Env3, 
  #   gamma_list=gamma_list_new,
  #   UB_list=UB_list_new,
  #   family="gaussian", link="identity", progbar=progbar
  # )
  
  
  
  

  # Env3_temp           <- disp_env_out$Env_out
  gamma_list_new <- disp_env_out$gamma_list
  UB_list_new    <- disp_env_out$UB_list
  low            <- gamma_list_new$disp_lower
  upp            <- gamma_list_new$disp_upper  
  diagnostics     <- disp_env_out$diagnostics
  
  ## USE reordered lg_prob factor
  UB_list_new$lg_prob_factor=Env3$lg_prob_factor
  
  

    
    cat("[DEBUG] disp_lower =", low,
      " disp_upper =", upp, "\n")
  
  cat("[DEBUG] Calling .rindep_norm_gamma_reg_std_V4_cpp \n")
  
  
  
  sim_temp=.rindep_norm_gamma_reg_std_cpp (n=n, y=y, x=x2, 
                                              mu=mu2, ## Should be zero vector
                                              P=P2,  ## Part of prior shifted to the Likelihood
                                              alpha=alpha, wt,
                                              f2=f2, Envelope=Env3, 
                                              gamma_list=gamma_list_new,
                                              UB_list=UB_list_new,
                                              family="gaussian",link="identity", progbar =progbar)

  cat("[DEBUG] Exiting .rindep_norm_gamma_reg_std_V4_cpp \n")
  

  print(paste("Interactive status:", interactive()))
  


  #proc.time()-ptm
  
  
  beta_out=sim_temp$beta_out
  disp_out=sim_temp$disp_out
  iters_out=sim_temp$iters_out
  weight_out=sim_temp$weight_out
  
  out=L2Inv%*%L3Inv%*%t(beta_out)
  
  for(i in 1:n){
    out[,i]=out[,i]+mu
  }
  
  
  famfunc=glmbfamfunc(gaussian())  
  f1=famfunc$f1
  
  outlist=list(
    coefficients=t(out), 
    coef.mode=betastar,  ## For now, use the conditional mode (not universal)
    dispersion=disp_out,
    ## For now, name items in list like this-eventually make format/names
    ## consistent with true prior (current names needed by summary function)
    Prior=list(mean=mu,Sigma=Sigma,shape=shape,rate=rate,Precision=solve(Sigma)), 
    family=gaussian(),
    prior.weights=wt,
    y=y,
    x=x,
    call=call,
    famfunc=famfunc,
    iters=iters_out,
    Envelope=NULL,
    loglike=NULL,
    weight_out=weight_out,
    sim_bounds=list(low=low,upp=upp)
    #,test_out=test_out
  )
  
  colnames(outlist$coefficients)<-colnames(x)
  outlist$offset2<-offset2
  class(outlist)<-c(outlist$class,"rglmb")
  
  return(outlist)  
  
  
}




## This function is used by the above (not sure why Neg_logLik is not working)
## Could be because it is not exported - replace with Neg_logLik
#' @rdname Neg_logLik
#' @export 

Neg_logLik2<-function(b, y, x, alpha, wt,family){
  
  ## Add required checks on other inputs at the top
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  okfamilies <- c("gaussian","poisson","binomial","quasipoisson","quasibinomial","Gamma")
  if(family$family %in% okfamilies){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-c("log")		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-c("logit","probit","cloglog")		
    if(family$family=="Gamma") oklinks<-c("log")		
    if(family$link %in% oklinks){
      
      ## This may be the R version of these files so may not be using the efficiency of C++
      ## This may be safer
      
      famfunc<-glmbfamfunc(family)
      f1<-famfunc$f1
      f2<-famfunc$f2
      f3<-famfunc$f3
      #      f5<-famfunc$f5
      #      f6<-famfunc$f6
    }
    else{
      stop(gettextf("link \"%s\" not available for selected family; available links are %s", 
                    family$link , paste(sQuote(oklinks), collapse = ", ")), 
           domain = NA)
      
    }	
    
  }		
  else {
    stop(gettextf("family \"%s\" not available in glmb; available families are %s", 
                  family$family , paste(sQuote(okfamilies), collapse = ", ")), 
         domain = NA)
  }
  
  return(f1(b, y, x, alpha, wt)) 
  
}



#' @family simfuncs 
#' @example inst/examples/Ex_rnorm_gamma_reg.R
#' @usage rNormal_Gamma_reg(n, y, x, prior_list, offset = NULL, weights = 1, family = gaussian(),
#'                   Gridtype = 2,n_envopt = NULL, 
#'                   use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE)
#' @export 
#' @rdname simfuncs
#' @order 3


rNormal_Gamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                            Gridtype=2,n_envopt = NULL,
                            use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE
){
  
  ## Added for consistency with earlier verion of function
  
  offset2=offset
  wt=weights
  
  ## Below code used precision matrix (not Sigma)
  ## Code checks for the presence of P in the prior
  ## if not present, it imputes Precision by inverting the Sigma matrix
  
  if(missing(prior_list)) stop("Prior Specification Missing")
  if(!missing(prior_list)){
    if(!is.null(prior_list$mu)) mu=prior_list$mu
    if (!is.null(prior_list$Sigma)) {
      Sigma <- prior_list$Sigma
      if (!isSymmetric(Sigma)) 
        stop("matrix Sigma must be symmetric")
    }
    
    if (!is.null(prior_list$P)) {
      P <- prior_list$P
      if (!isSymmetric(P)) 
        stop("matrix P must be symmetric")
    }
    
    if (is.null(prior_list$P)) {
      P <- solve(prior_list$Sigma)
      # enforce symmetry on the computed precision
      P <- 0.5 * (P + t(P))
    }
    if(!is.null(prior_list$dispersion)) dispersion=prior_list$dispersion
    else dispersion=NULL
    if(!is.null(prior_list$shape)) shape=prior_list$shape
    else shape=NULL
    if(!is.null(prior_list$rate)) rate=prior_list$rate
    else rate=NULL
  }
  
  if(is.numeric(n)==FALSE||is.numeric(y)==FALSE||is.numeric(x)==FALSE||
     is.numeric(mu)==FALSE||is.numeric(P)==FALSE) stop("non-numeric argument to numeric function")
  
  x <- as.matrix(x)
  mu<-as.matrix(as.vector(mu))
  P<-as.matrix(P)    
  
  ## Allow function to be called without offset2
  
  if(length(n)>1) n<-length(n)	   
  nobs <- NROW(y)
  nvars <- ncol(x)
  
  if(is.null(offset2)) offset2=rep(0,nobs)
  nvars2<-length(mu)	
  if(!nvars==nvars2) stop("incompatible dimensions")
  if (!all(dim(P) == c(nvars2, nvars2))) 
    stop("incompatible dimensions")
  if(!isSymmetric(P))stop("matrix P must be symmetric")
  if(length(wt)==1) wt=rep(wt,nobs)
  nobs2=NROW(wt)
  nobs3=NROW(x)
  nobs4=NROW(offset2)
  if(nobs2!=nobs) stop("weighting vector must have same number of elements as y")
  if(nobs3!=nobs) stop("matrix X must have same number of rows as y")
  if(nobs4!=nobs) stop("offset vector must have same number of rows as y")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(P, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'P' is not positive definite")
  
  ## Should add dimension checks here  
  ## Should move core part of rmultireg inside this code to avoid call
  
  famfunc=glmbfamfunc(gaussian())  
  f1=famfunc$f1
  
  fit=rNormal_reg.wfit(x,y,P,mu, w=wt,offset=offset2,method="qr",tol=1e-7,singular.ok=TRUE)
  Btilde=fit$Btilde
  IR=fit$IR
  S=fit$S
  k=fit$k
  m=1  ## Used to fix dimension to 1 --> Will eventually do away with
  
  out1<-matrix(0,nrow=n,ncol=k)
  
  if(m==1){  out2=matrix(0,nrow=n,ncol=1)  }
  if(m>1){  out2<-vector("list", n)}   # Code if we had multiple columns of data (should perhaps have m items)
  
  a_prior=shape     ## Should be relationship to shape in Wishart  
  b_prior=rate  ## Should be relationship to scale in Wishart (could also be V/2)
  
  a_post=a_prior+(nobs/2) # Posterior Shape parameter
  b_post=b_prior+0.5*S # Posterior rate  (S is scaled differently than V?)
  
  out1<-matrix(0,nrow=n,ncol=k)
  dispersion=1/rgamma(n=n,shape=a_post,rate=b_post)
  out2<-matrix(dispersion,nrow=n,ncol=1)
  
  for(i in 1:n){out1[i,1:k]<- t(Btilde + IR%*%matrix(rnorm(m*k),ncol=m)*sqrt(dispersion[i])) }
  
  draws<-matrix(1,n)
  ##LL<-matrix(1,n)
  
  #for(i in 1:n){
  ## This function should return the negative log-likelihood 
  #  LL[i]=f1(b=out1[i,],y=y,x=x,alpha=offset2,wt=wt/out2[i])	
  #}
  
  outlist<-list(
    coefficients=out1
    ,coef.mode=Btilde,
    dispersion=dispersion,
    offset=offset,
    Prior=list(mean=as.numeric(mu),Precision=P),
    prior.weights=wt,
    y=y,
    x=x,
    fit=fit,
    famfunc=famfunc,
    iters=draws,
    Envelope=NULL
    # , loglike=LL
  )
  
  
  colnames(outlist$coefficients)<-colnames(x)
  
  outlist$call<-match.call()
  
  class(outlist)<-c(outlist$class,"rglmb")
  
  return(outlist)
  
  
}





#' @family simfuncs 
#' @example inst/examples/Ex_rnorm_gamma_reg.R
#' @usage rNormal_reg(n, y, x, prior_list, offset = NULL, weights = 1,
#'             family = gaussian(), Gridtype = 2, n_envopt = NULL,
#'             use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE)
#' @export 
#' @rdname simfuncs
#' @order 2



rNormal_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                      Gridtype=2,n_envopt = NULL,
                      use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE){
  
  ## Added for consistency with earlier verion of function
  ## Useful to copy offset and weight and then to modify to non-null as needed  
  
  offset2=offset
  wt=weights
  
  ## Missing control variables (add option to pass these)
  ## Setting for default Gridtype might be important
  
  #Gridtype=2
  
  
  ## Below code used precision matrix (not Sigma)
  ## Code checks for the presence of P in the prior
  ## if not present, it imputes Precision by inverting the Sigma matrix
  
  if(missing(prior_list)) stop("Prior Specification Missing")
  if(!missing(prior_list)){
    if(!is.null(prior_list$mu)) mu=prior_list$mu
    if(!is.null(prior_list$Sigma)) Sigma=prior_list$Sigma
    if(!is.null(prior_list$P)) P=prior_list$P
    if(is.null(prior_list$P)) P=(solve(prior_list$Sigma)+t(solve(prior_list$Sigma)))/2
    if(!is.null(prior_list$dispersion)) dispersion=prior_list$dispersion
    else dispersion=NULL
    if(!is.null(prior_list$shape)) shape=prior_list$shape
    else shape=NULL
    if(!is.null(prior_list$rate)) rate=prior_list$rate
    else rate=NULL
  }
  
  
  if(is.numeric(n)==FALSE||is.numeric(y)==FALSE||is.numeric(x)==FALSE||
     is.numeric(mu)==FALSE||is.numeric(P)==FALSE) stop("non-numeric argument to numeric function")
  
  # normalize n_envopt
  if (is.null(n_envopt)) n_envopt <- n
  n_envopt <- as.integer(n_envopt)
  
  
  x <- as.matrix(x)
  mu<-as.matrix(as.vector(mu))
  P<-as.matrix(P)    
  
  ## Start value should be contingent on the family and link
  
  start <- mu
  
  ## Allow function to be called without offset2
  
  if(length(n)>1) n<-length(n)	   
  nobs <- NROW(y)
  nvars <- ncol(x)
  
  if(is.null(offset2)) offset2=rep(0,nobs)
  nvars2<-length(mu)	
  if(!nvars==nvars2) stop("incompatible dimensions")
  if (!all(dim(P) == c(nvars2, nvars2))) 
    stop("incompatible dimensions")
  if(!isSymmetric(P))stop("matrix P must be symmetric")
  if(length(wt)==1) wt=rep(wt,nobs)
  nobs2=NROW(wt)
  nobs3=NROW(x)
  nobs4=NROW(offset2)
  if(nobs2!=nobs) stop("weighting vector must have same number of elements as y")
  if(nobs3!=nobs) stop("matrix X must have same number of rows as y")
  if(nobs4!=nobs) stop("offset vector must have same number of rows as y")
  
  tol<- 1e-06 # Link this to Magnitude of P	
  eS <- eigen(P, symmetric = TRUE,only.values = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'P' is not positive definite")
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  okfamilies <- c("gaussian","poisson","binomial","quasipoisson","quasibinomial","Gamma")
  if(family$family %in% okfamilies){
    if(family$family=="gaussian") oklinks<-c("identity")
    if(family$family=="poisson"||family$family=="quasipoisson") oklinks<-c("log")		
    if(family$family=="binomial"||family$family=="quasibinomial") oklinks<-c("logit","probit","cloglog")		
    if(family$family=="Gamma") oklinks<-c("log")		
    if(family$link %in% oklinks){
      
      famfunc<-glmbfamfunc(family)
      f1<-famfunc$f1
      f2<-famfunc$f2
      f3<-famfunc$f3
      #      f5<-famfunc$f5
      #      f6<-famfunc$f6
    }
    else{
      stop(gettextf("link \"%s\" not available for selected family; available links are %s", 
                    family$link , paste(sQuote(oklinks), collapse = ", ")), 
           domain = NA)
      
    }	
    
  }		
  else {
    stop(gettextf("family \"%s\" not available in glmb; available families are %s", 
                  family$family , paste(sQuote(okfamilies), collapse = ", ")), 
         domain = NA)
  }
  
  
  if(family$family=="gaussian"){ 
    outlist<-.rnorm_reg_cpp(n=n,y=y,x=x,mu=mu,P=P,offset=offset2,wt=wt,dispersion=dispersion,
                            ##                      famfunc=famfunc,f1=f1,
                            f2=f2,f3=f3,start=mu)
    class(outlist$fit)="lm"
    
  }
  else{
    if(is.null(dispersion)){dispersion2=1}
    else{dispersion2=dispersion}
    
    #  stop("Inputs to function above")
    outlist<-.rnnorm_reg_cpp(n=n,y=y,x=x,mu=mu,P=P,offset=offset2,wt=wt,
                             dispersion=dispersion2,
                             ##famfunc=famfunc,f1=f1,
                             f2=f2,f3=f3,
                             start=start,family=family$family,link=family$link,Gridtype=Gridtype,
                             n_envopt = n_envopt,       # pass through
                             use_parallel = use_parallel,
                             use_opencl = use_opencl,
                             verbose = verbose)
    
    
    betastar=outlist$coef.mode  # Posterior mode from optim
    x=outlist$x
    y=outlist$y
    #offset=glmb.D93$offset   # not present in the output --> For now set to 0 vector
    #offset=offset2   # Should return this from lower level functions
    weights=outlist$prior.weights
    
    
    if(family$family=="quasipoisson"||family$family=="quasibinomial"){
      
      
      linkinv<-family$linkinv
      
      ## Compute dispersion and then rerun
      disp_temp=rep(0,n)
      m=length(y)  
      k=ncol(x)    
      res_temp=matrix(0,nrow=n,ncol=m)
      fit_temp=x%*%t(outlist$coefficients)
      for(l in 1:n){
        #fit_temp[1:m,l]=exp(offset2+fit_temp[1:m,l])
        fit_temp[1:m,l]=linkinv(offset2+fit_temp[1:m,l])
        res_temp[l,1:m]=(y-fit_temp[1:m,l])
        disp_temp[l]=(1/(m-k))*sum(res_temp[l,1:m]^2*wt/fit_temp[1:m,l])
        
      }
      
      #  stop("Inputs to function above")
      # Rerun model with updated dispersion    
      
      outlist<-.rnnorm_reg_cpp(n=n,y=y,x=x,mu=mu,P=P,offset=offset2,
                               #                             wt=wt/mean(disp_temp),
                               wt=wt,
                               dispersion=mean(disp_temp),
                               #                              dispersion=dispersion2,
                               ##famfunc=famfunc,f1=f1,
                               f2=f2,f3=f3,
                               start=start,family=family$family,link=family$link,Gridtype=Gridtype,
                               n_envopt = n_envopt,       # pass through
                               use_parallel = use_parallel,
                               use_opencl = use_opencl,
                               verbose = verbose)
      
      outlist$call <- match.call()  # overwrite with the rNormal_reg call
      outlist$dispersion=mean(disp_temp)
      
      #print("Mean_residuals")
      #print(colMeans(res_temp))
      #print("Weights")
      #print(wt)
      #print("mean_fit_temp")
      #print(rowMeans(fit_temp))
      #print("Mean Dispersion")
      #print(mean(disp_temp))
      
      # Update the dispersion
      #outlist$dispersion=disp_temp
    }
    
    if(family$family=="quasibinomial"){
      
      print("Hello quasibinomial")
    }
    
    
    ## get influence info for original model
    outlist$fit=glmb.wfit(x,y,weights,offset=offset2,family=family,Bbar=mu,P,betastar)
    
    
  }
  
  
  colnames(outlist$coefficients)<-colnames(x)
  
  # include family in final list
  
  rglmb_df=as.data.frame(cbind(y,x))
  rglmb_f=DF2formula(rglmb_df)
  rglmb_mf=model.frame(rglmb_f,rglmb_df)
  
  outlist$family=family
  outlist$famfunc=famfunc
  outlist$call<-match.call()
  outlist$offset2<-offset2
  outlist$formula<-rglmb_f
  outlist$model<-rglmb_mf
  outlist$data<-rglmb_df
  
  class(outlist)<-c(outlist$class,c("rglmb","glmb","glm","lm"))
  outlist
  
}








################################## Utility functions used by the above  #################

p_inv_gamma<-function(dispersion,shape,rate){
  1-pgamma(1/dispersion,shape=shape,rate=rate)
}

q_inv_gamma<-function(p,shape,rate,disp_upper,disp_lower){
  p_upp=p_inv_gamma(disp_upper,shape=shape,rate=rate)
  p_low=p_inv_gamma(disp_lower,shape=shape,rate=rate)
  p1=p_low+p*(p_upp-p_low)
  p2=1-p1
  1/qgamma(p2,shape,rate)
}

r_invgamma<-function(n,shape,rate,disp_upper,disp_lower){
  p=runif(n)
  q_inv_gamma(p=p,shape=shape,rate=rate,disp_upper=disp_upper,disp_lower)
}

