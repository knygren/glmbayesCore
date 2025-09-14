#' Simulation Functions for Bayesian Generalized Linear Models
#'
#' Simulation functions provide a unified interface for generating posterior samples from Bayesian GLMs. These functions are typically used within model fitting routines such as \code{\link{glmb}}, \code{\link{lmb}}, and \code{\link{rglmb}}, and support block Gibbs sampling and other simulation-based inference techniques.
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
#' @param use_parallel Logical. Whether to use parallel processing.
#' @param use_opencl Logical. Whether to use OpenCL acceleration.
#' @param verbose Logical. Whether to print progress messages.
#' @param digits Number of significant digits to use for printed output.
#' @param \ldots Additional arguments passed to or from other methods.
#'
#' @return
#' - \code{simfunction()}: An object of class \code{"simfunction"} containing:
#'   \itemize{
#'     \item{\code{name}}{Character string with the name of the simulation function}
#'     \item{\code{call}}{The matched call used to generate the simulation}
#'     \item{\code{args}}{A named list of arguments passed to the simulation function}
#'   }
#'
#' - \code{rGamma_reg()}: An object of class \code{"rGamma_reg"} containing:
#'   \itemize{
#'     \item{\code{coefficients}}{A 1 × p matrix of assumed regression coefficients}
#'     \item{\code{coef.mode}}{Currently \code{NULL}; reserved for future use}
#'     \item{\code{dispersion}}{A vector of simulated dispersion values}
#'     \item{\code{Prior}}{A list with prior parameters: \code{shape} and \code{rate}}
#'     \item{\code{prior.weights}}{Vector of prior weights used in the simulation}
#'     \item{\code{y}}{The response vector}
#'     \item{\code{x}}{The design matrix}
#'     \item{\code{famfunc}}{A processed family object used internally}
#'     \item{\code{iters}}{A vector indicating the number of iterations per sample (typically \code{rep(1, n)})}
#'     \item{\code{Envelope}}{Currently \code{NULL}; reserved for envelope diagnostics}
#'   }
#'   
#' @details
#' The \code{simfunction()} generic extracts metadata from simulation objects, including the function name, call, and arguments used. This is useful for introspection, reproducibility, and diagnostics.
#'
#' The lower-level simulation functions such as \code{rGamma_reg()} generate iid samples from posterior distributions for specific model components. These functions are used internally by \code{pfamily} constructors and model fitting routines.
#'
#' ## Simulation Functions
#'
#' - **\code{rGamma_reg()}**: Simulates dispersion parameters for Gaussian and Gamma families using either standard gamma sampling or accept-reject methods based on likelihood subgradients \insertCite{Nygren2006}{glmbayes}.
#'
#' - **\code{rNormal_reg()}**, **\code{rNormal_Gamma_reg()}**: Simulate regression coefficients and dispersion jointly or independently under Normal-Gamma priors.
#'
#' Each simulation function returns a structured object containing:
#' - Simulated values (e.g., coefficients, dispersion)
#' - Prior specification
#' - Model inputs (e.g., \code{y}, \code{x}, \code{weights})
#' - Metadata for reproducibility
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
#' @example inst/examples/Ex_rglmb_dispersion.R
#' @importFrom Rdpack reprompt
#' @rdname simfuncs
#' @order 1
#' @export
#' 
#' 


simfunction <- function(object, ...) {
  UseMethod("simfunction")
}



#' @export
#' @method simfunction default
#' @rdname simfuncs
#' @order 2
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
#' @order 3
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
#' @export 
#' @rdname simfuncs
#' @order 4
#' @export



rGamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                     Gridtype=2,
                     use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE
){
  
  call <- match.call()
  
  ## Renaming for consistency with earlier version
  
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
    
    mu1<-t(exp(alpha+x%*%b))
    
    testfunc<-function(v,wt){  
      -sum(lgamma(wt*v)+0.5*log(wt*v)+wt*v-wt*v*log(wt*v))
    }
    
    
    shape2=shape + 0.5 *n1
    rate1=rate +sum(wt*((y/mu1)-log(y/mu1)-1))
    
    vstar1<-shape2/rate1
    
    vout<-function(v){
      vstar1-(v/rate1)*sum((wt*digamma(wt*v) -wt*log(wt*v) + 0.5/v) )  
    }
    
    # Initialize vstar2
    vstar<-vstar1
    
    ## Optimize vstar?
    for(j in 1:20){
      vstar<-vout(vstar)
    }
    
    testbar<-testfunc(vstar,wt)
    cbar<--sum((wt*digamma(wt*vstar) -wt*log(wt*vstar) + 0.5/vstar))
    
    
    
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
#' @order 5
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
#' @order 6
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

