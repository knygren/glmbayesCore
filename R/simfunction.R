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
#' @param use_parallel Logical. Whether to use parallel processing.
#' @param use_opencl Logical. Whether to use OpenCL acceleration.
#' @param verbose Logical. Whether to print progress messages.
#' @param digits Number of significant digits to use for printed output.
#' @param max_disp_perc parameter currently used to control upper bound in accept-reject procedure 
#' @param progbar Logical. Whether to display a progress base during simulation.
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
#' @details The low-level simulation functions **\code{rNormal_reg()}**, **\code{rNormal_Gamma_reg()}**, 
#' **\code{rindependent_norm_gamma_reg()}**, and **\code{rGamma_reg()}** generate iid samples from posterior 
#' distributions for specific model components. These model functions are used internally by the functions
#' **\code{rglmb()}** and **\code{rlmb()}** to generate samples.  
#'  
#' The \code{simfunction()} generic extracts metadata from simulation objects, including the function name, call, and arguments used. This is useful for introspection, reproducibility, and diagnostics.
#'
#' The lower-level simulation functions such as \code{rGamma_reg()} generate iid samples from posterior distributions for specific model components. 
#' These functions are used internally by \code{pfamily} constructors and model fitting routines.
#'
#' ## Simulation Functions
#'
#'
#' - **\code{rNormal_reg()}**: Produces iid draws for regression coefficients in models with Multivariate normal priors and Log-concave likelihood functions. For gaussian likelihood functions,
#'  these are conjugate prior distributions and standard simulation procedures for multivariate normal distributoosn 
#'  are utilized. For all other families/link functions, the Likelihood subgrandient approach of \insertCite{Nygren2006}{glmbayes} is used to generate iid samples.  
#' 
#' - **\code{rNormal_Gamma_reg()}**: Produces iid draws for regression coefficients and the dispersion parameter in models 
#' with Normal-Gamma priors and gaussian likelihood functions where this is a conjugate prior distribution. Standard simulation procedures
#' for gamma distributions and Multivariate normal distributions are utilized.
#'
#' - **\code{rindependent_norm_gamma_reg()}**: Produces iid draws for regression coefficients and the dispersion parameter in models 
#' with independent Normal and Gamma priors for the two. This is a non-conjugate specification but can still be sampled using accept reject procedures
#' based on an enveloping approach (see Vignette covering this topic).
#' 
#' - **\code{rGamma_reg()}**: Simulates dispersion parameters for Gaussian and Gamma families using either standard gamma sampling or accept-reject methods based on likelihood subgradients \insertCite{Nygren2006}{glmbayes}.
#'
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
#'            Gridtype = 2, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE)
#' @export 
#' @rdname simfuncs
#' @order 5
#' @export



rGamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                     Gridtype=2,
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
#'                              family = gaussian(), Gridtype = 2, use_parallel = TRUE,
#'                              use_opencl = FALSE, verbose = FALSE, max_disp_perc = 0.99,
#'                              progbar = TRUE)
#' @export 
#' @rdname simfuncs
#' @order 4



rindependent_norm_gamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                                      Gridtype=2,
                                      use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE,
                                      max_disp_perc=0.99,
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
  mu2=Standard_Mod$mu2
  P2=Standard_Mod$P2
  L2Inv=Standard_Mod$L2Inv
  L3Inv=Standard_Mod$L3Inv
  
  
  
  ## Step 5: Build Model
  
  ## Build initial Envelope based on the optimized values
  
  ## Note, use Gridtype =4 here temporarily (Single Likelihood subgradient)
  
  Gridtype=as.integer(3)
  
  ## Pull the initial Envelope based on optimized values above
  
  
  Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
                     as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),
                     as.vector(wt2),
                     family="gaussian",link="identity",
                     Gridtype=Gridtype, n=as.integer(n),
                     sortgrid=TRUE)
  
  
  ###  Call new function to build shared envelope
  
  
  disp_env_out <- EnvelopeDispersionBuild(
    Env        = Env2,
    Shape      = shape,
    Rate       = rate,
    P          = P2,
    y=y,
    x          = x2,
    alpha      =as.vector(alpha),
    n_obs      = n_obs,
    RSS_post   = RSS_Post2,
    RSS_ML     =RSS_ML,
    max_disp_perc = max_disp_perc
  )
  
  
##  Env2$PLSD       <- disp_env_out$Env_out$PLSD  # updated mixture weights
##  UB_list         <- disp_env_out$UB_list
  diagnostics     <- disp_env_out$diagnostics
  
  # Optional: extract gamma parameters if needed
##  shape3          <- disp_env_out$diagnostics$shape3
##  rate3           <- disp_env_out$diagnostics$rate3
  
  Env3           <- disp_env_out$Env_out
  gamma_list_new <- disp_env_out$gamma_list
  UB_list_new    <- disp_env_out$UB_list
  low            <- gamma_list_new$disp_lower
  upp            <- gamma_list_new$disp_upper  

  

  
  
  sim_temp=.rindep_norm_gamma_reg_std_V4_cpp (n=n, y=y, x=x2, mu=mu2, P=P2, alpha=alpha, wt,
                                              f2=f2, Envelope=Env3, 
                                              gamma_list=gamma_list_new,
                                              UB_list=UB_list_new,
                                              family="gaussian",link="identity", progbar =progbar)
  
  # 
  # ## Extract outputs
  # beta_out  <- sim_temp$beta_out
  # disp_out  <- sim_temp$disp_out
  # iters_out <- sim_temp$iters_out
  # 
  # 
  # 
  #     
  # #### End of Standardization - Steps below might change  
  # 
  # ## Update Gamma parameters [Should this be just RSS_Post2?]
  # 
  # shape2= shape + n_obs/2
  # rate2 =rate + (RSS_ML/2)
  # #rate3 =rate + (RSS2_post/2)
  # rate3 =rate + (RSS_Post2/2)
  # 
  # 
  # 
  # #################  this Block Does evaluations at lower and upper bounds   #################
  # ## Model effectively implies a restricted prior
  # ## Extract envelope information --> gs is grid size 
  # 
  # 
  # cbars=Env2$cbars
  # 
  # 
  #   
  # gs=nrow(Env2$cbars)
  # logP1=Env2$logP
  # 
  #   New_LL=c(1:gs)
  # New_logP2=c(1:gs)
  # 
  ## Step 6: Set dispersion bounds
  
  ## Set upper and lower bounds for dispersion.
  
  # upp=1/qgamma(c(1-max_disp_perc),shape2,rate3)
  # low=1/qgamma(c(max_disp_perc),shape2,rate3)

    #wt_upp=wt/rep(upp,length(y))
  #wt_low=wt/rep(low,length(y))
  
  #theta_upp=Inv_f3_gaussian(t(Env2$cbars), y, as.matrix(x2),as.matrix(mu2,ncol=1), as.matrix(P2), 
  #                          as.vector(alpha), as.vector(wt_upp))
  
  #theta_low=Inv_f3_gaussian(t(Env2$cbars), y, as.matrix(x2),as.matrix(mu2,ncol=1), as.matrix(P2), 
  #                          as.vector(alpha), as.vector(wt_low))
  
  
  #New_LL_upp=c(1:gs)
  #New_LL_low=c(1:gs)
  #New_LL_Slope_test=c(1:gs)
  #New_LL_Slope_test2=c(1:gs)
  #New_LL_Slope_test3=c(1:gs)
  #New_LL_Slope_diff=c(1:gs)
  
  # thetabars=Env2$thetabars
  # thetabar_const_base=thetabar_const(P2,cbars,Env2$thetabars)

    #thetabar_const_upp=thetabar_const(P2,cbars,theta_upp)
  #thetabar_const_low=thetabar_const(P2,cbars,theta_low)
  
  ###########################################################################################
  ### This part is currently used for QC - Move to inside function if needed at all
  
  #thetabar_const_matrix<-matrix(0,nrow=21,ncol=nrow(cbars))
  
  #for(i in 1:21){
  
  #  wt_temp=wt/rep(low+((i-1)/20)*(upp-low),length(y))
  #  theta_temp=Inv_f3_gaussian(t(Env2$cbars), y, as.matrix(x2),as.matrix(mu2,ncol=1), as.matrix(P2), 
  #                             as.vector(alpha), as.vector(wt_temp))
  
  #  thetabar_const_temp=thetabar_const(P2,cbars,theta_temp)
  #  thetabar_const_matrix[i,1:nrow(cbars)]=thetabar_const_temp
  
  
  #}
  
  ###########################################################################################
  
  ## This function EnvBuildLinBound() plays a critical role in your sampler: 
  ## it computes how each envelope face in the coefficient space shifts with respect to 
  ## the dispersion parameter \sigma^2. More precisely, it returns the derivative of each 
  ## face constant K_j with respect to \sigma^2, evaluated at the anchor point dispstar.
  
  
    
  # New_LL_Slope=EnvBuildLinBound(thetabars,cbars,y,x2,P2,alpha,dispstar)
  # 
  # 
  # ## Linear approximation of face constants
  # 
  # thetabar_const_upp_apprx=thetabar_const_base+(upp-dispstar)*New_LL_Slope
  # thetabar_const_low_apprx=thetabar_const_base+(low-dispstar)*New_LL_Slope
  # 
  # 
  # ## Summary statistics of slopes
  # 
  # m_New_LL_Slope=mean(New_LL_Slope)  
  # min_New_LL_Slope=min(New_LL_Slope)  
  # max_New_LL_Slope=max(New_LL_Slope)  
  # 
  # 
  # prob_factor<-c(1:gs)
  # min_log_accept<-c(1:gs)
  # 
  # ## Compute max face constants at endpoints
  # 
  # max_low=max(thetabar_const_low_apprx)
  # max_upp=max(thetabar_const_upp_apprx)
  # 
  # 
  # max_low=max_low+0*(max_upp-max_low)
  # 
  # 
  # 
  # ### Global upper lien over dispersion
  # 
  # max_low_mean=max_upp-m_New_LL_Slope*(upp-low)
  # old_slope=(max_upp-max_low)/(upp-low)
  # max_low=max_low_mean
  
  ## Per face mixture weight
  
  ## Adjust lopP2 for each face
  # 
  # 
  # for(j in 1:gs){
  #   cbars_temp=as.matrix(cbars[j,1:ncol(x)],ncol=1)
  #   New_logP2[j]=logP1[j]+0.5*t(cbars_temp)%*%cbars_temp
  #   
  #   prob_factor[j]=max(thetabar_const_upp_apprx[j]-max_upp,thetabar_const_low_apprx[j]-max_low)
  #   min_log_accept[j]=min(thetabar_const_upp_apprx[j]-max_upp,thetabar_const_low_apprx[j]-max_low)- prob_factor[j]  
  # }
  
  
  ## Normalize mixture weight
  
  # lg_prob_factor=prob_factor
  # prob_factor=exp(New_logP2+prob_factor)
  # 
  # prob_factor=prob_factor/sum(prob_factor)
  # 
  # 
  # ## Dispersion envelope geometry
  # new_slope=(max_upp-max_low)/(upp-low)
  # 
  # 
  # new_int=max_low-new_slope*low
  # b1=(upp-low)
  # c1=-log(upp/low)
  
  #dispstar= (-b1+ sqrt(b1^2-4*a1*c1))/(2*a1) # Point of tangency
  ## Outputs from this block
  
  ## 1) upp, low
  ## 2) log_P_diff
  ## 3) lm_log1,lm_log2
  ## 4) shape3
  ## 5) max_New_LL_UB
  ## 6) max_LL_log_disp
  
  ## Not currently from block: RSS_ML, rate2 
  
  #####################################  New Derivations ##################################
  
  # Verify that this calculation is correct
  
  
  ## adjust Gamma shape for dispersion proposal
  # 
  # dispstar=b1/(-c1)
  # 
  # shape_shift=New_LL_Slope*dispstar
  # 
  # 
  # lm_log2=new_slope*dispstar
  # lm_log1=new_int+new_slope*dispstar-new_slope*log(dispstar)
  # 
  # shape3=shape2-lm_log2
  # 
  # shape3_vector=shape2-shape_shift
  # max_LL_log_disp=lm_log1+lm_log2*log(upp) ## From above
  # 
  # 
  # 
  ## No longer used - can remove later 
  #log_P_diff_new=0*log_P_diff
  
  ########################################
  
  
  ## Store mixture weights in Envelope
  # Env3=Env2
  # Env3$PLSD=prob_factor
  
  #  print("Starting Candidate Sampling using sample function")
  #  cand1=sample(x=gs,size=n*413,replace=TRUE,prob=prob_factor)
  
  #print(cand1)
  #print(length(cand1))
  
  #  print("Finished Sampling using Sample function")
  
  ### Set constants used during sampling
  
  # gamma_list_new=list(shape3=shape3,rate2=rate2,disp_upper=upp,disp_lower=low)
  # 
  # UB_list_new=list(RSS_ML=RSS_ML,max_New_LL_UB=max_upp,
  #                  max_LL_log_disp=max_LL_log_disp,lm_log1=lm_log1,lm_log2=lm_log2, 
  #                  #log_P_diff=log_P_diff_new,
  #                  lg_prob_factor=lg_prob_factor,lmc1=new_int,lmc2=new_slope)
  # #                   ,cand=cand1)
  # 
  
  
  ## log_P_Diff should no longer be used in the same way!
  
  print(paste("Interactive status:", interactive()))
  
  ##  ptm <- proc.time()

  ## Run simulation
    
  # sim_temp=.rindep_norm_gamma_reg_std_V4_cpp (n=n, y=y, x=x2, mu=mu2, P=P2, alpha=alpha, wt,
  #                                             f2=f2, Envelope=Env3, 
  #                                             gamma_list=gamma_list_new,
  #                                             UB_list=UB_list_new,
  #                                             family="gaussian",link="identity", progbar =progbar)
  
  
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




EnvBuildLinBound<-function(thetabars,cbars,y,x2,P2,alpha,dispstar){
  
  gs=nrow(cbars)
  n_vars=ncol(cbars)
  
  New_LL_Slope_test2=c(1:gs)
  New_LL_Slope_test3=c(1:gs)
  
  for(j in 1:gs){
    
    cbars_temp=as.matrix(cbars[j,1:n_vars],ncol=1)
    thetabars_temp=as.matrix(thetabars[j,1:n_vars],ncol=1)
    New_LL_Slope_test2[j]=(-t(thetabars_temp)%*%P2+t(cbars_temp))%*%solve(t(x2)%*%x2+dispstar*P2)%*%cbars_temp
    
    H1=-solve(t(x2)%*%x2+dispstar*P2)%*%P2%*%solve(t(x2)%*%x2+dispstar*P2)
    New_LL_Slope_test3[j]=New_LL_Slope_test2[j]+(-t(thetabars_temp)%*%P2+t(cbars_temp))%*%H1%*%(t(x2)%*%(y-alpha)+dispstar*cbars_temp)
    
  }
  
  return(New_LL_Slope_test3)
  
}


thetabar_const<-function(P,cbars,thetabar){
  
  gs=nrow(cbars)
  thetaconst=c(1:gs)
  n_var=nrow(P)
  
  for(j in 1:gs){
    theta_temp=as.matrix(thetabar[j,1:n_var],ncol=1)
    cbars_temp=as.matrix(cbars[j,1:n_var],ncol=1)
    thetaconst[j]=-0.5*t(theta_temp)%*%P%*%theta_temp+t(cbars_temp)%*% theta_temp
    
  }
  
  return(thetaconst)
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
#'                   Gridtype = 2, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE)
#' @export 
#' @rdname simfuncs
#' @order 3


rNormal_Gamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                            Gridtype=2,
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
    if(!is.null(prior_list$Sigma)) Sigma=prior_list$Sigma
    if(!is.null(prior_list$P)) P=prior_list$P
    if(is.null(prior_list$P)) P=solve(prior_list$Sigma)
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
#' @usage rNormal_reg(n, y, x, prior_list, offset = NULL, weights = 1, family = gaussian(),
#'             Gridtype = 2, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE)
#' @export 
#' @rdname simfuncs
#' @order 2



rNormal_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
                      Gridtype=2,
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

