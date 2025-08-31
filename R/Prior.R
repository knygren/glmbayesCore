#' Setup Prior Objects
#' 
#' Sets up the structure for the Prior mean and Variance Matrices using information from a classical model.
#' @param na.action how \code{NAs} are treated. The default is first, any \code{\link{na.action}} attribute of 
#' data, second a \code{na.action} setting of \link{options}, and third \code{na.fail} if that is unset. 
#' The \code{factory-fresh} default is \code{na.omit}. Another possible value is \code{NULL}.
#' @param family a description of the error distribution and linke function to be used in the model.
#' @param pwt Weight on the prior relative to the likelihood function at the the maximum likelihood estimate.
#' @param n_prior Optional argument with number of prior observations (either a scalar or a vector). When provided, this used together with the number of likelihood observations to compute the pwt.
#' @param intercept_source Specifies the method through which the prior mean for the intercept term is set. Options are based on the null intercept only model (null_model) or full_models. The default is the null model which is safer if variables are not centered. 
#' @param effects_source Specifies the method through which the prior means for the effects terms are set. Options are null_effects (prior means set to zero) or full_model (effect means set to match maximum likelihood estimates).  
#' @inheritParams stats::model.frame
#' @return A list with items related to the prior.
#' \item{mu}{A prior mean vector}
#' \item{Sigma}{A prior variance-covariance matrix}
#' \item{dispersion}{Empirical bayes estimate for the dispersion (gaussian model only)}
#' \item{model}{The model frame from \code{object} if it exists}
#' \item{x}{The design matrix from \code{object} if it exists}
#' @family prior


#' @example inst/examples/Ex_Prior_Setup.R
#' @export
#' @rdname Prior_Setup
#' @order 1

## Note arguments outside of first two are currently not used

Prior_Setup<-function(formula,data=NULL,family=gaussian(),pwt=0.01 ,
                      n_prior=NULL,
                      intercept_source = c("null_model", "full_model"),
                      effects_source = c("null_effects", "full_model"),
                      subset = NULL, na.action = na.fail, 
                         drop.unused.levels = FALSE, xlev = NULL, ...){

  call <- match.call()  
  #mf<-model.frame(formula,data,subset=subset,na.action=na.action,
  #                drop.unused.levels=drop.unused.levels,xlev=xlev)
  
  intercept_source <- match.arg(intercept_source)
  effects_source <- match.arg(effects_source)
  
  
  
  
  
  mf<-model.frame(formula,data)
  x<-model.matrix(formula,mf)
  
  nvar=ncol(x)
  
  
  
  ## Make sure the *columns* of x are named correctly:
  
  
  ## validate pwt  
  if (!is.numeric(pwt) || any(is.na(pwt))) {  
    stop("pwt must be numeric and non NA, either length 1 or length ", nvar)  
  }  
  if (! (length(pwt) %in% c(1, nvar)) ) {  
    stop("pwt must have length 1 or length(coef) = ", nvar,  
         "; you supplied length ", length(pwt))  
  }  
  if (any(pwt <= 0 | pwt >= 1)) {  
    stop("All elements of pwt must lie strictly between 0 and 1; you supplied:",  
         paste0(round(pwt, 3), collapse = ", "))  
  }  
  
  
  
  var_names <- colnames(x)
  colnames(x) <- var_names
  #nvar=length(object$coefficients)
  mu=matrix(0,nrow=nvar,ncol=1, 
            dimnames = list(var_names, "mu")
            )

  glm_full=glm(formula, family = family,data=data)
  V0 <- vcov(glm_full)
  
  glm_summary=summary(glm_full)
  
  n_likelihood <- glm_summary$df.residual + glm_summary$df[1]  # residual df + model rank
  
  # Override pwt if n_prior is provided
  if (!is.null(n_prior)) {
    if (!is.numeric(n_prior) || length(n_prior) != 1 || n_prior <= 0) {
      stop("n_prior must be a single positive numeric value")
    }
    pwt <- n_prior / (n_prior + n_likelihood)
    message("Computed pwt = ", round(pwt, 4),
            " from n_prior = ", n_prior,
            " and n_likelihood = ", n_likelihood)
  }
  
  # Compute n_prior if not supplied and pwt is scalar
  if (is.null(n_prior) && length(pwt) == 1L) {
    n_prior <- (pwt/(1-pwt)) * n_likelihood
  #  message("Computed n_prior = ", round(n_prior, 4),
  #          " from pwt = ", round(pwt, 4),
  #          " and n_likelihood = ", n_likelihood)
  }
  

  
  ## conditional dispersion
  if (family$family == "gaussian") {
    ## summary(glm) gives you the MLE σ̂² via $dispersion
    dispersion <- summary(glm_full)$dispersion
  } else {
    dispersion <- NULL
  }
  
  
  

    if (!is.matrix(V0) || nrow(V0) != ncol(V0)) {
    stop("vcov(glm_full) (V0) must be a square matrix.")
  }
  if (anyNA(V0)) {
    stop("vcov(glm_full) (V0) contains missing values.")
  }
  
  # 2. symmetry (up to numerical tolerance)
  if (!isSymmetric(V0, tol = sqrt(.Machine$double.eps))) {
    stop("vcov(glm_full) (V0) is not symmetric.")
  }
  
  # 3. positive-definiteness via Cholesky
  pd_try <- try(chol(V0), silent = TRUE)
  if (inherits(pd_try, "try-error")) {
    stop(
      "Variance-covariance matrix V0 is not positive-definite.\n",
      "This usually means the classical GLM is rank-deficient."
    )
  }
  
  
  
  if(var_names[1]=='(Intercept)'){
    ##lm_out=lm(formula,data=mf,y=TRUE)
    ##y=lm_out$y
    ##mu[1,1]=mean(y)
    
    f<-formula
    
    lhs<-f[[2]]
    intercept_only<-as.formula(paste(deparse(lhs),"~1") ,env=environment(f))
    

    glm_null=update(glm_full,formula=intercept_only)
    
    chosen_int <- switch(intercept_source,
                         null_model = coef(glm_null)[1],
                         full_model = coef(glm_full)[1])
    mu[1, 1] <- chosen_int
    
    
    
    
  } 
  
  
  
  
  # 5) effects prior means
  if (nvar > 1) {
    effect_names <- var_names[-1]
    if (effects_source == "full_model") {
      coefs <- coef(glm_full)[effect_names]
      mu[effect_names, 1] <- coefs
    }
    # else null_effects leaves mu[...] as zero
  }
  
  
  
#  Sigma=as.matrix(diag(nvar))
 
#  Sigma=(1-pwt)/pwt*V0
   
  ## build prior covariance  
  if (length(pwt) == 1L) {  
    ## full matrix prior  
    Sigma <- ((1 - pwt) / pwt) * V0  
  }
  else {  
    scale_vec <- sqrt((1 - pwt) / pwt)
    scale_mat <- outer(scale_vec, scale_vec)
    Sigma <- V0 * scale_mat
  }  
  
  
  
  rownames(mu)=var_names
  colnames(mu)=c("mu")
  rownames(Sigma)=var_names
  colnames(Sigma)=var_names
  
  prior_list <- list(
    mu = mu,
    Sigma = Sigma,
    dispersion = dispersion,
    model = mf,
    x = x,
    call=call,
    PriorSettings = list(
      pwt = pwt,
      n_prior = n_prior,
      intercept_source = intercept_source,
      effects_source = effects_source,
      n_likelihood=n_likelihood
    )
  )
  
  class(prior_list) <- "PriorSetup"
  return(prior_list)
  
}

#' @export
#' @method print PriorSetup

print.PriorSetup <- function(x, ...) {

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    
  settings <- x$PriorSettings
  
  if (!is.null(settings$pwt) && length(settings$pwt) == 1L) {
    g <- (1 - settings$pwt) / settings$pwt
    cat("Setting up a Zellner g-type prior: \n")
    cat("  pwt =", round(settings$pwt, 4), "\n")
    cat("  g   = (1 - pwt)/pwt =", round(g, 4), "\n\n")
  }
  
  if (!is.null(settings$n_prior) && !is.null(settings$n_likelihood)) {
    cat("Note: n_prior was computed as (pwt / (1 - pwt)) * n_likelihood: \n")
    cat("  n_prior      =", round(settings$n_prior, 4), "\n")
    cat("  n_likelihood =", round(settings$n_likelihood, 4), "\n\n")
  }
  

  if (!is.null(settings$pwt) && length(settings$pwt) > 1L) {
    cat("Note: Differential prior weights (pwt) were specified per coefficient.\n\n")
  }
  
  
  cat("Prior Setup Summary\n")
  cat("====================\n\n")
  
  # Check for Zellner g-prior structure
  Sigma <- x$Sigma
  mu <- x$mu
  var_names <- rownames(mu)
  nvar <- length(var_names)
  
  # Extract diagonal SDs
  prior_sd <- sqrt(diag(Sigma))
  
  # Always compute prior correlation matrix
  prior_cor <- cov2cor(Sigma)
  
  # Extract pwt vector for display
  if (!is.null(settings$pwt)) {
    if (length(settings$pwt) == 1L) {
      pwt_vec <- rep(settings$pwt, nvar)
    } else if (length(settings$pwt) == nvar) {
      pwt_vec <- settings$pwt
    } else {
      warning("Length of pwt does not match number of variables; skipping pwt column.")
      pwt_vec <- rep(NA_real_, nvar)
    }
  } else {
    pwt_vec <- rep(NA_real_, nvar)
  }
  

  # Compute 95% intervals
  z <- qnorm(0.975)
  lower <- mu[, 1] - z * prior_sd
  upper <- mu[, 1] + z * prior_sd
  

  # Build output table

  out <- data.frame(
    Prior.Mean = round(mu[, 1], 6),
    Prior.SD   = round(prior_sd, 6),
    CI.Lower   = round(lower, 6),
    CI.Upper   = round(upper, 6),
    pwt        = round(pwt_vec, 6)
  )
  
  # Print table
  cat("Prior Estimates with 95% Confidence Intervals\n")
  print(out)
  
  if (nvar <= 10) {
    cat("\nPrior Correlation Matrix\n")
    print(round(prior_cor, 4))
  }
  
  # Optional: print dispersion
  if (!is.null(x$dispersion)) {
    cat("\nConditional Dispersion (Gaussian family): ", round(x$dispersion, 4), "\n")
  }
  
  invisible(x)
}



#' Checks for Prior-data conflicts
#'
#' Checks if the credible intervals for the prior overlaps with the implied confidence intervals from 
#' the classical model that comes from a call to the glm function 
#' @param level the confidence level at which the Prior-data conflict should be checked.
#' @inheritParams glmb
#' @return A vector where each item provided the ratio of the absolue value for the difference between the 
#' prior and maximum likelihood estimate divided by the length of the sum of half of the two intervals 
#' (where normality is assumed)
#' @family prior
#' @example inst/examples/Ex_Prior_Check.R
#' @export
#' @rdname Prior_Check
#' @order 1

Prior_Check<-function(formula,family,pfamily,level=0.95,data=NULL, weights, subset,na.action, 
                      start = NULL, etastart, mustart, offset ,control = list(...) , model = TRUE, 
                      method = "glm.fit",x = FALSE, y = TRUE, contrasts = NULL, ...){
  
  pf=pfamily
  prior_list=pfamily$prior_list
  
  ## For now, the below is really only correct for the dNormal pfamily
  
  mu=prior_list$mu
  Sigma=prior_list$Sigma
  
  
  object=glm(formula=formula,family=family,data=data)
  
  Like_est=object$coefficients
  Like_std=summary(object)$coefficients[,2]
  
  if(is.null(mu)){
    print("No Prior mean vector provided. Variables with needed Priors are:")
    print(names(Like_est))
    names(Like_est)    
    
  }
  
  
  if(level<=0.5) stop("level must be greater than 0.5")
  
  Sigma=as.matrix(Sigma)
  Prior_std=sqrt(diag(Sigma))
  
  print("Variables in the Model Are:")
  print(names(Like_est))
  std_dev_sum=qnorm(level)*(Prior_std+Like_std)
  
  abs_ratio=matrix(rep(0,length(Like_est),nrow=length(Like_est),ncol=1))
  abs_diff=abs(mu-Like_est)
  abs_ratio[1:length(Like_est),]=abs_diff/std_dev_sum
  
  rownames(abs_ratio)=names(Like_est)
  colnames(abs_ratio)=c("abs_ratio")
  max_abs_ratio=max(abs_ratio)
  
  if(max_abs_ratio>1) {
    print("At least one of the maximum likelihood estimates appears to be inconsistent with the prior")
  }
  
  else{
    print("The maximum likelihood estimates for all coefficients appear to be roughly consistent with the prior.")
  }
  return(abs_ratio)
  
}




