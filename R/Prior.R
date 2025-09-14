#' Setup Prior Objects
#' 
#' Helper function to facilitate the Setup of Prior Distributions for glm models.
#' @param na.action how \code{NAs} are treated. The default is first, any \code{\link{na.action}} attribute of 
#' data, second a \code{na.action} setting of \link{options}, and third \code{na.fail} if that is unset. 
#' The \code{factory-fresh} default is \code{na.omit}. Another possible value is \code{NULL}.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param pwt Weight on the prior relative to the likelihood function at the the maximum likelihood estimate. If n_prior is provided, it is calculated as pwt=n_prior/(n_prior+n_likelihood) where n_likelihood is the number of observations for the likelihood function.
#' @param n_prior Optional argument with number of prior observations (either a scalar or a vector). When provided, this is used together with the number of likelihood observations to compute pwt. If not provided but pwt is a scalar, it is computed as n_prior=n_likelihood*(pwt/(1-pwt)).
#' @param sd Optional vector argument with the prior standard deviations for the coefficients
#' @param intercept_source Specifies the method through which the prior mean for the intercept term is set. Options are based on the null intercept only model (null_model) or full_models. The default is the null model which is safer if variables are not centered. 
#' @param effects_source Specifies the method through which the prior means for the effects terms are set. Options are null_effects (prior means set to zero) or full_model (effect means set to match maximum likelihood estimates).  
#' @param mu Optional vector argument with the prior means for the coefficients
#' @inheritParams stats::glm
#' @inheritParams stats::model.frame
#' @details
#' `Prior_Setup()` initializes a structured set of prior parameters for generalized linear models (GLMs), supporting both Gaussian and non-Gaussian families.
#'  It is designed to provide the full set of inputs required for multiple prior specifications (referred to as "pfamilies"), including conjugate normal priors, 
#'  Normal-Gamma priors, and independent normal-gamma priors.
#'
#' The function returns a list containing:
#' * `mu`: the prior mean vector
#' * `Sigma`: the prior variance-covariance matrix
#' * `dispersion`: the estimated dispersion (for Gaussian models)
#' * `shape`: the shape parameter for Normal-Gamma priors (if applicable)
#' * `rate`: the rate parameter for Normal-Gamma priors (if applicable)
#' * `model`: the model frame used to construct the design matrix
#' * `x`: the design matrix
#' * `call`: the matched call to `Prior_Setup()`
#' * `PriorSettings`: a list of metadata including:
#'   - `pwt`: prior weight (scalar or vector)
#'   - `n_prior`: effective prior sample size
#'   - `n_likelihood`: effective likelihood sample size
#'   - `intercept_source`: method used to set the prior mean for the intercept
#'   - `effects_source`: method used to set the prior mean for the effects
#'
#' ### Inputs to the function
#'
#' The inputs to `Prior_Setup()` fall into three conceptual categories:
#'
#' **1. Model specification**
#' * `formula`: defines the structure of the GLM, including response and predictors.
#' * `family`: specifies the error distribution and link function (e.g., `gaussian()`, `binomial()`).
#' * `data`: optional data frame used to evaluate the formula and extract variables.
#'
#' **2. Prior variance-covariance specification**
#' * `pwt`: prior weight relative to the likelihood. If scalar, used to compute Zellner's g-prior.
#' * `n_prior`: optional effective prior sample size. If provided, overrides `pwt`.
#' * `sd`: optional vector of prior standard deviations. If provided, used to compute `pwt`.
#'
#' **3. Prior mean specification**
#' * `intercept_source`: method for setting the prior mean of the intercept (`"null_model"` or `"full_model"`).
#' * `effects_source`: method for setting the prior mean of the effects (`"null_effects"` or `"full_model"`).
#' * `mu`: optional user-specified prior mean vector. Overrides other centering logic if provided.
#'
#' ### Mathematical structure and interpretation
#'
#' When `pwt` is a scalar, the prior covariance is scaled from the likelihood covariance using a Zellner g-prior:
#' \deqn{\Sigma_{Prior} = \frac{1 - pwt}{pwt} \cdot V_{MLE}}
#' where \eqn{V_{MLE}} is the variance–covariance matrix of the maximum likelihood estimator.
#'
#' If the likelihood covariance is not full rank, the function aborts with an error, as a g-prior cannot be constructed.
#'  While Bayesian models can still be estimated in such cases, users should proceed with caution.
#'  
#' The corresponding prior precision is:
#' \deqn{P_{Prior} = \frac{pwt}{1 - pwt} \cdot P_{MLE}}
#' and the posterior mean (in the Gaussian case) simplifies to:
#' \deqn{\mu_{Post} = (1 - pwt) \cdot \widehat{\beta}_{MLE} + pwt \cdot \mu_{Prior}}
#'
#' Default prior centering reflects classical reference structures:
#' \itemize{
#'   \item The intercept prior is centered at the null model estimate (`intercept_source = "null_model"`), consistent with how R\eqn{^2}, F-tests, and t-statistics are defined.
#'   \item The effect priors are centered at zero (`effects_source = "null_effects"`), supporting intuitive posterior summaries and direct interpretation of tail probabilities.
#' }
#' 
#' These defaults ensure that posterior means remain close to classical estimates when `pwt` is small, 
#' while allowing additional flexible prior specifications through `mu`, `sd`, and `n_prior`. 
#' 
#' If `n_prior` is not provided, it is derived from `pwt` and the effective likelihood sample size `n_likelihood` using:
#' \deqn{n_{\mathrm{prior}} = \frac{pwt}{1 - pwt} \cdot n_{\mathrm{likelihood}}}
#' where `n_likelihood` is the number of observations in the model.
#'
#' Likewise, if `n_prior` is provided, `pwt` is computed as:
#' \deqn{pwt = \frac{n_{\mathrm{prior}}}{n_{\mathrm{prior}} + n_{\mathrm{likelihood}}}}
#'
#'#' When applicable, `Prior_Setup()` computes the shape and rate parameters for a Gamma prior on the residual precision (inverse variance), used in 
#'compound prior families such as `dNormal_Gamma()` , `dIndependent_Normal_Gamma()`, and `dGamma()`. These are derived from the effective prior sample size and the estimated dispersion:
#'
#' \deqn{\text{shape} = \frac{n_{\mathrm{prior}}}{2}}
#' \deqn{\text{rate} = \text{shape} \cdot \text{dispersion} = \frac{n_{\mathrm{prior}}}{2} \cdot \text{dispersion}}
#'    
#' where RSS is the residual sum of squares from the likelihood model.
#' 
#' The posterior shape and rate parameters under a Normal-Gamma model are then:
#' \deqn{\text{shape}_{\mathrm{post}} = \text{shape} + \frac{n_{\mathrm{likelihood}}}{2} = \frac{n_{\mathrm{prior}} + n_{\mathrm{likelihood}}}{2}}
#' \deqn{\text{rate}_{\mathrm{post}} = \text{rate} + \frac{1}{2} \cdot \text{RSS} =  \frac{n_{\mathrm{prior}}+n_{\mathrm{likelihood}} - k}{2}  \cdot dispersion}
#' 
#' This structure allows the prior to contribute pseudo-observations to the residual precision estimate, enabling adaptive shrinkage and hierarchical 
#' regularization — especially valuable in small-sample or high-dimensional settings.
#' especially in small-sample or high-dimensional settings.
#'
#'
#' @return A list with items related to the prior.
#' \item{mu}{A prior mean vector}
#' \item{Sigma}{A prior variance-covariance matrix}
#' \item{dispersion}{Empirical bayes estimate for the dispersion (gaussian model only)}
#' \item{shape}{Derived prior shape parameter (gaussian model only). Defaults to n_prior/2 where n_prior is derived from pwt if not provided}
#' \item{rate}{Derived prior rate parameter (gaussian model only). Defaults to (n_prior*dispersion)/2 where n_prior is derived from pwt if not provided}
#' \item{model}{The model frame from \code{object} if it exists}
#' \item{x}{The design matrix from \code{object} if it exists}
#' \item{PriorSettings}{A list containing prior configuration details}
#' @family prior
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_Prior_Setup.R
#' @export
#' @rdname Prior_Setup
#' @order 1

## Note arguments outside of first two are currently not used

Prior_Setup <- function(
    formula,
    family      = gaussian(),
    data=NULL,
    weights=NULL,
    subset=NULL,
    na.action   = na.fail,
    offset=NULL,
    contrasts   = NULL,
    pwt         = 0.01,
    n_prior     = NULL,
    sd          = NULL,
    intercept_source = c("null_model", "full_model"),
    effects_source   = c("null_effects",  "full_model"),
    mu          = NULL,  ...
  ) 
  
  {

  call <- match.call()  
  intercept_source <- match.arg(intercept_source)
  effects_source <- match.arg(effects_source)
  
  #mf<-model.frame(formula,data,subset=subset,na.action=na.action,
  #                drop.unused.levels=drop.unused.levels,xlev=xlev)
  
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  if (missing(data))   data <- environment(formula)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
##  mf<-model.frame(formula,data)
  
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  
  if (family$family == "binomial"
      && is.numeric(Y) && is.vector(Y)
      && all(Y >= 0 & Y <= 1)    # use <= instead of ≤
      && is.null(weights)) {
    warning(
      "You supplied a proportion response (0 <= y <= 1) to a binomial family\n",
      "without `weights`. Each case will be treated as a single trial (n=1).\n",
      "If you meant to model counts, either use `cbind(success, failure)`\n",
      "or supply `weights =` the number of trials."
    )
  }
  
  
##  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(, NROW(Y), 0L)
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, ...) else matrix(, NROW(Y), 0L)
  
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) 
    stop("negative weights not allowed")
  
  
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                    length(offset), NROW(Y)), domain = NA)
  }
  
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  
  
  x<-X  
##  x<-model.matrix(formula,mf)
  
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
#  mu=matrix(0,nrow=nvar,ncol=1, 
#            dimnames = list(var_names, "mu")
#            )

  mu_internal <- matrix(0, nrow = nvar, ncol = 1, dimnames = list(var_names, "mu"))
 
  
  
##  glm_full=glm(formula, family = family,data=data)
  
##  glm_full <- glm(
##    formula = formula,
##    family = family,
##    data = data,
##    data = mf,
##    weights = weights,
##    offset = offset,
##    subset = subset,
##    na.action = na.action,
##    ...
##  )
  
  glm_full <- glm.fit(
    x       = X,
    y       = Y,
    weights = weights,
    offset  = offset,
    family  = family
    ,control = glm.control(...)
  )
  
  #glm.fit(x, y, weights = rep.int(1, nobs),
  #        start = NULL, etastart = NULL, mustart = NULL,
  #        offset = rep.int(0, nobs), family = gaussian(),
  #        control = list(), intercept = TRUE, singular.ok = TRUE)
  
  
  glm_full$call      <- call
  glm_full$formula   <- formula
  glm_full$terms     <- mt
  glm_full$data      <- mf
  glm_full$offset    <- offset
  glm_full$contrasts <- attr(X, "contrasts")
  glm_full$xlevels   <- .getXlevels(mt, mf)
  class(glm_full)    <- c("glm", "lm")
  
##    glm_full=glm(formula, family = family,data=data)
  
  V0 <- vcov(glm_full)
  
  glm_summary=summary(glm_full)
  
  n_likelihood <- glm_summary$df.residual + glm_summary$df[1]  # residual df + model rank

  # If sd is provided, use it to compute pwt
if (!is.null(sd)) {
  if (!is.numeric(sd) || any(is.na(sd))) {
    stop("sd must be a numeric vector with no missing values.")
  }
  if (length(sd) != nvar) {
    stop("Length of sd must match number of coefficients (", nvar, ").")
  }

  # Compute pwt from sd and V0
  V0_diag <- diag(V0)
  if (any(V0_diag <= 0)) {
    stop("Diagonal entries of V0 must be positive to compute pwt from sd.")
  }

  pwt <- V0_diag / (V0_diag + sd^2)
  message("Computed pwt from user-specified prior standard deviations (sd).")
}
  
    
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
  
  ## Compute shape and rate if n_prior is not null and n_prior is scalar
  if (!is.null(n_prior)&& length(n_prior) == 1L&&!is.null(dispersion) ) {
    shape<-n_prior/2
    rate<- dispersion*shape
  }
  
  else{
    shape<-NULL
    rate <-NULL    
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
  
  
  
##  if(var_names[1]=='(Intercept)'){
    ##lm_out=lm(formula,data=mf,y=TRUE)
    ##y=lm_out$y
    ##mu[1,1]=mean(y)
    
##    f<-formula
    
##    lhs<-f[[2]]
##    intercept_only<-as.formula(paste(deparse(lhs),"~1") ,env=environment(f))
    

##    glm_null=update(glm_full,formula=intercept_only)
    
##    chosen_int <- switch(intercept_source,
##                         null_model = coef(glm_null)[1],
##                         full_model = coef(glm_full)[1])
##    mu_internal[1, 1] <- chosen_int
  
##  } 
  
  if (var_names[1] == "(Intercept)") {
    # build 1-column design matrix for intercept only
    X0 <- matrix(1, nrow = NROW(Y), ncol = 1,
                 dimnames = list(NULL, "(Intercept)"))
    
    # fit intercept-only model via glm.fit()
    fit0 <- glm.fit(
      x       = X0,
      y       = Y,
      weights = weights,
      offset  = offset,
      family  = family,
      control = glm.control(...)
    )
    
    # pick the intercept from null or full model
    chosen_int <- switch(
      intercept_source,
      null_model = fit0$coefficients[1],
      full_model = glm_full$coefficients[1]
    )
    
    mu_internal[1, 1] <- chosen_int
  }
  
  
  
  
  # 5) effects prior means
  if (nvar > 1) {
    effect_names <- var_names[-1]
    if (effects_source == "full_model") {
      coefs <- coef(glm_full)[effect_names]
      mu_internal[effect_names, 1] <- coefs
    }
    # else null_effects leaves mu[...] as zero
  }
  

  # Validate user-supplied mu if provided
  if (!is.null(mu)) {
    if (!is.numeric(mu)) {
      stop("mu must be numeric.")
    }
    if (is.vector(mu)) {
      if (length(mu) != nvar) {
        stop("Length of mu vector must match number of coefficients (", nvar, ").")
      }
      mu <- matrix(mu, ncol = 1, dimnames = list(var_names, "mu"))
    } else if (is.matrix(mu)) {
      if (!all(dim(mu) == c(nvar, 1))) {
        stop("mu matrix must have dimensions [", nvar, ", 1].")
      }
      rownames(mu) <- var_names
      colnames(mu) <- "mu"
    } else {
      stop("mu must be either a numeric vector or a matrix.")
    }
    message("Using user-specified prior mean vector (mu).")
  } else {
    mu <- mu_internal
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
    shape=shape,
    rate=rate,
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
  invisible(print(out))
  
  if (nvar <= 10) {
    cat("\nPrior Correlation Matrix\n")
    invisible(print(round(prior_cor, 4)))
  }
  
  # Optional: print dispersion
  if (!is.null(x$dispersion)) {
    cat("\nConditional Dispersion (Gaussian family): ", round(x$dispersion, 4), "\n\n")
  }
  
  if (!is.null(x$shape) && !is.null(x$rate)) {
    cat("Gamma Prior on Residual Precision:\n")
    cat("  shape =", round(x$shape, 4), "\n")
    cat("  rate  =", round(x$rate, 4), "\n")
    cat("  Expected precision (inverse variance) =", round(x$shape / x$rate, 6), "which implies 1/Expected precision  =", round(x$rate/x$shape , 6), "\n\n"  )
    cat("  Applicable to gaussian models with compound pfamilies (e.g., dNormal_Gamma, dIndependent_Normal_Gamma),\n")
    cat("  as well as for Gamma regression, quasipoisson, and quasibinomial models.\n\n")
  }
  
  if (is.null(x$shape) && is.null(x$rate) && !is.null(x$dispersion)) {
    cat("Note: Gaussian family detected, but shape/rate parameters were not computed.\n")
    cat("This may occur if n_prior is not scalar.\n\n")
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




