#' Simulation Functions for Bayesian Generalized Linear Models
#'
#' @description
#' Simulation functions provide a unified interface for generating posterior samples from Bayesian GLMs.
#' These functions are typically used within model fitting routines such as \code{\link{rglmb}} and \code{\link{rlmb}}, and
#' are also suitable for use in Block Gibbs sampling and other simulation-based inference techniques.
#'
#' @name simfuncs
#' @param object A fitted model object containing a \code{pfamily} component. The generic function \code{simfunction()} accesses the simulation metadata stored within such objects.
#' @param x An object of class \code{"simfunction"} or \code{"rGamma_reg"} to be printed.
#' @param n Number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y A vector of observations of length \code{m}.
#' @param x A design matrix of dimension \code{m * p}.
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
#'       \item{\code{coefficients}}{Matrix (\code{n * p}) of simulated regression coefficients, with column names from \code{x}.}
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
#'   \item{\code{rNormalGamma_reg()}}{A list with class \code{"rglmb"} containing:
#'     \describe{
#'       \item{\code{coefficients}}{Matrix (\code{n * p}) of simulated regression coefficients; row \code{i} equals \code{Btilde + IR \%*\% rnorm(p) * sqrt(dispersion[i])}. Column names are set to \code{colnames(x)}.}
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
#'       \item{\code{call}}{Matched call to \code{rNormalGamma_reg()}.}
#'     }
#'   }
#'
#'   \item{\code{rindepNormalGamma_reg()}}{A list with class \code{"rglmb"} containing:
#'     \describe{
#'       \item{\code{coefficients}}{Matrix (\code{n * p}) of simulated regression coefficients, back-transformed to the original scale; column names set to \code{colnames(x)}.}
#'       \item{\code{coef.mode}}{Vector with the conditional posterior mode used for envelope anchoring (from the Gaussian fit).}
#'       \item{\code{dispersion}}{Numeric vector of length \code{n} with simulated dispersion draws.}
#'       \item{\code{Prior}}{List with prior components: \code{mean} (prior mean \code{mu}), \code{Sigma} (prior covariance), \code{shape} and \code{rate} (Gamma prior for dispersion), \code{Precision} (\code{solve(Sigma)}).}
#'       \item{\code{family}}{The \code{gaussian()} family object.}
#'       \item{\code{prior.weights}}{Vector of prior weights used in the simulation.}
#'       \item{\code{y}}{Response vector.}
#'       \item{\code{x}}{Design matrix.}
#'       \item{\code{call}}{Matched call to \code{rindepNormalGamma_reg()}.}
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
#'       \item{\code{coefficients}}{A 1 * p matrix of assumed regression coefficients.}
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
#' @details The low-level simulation functions **\code{rNormal_reg()}**, **\code{rNormalGamma_reg()}**, 
#' **\code{rindepNormalGamma_reg()}**, and **\code{rGamma_reg()}** generate iid samples from posterior 
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
#' - **\code{rNormalGamma_reg()}**: Produces iid draws for regression coefficients and the
#'   dispersion parameter in models with Normal-Gamma priors and Gaussian likelihoods, where
#'   this is a conjugate prior distribution. Standard simulation procedures for gamma and
#'   multivariate normal distributions are utilized
#'   \insertCite{Raiffa1961,LindleySmith1972}{glmbayes}.
#'
#' - **\code{rindepNormalGamma_reg()}**: Produces iid draws for regression coefficients
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
#' \code{\link{rNormal_reg}}, \code{\link{rNormalGamma_reg}}, \code{\link{rGamma_reg}} for individual simulation functions.
#'
#' \code{\link{EnvelopeBuild}} for envelope construction methods used in likelihood subgradient sampling.
#'
#' @usage simfunction(object, ...)
#' @export
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



rGamma_reg <- function(
    n,
    y,
    x,
    prior_list,
    offset = NULL,
    weights = 1,
    family = gaussian(),
    Gridtype = 2,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE
) {
  call <- match.call()
  
  ## Argument renaming and prior
  wt    <- weights
  alpha <- offset
  
  ## Basic checks on y and x
  n1 <- length(y)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (nrow(x) != n1)
    stop("Number of rows in x must match length of y.")
  
  ## 1) Validate and normalize offset (alpha)
  if (is.null(alpha)) {
    alpha <- rep(0, n1)
  } else {
    if (!is.numeric(alpha))
      stop("offset (alpha) must be numeric if supplied.")
    if (length(alpha) == 1L) {
      alpha <- rep(alpha, n1)
    } else if (length(alpha) != n1) {
      stop("offset (alpha) must be scalar or have length equal to length(y).")
    }
  }
  
  ## 2) Validate and normalize weights (wt)
  if (!is.numeric(wt))
    stop("weights must be numeric.")
  
  if (length(wt) == 1L) {
    wt <- rep(wt, n1)
  } else if (length(wt) != n1) {
    stop("weights must be either a scalar or have length equal to length(y).")
  }
  
  if (any(wt < 0))
    stop("weights must be nonnegative.")
  
  b     <- prior_list$beta
  shape <- prior_list$shape
  rate  <- prior_list$rate
  
#  if (shape <= 1) {
#    stop("The inverse-gamma prior on dispersion has no finite mean unless shape > 1. 
#       Please supply a prior_list$shape value greater than 1.")
#  }

    
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
  
  ## Family handling
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  okfamilies <- c("gaussian", "Gamma")
  if (family$family %in% okfamilies) {
    if (family$family == "gaussian") oklinks <- c("identity")
    if (family$family == "Gamma")    oklinks <- c("log")
    if (!(family$link %in% oklinks)) {
      stop(gettextf(
        "link \"%s\" not available for selected family; available links are %s",
        family$link, paste(sQuote(oklinks), collapse = ", ")
      ), domain = NA)
    }
  } else {
    stop(gettextf(
      "family \"%s\" not available in glmbdisp; available families are %s",
      family$family, paste(sQuote(okfamilies), collapse = ", ")
    ), domain = NA)
  }
  
  ## ----------------------
  ## Gaussian case (updated)
  ## ----------------------
  if (family$family == "gaussian") {
    
    
    sim<-  .rGammaGaussian_cpp(
      n, y, x, b, wt, alpha, shape, rate,
      disp_lower, disp_upper, verbose = verbose
    )
    
    # Validate output
    if (!is.list(sim) || is.null(sim$dispersion) || is.null(sim$draws)) {
      stop("C++ .rGammaGaussian_cpp returned an invalid structure.")
    }
    
    out  <- sim$dispersion
    draws <- sim$draws
    
  }
  
  ## ----------------------
  ## Gamma case (MODIFIED)
  ## ----------------------
  if (family$family == "Gamma") {
    
    # Call C++ sampler — do NOT hide errors
    sim <- .rGammaGamma_cpp(
      n, y, x, b, wt, alpha, shape, rate,
      max_disp_perc, disp_lower, disp_upper,
      verbose = verbose
    )
    
    # Validate output
    if (!is.list(sim) || is.null(sim$dispersion) || is.null(sim$draws)) {
      stop("C++ .rGammaGamma_cpp returned an invalid structure.")
    }
    
    out  <- sim$dispersion
    draws <- sim$draws

  }
  
  outlist <- list(
    coefficients   = matrix(b, nrow = 1, ncol = length(b)),
    coef.mode      = NULL,
    dispersion     = out,
    Prior          = list(shape = shape, rate = rate),
    prior.weights  = wt,
    y              = y,
    x              = x,
    famfunc        = glmbfamfunc(family),
    iters          = draws,
    Envelope       = NULL
  )
  
  outlist$call <- call
  class(outlist) <- c(outlist$class, "rGamma_reg")
  
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





#' @family simfuncs 
#' @example inst/examples/Ex_rindep_norm_gamma_reg.R
#' @usage rindepNormalGamma_reg(n, y, x, prior_list, offset = NULL, weights = 1,
#'                              family = gaussian(), Gridtype = 2,n_envopt = NULL,
#'                               use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, 
#'                              progbar = TRUE)
#' @export 
#' @rdname simfuncs
#' @order 4



rindepNormalGamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
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
  

  
  # Reconstruct P from Sigma (needed by the temporary core)
  R <- chol(Sigma)
  P <- chol2inv(R)
  P <- 0.5 * (P + t(P))
  
  
  ##########################  BEGIN *.CPP  MIGRATION   #########################################################

  ## --- NEW: Normalize and validate inputs before calling Rcore ---
  
  # Coerce basic types
  y  <- as.numeric(y)
  x  <- as.matrix(x)
  mu <- as.numeric(mu)
  wt <- as.numeric(wt)
  
  # Normalize offset
  if (is.null(offset2)) offset2 <- rep(0, length(y))
  offset2 <- as.numeric(offset2)
  
  # Normalize weights
  if (length(wt) == 1L) wt <- rep(wt, length(y))
  stopifnot(length(wt) == length(y))
  
  # Dimension checks
  stopifnot(nrow(x) == length(y))
  stopifnot(length(mu) == ncol(x))
  
  # Reconstruct P from Sigma and enforce SPD
  R    <- chol(Sigma)
  Pinv <- chol2inv(R)
  P    <- 0.5 * (Pinv + t(Pinv))
  
  stopifnot(isSymmetric(P))
  
  tol <- 1e-6
  ev  <- eigen(P, symmetric = TRUE)$values
  stopifnot(all(ev >= -tol * abs(ev[1L])))
  
  # dispersion must be numeric scalar or NULL
  if (!is.null(dispersion)) {
    dispersion <- as.numeric(dispersion)
    stopifnot(length(dispersion) == 1L, is.finite(dispersion))
  }
  
  if (is.null(n_envopt)) n_envopt <- n
  n_envopt <- as.integer(n_envopt)


  
  core_out <- .rIndepNormalGammaReg_cpp(
    n,
    y,
    x,
    mu,
    P,
    offset2,
    wt,
    shape,
    rate,
    max_disp_perc,
    disp_lower,
    disp_upper,
    Gridtype,
    n_envopt,
    use_parallel,
    use_opencl,
    verbose,
    progbar
  )

  out        <- core_out$out
  betastar   <- core_out$betastar
  disp_out   <- core_out$disp_out
  iters_out  <- core_out$iters_out
  weight_out <- core_out$weight_out
  low        <- core_out$low
  upp        <- core_out$upp
  
  famfunc=glmbfamfunc(gaussian())  
  f1=famfunc$f1
  
  R <- chol(Sigma)
  Prec <- chol2inv(R)
  Prec <- 0.5 * (Prec + t(Prec))   # enforce symmetry
  

  outlist=list(
    coefficients=t(out), 
    coef.mode=betastar,  ## For now, use the conditional mode (not universal)
    dispersion=disp_out,
    ## For now, name items in list like this-eventually make format/names
    ## consistent with true prior (current names needed by summary function)
    Prior=list(mean=mu,Sigma=Sigma,shape=shape,rate=rate,Precision=Prec), 
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


#' @family simfuncs 
#' @example inst/examples/Ex_rnorm_gamma_reg.R
#' @usage rNormalGamma_reg(n, y, x, prior_list, offset = NULL, weights = 1, family = gaussian(),
#'                   Gridtype = 2,n_envopt = NULL, 
#'                   use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE)
#' @export 
#' @rdname simfuncs
#' @order 3


rNormalGamma_reg<-function(n,y,x,prior_list,offset=NULL,weights=1,family=gaussian(),
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
      R <- chol(prior_list$Sigma)
      P <- chol2inv(R)
      P <- 0.5 * (P + t(P))   # enforce symmetry
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
  if(any(wt < 0)) stop("weights must be non-negative")
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
  

  ## Add Call to new function here or after famfunc / f1 is set (not sure if f1 is actually used)
  
  sim <- .rNormalGammaReg_cpp(
    n          = n,
    y          = y,
    x          = x,
    mu         = mu,
    P          = P,
    offset     = offset2,
    wt         = wt,
    shape      = shape,
    rate       = rate,
    max_disp_perc = NULL,
    disp_lower = NULL,
    disp_upper = NULL,
    verbose    = verbose
  )
  
  # Extract components returned by C++
  out1        <- sim$coefficients
  Btilde      <- sim$coef.mode
  dispersion  <- sim$dispersion
  fit         <- sim$fit
  famfunc     <- sim$famfunc
  draws       <- sim$iters
  
  # Build final outlist (mirrors original structure)
  outlist <- list(
    coefficients   = out1,
    coef.mode      = Btilde,
    dispersion     = dispersion,
    offset         = offset,
    Prior          = list(mean = as.numeric(mu), Precision = P),
    prior.weights  = wt,
    y              = y,
    x              = x,
    fit            = fit,
    famfunc        = famfunc,
    iters          = draws,
    Envelope       = NULL
  )
  
  colnames(outlist$coefficients) <- colnames(x)
  
  outlist$call <- match.call()
  class(outlist) <- c(outlist$class, "rglmb")
  
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

    if (is.null(prior_list$P)) {
      R <- chol(prior_list$Sigma)
      Pinv <- chol2inv(R)
      P <- 0.5 * (Pinv + t(Pinv))   # enforce symmetry
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
    outlist<-.rNormalReg_cpp(n=n,y=y,x=x,mu=mu,P=P,offset=offset2,wt=wt,dispersion=dispersion,
                            ##                      famfunc=famfunc,f1=f1,
                            f2=f2,f3=f3,start=mu)
    class(outlist$fit)="lm"
    
  }
  else{
    if(is.null(dispersion)){dispersion2=1}
    else{dispersion2=dispersion}
    
    #  stop("Inputs to function above")
    outlist<-.rNormalGLM_cpp(n=n,y=y,x=x,mu=mu,P=P,offset=offset2,wt=wt,
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
      
      outlist<-.rNormalGLM_cpp(n=n,y=y,x=x,mu=mu,P=P,offset=offset2,
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



# Helpers --------------------------------------------------------------------

################################## Utility functions used by the above  #################



# Helper: log(exp(b) - exp(a)) for a < b, in a numerically stable way
logdiffexp <- function(a, b) {
  # assumes a < b (strict), returns log(exp(b) - exp(a))
  # exp(b) - exp(a) = exp(b) * (1 - exp(a - b))
  # so log(...) = b + log(1 - exp(a - b))
  b + log(-expm1(a - b))
}



# .rnnorm_reg_std_cpp -->rNormalGLM_std_cpp
# .rnorm_reg_cpp --> rNormalReg_cpp
# .rindep_norm_gamma_reg_cpp --> rIndepNormalGammaReg_cpp
# .rindep_norm_gamma_reg_std_cpp -->rIndepNormalGammaReg_std_cpp
# .rindep_norm_gamma_reg_std_parallel_cpp --> rIndepNormalGammaReg_std_parallel_cpp
# .rnnorm_reg_cpp --> .rNormalGLM_cpp
