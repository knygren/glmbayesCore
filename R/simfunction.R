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
#'   \item{\code{rNormal_Gamma_reg()}}{A list with class \code{"rglmb"} containing:
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
#'       \item{\code{call}}{Matched call to \code{rNormal_Gamma_reg()}.}
#'     }
#'   }
#'
#'   \item{\code{rindependent_norm_gamma_reg()}}{A list with class \code{"rglmb"} containing:
#'     \describe{
#'       \item{\code{coefficients}}{Matrix (\code{n * p}) of simulated regression coefficients, back-transformed to the original scale; column names set to \code{colnames(x)}.}
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
#' @keywords internal
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



# Helper: log(exp(b) - exp(a)) for a < b, in a numerically stable way
logdiffexp <- function(a, b) {
  # assumes a < b (strict), returns log(exp(b) - exp(a))
  # exp(b) - exp(a) = exp(b) * (1 - exp(a - b))
  # so log(...) = b + log(1 - exp(a - b))
  b + log(-expm1(a - b))
}

# Constrained truncated gamma sampler on precision scale (numerically robust)
ctrgamma <- function(nn, shape, rate, lower_prec = NULL, upper_prec = NULL) {
  # ---- Case 0: no truncation ----
  if (is.null(lower_prec) && is.null(upper_prec)) {
    return(stats::rgamma(nn, shape = shape, rate = rate))
  }
  
  # ---- Extract numeric bounds (support-aware) ----
  # Gamma support is [0, infinity), so clamp at 0 if lower is NULL
  L <- if (!is.null(lower_prec)) lower_prec else 0
  U <- if (!is.null(upper_prec)) upper_prec else Inf
  
  # ---- Exact numeric degeneracy in the bounds ----
  if (is.finite(L) && is.finite(U) && L == U) {
    # Interval truly collapsed
    return(rep(L, nn))
  }
  
  # ---- One-sided truncation: lower only ----
  if (!is.null(lower_prec) && is.null(upper_prec)) {
    # Sample from Gamma conditioned on X >= L.
    # Work on the upper-tail CDF: T(x) = P(X >= x) = 1 - F(x)
    log_tail_L <- stats::pgamma(L, shape = shape, rate = rate,
                                lower.tail = FALSE, log.p = TRUE)
    # U ~ Unif(0,1), target tail probability: T_star = U * T(L)
    u      <- stats::runif(nn)
    log_u  <- log(u)
    log_T  <- log_u + log_tail_L  # log T_star
    
    # Invert using qgamma on the upper tail with log.p
    out <- stats::qgamma(p = log_T, shape = shape, rate = rate,
                         lower.tail = FALSE, log.p = TRUE)
    return(out)
  }
  
  # ---- One-sided truncation: upper only ----
  if (is.null(lower_prec) && !is.null(upper_prec)) {
    # Sample from Gamma conditioned on X <= U.
    # Work on the lower-tail CDF: F(x) = P(X <= x)
    log_F_U <- stats::pgamma(U, shape = shape, rate = rate,
                             lower.tail = TRUE, log.p = TRUE)
    
    u      <- stats::runif(nn)
    log_u  <- log(u)
    log_F  <- log_u + log_F_U  # log F_star
    
    out <- stats::qgamma(p = log_F, shape = shape, rate = rate,
                         lower.tail = TRUE, log.p = TRUE)
    return(out)
  }
  
  # ---- Two-sided truncation: L and U both finite ----
  # Now both lower_prec and upper_prec are non-NULL, and L < U (checked above).
  
  # Work on the lower-tail CDF: F(x) = P(X <= x)
  log_F_L <- stats::pgamma(L, shape = shape, rate = rate,
                           lower.tail = TRUE, log.p = TRUE)
  log_F_U <- stats::pgamma(U, shape = shape, rate = rate,
                           lower.tail = TRUE, log.p = TRUE)
  
  # Enforce ordering (should normally hold, but guard anyway)
  if (log_F_U < log_F_L) {
    tmp      <- log_F_L
    log_F_L  <- log_F_U
    log_F_U  <- tmp
    tmp      <- L
    L        <- U
    U        <- tmp
  }
  
  # If pgamma collapses the CDF in floating point, the conditional law
  # is essentially degenerate at L. We do not use an arbitrary epsilon;
  # we only test exact equality on the log scale.
  if (log_F_U == log_F_L) {
    return(rep(L, nn))
  }
  
  # Mass of the interval (L, U] on the CDF scale, in log:
  # delta = F(U) - F(L), log_delta = logdiffexp(log_F_L, log_F_U)
  log_delta <- logdiffexp(log_F_L, log_F_U)
  
  # U ~ Unif(0,1), target unconditional CDF value:
  # F_star = F(L) + U * (F(U) - F(L))
  # Work in log-space via log-sum-exp of two terms:
  #   term1 = log_F_L
  #   term2 = log(U) + log_delta
  u      <- stats::runif(nn)
  log_u  <- log(u)
  t1     <- log_F_L
  t2     <- log_u + log_delta
  
  m      <- pmax(t1, t2)
  log_F  <- m + log(exp(t1 - m) + exp(t2 - m))  # log(F_star)
  
  # Invert using qgamma on the lower tail with log.p
  out <- stats::qgamma(p = log_F, shape = shape, rate = rate,
                       lower.tail = TRUE, log.p = TRUE)
  out
}


## Constrained truncated gamma sampler on precision scale
# ctrgamma <- function(nn, shape, rate, lower_prec = NULL, upper_prec = NULL) {
#   
#   # Case: no truncation at all
#   if (is.null(lower_prec) && is.null(upper_prec)) {
#     return(stats::rgamma(nn, shape = shape, rate = rate))
#   }
#   
#   # Compute CDF bounds
#   a <- if (!is.null(lower_prec)) stats::pgamma(lower_prec, shape = shape, rate = rate) else 0
#   b <- if (!is.null(upper_prec)) stats::pgamma(upper_prec, shape = shape, rate = rate) else 1
#   
#   # Extract numeric bounds (for uniform fallback)
#   L <- if (!is.null(lower_prec)) lower_prec else -Inf
#   U <- if (!is.null(upper_prec)) upper_prec else  Inf
#   
#   # ---- Case 3: exact numeric degeneracy ----
#   if (!is.null(lower_prec) && !is.null(upper_prec) && lower_prec == upper_prec) {
#     return(rep(lower_prec, nn))
#   }
#   
#   # ---- Case 2: CDF interval collapsed but numeric interval valid ----
#   # (i.e., proposal has zero mass between L and U in floating point)
#   if (abs(b - a) < 1e-15) {
#     # Draw uniformly on [L, U]
#     u <- stats::runif(nn)
#     return(L + u * (U - L))
#   }
#   
#   # ---- Case 1: normal truncated Gamma ----
#   u <- stats::runif(nn, min = a, max = b)
#   stats::qgamma(u, shape = shape, rate = rate)
# }




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
    ## Center y by offset
    y1 <- as.numeric(y) - alpha
    
    ## Fixed linear predictor and residuals
    xb   <- drop(x %*% b)
    res  <- y1 - xb
    SS   <- res * res      # unweighted squared residuals
    
    ## Weighted posterior parameters for precision tau = 1/phi
    ## a1 = shape + 0.5 * sum(w_i)
    ## b1 = rate  + 0.5 * sum(w_i * res_i^2)
    sum_wt <- sum(wt)
    a1     <- shape + sum_wt / 2
    b1     <- rate  + sum(wt * SS) / 2

    draws <- integer(n)+1
    
        
    if (is.null(disp_lower) && is.null(disp_upper)) {
      ## Unconstrained sampler: sample precision then invert to dispersion
      out <- 1 / rgamma(n, shape = a1, rate = b1)
    } else {
      ## Constrained sampler on precision
      prec <- ctrgamma(
        n,
        shape = a1,
        rate  = b1,
        lower_prec = if (!is.null(disp_upper)) 1 / disp_upper else NULL,
        upper_prec = if (!is.null(disp_lower)) 1 / disp_lower else NULL
      )
      out <- 1 / prec
    }
  }
  
  ## ----------------------
  ## Gamma case (MODIFIED)
  ## ----------------------
  if (family$family == "Gamma") {
    
    ## Compute mu1 using the fixed coefficients
    mu1 <- t(exp(alpha + x %*% b))
    
    ## To build effective envelopes for the Gamma family with a Gamma prior
    ## for the precision, it is important to be able to center the proposal
    ## close to the center of the posterior. To enable this, we proceed in two steps
    ## 
    ## Step 1: We build an equivalent "proxy model" with the following three key properties
    ##
    ## (i)  shape_prop = shape_prior + (sum_i wt_i) / 2
    ## (ii) ll*(v)     = ll(v) - (sum_i wt_i / 2) * log v
    ## (iii) ll*(.) still retains the concavity of the log-likelihood
    ## 
    ## Step 2: We pick one or more tangency points and for each component of the grid
    ##         (if more than one), we shift a linear term from the new proxy model
    ##         to the rate. In the case of a single tangency, it is best to use a
    ##         tangency point at the posterior mean. Generally denote by cbar(v)
    ##         a gradient for -ll(v) at v. Then we build a proposal and a bounding
    ##         function log_h(.) as follows:
    ##
    ## (i)  A Gamma proposal with
    ##      (a) shape_prop = shape_prior + (sum_i wt_i) / 2
    ##      (b) rate_prop  = rate_prior - cbar(v_tangent)
    ##
    ##      where -cbar(v_tangent) =
    ##          sum_i wt_i [log(v_tangent) + 1 - log(mu_i) + log(y_i)
    ##                       - digamma(v_tangent) - y_i / mu_i]
    ##
    ##       Note: Remove this from cbar - Moves proposal too far from the posterior 
    ##                       
    ##          - (sum_i wt_i) / (2 v_tangent)
    ##
    ## (ii) A bounding function log_h(.) defined by
    ##
    ##      log_h(v) = ll(v) - [ll(v_tangent) + ll'(v_tangent)(v - v_tangent)] -[(sum_i wt_i) / 2] *(log v/v_min)]
    ##               = ll(v) - [ll*(v_tangent) - cbar(v_tangent)(v - v_tangent)]-[(sum_i wt_i) / 2] *(log v/v_min)]
    ##      where cbar(v_tangent) = -ll'(v_tangent).
    ##
    ##  To implement the sampler, we proceed as follows:
    ##
    ##  (i)   Find the posterior mode (v_star)
    ##  (ii)  Use curvature at v_star to approximate the posterior mean v_mean
    ##  (iii) Build the proposal and bounding functions using v_tangent = v_mean
    ##  (iv)  Sample as follows:
    ##        (a) Generate a sample from the proposal Gamma(shape_prop, rate_prop)
    ##        (b) Draw u ~ Uniform(0,1)
    ##        (c) Accept if lg_h(cand) - log(u) > 0
    
    ## ----------------------------------------------------------------------
    ## 1. Setup: means, weights, and convenience quantities
    ## ----------------------------------------------------------------------
    
    mu_vec <- as.numeric(mu1)
    W_sum  <- sum(wt)
    
    ## 2. Original log-likelihood and proxy log-likelihood ll*(v)
    
    loglik <- function(v) {
      sum(
        wt * (
          v * log(v) -
            v * log(mu_vec) +
            (v - 1) * log(y) -
            lgamma(v) -
            v * y / mu_vec
        )
      )
    }
    
    ## ll*(v) = ll(v) - 0.5 * sum_i wt_i * log(v)
    loglik_star <- function(v) {
      loglik(v) - 0.5 * W_sum * log(v)
    }
    
    ## 3. Prior: v ~ Gamma(shape, rate)
    ##    rate1: deviance-based piece used below in an approximate mode iteration
    
    rate1 <- rate + sum(wt * ((y / mu1) - log(y / mu1) - 1))
    
    ## ----------------------------------------------------------------------
    ## 4. Posterior mode finder (v_star) via a Gamma-surrogate fixed-point step
    ## ----------------------------------------------------------------------
    ## We use a Gamma-like surrogate with "shape2_mode" and rate1 to locate
    ## an approximate posterior mode v_star, then refine the posterior mean
    ## using curvature at v_star.
    
    shape2_mode <- shape
    vstar1      <- shape2_mode / rate1
    
    vout_mode <- function(v) {
      vstar1 - (v / rate1) *
        sum(wt * digamma(wt * v) - wt * log(wt * v))
    }
    
    vstar <- vstar1
    for (j in 1:20) {
      vstar <- vout_mode(vstar)
    }
    v_mode <- vstar
    
    ## 5. Exact second derivative f''(v) of the log-posterior for curvature
    ##    f(v) = loglik(v) + (shape - 1) log v - rate v
    
    f2_exact <- function(v) {
      term1 <- W_sum * (1 / v)       # sum_i w_i * 1/v   from log v terms
      term2 <- -W_sum * trigamma(v)  # - sum_i w_i * psi1(v) from -lgamma(v)
      term3 <- -(shape - 1) / (v^2)  # prior curvature
      term1 + term2 + term3
    }
    
    f2_vstar  <- f2_exact(vstar)
    
    ## 6. Gamma approximation to posterior: v | y ~ approx Gamma(alpha_bar, beta_bar)
    
    alpha_bar <- 1 - f2_vstar * vstar^2
    beta_bar  <- -f2_vstar * vstar
    
    v_mean <- alpha_bar / beta_bar
    
    v_min <- qgamma(1-max_disp_perc, shape = alpha_bar, rate = beta_bar)
    v_max <- qgamma(max_disp_perc, shape = alpha_bar, rate = beta_bar)

    
    ## New: extract optional low/upp from prior_list when available
    if (!is.null(prior_list$disp_lower))  disp_lower <- prior_list$disp_lower  else disp_lower <- 1/v_max
    if (!is.null(prior_list$disp_upper))  disp_upper <- prior_list$disp_upper  else disp_upper <- 1/v_min
    
    # After computing disp_lower and disp_upper (including defaults)
    if (disp_lower >= disp_upper) {
      stop("Final dispersion bounds invalid: disp_lower must be < disp_upper.")
    }
    
    v_min=1/disp_upper
    v_max=1/disp_lower
    
        
    ## Optional diagnostics on the surrogate posterior for v
    
#    cat("\nApproximate posterior Gamma surrogate (precision v):\n")
#    cat("  shape_bar          =", alpha_bar, "\n")
#    cat("  rate_bar           =", beta_bar,  "\n")
#    cat("  v_mode (precision) =", v_mode,    "\n")
#    cat("  v_mean (precision) =", v_mean,    "\n\n")
#    cat("  v_min (precision)  =", v_min,     "\n\n")
    
    ## 7. Tangency point at approximate posterior mean (A6-optimal anchor)
    
    v_tangent <- v_mean
    
    ## 8. Proxy log-likelihood at the tangent point
    
    ll_star_v0 <- loglik_star(v_tangent)
    ll_v0      <- loglik(v_tangent)
    
#    cat("  ll_v0        =", ll_v0,      "\n")
#    cat("  ll_star_v0   =", ll_star_v0, "\n")
    
    ## 9. Gradient of -ll*(v) at v_tangent: cbar(v_tangent)
    ##
    ##    ll'(v) = sum_i w_i [log v + 1 - log mu_i + log y_i - digamma(v) - y_i / mu_i]
    ##    ll*'(v) = ll'(v) - 0.5 * sum_i wt_i * 1/v
    
    ell_prime_at <- function(v) {
      sum(
        wt * (
          log(v) + 1 -
            log(mu_vec) +
            log(y) -
            digamma(v) -
            y / mu_vec
        )
      )
    }
    
    ## Temporarily remove adjustment term
    
    ell_star_prime_at <- function(v) {
      ell_prime_at(v) #- 0.5 * W_sum / v
    }
    
#    cat("  ell_prime_tangent      =", ell_prime_at(v_tangent),      "\n")
#    cat("  ell_star_shape_shift   =", -0.5 * W_sum / v_tangent,     "\n")
#    cat("  ell_star_prime_tangent =", ell_star_prime_at(v_tangent), "\n")
    
    ## cbar is the gradient of -ll*(v) at v_tangent
    cbar <- -ell_star_prime_at(v_tangent)
    
    ## 10. Proposal Gamma(shape_prop, rate_prop)
    ##     (i)  shape_prop = shape + 0.5 * sum_i wt_i
    ##     (ii) rate_prop  = rate - cbar(v_tangent)
    
    shape_prop <- shape + 0.5 * W_sum
    rate_prop  <- rate - cbar
    
    if (rate_prop <= 0) {
      stop("Gamma proposal rate_prop <= 0; check v_tangent and curvature diagnostics.")
    }
    
    
    #    cat("  shape_Proposal =", shape_prop, "\n")
    #    cat("  rate_Proposal  =", rate_prop,  "\n")
    
    ## Proposal distribution diagnostics
    prop_mean <- shape_prop / rate_prop
    prop_mode <- if (shape_prop > 1) (shape_prop - 1) / rate_prop else NA
    
#    cat("  Proposal mean (precision) =", prop_mean, "\n")
#    cat("  Proposal mode (precision) =", prop_mode, "\n\n")
    
    ## 11. log_h(v) based on the concave ll*(v)
    ##
    ## Concavity of ll* implies:
    ##   ll*(v) <= ll*(v0) + ll*'(v0)(v - v0)
    ## With cbar = -ll*'(v0), the remainder is:
    ##   log h(v) = ll*(v) - (ll*(v0) - cbar (v - v0) ) <= 0.
    
    lg_h <- function(v) {
      loglik_star(v) - (ll_star_v0 - cbar * (v - v_tangent))-  0.5 * W_sum * log(v / v_min)
    }
    
    ## 12. Rejection sampling for precision v
    
    out   <- numeric(n)
    test  <- numeric(n)
    a     <- numeric(n)
    draws <- integer(n)
    
    for (i in 1:n) {
      draws[i] <- 1
      while (a[i] == 0) {
        
        if (is.null(disp_lower) && is.null(disp_upper)) {
#          cand_prec <- rgamma(1, shape = shape_prop, rate = rate_prop)
          cand_prec <- ctrgamma(1, shape = shape_prop, rate = rate_prop,lower_prec=v_min,upper_prec=v_max)
          
          
        } else {
          cand_prec <- ctrgamma(
            1,
            shape = shape_prop,
            rate  = rate_prop,
            lower_prec = if (!is.null(disp_upper)) 1 / disp_upper else NULL,
            upper_prec = if (!is.null(disp_lower)) 1 / disp_lower else NULL
          )
        }
        out[i] <- cand_prec
        
        ## Accept–reject based on log_h with additional v_min penalty:
        ##   - (W_sum / 2) * log(v / v_min)
        test[i] <- lg_h(out[i]) - log(runif(1, 0, 1))

#          0.5 * W_sum * log(out[i] / v_min) -
        
        if (test[i] > 0) {
          a[i] <- 1
        } else {
          draws[i] <- draws[i] + 1
        }
      }
    }
    
    ## Convert precision to dispersion
    out <- 1 / out
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
#' @keywords internal
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
  
  shape <- object$Prior$shape
  rate  <- object$Prior$rate
  
  # Prior mean of dispersion
  if (shape > 1) {
    prior_mean_disp <- rate / (shape - 1)
  } else {
    prior_mean_disp <- NA
  }
  
  # Prior SD of dispersion
  if (shape > 2) {
    prior_sd_disp <- rate / ((shape - 1) * sqrt(shape - 2))
  } else {
    prior_sd_disp <- NA
  }
  
  Tab1 <- cbind(
    "Prior.Mean" = prior_mean_disp,
    "Prior.Sd"   = prior_sd_disp,
    "Approx.Prior.wt" = Priorwt
  )
  
  rownames(Tab1) <- "dispersion"
  colnames(Tab1) <- c("Prior.Mean", "Prior.Sd", "Approx.Prior.wt")
  
  
  
  
#  Tab1<-cbind("Prior.Mean"=object$Prior$shape/object$Prior$rate,"Prior.Sd"=sqrt(object$Prior$shape)/object$Prior$rate
#              ,"Approx.Prior.wt"=Priorwt  )

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
  
  colnames(Tab1) <- c("Prior.Mean", "Prior.Sd", "Approx.Prior.wt")
  colnames(TAB) <- c("Post.Mean", "Post.Sd", "MC.Error", "Pr(tail)")
  
  res<-list(call=object$call,
            n=n,
            coefficients1=Tab1,
            coefficients=TAB,
            Percentiles=TAB2
  )
  
  # Reuse summary.rglmb class
  
  class(res) <- "summary.rGamma_reg"  
  res
  
}


#' @export
#' @keywords internal
#' @rdname simfuncs
#' @order 9
#' @method print summary.rGamma_reg


print.summary.rGamma_reg <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  
  ## --- Call ---
  cat("Call\n")
  print(x$call)
  cat("\n")
  
  ## --- Prior Estimates ---
  cat("Prior Estimates with Standard Deviations\n\n")
  print(round(x$coefficients1, digits))
  cat("\n")
  
  ## --- Posterior Estimates ---
  cat("Bayesian Estimates Based on", x$n, "iid draws\n\n")
  
  # Extract posterior table
  TAB <- round(x$coefficients, digits)
  
  # Compute significance stars
  pvals <- x$coefficients[, "Pr(tail)"]
  stars <- ifelse(pvals < 0.001, "***",
                  ifelse(pvals < 0.01,  "**",
                         ifelse(pvals < 0.05,  "*",
                                ifelse(pvals < 0.1,   ".", " "))))
  
  # Build final table with stars
  TAB2 <- cbind(TAB, Signif = stars)
  
  print(TAB2)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  ## --- Percentiles ---
  cat("Distribution Percentiles\n\n")
  print(round(x$Percentiles, digits))
  cat("\n")
  
  ## --- Dispersion summary ---
  disp.mean <- x$coefficients["dispersion", "Post.Mean"]
  cat("Expected Mean dispersion:", round(disp.mean, digits), "\n")
  cat("Sq.root of Expected Mean dispersion:", round(sqrt(disp.mean), digits), "\n\n")
  
  invisible(x)
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
  
  if (verbose) {
    start_anchor <- as.numeric(Sys.time())
    cat("[DispersionAnchoring] >>> Entering iterative dispersion anchoring at",
        format(Sys.time(), "%H:%M:%S"), "<<<\n")
  }
  
  
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
  
  
  if (verbose) {
    end_anchor <- as.numeric(Sys.time())
    elapsed <- end_anchor - start_anchor
    h <- as.integer(elapsed / 3600)
    m <- as.integer((elapsed - h*3600) / 60)
    s <- as.integer(elapsed - h*3600 - m*60)
    
    cat("[DispersionAnchoring] >>> Exiting iterative dispersion anchoring at",
        format(Sys.time(), "%H:%M:%S"), "<<<\n")
    cat("[DispersionAnchoring] Dispersion anchoring completed in:",
        h, "h ", m, "m ", s, "s.\n")
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

  R <- chol(Sigma)
  P <- chol2inv(R)
  P <- 0.5 * (P + t(P))   # enforce symmetry
  
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
  
  if (verbose) {
    start_optim <- as.numeric(Sys.time())
    cat("[PosteriorMode] >>> Entering optim() call at",
        format(Sys.time(), "%H:%M:%S"), "<<<\n")
  }
  
  
  opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x2),mu=as.vector(mu2),
                P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
                method="BFGS",hessian=TRUE
  )
  
  bstar=opt_out$par  ## Posterior mode for adjusted model
  
  ## Temporarily use bstar as posterior mode
  
  
  #  bstar
  #  bstar+as.vector(mu)  # mode for actual model
  A1=opt_out$hessian # Approximate Precision at mode
  
  
  if (verbose) {
    end_optim <- as.numeric(Sys.time())
    elapsed <- end_optim - start_optim
    h <- as.integer(elapsed / 3600)
    m <- as.integer((elapsed - h*3600) / 60)
    s <- as.integer(elapsed - h*3600 - m*60)
    
    cat("[PosteriorMode] >>> Exiting optim() call at",
        format(Sys.time(), "%H:%M:%S"), "<<<\n")
    cat("[PosteriorMode] optim() completed in:",
        h, "h ", m, "m ", s, "s.\n")
  }
  
  
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
  
  if (verbose) {
    start_envbuild <- as.numeric(Sys.time())
    cat("[EnvelopeBuild] >>> Entering EnvelopeBuild at",
        format(Sys.time(), "%H:%M:%S"), "<<<\n")
  }
  
  
  Env2=EnvelopeBuild(as.vector(bstar2), as.matrix(A),y, as.matrix(x2),
                     as.matrix(mu2,ncol=1),as.matrix(P2),as.vector(alpha),
                     as.vector(wt2),
                     family="gaussian",link="identity",
                     Gridtype=Gridtype, n=as.integer(n),
                     n_envopt=n_envopt,
                     sortgrid=TRUE,use_opencl = use_opencl,verbose = verbose)
  
  if (verbose) {
    end_envbuild <- as.numeric(Sys.time())
    elapsed <- end_envbuild - start_envbuild
    h <- as.integer(elapsed / 3600)
    m <- as.integer((elapsed - h*3600) / 60)
    s <- as.integer(elapsed - h*3600 - m*60)
    
    cat("[EnvelopeBuild] >>> Exiting EnvelopeBuild at",
        format(Sys.time(), "%H:%M:%S"), "<<<\n")
    cat("[EnvelopeBuild] EnvelopeBuild completed in:",
        h, "h ", m, "m ", s, "s.\n")
  }
  
  
  
  ###  Call new function to build shared envelope

  if (verbose) {

    start_dispbuild <- as.numeric(Sys.time())
    
    cat("[EnvelopeDispersionBuild] >>> Entering EnvelopeDispersionBuild at",
        format(Sys.time(), "%H:%M:%S"), "<<<\n")
  }
  
  
  
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
    verbose    = verbose,   # optional, matches the C++ signature
    use_parallel=use_parallel
  )
  
  
  if (verbose) {
    end_dispbuild <- as.numeric(Sys.time())
    elapsed <- end_dispbuild - start_dispbuild
    h <- as.integer(elapsed / 3600)
    m <- as.integer((elapsed - h*3600) / 60)
    s <- as.integer(elapsed - h*3600 - m*60)
    
    cat("[EnvelopeDispersionBuild] >>> Exiting EnvelopeDispersionBuild at",
        format(Sys.time(), "%H:%M:%S"), "<<<\n")
    cat("[EnvelopeDispersionBuild] EnvelopeDispersionBuild completed in:",
        h, "h ", m, "m ", s, "s.\n")
  }
  
  
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
    lg_prob_factor = disp_env_out$UB_list$lg_prob_factor,
    UB2min = disp_env_out$UB_list$UB2min
  )

  

  # Env3_temp           <- disp_env_out$Env_out
  gamma_list_new <- disp_env_out$gamma_list
  UB_list_new    <- disp_env_out$UB_list
  low            <- gamma_list_new$disp_lower
  upp            <- gamma_list_new$disp_upper  
  diagnostics     <- disp_env_out$diagnostics
  
  ## USE reordered lg_prob factor and UB2min
  UB_list_new$lg_prob_factor=Env3$lg_prob_factor
  UB_list_new$UB2min=Env3$UB2min
  

  

  
  if (!use_parallel || n == 1) {
    if (verbose) {
      start_sim <- as.numeric(Sys.time())
      
      cat("[Simulation] >>> Calling .rindep_norm_gamma_reg_std_cpp at",
          format(Sys.time(), "%H:%M:%S"), "<<<\n")
    }
    
    sim_temp <- .rindep_norm_gamma_reg_std_cpp(
      n = n, y = y, x = x2,
      mu = mu2,  # Should be zero vector
      P = P2,    # Part of prior shifted to the Likelihood
      alpha = alpha, wt = wt,
      f2 = f2, Envelope = Env3,
      gamma_list = gamma_list_new,
      UB_list = UB_list_new,
      family = "gaussian", link = "identity",
      progbar = progbar,
      verbose=verbose
    )
    
    if (verbose) {
      end_sim <- as.numeric(Sys.time())
      elapsed <- end_sim - start_sim
      h <- as.integer(elapsed / 3600)
      m <- as.integer((elapsed - h*3600) / 60)
      s <- as.integer(elapsed - h*3600 - m*60)
      

      cat("[Simulation] >>> Exiting .rindep_norm_gamma_reg_std_cpp at",
          format(Sys.time(), "%H:%M:%S"), "<<<\n")
      cat("[Simulation] Simulation completed in: ",
          h, "h ", m, "m ", s, "s.\n")
      
            
    }
    
  } else {
    
    if (verbose) {
      start_sim <- as.numeric(Sys.time())
      
      cat("[Simulation] >>> Calling .rindep_norm_gamma_reg_std_parallel_cpp at",
          format(Sys.time(), "%H:%M:%S"), "<<<\n")
    }
    
    sim_temp <- .rindep_norm_gamma_reg_std_parallel_cpp(
      n = n, y = y, x = x2,
      mu = mu2,  # Should be zero vector
      P = P2,    # Part of prior shifted to the Likelihood
      alpha = alpha, wt = wt,
      f2 = f2, Envelope = Env3,
      gamma_list = gamma_list_new,
      UB_list = UB_list_new,
      family = "gaussian", link = "identity",
      progbar = progbar,
      verbose=verbose
    )
    
    if (verbose) {
      end_sim <- as.numeric(Sys.time())
      elapsed <- end_sim - start_sim
      h <- as.integer(elapsed / 3600)
      m <- as.integer((elapsed - h*3600) / 60)
      s <- as.integer(elapsed - h*3600 - m*60)
      
    

      cat("[Simulation] >>> Exiting .rindep_norm_gamma_reg_std_parallel_cpp at",
          format(Sys.time(), "%H:%M:%S"), "<<<\n")
      cat("[Simulation] Simulation completed in: ",
          h, "h ", m, "m ", s, "s.\n")
      
        
    }
  }
  
  
##  print(paste("Interactive status:", interactive()))
  


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
  
  # Effective sample size for Gaussian likelihood
  sum_wt <- sum(wt)  # wt has already been expanded to length nobs
  
  a_prior=shape     ## Should be relationship to shape in Wishart  
  b_prior=rate  ## Should be relationship to scale in Wishart (could also be V/2)
  
  # Posterior Gamma parameters for precision tau = 1/phi
  a_post <- a_prior + sum_wt / 2      # instead of nobs / 2
  b_post <- b_prior + 0.5 * S         # S is weighted RSS from rNormal_reg.wfit
  
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

