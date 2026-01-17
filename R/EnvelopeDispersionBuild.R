#' @name dispenvelopes 
#' @title 
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
#' @param cache       Precomputed set of outputs passed downstream
#' @param cbars_j     gradient vector for standarized negative log-likelihood function
#' @param wt          weight vector
#' @param rss_min_global Minimum RSS across the face tangencies over the range defined by disp_lower and disp_upper
#' @param dispersion  Dispersion parameter
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
#' @rdname dispenvelopes
#' @order 1


EnvelopeDispersionBuild <- function(Env, Shape, Rate, P, y, x, alpha, n_obs, RSS_post, RSS_ML,
                                    mu, wt, max_disp_perc = 0.99, disp_lower = NULL, disp_upper = NULL,
                                    verbose = FALSE, use_parallel = TRUE) {
  .Call(`_glmbayes_EnvelopeDispersionBuild_cpp`,
        Env, Shape, Rate, P, y, x, alpha, n_obs,
        RSS_post, RSS_ML, mu, wt,
        max_disp_perc, disp_lower, disp_upper,
        verbose, use_parallel)
}






#' @noRd

EnvBuildLinBound<-function(thetabars,cbars,y,x2,P2,alpha,dispstar){
  
  # gs=nrow(cbars)
  # n_vars=ncol(cbars)
  # 
  # New_LL_Slope_test2=c(1:gs)
  # New_LL_Slope_test3=c(1:gs)
  
  XtX   <- crossprod(x2)
  rhs   <- crossprod(x2, y - alpha)
  M     <- XtX + dispstar * P2
  # Replace solve(M) with Cholesky inversion
  R    <- chol(M)
  Minv <- chol2inv(R)
  Minv <- 0.5 * (Minv + t(Minv))   # enforce symmetry
  H1    <- -Minv %*% P2 %*% Minv
  
  V <- -(thetabars %*% P2) + cbars          # gs x p
  Minv_cbars <- t(Minv %*% t(cbars))        # gs x p
  term1 <- rowSums(V * Minv_cbars)
  
  # replicate rhs across gs columns
  rhs_mat <- matrix(rhs, nrow = length(rhs), ncol = nrow(cbars))
  H1_rhs  <- t(H1 %*% (rhs_mat + dispstar * t(cbars)))  # gs * p
  
  term2 <- rowSums(V * H1_rhs)
  
  New_LL_Slope <- term1 + term2
  
  return(New_LL_Slope)

}

#' @noRd

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


#' @noRd

Inv_f3_with_disp <- function(cache, dispersion, cbars_small) {
  .Call(`_glmbayes_Inv_f3_with_disp`, cache, dispersion, cbars_small)
}

#' @usage UB2(dispersion, cache, cbars_j, y, x, alpha, wt, rss_min_global)
#' @export
#' @rdname dispenvelopes
#' @order 5


UB2 <- function(dispersion, cache, cbars_j, y, x, alpha, wt, rss_min_global) {
  .Call(`_glmbayes_UB2`, dispersion, cache, cbars_j, y, x, alpha, wt, rss_min_global)
}

#' @usage rss_face_at_disp(dispersion, cache, cbars_j, y, x, alpha, wt)
#' @export
#' @rdname dispenvelopes
#' @order 6

rss_face_at_disp<- function(dispersion, cache, cbars_j, y, x, alpha, wt) {
  .Call(`_glmbayes_rss_face_at_disp`, dispersion, cache, cbars_j, y, x, alpha, wt)
}



#' @usage drss_ddisp(dispersion, cache, cbars_j, y, x, alpha, wt)
#' @export
#' @rdname dispenvelopes
#' @order 7

drss_ddisp <- function(dispersion, cache, cbars_j, y, x, alpha, wt) {
  .Call(`_glmbayes_drss_ddisp`, dispersion, cache, cbars_j, y, x, alpha, wt)
}

#' @noRd
#' @order 11

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
          fn     = rss_face_at_disp,
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
#' @order 8

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
#' @order 12


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
          fn     = UB2,                 # UB2 objective function
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
#' @order 9

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