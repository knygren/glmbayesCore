#' Builds Dispersion-Aware Envelope for Simulation
#'
#' Constructs a dispersion-aware envelope for simulation in Gaussian models with uncertain variance.
#' This function extrapolates the coefficient envelope across a high-probability interval for the
#' dispersion parameter \code{sigma^2}, and builds a global upper bound for the log-posterior remainder.
#' It also computes mixture weights for envelope faces and adjusts the Gamma proposal for precision.
#'
#' The envelope is constructed using the slopes of the face constants with respect to dispersion,
#' evaluated at an anchor point. The resulting structure supports exact i.i.d. sampling via
#' accept–reject correction.
#'
#' The procedure follows these steps:
#' \enumerate{
#'   \item Compute posterior Gamma parameters for precision using the prior and posterior-predictive RSS.
#'   \item Determine a central credible interval for dispersion using \code{max_disp_perc}.
#'   \item Evaluate the slope of each envelope face constant with respect to dispersion.
#'   \item Linearly extrapolate each face constant to the upper and lower bounds of the dispersion interval.
#'   \item Construct a global upper bound for the log-posterior remainder across dispersion.
#'   \item Compute mixture weights for each face based on tightness at the interval endpoints.
#'   \item Adjust the Gamma shape parameter to match the envelope slope at the anchor point.
#' }
#'
#' @param Env        Envelope object from \code{\link{EnvelopeBuild}}, containing tangency points and gradients
#' @param Shape      Prior shape parameter for precision \code{v = 1 / sigma^2}
#' @param Rate       Prior rate parameter for precision
#' @param P          Prior precision matrix for coefficients
#' @param y          Design matrix
#' @param x          Design matrix
#' @param alpha          Design matrix
#' @param n_obs      Number of observations
#' @param RSS_post   Posterior-predictive residual sum of squares
#' @param RSS_ML   Residual sum of squares associated with MLE estimate
#' @param max_disp_perc Truncation level for dispersion (default 0.99)
#' @param verbose Option to have verbose output
#'
#' @return A list with elements:
#'   \item{Env_out}{Updated envelope object with dispersion-aware mixture weights}
#'   \item{gamma_list}{List of parameters for truncated Gamma proposal}
#'   \item{UB_list}{Constants for evaluating the dispersion envelope}
#'   \item{diagnostics}{Optional diagnostics including slopes and bounds}
#'
#' @details
#' This function is designed to complement \code{\link{EnvelopeBuild}} for Gaussian models
#' with Normal-Gamma priors. It enables exact sampling of both coefficients and dispersion
#' by constructing a joint envelope that respects posterior curvature in both dimensions.
#'
#' The dispersion anchor point is chosen as the log-scale center of the credible interval,
#' and the Gamma proposal is tilted to match the envelope slope at this point.
#'
#' @seealso \code{\link{EnvelopeBuild}}, \code{\link{glmb}}, \code{\link{glmbfamfunc}}
#' @export



EnvelopeDispersionBuild <- function(
    Env, Shape, Rate, P, y, x, alpha,
    n_obs, RSS_post, RSS_ML,
    max_disp_perc = 0.99,
    verbose = FALSE
) {

  # Step 1: Posterior Gamma parameters (precision prior)
  shape2 <- Shape + n_obs / 2
  rate3  <- Rate  + RSS_post / 2
  
  # Step 2: Dispersion bounds (on sigma^2)
  low <- 1 / qgamma(max_disp_perc,     shape2, rate3)
  upp <- 1 / qgamma(1 - max_disp_perc, shape2, rate3)
  
  
  
  
  # Step 3: Extract envelope faces
  cbars     <- Env$cbars
  thetabars <- Env$thetabars
  logP1     <- Env$logP
  gs        <- nrow(cbars)
  
  # Step 4: Base face constants
  thetabar_const_base <- thetabar_const(P, cbars, thetabars)
  
  # Step 5: initial anchor (posterior mean, optional)
  dispstar <- rate3 / (shape2 - 1)
  
  
  
  # Step 6: Face slopes at dispstar
  New_LL_Slope <- EnvBuildLinBound(thetabars, cbars, y, x, P, alpha, dispstar)
  
  # Step 7: Linear extrapolation of face constants to bounds
  thetabar_const_upp_apprx <- thetabar_const_base + (upp - dispstar) * New_LL_Slope
  thetabar_const_low_apprx <- thetabar_const_base + (low - dispstar) * New_LL_Slope
  
  
  # Step 8: Global upper line geometry (match original mean-slope correction)
  
  ## Compute max face constants at endpoints
  max_low = max(thetabar_const_low_apprx)
  max_upp = max(thetabar_const_upp_apprx)
  
  ## This line is a no‑op in the original, but keep it for parity
  max_low = max_low + 0 * (max_upp - max_low)
  

  
  ### Global upper line over dispersion
  m_New_LL_Slope = mean(New_LL_Slope)
  max_low_mean   = max_upp - m_New_LL_Slope * (upp - low)
  old_slope      = (max_upp - max_low) / (upp - low)
  max_low        = max_low_mean
  
  new_slope = (max_upp - max_low) / (upp - low)
  new_int   = max_low - new_slope * low

  # Step 9a: Dispersion anchor (exactly as in original: b1/(-c1))
  b1 <- (upp - low)
  c1 <- -log(upp / low)
  dispstar <- b1 / (-c1)  # equivalently (upp - low)/log(upp/low)
  
    

  # Step 9: Mixture weights per face (match original)
  New_logP2  <- numeric(gs)
  prob_factor <- numeric(gs)
  for (j in 1:gs) {
    cbars_temp <- matrix(cbars[j, ], ncol = 1)
    New_logP2[j] <- logP1[j] + 0.5 * t(cbars_temp) %*% cbars_temp
    prob_factor[j] <- max(
      thetabar_const_upp_apprx[j] - max_upp,
      thetabar_const_low_apprx[j] - max_low
    )
  }
  lg_prob_factor <- prob_factor
  prob_factor <- exp(New_logP2 + prob_factor)
  prob_factor <- prob_factor / sum(prob_factor)
  
  # Step 10: Envelope constants for dispersion and gamma tilt
  lm_log2 <- new_slope * dispstar
  lm_log1 <- new_int + new_slope * dispstar - new_slope * log(dispstar)
  shape3  <- shape2 - lm_log2
  
  # Step 11: Package outputs
  Env_out <- Env
  Env_out$PLSD <- prob_factor
  

  
  gamma_list <- list(
    shape3      = shape3,
    rate2       = Rate + RSS_ML/2,   # matches original definition
    disp_upper  = upp,
    disp_lower  = low
  )
  
  UB_list <- list(
    RSS_ML          = RSS_ML,        # not RSS_post
    max_New_LL_UB   = max_upp,
    max_LL_log_disp = lm_log1 + lm_log2 * log(upp),
    lm_log1         = lm_log1,
    lm_log2         = lm_log2,
    lg_prob_factor  = lg_prob_factor,
    lmc1            = new_int,
    lmc2            = new_slope
  )
  
  
  diagnostics <- list(
    dispstar         = dispstar,
    New_LL_Slope     = New_LL_Slope,
    shape2           = shape2,
    rate3            = rate3,
    shape3           = shape3,
    max_low          = max_low,
    max_upp          = max_upp,
    new_slope        = new_slope,
    new_int          = new_int,
    prob_factor      = prob_factor
  )
  
  if (verbose) {
    cat("EnvelopeDispersionBuild diagnostics:\n")
    cat("  dispstar      =", dispstar, "\n")
    cat("  new_slope     =", new_slope, "\n")
    cat("  new_int       =", new_int, "\n")
    cat("  lm_log1       =", lm_log1, "\n")
    cat("  lm_log2       =", lm_log2, "\n")
    cat("  shape3        =", shape3, "\n")
    cat("  max_low       =", max_low, "\n")
    cat("  max_upp       =", max_upp, "\n")
    cat("  prob_factor   =", paste(round(prob_factor, 4), collapse = ", "), "\n")
  }
  
  list(
    Env_out     = Env_out,
    gamma_list  = gamma_list,
    UB_list     = UB_list,
    diagnostics = diagnostics
  )
}