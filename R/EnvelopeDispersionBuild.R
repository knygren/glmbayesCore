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
#'         Using the prior and posterior-predictive RSS, set
#'
#'         \deqn{ \mathrm{shape2} = \mathrm{Shape} + n_{\mathrm{obs}}/2, \quad
#'                \mathrm{rate3}  = \mathrm{Rate} + \mathrm{RSS}_{\mathrm{post}}/2. }
#'
#'         These parameterize the posterior precision \eqn{v = 1/\sigma^2 \sim \mathrm{Gamma}(\mathrm{shape2}, \mathrm{rate3})}.
#'
#'   \item \strong{Central credible interval for dispersion (low, upp).}
#'         Choose a central mass level \code{max_disp_perc} (e.g., 0.99) for precision, then invert the corresponding
#'         Gamma quantiles to dispersion:
#'
#'         \deqn{ \mathrm{low} = 1 / Q_{\Gamma}(\mathrm{max\_disp\_perc}; \mathrm{shape2}, \mathrm{rate3}), \quad
#'                \mathrm{upp} = 1 / Q_{\Gamma}(1 - \mathrm{max\_disp\_perc}; \mathrm{shape2}, \mathrm{rate3}). }
#'
#'         The interval \eqn{[\mathrm{low}, \mathrm{upp}]} is the domain over which all envelopes must dominate.
#'
#'   \item \strong{Face slopes at an anchor (dispstar).}
#'         \deqn{ \mathrm{dispstar} = \mathrm{rate3} / (\mathrm{shape2} - 1) \quad \text{(posterior mean of } \sigma^2\text{)}. }
#'
#'         Compute \eqn{\mathrm{New\_LL\_Slope}_j} for each face \eqn{j}.
#'
#'   \item \strong{Linear extrapolation of face constants.}
#'         \deqn{ \theta^{\mathrm{low}}_{\bar{j}} = \theta^{\mathrm{base}}_{\bar{j}} + (\mathrm{low} - \mathrm{dispstar}) \cdot \mathrm{New\_LL\_Slope}_j, }
#'         \deqn{ \theta^{\mathrm{upp}}_{\bar{j}} = \theta^{\mathrm{base}}_{\bar{j}} + (\mathrm{upp} - \mathrm{dispstar}) \cdot \mathrm{New\_LL\_Slope}_j. }
#'
#'   \item \strong{Global upper line and endpoint maxima.}
#'         \deqn{ \mathrm{max\_low} = \max_j \theta^{\mathrm{low}}_{\bar{j}}, \quad
#'                \mathrm{max\_upp} = \max_j \theta^{\mathrm{upp}}_{\bar{j}}. }
#'
#'         \deqn{ \mathrm{new\_slope} = (\mathrm{max\_upp} - \mathrm{max\_low}) / (\mathrm{upp} - \mathrm{low}), \quad
#'                \mathrm{new\_int}   = \mathrm{max\_low} - \mathrm{new\_slope} \cdot \mathrm{low}. }
#'
#'   \item \strong{Face slack and mixture weights.}
#'         \deqn{ \mathrm{lg\_prob\_factor}_j =
#'                \max\big(\theta^{\mathrm{upp}}_{\bar{j}} - \mathrm{max\_upp},\;
#'                          \theta^{\mathrm{low}}_{\bar{j}} - \mathrm{max\_low}\big). }
#'
#'         Combine with \eqn{\mathrm{New\_logP2}_j = \mathrm{logP}_j + \tfrac{1}{2}\|\bar{c}_j\|^2}
#'         to form mixture weights
#'         \eqn{\mathrm{PLSD}_j \propto \exp\!\bigl(\mathrm{New\_logP2}_j + \mathrm{lg\_prob\_factor}_j\bigr)}.
#'
#'   \item \strong{Gamma tilt and dispersion-axis envelope.}
#'         \deqn{ \mathrm{dispstar} = (\mathrm{upp} - \mathrm{low}) / \log(\mathrm{upp}/\mathrm{low}). }
#'
#'         \deqn{ \mathrm{lm\_log2} = \mathrm{new\_slope} \cdot \mathrm{dispstar}, \quad
#'                \mathrm{lm\_log1} = \mathrm{new\_int} + \mathrm{new\_slope} \cdot \mathrm{dispstar}
#'                                   - \mathrm{new\_slope} \cdot \log(\mathrm{dispstar}). }
#'
#'         Tilt the Gamma proposal via \eqn{\mathrm{shape3} = \mathrm{shape2} - \mathrm{lm\_log2}}.
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
#' @param disp_lower lower bound truncation for dispersion 
#' @param disp_upper upper bound truncation for dispersion
#' @param verbose Option to have verbose output
#' @param thetabars   var1
#' @param cbars       var2
#' @param thetabars   var3
#' @param x2          var4
#' @param P2          var5
#' @param dispstar    var6
#' @param thetabar    var7
#' @param cache       var8
#' @param cbars_small var9  
#' @param cbars_j     var10
#' @param wt          var11
#' @param rss_min_global var12
#' @param dispersion  var13
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
#'     Placeholder: explain how lg_prob_factor, lmc1, and lmc2 are derived and why UB3A ≥ 0.
#'   }
#'
#'   \item{UB3B (dispersion-axis envelope surplus)}{
#'     Placeholder: explain how lm_log1, lm_log2, and max_New_LL_UB are used and why UB3B ≥ 0.
#'   }
#'
#' }
#'
#' Together, these components define
#' \deqn{ test = test1 - UB2 - UB3A - UB3B, }
#' with \eqn{test1 \le 0} and each UB term \eqn{\ge 0}, ensuring the accept/reject
#' procedure is valid and unbiased.
#' @seealso \code{\link{EnvelopeBuild}}, \code{\link{glmb}}, \code{\link{glmbfamfunc}}
#' @export
#' @rdname dispenvelopes
#' @order 1



EnvelopeDispersionBuild <- function(
    Env, Shape, Rate, P, y, x, alpha,
    n_obs, RSS_post, RSS_ML,
    max_disp_perc = 0.99,
    disp_lower = NULL, disp_upper = NULL,
    verbose = FALSE
) {

  # Step 1: Posterior Gamma parameters (precision prior)
  shape2 <- Shape + n_obs / 2
  rate3  <- Rate  + RSS_post / 2
  

  # Step 2: Dispersion bounds (on sigma^2)
  if (is.null(disp_lower) || is.null(disp_upper)) {
    low <- 1 / qgamma(max_disp_perc,     shape2, rate3)
    upp <- 1 / qgamma(1 - max_disp_perc, shape2, rate3)
  } else {
    # Validation checks
    if (!is.numeric(disp_lower) || !is.numeric(disp_upper)) {
      stop("Arguments 'disp_lower' and 'disp_upper' must be numeric.")
    }
    if (disp_lower <= 0 || disp_upper <= 0) {
      stop("Arguments 'disp_lower' and 'disp_upper' must be positive.")
    }
    if (disp_upper <= disp_lower) {
      stop("Argument 'disp_upper' must be strictly greater than 'disp_lower'.")
    }
    low=disp_lower
    upp=disp_upper
  }
  
  
  #low <- 1 / qgamma(max_disp_perc,     shape2, rate3)
  #upp <- 1 / qgamma(1 - max_disp_perc, shape2, rate3)
  
  
  
  
  # Step 3: Extract envelope faces
  cbars     <- Env$cbars
  thetabars <- Env$thetabars
  logP1     <- Env$logP
  gs        <- nrow(cbars)
  
  # Step 4: Base face constants
  # thetabar at dispstar
  
  thetabar_const_base <- thetabar_const(P, cbars, thetabars)
  
  # Step 5: initial anchor (posterior mean, optional)
  ## Is there a potential issue with this dispstar? should this be rate 2 to be consistent with how cbars and 
  ## rate3=Rate+RSS_Post - Outside the function, 
  ##     rate2=rate + RSS_Post2/2   
  ## in the call to Envelpe dispersionbuild,    RSS_post   = RSS_Post2, so rate3=rate2 (outside function) I think
  #dispstar=rate2/(shape2-1)  
  
  
    
  dispstar <- rate3 / (shape2 - 1)
  
  
  cat("[DEBUG] Entering EnvBuildLinBound  \n")
  
  # Step 6: Face slopes at dispstar
  New_LL_Slope <- EnvBuildLinBound(thetabars, cbars, y, x, P, alpha, dispstar)
  
  cat("[DEBUG] Exiting EnvBuildLinBound  \n")

  ## These lines bound the 
  
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
  
    
  cat("[DEBUG] Entering Loop inside EnvelopeDispersionBuild  \n")
  

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
  
  cat("[DEBUG] Exited Loop inside EnvelopeDispersionBuild  \n")
  
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



#' @export 
#' @rdname dispenvelopes
#' @order 2

EnvBuildLinBound<-function(thetabars,cbars,y,x2,P2,alpha,dispstar){
  
  # gs=nrow(cbars)
  # n_vars=ncol(cbars)
  # 
  # New_LL_Slope_test2=c(1:gs)
  # New_LL_Slope_test3=c(1:gs)
  
  XtX   <- crossprod(x2)
  rhs   <- crossprod(x2, y - alpha)
  M     <- XtX + dispstar * P2
  Minv  <- solve(M)
  H1    <- -Minv %*% P2 %*% Minv
  
  V <- -(thetabars %*% P2) + cbars          # gs x p
  Minv_cbars <- t(Minv %*% t(cbars))        # gs x p
  term1 <- rowSums(V * Minv_cbars)
  
  # replicate rhs across gs columns
  rhs_mat <- matrix(rhs, nrow = length(rhs), ncol = nrow(cbars))
  H1_rhs  <- t(H1 %*% (rhs_mat + dispstar * t(cbars)))  # gs × p
  
  term2 <- rowSums(V * H1_rhs)
  
  New_LL_Slope <- term1 + term2
  
  return(New_LL_Slope)
  # for(j in 1:gs){
  #   
  #   cbars_temp=as.matrix(cbars[j,1:n_vars],ncol=1)
  #   thetabars_temp=as.matrix(thetabars[j,1:n_vars],ncol=1)
  #   New_LL_Slope_test2[j]=(-t(thetabars_temp)%*%P2+t(cbars_temp))%*%solve(t(x2)%*%x2+dispstar*P2)%*%cbars_temp
  #   
  #   H1=-solve(t(x2)%*%x2+dispstar*P2)%*%P2%*%solve(t(x2)%*%x2+dispstar*P2)
  #   New_LL_Slope_test3[j]=New_LL_Slope_test2[j]+(-t(thetabars_temp)%*%P2+t(cbars_temp))%*%H1%*%(t(x2)%*%(y-alpha)+dispstar*cbars_temp)
  #   
  # }
  # 
  # return(New_LL_Slope_test3)
  
}

#' @export 
#' @rdname dispenvelopes
#' @order 3




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



#' @export
#' @rdname dispenvelopes
#' @order 4

Inv_f3_with_disp <- function(cache, dispersion, cbars_small) {
  .Call(`_glmbayes_Inv_f3_with_disp`, cache, dispersion, cbars_small)
}

#' @export
#' @rdname dispenvelopes
#' @order 5


UB2 <- function(dispersion, cache, cbars_j, y, x, alpha, wt, rss_min_global) {
  .Call(`_glmbayes_UB2`, dispersion, cache, cbars_j, y, x, alpha, wt, rss_min_global)
}

#' @export
#' @rdname dispenvelopes
#' @order 6

rss_face_at_disp_export <- function(dispersion, cache, cbars_j, y, x, alpha, wt) {
  .Call(`_glmbayes_rss_face_at_disp_export`, dispersion, cache, cbars_j, y, x, alpha, wt)
}
