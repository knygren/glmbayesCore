
#' Directional Tail Diagnostic
#'
#' Computes the directional tail probability based on posterior draws and prior mean,
#' using whitening transformation and projection onto the direction of disagreement.
#' This diagnostic identifies directional disagreement between posterior and prior,
#' and is especially useful for visualizing rejection regions in whitened space.
#'
#' @param fit A fitted model object of class 'glmb' or 'lmb'
#' @return An object of class 'directional_tail' containing:
#'   \item{mahalanobis_shift}{Measures the standardized Mahalanobis distance between the posterior and prior means, 
#' using posterior precision for scaling. In the Gaussian case, this directly determines the directional tail probability via Phi(-||w||).}
#'   \item{p_directional}{Directional tail probability (proportion of draws in the direction of disagreement)}
#'   \item{delta}{Mean shift in whitened space}
#'   \item{draws}{List containing whitened draws, raw draws, and tail flags}
#' @details Whitening is performed using the posterior precision matrix.
#' The direction vector is computed as the mean shift in whitened space.
#' Tail probability is the proportion of draws with negative projection onto this direction.
#' @example inst/examples/Ex_directional_tail.R
#' @export
#' @keywords diagnostic geometry Bayesian




directional_tail <- function(fit) {
  main_class <- class(fit)[1]
  
  ## 1) Extract posterior draws and prior mean
  B      <- as.matrix(fit$coefficients)       # M * p
  mu0    <- as.numeric(fit$Prior$mean)        # length p
  V0     <- fit$Prior$Variance                # p * p
  Sigma  <- as.matrix(vcov(fit))              # p * p
  
  ## 2) Compute likelihood precision
  if (main_class == "glmb") {
    f7     <- fit$famfunc$f7
    bstar  <- fit$coef.mode
##    bstar <- matrix(bstar, ncol = 1)
    y      <- fit$y
    x      <- fit$x
    P0     <- solve(V0)
    alpha  <- tryCatch(fit$offset, error = function(e) fit$glm$offset)
    if (is.null(alpha)) alpha <- 0
    wt     <- fit$prior.weights
    

    
    Prec_lik <- f7(bstar, y, x, mu0, P0, alpha, wt)
  } else if (main_class == "lmb") {
    Prec_lik <- solve(vcov(fit$lm))
  } else {
    stop("Unsupported model class")
  }
  
  ## 3) Combine with prior precision
  Prec_post <- Prec_lik + chol2inv(chol(V0))

  ## 4) Whiten draws relative to prior mean
  P <- Prec_post
  eig <- eigen(P, symmetric = TRUE)
  D_eig <- diag(sqrt(eig$values), nrow = length(eig$values))
  W      <- eig$vectors %*% D_eig %*% t(eig$vectors)
  Z      <- t(W %*% t(sweep(B, 2, mu0)))       # M * p
  
  
  ## 5) Direction vector from prior to posterior (in whitened space)
  delta  <- colMeans(Z)
  w      <- delta
  
  
  ## 6) Project draws and compute tail probability
  proj   <- as.vector(Z %*% w)
  flag   <- proj <= 0
  p_tail <- mean(flag)
  
  ## 7) Return simplified output
  out <- list(
    mahalanobis_shift = sqrt(sum(delta^2)),
    p_directional   = p_tail,
    delta           = delta,
    draws = list(
      is_tail = flag,
      Z       = Z,
      B       = B
    )
  )
  
  
  class(out) <- "directional_tail"
  return(out)
}

#' @export

print.directional_tail <- function(x, ...) {
  # Determine number of draws
  n_draws <- nrow(x$draws$Z)
  
  cat("Bayesian Estimates Based on", n_draws, "iid draws\n")
  cat("--------------------------------------------------\n")
  cat("Standardized Prior-Posterior Mahalanobis distance, and\n")
  cat("associated tail probability (P[delta^T * Z <= 0]):\n\n")
  cat("  Mahalanobis Prior-Posterior Distance    :", formatC(x$mahalanobis_shift, digits = 4, format = "f"), "\n")
  cat("  Associated Directional Tail Probability :", formatC(x$p_directional, digits = 4, format = "f"), "\n\n")
  
  # Compute MC error for delta on the fly
  Z <- x$draws$Z
  mcse_delta <- apply(Z, 2, function(zj) sd(zj) / sqrt(length(zj)))
  
  cat("Standardized Prior-posterior shifts mcses: \n")
  
    # Create and print table of delta and MCSE
  delta_table <- data.frame(
    Delta = round(x$delta, 4),
    MCSE  = round(mcse_delta, 4),
    row.names = colnames(Z)
  )
  print(delta_table)
  
  cat("\nPosterior draws available:\n")
  cat("  - `x$draws$Z` for whitened draws\n")
  cat("  - `x$draws$B` for raw coefficient draws\n")
  cat("  - `x$draws$is_tail` to flag tail draws\n")
  
  invisible(x)
}