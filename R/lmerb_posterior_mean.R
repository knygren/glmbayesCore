#' Joint posterior mean of the two-block Gaussian model
#'
#' Finds the joint posterior mean (= joint mode, since the posterior is exactly
#' multivariate normal when variance components are fixed) of the two-block
#' model sampled by \code{\link{two_block_rNormal_reg}}, using an \emph{iterated
#' conditional means} (ICM) algorithm.
#'
#' @details
#' \strong{Algorithm.}
#' For any jointly Gaussian distribution the conditional mean of each block is
#' an affine function of the other block's value.  ICM alternates between the
#' two closed-form conditional mean updates:
#'
#' \describe{
#'   \item{Block 1 mean}{For each group \eqn{j}:
#'     \deqn{
#'       E[b_j \mid \gamma] =
#'       \bigl(Z_j^\top Z_j / \sigma^2 + P_b\bigr)^{-1}
#'       \bigl(Z_j^\top y_j / \sigma^2 + P_b \,\mu_j(\gamma)\bigr)
#'     }
#'     where \eqn{P_b = \Sigma_b^{-1}} and
#'     \eqn{\mu_j(\gamma)} is the Block 2 prior mean from
#'     \code{\link{build_mu_all}}.
#'   }
#'   \item{Block 2 mean}{For each RE component \eqn{k}:
#'     \deqn{
#'       E[\gamma_k \mid b_k] =
#'       \bigl(X_k^\top X_k / \tau^2_k + P_{\gamma_k}\bigr)^{-1}
#'       \bigl(X_k^\top b_k / \tau^2_k + P_{\gamma_k} \mu_{\gamma_k}\bigr)
#'     }
#'     where \eqn{b_k} is the \eqn{k}-th column of the current Block 1 mean
#'     matrix, \eqn{X_k = } \code{design$X_hyper[[k]]},
#'     \eqn{\tau^2_k = } \code{dispersion_fixef}, and
#'     \eqn{P_{\gamma_k} = \Sigma_{\gamma_k}^{-1}} from \code{prior_list}.
#'   }
#' }
#'
#' @param design Design list with \code{y}, \code{Z}, \code{groups},
#'   \code{X_hyper}, and \code{re_coef_names}.
#' @param measurement_prior_list List with \code{dispersion_ranef}
#'   (Gaussian only), \code{Sigma_ranef}, and \code{prior_list}. Each \code{prior_list[[k]]}
#'   must contain \code{mu_fixef}, \code{Sigma_fixef}, and
#'   \code{dispersion_fixef}.
#' @param tol Convergence tolerance on the \eqn{\ell_\infty} change in
#'   \code{fixef} between successive iterations.  Default \code{1e-10}.
#' @param maxit Maximum number of ICM iterations.  Default \code{200L}.
#'
#' @return A list with components \code{fixef}, \code{b_mean}, \code{converged},
#'   \code{iterations}, and \code{delta}.
#'
#' @seealso \code{\link{two_block_rNormal_reg}}, \code{\link{build_mu_all}}
#' @export
lmerb_posterior_mean <- function(design,
                                 measurement_prior_list,
                                 tol   = 1e-10,
                                 maxit = 200L) {

  .lmerb_validate_design(design)
  .lmerb_validate_measurement_prior_list(measurement_prior_list)

  if (is.null(design$y) || is.null(design$Z)) {
    stop("'design' must contain 'y' and 'Z'.", call. = FALSE)
  }

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)
  J            <- length(group_levels)
  p_re         <- length(re_names)
  g_chr        <- as.character(design$groups)

  sigma2 <- measurement_prior_list$dispersion_ranef
  if (is.null(sigma2)) {
    stop(
      "'measurement_prior_list' must contain 'dispersion_ranef' ",
      "for lmerb_posterior_mean().",
      call. = FALSE
    )
  }
  P_b    <- solve(measurement_prior_list$Sigma_ranef)

  P_gamma  <- stats::setNames(
    lapply(re_names, function(k) {
      solve(measurement_prior_list$prior_list[[k]]$Sigma_fixef)
    }),
    re_names
  )
  mu_gamma <- stats::setNames(
    lapply(re_names, function(k) measurement_prior_list$prior_list[[k]]$mu_fixef),
    re_names
  )
  tau2 <- stats::setNames(
    lapply(re_names, function(k) measurement_prior_list$prior_list[[k]]$dispersion_fixef),
    re_names
  )

  ZtZ_scaled <- vector("list", J)
  Zty_scaled <- vector("list", J)
  names(ZtZ_scaled) <- names(Zty_scaled) <- group_levels

  for (lev in group_levels) {
    rows <- which(g_chr == lev)
    Z_j  <- design$Z[rows, , drop = FALSE]
    y_j  <- design$y[rows]
    ZtZ_scaled[[lev]] <- crossprod(Z_j) / sigma2
    Zty_scaled[[lev]] <- crossprod(Z_j, y_j) / sigma2
  }

  fixef <- lapply(measurement_prior_list$prior_list, `[[`, "mu_fixef")
  names(fixef) <- re_names

  b_mean <- matrix(
    0.0, nrow = J, ncol = p_re,
    dimnames = list(group_levels, re_names)
  )

  converged <- FALSE
  delta     <- NA_real_

  for (iter in seq_len(maxit)) {

    mu_all <- as.matrix(build_mu_all(design, fixef)$mu_all)

    for (jj in seq_len(J)) {
      lev      <- group_levels[jj]
      mu_j     <- mu_all[, jj]
      post_P_j <- ZtZ_scaled[[lev]] + P_b
      post_v_j <- Zty_scaled[[lev]] + P_b %*% mu_j
      b_mean[jj, ] <- solve(post_P_j, post_v_j)
    }

    fixef_new <- vector("list", p_re)
    names(fixef_new) <- re_names

    for (k in re_names) {
      X_k      <- design$X_hyper[[k]]
      b_k      <- b_mean[, k]
      tau2_k   <- tau2[[k]]
      P_gam_k  <- P_gamma[[k]]
      mu_gam_k <- mu_gamma[[k]]

      post_P_k  <- crossprod(X_k) / tau2_k + P_gam_k
      post_v_k  <- crossprod(X_k, b_k) / tau2_k + P_gam_k %*% mu_gam_k
      gam_k <- as.vector(solve(post_P_k, post_v_k))
      names(gam_k) <- colnames(X_k)
      fixef_new[[k]] <- gam_k
    }

    delta <- max(vapply(re_names, function(k) {
      max(abs(fixef_new[[k]] - fixef[[k]]))
    }, numeric(1L)))

    fixef <- fixef_new

    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      "lmerb_posterior_mean() did not converge in ", maxit, " iterations ",
      "(final delta = ", signif(delta, 3L), "). ",
      "Consider increasing 'maxit' or checking model identifiability.",
      call. = FALSE
    )
  }

  list(
    fixef      = fixef,
    b_mean     = b_mean,
    converged  = converged,
    iterations = iter,
    delta      = delta
  )
}

#' @noRd
.lmerb_validate_measurement_prior_list <- function(mpl) {
  if (!is.list(mpl)) {
    stop("'measurement_prior_list' must be a list.", call. = FALSE)
  }
  for (nm in c("Sigma_ranef", "prior_list")) {
    if (is.null(mpl[[nm]])) {
      stop(
        "'measurement_prior_list' must contain '", nm, "'.",
        call. = FALSE
      )
    }
  }
  invisible(mpl)
}
