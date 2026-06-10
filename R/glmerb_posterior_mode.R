#' Joint posterior mode of the two-block GLMM (Gaussian case)
#'
#' Finds the joint posterior mode of the two-block model using an \emph{iterated
#' conditional modes} (ICM) algorithm.  When the response model is Gaussian and
#' variance components are fixed, the joint posterior is multivariate normal so
#' the posterior mode equals the posterior mean; in that case
#' \code{glmerb_posterior_mode()} returns the same result as
#' \code{\link{lmerb_posterior_mean}}.
#'
#' @details
#' \strong{Algorithm.}
#' For any jointly Gaussian distribution the conditional mean of each block is
#' an affine function of the other block's value.  ICM alternates between the
#' two closed-form conditional mean updates:
#'
#' \describe{
#'   \item{Block 1 mode}{For each group \eqn{j}, the conditional posterior mode
#'     of \eqn{b_j} given \eqn{\gamma} is obtained via \code{\link{rglmb}} with
#'     \code{n = 1L} and a \code{\link{dNormal}} prior with mean \eqn{\mu_j}
#'     and covariance \eqn{\Sigma_b}.
#'     The mode is read from \code{coef.mode}.  When \code{family = gaussian()},
#'     this matches the closed-form normal update used in
#'     \code{\link{lmerb_posterior_mean}}.}
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
#' @param family A \code{\link[stats]{family}} object. Defaults to
#'   \code{gaussian()}. Reserved for future non-Gaussian response models.
#' @param measurement_prior_list List with \code{Sigma_ranef} and
#'   \code{prior_list}.  \code{dispersion_ranef} (\eqn{\sigma^2}) is required
#'   for \code{family = gaussian()} and omitted otherwise.  Each
#'   \code{prior_list[[k]]} must contain \code{mu_fixef},
#'   \code{Sigma_fixef}, and \code{dispersion_fixef}.
#' @param tol Convergence tolerance on the \eqn{\ell_\infty} change in
#'   \code{fixef} between successive iterations.  Default \code{1e-10}.
#' @param maxit Maximum number of ICM iterations.  Default \code{200L}.
#'
#' @return A list with components \code{fixef}, \code{b_mean}, \code{converged},
#'   \code{iterations}, and \code{delta}.
#'
#' @seealso \code{\link{lmerb_posterior_mean}}, \code{\link{rglmb}},
#'   \code{\link{two_block_rNormal_reg}}, \code{\link{build_mu_all}}
#' @export
glmerb_posterior_mode <- function(design,
                                  family = gaussian(),
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

  sigma2   <- measurement_prior_list$dispersion_ranef
  Sigma_b  <- measurement_prior_list$Sigma_ranef
  is_gaussian <- identical(family$family, "gaussian")
  if (is_gaussian && is.null(sigma2)) {
    stop(
      "'measurement_prior_list' must contain 'dispersion_ranef' ",
      "when family = gaussian().",
      call. = FALSE
    )
  }
  if (!is_gaussian && !is.null(sigma2)) {
    stop(
      "'measurement_prior_list$dispersion_ranef' must be NULL ",
      "for non-Gaussian families.",
      call. = FALSE
    )
  }

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
      lev  <- group_levels[jj]
      rows <- which(g_chr == lev)
      y_j  <- design$y[rows]
      Z_j  <- design$Z[rows, , drop = FALSE]
      mu_j <- mu_all[, jj]

      pf_j <- if (is.null(sigma2)) {
        dNormal(mu = mu_j, Sigma = Sigma_b)
      } else {
        dNormal(mu = mu_j, Sigma = Sigma_b, dispersion = sigma2)
      }
      fit_j <- rglmb(
        n       = 1L,
        y       = y_j,
        x       = Z_j,
        family  = family,
        pfamily = pf_j,
        verbose = FALSE
      )
      b_mean[jj, ] <- fit_j$coef.mode
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
      "glmerb_posterior_mode() did not converge in ", maxit, " iterations ",
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
