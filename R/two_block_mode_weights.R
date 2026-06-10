## Likelihood precision (IRLS/Fisher weights) at the random-effects posterior
## mode, for the local-Gaussian heuristic extension of two_block_rate() to
## non-Gaussian families.
##
## The per-observation weight is the (expected) negative second derivative of
## that observation's log-likelihood with respect to its own linear predictor:
##
##   w_i = wt_i * mu'(eta_i)^2 / (V(mu_i) * dispersion)
##
## evaluated at eta_i = offset_i + Z_i b_mode[j(i), ].  For canonical links
## (gaussian-identity, poisson-log, binomial-logit) this is the exact observed
## Hessian; otherwise it is the expected (Fisher) information.  The weights
## are computed generically from the stats family object (mu.eta, variance),
## which reproduces the family-specific glmbfamfunc()$f7 curvature weights on
## the branches where f7 is correct and avoids its copy-pasted logistic
## weights in the probit / cloglog / Gamma-log branches.  The matrix-level
## precision with respect to b_j is the assembly B_j^lik = Z_j' W_j Z_j.

#' Likelihood precision weights at the posterior mode
#'
#' Evaluates the per-observation likelihood precisions (IRLS/Fisher weights)
#' of a GLMM Block 1 model at a supplied random-effects value -- typically
#' the joint posterior mode from \code{\link{glmerb_posterior_mode}} -- and
#' assembles the corresponding per-group likelihood precision blocks
#' \eqn{Z_j^\top W_j Z_j}.
#'
#' The weight for observation \eqn{i} is
#' \deqn{w_i = \frac{wt_i \, \mu'(\eta_i)^2}{V(\mu_i)\,\phi}, \qquad
#'       \eta_i = \mathrm{offset}_i + z_i^\top b_{j(i)},}
#' the (expected) negative curvature of its log-likelihood with respect to
#' its own linear predictor.  For canonical links (gaussian-identity,
#' poisson-log, binomial-logit) this equals the exact observed Hessian at
#' \code{b_mode}; for other links it is the expected (Fisher) information.
#' Examples: gaussian \eqn{w_i = wt_i/\phi}; poisson-log
#' \eqn{w_i = wt_i\,\lambda_i}; binomial-logit \eqn{w_i = wt_i\,p_i(1-p_i)};
#' binomial-probit \eqn{w_i = wt_i\,\phi(\eta_i)^2/[p_i(1-p_i)]}.
#'
#' The returned \code{weights} vector is designed to be passed directly to
#' \code{\link{two_block_rate}(weights = )} to obtain the local-Gaussian
#' heuristic convergence rate for non-Gaussian families.  Note that for
#' non-Gaussian responses the joint posterior is not normal, so rates and TV
#' bounds derived from these weights are a heuristic approximation, not a
#' theorem.
#'
#' @param x Level-1 RE design matrix \code{Z} (\code{l2 x p_re}), as passed
#'   to \code{\link{two_block_rNormal_reg}}.
#' @param block Grouping factor or block partition of length \code{l2}.
#' @param b_mode \code{J x p_re} matrix of random-effect values at which to
#'   evaluate the curvature (rows aligned to \code{group_levels}), e.g.
#'   \code{glmerb_posterior_mode(...)$b_mean}.
#' @param family Response \code{\link[stats]{family}} object (default
#'   \code{gaussian()}).
#' @param wt Prior weights (length 1 or \code{l2}; e.g. binomial trial
#'   counts).  Default 1.
#' @param offset Linear-predictor offset (length 1 or \code{l2}).  Default 0.
#' @param dispersion Dispersion \eqn{\phi}.  Required for \code{gaussian()}
#'   and \code{Gamma} families; defaults to 1 for poisson and binomial.
#' @param group_levels Character vector defining group order (default
#'   \code{levels(block)}); must match the rows of \code{b_mode}.
#' @return Object of class \code{"two_block_mode_weights"}: a list with
#'   \describe{
#'     \item{\code{weights}}{Length-\code{l2} per-observation precisions
#'       \eqn{w_i}; pass as \code{two_block_rate(weights = )}.}
#'     \item{\code{eta}, \code{mu}}{Linear predictor and fitted means at
#'       \code{b_mode}.}
#'     \item{\code{B_lik}}{Named list (per group) of \code{p_re x p_re}
#'       likelihood precision blocks \eqn{Z_j^\top W_j Z_j} (the likelihood
#'       part of \eqn{B_j} before adding the Block 1 prior precision).}
#'     \item{\code{info_total}}{\code{p_re x p_re} total \eqn{Z^\top W Z}.}
#'     \item{\code{family}, \code{link}, \code{dispersion},
#'       \code{group_levels}, \code{b_mode}, \code{call}}{Inputs echoed for
#'       reference.}
#'   }
#' @seealso \code{\link{two_block_rate}}, \code{\link{glmerb_posterior_mode}},
#'   \code{\link{glmbfamfunc}}
#' @family simfuncs
#' @export
two_block_mode_weights <- function(x,
                                   block,
                                   b_mode,
                                   family = gaussian(),
                                   wt = 1,
                                   offset = 0,
                                   dispersion = NULL,
                                   group_levels = levels(block)) {
  cl <- match.call()

  x <- as.matrix(x)
  l2 <- nrow(x)
  p_re <- ncol(x)
  re_names <- colnames(x)
  if (is.null(re_names) || length(re_names) != p_re) {
    re_names <- paste0("RE", seq_len(p_re))
    colnames(x) <- re_names
  }

  family <- .two_block_normalize_family(family)

  ## --- dispersion ------------------------------------------------------------
  needs_disp <- family$family %in% c("gaussian", "Gamma")
  if (is.null(dispersion)) {
    if (needs_disp) {
      stop(
        "'dispersion' is required for family \"", family$family, "\".",
        call. = FALSE
      )
    }
    dispersion <- 1
  }
  dispersion <- as.numeric(dispersion)[1L]
  if (!is.finite(dispersion) || dispersion <= 0) {
    stop("'dispersion' must be a single positive number.", call. = FALSE)
  }

  ## --- group partition (same path as two_block_rate) --------------------------
  block_info <- normalize_block(block, l2)
  if (is.null(group_levels)) {
    group_levels <- block_info$ids
  }
  group_levels <- as.character(group_levels)
  idx_map <- match(group_levels, as.character(block_info$ids))
  if (anyNA(idx_map)) {
    stop(
      "group_levels not found in block ids: ",
      paste(group_levels[is.na(idx_map)], collapse = ", "),
      call. = FALSE
    )
  }
  row_idx <- block_info$rows[idx_map]
  J <- length(group_levels)

  ## --- b_mode aligned to group_levels ----------------------------------------
  b_mode <- as.matrix(b_mode)
  if (ncol(b_mode) != p_re) {
    stop("ncol(b_mode) must equal ncol(x) = ", p_re, ".", call. = FALSE)
  }
  rn <- rownames(b_mode)
  if (!is.null(rn)) {
    if (!all(group_levels %in% rn)) {
      stop("b_mode is missing rows for some group levels.", call. = FALSE)
    }
    b_mode <- b_mode[group_levels, , drop = FALSE]
  } else if (nrow(b_mode) != J) {
    stop("nrow(b_mode) must equal length(group_levels) = ", J, ".",
         call. = FALSE)
  }

  ## --- wt / offset -------------------------------------------------------------
  wt <- as.numeric(wt)
  if (length(wt) == 1L) wt <- rep(wt, l2)
  if (length(wt) != l2) {
    stop("length(wt) must be 1 or nrow(x).", call. = FALSE)
  }
  if (any(!is.finite(wt)) || any(wt < 0)) {
    stop("'wt' must be finite and non-negative.", call. = FALSE)
  }
  offset <- as.numeric(offset)
  if (length(offset) == 1L) offset <- rep(offset, l2)
  if (length(offset) != l2) {
    stop("length(offset) must be 1 or nrow(x).", call. = FALSE)
  }

  ## --- curvature weights at b_mode --------------------------------------------
  ## eta_i = offset_i + z_i' b_mode[j(i), ], evaluated group by group.
  eta <- numeric(l2)
  for (j in seq_len(J)) {
    rows <- row_idx[[j]]
    eta[rows] <- offset[rows] +
      as.vector(x[rows, , drop = FALSE] %*% b_mode[j, ])
  }
  mu <- family$linkinv(eta)
  w <- wt * family$mu.eta(eta)^2 / (family$variance(mu) * dispersion)
  if (any(!is.finite(w))) {
    stop(
      "non-finite curvature weights at b_mode (fitted means on the ",
      "boundary of the parameter space?).",
      call. = FALSE
    )
  }

  ## --- per-group and total likelihood precision blocks ------------------------
  B_lik <- vector("list", J)
  info_total <- matrix(0, p_re, p_re, dimnames = list(re_names, re_names))
  for (j in seq_len(J)) {
    rows <- row_idx[[j]]
    Z_j <- x[rows, , drop = FALSE]
    B_j <- crossprod(Z_j, Z_j * w[rows])
    B_j <- 0.5 * (B_j + t(B_j))
    dimnames(B_j) <- list(re_names, re_names)
    B_lik[[j]] <- B_j
    info_total <- info_total + B_j
  }
  names(B_lik) <- group_levels

  structure(
    list(
      weights = w,
      eta = eta,
      mu = mu,
      B_lik = B_lik,
      info_total = info_total,
      family = family$family,
      link = family$link,
      dispersion = dispersion,
      group_levels = group_levels,
      b_mode = b_mode,
      call = cl
    ),
    class = "two_block_mode_weights"
  )
}

#' Print method for two_block_mode_weights objects
#'
#' @param x Object of class \code{"two_block_mode_weights"}.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.two_block_mode_weights <- function(x, ...) {
  w <- x$weights
  cat("Likelihood precision weights at posterior mode\n")
  cat(sprintf("  family: %s(link = %s),  dispersion = %g\n",
              x$family, x$link, x$dispersion))
  cat(sprintf("  n_obs = %d,  groups J = %d,  p_re = %d\n",
              length(w), length(x$group_levels), ncol(x$b_mode)))
  cat(sprintf("  weights: min = %.4g, median = %.4g, max = %.4g\n",
              min(w), stats::median(w), max(w)))
  if (!identical(x$family, "gaussian")) {
    cat("  (non-Gaussian: local-Gaussian heuristic input for two_block_rate)\n")
  }
  invisible(x)
}
