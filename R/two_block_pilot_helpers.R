#' Print pilot-stage diagnostics between pilot and main sampling (UB path)
#' @noRd
.two_block_print_pilot_stage_diagnostics <- function(
    n_pilot,
    n_main,
    pilot_ub,
    rate_calibration,
    m_convergence_used) {

  if (!is.null(pilot_ub) && !is.null(rate_calibration)) {
    .fmt_eigs <- function(ev) paste(sprintf("%.4f", ev), collapse = ", ")
    cat(sprintf(
      "--- glmerb: post-pilot convergence bounds (%d pilot draws) ---\n    ML estimate (local-Gaussian at mode):    lambda* = %.4f, m_min = %d, eigenvalues = [%s]\n    Pilot upper bound (per-eig max, #%d/%d): lambda* = %.4f, m_min = %d, eigenvalues = [%s]\n    => using m_convergence = %d ---\n\n",
      n_pilot,
      rate_calibration$lambda_star,
      rate_calibration$m_min,
      .fmt_eigs(rate_calibration$eigenvalues),
      pilot_ub$i_max_rate, n_pilot,
      pilot_ub$rate_upper$lambda_star,
      pilot_ub$m_min_upper,
      .fmt_eigs(pilot_ub$max_eigenvalues),
      m_convergence_used
    ))
  }

  cat(sprintf(
    "--- glmerb: pilot complete; main stage (%d independent chains from pilot mean; m_convergence = %d) ---\n\n",
    n_main, m_convergence_used
  ))
  flush.console()
}

#' Rename engine fields to public group/pop (\code{groupef}/\code{popef}) names
#'
#' Internal engine keys (\code{fixef_draws}, \code{coefficients}, \ldots) are
#' mapped to package notation: Block~1 \eqn{\beta_j} as \code{groupef*},
#' Block~2 \eqn{\gamma} as \code{popef*}. Unused engine leftovers are dropped.
#' @noRd
.two_block_as_staged_names <- function(x, fixef_mode, fixef_init) {
  renames <- c(
    fixef_draws            = "popef",
    coefficients           = "groupef",
    dispersion_fixef_draws = "popef.dispersion",
    iters_fixef_draws      = "popef.iters",
    iters_ranef_draws      = "groupef.iters"
  )
  for (old_nm in names(renames)) {
    if (!is.null(x[[old_nm]])) {
      x[[renames[[old_nm]]]] <- x[[old_nm]]
      x[[old_nm]] <- NULL
    }
  }
  drop <- c(
    "fixef_last", "b_last", "mu_all_last", "groupef.names",
    "group_levels", "fixef_start", "sampling", "m_convergence",
    "fixef.last", "coef.last", "fixef.mu", "coef.names"
  )
  for (nm in drop) {
    x[[nm]] <- NULL
  }
  x$popef.mode <- fixef_mode
  x$popef.init <- fixef_init
  x
}

#' Hotelling chi-squared test: pilot fixef mean vs start
#' @noRd
.two_block_pilot_chisq_test <- function(fixef_draws, re_names, fixef_mode,
                                          n_pilot) {
  X_pilot <- do.call(cbind, lapply(re_names, function(k) fixef_draws[[k]]))
  pnames <- unlist(lapply(re_names, function(k) {
    paste0(k, "::", colnames(fixef_draws[[k]]))
  }))
  colnames(X_pilot) <- pnames
  mu_pilot <- colMeans(X_pilot)
  mode_vec <- unlist(lapply(re_names, function(k) fixef_mode[[k]]))
  names(mode_vec) <- pnames
  d_pm <- mu_pilot - mode_vec
  S_pilot <- stats::cov(X_pilot)
  p_dim2 <- ncol(X_pilot)
  ridge <- 1e-8 * sum(diag(S_pilot)) / p_dim2
  if (!is.finite(ridge) || ridge <= 0) {
    ridge <- 1e-8
  }
  S_inv <- tryCatch(
    solve(S_pilot),
    error = function(e) {
      solve(S_pilot + diag(ridge, p_dim2))
    }
  )
  Q_pm <- as.numeric(n_pilot * t(d_pm) %*% S_inv %*% d_pm)
  p_pm <- stats::pchisq(Q_pm, df = p_dim2, lower.tail = FALSE)
  list(
    Q       = Q_pm,
    df      = p_dim2,
    p_value = p_pm,
    n_pilot = n_pilot
  )
}

#' Column means of fixef draws with names copied from fixef mode
#' @noRd
.two_block_fixef_colmeans <- function(fixef_draws, re_names, fixef_mode) {
  out <- lapply(re_names, function(k) {
    v <- colMeans(fixef_draws[[k]])
    names(v) <- names(fixef_mode[[k]])
    v
  })
  stats::setNames(out, re_names)
}


#' Post-pilot eigenvalue upper bounds from v6 stacked \code{coefficients}
#' @noRd
.two_block_pilot_ub_from_coefficients <- function(
    pilot_coefficients,
    n_pilot,
    re_names,
    group_levels,
    group_name,
    x,
    group,
    x_hyper,
    prior_list,
    pfamily_list,
    family,
    tv_tol,
    dispersion = NULL,
    weights = NULL
) {
  n_grp        <- length(group_levels)
  re_col_names <- re_names
  grp_col_name <- group_name

  lambda_star_vec  <- numeric(n_pilot)
  max_eigenvalues  <- NULL
  rate_upper       <- NULL
  lambda_star_best <- -Inf
  i_max_rate       <- NA_integer_

  for (i in seq_len(n_pilot)) {
    rows_i   <- ((i - 1L) * n_grp + 1L):(i * n_grp)
    block_df <- pilot_coefficients[rows_i, , drop = FALSE]
    if (!is.null(grp_col_name) && grp_col_name %in% colnames(block_df)) {
      ord <- match(group_levels, block_df[[grp_col_name]])
      b_i <- as.matrix(block_df[ord, re_col_names, drop = FALSE])
    } else {
      b_i <- as.matrix(block_df[, re_col_names, drop = FALSE])
    }
    rownames(b_i) <- group_levels
    # Explicit 'weights' (e.g. per-group precision from disp_upper_group)
    # bypasses the mode-dependent curvature weights entirely; for the
    # gaussian() likelihood these are mode-independent anyway (constant
    # across b_i draws), so this is exact, not an approximation.
    weights_i <- if (!is.null(weights)) {
      weights
    } else {
      two_block_mode_weights(
        x            = x,
        group        = group,
        b_mode       = b_i,
        family       = family,
        dispersion   = dispersion,
        group_levels = group_levels
      )$weights
    }
    rate_i <- two_block_rate_from_pfamily_list(
      x                 = x,
      group             = group,
      x_hyper           = x_hyper,
      prior_list_block1 = prior_list,
      pfamily_list      = pfamily_list,
      weights           = weights_i,
      family            = family,
      group_levels      = group_levels,
      warn_slow         = FALSE
    )
    if (is.null(max_eigenvalues)) {
      max_eigenvalues <- rep(-Inf, length(rate_i$eigenvalues))
    }
    lambda_star_vec[i] <- rate_i$lambda_star
    max_eigenvalues    <- pmax(max_eigenvalues, rate_i$eigenvalues)
    if (rate_i$lambda_star > lambda_star_best) {
      lambda_star_best <- rate_i$lambda_star
      rate_upper       <- rate_i
      i_max_rate       <- i
    }
  }

  rate_upper_eig             <- rate_upper
  rate_upper_eig$eigenvalues <- max_eigenvalues
  rate_upper_eig$lambda_star <- max_eigenvalues[1L]

  m_min_upper <- .two_block_cap_inner_sweeps(
    two_block_l_for_tv(rate_upper_eig, tv_tol, method = "theorem3") + 1L
  )

  list(
    rate_upper      = rate_upper_eig,
    m_min_upper     = m_min_upper,
    lambda_star_vec = lambda_star_vec,
    i_max_rate      = i_max_rate,
    max_eigenvalues = max_eigenvalues
  )
}
