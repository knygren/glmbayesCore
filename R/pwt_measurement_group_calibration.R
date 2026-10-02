## Group-specific pwt_measurement calibration targeting a fixed ellipsoid
## violation rate (Section 16.6 exact/truncated criterion). Promoted from
## the investigation-only data-raw/_scratch_group_pwt_measurement_noncentral.R
## -- see inst/omega-ing-marginal-multivariate-t.md Section 9 and
## inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 16.6 for the derivation.

#' Exact truncated-Gamma moments of measurement precision (Section 16.6)
#' @noRd
.two_block_truncated_omega_moments <- function(shape, rate, omega_L, omega_U,
                                                dP0_floor = 1e-8) {
  dP <- function(k) {
    stats::pgamma(omega_U, shape = k, rate = rate) -
      stats::pgamma(omega_L, shape = k, rate = rate)
  }
  dP0 <- dP(shape)
  dP1 <- dP(shape + 1)
  dP2 <- dP(shape + 2)

  bad <- !is.finite(dP0) | dP0 < dP0_floor
  E_omega   <- rep(NA_real_, length(rate))
  Var_omega <- rep(NA_real_, length(rate))

  ok <- !bad
  if (any(ok)) {
    E_omega[ok]  <- (shape / rate[ok]) * (dP1[ok] / dP0[ok])
    E2           <- (shape * (shape + 1) / rate[ok]^2) * (dP2[ok] / dP0[ok])
    Var_omega[ok] <- E2 - E_omega[ok]^2
  }

  list(E_omega = E_omega, Var_omega = Var_omega, dP0 = dP0, ok = ok)
}

#' Per-group \code{pwt_measurement} search against the exact ellipsoid criterion
#'
#' For each group level, simulates the (\code{pwt_measurement}-invariant, per
#' Section 9.3 of \code{inst/omega-ing-marginal-multivariate-t.md})
#' distribution of \eqn{q_j = (\beta_j-\hat\beta_j^{\mathrm{OLS}})' D_j'D_j
#' (\beta_j-\hat\beta_j^{\mathrm{OLS}})} under the noncentral Gaussian full
#' conditional at \code{ing_prior_measurement_group}'s (pass-1)
#' \eqn{\hat\sigma_j^2}, then 1-D \code{uniroot()}s over the per-group prior
#' weight \eqn{w_j} against the EXACT, truncation-aware boundary \eqn{q_j >
#' E_t[\Omega_j]/\mathrm{Var}_t[\Omega_j]} (Section 16.6) to find the
#' smallest \eqn{w_j} driving the predicted violation rate down to
#' \code{alpha_target}. The result is floored at \code{floor_vec} (the
#' caller's already-resolved \code{pwt_measurement}) -- calibration only
#' ever sharpens (increases \eqn{w_j}), never loosens, the caller's own
#' prior weight.
#'
#' Keeps only the exact/truncated criterion; the scratch version's
#' untruncated and central-F-only (no-shrinkage) columns were investigation
#' contrasts only and are not reproduced here.
#'
#' @param fit_ref Reference fit (\code{merMod} or \code{glmmTMB}) used for
#'   \code{gamma_hat} (the Block~2 hyper-regression's fixed-effect point
#'   estimates, i.e. each group's shrinkage-target mean \eqn{\mu_j}).
#' @param design \code{model_setup()} design list (\code{D}, \code{y},
#'   \code{W}, \code{group}, \code{groupef.names}).
#' @param group_levels Character vector of group levels.
#' @param groupef.names Character vector of random-coefficient names.
#' @param group.Sigma Diagonal RE covariance matrix (Block~1's \eqn{\Psi}).
#' @param ing_prior_measurement_group Pass-1 (seed) per-group Block~1
#'   calibration list, as returned by
#'   \code{.lmebayes_calibrate_ing_prior_measurement_group()}; supplies each
#'   group's \code{sigma2_hat} and \code{max_disp_perc}.
#' @param alpha_target Target ellipsoid violation rate, a scalar in
#'   \eqn{(0, 1)} (e.g. \code{0.01}).
#' @param floor_vec Named length-\eqn{J} numeric vector, the per-group
#'   \code{pwt_measurement} floor (pass-1's resolved
#'   \code{meas_group$pwt_measurement}).
#' @param n_sim Monte Carlo draws per group.
#' @param seed RNG seed (calibration is deterministic given the reference
#'   fit/design, for reproducibility across calls).
#' @return \code{list(pwt_measurement = <named length-J vector>, table =
#'   <data.frame>)}.
#' @noRd
.lmebayes_calibrate_pwt_measurement_group <- function(
    fit_ref, design, group_levels, groupef.names, group.Sigma,
    ing_prior_measurement_group, alpha_target, floor_vec,
    n_sim = 200000L, seed = 1L
) {
  p_re <- length(groupef.names)

  fe <- .lmebayes_reference_fixef(fit_ref)
  ## diag(x) with a length-1 numeric x builds an x-by-x identity matrix, not
  ## a 1x1 diagonal matrix containing x -- must pass nrow explicitly to get
  ## the intended 1x1 result when p_re == 1 (a single RE component, e.g.
  ## (1 | group) with no other random slopes).
  Psi_inv <- diag(1 / diag(group.Sigma), nrow = p_re)
  dimnames(Psi_inv) <- dimnames(group.Sigma)

  ## Map hyper-predictor 'col' of RE component 'k' to fe's name: a main
  ## effect for intercept-associated predictors, but the observation-level
  ## interaction term for a non-intercept RE component's own hyper-
  ## predictors -- mirrors Prior_Setup_GLMM()'s own private fe_name_for().
  fe_name_for <- function(k, col) {
    if (identical(k, "(Intercept)") && identical(col, "(Intercept)")) {
      return("(Intercept)")
    }
    if (identical(col, "(Intercept)")) {
      return(if (k %in% names(fe)) k else NA_character_)
    }
    if (identical(k, "(Intercept)")) {
      return(if (col %in% names(fe)) col else NA_character_)
    }
    cand <- c(paste0(col, ":", k), paste0(k, ":", col))
    hit  <- cand[cand %in% names(fe)]
    if (length(hit)) hit[1L] else NA_character_
  }

  gamma_hat_by_component <- stats::setNames(
    lapply(groupef.names, function(k) {
      Wk   <- design$W[[k]]
      cols <- colnames(Wk)
      fe_nm <- vapply(cols, fe_name_for, character(1L), k = k)
      if (anyNA(fe_nm)) {
        stop(
          "'alpha_target_measurement' calibration could not map hyper-",
          "predictor(s) for '", k, "' to the reference fit's fixed effects.",
          call. = FALSE
        )
      }
      stats::setNames(unname(fe[fe_nm]), cols)
    }),
    groupef.names
  )

  set.seed(seed)
  group_chr <- as.character(design$group)

  rows <- lapply(group_levels, function(lev) {
    idx <- which(group_chr == lev)
    D_j <- design$D[idx, groupef.names, drop = FALSE]
    y_j <- design$y[idx]
    n_j <- length(idx)
    DtD <- crossprod(D_j)
    beta_ols <- as.vector(solve(DtD, crossprod(D_j, y_j)))
    RSS_ols  <- sum((y_j - D_j %*% beta_ols)^2)

    ing_j         <- ing_prior_measurement_group[[lev]]
    sigma2_hat    <- ing_j$sigma2_hat
    Omega_hat_j   <- 1 / sigma2_hat
    max_disp_perc <- ing_j$max_disp_perc

    mu_j <- vapply(groupef.names, function(k) {
      Wk <- design$W[[k]]
      as.numeric(Wk[lev, , drop = FALSE] %*% gamma_hat_by_component[[k]])
    }, numeric(1))

    Sigma_j_inv <- Omega_hat_j * DtD + Psi_inv
    Sigma_j     <- solve(Sigma_j_inv)
    beta_bar_j  <- as.vector(
      Sigma_j %*% (Omega_hat_j * crossprod(D_j, y_j) + Psi_inv %*% mu_j)
    )

    ## ONE simulation of q_j's distribution -- invariant to pwt_measurement
    ## (Section 9.3): reused across every candidate w below.
    L     <- chol(Sigma_j)
    Z     <- matrix(stats::rnorm(n_sim * p_re), nrow = n_sim)
    draws <- sweep(Z %*% L, 2, beta_bar_j, "+")
    diffs <- sweep(draws, 2, beta_ols, "-")
    q_draws   <- rowSums((diffs %*% DtD) * diffs)
    ete_draws <- q_draws + RSS_ols

    ## Exact criterion (Section 16.6) as a function of candidate w: a0(w)/
    ## r0(w) mean-match sigma2_hat (fixed, pwt_measurement-invariant) while
    ## n_prior(w) = w/(1-w)*n_j moves them together; disp_lower(w)/
    ## disp_upper(w) are recomputed from the POSTERIOR-shape window (a0(w) +
    ## n_j/2, mean-matched at the same sigma2_hat), matching
    ## Prior_Setup_GLMM()'s own .lmebayes_calibrate_ing_prior_measurement_group()
    ## construction -- so the window this search certifies w against is the
    ## same one that will actually be shipped to the sampler once w_final is
    ## fed back into Prior_Setup_GLMM().
    pct_outside_for_w <- function(w) {
      n_prior_w <- w / (1 - w) * n_j
      a0_w <- (n_prior_w + p_re + 1) / 2
      r0_w <- sigma2_hat * (a0_w - 1)
      shape_post_w <- a0_w + n_j / 2
      rate_post_w  <- sigma2_hat * (shape_post_w - 1)
      win_w <- .lmebayes_ing_prior_quantile_window(
        shape_post_w, rate_post_w, max_disp_perc
      )
      omega_L_w <- 1 / win_w$disp_upper
      omega_U_w <- 1 / win_w$disp_lower
      shape_j_w <- a0_w + n_j / 2
      rate_j_w  <- r0_w + 0.5 * ete_draws
      mom_w <- .two_block_truncated_omega_moments(
        shape_j_w, rate_j_w, omega_L_w, omega_U_w
      )
      n_ok <- sum(mom_w$ok)
      if (n_ok == 0L) {
        return(1) ## conservative: can't certify, treat as all-violating
      }
      sum(mom_w$ok & (q_draws > mom_w$E_omega / mom_w$Var_omega)) / n_ok
    }

    floor_j <- unname(floor_vec[[lev]])
    pct_outside_before <- pct_outside_for_w(floor_j)

    w_lo <- 1e-6
    w_hi <- 0.5 - 1e-6
    pct_lo <- pct_outside_for_w(w_lo)
    pct_hi <- pct_outside_for_w(w_hi)

    if (pct_lo <= alpha_target) {
      ## Already below target with almost no pseudo-observations.
      w_star <- 0
      clipped_at_ceiling <- FALSE
    } else if (pct_hi > alpha_target) {
      ## Even the package's 0.5 per-group ceiling isn't enough.
      w_star <- 0.5
      clipped_at_ceiling <- TRUE
    } else {
      g <- function(w) pct_outside_for_w(w) - alpha_target
      w_star <- stats::uniroot(g, lower = w_lo, upper = w_hi, tol = 1e-4)$root
      clipped_at_ceiling <- FALSE
    }

    w_final     <- max(w_star, floor_j)
    floor_binds <- floor_j > w_star
    pct_outside_after <- pct_outside_for_w(w_final)

    data.frame(
      group               = lev,
      n_j                 = n_j,
      pwt_floor           = floor_j,
      pct_outside_before  = 100 * pct_outside_before,
      w_star              = w_star,
      clipped_at_ceiling  = clipped_at_ceiling,
      w_final             = w_final,
      floor_binds         = floor_binds,
      pct_outside_after   = 100 * pct_outside_after,
      stringsAsFactors    = FALSE
    )
  })

  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  pwt <- stats::setNames(tab$w_final, tab$group)

  list(pwt_measurement = pwt, table = tab)
}
