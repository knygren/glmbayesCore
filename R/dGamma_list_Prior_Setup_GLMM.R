#' @rdname dGamma_list
#' @order 2
#' @param max_disp_perc_measurement \code{NULL} (default), a scalar in
#'   \eqn{(0.5, 1)}, or a named/positional length-\eqn{J} vector (one value
#'   per group level). \code{NULL} reuses each group's stored
#'   \code{disp_lower}/\code{disp_upper}. A scalar or vector recomputes
#'   those bounds per group.
#' @details
#' For \code{\link{Prior_Setup_GLMM}} with a per-group \code{dispformula},
#' each group's \code{dGamma()} uses \code{shape_ING}/\code{rate} from
#' \code{object$group.ing_prior} (Part VI calibration; see
#' \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md}). Truncation bounds are
#' posterior-shape quantiles at that group's \code{max_disp_perc}.
#' @export
#' @method dGamma_list Prior_Setup_GLMM
dGamma_list.Prior_Setup_GLMM <- function(
    object,
    max_disp_perc_measurement = NULL,
    ...
) {
  if (!identical(object$family$family, "gaussian")) {
    stop(
      "dGamma_list() for Prior_Setup_GLMM requires family = gaussian().",
      call. = FALSE
    )
  }

  ing_grp <- object$group.ing_prior
  if (is.null(ing_grp) || !.lmebayes_ing_prior_is_grouped(ing_grp)) {
    grp_nm <- object$design$group_name
    stop(
      "object's group.ing_prior is not per-group; call Prior_Setup_GLMM(",
      "..., dispformula = ~", if (!is.null(grp_nm)) grp_nm else "<group_name>",
      ") on a Gaussian model to calibrate per-group measurement-dispersion ",
      "priors (dispformula = ~1, the default, skips this calibration).",
      call. = FALSE
    )
  }

  group_levels <- names(ing_grp)
  if (is.null(group_levels)) {
    group_levels <- levels(object$design$group)
    names(ing_grp) <- group_levels
  }

  ## NULL: reuse each group's own stored max_disp_perc (no override, no
  ## recompute). Supplied: resolve to a per-group vector (scalar recycled,
  ## or a length-J vector) and recompute per group below.
  mdp_override <- if (!is.null(max_disp_perc_measurement)) {
    .lmebayes_expand_scalar_or_vector(
      max_disp_perc_measurement, group_levels, "max_disp_perc_measurement"
    )
  } else {
    NULL
  }

  diag_rows <- vector("list", length(group_levels))
  out <- stats::setNames(
    lapply(seq_along(group_levels), function(i) {
      lev <- group_levels[[i]]
      g <- ing_grp[[lev]]

      ## The bounds stored on g were computed at g$max_disp_perc; if the
      ## caller asks for a different value for this group, recompute fresh
      ## quantiles of Gamma(shape_ING + n_j/2, rate_post) -- the posterior-
      ## shape window (see Prior_Setup_GLMM()'s group.ing_prior
      ## docs), NOT g$shape_ING/g$rate directly (those are the prior fed to
      ## the sampler, unchanged) -- rather than reusing the stored bounds.
      mdp_lev <- if (!is.null(mdp_override)) {
        unname(mdp_override[[lev]])
      } else {
        g$max_disp_perc
      }
      recompute <- !isTRUE(all.equal(mdp_lev, g$max_disp_perc))

      if (recompute || is.null(g$disp_lower) || is.null(g$disp_upper)) {
        shape_post_lev <- g$shape_ING + g$n_j / 2
        rate_post_lev  <- g$sigma2_hat * (shape_post_lev - 1)
        win <- .lmebayes_ing_prior_quantile_window(
          shape_post_lev, rate_post_lev, mdp_lev
        )
        disp_lower <- win$disp_lower
        disp_upper <- win$disp_upper
      } else {
        disp_lower <- g$disp_lower
        disp_upper <- g$disp_upper
      }

      diag_rows[[i]] <<- data.frame(
        group         = lev,
        n_j           = g$n_j,
        sigma2_hat    = unname(g$sigma2_hat),
        shape_ING     = unname(g$shape_ING),
        rate          = unname(g$rate),
        max_disp_perc = mdp_lev,
        disp_lower    = disp_lower,
        disp_upper    = disp_upper,
        stringsAsFactors = FALSE
      )

      dGamma(
        shape          = g$shape_ING,
        rate           = g$rate,
        beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
        Inv_Dispersion = TRUE,
        max_disp_perc  = mdp_lev,
        disp_lower     = disp_lower,
        disp_upper     = disp_upper
      )
    }),
    group_levels
  )

  window_diagnostics <- do.call(rbind, diag_rows)
  rownames(window_diagnostics) <- NULL
  attr(out, "window_diagnostics")  <- window_diagnostics
  ## Carried forward so lmerb()/glmerb() can reuse this glmmTMB reference fit
  ## as their diagnostic-only group.dispersion.fit instead of re-fitting it; see
  ## lmebayes:::.lmebayes_fit_glmmtmb_dispersion().
  attr(out, "group.dispersion.fit")      <- object$group.dispersion.fit
  attr(out, "calibration_source")  <- object$calibration_source

  ## Flattened view of the SAME four per-group quantities already collected
  ## into window_diagnostics above (shape_ING/rate/disp_lower/disp_upper),
  ## reshaped as named-by-group-level vectors instead of data-frame columns
  ## -- ready to pass straight through as 'prior_list' to
  ## rLMMindepNormalGamma_reg_known_vcov()/_estimated_vcov() (which now
  ## accept this list argument-name-for-argument-name; note 'mu'/'Sigma' are
  ## deliberately NOT included here: 'mu' is never read by those routes --
  ## each group's mean is always W_j %*% gamma, recomputed every sweep -- and
  ## 'Sigma' is the Block~2-derived random-effects covariance, which those
  ## routes now always derive internally from 'pfamily_list' as
  ## solve(P) rather than accepting it here, so it can never drift out of
  ## sync with pfamily_list). All four vectors share identical names/order
  ## (= group_levels) by construction: they come from the same per-group
  ## loop above, never independently subset afterward.
  attr(out, "measurement_prior_group") <- list(
    shape_group      = stats::setNames(window_diagnostics$shape_ING,  window_diagnostics$group),
    rate_group       = stats::setNames(window_diagnostics$rate,       window_diagnostics$group),
    disp_lower_group = stats::setNames(window_diagnostics$disp_lower, window_diagnostics$group),
    disp_upper_group = stats::setNames(window_diagnostics$disp_upper, window_diagnostics$group)
  )

  class(out) <- c("dGamma_list", "list")
  out
}
