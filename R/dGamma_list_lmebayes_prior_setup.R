#' Build per-group dGamma priors from a Prior_Setup_lmebayes object
#'
#' Converts the per-group Block~1 measurement-dispersion calibration stored
#' in a \code{\link{Prior_Setup_lmebayes}} object into a named list of
#' \code{\link{dGamma}} \code{pfamily} objects, one per group level.
#'
#' Prior density (\code{shape_ING}, \code{rate}) comes from
#' \code{object$ing_prior_measurement_group}, calibrated in
#' \code{Prior_Setup_lmebayes()} via \code{\link{compute_gaussian_prior}} with
#' shared population \code{sd_tau}.  As of the Chapter A12 Part~VI switch,
#' this is the "full marginal" calibration (\eqn{\mu_j \to \bar\mu}, also
#' integrating out \code{fixef} against Block~2's own \code{dNormal} prior),
#' falling back to the Part~I calibration (\eqn{S_{\mathrm{marg}}}, A12 3.3.4)
#' when Part~VI's calibration fails; see
#' \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} Part VI and
#' \code{object$ing_prior_measurement_group_classical} (the Part~I
#' calibration, always computed, kept for comparison/fallback). Truncation
#' bounds use an approximate posterior at
#' \eqn{n_{\mathrm{combined},j} = n_{\mathrm{prior},j} + n_j}, mean-matched at
#' one of two per-group centers (\code{disp_center}):
#' \describe{
#'   \item{\code{"sigma2_hat"} (default)}{\code{disp_lower} is OLS-anchored at
#'     \code{sigma2_hat_j}; \code{disp_upper} widens by the BLUP/OLS residual
#'     RSS inflation ratio when \code{disp_upper_anchor = "blup"} (default).}
#'   \item{\code{"dispersion2"}}{Both bounds are mean-matched, symmetrically,
#'     at a dispersion estimate that integrates over the random effect's own
#'     posterior uncertainty (an \code{EnvelopeCentering()}-style trace
#'     correction; \code{disp_upper_anchor} is ignored).  Tends to produce
#'     narrower, better-centered upper tails than \code{"sigma2_hat"} for
#'     groups with large BLUP/OLS inflation.}
#' }
#' See \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} (Parts II and III) for
#' the full derivation of both constructions and their empirical validation.
#'
#' The returned list carries attribute \code{"window_diagnostics"}: a data
#' frame with cross-percentiles of the asymmetric truncation window and flag
#' \code{asymmetric_window} when
#' \eqn{R_{\mathrm{lo}} = \mathrm{lo\_pct\_BLUP}/\mathrm{lo\_pct\_OLS} <}
#' \code{asymmetric_R_lo} (default \code{0.25}) or
#' \eqn{R_{\mathrm{hi}} = \mathrm{hi\_pct\_OLS}/\mathrm{hi\_pct\_BLUP} >}
#' \code{asymmetric_R_hi} (default \code{4}).  A \code{\link[base]{warning}}
#' is emitted for flagged groups when \code{warn_asymmetric = TRUE}.
#'
#' @param object An object of class \code{"lmebayes_prior_setup"} as returned
#'   by \code{\link{Prior_Setup_lmebayes}} (Gaussian models only).
#' @param max_disp_perc Scalar in \eqn{(0.5, 1)}; defaults to
#'   \code{object$max_disp_perc}.
#' @param disp_center Character: \code{"sigma2_hat"} (default) mean-matches
#'   the truncation window at each group's calibrated point estimate
#'   (\code{disp_upper_anchor} controls upper-tail widening); \code{"dispersion2"}
#'   mean-matches symmetrically at an envelope-centering dispersion estimate
#'   that integrates over random-effect uncertainty instead (Part III of
#'   \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md}; ignores
#'   \code{disp_upper_anchor}).
#' @param disp_upper_anchor Character: \code{"blup"} (default) scales the
#'   upper-bound rate by \eqn{\mathrm{RSS}_{\mathrm{blup}}/\mathrm{RSS}_{\mathrm{ols}}}
#'   per group; \code{"symmetric"} uses the same rate as the lower bound.
#'   Ignored when \code{disp_center = "dispersion2"}.
#' @param n_rss_iter Positive integer; number of fixed-point iterations for
#'   the \code{disp_center = "dispersion2"} envelope-centering estimate
#'   (default \code{10L}). Ignored when \code{disp_center = "sigma2_hat"}.
#' @param warn_asymmetric If \code{TRUE} (default), warn when any group's
#'   truncation window is flagged as asymmetric.  Defaults to option
#'   \code{glmbayesCore.dgamma_window_warn} when \code{NULL}.
#' @param print_asymmetric If \code{TRUE}, print a table of flagged groups
#'   after the warning.  Defaults to option
#'   \code{glmbayesCore.dgamma_window_print} when \code{NULL}.
#' @param asymmetric_R_lo Flag when \code{R_lo = lo_pct_BLUP/lo_pct_OLS} is
#'   below this threshold (default \code{0.25}, or option
#'   \code{glmbayesCore.dgamma_window_R_lo}).
#' @param asymmetric_R_hi Flag when \code{R_hi = hi_pct_OLS/hi_pct_BLUP} exceeds
#'   this threshold (default \code{4}, or option
#'   \code{glmbayesCore.dgamma_window_R_hi}).
#' @param ... Currently ignored.
#'
#' @return A named list of \code{"pfamily"} objects keyed by group levels,
#'   suitable for \code{lmerb(..., dispersion_ranef = dGamma_list(ps))}, with
#'   attributes \code{"window_diagnostics"} (data frame, one row per group),
#'   \code{"dispersion_fit"} (\code{object$dispersion_fit}, the \code{glmmTMB}
#'   reference fit used for calibration -- \code{lmerb()}/\code{glmerb()}
#'   reuse it as their own \code{dispersion_fit} instead of re-fitting
#'   \code{glmmTMB} when \code{dispersion_ranef} carries this attribute), and
#'   \code{"calibration_source"} (\code{object$calibration_source}).
#'
#' @seealso \code{\link{Prior_Setup_lmebayes}}, \code{\link{dGamma_list}},
#'   \code{\link{dGamma}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("bayesrules", quietly = TRUE)) {
#'   data(big_word_club, package = "bayesrules")
#'   dat <- big_word_club
#'   dat$school_id <- factor(dat$school_id)
#'   dat <- subset(dat, !is.na(score_ppvt))
#'
#'   ps <- Prior_Setup_lmebayes(
#'     score_ppvt ~ private_school + (1 | school_id),
#'     data = dat,
#'     pwt_measurement = 0.01,
#'     dispformula = ~school_id
#'   )
#'   disp_pf <- dGamma_list(ps)
#'   print(disp_pf[[1L]])
#' }
#' }
#'
#' @export
#' @method dGamma_list lmebayes_prior_setup
dGamma_list.lmebayes_prior_setup <- function(
    object,
    max_disp_perc = NULL,
    disp_center = c("sigma2_hat", "dispersion2"),
    disp_upper_anchor = c("blup", "symmetric"),
    warn_asymmetric = NULL,
    print_asymmetric = NULL,
    asymmetric_R_lo = NULL,
    asymmetric_R_hi = NULL,
    n_rss_iter = 10L,
    ...
) {
  disp_center       <- match.arg(disp_center)
  disp_upper_anchor <- match.arg(disp_upper_anchor)
  if (!is.numeric(n_rss_iter) || length(n_rss_iter) != 1L ||
      is.na(n_rss_iter) || n_rss_iter < 1L) {
    stop("'n_rss_iter' must be a positive integer.", call. = FALSE)
  }
  n_rss_iter <- as.integer(n_rss_iter)

  if (!identical(object$family$family, "gaussian")) {
    stop(
      "dGamma_list() for lmebayes_prior_setup requires family = gaussian().",
      call. = FALSE
    )
  }

  ing_grp <- object$ing_prior_measurement_group
  if (is.null(ing_grp)) {
    grp_nm <- object$design$group_name
    stop(
      "object has no ing_prior_measurement_group; call Prior_Setup_lmebayes(",
      "..., dispformula = ~", if (!is.null(grp_nm)) grp_nm else "<group_name>",
      ") on a Gaussian model to calibrate per-group measurement-dispersion ",
      "priors (dispformula = ~1, the default, skips this calibration).",
      call. = FALSE
    )
  }

  if (is.null(max_disp_perc)) {
    max_disp_perc <- object$max_disp_perc
  }
  if (is.null(max_disp_perc)) {
    max_disp_perc <- 0.99
  }
  if (!is.numeric(max_disp_perc) || length(max_disp_perc) != 1L ||
      is.na(max_disp_perc) || max_disp_perc <= 0.5 || max_disp_perc >= 1) {
    stop("'max_disp_perc' must be a scalar in (0.5, 1).", call. = FALSE)
  }

  if (is.null(warn_asymmetric)) {
    warn_asymmetric <- getOption("glmbayesCore.dgamma_window_warn", TRUE)
  }
  if (is.null(print_asymmetric)) {
    print_asymmetric <- getOption("glmbayesCore.dgamma_window_print", TRUE)
  }
  if (is.null(asymmetric_R_lo)) {
    asymmetric_R_lo <- getOption("glmbayesCore.dgamma_window_R_lo", 0.25)
  }
  if (is.null(asymmetric_R_hi)) {
    asymmetric_R_hi <- getOption("glmbayesCore.dgamma_window_R_hi", 4)
  }
  if (!is.numeric(asymmetric_R_lo) || length(asymmetric_R_lo) != 1L ||
      !is.finite(asymmetric_R_lo) || asymmetric_R_lo <= 0 || asymmetric_R_lo >= 1) {
    stop("'asymmetric_R_lo' must be a single value in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(asymmetric_R_hi) || length(asymmetric_R_hi) != 1L ||
      !is.finite(asymmetric_R_hi) || asymmetric_R_hi <= 1) {
    stop("'asymmetric_R_hi' must be a single value > 1.", call. = FALSE)
  }

  group_levels <- names(ing_grp)
  if (is.null(group_levels)) {
    group_levels <- levels(object$design$groups)
    names(ing_grp) <- group_levels
  }

  blup_inflation <- if (disp_center == "sigma2_hat" && disp_upper_anchor == "blup") {
    .lmebayes_group_blup_rss_inflation(
      data          = object$data,
      block_formula = object$block_formula,
      fit_ref       = object$fit_ref,
      groups        = object$design$groups,
      group_levels  = group_levels,
      group_name    = object$design$group_name
    )
  } else {
    stats::setNames(rep(1, length(group_levels)), group_levels)
  }

  ## Part III (inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md): mean-match
  ## symmetrically at an envelope-centering dispersion estimate that
  ## integrates over b_j's own posterior uncertainty, instead of
  ## BLUP-inflating sigma2_hat. Computed once per group up front (like
  ## blup_inflation above), not inside the per-group lapply() below.
  dispersion2_vals <- if (disp_center == "dispersion2") {
    .lmebayes_group_dispersion2_envelope_centering(
      data              = object$data,
      block_formula     = object$block_formula,
      Sigma_ranef       = object$Sigma_ranef,
      groups            = object$design$groups,
      group_levels      = group_levels,
      intercept_source  = object$intercept_source,
      effects_source    = object$effects_source,
      n_rss_iter        = n_rss_iter
    )
  } else {
    NULL
  }

  diag_rows <- vector("list", length(group_levels))
  out <- stats::setNames(
    lapply(seq_along(group_levels), function(i) {
      lev <- group_levels[[i]]
      g <- ing_grp[[lev]]
      p_re <- g$p_re
      n_combined <- g$n_combined
      sigma2_hat <- g$sigma2_hat

      shape_w <- (n_combined + 1) / 2 + p_re / 2

      if (disp_center == "dispersion2") {
        rate_w <- dispersion2_vals[[lev]] * (n_combined + p_re - 1) / 2
        rate_u <- rate_w
      } else {
        rate_w <- sigma2_hat * (n_combined + p_re - 1) / 2
        rate_u <- if (disp_upper_anchor == "blup") {
          rate_w * blup_inflation[[lev]]
        } else {
          rate_w
        }
      }

      if (rate_u < rate_w - sqrt(.Machine$double.eps) * max(rate_w, 1)) {
        stop(
          "Group '", lev, "': upper-bound rate (", rate_u,
          ") is below lower-bound rate (", rate_w, ").",
          call. = FALSE
        )
      }

      blup_infl_val <- if (disp_center == "dispersion2") {
        NA_real_
      } else {
        unname(blup_inflation[[lev]])
      }
      dispersion2_val <- if (disp_center == "dispersion2") {
        unname(dispersion2_vals[[lev]])
      } else {
        NA_real_
      }

      xwin <- .lmebayes_dgamma_window_cross_percentiles(
        shape         = shape_w,
        rate_w        = rate_w,
        rate_u        = rate_u,
        max_disp_perc = max_disp_perc,
        blup_infl     = blup_infl_val,
        sigma2_hat    = unname(sigma2_hat)
      )

      diag_rows[[i]] <<- data.frame(
        group              = lev,
        n_j                = g$n_j,
        sigma2_hat         = unname(sigma2_hat),
        dispersion2        = dispersion2_val,
        blup_infl          = blup_infl_val,
        disp_lower         = xwin$disp_lower,
        disp_upper         = xwin$disp_upper,
        lo_pct_OLS         = xwin$lo_pct_OLS,
        lo_pct_BLUP        = xwin$lo_pct_BLUP,
        hi_pct_BLUP        = xwin$hi_pct_BLUP,
        hi_pct_OLS         = xwin$hi_pct_OLS,
        R_lo               = xwin$R_lo,
        R_hi               = xwin$R_hi,
        asymmetric_window  = .lmebayes_dgamma_window_asymmetric_flag(
          xwin$R_lo, xwin$R_hi,
          asymmetric_R_lo = asymmetric_R_lo,
          asymmetric_R_hi = asymmetric_R_hi
        ),
        stringsAsFactors   = FALSE
      )

      dGamma(
        shape          = g$shape_ING,
        rate           = g$rate,
        beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
        Inv_Dispersion = TRUE,
        max_disp_perc  = max_disp_perc,
        disp_lower     = xwin$disp_lower,
        disp_upper     = xwin$disp_upper
      )
    }),
    group_levels
  )

  window_diagnostics <- do.call(rbind, diag_rows)
  rownames(window_diagnostics) <- NULL
  attr(out, "window_diagnostics")  <- window_diagnostics
  ## Carried forward so lmerb()/glmerb() can reuse this glmmTMB reference fit
  ## as their diagnostic-only dispersion_fit instead of re-fitting it; see
  ## lmebayes:::.lmebayes_fit_glmmtmb_dispersion().
  attr(out, "dispersion_fit")      <- object$dispersion_fit
  attr(out, "calibration_source")  <- object$calibration_source

  if (isTRUE(warn_asymmetric)) {
    .lmebayes_warn_dgamma_window_asymmetry(
      window_diagnostics,
      asymmetric_R_lo = asymmetric_R_lo,
      asymmetric_R_hi   = asymmetric_R_hi,
      print_table       = isTRUE(print_asymmetric)
    )
  }

  out
}
