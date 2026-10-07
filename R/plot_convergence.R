#' Combined variance-convergence plot for Block~2 fixed effects (Claim 3)
#'
#' Tracks convergence of the Block~2 fixed-effects (hyperparameter) posterior
#' variance across inner Gibbs sweeps as a single combined chart (all traces
#' overlaid on one axis, following \code{lmebayes::Ex_16}'s style), per the
#' two-block Gibbs ergodicity result (Claim~3 of the package's ergodicity
#' reference, \code{inst/BLOCK_GIBBS_ERGODICITY.md}): after \code{l} sweeps,
#' the cross-chain covariance \eqn{\Sigma^{(l)}_{11}} satisfies
#' \deqn{
#'   \Sigma_{11}^{-1/2} \Sigma^{(l)}_{11} \Sigma_{11}^{-1/2} = I - A^{2l}
#' }
#' for a convergent (symmetric, eigenvalues in \eqn{[0, 1)}) matrix \eqn{A}.
#' \code{Var(l) / Var_ref} ratios are therefore bounded above by 1 and
#' increase toward it as \code{l} grows.
#'
#' @param hist Object of class \code{"two_block_sweep_history"} (see
#'   \code{\link{print.two_block_sweep_history}}). For \code{whitened = TRUE},
#'   \code{hist} must also carry \code{cov_by_sweep}/\code{coef_index}, which
#'   are only populated by sweeps-outer/chains-inner engines -- i.e. every
#'   engine driven by \code{\link{rGLMM_sweep}}: the \verb{_run_with_pilot()}
#'   family (\code{rLMMNormal_reg_estimated_vcov()},
#'   \code{rLMMindepNormalGamma_reg_known_vcov()},
#'   \code{rLMMindepNormalGamma_reg_estimated_vcov()}) and, as of the
#'   \code{rGLMM_sweep()}-based rewrite of its engine,
#'   \code{rLMMNormal_reg_known_vcov(sim_method = "TWO_BLOCK_GIBBS")} as well.
#' @param coef_focus Optional list of length-2 character vectors
#'   \code{c(re_component, covariate)} to restrict the plot to; default
#'   (\code{NULL}) uses every coefficient tracked in \code{hist}.
#' @param component \code{"fixef"} (default) plots the Block~2
#'   fixed-effects/hyperparameter series described above. \code{"precision"}
#'   instead plots the estimated random-effects \emph{precision}
#'   (\eqn{1/\tau^2_k}, i.e. the diagonal of \eqn{\Sigma_{ranef}^{-1}}) for
#'   every ING (\code{dIndependent_Normal_Gamma}) Block~2 component --
#'   requires \code{hist$disp_table} (populated whenever at least one
#'   Block~2 component has an estimated, not fixed, variance -- see
#'   \code{fit$any_non_normal}); errors if that is empty (e.g. a
#'   known-vcov fit) or if \code{whitened = TRUE} (only the cross-chain
#'   \emph{variance} of each component's precision is captured so far, not
#'   its covariance with other components/fixed effects, so whitening is
#'   not yet possible -- see \code{inst/BLOCK_GIBBS_ERGODICITY.md} for the
#'   deferred extension). Precision is plotted (not the variance
#'   \eqn{\tau^2_k} itself) because \eqn{E[1/\tau^2_k]} stays well-defined
#'   under a weak/vague prior where \eqn{E[\tau^2_k]} need not be. This is
#'   a distinct quantity from the observation-level dispersion estimated by
#'   \code{rLMMindepNormalGamma_reg_*()} (despite both being called
#'   "dispersion" internally) -- \code{component = "precision"} is always
#'   about the RE variance-covariance, never the observation noise
#'   variance. \code{design}/\code{measurement_prior_list} are ignored for
#'   \code{component = "precision"}: no exact reference exists yet for an
#'   estimated \eqn{\Sigma_{ranef}}, so the denominator is always the
#'   empirical last-sweep value.
#' @param design,measurement_prior_list When both are supplied, the *exact*
#'   reference covariance is resolved automatically via
#'   \code{\link{lmerb_posterior_covariance}} and used as the denominator --
#'   this is not an opt-in toggle, it is used whenever the ingredients for it
#'   exist. Only valid when dispersion and the RE variance-covariance are
#'   both fixed (not sampled). When omitted (the only option for estimated
#'   dispersion or estimated-vcov models), the denominator falls back to the
#'   empirical last-sweep cross-chain variance/covariance ("Var_final").
#' @param whitened When \code{FALSE} (default), plots \code{Var(l) / Var_ref}
#'   per named coefficient in \code{coef_focus}. When \code{TRUE}, whitens
#'   each sweep's full cross-chain covariance by the reference covariance
#'   (\eqn{\Sigma_{ref}^{-1/2} \mathrm{Cov}(l) \Sigma_{ref}^{-1/2}}) and plots
#'   its eigenvalues -- basis-invariant, so no eigenvectors of the
#'   theoretical \eqn{A} are ever needed -- labeled \code{var1, var2, ...}
#'   instead of named coefficients (per Claim~3, these eigen-directions mix
#'   the original coefficients together, so named labels no longer apply).
#' @param engine \code{"base"} (default) or \code{"ggplot"} (requires
#'   \pkg{ggplot2}); both draw one combined chart with all traces overlaid.
#' @param n_chains Optional integer, the number of independent chains each
#'   plotted \code{Var(l)} (and \code{Var_ref}, when it is the empirical
#'   \code{Var_final} rather than an exact reference) was computed across --
#'   i.e. the same \code{n_chains}/\code{n_pilot}/\code{n} passed to the
#'   sweep-generating call (see Details). When supplied, draws a single pair
#'   of horizontal dotted reference lines: a naive \code{conf_level} band
#'   around 1 for \code{Ratio(l)}, under the null that the true ratio at that
#'   sweep already equals 1 (unrealistic at small \code{l}, but useful as a
#'   ruler for how much of the departure from 1 is just Monte Carlo noise).
#'   The band's distribution depends on whether \code{Var_ref} is exact or
#'   empirical (see \code{design}/\code{measurement_prior_list} above):
#'   \itemize{
#'     \item \code{ref_source = "exact"}: only \code{Var(l)} (the numerator)
#'       has sampling error, so
#'       \eqn{\mathrm{Ratio}(l) \sim \chi^2_{n_{chains}-1} / (n_{chains}-1)}.
#'     \item \code{ref_source = "empirical (Var_final)"}: \code{Var_ref}
#'       itself is a sample variance/covariance across the same
#'       \code{n_chains} chains (just at the final sweep), so both
#'       numerator and denominator carry sampling error and the classic
#'       two-sample variance-ratio result applies instead:
#'       \eqn{\mathrm{Ratio}(l) \sim F(n_{chains}-1,\, n_{chains}-1)}. This
#'       is *wider* than the chi-squared band, and reduces to it exactly as
#'       \code{Var_final}'s own degrees of freedom go to infinity
#'       (\code{stats::qf(p, df1, Inf) == stats::qchisq(p, df1) / df1}), so
#'       both cases are computed via a single \code{stats::qf()} call with
#'       \code{df2 = Inf} for the exact case.
#'   }
#'   Not series-specific (a single band, not one per coefficient/eigenvalue)
#'   since it depends only on \code{n_chains}/\code{ref_source}. For
#'   \code{whitened = TRUE} this is an additional approximation regardless of
#'   \code{ref_source}: whitened eigenvalues are eigenvalues of a
#'   Wishart-distributed matrix, not marginally chi-squared/F, and tend to
#'   spread out more than this band implies (eigenvalue repulsion), so treat
#'   it as an illustrative ruler, not an exact whitened confidence band. Also
#'   note that \code{Var(l)} and \code{Var_final} come from the *same*
#'   \code{n_chains} chains at two different sweep indices, not independent
#'   samples the way a textbook two-sample F-test assumes -- another reason
#'   this is a naive/illustrative band rather than an exact one. Omit (the
#'   default, \code{NULL}) to skip the band entirely.
#' @param conf_level Confidence level for the \code{n_chains} band (default
#'   \code{0.95}); ignored when \code{n_chains} is \code{NULL}.
#' @param stage_label Character label for the title; defaults to \code{hist$stage}.
#' @param split \code{"auto"} (default) or \code{"none"}. \code{"auto"} draws
#'   one \emph{separate} chart -- its own plot page/window, never stacked via
#'   \code{\link[graphics]{layout}}/\code{par(mfrow = ...)} on top of another
#'   -- per group: for \code{whitened = FALSE}, one chart for every
#'   coefficient whose \code{re_component} is \code{"(Intercept)"}
#'   ("Intercept predictors") and one for every other \code{re_component}
#'   combined ("Slope predictors"), whichever of the two is non-empty; for
#'   \code{whitened = TRUE}, consecutive batches of at most
#'   \code{max_whitened} eigenvalues each. \code{"none"} restores the
#'   single-combined-chart behaviour (every tracked coefficient/eigenvalue
#'   overlaid on one panel), same as before \code{split} existed.
#' @param max_whitened Maximum number of eigenvalue series per chart when
#'   \code{whitened = TRUE} and \code{split = "auto"} (default \code{4});
#'   ignored otherwise.
#' @return \code{hist} invisibly.
#' @details
#' ## What is plotted
#'
#' The x-axis is the inner Gibbs sweep index \eqn{l}, not the stored-draw
#' index: these charts diagnose convergence \emph{within} the burn-in that
#' precedes each stored draw.  Each point is a \strong{cross-chain} statistic
#' at sweep \eqn{l} --- the variance across the \code{n_chains} independent
#' replicate chains, all of which are at sweep \eqn{l} simultaneously --- not
#' a running statistic along a single chain.  This is why the sweep-outer
#' loop order matters: it is what makes "the sampler's distribution after
#' \eqn{l} sweeps" an observable quantity.
#'
#' The y-axis is the ratio \code{Var(l) / Var_ref}.  Under Claim 3 that ratio
#' is bounded above by 1 and rises toward it, so a converged run shows every
#' trace flattening just below 1.  Traces that are still climbing at the last
#' sweep indicate the sweep budget is too small; traces that overshoot 1
#' substantially are Monte Carlo noise, and \code{n_chains} draws the band
#' that says how much overshoot is expected.
#'
#' \code{Var_ref} is the exact posterior variance when \code{design} and
#' \code{measurement_prior_list} are both supplied (only possible when
#' \eqn{\sigma^2} and the RE variance components are all fixed), and
#' otherwise the empirical last-sweep variance.  The empirical fallback
#' cannot detect a run that has not converged \emph{at all}, since the
#' denominator is then drawn from the same unconverged distribution; prefer
#' the exact reference whenever the model allows it.
#'
#' ## Reading it alongside the mean chart
#'
#' Location and spread converge at different speeds, and the two plotting
#' functions separate them: \code{\link{plot_mean_convergence}} tracks the
#' mean, \code{plot_var_convergence} the variance.  The variance converges
#' at twice the rate, so the mean chart settles later and is usually the one
#' that decides whether the sweep budget is adequate.  A run whose variance
#' ratios sit at 1 while its standardized mean deviations still drift has
#' converged in spread but not in location.
#'
#' Because every chain starts at the posterior mean, the mean chart is flat
#' from sweep 1 on well-behaved Gaussian fits, leaving the variance chart to
#' carry the information.  The lmebayesCore vignette Chapter-B02 derives
#' the two rates.
#'
#' ## When there is nothing to plot
#'
#' \code{sweep_history} is recorded only when all chains advance together,
#' one sweep at a time, so that the cross-chain spread at sweep \eqn{l} is
#' meaningful.  Two situations produce no history:
#'
#' \itemize{
#'   \item The exact route (fixed \eqn{\sigma^2}, all components
#'     \code{dNormal}, \code{sim_method = "DEFAULT"}) runs no sweeps at all,
#'     because the draws are already exact.  There is nothing to diagnose.
#'     Passing \code{sim_method = "TWO_BLOCK_GIBBS"} forces the iterative
#'     path on the same target if you want a history anyway.
#'   \item Some Gaussian paths run sweeps but do not export per-sweep
#'     cross-chain statistics; \code{fit$sweep_history} is then \code{NULL}.
#' }
#'
#' \code{whitened = TRUE} additionally requires \code{cov_by_sweep} and
#' \code{coef_index}, which only the sweep-outer engines above populate, and
#' \code{component = "precision"} additionally requires \code{disp_table},
#' which requires at least one ING Block 2 component.
#'
#' ## Supplying n_chains
#'
#' \code{n_chains} is not stored on \code{hist} itself, so it must be
#' supplied from the fitted object that produced it: for the main-stage
#' history (\code{fit$sweep_history}) use \code{nrow(fit$popef[[k]])} for
#' any RE component \code{k} (all have the same number of rows, the \code{n}
#' passed to the \verb{_run_with_pilot()} call); for the pilot-stage history
#' (\code{fit$pilot$draws$sweep_history}) use \code{fit$pilot$n}.
#' @param ... Passed to methods (e.g. \code{coef_focus}, \code{design},
#'   \code{whitened}, \code{n_chains}, ...; see below). Ignored by
#'   \code{\link{plot_var_convergence.default}} beyond the documented
#'   formals.
#' @seealso \code{\link{plot_mean_convergence}},
#'   \code{\link{lmerb_posterior_covariance}}
#' @export
plot_var_convergence <- function(hist, ...) {
  UseMethod("plot_var_convergence")
}

#' @rdname plot_var_convergence
#' @method plot_var_convergence default
#' @export
plot_var_convergence.default <- function(
    hist,
    coef_focus = NULL,
    component = c("fixef", "precision"),
    design = NULL,
    measurement_prior_list = NULL,
    whitened = FALSE,
    engine = c("base", "ggplot"),
    n_chains = NULL,
    conf_level = 0.95,
    stage_label = hist$stage,
    split = c("auto", "none"),
    max_whitened = 4L,
    ...
) {
  if (!inherits(hist, "two_block_sweep_history")) {
    stop(
      "'hist' must be a two_block_sweep_history object ",
      "(e.g. fit$sweep_history$main).",
      call. = FALSE
    )
  }
  component <- match.arg(component)
  engine <- match.arg(engine)
  split  <- match.arg(split)
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) {
    stage_label <- if (!is.null(hist$stage)) hist$stage else "stage"
  }
  if (identical(engine, "ggplot") && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("'engine = \"ggplot\"' requires the ggplot2 package.", call. = FALSE)
  }
  if (identical(component, "precision")) {
    if (is.null(hist$disp_table) || !nrow(hist$disp_table)) {
      stop(
        "'hist' has no RE-precision rows ('disp_table') -- component = ",
        "\"precision\" requires a sweep history from an estimated-vcov ",
        "engine with at least one ING Block~2 component (e.g. ",
        "rLMMNormal_reg_estimated_vcov(), rLMMindepNormalGamma_reg_estimated_vcov()).",
        call. = FALSE
      )
    }
    if (isTRUE(whitened)) {
      stop(
        "'whitened = TRUE' is not yet supported for component = \"precision\" ",
        "-- only the cross-chain variance of each RE component's precision ",
        "is captured so far, not its cross-chain covariance with other ",
        "components/fixed effects.",
        call. = FALSE
      )
    }
  }

  sh_tab <- if (identical(component, "precision")) hist$disp_table else hist$table
  sh_tab$re_component <- as.character(sh_tab$re_component)
  sh_tab$covariate <- as.character(sh_tab$covariate)
  sh_sweeps <- sh_tab[sh_tab$sweep > 0L, , drop = FALSE]
  if (!nrow(sh_sweeps)) {
    warning("No sweep rows in sweep history for stage ", stage_label, call. = FALSE)
    return(invisible(hist))
  }

  ## No exact reference exists yet for the RE-precision series (no analogue
  ## of lmerb_posterior_covariance() for an estimated Sigma_ranef), so
  ## design/measurement_prior_list are only consulted for component = "fixef".
  exact_ref <- NULL
  ref_source <- "empirical (Var_final)"
  if (identical(component, "fixef") && !is.null(design) && !is.null(measurement_prior_list)) {
    exact_ref <- lmerb_posterior_covariance(design, measurement_prior_list)
    ref_source <- "exact"
  }
  ## Exact reference => only Var(l) has sampling error (chi-sq/df); empirical
  ## Var_final reference => both numerator and denominator are sample
  ## variances across the same n_chains chains (two-sample F). The former is
  ## the df2 = Inf limit of the latter, so one qf() call covers both.
  band <- .convergence_var_band(n_chains, conf_level, identical(ref_source, "exact"))

  metric_label <- if (identical(component, "precision")) "RE precision (1/tau^2)" else "cross-chain variance"
  if (isTRUE(whitened)) {
    var_ratio_full <- .convergence_var_whitened_series(hist, coef_focus, exact_ref)
    ylab <- "Whitened variance ratio (eigenvalue)"
    sub <- sprintf(
      "(whitened cross-chain covariance eigenvalues; reference = %s; Claim 3: I - A^(2l))",
      ref_source
    )
    groups <- .convergence_whitened_groups(var_ratio_full, split, max_whitened)
  } else {
    ylab <- if (identical(ref_source, "exact")) {
      "Var / Var_exact ratio"
    } else {
      "Var / Var_final ratio"
    }
    sub <- sprintf("(%s from sweep history; reference = %s)", metric_label, ref_source)
    groups <- .convergence_named_groups(
      sh_sweeps, coef_focus, exact_ref, .convergence_var_named_series, split
    )
  }
  if (!is.null(band)) {
    ref_note <- if (identical(ref_source, "exact")) {
      "chi-sq"
    } else {
      "F, Var_ref also estimated"
    }
    sub <- paste0(
      sub,
      sprintf(
        "\ndotted: naive %.0f%% band around 1 (%s; n_chains = %d%s)",
        conf_level * 100, ref_note, n_chains,
        if (isTRUE(whitened)) "; approx. for eigenvalues" else ""
      )
    )
  }

  .convergence_render_groups(
    groups, stage_label, ylab, sub, engine, band, ref_line = 1, floor_at_ref = TRUE
  )

  invisible(hist)
}

#' Combined mean-convergence plot for Block~2 fixed effects (Claim 1)
#'
#' Tracks convergence of the Block~2 fixed-effects (hyperparameter) posterior
#' \emph{mean} across inner Gibbs sweeps as a single combined chart (all
#' traces overlaid on one axis), the mean-side counterpart of
#' \code{\link{plot_var_convergence}}'s variance-side Claim~3 chart. Per the
#' two-block Gibbs ergodicity result's Claim~1, the cross-chain mean bias
#' \eqn{E[\gamma^{(l)}] - \gamma^*} decays like \eqn{A^l} (linear in \code{l},
#' vs. Claim~3's \eqn{A^{2l}} for variance). Each coefficient's per-sweep
#' deviation from its reference mean is standardized by that same sweep's
#' cross-chain SD, \eqn{Z(l) = (\mathrm{mean}(l) - \mathrm{mean}_{ref}) /
#' \mathrm{sd}(l)}, so coefficients on very different scales can be overlaid
#' on one chart and \code{Z(l)} is expected to shrink toward 0 as \code{l}
#' grows.
#'
#' @param hist,coef_focus,component,engine,n_chains,conf_level,split,max_whitened
#'   Same meaning as in \code{\link{plot_var_convergence}}.
#' @param design,measurement_prior_list When both are supplied, the exact
#'   reference \emph{mean} is resolved automatically via
#'   \code{\link{lmerb_posterior_mean}} and, for \code{whitened = TRUE}, the
#'   exact reference \emph{covariance} (via
#'   \code{\link{lmerb_posterior_covariance}}) needed to whiten the mean
#'   deviation vector; see \code{\link{plot_var_convergence}} for when these
#'   are available.
#' @param stage_label Defaults to the resolved \code{hist$stage}; pass
#'   explicitly to override.
#' @param whitened When \code{FALSE} (default), plots \code{Z(l)} per named
#'   coefficient in \code{coef_focus}. When \code{TRUE}, whitens each sweep's
#'   cross-chain mean deviation vector by the reference covariance
#'   (\eqn{\Sigma_{ref}^{-1/2} (\mathrm{mean}(l) - \mathrm{mean}_{ref})}) and
#'   plots its components -- basis-invariant -- labeled \code{z1, z2, ...}
#'   instead of named coefficients.
#' @param ... Passed to methods; see \code{\link{plot_mean_convergence.default}}.
#' @return \code{hist} invisibly.
#' @details
#' ## What is plotted
#'
#' The x-axis is the inner Gibbs sweep index \eqn{l}, and each point is a
#' \strong{cross-chain} statistic computed over the \code{n_chains}
#' independent replicate chains, all of which sit at sweep \eqn{l}
#' simultaneously.  Nothing here is a running statistic along a single
#' chain.
#'
#' The plotted quantity is the standardized deviation
#' \deqn{Z(l) = \frac{\mathrm{mean}(l) - \mathrm{mean}_{ref}}{\mathrm{sd}(l)},}
#' which divides out the scale of each coefficient so that an intercept
#' measured in hundreds and a slope measured in tenths can share one axis.
#' A converged run shows every trace scattered around 0 with no trend; a
#' systematic drift toward or away from 0 means the sweep budget is too
#' small.  Supplying \code{n_chains} adds a naive band that separates
#' Monte Carlo noise from genuine bias.
#'
#' \code{mean_ref} is the exact posterior mean when \code{design} and
#' \code{measurement_prior_list} are supplied (only possible when
#' \eqn{\sigma^2} and the RE variance components are all fixed), and
#' otherwise the empirical last-sweep mean.  As with
#' \code{\link{plot_var_convergence}}, the empirical fallback cannot detect a
#' run that never converged, because the reference is drawn from the same
#' unconverged distribution.
#'
#' ## Reading it alongside the variance chart
#'
#' This chart is the location-side counterpart of
#' \code{\link{plot_var_convergence}}, which tracks spread.  The variance
#' converges at twice the rate of the mean, so this is the slower --- and
#' therefore the more informative --- of the two when deciding whether
#' \code{m_convergence} is large enough.  A fit whose variance ratios have
#' settled at 1 while these traces still drift has converged in spread but
#' not in location.
#'
#' Because every replicate chain starts at the posterior mean, the bias
#' begins at zero on Gaussian fits with fixed variance components: expect a
#' flat chart from sweep 1 there.  Read a non-flat one as evidence that
#' something upstream --- an estimated variance component, or a start at the
#' mode rather than the mean --- has moved the starting point off the target
#' mean.  The lmebayesCore vignette Chapter-B02 derives
#' the two rates.
#'
#' ## When there is nothing to plot
#'
#' As for \code{\link{plot_var_convergence}}: the exact route runs no sweeps
#' to record, some Gaussian paths do not export per-sweep cross-chain
#' statistics, and \code{whitened = TRUE} needs slots that only the
#' sweep-recording engines populate.  See that topic, which also explains
#' how to obtain \code{n_chains} from a fitted object.
#' @seealso \code{\link{plot_var_convergence}},
#'   \code{\link{lmerb_posterior_mean}}
#' @export
plot_mean_convergence <- function(hist, ...) {
  UseMethod("plot_mean_convergence")
}

#' @rdname plot_mean_convergence
#' @method plot_mean_convergence default
#' @export
plot_mean_convergence.default <- function(
    hist,
    coef_focus = NULL,
    component = c("fixef", "precision"),
    design = NULL,
    measurement_prior_list = NULL,
    whitened = FALSE,
    engine = c("base", "ggplot"),
    n_chains = NULL,
    conf_level = 0.95,
    stage_label = hist$stage,
    split = c("auto", "none"),
    max_whitened = 4L,
    ...
) {
  if (!inherits(hist, "two_block_sweep_history")) {
    stop(
      "'hist' must be a two_block_sweep_history object ",
      "(e.g. fit$sweep_history$main).",
      call. = FALSE
    )
  }
  component <- match.arg(component)
  engine <- match.arg(engine)
  split  <- match.arg(split)
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) {
    stage_label <- if (!is.null(hist$stage)) hist$stage else "stage"
  }
  if (identical(engine, "ggplot") && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("'engine = \"ggplot\"' requires the ggplot2 package.", call. = FALSE)
  }
  if (identical(component, "precision")) {
    if (is.null(hist$disp_table) || !nrow(hist$disp_table)) {
      stop(
        "'hist' has no RE-precision rows ('disp_table') -- component = ",
        "\"precision\" requires a sweep history from an estimated-vcov ",
        "engine with at least one ING Block~2 component (e.g. ",
        "rLMMNormal_reg_estimated_vcov(), rLMMindepNormalGamma_reg_estimated_vcov()).",
        call. = FALSE
      )
    }
    if (isTRUE(whitened)) {
      stop(
        "'whitened = TRUE' is not yet supported for component = \"precision\" ",
        "-- only the cross-chain variance of each RE component's precision ",
        "is captured so far, not its cross-chain covariance with other ",
        "components/fixed effects.",
        call. = FALSE
      )
    }
  }

  sh_tab <- if (identical(component, "precision")) hist$disp_table else hist$table
  sh_tab$re_component <- as.character(sh_tab$re_component)
  sh_tab$covariate <- as.character(sh_tab$covariate)
  sh_sweeps <- sh_tab[sh_tab$sweep > 0L, , drop = FALSE]
  if (!nrow(sh_sweeps)) {
    warning("No sweep rows in sweep history for stage ", stage_label, call. = FALSE)
    return(invisible(hist))
  }

  ## No exact reference exists yet for the RE-precision series (no analogue
  ## of lmerb_posterior_mean() for an estimated Sigma_ranef), so
  ## design/measurement_prior_list are only consulted for component = "fixef".
  exact_ref <- NULL
  ref_source <- "empirical (mean_final)"
  if (identical(component, "fixef") && !is.null(design) && !is.null(measurement_prior_list)) {
    pm <- lmerb_posterior_mean(design, measurement_prior_list)
    exact_ref <- list(
      mean = .convergence_flatten_fixef(pm$fixef),
      cov  = lmerb_posterior_covariance(design, measurement_prior_list)
    )
    ref_source <- "exact"
  }
  ## Exact reference => mean_ref is a fixed/known constant, so the one-sample
  ## t-statistic has n_chains - 1 df (Var(l) is estimated) but a fixed mean;
  ## empirical mean_final reference => mean_ref is itself a sample mean from
  ## the same n_chains chains, so both the sweep-l and reference sampling
  ## error are naively combined (same "df2 = Inf limit" trick as the
  ## variance band via qt(df = Inf) == qnorm()).
  band <- .convergence_mean_band(n_chains, conf_level, identical(ref_source, "exact"))

  metric_label <- if (identical(component, "precision")) {
    "RE precision (1/tau^2) cross-chain mean deviation"
  } else {
    "cross-chain mean deviation"
  }
  if (isTRUE(whitened)) {
    z_series_full <- .convergence_mean_whitened_series(hist, coef_focus, exact_ref)
    ylab <- "Whitened mean deviation (z-score)"
    sub <- sprintf(
      "(whitened cross-chain mean deviation; reference = %s; Claim 1: A^l bias decay)",
      ref_source
    )
    groups <- .convergence_whitened_groups(z_series_full, split, max_whitened)
  } else {
    exact_ref_mean <- if (!is.null(exact_ref)) exact_ref$mean else NULL
    ylab <- if (identical(ref_source, "exact")) {
      "(mean - mean_exact) / SD"
    } else {
      "(mean - mean_final) / SD"
    }
    sub <- sprintf(
      "(%s, standardized by cross-chain SD; reference = %s)",
      metric_label, ref_source
    )
    groups <- .convergence_named_groups(
      sh_sweeps, coef_focus, exact_ref_mean, .convergence_mean_named_series, split
    )
  }
  if (!is.null(band)) {
    ref_note <- if (identical(ref_source, "exact")) {
      "normal"
    } else {
      "t, mean_ref also estimated"
    }
    sub <- paste0(
      sub,
      sprintf(
        "\ndotted: naive %.0f%% band around 0 (%s; n_chains = %d%s)",
        conf_level * 100, ref_note, n_chains,
        if (isTRUE(whitened)) "; approx. for whitened components" else ""
      )
    )
  }

  .convergence_render_groups(
    groups, stage_label, ylab, sub, engine, band, ref_line = 0, floor_at_ref = FALSE
  )

  invisible(hist)
}

#' Naive band around 1 for a variance-ratio estimated from \code{n_chains}
#' draws, against either an exact or an empirical reference variance
#'
#' Under the (naive) null that the true ratio at a sweep already equals 1:
#' \itemize{
#'   \item \code{exact_ref = TRUE} (\code{Var_ref} is a fixed/known number,
#'     e.g. \code{lmerb_posterior_covariance()}'s \code{M^{-1}}): only the
#'     numerator \code{Var(l)} has sampling error, so
#'     \eqn{(n_{chains}-1) \cdot \mathrm{Ratio}(l) \sim \chi^2_{n_{chains}-1}},
#'     i.e. \code{Ratio(l)} itself is distributed as
#'     \eqn{\chi^2_{n_{chains}-1} / (n_{chains}-1)}.
#'   \item \code{exact_ref = FALSE} (\code{Var_ref} is the empirical
#'     \code{Var_final}, itself a sample variance across the same
#'     \code{n_chains} chains, just at the last sweep): both numerator and
#'     denominator are sample variances of the same true (co)variance, so
#'     the classic two-sample variance-ratio result applies instead:
#'     \eqn{\mathrm{Ratio}(l) \sim F(n_{chains}-1,\, n_{chains}-1)}.
#' }
#' These are computed via a single \code{stats::qf()} call: the exact case
#' is simply \code{df2 = Inf}, since
#' \code{stats::qf(p, df1, Inf) == stats::qchisq(p, df1) / df1} exactly (a
#' known limiting identity -- a fixed/known reference variance behaves like
#' a denominator sample variance with infinite degrees of freedom). Both are
#' exact for the named-coefficient (non-whitened) series, and only an
#' approximation for whitened eigenvalues (see
#' \code{\link{plot_var_convergence}}'s \code{n_chains} docs).
#' @noRd
.convergence_var_band <- function(n_chains, conf_level, exact_ref) {
  if (is.null(n_chains)) {
    return(NULL)
  }
  n_chains <- as.integer(n_chains)[1L]
  if (is.na(n_chains) || n_chains <= 1L) {
    stop("'n_chains' must be a single integer > 1.", call. = FALSE)
  }
  conf_level <- as.numeric(conf_level)[1L]
  if (is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single value strictly between 0 and 1.", call. = FALSE)
  }
  df1 <- n_chains - 1L
  df2 <- if (isTRUE(exact_ref)) Inf else df1
  alpha <- 1 - conf_level
  stats::qf(c(alpha / 2, 1 - alpha / 2), df1 = df1, df2 = df2)
}

#' Naive band around 0 for a standardized mean deviation \code{Z(l)}
#' estimated from \code{n_chains} draws, against either an exact or an
#' empirical reference mean
#'
#' Companion to \code{\link{.convergence_var_band}} for
#' \code{\link{plot_mean_convergence}}'s \code{Z(l) = (mean(l) - mean_ref) /
#' sd(l)} statistic: under the (naive) null that the true mean at sweep
#' \code{l} already equals \code{mean_ref}, \code{sqrt(n_chains) * Z(l)} is a
#' one-sample t-statistic on \code{n_chains - 1} df when \code{mean_ref} is a
#' fixed/known constant (\code{exact_ref = TRUE}), or an approximate
#' (naive) t-statistic on the same df when \code{mean_ref} is itself a
#' sample mean from the same chains (\code{exact_ref = FALSE}, e.g.
#' \code{mean_final}) -- mirroring \code{.convergence_var_band()}'s
#' \code{df2 = Inf} trick, \code{stats::qt(p, df = Inf) == stats::qnorm(p)}
#' exactly, so a single \code{stats::qt()} call covers both cases.
#' @noRd
.convergence_mean_band <- function(n_chains, conf_level, exact_ref) {
  if (is.null(n_chains)) {
    return(NULL)
  }
  n_chains <- as.integer(n_chains)[1L]
  if (is.na(n_chains) || n_chains <= 1L) {
    stop("'n_chains' must be a single integer > 1.", call. = FALSE)
  }
  conf_level <- as.numeric(conf_level)[1L]
  if (is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single value strictly between 0 and 1.", call. = FALSE)
  }
  df <- if (isTRUE(exact_ref)) Inf else n_chains - 1L
  alpha <- 1 - conf_level
  q <- stats::qt(1 - alpha / 2, df = df)
  c(-q, q) / sqrt(n_chains)
}

#' Flatten \code{lmerb_posterior_mean()$fixef} into a \code{"re_component |
#' covariate"}-named vector matching \code{\link{lmerb_posterior_covariance}}'s
#' \code{lbl} convention
#' @noRd
.convergence_flatten_fixef <- function(fixef_list) {
  unlist(lapply(names(fixef_list), function(k) {
    v <- fixef_list[[k]]
    stats::setNames(as.numeric(v), paste(k, names(v), sep = " | "))
  }))
}

#' Resolve \code{coef_focus} into a (re_component, covariate) key data frame
#' @noRd
.convergence_resolve_keys <- function(coef_focus, available) {
  if (is.null(coef_focus)) {
    keys <- unique(available[, c("re_component", "covariate")])
    rownames(keys) <- NULL
    return(keys)
  }
  if (!is.list(coef_focus) || !length(coef_focus)) {
    stop("'coef_focus' must be a non-empty list of c(re_component, covariate) pairs.",
         call. = FALSE)
  }
  keys <- do.call(rbind, lapply(coef_focus, function(cc) {
    if (length(cc) < 2L) {
      stop("Each element of 'coef_focus' must be c(re_component, covariate).",
           call. = FALSE)
    }
    data.frame(
      re_component = as.character(cc[1L]),
      covariate    = as.character(cc[2L]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(keys) <- NULL
  keys
}

#' \code{"re_component | covariate"} labels for a key data frame
#' @noRd
.convergence_labels <- function(keys) {
  paste(keys$re_component, keys$covariate, sep = " | ")
}

#' Check that every label in \code{lbl} is a row/col name of \code{mat}
#' @noRd
.convergence_check_labels <- function(mat, lbl, what) {
  missing <- setdiff(lbl, rownames(mat))
  if (length(missing)) {
    stop(
      what, " has no entry for: ", paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Named-coefficient \code{Var(l) / Var_ref} series (non-whitened mode)
#' @noRd
.convergence_var_named_series <- function(sh_sweeps, coef_focus, exact_ref) {
  keys <- .convergence_resolve_keys(coef_focus, sh_sweeps)
  lbl  <- .convergence_labels(keys)
  if (!is.null(exact_ref)) {
    .convergence_check_labels(exact_ref, lbl, "Exact reference covariance")
  }

  var_list <- stats::setNames(
    lapply(seq_len(nrow(keys)), function(i) {
      sub <- sh_sweeps[
        sh_sweeps$re_component == keys$re_component[i] &
          sh_sweeps$covariate == keys$covariate[i],
        ,
        drop = FALSE
      ]
      if (!nrow(sub)) {
        stop("No sweep rows for ", lbl[i], call. = FALSE)
      }
      sub <- sub[order(sub$sweep), , drop = FALSE]
      stats::setNames(sub$sd^2, sub$sweep)
    }),
    lbl
  )

  var_ref <- stats::setNames(vapply(seq_along(var_list), function(i) {
    if (!is.null(exact_ref)) {
      exact_ref[lbl[[i]], lbl[[i]]]
    } else {
      v <- var_list[[i]]
      m_final <- max(as.integer(names(v)))
      as.numeric(v[as.character(m_final)])
    }
  }, numeric(1L)), lbl)

  stats::setNames(
    lapply(seq_along(var_list), function(i) var_list[[i]] / var_ref[[i]]),
    lbl
  )
}

#' Named-coefficient \code{Z(l) = (mean(l) - mean_ref) / sd(l)} series
#' (non-whitened mode)
#' @noRd
.convergence_mean_named_series <- function(sh_sweeps, coef_focus, exact_ref) {
  keys <- .convergence_resolve_keys(coef_focus, sh_sweeps)
  lbl  <- .convergence_labels(keys)
  if (!is.null(exact_ref)) {
    missing <- setdiff(lbl, names(exact_ref))
    if (length(missing)) {
      stop(
        "Exact reference mean has no entry for: ", paste(missing, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  stats::setNames(
    lapply(seq_len(nrow(keys)), function(i) {
      sub <- sh_sweeps[
        sh_sweeps$re_component == keys$re_component[i] &
          sh_sweeps$covariate == keys$covariate[i],
        ,
        drop = FALSE
      ]
      if (!nrow(sub)) {
        stop("No sweep rows for ", lbl[i], call. = FALSE)
      }
      sub <- sub[order(sub$sweep), , drop = FALSE]
      ref_i <- if (!is.null(exact_ref)) {
        exact_ref[[lbl[i]]]
      } else {
        sub$mean[[which.max(sub$sweep)]]
      }
      sd_i <- sub$sd
      sd_i[!is.finite(sd_i) | sd_i <= 0] <- NA_real_
      stats::setNames((sub$mean - ref_i) / sd_i, sub$sweep)
    }),
    lbl
  )
}

#' Whitened \code{Var(l) / Var_ref} eigenvalue series (whitened mode)
#' @noRd
.convergence_var_whitened_series <- function(hist, coef_focus, exact_ref) {
  if (is.null(hist$cov_by_sweep) || !length(hist$cov_by_sweep)) {
    stop(
      "'hist' has no per-sweep covariance ('cov_by_sweep') -- whitened = TRUE ",
      "requires a sweep history built by a sweeps-outer/chains-inner engine ",
      "(e.g. rLMMNormal_reg_estimated_vcov(), rLMMindepNormalGamma_reg_known_vcov(), ",
      "rLMMindepNormalGamma_reg_estimated_vcov()).",
      call. = FALSE
    )
  }
  coef_index <- hist$coef_index
  keys <- .convergence_resolve_keys(coef_focus, coef_index)
  lbl  <- .convergence_labels(keys)

  n_sweeps <- length(hist$cov_by_sweep)
  .convergence_check_labels(hist$cov_by_sweep[[n_sweeps]], lbl, "Sweep-history covariance")

  Sigma_ref <- if (!is.null(exact_ref)) {
    .convergence_check_labels(exact_ref, lbl, "Exact reference covariance")
    exact_ref[lbl, lbl, drop = FALSE]
  } else {
    hist$cov_by_sweep[[n_sweeps]][lbl, lbl, drop = FALSE]
  }

  eg <- eigen(Sigma_ref, symmetric = TRUE)
  d_inv_sqrt <- ifelse(eg$values > .Machine$double.eps, 1 / sqrt(eg$values), 0)
  Sigma_ref_inv_sqrt <- eg$vectors %*% (d_inv_sqrt * t(eg$vectors))

  P <- length(lbl)
  eig_by_sweep <- matrix(NA_real_, nrow = n_sweeps, ncol = P)
  for (m in seq_len(n_sweeps)) {
    Cov_l <- hist$cov_by_sweep[[m]][lbl, lbl, drop = FALSE]
    W_l <- Sigma_ref_inv_sqrt %*% Cov_l %*% Sigma_ref_inv_sqrt
    W_l <- (W_l + t(W_l)) / 2
    eig_by_sweep[m, ] <- eigen(W_l, symmetric = TRUE, only.values = TRUE)$values
  }

  var_names <- paste0("var", seq_len(P))
  stats::setNames(
    lapply(seq_len(P), function(j) {
      stats::setNames(eig_by_sweep[, j], seq_len(n_sweeps))
    }),
    var_names
  )
}

#' Whitened mean-deviation series (whitened mode)
#'
#' Whitens each sweep's cross-chain mean-deviation \emph{vector} (not a
#' per-sweep matrix, unlike \code{\link{.convergence_var_whitened_series}})
#' by the reference covariance's inverse square root, computed once:
#' \eqn{z^{(l)} = \Sigma_{ref}^{-1/2} (\mathrm{mean}(l) - \mathrm{mean}_{ref})}.
#' @noRd
.convergence_mean_whitened_series <- function(hist, coef_focus, exact_ref) {
  if (is.null(hist$cov_by_sweep) || !length(hist$cov_by_sweep)) {
    stop(
      "'hist' has no per-sweep covariance ('cov_by_sweep') -- whitened = TRUE ",
      "requires a sweep history built by a sweeps-outer/chains-inner engine ",
      "(e.g. rLMMNormal_reg_estimated_vcov(), rLMMindepNormalGamma_reg_known_vcov(), ",
      "rLMMindepNormalGamma_reg_estimated_vcov()).",
      call. = FALSE
    )
  }
  coef_index <- hist$coef_index
  keys <- .convergence_resolve_keys(coef_focus, coef_index)
  lbl  <- .convergence_labels(keys)

  n_sweeps <- length(hist$cov_by_sweep)
  .convergence_check_labels(hist$cov_by_sweep[[n_sweeps]], lbl, "Sweep-history covariance")

  sh_tab <- hist$table
  sh_tab$re_component <- as.character(sh_tab$re_component)
  sh_tab$covariate    <- as.character(sh_tab$covariate)
  sh_sweeps <- sh_tab[sh_tab$sweep > 0L, , drop = FALSE]

  Sigma_ref <- if (!is.null(exact_ref)) {
    .convergence_check_labels(exact_ref$cov, lbl, "Exact reference covariance")
    exact_ref$cov[lbl, lbl, drop = FALSE]
  } else {
    hist$cov_by_sweep[[n_sweeps]][lbl, lbl, drop = FALSE]
  }

  mean_ref <- if (!is.null(exact_ref)) {
    missing <- setdiff(lbl, names(exact_ref$mean))
    if (length(missing)) {
      stop(
        "Exact reference mean has no entry for: ", paste(missing, collapse = ", "), ".",
        call. = FALSE
      )
    }
    as.numeric(exact_ref$mean[lbl])
  } else {
    vapply(seq_len(nrow(keys)), function(i) {
      sub <- sh_sweeps[
        sh_sweeps$re_component == keys$re_component[i] &
          sh_sweeps$covariate == keys$covariate[i],
        ,
        drop = FALSE
      ]
      sub$mean[[which.max(sub$sweep)]]
    }, numeric(1L))
  }

  eg <- eigen(Sigma_ref, symmetric = TRUE)
  d_inv_sqrt <- ifelse(eg$values > .Machine$double.eps, 1 / sqrt(eg$values), 0)
  Sigma_ref_inv_sqrt <- eg$vectors %*% (d_inv_sqrt * t(eg$vectors))

  P <- length(lbl)
  mean_by_sweep <- do.call(rbind, lapply(seq_len(n_sweeps), function(m) {
    vapply(seq_len(nrow(keys)), function(i) {
      v <- sh_sweeps$mean[
        sh_sweeps$re_component == keys$re_component[i] &
          sh_sweeps$covariate == keys$covariate[i] &
          sh_sweeps$sweep == m
      ]
      if (length(v) == 1L) v else NA_real_
    }, numeric(1L))
  }))

  z_by_sweep <- matrix(NA_real_, nrow = n_sweeps, ncol = P)
  for (m in seq_len(n_sweeps)) {
    dev <- as.numeric(mean_by_sweep[m, ] - mean_ref)
    z_by_sweep[m, ] <- as.numeric(Sigma_ref_inv_sqrt %*% dev)
  }

  z_names <- paste0("z", seq_len(P))
  stats::setNames(
    lapply(seq_len(P), function(j) {
      stats::setNames(z_by_sweep[, j], seq_len(n_sweeps))
    }),
    z_names
  )
}

#' Split resolved \code{(re_component, covariate)} keys into per-plot groups
#'
#' Two-block Gibbs models with a random intercept plus one or more random
#' slopes naturally split into (a) the intercept's own hyper-predictors and
#' (b) every slope's hyper-predictors combined -- keeping each chart's
#' traces/legend legible instead of cramming every coefficient (e.g. 7+
#' series) onto one panel. Returns a named list of \code{coef_focus}-style
#' lists (one \code{c(re_component, covariate)} pair per row of that
#' group's keys); the name is the group's plot-title suffix. When every key
#' already shares the same intercept/non-intercept status (so there is only
#' one non-empty group), returns a single unnamed (\code{""}) group,
#' signalling "no real split" to the caller.
#' @noRd
.convergence_split_coef_focus <- function(keys) {
  to_coef_focus <- function(k) {
    lapply(seq_len(nrow(k)), function(i) c(k$re_component[i], k$covariate[i]))
  }
  is_icpt <- keys$re_component == "(Intercept)"
  groups <- list()
  if (any(is_icpt)) {
    groups[["Intercept predictors"]] <- to_coef_focus(keys[is_icpt, , drop = FALSE])
  }
  if (any(!is_icpt)) {
    groups[["Slope predictors"]] <- to_coef_focus(keys[!is_icpt, , drop = FALSE])
  }
  if (length(groups) <= 1L) {
    return(stats::setNames(list(to_coef_focus(keys)), ""))
  }
  groups
}

#' Chunk a whitened series list (\code{var1, var2, ...}/\code{z1, z2, ...})
#' into consecutive groups of at most \code{max_per_plot} components
#'
#' Unlike \code{\link{.convergence_split_coef_focus}}, splitting happens
#' \emph{after} the full eigen-decomposition (over every tracked coefficient
#' jointly, as Claim~3/Claim~1 require) rather than by re-running it on a
#' subset -- this only chunks the already-computed eigenvalue/whitened-
#' component series for legibility, it never changes their values. Returns
#' \code{list(series)} unnamed (\code{""}, "no real split") when
#' \code{length(series) <= max_per_plot}.
#' @noRd
.convergence_split_whitened <- function(series, max_per_plot) {
  p <- length(series)
  if (!p) {
    return(list())
  }
  if (p <= max_per_plot) {
    return(stats::setNames(list(series), ""))
  }
  idx <- split(seq_len(p), ceiling(seq_len(p) / max_per_plot))
  groups <- lapply(idx, function(ix) series[ix])
  labels <- vapply(idx, function(ix) {
    if (length(ix) > 1L) {
      paste0(names(series)[ix[1L]], "-", names(series)[ix[length(ix)]])
    } else {
      names(series)[ix]
    }
  }, character(1L))
  stats::setNames(groups, labels)
}

#' Build the named-coefficient series for every split group (shared by
#' \code{\link{plot_var_convergence.default}}/
#' \code{\link{plot_mean_convergence.default}})
#'
#' \code{split = "none"} skips grouping entirely and calls \code{series_fn}
#' once on the original (unsplit) \code{coef_focus}, i.e. the pre-\code{split}
#' single-combined-chart behaviour.
#' @noRd
.convergence_named_groups <- function(sh_sweeps, coef_focus, exact_ref, series_fn, split) {
  if (identical(split, "none")) {
    return(stats::setNames(list(series_fn(sh_sweeps, coef_focus, exact_ref)), ""))
  }
  keys <- .convergence_resolve_keys(coef_focus, sh_sweeps)
  coef_focus_groups <- .convergence_split_coef_focus(keys)
  stats::setNames(
    lapply(coef_focus_groups, function(cf) series_fn(sh_sweeps, cf, exact_ref)),
    names(coef_focus_groups)
  )
}

#' Split an already-computed whitened series list into per-plot groups
#' (shared by \code{\link{plot_var_convergence.default}}/
#' \code{\link{plot_mean_convergence.default}})
#' @noRd
.convergence_whitened_groups <- function(series_full, split, max_whitened) {
  if (identical(split, "none")) {
    return(stats::setNames(list(series_full), ""))
  }
  .convergence_split_whitened(series_full, max_whitened)
}

#' Render every split group as its own separate combined chart
#'
#' Each group gets its own \code{\link{.convergence_render}} call in
#' sequence -- never wrapped in a shared \code{\link[graphics]{layout}}/
#' \code{par(mfrow = ...)} across groups -- so on an interactive/on-screen
#' device each group's chart lands on its own plot page (never stacked
#' underneath another group's chart on the same page), and on a multi-page
#' file device (e.g. \code{\link[grDevices]{pdf}}) each group becomes its
#' own page.
#' @noRd
.convergence_render_groups <- function(groups, stage_label, ylab, sub, engine, band,
                                        ref_line, floor_at_ref) {
  for (g in names(groups)) {
    stage_label_g <- if (nzchar(g)) paste0(stage_label, " -- ", g) else stage_label
    cat(sprintf("\n=== %s sweep history (%s; %s) ===\n\n", stage_label_g, ylab, engine))
    .convergence_render(
      groups[[g]], stage_label_g, ylab, sub, engine, band, ref_line, floor_at_ref
    )
  }
  invisible(NULL)
}

#' Render one combined convergence chart (base or ggplot)
#' @noRd
.convergence_render <- function(series, stage_label, ylab, sub, engine, band = NULL,
                                 ref_line = 1, floor_at_ref = TRUE) {
  if (!length(series)) {
    warning("No convergence series to plot for stage ", stage_label, call. = FALSE)
    return(invisible(NULL))
  }
  if (identical(engine, "base")) {
    .convergence_plot_base(series, stage_label, ylab, sub, band, ref_line, floor_at_ref)
  } else {
    .convergence_plot_ggplot(series, stage_label, ylab, sub, band, ref_line)
  }
}

#' Recycled palette of maximally-distinguishable point shapes
#'
#' \code{graphics}' \code{pch} and \code{ggplot2}'s \code{shape} aesthetic
#' both accept the same integer codes, so this one vector drives point
#' shapes in both \code{\link{.convergence_plot_base}} and
#' \code{\link{.convergence_plot_ggplot}}. Colour alone (via
#' \code{grDevices::hcl.colors()}) can be hard to tell apart across many
#' similarly-hued series -- shape gives a second, colour-blind-friendly
#' channel to distinguish traces/legend entries. Recycles (with a
#' \code{warning()}) if there are more series than shapes.
#' @noRd
.convergence_shapes <- function(n) {
  ## filled circle, triangle, square, diamond, plus, cross, asterisk, open
  ## circle, open triangle, open square, open diamond, inverted triangle.
  shapes <- c(16L, 17L, 15L, 18L, 3L, 4L, 8L, 1L, 2L, 0L, 5L, 6L)
  rep(shapes, length.out = n)
}

#' Estimate how many legend columns fit the current device width
#'
#' \code{graphics::legend(..., ncol = k)} does not auto-wrap based on
#' available width -- if the labels don't actually fit \code{k} across, the
#' extra columns silently spill off the left/right of the device/panel
#' instead of wrapping to more rows (this is exactly what happened with a
#' hardcoded \code{ncol = 6}: long \code{"re_component | covariate"} labels,
#' e.g. \code{"distracted_a1 | free_reduced_lunch"}, don't fit six-across on
#' an ordinary device, and the overflow columns ended up positioned outside
#' the plot area). Estimates the widest label's rendered width (in inches,
#' via \code{graphics::strwidth()}, so it reflects the actual font/cex/
#' device), adds a fixed allowance for the point symbol and inter-column
#' gap, and divides that into the current device width to get a column
#' count that actually fits.
#' @noRd
.convergence_legend_ncol <- function(labels, cex = 0.8, max_ncol = 6L) {
  label_w_in <- max(graphics::strwidth(labels, units = "inches", cex = cex))
  ## Allowance for the pch point plus legend()'s own internal symbol/text
  ## gap (a generous fixed estimate; legend() itself does not expose this).
  col_w_in   <- label_w_in + 0.35
  avail_w_in <- graphics::par("din")[1L] * 0.94
  ncol_fit   <- max(1L, floor(avail_w_in / col_w_in))
  max(1L, min(max_ncol, ncol_fit, length(labels)))
}

#' Base-R combined convergence chart (all series on one panel + legend below)
#'
#' The legend lives in its own \code{\link[graphics]{layout}} panel below the
#' main plot, not hand-placed via \code{par("usr")} data-range arithmetic --
#' the previous approach positioned the legend at a fixed *fraction of the
#' plot's own y data range* below \code{usr[3]}, which (for any data range
#' larger than a couple of units) lands well outside the actual figure/device
#' area and is silently clipped, i.e. never drawn at all. A dedicated layout
#' panel is robust to the data range, the device size, and the number of
#' series (wraps to multiple legend rows via \code{ncol}).
#' @param band Optional length-2 numeric vector \code{c(lower, upper)} from
#'   \code{\link{.convergence_var_band}}/\code{\link{.convergence_mean_band}},
#'   drawn as a pair of horizontal dotted reference lines (not a per-series/
#'   legend entry).
#' @param ref_line Horizontal dashed reference line (\code{1} for the
#'   variance-ratio chart, \code{0} for the mean-deviation chart).
#' @param floor_at_ref When \code{TRUE} (variance ratios, which cannot be
#'   negative), the y-axis floor is fixed at \code{0}; when \code{FALSE}
#'   (mean deviations, which can be positive or negative), the y-axis is a
#'   padded, symmetric-ish range around the data/band/\code{ref_line}.
#' @noRd
.convergence_plot_base <- function(series, stage_label, ylab, sub, band = NULL,
                                    ref_line = 1, floor_at_ref = TRUE) {
  sweeps <- sort(unique(as.integer(unlist(lapply(series, names)))))
  if (!length(sweeps)) {
    warning("No sweep rows to plot for stage ", stage_label, call. = FALSE)
    return(invisible(NULL))
  }
  if (isTRUE(floor_at_ref)) {
    y_top <- max(ref_line, unlist(series), band, na.rm = TRUE)
    if (!is.finite(y_top) || y_top <= 0) {
      y_top <- ref_line
    }
    ylim <- c(0, y_top * 1.05)
  } else {
    y_rng <- range(c(ref_line, unlist(series), band), na.rm = TRUE, finite = TRUE)
    pad <- diff(y_rng) * 0.08
    if (!is.finite(pad) || pad <= 0) {
      pad <- 0.1
    }
    ylim <- y_rng + c(-pad, pad)
  }

  n_series <- length(series)
  cols   <- grDevices::hcl.colors(n_series, palette = "Dark 3")
  shapes <- .convergence_shapes(n_series)
  y_first <- series[[1L]][match(sweeps, names(series[[1L]]))]

  has_legend <- n_series > 1L
  legend_ncol <- if (has_legend) {
    .convergence_legend_ncol(names(series), cex = 0.8)
  } else {
    1L
  }
  n_legend_rows <- if (has_legend) ceiling(n_series / legend_ncol) else 0L

  ## 'sub' may carry multiple '\n'-separated lines (e.g. the naive-band
  ## note appended below the base caption) -- mtext() does not honour
  ## embedded newlines the way text() does, so split and place each line at
  ## its own 'line=' offset, sizing the bottom margin to fit them all.
  sub_lines <- if (nzchar(sub)) strsplit(sub, "\n", fixed = TRUE)[[1L]] else character(0L)
  bottom_mar <- 3.2 + 1.3 * max(1L, length(sub_lines))

  ## Only 'mar' (and the layout panel structure, reset separately via
  ## layout(1L)) are ever changed below -- restore just that, not a full
  ## par(no.readonly = TRUE) snapshot. Restoring derived/device-dependent
  ## entries from such a snapshot (e.g. 'pin') alongside a layout() panel
  ## change can raise "invalid value specified for graphical parameter
  ## 'pin'" if the device's plotting-region geometry at restore time no
  ## longer matches what was captured (mirrors the plain par(mfrow = ...)/
  ## local par(op) restore pattern.
  old_mar <- graphics::par("mar")
  on.exit({
    graphics::layout(1L)
    graphics::par(mar = old_mar)
  }, add = TRUE)

  if (has_legend) {
    ## Reserve a short bottom panel sized to the number of legend rows
    ## needed (each row ~0.09 of the plot panel's height), so the legend
    ## always has room regardless of how many series are being plotted.
    graphics::layout(matrix(1:2, nrow = 2L), heights = c(1, 0.12 + 0.09 * n_legend_rows))
    graphics::par(mar = c(bottom_mar, 4, 3, 1))
  } else {
    graphics::par(mar = c(bottom_mar, 4, 3, 1))
  }

  graphics::plot(
    sweeps, y_first,
    type = "n",
    xlab = "Inner sweep",
    ylab = ylab,
    ylim = ylim,
    main = paste0(stage_label, ": ", ylab, " vs sweep")
  )
  graphics::grid()
  graphics::abline(h = ref_line, lty = 2, col = "gray40")
  if (!is.null(band)) {
    graphics::abline(h = band, lty = 3, col = "gray50")
  }
  for (i in seq_len(n_series)) {
    v <- series[[i]]
    x <- as.integer(names(v))
    graphics::lines(x, v, type = "b", pch = shapes[i], col = cols[i])
  }
  for (i in seq_along(sub_lines)) {
    graphics::mtext(sub_lines[i], side = 1, line = 2.2 + 1.3 * i, cex = 0.8)
  }

  if (has_legend) {
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::legend(
      "center", legend = names(series), col = cols, pch = shapes,
      ncol = legend_ncol, bty = "n", cex = 0.8
    )
  }
  invisible(series)
}

#' ggplot2 combined convergence chart (all series on one panel)
#'
#' Shape (via \code{\link{.convergence_shapes}}, the same palette
#' \code{\link{.convergence_plot_base}} uses) is mapped alongside colour,
#' so series stay distinguishable even where the \code{hcl.colors()} hues
#' are close. \code{coef} is kept as a factor in \code{names(series)}'s
#' own order (not ggplot's default alphabetical order) so its levels line up
#' 1:1 with \code{shapes}.
#' @param band Optional length-2 numeric vector \code{c(lower, upper)} from
#'   \code{\link{.convergence_var_band}}/\code{\link{.convergence_mean_band}},
#'   drawn as a pair of horizontal dotted reference lines (not a per-series/
#'   legend entry).
#' @param ref_line Horizontal dashed reference line (\code{1} for the
#'   variance-ratio chart, \code{0} for the mean-deviation chart).
#' @noRd
.convergence_plot_ggplot <- function(series, stage_label, ylab, sub, band = NULL, ref_line = 1) {
  coef_levels <- names(series)
  shapes <- .convergence_shapes(length(coef_levels))

  df <- do.call(rbind, lapply(coef_levels, function(nm) {
    v <- series[[nm]]
    data.frame(
      sweep = as.integer(names(v)), value = as.numeric(v), coef = nm,
      stringsAsFactors = FALSE
    )
  }))
  df$coef <- factor(df$coef, levels = coef_levels)

  ## Columns referenced via !!as.name(...) (injected symbols), not bare names,
  ## so aes() sees ordinary symbols resolved against 'df' at plot-build time --
  ## no "no visible binding for global variable" NOTE, and no need for a
  ## bare (unqualified) '.data' pronoun, which only resolves correctly when
  ## written literally inside aes() (see .convergence_plot_ggplot() docs).
  aes_sweep <- as.name("sweep")
  aes_value <- as.name("value")
  aes_coef  <- as.name("coef")
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = !!aes_sweep, y = !!aes_value,
      group = !!aes_coef, colour = !!aes_coef, shape = !!aes_coef
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = ref_line, linetype = "dashed", colour = "gray40") +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::labs(
      x = "Inner sweep", y = ylab,
      title = paste0(stage_label, ": ", ylab, " vs sweep"),
      subtitle = sub, colour = NULL, shape = NULL
    )
  if (!is.null(band)) {
    p <- p + ggplot2::geom_hline(yintercept = band, linetype = "dotted", colour = "gray50")
  }
  print(p)
  invisible(series)
}
