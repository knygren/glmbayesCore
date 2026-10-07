## Fit-object S3 methods for plot_var_convergence()/plot_mean_convergence(),
## so callers can plot straight from an rlmerb() or rLMM*_reg() fit
## instead of assembling 'hist'/'design'/'measurement_prior_list' by hand.
## All the class-specific logic lives in the normalizer helpers below;
## every fit class's method is the same thin wrapper. The rglmerb() /
## rGLMM_reg() methods stay in lmebayesCore.

#' Resolve the sweep-history object, \code{n_chains}, and (when eligible) the
#' exact reference ingredients from a fit object
#'
#' Shared by every \code{plot_var_convergence.*}/\code{plot_mean_convergence.*}
#' fit-object method. \code{"main"} reads
#' \code{fit$sweep_history}/\code{fit$n} (the \code{n} chains passed to the
#' sampler); \code{"pilot"} reads \code{fit$pilot$draws$sweep_history}/
#' \code{fit$pilot$n}.
#' @noRd
.lmebayes_convergence_inputs <- function(fit, stage = c("main", "pilot")) {
  stage <- match.arg(stage)
  hist <- if (identical(stage, "pilot")) {
    fit$pilot$draws$sweep_history
  } else {
    fit$sweep_history
  }
  if (is.null(hist)) {
    hint <- if (identical(stage, "main")) {
      paste0(
        "This is expected for rLMMNormal_reg_known_vcov(sim_method = \"DEFAULT\")",
        "/rlmerb(sim_method = \"DEFAULT\") exact-iid fits (no Gibbs sweeps, ",
        "nothing to plot convergence for) -- pass sim_method = \"TWO_BLOCK_GIBBS\" ",
        "for a sweep history."
      )
    } else {
      "No pilot stage ran for this fit (e.g. n_pilot = 0, or a route that never pilots)."
    }
    stop(
      "No '", stage, "' sweep_history on this fit (class ",
      class(fit)[1L], "). ", hint,
      call. = FALSE
    )
  }
  list(
    hist                   = hist,
    design                 = fit$design,
    measurement_prior_list = .lmebayes_convergence_measurement_prior_list(fit),
    n_chains               = if (identical(stage, "pilot")) fit$pilot$n else fit$n
  )
}

#' Is an *exact* Claim~1/Claim~3 reference available for this fit?
#'
#' Exact references (\code{\link{lmerb_posterior_mean}}/
#' \code{\link{lmerb_posterior_covariance}}) require a Gaussian response with
#' *both* the observation dispersion and the Block~2 hyperparameter vcov
#' fixed (not sampled): \code{rLMMNormal_reg_known_vcov()} (and
#' \code{rlmerb()}/\code{rglmerb()} routed there) qualify;
#' \code{rLMMNormal_reg_estimated_vcov()} (vcov estimated),
#' \code{rLMMindepNormalGamma_reg_*()} (dispersion always estimated per
#' group), and \code{rGLMM_reg_*()} (non-Gaussian in every real
#' \code{glmerb()} call path -- \code{gaussian()} routes to
#' \code{rLMMNormal_reg_*()}/\code{rLMMindepNormalGamma_reg_*()} instead) do
#' not. \code{rlmerb()}/\code{rglmerb()} lose the informative
#' \code{result_class} (their final \code{class()} is always just
#' \code{"rlmerb"}/\code{"rglmerb"}), so they are resolved from
#' \code{fit$Prior$dispersion_mode} instead of \code{inherits()}.
#' @noRd
.lmebayes_convergence_exact_ref_ok <- function(fit) {
  if (inherits(fit, "rGLMM_reg")) {
    return(FALSE)
  }
  fam <- fit$family
  if (is.null(fam) || !identical(fam$family, "gaussian")) {
    return(FALSE)
  }
  disp_ok <- if (!is.null(fit$Prior$dispersion_mode)) {
    fit$Prior$dispersion_mode %in% c("fixed", "fixed_vector")
  } else {
    !inherits(fit, "rLMMindepNormalGamma_reg")
  }
  isTRUE(disp_ok) && !isTRUE(fit$any_non_normal)
}

#' Build \code{measurement_prior_list} from a fit's own \code{dispprior_list}/
#' \code{prior_list}/\code{pfamily_list}/\code{design}/\code{family}, or
#' \code{NULL} when no exact reference is available (see
#' \code{\link{.lmebayes_convergence_exact_ref_ok}})
#'
#' LMM routes store the Block~1 measurement prior as
#' \code{fit$dispprior_list} (with fallback to legacy \code{fit$prior_list}).
#' The exact-iid engine (\code{sim_method = "DEFAULT"}) does not store its
#' internally-derived Block~1 precision back onto that object's \code{P}, so
#' it is recomputed here via \code{.rLMM_P_from_pfamily_list()} whenever
#' missing (deterministic given \code{pfamily_list}, so this always agrees
#' with whatever precision the sampler actually used).
#' @noRd
.lmebayes_convergence_measurement_prior_list <- function(fit) {
  if (!.lmebayes_convergence_exact_ref_ok(fit)) {
    return(NULL)
  }
  design   <- fit$design
  re_names <- design$groupef.names
  prior_list_block1 <- fit$dispprior_list
  if (is.null(prior_list_block1)) {
    prior_list_block1 <- fit$prior_list
  }
  if (is.null(prior_list_block1$P) && is.null(prior_list_block1$Sigma)) {
    prior_list_block1$P <- .rLMM_P_from_pfamily_list(fit$pfamily_list, re_names)
  }
  .two_block_measurement_prior_list(
    prior_list_block1 = prior_list_block1,
    pfamily_list       = fit$pfamily_list,
    re_names           = re_names,
    x_hyper            = design$W,
    family             = fit$family
  )
}

#' Shared body for every \code{plot_var_convergence.*}/\code{plot_mean_convergence.*}
#' fit-object method: resolve inputs from \code{fit} then delegate to
#' \code{default_fn} (one of \code{\link{plot_var_convergence.default}}/
#' \code{\link{plot_mean_convergence.default}})
#'
#' \code{n_chains} is only auto-filled from the fit (\code{inputs$n_chains})
#' when the caller omits it entirely (\code{missing(n_chains)}); explicitly
#' passing \code{n_chains = NULL} still means "no band" (forwarded as-is to
#' \code{default_fn}), matching \code{\link{plot_var_convergence.default}}'s
#' own convention. \code{...} (e.g. \code{split}, \code{max_whitened}) is
#' forwarded to \code{default_fn} as-is, so fit-object methods need not
#' repeat every \code{.default}-only formal in their own signature.
#' @noRd
.lmebayes_convergence_from_fit <- function(
    fit, default_fn, coef_focus, whitened, engine, n_chains, has_n_chains,
    conf_level, stage, stage_label, ...
) {
  inputs <- .lmebayes_convergence_inputs(fit, stage = stage)
  if (!has_n_chains) {
    n_chains <- inputs$n_chains
  }
  if (is.null(stage_label)) {
    stage_label <- inputs$hist$stage
  }
  default_fn(
    inputs$hist,
    coef_focus             = coef_focus,
    design                 = inputs$design,
    measurement_prior_list = inputs$measurement_prior_list,
    whitened               = whitened,
    engine                 = engine,
    n_chains               = n_chains,
    conf_level             = conf_level,
    stage_label            = stage_label,
    ...
  )
}

#' @param stage \code{"main"} (default) or \code{"pilot"} -- which sweep
#'   history to plot (\code{fit$sweep_history} or
#'   \code{fit$pilot$draws$sweep_history}). Ignored by
#'   \code{\link{plot_var_convergence.default}}, which always takes
#'   \code{hist} literally.
#' @param stage_label Defaults to the resolved \code{hist$stage}; pass
#'   explicitly to override.
#' @rdname plot_var_convergence
#' @method plot_var_convergence rLMMNormal_reg
#' @export
plot_var_convergence.rLMMNormal_reg <- function(
    hist,
    coef_focus = NULL,
    component = c("fixef", "precision"),
    whitened = FALSE,
    engine = c("base", "ggplot"),
    n_chains,
    conf_level = 0.95,
    stage = c("main", "pilot"),
    stage_label = NULL,
    ...
) {
  component <- match.arg(component)
  engine <- match.arg(engine)
  stage  <- match.arg(stage)
  .lmebayes_convergence_from_fit(
    hist, plot_var_convergence.default,
    coef_focus = coef_focus, component = component, whitened = whitened, engine = engine,
    n_chains = if (missing(n_chains)) NULL else n_chains,
    has_n_chains = !missing(n_chains),
    conf_level = conf_level, stage = stage, stage_label = stage_label,
    ...
  )
}

#' @rdname plot_var_convergence
#' @method plot_var_convergence rLMMindepNormalGamma_reg
#' @export
plot_var_convergence.rLMMindepNormalGamma_reg <- plot_var_convergence.rLMMNormal_reg

#' @rdname plot_var_convergence
#' @method plot_var_convergence rlmerb
#' @export
plot_var_convergence.rlmerb <- plot_var_convergence.rLMMNormal_reg

#' @param stage \code{"main"} (default) or \code{"pilot"}; see
#'   \code{\link{plot_var_convergence.rLMMNormal_reg}}.
#' @param stage_label Defaults to the resolved \code{hist$stage}; pass
#'   explicitly to override.
#' @rdname plot_mean_convergence
#' @method plot_mean_convergence rLMMNormal_reg
#' @export
plot_mean_convergence.rLMMNormal_reg <- function(
    hist,
    coef_focus = NULL,
    component = c("fixef", "precision"),
    whitened = FALSE,
    engine = c("base", "ggplot"),
    n_chains,
    conf_level = 0.95,
    stage = c("main", "pilot"),
    stage_label = NULL,
    ...
) {
  component <- match.arg(component)
  engine <- match.arg(engine)
  stage  <- match.arg(stage)
  .lmebayes_convergence_from_fit(
    hist, plot_mean_convergence.default,
    coef_focus = coef_focus, component = component, whitened = whitened, engine = engine,
    n_chains = if (missing(n_chains)) NULL else n_chains,
    has_n_chains = !missing(n_chains),
    conf_level = conf_level, stage = stage, stage_label = stage_label,
    ...
  )
}

#' @rdname plot_mean_convergence
#' @method plot_mean_convergence rLMMindepNormalGamma_reg
#' @export
plot_mean_convergence.rLMMindepNormalGamma_reg <- plot_mean_convergence.rLMMNormal_reg

#' @rdname plot_mean_convergence
#' @method plot_mean_convergence rlmerb
#' @export
plot_mean_convergence.rlmerb <- plot_mean_convergence.rLMMNormal_reg
