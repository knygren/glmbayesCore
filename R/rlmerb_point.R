## Point-estimate (simulate = FALSE) front-door helpers for rlmerb / rglmerb.

#' @noRd
.rlmerb_prepare_prior <- function(
    design,
    pfamily_list,
    dispprior_list = NULL,
    family,
    fn_name,
    offset = NULL,
    weights = 1,
    offset_missing = FALSE,
    weights_missing = FALSE
) {
  if (!inherits(design, "model_setup")) {
    stop("'design' must be a model_setup object.", call. = FALSE)
  }
  if (is.null(pfamily_list)) {
    stop(
      "'pfamily_list' is required for ", fn_name, "() (named Block~2 pfamily ",
      "list, same as rLMMNormal_reg() / rGLMM_reg()).",
      call. = FALSE
    )
  }

  group.dispersion <- .lmebayes_dispprior_list_as_group_dispersion(dispprior_list)
  prior <- .priors_from_pfamily_list(
    pfamily_list     = pfamily_list,
    group.dispersion = group.dispersion,
    design           = design,
    family           = family,
    fn_name          = fn_name
  )
  disp_info <- list(
    mode                  = prior$dispersion_mode,
    dispersion_fix        = prior$group.dispersion,
    dispersion_prior_list = prior$dispersion_prior_list,
    dispersion_pfamily    = prior$dispersion_pfamily,
    window_diagnostics    = prior$window_diagnostics
  )
  ow <- .lmebayes_resolve_offset_weights(
    offset, weights, design,
    offset_missing = offset_missing,
    weights_missing = weights_missing
  )
  block1_prior <- .lmebayes_block1_prior_list(
    prior,
    group.dispersion = if (identical(family$family, "gaussian")) {
      disp_info$dispersion_fix
    } else {
      NULL
    }
  )
  list(
    prior        = prior,
    disp_info    = disp_info,
    block1_prior = block1_prior,
    ow           = ow,
    re_names     = design$groupef.names
  )
}

#' @noRd
.rlmerb_thin_Prior <- function(prior, disp_info, block1_prior) {
  list(
    block1_prior          = block1_prior,
    pfamily_list          = prior$pfamily_list,
    group.dispersion      = disp_info$dispersion_fix,
    dispersion_mode       = disp_info$mode,
    dispersion_pfamily    = disp_info$dispersion_pfamily,
    dispersion_prior_list = disp_info$dispersion_prior_list
  )
}

#' ICM / exact posterior location at fixed variance-component plug-ins
#' @noRd
.rlmerb_icm_at_fixed_vc <- function(
    design,
    prior,
    family,
    tol = 1e-10,
    maxit = 200L
) {
  measurement_prior_list <- list(
    group.Sigma      = prior$group.Sigma,
    pop.prior_list   = prior$pop.prior_list,
    group.dispersion = prior$group.dispersion
  )
  re_names <- design$groupef.names
  is_gaussian <- identical(family$family, "gaussian")
  if (!is_gaussian) {
    stop(
      "rlmerb() only fits gaussian() responses.",
      call. = FALSE
    )
  }
  pm <- lmerb_posterior_mean(
    design                 = design,
    measurement_prior_list = measurement_prior_list,
    tol                    = tol,
    maxit                  = maxit
  )
  icm_label <- "ICM mean"
  fixef_init <- lapply(prior$pop.prior_list, `[[`, "mu")
  names(fixef_init) <- re_names
  tau2_mode <- stats::setNames(
    vapply(re_names, function(k) prior$pop.prior_list[[k]]$dispersion, numeric(1)),
    re_names
  )
  list(
    fixef      = pm$fixef,
    fixef_init = fixef_init,
    b_mean     = pm$b_mean,
    icm_info   = list(
      converged  = pm$converged,
      iterations = pm$iterations,
      delta      = pm$delta
    ),
    icm_label  = icm_label,
    joint_mode = FALSE,
    sigma2     = prior$group.dispersion,
    tau2       = tau2_mode
  )
}

#' Point estimates for Bayesian LMMs (\code{rlmerb(..., simulate = FALSE)})
#'
#' Computes the exact joint Gaussian posterior mean (at fixed variance-component
#' plug-ins from \code{pfamily_list} / \code{dispprior_list}) without drawing
#' MCMC samples. Called from \code{\link{rlmerb}} when \code{simulate = FALSE};
#' not exported.
#'
#' @inheritParams rlmerb
#' @param offset_missing,weights_missing Logical flags from the parent
#'   \code{rlmerb()} call so design inheritance still works when those
#'   formals are forwarded. Default \code{NULL} means use \code{missing()}.
#' @return An object of class \code{c("rlmerb", "list")} with \code{popef.mode},
#'   \code{groupef.mode}, plug-in dispersion fields, full \code{prior}, thin
#'   \code{Prior}, and \code{NULL} draw / pilot / convergence slots.
#' @seealso \code{\link{rlmerb}}, \code{\link{lmerb_posterior_mean}}
#' @keywords internal
rlmerb_point <- function(
    design,
    pfamily_list,
    dispprior_list,
    offset = NULL,
    weights = 1,
    verbose = TRUE,
    print_icm_table = TRUE,
    offset_missing = NULL,
    weights_missing = NULL
) {
  cl <- match.call()
  if (is.null(offset_missing)) offset_missing <- missing(offset)
  if (is.null(weights_missing)) weights_missing <- missing(weights)

  prep <- .rlmerb_prepare_prior(
    design          = design,
    pfamily_list    = pfamily_list,
    dispprior_list  = dispprior_list,
    family          = stats::gaussian(),
    fn_name         = "rlmerb",
    offset          = offset,
    weights         = weights,
    offset_missing  = offset_missing,
    weights_missing = weights_missing
  )
  prior <- prep$prior
  disp_info <- prep$disp_info
  re_names <- prep$re_names

  icm <- .rlmerb_icm_at_fixed_vc(
    design = design,
    prior  = prior,
    family = stats::gaussian()
  )

  if (isTRUE(print_icm_table)) {
    icm_lbl <- .lmebayes_block2_icm_labels(prior, stats::gaussian())
    .lmebayes_print_icm_fixef_table(
      prior_list = prior$pop.prior_list,
      re_names   = re_names,
      fixef_icm  = icm$fixef,
      icm_info   = icm$icm_info,
      ref_label  = icm_lbl$ref_label,
      icm_label  = icm$icm_label,
      conv_label = icm_lbl$conv_label,
      header     = "--- rlmerb: population effects ---",
      verbose    = verbose
    )
  }

  n_obs <- length(design$y)
  off <- prep$ow$offset
  wts <- prep$ow$weights
  offset2 <- if (is.null(off)) {
    rep.int(0, n_obs)
  } else if (length(off) == 1L) {
    rep.int(as.numeric(off), n_obs)
  } else {
    as.numeric(off)
  }
  prior.weights <- if (length(wts) == 1L) {
    rep.int(as.numeric(wts), n_obs)
  } else {
    as.numeric(wts)
  }

  out <- list(
    call                  = cl,
    popef.mode            = icm$fixef,
    popef.init            = icm$fixef_init,
    groupef.mode          = icm$b_mean,
    popef.means           = NULL,
    popef                 = NULL,
    groupef               = NULL,
    popef.dispersion      = NULL,
    popef.dispersion.mean = NULL,
    popef.iters           = NULL,
    popef.iters.mean      = NULL,
    groupef.iters         = NULL,
    groupef.iters.mean    = NULL,
    group.dispersion      = prior$group.dispersion,
    group.dispersion.mean = prior$group.dispersion,
    group.dispersion.mode = prior$group.dispersion,
    group.dispersion.iters = NULL,
    group.dispersion.iters.mean = NULL,
    tau2.mode             = icm$tau2,
    joint_mode            = icm$joint_mode,
    icm_info              = icm$icm_info,
    m_convergence         = NULL,
    convergence           = NULL,
    convergence_info      = NULL,
    pilot                 = NULL,
    sweep_history         = NULL,
    prior                 = prior,
    Prior                 = .rlmerb_thin_Prior(prior, disp_info, prep$block1_prior),
    design                = design,
    offset                = off,
    offset2               = offset2,
    prior.weights         = prior.weights
  )
  class(out) <- c("rlmerb", "list")
  out
}
