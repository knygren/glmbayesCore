# Prior pack, route arguments, and ICM printing for rlmerb().
# Copied from lmebayesCore; the originals stay there for rglmerb().

## Internal helpers for matrix-level mixed-model samplers (rlmerb / rglmerb).
## Called by rlmerb(). lmerb() in lmebayes still reaches the lmebayesCore
## copies of the dispersion helpers with :::.

#' Resolve \code{offset}/\code{weights} for \code{rlmerb}/\code{rglmerb}.
#'
#' Formals match the matrix engines (\code{offset = NULL}, \code{weights = 1}).
#' When a formal is \emph{missing} (caller omitted it), inherit
#' \code{design$offset}/\code{design$weights} when present, else the engine
#' default. An explicitly passed value always wins (including
#' \code{weights = 1} or \code{offset = NULL} that clears a design value).
#' @noRd
.lmebayes_resolve_offset_weights <- function(
    offset,
    weights,
    design,
    offset_missing = FALSE,
    weights_missing = FALSE
) {
  if (isTRUE(offset_missing)) {
    offset <- if (!is.null(design$offset)) design$offset else NULL
  }
  if (isTRUE(weights_missing)) {
    weights <- if (!is.null(design$weights)) design$weights else 1
  } else if (is.null(weights)) {
    weights <- 1
  }
  list(offset = offset, weights = weights)
}


#' Map a matrix-engine \code{dispprior_list} to the \code{group.dispersion}
#' shapes accepted by \code{.lmebayes_resolve_group_dispersion()}.
#'
#' Engine fixed form \code{list(dispersion = ...)} unwraps to the numeric
#' value; \code{list(ddef = TRUE)} (no observation dispersion) becomes
#' \code{NULL}. Scalars, \code{dGamma()} / \code{dGamma_list()}, and other
#' list forms pass through unchanged.
#' @noRd
.lmebayes_dispprior_list_as_group_dispersion <- function(dispprior_list) {
  if (is.null(dispprior_list)) {
    return(NULL)
  }
  if (inherits(dispprior_list, "pfamily")) {
    return(dispprior_list)
  }
  if (!is.list(dispprior_list)) {
    return(dispprior_list)
  }
  nms <- names(dispprior_list)
  if (is.null(nms)) {
    return(dispprior_list)
  }
  ## Engine fixed-dispersion / GLMM "no dispersion" containers.
  if ("dispersion" %in% nms &&
      !any(c("shape", "rate", "shape_group", "rate_group") %in% nms)) {
    return(dispprior_list$dispersion)
  }
  if (isTRUE(dispprior_list$ddef) && is.null(dispprior_list$dispersion) &&
      !any(c("shape", "rate", "shape_group", "rate_group") %in% nms)) {
    return(NULL)
  }
  dispprior_list
}

#' @noRd
.lmebayes_resolve_group_dispersion <- function(
    group.dispersion,
    family,
    design = NULL,
    fn_name = "lmerb"
) {
  has_dispersion <- family$family %in%
    c("gaussian", "Gamma", "quasipoisson", "quasibinomial")

  if (!has_dispersion) {
    if (!is.null(group.dispersion)) {
      stop(
        "'group.dispersion' must be NULL for family = ", family$family,
        "() (no observation-level dispersion).",
        call. = FALSE
      )
    }
    return(list(
      mode                   = "none",
      dispersion_fix         = NULL,
      dispersion_prior_list  = NULL,
      dispersion_pfamily     = NULL
    ))
  }

  if (is.list(group.dispersion) && !inherits(group.dispersion, "pfamily")) {
    return(.lmebayes_resolve_group_dispersion_group_list(
      group.dispersion = group.dispersion,
      design           = design,
      fn_name          = fn_name
    ))
  }

  if (inherits(group.dispersion, "pfamily")) {
    if (!identical(group.dispersion$pfamily, "dGamma")) {
      stop(
        fn_name, "(): 'group.dispersion' pfamily must be dGamma(); got ",
        group.dispersion$pfamily, ". RE priors belong in 'pfamily_list'.",
        call. = FALSE
      )
    }
    pl <- group.dispersion$prior_list
    if (!isTRUE(pl$Inv_Dispersion)) {
      stop(
        fn_name, "(): dGamma() observation-dispersion prior requires ",
        "Inv_Dispersion = TRUE.",
        call. = FALSE
      )
    }
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
      stop(
        fn_name, "(): dGamma() group.dispersion prior_list requires positive ",
        "'shape' and 'rate'.",
        call. = FALSE
      )
    }
    return(list(
      mode                  = "gamma",
      dispersion_fix        = shape / rate,
      dispersion_prior_list = pl,
      dispersion_pfamily    = group.dispersion
    ))
  }

  if (is.numeric(group.dispersion) && length(group.dispersion) > 1L) {
    return(.lmebayes_resolve_group_dispersion_fixed_vector(
      group.dispersion = group.dispersion,
      design           = design,
      fn_name          = fn_name
    ))
  }

  if (is.null(group.dispersion) || !is.numeric(group.dispersion) ||
      length(group.dispersion) != 1L || !is.finite(group.dispersion) ||
      group.dispersion <= 0) {
    stop(
      "'group.dispersion' must be a single positive number, a named numeric ",
      "vector of positive per-group values, a dGamma() pfamily, or a named ",
      "list of dGamma() pfamilies (one per group) for family = ",
      family$family, "().",
      call. = FALSE
    )
  }
  list(
    mode                  = "fixed",
    dispersion_fix        = as.numeric(group.dispersion),
    dispersion_prior_list = NULL,
    dispersion_pfamily    = NULL
  )
}

#' Resolve a fixed per-group numeric vector for \code{group.dispersion}
#'
#' Fourth \code{group.dispersion} option alongside a fixed scalar, a single
#' (pooled) \code{dGamma()}, and a named list of \code{dGamma()} pfamilies: a
#' plain named numeric vector with one known, fixed dispersion value per
#' group. Unlike \code{dGamma_list(...)} (\code{"gamma_list"}), this is not a
#' prior to be sampled -- each group's \eqn{\sigma^2_j} is treated as known
#' for the duration of sampling, exactly like the pooled \code{"fixed"} mode
#' but allowed to vary by group. No \code{glmmTMB} reference fit or full-rank
#' requirement applies (there is no ING accept/reject envelope to build).
#' @noRd
.lmebayes_resolve_group_dispersion_fixed_vector <- function(
    group.dispersion,
    design,
    fn_name = "lmerb"
) {
  if (is.null(design) || is.null(design$group)) {
    stop(
      fn_name, "(): a named numeric vector for 'group.dispersion' requires ",
      "'design' with grouping information.",
      call. = FALSE
    )
  }
  group_levels <- levels(design$group)
  J <- length(group_levels)

  if (length(group.dispersion) != J) {
    stop(
      fn_name, "(): 'group.dispersion' is a vector of length ",
      length(group.dispersion), " but there are ", J, " group level(s) (",
      paste(group_levels, collapse = ", "), "). Supply exactly one fixed ",
      "dispersion value per group.",
      call. = FALSE
    )
  }
  nms <- names(group.dispersion)
  if (is.null(nms) || any(!nzchar(nms)) || !setequal(nms, group_levels)) {
    stop(
      fn_name, "(): names(group.dispersion) must match the group levels (",
      paste(group_levels, collapse = ", "), ") exactly.",
      call. = FALSE
    )
  }
  group.dispersion <- stats::setNames(
    as.numeric(group.dispersion[group_levels]),
    group_levels
  )

  if (any(!is.finite(group.dispersion)) || any(group.dispersion <= 0)) {
    stop(
      fn_name, "(): 'group.dispersion' values must all be positive and ",
      "finite.",
      call. = FALSE
    )
  }

  list(
    mode                  = "fixed_vector",
    dispersion_fix        = group.dispersion,
    dispersion_prior_list = NULL,
    dispersion_pfamily    = NULL
  )
}

#' Resolve a per-group list of \code{dGamma()} pfamilies for \code{group.dispersion}
#'
#' Third \code{group.dispersion} option alongside a fixed scalar and a single
#' (pooled) \code{dGamma()}: a named list with one \code{dGamma()} pfamily per
#' group level. Each entry keeps its own \code{shape}/\code{rate}/\code{disp_lower}/
#' \code{disp_upper} -- there is no requirement that groups share the same
#' \code{shape}/\code{rate} (\code{Prior_Setup_GLMM()} may choose to build
#' them from a shared hyperprior, but the engine itself is agnostic).
#' Requires every group to be full column rank (rank-deficient groups are not
#' yet supported for this option).
#' @noRd
.lmebayes_resolve_group_dispersion_group_list <- function(
    group.dispersion,
    design,
    fn_name = "lmerb"
) {
  if (is.null(design) || is.null(design$group)) {
    stop(
      fn_name, "(): a list of dGamma() priors for 'group.dispersion' requires ",
      "'design' with grouping information.",
      call. = FALSE
    )
  }
  group_levels <- levels(design$group)
  J <- length(group_levels)

  if (length(group.dispersion) != J) {
    stop(
      fn_name, "(): 'group.dispersion' is a list of length ",
      length(group.dispersion), " but there are ", J, " group level(s) (",
      paste(group_levels, collapse = ", "), "). Supply exactly one dGamma() ",
      "pfamily per group.",
      call. = FALSE
    )
  }
  nms <- names(group.dispersion)
  if (is.null(nms) || any(!nzchar(nms)) || !setequal(nms, group_levels)) {
    stop(
      fn_name, "(): names(group.dispersion) must match the group levels (",
      paste(group_levels, collapse = ", "), ") exactly.",
      call. = FALSE
    )
  }
  window_diagnostics <- attr(group.dispersion, "window_diagnostics")
  group.dispersion <- group.dispersion[group_levels]

  if (is.null(design$groupef.rank) || !all(design$groupef.rank[group_levels])) {
    deficient <- if (!is.null(design$groupef.rank)) {
      group_levels[!design$groupef.rank[group_levels]]
    } else {
      group_levels
    }
    stop(
      fn_name, "(): a list of per-group dGamma() dispersion priors currently ",
      "requires every group to be full column rank; rank-deficient group(s): ",
      paste(deficient, collapse = ", "), ". Use a single dGamma() or a fixed ",
      "scalar 'group.dispersion' for models with rank-deficient groups.",
      call. = FALSE
    )
  }

  shape_group      <- stats::setNames(numeric(J), group_levels)
  rate_group       <- stats::setNames(numeric(J), group_levels)
  disp_lower_group <- stats::setNames(numeric(J), group_levels)
  disp_upper_group <- stats::setNames(numeric(J), group_levels)

  for (lev in group_levels) {
    pf <- group.dispersion[[lev]]
    if (!inherits(pf, "pfamily") || !identical(pf$pfamily, "dGamma")) {
      stop(
        fn_name, "(): group.dispersion[[\"", lev, "\"]] must be a dGamma() ",
        "pfamily.",
        call. = FALSE
      )
    }
    pl <- pf$prior_list
    if (!isTRUE(pl$Inv_Dispersion)) {
      stop(
        fn_name, "(): group.dispersion[[\"", lev, "\"]] dGamma() prior ",
        "requires Inv_Dispersion = TRUE.",
        call. = FALSE
      )
    }
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
      stop(
        fn_name, "(): group.dispersion[[\"", lev, "\"]] must have positive, ",
        "finite 'shape' and 'rate'.",
        call. = FALSE
      )
    }
    lo <- pl$disp_lower
    hi <- pl$disp_upper
    if (is.null(lo) || is.null(hi) ||
        !is.numeric(lo) || !is.numeric(hi) ||
        length(lo) != 1L || length(hi) != 1L ||
        !is.finite(lo) || !is.finite(hi) ||
        lo <= 0 || hi <= lo) {
      stop(
        fn_name, "(): group.dispersion[[\"", lev, "\"]] must supply finite ",
        "'disp_lower' and 'disp_upper' with 0 < disp_lower < disp_upper -- a ",
        "list of dGamma() priors requires explicit per-group truncation bounds.",
        call. = FALSE
      )
    }
    shape_group[[lev]]      <- shape
    rate_group[[lev]]       <- rate
    disp_lower_group[[lev]] <- as.numeric(lo)
    disp_upper_group[[lev]] <- as.numeric(hi)
  }

  list(
    mode                  = "gamma_list",
    dispersion_fix        = mean(shape_group / rate_group),
    dispersion_prior_list = list(
      shape_group      = shape_group,
      rate_group       = rate_group,
      disp_lower_group = disp_lower_group,
      disp_upper_group = disp_upper_group
    ),
    dispersion_pfamily    = group.dispersion,
    window_diagnostics    = window_diagnostics
  )
}

.priors_from_pfamily_list <- function(pfamily_list,
                                      group.dispersion,
                                      design,
                                      family,
                                      fn_name = "lmerb") {

  re_names <- design$groupef.names
  p_re     <- length(re_names)

  ## --- group.dispersion (Block 1 measurement dispersion) -------------------
  disp_res <- .lmebayes_resolve_group_dispersion(
    group.dispersion = group.dispersion,
    family           = family,
    design           = design,
    fn_name          = fn_name
  )
  group.dispersion <- disp_res$dispersion_fix

  ## --- pfamily_list ---------------------------------------------------------
  if (!is.list(pfamily_list) || length(pfamily_list) != p_re) {
    stop(
      "'pfamily_list' must be a list with one pfamily per random-effect ",
      "component (", p_re, " expected: ", paste(re_names, collapse = ", "),
      "). Build it with pfamily_list(Prior_Setup_GLMM(...)).",
      call. = FALSE
    )
  }
  if (is.null(names(pfamily_list)) || !setequal(names(pfamily_list), re_names)) {
    stop(
      "Names of 'pfamily_list' must match the random-effect coefficient ",
      "names: ", paste(re_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  pfamily_list <- pfamily_list[re_names]

  prior_list <- stats::setNames(vector("list", p_re), re_names)
  tau2   <- stats::setNames(numeric(p_re), re_names)
  ptypes <- stats::setNames(character(p_re), re_names)

  for (k in re_names) {
    pf <- pfamily_list[[k]]
    if (!inherits(pf, "pfamily")) {
      stop("pfamily_list[[\"", k, "\"]] must be a pfamily object.",
           call. = FALSE)
    }
    if (!pf$pfamily %in% c("dNormal", "dIndependent_Normal_Gamma")) {
      stop(
        fn_name, "() supports only dNormal and dIndependent_Normal_Gamma ",
        "pfamilies in 'pfamily_list'; component \"", k, "\" is ",
        pf$pfamily, ".",
        call. = FALSE
      )
    }
    ptypes[[k]] <- pf$pfamily

    par_names <- colnames(design$W[[k]])
    q_k <- length(par_names)

    mu_k <- as.numeric(pf$prior_list$mu)
    if (length(mu_k) != q_k) {
      stop(
        "pfamily_list[[\"", k, "\"]]$prior_list$mu has length ",
        length(mu_k), " but the hyper design has ", q_k, " column(s): ",
        paste(par_names, collapse = ", "), ".",
        call. = FALSE
      )
    }
    mu_nms <- rownames(pf$prior_list$mu)
    if (!is.null(mu_nms) && all(nzchar(mu_nms))) {
      if (!setequal(mu_nms, par_names)) {
        stop(
          "Parameter names of pfamily_list[[\"", k, "\"]] (",
          paste(mu_nms, collapse = ", "), ") do not match the hyper design ",
          "columns (", paste(par_names, collapse = ", "), ").",
          call. = FALSE
        )
      }
      ord <- match(par_names, mu_nms)
      mu_k <- mu_k[ord]
      Sigma_k <- as.matrix(pf$prior_list$Sigma)[ord, ord, drop = FALSE]
    } else {
      Sigma_k <- as.matrix(pf$prior_list$Sigma)
    }
    names(mu_k) <- par_names
    dimnames(Sigma_k) <- list(par_names, par_names)

    ## Keep the pfamily object itself aligned with the hyper-design column
    ## order: it is passed straight to the v2 sampler as the Block 2 source
    ## of truth, so its mu/Sigma must match x_hyper[[k]].
    pfamily_list[[k]]$prior_list$mu <-
      matrix(mu_k, ncol = 1L, dimnames = list(par_names, NULL))
    pfamily_list[[k]]$prior_list$Sigma <- Sigma_k

    if (identical(pf$pfamily, "dNormal")) {
      d_k <- pf$prior_list$dispersion
      if (isTRUE(pf$prior_list$ddef)) {
        warning(
          fn_name, ": pfamily_list[[\"", k, "\"]] uses the default ",
          "dispersion = 1 (none was supplied to dNormal()); the Block 1 ",
          "random-effect variance tau^2 for \"", k, "\" is therefore 1.",
          call. = FALSE
        )
      }
    } else {
      ## ING: disp_lower/disp_upper fix the truncation window and lambda*
      ## calibration only; ICM plug-in tau^2 comes from the pfamily spec.
      d_k <- pf$prior_list$disp_lower
      if (is.null(d_k) || !is.numeric(d_k) || length(d_k) != 1L ||
          !is.finite(d_k) || d_k <= 0) {
        stop(
          fn_name, "(): pfamily_list[[\"", k, "\"]] is ",
          "dIndependent_Normal_Gamma and must supply a positive scalar ",
          "'disp_lower' (lower dispersion truncation) for lambda* calibration.",
          call. = FALSE
        )
      }
      u_k <- pf$prior_list$disp_upper
      if (is.null(u_k) || !is.numeric(u_k) || length(u_k) != 1L ||
          !is.finite(u_k) || u_k <= as.numeric(d_k)) {
        stop(
          fn_name, "(): pfamily_list[[\"", k, "\"]] is ",
          "dIndependent_Normal_Gamma and must supply a finite scalar ",
          "'disp_upper' > 'disp_lower' (upper dispersion truncation), so ",
          "the tau^2 truncation window is fixed across Gibbs sweeps. ",
          "pfamily_list(Prior_Setup_GLMM(...)) sets both bounds to the ",
          "0.01/0.99 prior dispersion quantiles by default.",
          call. = FALSE
        )
      }
    }

    tau2_k <- .two_block_tau2_plug_in_from_pfamily(pf)
    tau2[[k]] <- tau2_k
    prior_list[[k]] <- list(
      mu         = mu_k,
      Sigma      = Sigma_k,
      dispersion = tau2_k
    )
  }

  group.Sigma <- diag(unname(tau2), nrow = p_re, ncol = p_re)
  dimnames(group.Sigma) <- list(re_names, re_names)

  list(
    pfamily_list          = pfamily_list,
    group.dispersion      = group.dispersion,
    dispersion_mode       = disp_res$mode,
    dispersion_pfamily    = disp_res$dispersion_pfamily,
    dispersion_prior_list = disp_res$dispersion_prior_list,
    window_diagnostics    = disp_res$window_diagnostics,
    group.Sigma           = group.Sigma,
    pop.prior_list        = prior_list,
    ptypes         = ptypes,
    any_non_normal = any(ptypes != "dNormal")
  )
}

.lmebayes_ing_measurement_prior_list <- function(prior, disp_info, design) {
  re_names <- design$groupef.names
  p_re     <- length(re_names)
  pl       <- disp_info$dispersion_prior_list
  if (is.null(pl$shape) || is.null(pl$rate)) {
    stop(
      "dGamma() group.dispersion prior_list must contain 'shape' and 'rate'.",
      call. = FALSE
    )
  }
  mu <- matrix(
    0,
    nrow = p_re,
    ncol = 1L,
    dimnames = list(re_names, NULL)
  )
  Sigma <- as.matrix(prior$group.Sigma)
  if (nrow(Sigma) != p_re || ncol(Sigma) != p_re) {
    stop(
      "prior$group.Sigma must be ", p_re, " x ", p_re, ".",
      call. = FALSE
    )
  }
  out <- list(
    mu            = mu,
    Sigma         = Sigma,
    shape         = pl$shape,
    rate          = pl$rate,
    max_disp_perc = if (!is.null(pl$max_disp_perc)) pl$max_disp_perc else 0.99
  )
  if (!is.null(pl$disp_lower)) out$disp_lower <- pl$disp_lower
  if (!is.null(pl$disp_upper)) out$disp_upper <- pl$disp_upper
  if (is.null(out$disp_lower) || is.null(out$disp_upper)) {
    win <- .lmebayes_ing_prior_quantile_window(
      shape = as.numeric(pl$shape[1L]),
      rate  = as.numeric(pl$rate[1L])
    )
    if (is.null(out$disp_lower)) out$disp_lower <- win$disp_lower
    if (is.null(out$disp_upper)) out$disp_upper <- win$disp_upper
  }
  if (out$disp_upper <= out$disp_lower) {
    stop(
      "dGamma() measurement prior: implied disp_upper must exceed disp_lower.",
      call. = FALSE
    )
  }
  out
}

#' Build per-group ING Block~1 measurement \code{prior_list} for lmebayes glue
#'
#' Third \code{group.dispersion} option: a named list of \code{dGamma()}
#' pfamilies, one per group. Unlike \code{.lmebayes_ing_measurement_prior_list()}
#' (single pooled \code{shape}/\code{rate}), each group keeps its own
#' \code{shape}/\code{rate}/\code{disp_lower}/\code{disp_upper}.
#' @noRd
.lmebayes_ing_measurement_prior_list_group <- function(prior, disp_info, design) {
  re_names <- design$groupef.names
  p_re     <- length(re_names)
  pl       <- disp_info$dispersion_prior_list
  req <- c("shape_group", "rate_group", "disp_lower_group", "disp_upper_group")
  miss <- req[!req %in% names(pl)]
  if (length(miss)) {
    stop(
      "Internal error: per-group dGamma() group.dispersion prior_list must ",
      "contain ", paste(req, collapse = ", "), ".",
      call. = FALSE
    )
  }
  mu <- matrix(
    0,
    nrow = p_re,
    ncol = 1L,
    dimnames = list(re_names, NULL)
  )
  Sigma <- as.matrix(prior$group.Sigma)
  if (nrow(Sigma) != p_re || ncol(Sigma) != p_re) {
    stop(
      "prior$group.Sigma must be ", p_re, " x ", p_re, ".",
      call. = FALSE
    )
  }
  list(
    mu               = mu,
    Sigma            = Sigma,
    shape_group      = pl$shape_group,
    rate_group       = pl$rate_group,
    disp_lower_group = pl$disp_lower_group,
    disp_upper_group = pl$disp_upper_group,
    max_disp_perc    = if (!is.null(pl$max_disp_perc)) pl$max_disp_perc else 0.99
  )
}

#' Assemble the flat argument list for the routed LMM \code{rLMM_reg} export
#'
#' @description
#' Shared front-door helper that turns \code{design} + the \code{prior}
#' returned by \code{.priors_from_pfamily_list} + \code{disp_info}
#' (from \code{.lmebayes_resolve_group_dispersion()}) into the flat,
#' named argument list passed via \code{do.call()} to whichever
#' \code{rLMM_reg} export \code{.lmebayes_reg_route_fn()} resolves to. There
#' are exactly four possible targets, chosen by \code{disp_info$mode} and
#' \code{prior$any_non_normal}: \code{\link{rLMMNormal_reg_known_vcov}},
#' \code{rLMMNormal_reg_estimated_vcov}, \code{rLMMindepNormalGamma_reg_known_vcov},
#' or \code{rLMMindepNormalGamma_reg_estimated_vcov} (the actual dispatch
#' happens one level up, in \code{.lmebayes_run_lmm_engine()} via
#' \code{.lmebayes_reg_route_key()}/\code{.lmebayes_reg_route_fn()} --
#' this function only builds the argument list, it does not pick the route).
#'
#' @details
#' Concretely, this function does three things:
#' \enumerate{
#'   \item Builds a base argument list that is always the same shape:
#'     the response/design matrices and grouping structure taken from
#'     \code{design}, \code{prior$pfamily_list} unchanged (the routed export
#'     derives its own Block~2 random-effect prior precision from
#'     \code{pfamily_list} internally, so it is not built here), and
#'     the naming/control arguments (\code{tv_tol}, \code{progbar},
#'     \code{verbose}), with \code{group_name} attached to \code{group} as
#'     an attribute rather than a separate argument.
#'   \item Builds \code{args$dispprior_list} (the Block~1/measurement prior),
#'     whose \emph{shape} depends on \code{disp_info$mode}: a plain known
#'     \code{dispersion} value for \code{"none"}/\code{"fixed"}/
#'     \code{"fixed_vector"}, or a shared/per-group ING measurement prior
#'     (built by \code{.lmebayes_ing_measurement_prior_list()}/
#'     \code{.lmebayes_ing_measurement_prior_list_group()}) for
#'     \code{"gamma"}/\code{"gamma_list"}.
#'   \item Conditionally adds route-specific controls: pilot-stage controls
#'     (\code{gap_tol}, \code{mode_gap_max}, \code{diag_sweeps},
#'     \code{stage_verbose}) when \code{prior$any_non_normal} is
#'     \code{TRUE} (an ING route, which needs a pilot chain), or
#'     \code{sim_method} when the route is fixed \emph{and} known-vcov
#'     (the only route with more than one sampling engine).
#' }
#' It currently combines all of the above (base argument assembly,
#' route-specific \code{dispprior_list} branching, and conditional pilot-stage/
#' \code{sim_method} fields) into a single flat list and is a refactor
#' candidate; its argument list and return shape are not yet considered
#' stable.
#'
#' @param n Number of stored draws.
#' @param design A \code{model_setup} object.
#' @param prior The list returned by \code{.priors_from_pfamily_list}.
#' @param disp_info Resolved dispersion info, as returned by
#'   \code{.lmebayes_resolve_group_dispersion()} (must supply
#'   \code{mode} and \code{dispersion_fix}).
#' @param tv_tol Total-variation tolerance passed through to the routed
#'   export.
#' @param progbar,verbose Passed through to the routed export.
#' @param gap_tol,mode_gap_max,diag_sweeps Pilot-stage controls, added to the
#'   argument list only when \code{prior$any_non_normal} is \code{TRUE}
#'   (i.e. an ING route).
#' @param sim_method Sampling engine (\code{"DEFAULT"} or
#'   \code{"TWO_BLOCK_GIBBS"}), added to the argument list only on the
#'   fixed+known-vcov (\code{rLMMNormal_reg_known_vcov}) route.
#' @return A flat named list (\code{args}) suitable for
#'   \code{do.call(route$export_fn, args)}, with elements:
#'   \describe{
#'     \item{\code{n}}{Number of stored draws, unchanged from the \code{n}
#'       argument.}
#'     \item{\code{y}}{\code{design$y}, the response vector.}
#'     \item{\code{D}}{\code{design$D}, the level-1 (\eqn{l_2 \times p_{re}})
#'       random-effect design matrix.}
#'     \item{\code{group}}{\code{design$group}, the grouping factor, with
#'       \code{attr(group, "group_name")} set to \code{design$group_name}
#'       (the routed export has no \code{group_name} formal and resolves it
#'       from this attribute, since \code{design$group} is never a bare
#'       variable at the routed export's call site).}
#'     \item{\code{W}}{\code{design$W}, the named list of
#'       group-level hyper-design matrices (one per random-effect
#'       coefficient).}
#'     \item{\code{pfamily_list}}{\code{prior$pfamily_list} unchanged. The
#'       routed export derives its own Block~2 random-effect prior
#'       precision from this internally; it is not built or passed here.}
#'     \item{\code{dispprior_list}}{Always present; shape depends on
#'       \code{disp_info$mode}: for \code{"none"}/\code{"fixed"}/
#'       \code{"fixed_vector"} it is \code{list(dispersion = disp_info$dispersion_fix)}
#'       (a scalar or named per-group numeric vector -- a \emph{known}, fixed
#'       observation dispersion); for \code{"gamma"} it is the pooled ING
#'       measurement prior list from
#'       \code{.lmebayes_ing_measurement_prior_list()}; for
#'       \code{"gamma_list"} it is the per-group ING measurement prior list
#'       from \code{.lmebayes_ing_measurement_prior_list_group()}.}
#'     \item{\code{tv_tol}}{The \code{tv_tol} argument, unchanged.}
#'     \item{\code{progbar}, \code{verbose}}{The \code{progbar}/\code{verbose}
#'       arguments, unchanged.}
#'     \item{\code{gap_tol}, \code{mode_gap_max}, \code{diag_sweeps},
#'       \code{stage_verbose}}{Present \strong{only} when
#'       \code{prior$any_non_normal} is \code{TRUE} (an ING/estimated-vcov
#'       route that needs a pilot chain): the \code{gap_tol}/
#'       \code{mode_gap_max}/\code{diag_sweeps} arguments unchanged, and
#'       \code{stage_verbose} set to \code{verbose}. Absent for
#'       \code{"none"}/\code{"fixed"}/\code{"fixed_vector"} dispersion
#'       modes.}
#'     \item{\code{sim_method}}{Present \strong{only} when
#'       \code{prior$any_non_normal} is \code{FALSE} \emph{and}
#'       \code{disp_info$mode} is not \code{"gamma"}/\code{"gamma_list"}
#'       (i.e. only on the route that resolves to
#'       \code{\link{rLMMNormal_reg_known_vcov}}, the only export with a
#'       real \code{sim_method} dispatch between exact iid draws and
#'       two-block Gibbs sweeps): the \code{sim_method} argument, unchanged.
#'       Every other route ignores \code{sim_method} (if it even accepts
#'       it), so it is only forwarded here.}
#'   }
#' @keywords internal
matrix_args_lmm <- function(
    n,
    design,
    prior,
    disp_info,
    tv_tol        = 0.01,
    progbar       = TRUE,
    verbose       = FALSE,
    gap_tol       = 0.0196,
    mode_gap_max  = 1.0,
    diag_sweeps   = FALSE,
    sim_method    = "DEFAULT",
    offset        = NULL,
    weights       = 1,
    offset_missing = FALSE,
    weights_missing = FALSE
) {
  ## The routed export has no 'group_name' formal: attach it to 'group'
  ## itself (design$group is never a bare variable here, so the export's
  ## substitute()-based fallback could not resolve it anyway).
  grp <- design$group
  attr(grp, "group_name") <- design$group_name
  ow <- .lmebayes_resolve_offset_weights(
    offset, weights, design,
    offset_missing = offset_missing,
    weights_missing = weights_missing
  )

  args <- list(
    n             = n,
    y             = design$y,
    D             = design$D,
    group         = grp,
    W             = design$W,
    pfamily_list  = prior$pfamily_list,
    offset        = ow$offset,
    weights       = ow$weights,
    tv_tol        = tv_tol,
    progbar       = progbar,
    verbose       = verbose
  )

  if (identical(disp_info$mode, "gamma")) {
    args$dispprior_list <- .lmebayes_ing_measurement_prior_list(
      prior     = prior,
      disp_info = disp_info,
      design    = design
    )
  } else if (identical(disp_info$mode, "gamma_list")) {
    args$dispprior_list <- .lmebayes_ing_measurement_prior_list_group(
      prior     = prior,
      disp_info = disp_info,
      design    = design
    )
  } else {
    args$dispprior_list <- list(dispersion = disp_info$dispersion_fix)
  }

  if (isTRUE(prior$any_non_normal)) {
    args$gap_tol       <- gap_tol
    args$mode_gap_max  <- mode_gap_max
    args$diag_sweeps   <- diag_sweeps
    args$stage_verbose <- verbose
  } else if (!identical(disp_info$mode, "gamma") &&
             !identical(disp_info$mode, "gamma_list")) {
    ## lmm_fixed_known route only: rLMMNormal_reg_known_vcov() is the only
    ## export with a real sim_method dispatch (exact iid vs. two-block
    ## Gibbs). Every other route accepts the argument (if at all) as a
    ## no-op, so it is only forwarded here.
    args$sim_method <- sim_method
  }

  args
}

#' @noRd
.lmebayes_run_lmm_engine <- function(
    n,
    design,
    prior,
    disp_info,
    tv_tol        = 0.01,
    progbar       = TRUE,
    verbose       = FALSE,
    gap_tol             = 0.0196,
    mode_gap_max        = 1.0,
    diag_sweeps         = FALSE,
    sim_method          = "DEFAULT",
    offset              = NULL,
    weights             = 1,
    offset_missing      = FALSE,
    weights_missing     = FALSE
) {
  route_key <- .lmebayes_reg_route_key(
    family         = gaussian(),
    disp_mode      = disp_info$mode,
    any_non_normal = prior$any_non_normal
  )
  route <- .lmebayes_reg_route_fn(route_key)
  args  <- matrix_args_lmm(
    n             = n,
    design        = design,
    prior         = prior,
    disp_info     = disp_info,
    tv_tol        = tv_tol,
    progbar       = progbar,
    verbose       = verbose,
    gap_tol       = gap_tol,
    mode_gap_max  = mode_gap_max,
    diag_sweeps   = diag_sweeps,
    sim_method    = sim_method,
    offset        = offset,
    weights       = weights,
    offset_missing = offset_missing,
    weights_missing = weights_missing
  )
  out <- do.call(route$export_fn, args)
  if (is.null(out$convergence_info$sim_method_used)) {
    out$convergence_info$sim_method_used <- "TWO_BLOCK_GIBBS"
  }
  .lmebayes_attach_group_dispersion(out, disp_info)
}

#' Build a P/Sigma-free Block~1 prior list (\code{dispersion}/\code{ddef} only)
#'
#' The Block~1 random-effect prior precision is now always derived
#' internally (from \code{pfamily_list}) by the routed \code{rGLMM_reg}/
#' \code{rLMM_reg} exports, which reject a caller-supplied \code{P}/
#' \code{Sigma}; this helper therefore no longer computes or returns one
#' (previously \code{solve(measurement_prior_list$group.Sigma)}).
#' @noRd
.lmebayes_block1_prior_list <- function(
    measurement_prior_list,
    group.dispersion = NULL
) {
  dispersion <- if (!is.null(group.dispersion)) {
    group.dispersion
  } else {
    measurement_prior_list$group.dispersion
  }
  if (is.null(dispersion)) {
    list(ddef = TRUE)
  } else {
    list(dispersion = dispersion, ddef = FALSE)
  }
}

#' Attach \code{group.dispersion} / \code{group.dispersion.mean} from mode + draws.
#'
#' Observation residual variance \eqn{\sigma^2} / \eqn{\sigma^2_j} (not RE
#' \eqn{\tau^2}). Fixed modes return a scalar or named length-\code{J} vector;
#' \code{dGamma()} returns length-\code{n} draws; a per-group
#' \code{dGamma()} list returns an \code{n x J} matrix. Families without
#' observation-level dispersion get \code{NULL}. Clears legacy
#' \code{dispersion_ranef}/\code{sigma2} slots if present.
#' @noRd
.lmebayes_attach_group_dispersion <- function(out, disp_info) {
  mode <- disp_info$mode
  ## Drop legacy public names if an older path left them on.
  out$dispersion_ranef <- NULL
  out$dispersion_ranef.mean <- NULL
  out$sigma2 <- NULL
  out$sigma2.mean <- NULL

  if (identical(mode, "none")) {
    out$group.dispersion <- NULL
    out$group.dispersion.mean <- NULL
    return(out)
  }
  if (identical(mode, "fixed")) {
    val <- as.numeric(disp_info$dispersion_fix)
    out$group.dispersion <- val
    out$group.dispersion.mean <- val
    return(out)
  }
  if (identical(mode, "fixed_vector")) {
    val <- disp_info$dispersion_fix
    out$group.dispersion <- val
    out$group.dispersion.mean <- val
    return(out)
  }
  if (identical(mode, "gamma")) {
    dr <- out$group.dispersion
    if (is.null(dr)) {
      stop(
        "Internal error: dGamma measurement dispersion requires ",
        "'group.dispersion' draws on the sampler output.",
        call. = FALSE
      )
    }
    out$group.dispersion <- as.numeric(dr)
    out$group.dispersion.mean <- mean(out$group.dispersion)
    return(out)
  }
  if (identical(mode, "gamma_list")) {
    dr <- out$group.dispersion
    if (is.null(dr)) {
      stop(
        "Internal error: per-group dGamma() measurement dispersion requires ",
        "'group.dispersion' draws on the sampler output.",
        call. = FALSE
      )
    }
    out$group.dispersion <- as.matrix(dr)
    out$group.dispersion.mean <- colMeans(out$group.dispersion)
    return(out)
  }
  stop("Unknown dispersion mode: ", mode, call. = FALSE)
}

#' @noRd
.lmebayes_add_popef_summaries <- function(x) {
  if (!is.null(x$popef)) {
    x$popef.means <- lapply(x$popef, colMeans)
  }
  if (!is.null(x$popef.dispersion)) {
    x$popef.dispersion.mean <- colMeans(x$popef.dispersion)
  }
  if (!is.null(x$popef.iters) && !is.null(x$m_convergence)) {
    x$popef.iters.mean <- colMeans(x$popef.iters) / x$m_convergence
  }
  if (!is.null(x$groupef.iters) && !is.null(x$m_convergence)) {
    x$groupef.iters.mean <- mean(x$groupef.iters) / x$m_convergence
  }
  if (!is.null(x$group.dispersion.iters) && !is.null(x$m_convergence)) {
    x$group.dispersion.iters.mean <-
      colMeans(x$group.dispersion.iters) / x$m_convergence
  }
  x
}

.lmebayes_block2_icm_labels <- function(prior, family = gaussian()) {
  any_ing <- isTRUE(prior$any_non_normal)
  is_gauss <- is.null(family) || identical(family$family, "gaussian")
  ref_label <- "prior mean"
  if (any_ing) {
    icm_label   <- "gamma @ lmer tau2"
    icm_verbose <- "Block 2 start at lmer tau^2 plug-in"
    conv_label  <- "Plug-in fixed point"
  } else if (is_gauss) {
    icm_label   <- "ICM mean"
    icm_verbose <- "ICM posterior mean"
    conv_label  <- "ICM"
  } else {
    icm_label   <- "ICM mode"
    icm_verbose <- "ICM posterior mode"
    conv_label  <- "ICM"
  }
  list(
    ref_label   = ref_label,
    icm_label   = icm_label,
    icm_verbose = icm_verbose,
    conv_label  = conv_label
  )
}

#' @noRd
.lmebayes_print_icm_fixef_table <- function(
    prior_list,
    re_names,
    fixef_icm,
    icm_info,
    ref_label,
    icm_label,
    conv_label = "ICM",
    header,
    verbose
) {
  if (!isTRUE(verbose) || is.null(fixef_icm)) {
    return(invisible(NULL))
  }
  fixef_ref <- lapply(prior_list, `[[`, "mu")
  names(fixef_ref) <- re_names
  hdr <- sprintf("  %-18s  %-30s  %14s  %18s",
                 "RE component", "parameter", ref_label, icm_label)
  sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
  cat(header, "\n")
  cat(hdr, "\n")
  cat(sep, "\n")
  for (k in re_names) {
    nms_k  <- names(fixef_ref[[k]])
    ref_v  <- fixef_ref[[k]]
    icm_v  <- fixef_icm[[k]]
    for (nm in nms_k) {
      cat(sprintf("  %-18s  %-30s  %14.4f  %18.4f\n",
                  k, nm, ref_v[[nm]], icm_v[[nm]]))
    }
  }
  if (!is.null(icm_info)) {
    cat(sprintf("  (%s converged: %s, %d iter, delta = %.2e)\n\n",
                conv_label,
                icm_info$converged, icm_info$iterations, icm_info$delta))
  } else {
    cat("\n")
  }
  invisible(NULL)
}

#' @noRd
.lmebayes_print_fixef_init <- function(
    fixef_init,
    re_names,
    verbose,
    header = "--- main-stage fixef.init (pilot colMeans) ---"
) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  cat(header, "\n")
  for (k in re_names) {
    for (nm in names(fixef_init[[k]])) {
      cat(sprintf("  %-18s  %-30s  %12.4f\n",
                  k, nm, fixef_init[[k]][[nm]]))
    }
  }
  cat("\n")
  invisible(NULL)
}
