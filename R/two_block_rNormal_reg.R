#' @noRd
.two_block_normalize_family <- function(family) {
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    stop("'family' not recognized.", call. = FALSE)
  }

  okfamilies <- c(
    "gaussian", "poisson", "binomial",
    "quasipoisson", "quasibinomial", "Gamma"
  )
  if (!family$family %in% okfamilies) {
    stop(
      "family \"", family$family, "\" is not supported by two-block samplers.",
      call. = FALSE
    )
  }

  oklinks <- switch(
    family$family,
    gaussian = "identity",
    poisson = "log",
    quasipoisson = "log",
    binomial = c("logit", "probit", "cloglog"),
    quasibinomial = c("logit", "probit", "cloglog"),
    Gamma = "log",
    character(0)
  )
  if (!family$link %in% oklinks) {
    stop(
      "link \"", family$link, "\" not available for family \"",
      family$family, "\".",
      call. = FALSE
    )
  }

  family
}

#' @noRd
.two_block_validate_block1_prior <- function(prior_list_block1, family) {
  if (!is.list(prior_list_block1)) {
    stop("'prior_list_block1' must be a list.", call. = FALSE)
  }
  if (is.null(prior_list_block1$P) && is.null(prior_list_block1$Sigma)) {
    stop("prior_list_block1 must contain 'P' or 'Sigma'.", call. = FALSE)
  }

  ddef <- if ("ddef" %in% names(prior_list_block1)) {
    prior_list_block1$ddef
  } else {
    is.null(prior_list_block1$dispersion)
  }

  dispersion <- prior_list_block1$dispersion

  if (identical(family$family, "gaussian")) {
    if (is.null(dispersion)) {
      stop(
        "prior_list_block1 must contain 'dispersion' for gaussian() Block~1.",
        call. = FALSE
      )
    }
    if (isTRUE(ddef)) {
      stop(
        "For gaussian() Block~1, dNormal() requires an explicit dispersion.",
        call. = FALSE
      )
    }
  }

  if (family$family %in% c("gaussian", "Gamma") && isTRUE(ddef)) {
    stop(
      "For gaussian() and Gamma() models, dNormal() requires an explicit dispersion.",
      call. = FALSE
    )
  }

  list(ddef = ddef, dispersion = dispersion)
}

#' Resolve group_name from group's own attribute or the caller's expression
#'
#' @description
#' \code{group_name} cannot be derived from \code{group}'s \emph{value} (R
#' variable names are not part of an object's data). It is resolved from,
#' in order: (1) \code{attr(group, "group_name")}, set internally when
#' \code{group} is threaded through a call chain that no longer has a
#' \code{group_name} formal to carry it (e.g. \code{lmerb()} via
#' \code{do.call()}, or the private-helper calls into
#' \code{\link{rLMMNormal_joint_iid}});
#' (2) \code{substitute(group)} in the caller's own frame, when \code{group}
#' was passed as a bare variable.
#' @noRd
.lmebayes_resolve_group_name <- function(group, group_expr, fn_name) {
  gn <- attr(group, "group_name", exact = TRUE)
  if (!is.null(gn) && length(gn) == 1L && !is.na(gn) && nzchar(gn)) {
    return(as.character(gn))
  }
  if (is.symbol(group_expr)) {
    return(deparse(group_expr))
  }
  stop(
    fn_name, "(): cannot derive 'group_name': 'group' is neither a plain ",
    "variable nor carries a 'group_name' attribute. Pass 'group' as a bare ",
    "variable (e.g. group = school_id), or attach one via ",
    "attr(group, \"group_name\") <- \"school_id\" beforehand.",
    call. = FALSE
  )
}


#' Validate a Block 2 pfamily list
#'
#' Checks that each component is a supported \code{pfamily} object and that
#' the type-specific hyperparameters needed by the two-block driver are present.
#' When \code{J} is supplied (sampling context), ING components must carry
#' \emph{both} truncation bounds (\code{disp_lower}, \code{disp_upper}) so
#' the tau^2_k window is fixed across all inner Gibbs sweeps, and must also
#' satisfy the prior-vs-data balance guard \code{n_prior <= J}: the Block 2
#' hyper-regression has \code{J} group-level observations, and the ING
#' dispersion envelope caps its log-tilt at the data contribution
#' \code{n_w/2 = J/2} (Remark 4.1.3 of the ING vignette), which presumes a
#' likelihood-dominated prior.  Under the calibration
#' \code{shape = (n_prior + 1 + q_k)/2} this is
#' \code{2*shape - 1 - q_k <= J} (equivalently \code{pwt_disp <= 0.5}).
#'
#' @param pfamily_list Named list of \code{pfamily} objects.
#' @param re_names Character vector of RE component names (defines order).
#' @param J Number of groups in the Block 2 hyper-regression, or \code{NULL}
#'   to skip the sampling-specific prior-vs-data guard (calibration-only
#'   contexts such as \code{two_block_rate_from_pfamily_list}).
#' @return The validated list, reordered to match \code{re_names}.
#' @noRd
.two_block_validate_pfamily_list <- function(pfamily_list, re_names, J = NULL) {
  if (!is.list(pfamily_list)) {
    stop("'pfamily_list' must be a named list of pfamily objects.",
         call. = FALSE)
  }
  if (!setequal(names(pfamily_list), re_names)) {
    stop("names(pfamily_list) must match groupef.names.", call. = FALSE)
  }
  pfamily_list <- pfamily_list[re_names]

  supported <- c("dNormal", "dIndependent_Normal_Gamma")
  for (k in re_names) {
    pf <- pfamily_list[[k]]
    if (!inherits(pf, "pfamily") ||
        is.null(pf$pfamily) || is.null(pf$prior_list)) {
      stop(
        "pfamily_list[[\"", k, "\"]] must be a pfamily object ",
        "(e.g. dNormal() or dIndependent_Normal_Gamma()).",
        call. = FALSE
      )
    }
    if (!pf$pfamily %in% supported) {
      stop(
        "pfamily_list[[\"", k, "\"]]: unsupported pfamily \"", pf$pfamily,
        "\" (allowed: ", paste(supported, collapse = ", "), ").",
        call. = FALSE
      )
    }
    pl <- pf$prior_list
    if (identical(pf$pfamily, "dNormal")) {
      if (is.null(pl$dispersion) || isTRUE(pl$ddef)) {
        stop(
          "pfamily_list[[\"", k, "\"]]: dNormal() requires an explicit ",
          "dispersion for the Block 2 (gaussian) update.",
          call. = FALSE
        )
      }
    } else {
      if (is.null(pl$shape) || is.null(pl$rate) ||
          !is.numeric(pl$shape) || !is.numeric(pl$rate) ||
          pl$shape[1L] <= 0 || pl$rate[1L] <= 0) {
        stop(
          "pfamily_list[[\"", k, "\"]]: dIndependent_Normal_Gamma() requires ",
          "positive 'shape' and 'rate'.",
          call. = FALSE
        )
      }
      if (is.null(pl$disp_lower) || !is.numeric(pl$disp_lower) ||
          pl$disp_lower[1L] <= 0) {
        stop(
          "pfamily_list[[\"", k, "\"]]: dIndependent_Normal_Gamma() requires ",
          "a positive 'disp_lower' (initial/calibration tau^2_k plug-in).",
          call. = FALSE
        )
      }
      if (!is.null(J)) {
        ## Sampling context: both truncation bounds are required so the
        ## tau^2_k window is fixed across all inner Gibbs sweeps (with
        ## either bound missing the envelope would re-derive a surrogate
        ## posterior window per sweep, making the truncation
        ## state-dependent and the disp_lower-based rate bound only
        ## approximate).
        if (is.null(pl$disp_upper) || !is.numeric(pl$disp_upper) ||
            !is.finite(pl$disp_upper[1L]) ||
            pl$disp_upper[1L] <= pl$disp_lower[1L]) {
          stop(
            "pfamily_list[[\"", k, "\"]]: dIndependent_Normal_Gamma() ",
            "sampling requires a finite 'disp_upper' > 'disp_lower' so the ",
            "tau^2_k truncation window is fixed across Gibbs sweeps ",
            "(e.g. the 0.99 prior dispersion quantile ",
            "1/qgamma(0.01, shape, rate)).",
            call. = FALSE
          )
        }
        q_k <- length(pl$mu)
        .ing_stop_if_prior_exceeds_data(
          shape       = pl$shape[1L],
          p           = q_k,
          n_w         = J,
          detail      = paste0(
            "the Block 2 hyper-regression has only J = ", J, " groups"
          ),
          limit_label = "J",
          prefix      = paste0("pfamily_list[[\"", k, "\"]]: ")
        )
      }
    }
  }
  pfamily_list
}

#' Summarize Block~2 \code{pfamily} types for fixed vs estimated RE dispersion
#'
#' When every component is \code{dNormal}, each RE-scale dispersion
#' \eqn{\tau^2_k} is fixed and the random-effects covariance structure
#' (via \code{P}) is treated as known up to the hyper means \eqn{\gamma_k}.
#' When any component is \code{dIndependent_Normal_Gamma}, at least one
#' \eqn{\tau^2_k} is sampled and fed back into Block~1 each inner sweep.
#'
#' @param pfamily_list Named list of validated \code{pfamily} objects.
#' @return List with \code{ptypes}, \code{any_non_normal}, and \code{all_dNormal}.
#' @noRd
.two_block_summarize_pfamily_list <- function(pfamily_list) {
  ptypes <- vapply(
    pfamily_list,
    function(pf) pf$pfamily,
    character(1)
  )
  names(ptypes) <- names(pfamily_list)
  any_non_normal <- any(ptypes != "dNormal")
  list(
    ptypes         = ptypes,
    any_non_normal = any_non_normal,
    all_dNormal    = !any_non_normal
  )
}
