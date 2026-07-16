#' Matrix-level replicate-chain Gibbs engines for Bayesian LMMs
#'
#' @description
#' Gaussian two-block LMM samplers at matrix level (\code{y}, \code{Z},
#' \code{x_hyper}, \code{pfamily_list}). Each stored draw runs
#' \code{m_convergence} inner Gibbs sweeps (Block~1 random effects, then
#' Block~2 hyperparameters). There is no standalone \code{rLMM()} function;
#' formula-level fitting lives in **lmebayes** (\code{lmerb}/\code{rlmerb});
#' non-Gaussian models use the \code{\link{rGLMM_reg}} routes.
#'
#' @section Four route engines:
#' \describe{
#'   \item{\code{rLMMNormal_reg_known_vcov}}{
#'     Fixed observation \eqn{\sigma^2}; all Block~2 components \code{dNormal}
#'     (known \eqn{\tau^2_k}).}
#'   \item{\code{rLMMNormal_reg_estimated_vcov}}{
#'     Fixed \eqn{\sigma^2}; at least one ING Block~2 component.}
#'   \item{\code{rLMMindepNormalGamma_reg_known_vcov}}{
#'     Random \eqn{\sigma^2} (per-group ING Block~1); all Block~2 \code{dNormal}.}
#'   \item{\code{rLMMindepNormalGamma_reg_estimated_vcov}}{
#'     Random \eqn{\sigma^2}; at least one ING Block~2 component.}
#' }
#'
#' @section Dispatchers:
#' \code{\link{rLMMNormal_reg}} and \code{\link{rLMMindepNormalGamma_reg}}
#' validate inputs and delegate to the appropriate route (or, for the legacy
#' outer-loop \code{rLMMindepNormalGamma_reg}, draw \eqn{\sigma^2} then call
#' \code{rLMMNormal_reg}).
#'
#' @param n Number of stored draws. If \code{length(n) > 1}, the length is used.
#' @param y Response vector of length \code{nrow(x)}.
#' @param x Level-1 design matrix \code{Z} (\code{l2 x p_re}).
#' @param block Grouping factor or block partition of length \code{l2}.
#' @param x_hyper Named list of group-level design matrices (\code{J x q_k}),
#'   one per column of \code{x}.
#' @param P Random-effect prior precision matrix (\code{p_re x p_re}).
#' @param prior_list Block~1 prior: \code{list(dispersion = sigma2)} for fixed
#'   \eqn{\sigma^2} routes, \code{dGamma()} fields for legacy
#'   \code{rLMMindepNormalGamma_reg}, or shared ING measurement prior
#'   (\code{mu}, \code{Sigma}, \code{shape}, \code{rate}, \ldots) for ING routes
#'   (plug-in \eqn{\sigma^2 =} \code{shape/rate} for ICM/TV is derived internally).
#' @param pfamily_list Named list of Block~2 \code{pfamily} objects.
#' @param icm_tol,icm_maxit ICM convergence controls for the internal Block~2 start.
#' @param tv_tol Total-variation tolerance in \code{(0, 1)} for calibration.
#'   Inner Gibbs sweeps per stored draw (\code{m_convergence}) are derived from
#'   Theorem~3 at the ICM Block~2 start; pilot chain counts likewise.
#' @param re_coef_names Names for columns of \code{x}.
#' @param group_levels Character vector defining row order of Block~1 draws.
#' @param group_name Name for the grouping column in \code{coefficients}.
#' @param progbar Show a text progress bar during sampling.
#' @param verbose Print convergence calibration / ICM lines.
#' @param gap_tol,mode_gap_max,diag_sweeps,stage_verbose Pilot-stage controls for
#'   \code{rLMMNormal_reg_estimated_vcov} and ING estimated routes (see route docs).
#' @family simfuncs
#' @seealso \code{\link{rGLMM_reg}}, \code{\link{rindepNormalGamma_reg}}
#' @name rLMM_reg
NULL

#' Shared matrix-level validation for LMM replicate-chain engines
#' @noRd
.rLMM_validate_matrix_inputs <- function(
    n,
    y,
    x,
    x_hyper,
    tv_tol,
    re_coef_names,
    group_levels,
    group_name,
    block
) {
  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
      !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
    stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
  }

  y <- as.vector(y)
  x <- as.matrix(x)
  l2 <- nrow(x)
  if (length(y) != l2) {
    stop("length(y) must equal nrow(x).", call. = FALSE)
  }

  if (is.null(re_coef_names) || length(re_coef_names) != ncol(x)) {
    re_coef_names <- if (ncol(x) >= 1L) {
      cn <- colnames(x)
      if (is.null(cn) || length(cn) != ncol(x)) {
        paste0("RE", seq_len(ncol(x)))
      } else {
        cn
      }
    } else {
      stop("'x' must have at least one column.", call. = FALSE)
    }
  }
  colnames(x) <- re_coef_names
  re_names <- re_coef_names

  group_levels <- as.character(group_levels)
  if (length(group_levels) < 1L) {
    stop("'group_levels' must contain at least one level.", call. = FALSE)
  }

  if (is.null(group_name) || !nzchar(group_name)) {
    group_name <- tryCatch(
      deparse(substitute(block))[1L],
      error = function(e) "group"
    )
    if (!nzchar(group_name)) group_name <- "group"
  }

  if (!is.list(x_hyper) || is.data.frame(x_hyper)) {
    stop("'x_hyper' must be a list of design matrices.", call. = FALSE)
  }
  if (length(x_hyper) != length(re_names)) {
    stop("length(x_hyper) must equal ncol(x) = ", length(re_names), ".",
         call. = FALSE)
  }
  if (!setequal(names(x_hyper), re_names)) {
    x_hyper <- x_hyper[re_names]
  }

  list(
    n             = n,
    y             = y,
    x             = x,
    x_hyper       = x_hyper,
    tv_tol        = tv_tol,
    re_names      = re_names,
    group_levels  = group_levels,
    group_name    = group_name
  )
}

#' @noRd
.rLMM_validate_P <- function(P, p_re, fn_name = "rLMMNormal_reg") {
  P <- as.matrix(P)
  if (!is.matrix(P) || nrow(P) != p_re || ncol(P) != p_re) {
    stop(
      fn_name, "(): 'P' must be a ", p_re, " x ", p_re, " matrix.",
      call. = FALSE
    )
  }
  if (!isSymmetric(P)) {
    stop(fn_name, "(): 'P' must be symmetric.", call. = FALSE)
  }
  P
}

#' @noRd
.rLMM_validate_fixed_dispersion_prior_list <- function(
    prior_list,
    fn_name = "rLMMNormal_reg"
) {
  if (!is.list(prior_list) || is.null(prior_list$dispersion)) {
    stop(
      fn_name, "(): 'prior_list' must contain 'dispersion' (fixed sigma^2).",
      call. = FALSE
    )
  }
  d <- prior_list$dispersion
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d <= 0) {
    stop(
      fn_name, "(): 'prior_list$dispersion' must be a single positive number.",
      call. = FALSE
    )
  }
  extra <- setdiff(
    names(prior_list),
    c("dispersion", "ddef")
  )
  if (length(extra)) {
    stop(
      fn_name, "(): 'prior_list' must contain fixed dispersion only; ",
      "unexpected fields: ", paste(extra, collapse = ", "), ".",
      call. = FALSE
    )
  }
  as.numeric(d)
}

#' @noRd
.rLMM_validate_dGamma_dispersion_prior_list <- function(
    prior_list,
    fn_name = "rLMMindepNormalGamma_reg"
) {
  if (!is.list(prior_list)) {
    stop(fn_name, "(): 'prior_list' must be a list.", call. = FALSE)
  }
  req <- c("shape", "rate", "beta", "Inv_Dispersion")
  miss <- req[!req %in% names(prior_list)]
  if (length(miss)) {
    stop(
      fn_name, "(): 'prior_list' must contain ",
      paste(req, collapse = ", "), " (from dGamma()).",
      call. = FALSE
    )
  }
  if (!isTRUE(prior_list$Inv_Dispersion)) {
    stop(
      fn_name, "(): dGamma() observation-dispersion prior requires ",
      "Inv_Dispersion = TRUE.",
      call. = FALSE
    )
  }
  shape <- prior_list$shape
  rate  <- prior_list$rate
  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape) ||
      shape <= 0) {
    stop(fn_name, "(): 'prior_list$shape' must be a positive scalar.",
         call. = FALSE)
  }
  if (!is.numeric(rate) || length(rate) != 1L || !is.finite(rate) ||
      rate <= 0) {
    stop(fn_name, "(): 'prior_list$rate' must be a positive scalar.",
         call. = FALSE)
  }
  beta <- as.matrix(prior_list$beta)
  if (nrow(beta) != 1L || ncol(beta) != 1L) {
    stop(
      fn_name, "(): 'prior_list$beta' must be a 1 x 1 matrix for ",
      "observation-level dispersion.",
      call. = FALSE
    )
  }
  if (!is.null(prior_list$disp_lower) && !is.null(prior_list$disp_upper)) {
    if (prior_list$disp_upper <= prior_list$disp_lower) {
      stop(
        fn_name, "(): 'prior_list$disp_upper' must exceed 'disp_lower'.",
        call. = FALSE
      )
    }
  }
  prior_list
}

#' Plug-in observation \eqn{\sigma^2} from a dGamma / ING measurement \code{prior_list}
#' @noRd
.rLMM_dispersion_fix_from_prior_list <- function(
    prior_list,
    fn_name = "rLMM_reg"
) {
  if (!is.list(prior_list)) {
    stop(fn_name, "(): 'prior_list' must be a list.", call. = FALSE)
  }
  if (!is.null(prior_list$shape_group) || !is.null(prior_list$rate_group)) {
    shape_group <- as.numeric(prior_list$shape_group)
    rate_group  <- as.numeric(prior_list$rate_group)
    if (length(shape_group) < 1L || length(shape_group) != length(rate_group) ||
        any(!is.finite(shape_group)) || any(shape_group <= 0) ||
        any(!is.finite(rate_group)) || any(rate_group <= 0)) {
      stop(
        fn_name, "(): 'prior_list$shape_group' and 'prior_list$rate_group' ",
        "must be positive, finite, and of the same length.",
        call. = FALSE
      )
    }
    # Pooled reference (mean of per-group prior means); used only as an ICM
    # plug-in starting point, not as the sampler's actual per-group prior.
    return(mean(shape_group / rate_group))
  }
  if (is.null(prior_list$shape) || is.null(prior_list$rate)) {
    stop(
      fn_name, "(): 'prior_list' must contain 'shape' and 'rate' ",
      "(plug-in sigma^2 = shape / rate is derived internally).",
      call. = FALSE
    )
  }
  shape <- as.numeric(prior_list$shape[1L])
  rate  <- as.numeric(prior_list$rate[1L])
  if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    stop(
      fn_name, "(): 'prior_list$shape' and 'prior_list$rate' must be positive scalars.",
      call. = FALSE
    )
  }
  as.numeric(shape / rate)
}

#' Observation-level linear predictor from group random effects
#' @noRd
.rLMM_observation_mu <- function(x, block, b_mat, group_levels) {
  g_chr <- as.character(block)
  g_idx <- match(g_chr, group_levels)
  if (anyNA(g_idx)) {
    stop("block levels not found in 'group_levels'.", call. = FALSE)
  }
  b_obs <- b_mat[g_idx, , drop = FALSE]
  rowSums(x * b_obs)
}

#' Group-level random-effect matrix from \code{coefficients} output
#' @noRd
.rLMM_b_matrix_from_coefficients <- function(
    coef_df,
    re_names,
    group_levels,
    group_name
) {
  b_mat <- matrix(
    NA_real_,
    nrow = length(group_levels),
    ncol = length(re_names),
    dimnames = list(group_levels, re_names)
  )
  for (lev in group_levels) {
    rows <- coef_df[[group_name]] == lev
    if (!any(rows)) {
      stop("coefficients missing group level: ", lev, call. = FALSE)
    }
    hit <- which(rows)[1L]
    b_mat[lev, ] <- as.numeric(coef_df[hit, re_names, drop = TRUE])
  }
  b_mat
}

#' ICM start for matrix-level LMM engines
#'
#' Iterated conditional modes for Block~2 hyperparameters at fixed Block~1
#' dispersion; used internally by \code{rLMMNormal_reg_*_vcov}.
#' @noRd
.rLMM_icm_at_start <- function(
    y,
    x,
    block,
    x_hyper,
    prior_list_block1,
    pfamily_list,
    re_names,
    group_levels,
    group_name,
    icm_tol,
    icm_maxit,
    verbose,
    engine_label
) {
  design_icm <- list(
    y             = y,
    Z             = x,
    groups        = factor(block, levels = group_levels),
    X_hyper       = x_hyper,
    re_coef_names = re_names,
    group_name    = group_name
  )
  icm <- .two_block_icm_at_start(
    design       = design_icm,
    prior_list   = prior_list_block1,
    pfamily_list = pfamily_list,
    re_names     = re_names,
    family       = gaussian(),
    tol          = icm_tol,
    maxit        = icm_maxit
  )
  if (isTRUE(verbose)) {
    cat(sprintf(
      "  %s: ICM posterior mean (converged: %s, %d iter, delta = %.2e)\n\n",
      engine_label,
      icm$icm$converged,
      icm$icm$iterations,
      icm$icm$delta
    ))
  }
  icm
}

#' Labels and \code{convergence_info$method} for rate calibration plug-ins
#' @noRd
.rLMM_rate_calibration_meta <- function(
    any_non_normal,
    random_measurement = FALSE
) {
  if (isTRUE(random_measurement) && isTRUE(any_non_normal)) {
    list(
      method = "disp_upper_bound+disp_lower_bound",
      label  = paste0(
        "conservative: measurement disp_upper + ",
        "RE disp_lower plug-ins"
      )
    )
  } else if (isTRUE(random_measurement)) {
    list(
      method = "disp_upper_bound",
      label  = "conservative: measurement disp_upper plug-in"
    )
  } else if (isTRUE(any_non_normal)) {
    list(
      method = "disp_lower_bound",
      label  = "conservative: non-dNormal RE dispersion (disp_lower plug-in)"
    )
  } else {
    list(method = "exact", label = "exact")
  }
}

#' Calibrate inner Gibbs sweeps for matrix-level LMM engines
#'
#' Computes the two-block convergence rate at the chain start, then sets
#' \code{m_convergence} to at least \code{l_for_tv(tv_tol) + 1} (Theorem 3).
#' @noRd
.rLMM_calibrate_m_convergence <- function(
    x,
    block,
    x_hyper,
    prior_list_block1,
    pfamily_list,
    group_levels,
    tv_tol,
    any_non_normal,
    random_measurement = FALSE,
    engine_label,
    verbose
) {
  rate <- two_block_rate_from_pfamily_list(
    x                 = x,
    block             = block,
    x_hyper           = x_hyper,
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
    family            = gaussian(),
    group_levels      = group_levels
  )
  m_min <- .two_block_cap_inner_sweeps(
    two_block_l_for_tv(rate, tv_tol, method = "theorem3") + 1L
  )
  m_convergence <- m_min
  calib_meta <- .rLMM_rate_calibration_meta(
    any_non_normal     = any_non_normal,
    random_measurement = random_measurement
  )
  if (isTRUE(verbose)) {
    cat(sprintf(
      paste0(
        "--- %s: convergence calibration [%s]: lambda* = %.4f, ",
        "tv_tol = %g => m_min = %d, using m_convergence = %d ---\n\n"
      ),
      engine_label, calib_meta$label, rate$lambda_star, tv_tol, m_min,
      m_convergence
    ))
  }
  list(
    m_convergence     = m_convergence,
    convergence_info  = list(
      method        = calib_meta$method,
      tv_tol        = tv_tol,
      lambda_star   = rate$lambda_star,
      eigenvalues   = rate$eigenvalues,
      m_min         = m_min,
      m_convergence = m_convergence
    )
  )
}
#' Validate a shared ING measurement prior for per-group Block~1 updates
#'
#' Accepts either a single pooled \code{shape}/\code{rate}/\code{disp_lower}/
#' \code{disp_upper} (scalar) prior, or a per-group prior carrying
#' \code{shape_group}/\code{rate_group}/\code{disp_lower_group}/
#' \code{disp_upper_group} named vectors (one entry per \code{group_levels}).
#' The two are mutually exclusive; which one is present determines the
#' dispersion mode used by the Block~1 sweep (see
#' \code{.two_block_block1_ing_one_chain()}).
#' @noRd
.rLMM_validate_ing_measurement_prior_list <- function(
    prior_list,
    p_re,
    fn_name = "rLMMindepNormalGamma_reg",
    group_levels = NULL
) {
  if (!is.list(prior_list)) {
    stop(fn_name, "(): 'prior_list' must be a list.", call. = FALSE)
  }
  req_common <- c("mu", "Sigma")
  miss <- req_common[!req_common %in% names(prior_list)]
  if (length(miss)) {
    stop(
      fn_name, "(): 'prior_list' must contain ",
      paste(req_common, collapse = ", "),
      " (from dIndependent_Normal_Gamma()).",
      call. = FALSE
    )
  }
  mu <- as.matrix(prior_list$mu, ncol = 1L)
  Sigma <- as.matrix(prior_list$Sigma)
  if (nrow(mu) != p_re || ncol(mu) != 1L) {
    stop(
      fn_name, "(): 'prior_list$mu' must be a ", p_re, " x 1 matrix.",
      call. = FALSE
    )
  }
  if (nrow(Sigma) != p_re || ncol(Sigma) != p_re) {
    stop(
      fn_name, "(): 'prior_list$Sigma' must be ", p_re, " x ", p_re, ".",
      call. = FALSE
    )
  }

  is_group <- !is.null(prior_list$shape_group) || !is.null(prior_list$rate_group)
  if (is_group) {
    if (is.null(group_levels)) {
      stop(
        fn_name, "(): 'group_levels' is required to validate a per-group ",
        "measurement dispersion prior_list.",
        call. = FALSE
      )
    }
    req_grp <- c("shape_group", "rate_group", "disp_lower_group", "disp_upper_group")
    miss <- req_grp[!req_grp %in% names(prior_list)]
    if (length(miss)) {
      stop(
        fn_name, "(): per-group 'prior_list' must contain ",
        paste(req_grp, collapse = ", "), ".",
        call. = FALSE
      )
    }
    for (nm in req_grp) {
      v <- prior_list[[nm]]
      if (!is.numeric(v) || is.null(names(v)) ||
          !setequal(names(v), group_levels)) {
        stop(
          fn_name, "(): 'prior_list$", nm, "' must be a numeric vector ",
          "named by every group level (", paste(group_levels, collapse = ", "),
          ").",
          call. = FALSE
        )
      }
      if (any(!is.finite(v)) || any(v <= 0)) {
        stop(
          fn_name, "(): 'prior_list$", nm, "' must be positive and finite ",
          "for every group.",
          call. = FALSE
        )
      }
    }
    if (any(prior_list$disp_upper_group[group_levels] <=
            prior_list$disp_lower_group[group_levels])) {
      stop(
        fn_name, "(): 'prior_list$disp_upper_group' must exceed ",
        "'disp_lower_group' for every group.",
        call. = FALSE
      )
    }
    return(prior_list)
  }

  req <- c("shape", "rate")
  miss <- req[!req %in% names(prior_list)]
  if (length(miss)) {
    stop(
      fn_name, "(): 'prior_list' must contain ",
      paste(req, collapse = ", "),
      " (from dIndependent_Normal_Gamma()).",
      call. = FALSE
    )
  }
  shape <- prior_list$shape
  rate  <- prior_list$rate
  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape) ||
      shape <= 0) {
    stop(fn_name, "(): 'prior_list$shape' must be a positive scalar.",
         call. = FALSE)
  }
  if (!is.numeric(rate) || length(rate) != 1L || !is.finite(rate) ||
      rate <= 0) {
    stop(fn_name, "(): 'prior_list$rate' must be a positive scalar.",
         call. = FALSE)
  }
  if (!is.null(prior_list$disp_lower) && !is.null(prior_list$disp_upper)) {
    if (prior_list$disp_upper <= prior_list$disp_lower) {
      stop(
        fn_name, "(): 'prior_list$disp_upper' must exceed 'disp_lower'.",
        call. = FALSE
      )
    }
  }
  if (is.null(prior_list$disp_upper) || !is.numeric(prior_list$disp_upper) ||
      length(prior_list$disp_upper) != 1L || !is.finite(prior_list$disp_upper) ||
      prior_list$disp_upper <= 0) {
    stop(
      fn_name, "(): 'prior_list$disp_upper' is required (conservative ",
      "measurement-dispersion plug-in for lambda* calibration).",
      call. = FALSE
    )
  }
  prior_list
}

#' Block~1 Gaussian prior list (\code{P}, \code{dispersion})
#' @noRd
.rLMM_block1_prior_gaussian <- function(P, dispersion) {
  list(P = P, dispersion = dispersion, ddef = FALSE)
}

#' Conservative measurement \eqn{\sigma^2} plug-in for rate calibration
#' @noRd
.rLMM_measurement_disp_upper_for_rate <- function(ing_prior_list, fn_name) {
  d <- ing_prior_list$disp_upper
  if (is.null(d) || !is.numeric(d) || length(d) != 1L || !is.finite(d) ||
      d <= 0) {
    stop(
      fn_name, "(): 'prior_list$disp_upper' is required for lambda* ",
      "calibration when measurement dispersion is random.",
      call. = FALSE
    )
  }
  as.numeric(d)
}

#' Conservative measurement-dispersion inputs for lambda* calibration
#'
#' Scalar \code{shape}/\code{rate} priors plug in a single \code{disp_upper};
#' per-group priors instead build a per-observation precision \code{weights}
#' vector (\code{1 / disp_upper_group[group]}), so that groups with a wider
#' truncation window do not drag down the conservative rate bound for groups
#' with a tighter one. \code{dispersion_scalar} is a pooled fallback, used
#' only to satisfy generic \code{prior_list_block1} validation and as the
#' fixed \code{dispersion} plug-in when \code{weights} is not used downstream.
#' @noRd
.rLMM_measurement_rate_inputs <- function(
    ing_prior_list,
    block,
    group_levels,
    fn_name
) {
  if (!is.null(ing_prior_list$disp_upper_group)) {
    du <- ing_prior_list$disp_upper_group
    miss <- setdiff(group_levels, names(du))
    if (length(miss)) {
      stop(
        fn_name, "(): 'prior_list$disp_upper_group' is missing group ",
        "level(s): ", paste(miss, collapse = ", "), ".",
        call. = FALSE
      )
    }
    w <- 1 / as.numeric(du[as.character(block)])
    if (any(!is.finite(w))) {
      stop(
        fn_name, "(): non-finite per-observation weights derived from ",
        "'prior_list$disp_upper_group'.",
        call. = FALSE
      )
    }
    return(list(
      weights           = w,
      dispersion_scalar = mean(as.numeric(du))
    ))
  }
  list(
    weights           = NULL,
    dispersion_scalar = .rLMM_measurement_disp_upper_for_rate(ing_prior_list, fn_name)
  )
}

#' Truncated Gamma draw for ING measurement dispersion (prior-only path)
#' @noRd
.rLMM_ing_sample_sigma2 <- function(pl_j) {
  shape <- as.numeric(pl_j$shape[1L])
  rate  <- as.numeric(pl_j$rate[1L])
  lo    <- pl_j$disp_lower
  hi    <- pl_j$disp_upper
  if (is.null(lo)) {
    lo <- as.numeric(stats::qgamma(0.01, shape = shape, rate = rate))
  }
  if (is.null(hi)) {
    hi <- as.numeric(stats::qgamma(0.99, shape = shape, rate = rate))
  }
  F_lo <- stats::pgamma(lo, shape = shape, rate = rate)
  F_hi <- stats::pgamma(hi, shape = shape, rate = rate)
  stats::qgamma(
    F_lo + stats::runif(1L) * (F_hi - F_lo),
    shape = shape,
    rate  = rate
  )
}

#' Prior-only ING draw for one group when \code{Z_j} is rank-deficient
#' @noRd
.rLMM_ing_prior_draw_one_group <- function(pl_j, re_names) {
  sigma2 <- .rLMM_ing_sample_sigma2(pl_j)
  mu     <- as.numeric(pl_j$mu[, 1L])
  names(mu) <- re_names
  Sigma  <- as.matrix(pl_j$Sigma)
  L <- tryCatch(
    chol(Sigma),
    error = function(e) chol(Sigma + 1e-8 * diag(nrow(Sigma)))
  )
  b <- mu + sqrt(sigma2) * as.numeric(crossprod(L, stats::rnorm(length(mu))))
  names(b) <- re_names
  list(coefficients = b, dispersion = sigma2, iters = 1L)
}

#' Subset an ING \code{prior_list} to identifiable \code{Z_j} columns
#' Diagnostic-only guard for the disp_lower/disp_upper crash trace
#'
#' Fires \emph{only} on an invalid value (never on a legitimate \code{NULL},
#' which is a valid "unset" state upstream); tagged with the function name
#' where the invalid value was first observed, to bisect which hop in the
#' R -> C++ call chain first sees a bad value.
#' @noRd
.lmebayes_check_disp_bounds_or_stop <- function(disp_lower, disp_upper, fn_name) {
  if (is.null(disp_lower) || is.null(disp_upper)) {
    return(invisible(NULL))
  }
  if (!is.numeric(disp_lower) || !is.numeric(disp_upper) ||
      length(disp_lower) < 1L || length(disp_upper) < 1L ||
      !is.finite(disp_lower[[1L]]) || !is.finite(disp_upper[[1L]]) ||
      disp_lower[[1L]] <= 0 || disp_upper[[1L]] <= 0 ||
      disp_upper[[1L]] <= disp_lower[[1L]]) {
    stop(sprintf(
      "invalid disp_lower or disp_upper in function %s. disp_lower=%s, disp_upper=%s",
      fn_name, paste(disp_lower, collapse = ", "), paste(disp_upper, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(NULL)
}

#' @noRd
.rLMM_ing_prior_list_subset <- function(pl_j, keep, re_names) {
  keep <- as.integer(keep)
  re_keep <- re_names[keep]
  list(
    mu            = pl_j$mu[re_keep, , drop = FALSE],
    Sigma         = pl_j$Sigma[re_keep, re_keep, drop = FALSE],
    shape         = pl_j$shape,
    rate          = pl_j$rate,
    max_disp_perc = pl_j$max_disp_perc,
    disp_lower    = pl_j$disp_lower,
    disp_upper    = pl_j$disp_upper
  )
}

#' One school/group Block~1 ING draw (wrapper around \code{rindepNormalGamma_reg})
#'
#' Called once per factor level inside \code{.two_block_block1_ing_one_chain}.
#' Returns a **fixed contract** for the batch driver:
#' \code{list(coefficients = named b_j, dispersion = sigma2_j, iters = ...)}.
#'
#' Why not call \code{rindepNormalGamma_reg()} directly in the school loop?
#' \itemize{
#'   \item \code{rindepNormalGamma_reg} expects a full-rank design and returns an
#'     \code{rglmb}-style object (\code{coefficients} matrix, \code{Prior}, etc.).
#'     The sweep driver only needs one named \code{b_j} row plus scalar
#'     \code{sigma2_j} and envelope iteration count.
#'   \item Per-school \code{Z_j} can be rank-deficient (too few pupils, collinear
#'     RE columns). The ING envelope then fails or is undefined on the full
#'     \code{p_re} columns; this helper routes deficient schools to prior-only or
#'     QR-identifiable subsets before calling the sampler.
#'   \item Optional \code{full_rank} (from \code{design$re_rank}) skips repeated
#'     QR when rank was computed once at setup.
#' }
#'
#' Three paths (mutually exclusive after rank is resolved):
#' \describe{
#'   \item[Path A ? standard]{Full column rank and \code{n_j >= p_re}: one
#'     \code{rindepNormalGamma_reg(n = 1)} on all RE columns.}
#'   \item[Path B ? prior only]{\code{rank(Z_j) == 0}: no likelihood information;
#'     draw \code{(b_j, sigma2_j)} from the ING prior only
#'     (\code{.rLMM_ing_prior_draw_one_group}).}
#'   \item[Path C ? partial rank]{Some but not all columns identified: ING on the
#'     QR-identifiable columns, then impute the remaining coordinates from the
#'     prior given the drawn \code{sigma2_j}.}
#' }
#'
#' @noRd
.rLMM_ing_one_group_draw <- function(
    y_j,
    Z_j,
    pl_j,
    re_names,
    family,
    full_rank = NULL
) {
  p_re <- length(re_names)
  n_j  <- nrow(Z_j)

  .lmebayes_check_disp_bounds_or_stop(
    pl_j$disp_lower, pl_j$disp_upper, ".rLMM_ing_one_group_draw"
  )

  # Resolve effective rank of this school's design (cached or QR).
  rk   <- if (isTRUE(full_rank)) {
    p_re
  } else if (isFALSE(full_rank)) {
    0L
  } else {
    as.integer(Matrix::rankMatrix(Z_j, method = "qr")[1L])
  }

  # Path A: data identify all RE columns ? direct ING envelope on full Z_j.
  if (rk >= p_re && n_j >= p_re) {
    ing <- rindepNormalGamma_reg(
      n            = 1L,
      y            = y_j,
      x            = Z_j,
      prior_list   = pl_j,
      family       = family,
      progbar      = FALSE,
      verbose      = FALSE,
      use_parallel = FALSE
    )
    return(list(
      coefficients = ing$coefficients[1L, , drop = TRUE],
      dispersion   = as.numeric(ing$dispersion[1L]),
      iters        = as.numeric(ing$iters[1L])
    ))
  }

  # Path B: no identifiable column ? skip likelihood; prior draw only.
  if (rk < 1L) {
    return(.rLMM_ing_prior_draw_one_group(pl_j, re_names))
  }

  # Path C: rank 1 .. p_re-1 ? ING on identifiable columns, impute the rest.
  qr_j  <- qr(Z_j)
  keep  <- qr_j$pivot[seq_len(qr_j$rank)]
  pl_sub <- .rLMM_ing_prior_list_subset(pl_j, keep, re_names)
  .lmebayes_check_disp_bounds_or_stop(
    pl_sub$disp_lower, pl_sub$disp_upper, ".rLMM_ing_one_group_draw (Path C)"
  )
  ing <- rindepNormalGamma_reg(
    n            = 1L,
    y            = y_j,
    x            = Z_j[, keep, drop = FALSE],
    prior_list   = pl_sub,
    family       = family,
    progbar      = FALSE,
    verbose      = FALSE,
    use_parallel = FALSE
  )

  mu     <- as.numeric(pl_j$mu[, 1L])
  names(mu) <- re_names
  Sigma  <- as.matrix(pl_j$Sigma)
  sigma2 <- as.numeric(ing$dispersion[1L])
  b_full <- mu
  coef_sub <- ing$coefficients[1L, , drop = TRUE]
  b_full[keep] <- as.numeric(coef_sub)
  non <- setdiff(seq_len(p_re), keep)
  if (length(non)) {
    for (k in non) {
      b_full[k] <- mu[k] + sqrt(sigma2 * Sigma[k, k]) * stats::rnorm(1L)
    }
  }
  names(b_full) <- re_names

  list(
    coefficients = b_full,
    dispersion   = sigma2,
    iters        = as.numeric(ing$iters[1L])
  )
}

#' Build a group-level \code{prior_list} for \code{rindepNormalGamma_reg()}
#' @noRd
.rLMM_ing_prior_list_for_group <- function(
    ing_prior_list,
    mu_j,
    re_names
) {
  .lmebayes_check_disp_bounds_or_stop(
    ing_prior_list$disp_lower, ing_prior_list$disp_upper,
    ".rLMM_ing_prior_list_for_group (entry)"
  )
  mu_j <- as.numeric(mu_j)
  if (is.null(names(mu_j)) || !setequal(names(mu_j), re_names)) {
    names(mu_j) <- re_names
  }
  mu_mat <- matrix(mu_j[re_names], ncol = 1L,
                   dimnames = list(re_names, NULL))
  out <- list(
    mu            = mu_mat,
    Sigma         = ing_prior_list$Sigma,
    shape         = ing_prior_list$shape,
    rate          = ing_prior_list$rate,
    max_disp_perc = if (!is.null(ing_prior_list$max_disp_perc)) {
      ing_prior_list$max_disp_perc
    } else {
      0.99
    }
  )
  if (!is.null(ing_prior_list$disp_lower)) {
    out$disp_lower <- ing_prior_list$disp_lower
  }
  if (!is.null(ing_prior_list$disp_upper)) {
    out$disp_upper <- ing_prior_list$disp_upper
  }
  .lmebayes_check_disp_bounds_or_stop(
    out$disp_lower, out$disp_upper, ".rLMM_ing_prior_list_for_group (exit)"
  )
  out
}

#' Build Block~1 ING \code{prior_list} for one chain (mu_all from prep)
#' @noRd
.two_block_block1_ing_prior_list_one_chain <- function(mu_all, ing_prior_list) {
  prior_list <- list(
    mu            = mu_all,
    Sigma         = ing_prior_list$Sigma,
    shape         = ing_prior_list$shape,
    rate          = ing_prior_list$rate,
    max_disp_perc = if (!is.null(ing_prior_list$max_disp_perc)) {
      ing_prior_list$max_disp_perc
    } else {
      0.99
    }
  )
  if (!is.null(ing_prior_list$disp_lower)) {
    prior_list$disp_lower <- ing_prior_list$disp_lower
  }
  if (!is.null(ing_prior_list$disp_upper)) {
    prior_list$disp_upper <- ing_prior_list$disp_upper
  }
  prior_list
}

#' Map \code{BlockEnvelopeSim} block draws to the batch \code{b} matrix layout
#' @noRd
.two_block_block1_map_envelope_sim_to_b <- function(sim, re_names, group_levels, p_re) {
  ids <- vapply(sim$block_results, function(br) as.character(br$block_id), character(1))
  b_rows <- lapply(sim$block_results, function(br) {
    v <- as.numeric(br$beta[, 1L])
    if (length(v) != p_re) {
      stop(
        "BlockEnvelopeSim beta length (", length(v),
        ") must equal ncol(Z) / length(re_names) (", p_re, ").",
        call. = FALSE
      )
    }
    stats::setNames(v, re_names)
  })
  b_mat <- do.call(rbind, b_rows)
  rownames(b_mat) <- ids
  colnames(b_mat) <- re_names
  if (!all(group_levels %in% ids)) {
    stop(
      "BlockEnvelopeSim block ids do not cover all group levels.",
      call. = FALSE
    )
  }
  b_draw <- b_mat[group_levels, re_names, drop = FALSE]
  rownames(b_draw) <- group_levels
  b_draw
}

#' Block~1 one chain: block ING envelope update via \code{rIndepNormalGammaRegBlock}
#' (Centering ? Build ? DispersionBuild ? Sim in one C++ call).
#' @noRd
.two_block_block1_envelope_draw_one_chain <- function(
    y,
    Z,
    groups,
    prior_list,
    p_re,
    re_names,
    group_levels
) {
  nobs <- length(y)
  # Joint sampler (gs^k product-face table ? only feasible for very few groups;
  # for k groups with gs faces each the table has gs^k entries, exponential in k):
  # .rIndepNormalGammaRegBlock_cpp(
  #   n             = 1L,
  #   y             = y,
  #   x             = Z,
  #   block         = groups,
  #   prior_list    = prior_list,
  #   prior_lists   = NULL,
  #   offset        = rep(0, nobs),
  #   wt            = rep(1, nobs),
  #   p_re          = p_re,
  #   n_rss_iter    = 10L,
  #   Gridtype      = 3L,
  #   n_envopt      = -1L,
  #   RSS_ML        = NA_real_,
  #   use_parallel  = TRUE,
  #   use_opencl    = FALSE,
  #   progbar       = FALSE,
  #   verbose       = FALSE,
  #   group_levels  = group_levels,
  #   re_names      = re_names
  # )[c("b", "dispersion_ranef", "iters_mean")]
  .rIndepNormalGammaRegBlockInd_cpp(
    n             = 1L,
    y             = y,
    x             = Z,
    block         = groups,
    prior_list    = prior_list,
    prior_lists   = NULL,
    offset        = rep(0, nobs),
    wt            = rep(1, nobs),
    p_re          = p_re,
    n_rss_iter    = 10L,
    Gridtype      = 3L,
    n_envopt      = -1L,
    RSS_ML        = NA_real_,
    use_parallel  = TRUE,
    use_opencl    = FALSE,
    progbar       = FALSE,
    verbose       = FALSE,
    group_levels  = group_levels,
    re_names      = re_names
  )[c("b", "dispersion_ranef", "iters_mean")]
}

#' Block~1 one chain: independent per-group ING draws for a list of
#' per-group \code{dGamma()} measurement dispersion priors.
#'
#' Reuses \code{.rLMM_ing_one_group_draw()} (rank paths A/B/C) once per group,
#' each with its own \code{shape_group}/\code{rate_group}/\code{disp_lower_group}/
#' \code{disp_upper_group}. Unlike the shared-\eqn{\sigma^2} joint envelope
#' path, groups are conditionally independent given the current Block~2
#' fixef/tau2 (no joint acceptance across groups).
#' @noRd
.two_block_block1_ing_group_draw_one_chain <- function(
    y,
    Z,
    groups,
    ing_prior_list,
    mu_all,
    re_names,
    group_levels,
    family,
    full_rank = TRUE
) {
  p_re  <- length(re_names)
  g_chr <- as.character(groups)

  b_mat <- matrix(
    NA_real_,
    nrow = length(group_levels),
    ncol = p_re,
    dimnames = list(group_levels, re_names)
  )
  dispersion_ranef <- stats::setNames(numeric(length(group_levels)), group_levels)
  iters_j <- numeric(length(group_levels))

  for (j in seq_along(group_levels)) {
    lev  <- group_levels[j]
    rows <- which(g_chr == lev)

    pl_group_j <- list(
      Sigma         = ing_prior_list$Sigma,
      shape         = ing_prior_list$shape_group[[lev]],
      rate          = ing_prior_list$rate_group[[lev]],
      max_disp_perc = ing_prior_list$max_disp_perc,
      disp_lower    = ing_prior_list$disp_lower_group[[lev]],
      disp_upper    = ing_prior_list$disp_upper_group[[lev]]
    )
    .lmebayes_check_disp_bounds_or_stop(
      pl_group_j$disp_lower, pl_group_j$disp_upper,
      sprintf(".two_block_block1_ing_group_draw_one_chain (group '%s')", lev)
    )
    pl_j <- .rLMM_ing_prior_list_for_group(
      ing_prior_list = pl_group_j,
      mu_j           = mu_all[, lev],
      re_names       = re_names
    )
    draw_j <- .rLMM_ing_one_group_draw(
      y_j       = y[rows],
      Z_j       = Z[rows, , drop = FALSE],
      pl_j      = pl_j,
      re_names  = re_names,
      family    = family,
      full_rank = full_rank
    )

    b_mat[lev, ]            <- draw_j$coefficients[re_names]
    dispersion_ranef[[lev]] <- draw_j$dispersion
    iters_j[j]              <- draw_j$iters
  }

  list(
    b                = b_mat,
    dispersion_ranef = dispersion_ranef,
    iters_mean       = mean(iters_j),
    iters_group      = stats::setNames(iters_j, group_levels)
  )
}

#' Block~1 one chain: block ING envelope update for all groups
#'
#' Dispatches on the shape of \code{ing_prior_list}: a shared (pooled)
#' \code{shape}/\code{rate} routes through the joint C++ envelope
#' (\code{.two_block_block1_envelope_draw_one_chain()}, one shared
#' \eqn{\sigma^2} for all groups); a per-group \code{shape_group}/
#' \code{rate_group} (from a list of \code{dGamma()} pfamilies, one per
#' group) routes through \code{.two_block_block1_ing_group_draw_one_chain()}
#' instead (independent \eqn{\sigma^2_j} per group).
#' @noRd
.two_block_block1_ing_one_chain <- function(
    batch,
    i,
    design,
    block1_prior,
    ing_prior_list,
    family,
    ptypes
) {
  prep <- .two_block_block1_prep_one_chain(
    batch              = batch,
    i                  = i,
    design             = design,
    block1_prior       = block1_prior,
    ptypes             = ptypes,
    use_cpp_mu_all     = FALSE,
    use_cpp_prior_tau2 = FALSE
  )
  mu_all <- prep$mu_all
  p_re <- length(batch$re_names)
  re_names <- batch$re_names
  group_levels <- batch$group_levels

  y <- design$y
  Z <- as.matrix(design$Z)

  if (!is.null(ing_prior_list$shape_group)) {
    out <- .two_block_block1_ing_group_draw_one_chain(
      y              = y,
      Z              = Z,
      groups         = design$groups,
      ing_prior_list = ing_prior_list,
      mu_all         = mu_all,
      re_names       = re_names,
      group_levels   = group_levels,
      family         = family
    )
  } else {
    prior_list <- .two_block_block1_ing_prior_list_one_chain(
      mu_all, ing_prior_list
    )

    out <- .two_block_block1_envelope_draw_one_chain(
      y            = y,
      Z            = Z,
      groups       = design$groups,
      prior_list   = prior_list,
      p_re         = p_re,
      re_names     = re_names,
      group_levels = group_levels
    )
  }

  if (isTRUE(getOption("glmbayesCore.debug_block1_ing_levels", FALSE))) {
    cat(sprintf(
      "  [Block1 ING debug] chain %d/%d: envelope sim dispersion=%s\n",
      i, batch$n, paste(signif(out$dispersion_ranef, 6), collapse = ", ")
    ))
    utils::flush.console()
  }

  out
}

#' Block~1 batch: ING per-group updates for all replicate chains
#' @noRd
.two_block_block1_ing_all_chains <- function(
    n,
    fixef,
    tau2,
    b,
    iters_ranef,
    re_names,
    group_levels,
    design,
    block1_prior,
    ing_prior_list,
    family,
    ptypes,
    iters_ranef_group = NULL,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE
) {
  show_bar <- isTRUE(progbar) && n > 1L
  batch <- list(
    n            = n,
    fixef        = fixef,
    tau2         = tau2,
    re_names     = re_names,
    group_levels = group_levels,
    b            = .two_block_ensure_batch_b_dimnames(b, group_levels, re_names, n),
    iters_ranef  = iters_ranef + 0
  )
  group_mode <- !is.null(ing_prior_list$shape_group)
  dispersion_ranef <- if (group_mode) {
    matrix(
      NA_real_, nrow = n, ncol = length(group_levels),
      dimnames = list(NULL, group_levels)
    )
  } else {
    numeric(n)
  }
  if (group_mode) {
    if (is.null(iters_ranef_group)) {
      iters_ranef_group <- matrix(
        0,
        nrow = n,
        ncol = length(group_levels),
        dimnames = list(NULL, group_levels)
      )
    } else {
      iters_ranef_group <- iters_ranef_group + 0
    }
  }
  debug_b1 <- isTRUE(getOption("glmbayesCore.debug_block1_ing_levels", FALSE))

  for (i in seq_len(n)) {
    if (show_bar) .two_block_progress_bar(i, n, prefix = progbar_prefix)
    if (debug_b1 && nzchar(progbar_prefix)) {
      cat(progbar_prefix, "chain ", i, "/", n, " Block1 ING\n", sep = "")
      utils::flush.console()
    }
    out <- .two_block_block1_ing_one_chain(
      batch          = batch,
      i              = i,
      design         = design,
      block1_prior   = block1_prior,
      ing_prior_list = ing_prior_list,
      family         = family,
      ptypes         = ptypes
    )
    batch$b <- .two_block_batch_b_assign_slice(
      batch$b, i, out$b, use_cpp_b_slice = FALSE
    )
    batch$iters_ranef <- .two_block_batch_iters_ranef_add(
      batch$iters_ranef, i, out$iters_mean, use_cpp_iters_ranef_add = FALSE
    )
    if (group_mode) {
      dispersion_ranef[i, ] <- out$dispersion_ranef[group_levels]
      iters_ranef_group[i, ] <- iters_ranef_group[i, ] +
        out$iters_group[group_levels]
    } else {
      dispersion_ranef[i] <- out$dispersion_ranef
    }
  }
  if (show_bar) {
    .two_block_progress_bar_finish(newline = progbar_finish_newline)
  }

  list(
    b                = .two_block_ensure_batch_b_dimnames(
      batch$b, group_levels, re_names, n
    ),
    iters_ranef      = batch$iters_ranef,
    iters_ranef_group = if (group_mode) iters_ranef_group else NULL,
    dispersion_ranef = dispersion_ranef
  )
}

#' Pack ING sweep outputs (includes per-chain \code{dispersion_ranef} from envelope sim)
#' @noRd
.rGLMM_sweep_save_ing <- function(
    n,
    fixef,
    tau2,
    b,
    iters,
    iters_ranef,
    dispersion_ranef,
    re_names,
    group_levels,
    design,
    collect_block1 = TRUE,
    iters_ranef_group = NULL
) {
  out <- .rGLMM_sweep_save(
    n              = n,
    fixef          = fixef,
    tau2           = tau2,
    b              = b,
    iters          = iters,
    iters_ranef    = iters_ranef,
    re_names       = re_names,
    group_levels   = group_levels,
    design         = design,
    collect_block1 = collect_block1
  )
  out$dispersion_ranef <- dispersion_ranef
  if (!is.null(iters_ranef_group)) {
    out$iters_sigma2_draws <- iters_ranef_group
  }
  out
}

#' Two-block Gibbs sweep with per-group ING Block~1 measurement dispersion
#' @noRd
.rGLMM_sweep_ing_block1 <- function(
    n_chains,
    start_fixef,
    inner_sweeps,
    design,
    block1_prior,
    ing_prior_list,
    pfamily_list,
    family,
    re_names,
    group_levels,
    collect_block1 = TRUE,
    progbar        = FALSE,
    stage_label    = "",
    diag_sweeps    = FALSE,
    fixef_mode     = NULL,
    b_mode         = NULL,
    b_start        = NULL,
    ptypes         = NULL,
    tau2_start     = NULL,
    use_cpp_block2 = TRUE
) {
  if (is.null(ptypes)) {
    ptypes <- vapply(pfamily_list, function(pf) pf$pfamily, character(1))
    names(ptypes) <- re_names
  }

  if (is.null(tau2_start)) {
    tau2_start <- .two_block_tau2_start_from_pfamily(pfamily_list, re_names)
  } else {
    if (is.null(names(tau2_start)) || !setequal(names(tau2_start), re_names)) {
      stop("'tau2_start' must be a named vector with names(re_names).",
           call. = FALSE)
    }
    tau2_start <- as.numeric(tau2_start[re_names])
    names(tau2_start) <- re_names
  }
  if (is.null(b_start)) {
    if (is.null(b_mode)) {
      stop("'b_start' or 'b_mode' required for batch init.", call. = FALSE)
    }
    b_start <- b_mode
  }

  batch <- .rGLMM_sweep_initialize(
    n_chains     = n_chains,
    start_fixef  = start_fixef,
    b_start      = b_start,
    tau2_start   = tau2_start,
    re_names     = re_names,
    group_levels = group_levels
  )

  progbar_use <- isTRUE(progbar)
  sweep_stats <- vector("list", inner_sweeps)
  dispersion_ranef <- numeric(n_chains)
  group_mode <- !is.null(ing_prior_list$shape_group)
  iters_ranef_group <- if (group_mode) {
    matrix(
      0,
      nrow = n_chains,
      ncol = length(group_levels),
      dimnames = list(NULL, group_levels)
    )
  } else {
    NULL
  }

  for (m in seq_len(inner_sweeps)) {
    prefix_b1 <- if (progbar_use) {
      .two_block_progbar_prefix(stage_label, m, inner_sweeps, "Block1")
    } else {
      ""
    }
    prefix_b2 <- if (progbar_use) {
      .two_block_progbar_prefix(stage_label, m, inner_sweeps, "Block2")
    } else {
      ""
    }

    b1 <- .two_block_block1_ing_all_chains(
      n                      = batch$n,
      fixef                  = batch$fixef,
      tau2                   = batch$tau2,
      b                      = batch$b,
      iters_ranef            = batch$iters_ranef,
      re_names               = re_names,
      group_levels           = group_levels,
      design                 = design,
      block1_prior           = block1_prior,
      ing_prior_list         = ing_prior_list,
      family                 = family,
      ptypes                 = ptypes,
      iters_ranef_group      = iters_ranef_group,
      progbar                = progbar_use,
      progbar_prefix         = prefix_b1,
      progbar_finish_newline = FALSE
    )

    batch$b           <- b1$b
    batch$iters_ranef <- b1$iters_ranef
    dispersion_ranef  <- b1$dispersion_ranef
    if (group_mode) {
      iters_ranef_group <- b1$iters_ranef_group
    }

    b2 <- .two_block_block2_all_chains(
      n                      = batch$n,
      b                      = batch$b,
      fixef                  = batch$fixef,
      tau2                   = batch$tau2,
      iters                  = batch$iters,
      re_names               = re_names,
      group_levels           = group_levels,
      design                 = design,
      pfamily_list           = pfamily_list,
      ptypes                 = ptypes,
      use_cpp_block2         = use_cpp_block2,
      progbar                = progbar_use,
      progbar_prefix         = prefix_b2,
      progbar_finish_newline = (m == inner_sweeps)
    )
    batch$fixef <- b2$fixef
    batch$tau2  <- b2$tau2
    batch$iters <- b2$iters

    sweep_stats[[m]] <- .two_block_snapshot_fixef_stats(
      fixef    = batch$fixef,
      re_names = re_names
    )
    if (progbar_use && n_chains <= 1L) {
      prefix_sweep <- if (nzchar(stage_label)) {
        sprintf("[%s] sweep %d/%d: ", stage_label, m, inner_sweeps)
      } else {
        sprintf("sweep %d/%d: ", m, inner_sweeps)
      }
      .two_block_progress_bar(m, inner_sweeps, prefix = prefix_sweep)
      .two_block_progress_bar_finish(newline = (m == inner_sweeps))
    }
  }

  out <- .rGLMM_sweep_save_ing(
    n              = batch$n,
    fixef          = batch$fixef,
    tau2           = batch$tau2,
    b              = batch$b,
    iters          = batch$iters,
    iters_ranef    = batch$iters_ranef,
    iters_ranef_group = iters_ranef_group,
    dispersion_ranef = dispersion_ranef,
    re_names       = re_names,
    group_levels   = group_levels,
    design         = design,
    collect_block1 = collect_block1
  )
  out$sweep_history <- .two_block_build_sweep_history(
    stage_label = stage_label,
    sweep_stats = sweep_stats,
    fixef_mode  = fixef_mode,
    re_names    = re_names
  )
  if (isTRUE(diag_sweeps)) {
    print(out$sweep_history)
  }
  invisible(out)
}

#' Format ING sweep output for staged \code{fixef.*} naming
#' @noRd
.rLMM_format_ing_sweep_out <- function(
    sweep_out,
    n,
    re_names,
    group_levels,
    fixef_mode,
    fixef_init
) {
  staged <- .rLMM_format_v2_out(
    v2_out       = sweep_out,
    n            = n,
    re_names     = re_names,
    group_levels = group_levels,
    fixef_mode   = fixef_mode,
    fixef_init   = fixef_init
  )
  staged$dispersion_ranef      <- sweep_out$dispersion_ranef
  staged$dispersion_ranef.mean <- if (is.matrix(sweep_out$dispersion_ranef)) {
    colMeans(sweep_out$dispersion_ranef)
  } else {
    mean(sweep_out$dispersion_ranef)
  }
  if (!is.null(sweep_out$iters_sigma2_draws)) {
    staged$sigma2.iters <- sweep_out$iters_sigma2_draws
  }
  staged$sweep_history         <- sweep_out$sweep_history
  staged
}

#' ING measurement LMM: pilot then main via ING sweep-outer engine
#' @noRd
.rLMMIngNormal_reg_run_with_pilot <- function(
    inp,
    block,
    P,
    ing_prior_list,
    pfamily_list,
    pf_summary,
    icm_tol,
    icm_maxit,
    progbar,
    verbose,
    stage_verbose = FALSE,
    gap_tol       = 0.0196,
    mode_gap_max  = 1.0,
    diag_sweeps   = FALSE,
    any_non_normal = TRUE,
    random_measurement = TRUE,
    engine_label  = "rLMMindepNormalGamma_reg_estimated_vcov",
    result_class  = "rLMMindepNormalGamma_reg_estimated_vcov",
    cl
) {
  re_names         <- inp$re_names
  group_levels     <- inp$group_levels
  group_name       <- inp$group_name
  tv_tol           <- inp$tv_tol
  n                <- inp$n
  n_pilot_arg      <- NULL
  m_convergence_pilot <- NULL
  rate_calibration <- NULL
  collect_block1   <- TRUE
  family           <- gaussian()
  is_gaussian      <- TRUE
  ptypes           <- pf_summary$ptypes

  dispersion_fix <- .rLMM_dispersion_fix_from_prior_list(
    ing_prior_list, fn_name = engine_label
  )

  rate_inputs <- .rLMM_measurement_rate_inputs(
    ing_prior_list, block, group_levels, engine_label
  )
  prior_list_block1_icm <- .rLMM_block1_prior_gaussian(P, dispersion_fix)
  prior_list_block1_rate <- .rLMM_block1_prior_gaussian(P, rate_inputs$dispersion_scalar)
  .two_block_validate_block1_prior(prior_list_block1_icm, family = family)
  .two_block_validate_block1_prior(prior_list_block1_rate, family = family)

  calib_meta <- .rLMM_rate_calibration_meta(
    any_non_normal     = any_non_normal,
    random_measurement = random_measurement
  )

  gap_tol <- .two_block_validate_gap_tol(gap_tol)

  will_pilot <- .two_block_pilot_will_run(
    is_gaussian,
    n_pilot_arg,
    gap_tol,
    tv_tol,
    any_non_normal     = any_non_normal,
    random_measurement = random_measurement
  )
  run_pilot <- will_pilot
  run_ub    <- will_pilot && !is.null(tv_tol)

  if (run_pilot && is.null(m_convergence_pilot)) {
    m_convergence_pilot <- if (!is.null(tv_tol)) {
      NULL
    } else {
      10L
    }
  }

  icm <- .rLMM_icm_at_start(
    y                     = inp$y,
    x                     = inp$x,
    block                 = block,
    x_hyper               = inp$x_hyper,
    prior_list_block1     = prior_list_block1_icm,
    pfamily_list          = pfamily_list,
    re_names              = re_names,
    group_levels          = group_levels,
    group_name            = group_name,
    icm_tol               = icm_tol,
    icm_maxit             = icm_maxit,
    verbose               = verbose,
    engine_label          = engine_label
  )
  fixef_mode <- icm$start
  ranef_mode <- icm$b_start
  icm_info   <- icm$icm
  if (isTRUE(verbose)) {
    cat(sprintf(
      "  %s: Block 2 start at lmer tau^2 plug-in (converged: %s, %d iter, delta = %.2e)\n\n",
      engine_label,
      icm_info$converged,
      icm_info$iterations,
      icm_info$delta
    ))
  }
  b_start <- ranef_mode

  design <- list(
    y             = inp$y,
    Z             = inp$x,
    groups        = factor(block, levels = group_levels),
    X_hyper       = inp$x_hyper,
    re_coef_names = re_names,
    group_name    = group_name,
    re_rank       = .lmebayes_re_rank_from_Z(
      inp$x, block, group_levels = group_levels
    )
  )

  if (is.null(b_start)) {
    b_start <- matrix(
      0,
      nrow = length(group_levels),
      ncol = length(re_names),
      dimnames = list(group_levels, re_names)
    )
  }

  fixef_mode_ref <- fixef_mode
  b_mode_ref     <- b_start
  progbar_use    <- isTRUE(progbar) || isTRUE(verbose) || isTRUE(stage_verbose)

  rate <- two_block_rate_from_pfamily_list(
    x                 = inp$x,
    block             = block,
    x_hyper           = inp$x_hyper,
    prior_list_block1 = prior_list_block1_rate,
    pfamily_list      = pfamily_list,
    weights           = rate_inputs$weights,
    family            = family,
    group_levels      = group_levels
  )

  m_min <- NULL
  if (!is.null(tv_tol)) {
    m_min <- .two_block_cap_inner_sweeps(
    two_block_l_for_tv(rate, tv_tol, method = "theorem3") + 1L
  )
  }

  p_dim            <- sum(vapply(fixef_mode, length, integer(1L)))
  D_max            <- if (!is.null(mode_gap_max)) sqrt(p_dim) * mode_gap_max else 0
  m_pilot_from_gap <- NULL

  if (run_pilot && is.null(m_convergence_pilot) && !is.null(tv_tol)) {
    erf1_inv_tv <- stats::qnorm((tv_tol + 1) / 2) / sqrt(2)
    c_tol       <- erf1_inv_tv * 2 * sqrt(2)
    m_pilot_from_gap <- .two_block_m_pilot_from_gap(rate, D_max, c_tol, m_min)
    m_convergence_pilot <- m_pilot_from_gap
  }

  pilot_plan <- .two_block_resolve_pilot_plan(
    is_gaussian         = is_gaussian,
    n                   = n,
    n_pilot_arg         = n_pilot_arg,
    gap_tol             = gap_tol,
    tv_tol              = tv_tol,
    m_convergence_user  = NULL,
    m_convergence_pilot = m_convergence_pilot,
    rate                = rate,
    p_dim               = p_dim,
    m_min               = m_min,
    any_non_normal      = any_non_normal,
    random_measurement  = random_measurement
  )
  n_pilot        <- pilot_plan$n_pilot
  m_convergence  <- pilot_plan$m_convergence
  pilot_cost_opt <- pilot_plan$pilot_cost_opt
  run_pilot      <- n_pilot > 0L
  run_ub         <- run_pilot && !is.null(tv_tol)

  if (is.null(m_min) && !run_pilot) {
    m_convergence <- 10L
  }

  if (is.null(rate_calibration) && !is.null(tv_tol)) {
    rate_calibration <- list(
      lambda_star = rate$lambda_star,
      eigenvalues = rate$eigenvalues,
      m_min       = m_min
    )
  }

  calib_label <- calib_meta$label

  if (isTRUE(verbose) && !is.null(tv_tol)) {
    cat(sprintf(
      paste0(
        "--- %s: convergence calibration [%s]:\n",
        "    lambda* = %.4f, tv_tol = %g => m_min = %d (mode start), ",
        "main m_convergence = %d ---\n\n"
      ),
      engine_label, calib_label, rate$lambda_star, tv_tol, m_min, m_convergence
    ))
    if (run_pilot && !is.null(mode_gap_max) && !is.null(m_pilot_from_gap)) {
      cat(sprintf(
        paste0(
          "--- %s: pilot sweep calibration [mode_gap_max = %g SD/dim, p = %d, ",
          "D_max = %.4f]:\n    m_min = %d, lambda* = %.4f => ",
          "m_convergence_pilot = %d ---\n\n"
        ),
        engine_label, mode_gap_max, p_dim, D_max, m_min,
        rate$lambda_star, m_convergence_pilot
      ))
    }
    if (run_pilot) {
      .two_block_print_pilot_plan(
        pilot_plan          = pilot_plan,
        n                   = n,
        m_convergence_pilot = m_convergence_pilot,
        rate                = rate,
        tv_tol              = tv_tol,
        p                   = p_dim,
        verbose             = verbose
      )
    }
  }

  convergence_info <- list(
    method              = calib_meta$method,
    tv_tol              = tv_tol,
    gap_tol             = gap_tol,
    n_pilot             = n_pilot,
    n_pilot_source      = pilot_plan$n_pilot_source,
    n_pilot_gap_tol     = pilot_plan$n_pilot_gap_tol,
    lambda_star         = rate$lambda_star,
    eigenvalues         = rate$eigenvalues,
    m_min               = m_min,
    m_certificate       = pilot_plan$m_certificate,
    m_convergence       = m_convergence,
    m_convergence_pilot = if (run_pilot) m_convergence_pilot else NULL,
    mode_gap_max        = if (run_pilot) mode_gap_max else NULL,
    m_pilot_from_gap    = if (run_pilot) m_pilot_from_gap else NULL,
    pilot_cost_opt      = pilot_cost_opt,
    draw_engine         = "rGLMM_sweep_ing_block1_ind"
  )

  m_convergence_used <- m_convergence
  fixef_init         <- fixef_mode
  pilot_res          <- NULL
  pilot_chisq        <- NULL
  pilot_ub           <- NULL
  tau2_start_main    <- if (!is.null(icm) && !is.null(icm$tau2_start)) {
    icm$tau2_start
  } else {
    .two_block_tau2_start_from_pfamily(pfamily_list, re_names)
  }

  sweep_common <- list(
    design         = design,
    block1_prior   = prior_list_block1_icm,
    ing_prior_list = ing_prior_list,
    pfamily_list   = pfamily_list,
    family         = family,
    re_names       = re_names,
    group_levels   = group_levels,
    collect_block1 = collect_block1,
    progbar        = progbar_use,
    fixef_mode     = fixef_mode_ref,
    b_mode         = b_mode_ref,
    b_start        = b_mode_ref,
    ptypes         = ptypes,
    use_cpp_block2 = TRUE
  )

  if (run_pilot) {
    if (isTRUE(verbose)) {
      cat(sprintf(
        "--- %s [sweep-outer]: pilot stage (%d chains; m_convergence_pilot = %d) ---\n\n",
        engine_label, n_pilot, m_convergence_pilot
      ))
    }

    pilot_raw <- do.call(
      .rGLMM_sweep_ing_block1,
      c(
        list(
          n_chains     = n_pilot,
          start_fixef  = fixef_mode,
          inner_sweeps = m_convergence_pilot,
          stage_label  = "pilot",
          diag_sweeps  = isTRUE(diag_sweeps),
          tau2_start   = tau2_start_main
        ),
        sweep_common
      )
    )

    pilot_chisq <- .two_block_pilot_chisq_test(
      fixef_draws = pilot_raw$fixef_draws,
      re_names    = re_names,
      fixef_mode  = fixef_mode,
      n_pilot     = n_pilot
    )

    if (isTRUE(stage_verbose) || isTRUE(verbose)) {
      cat(sprintf(
        "--- %s: pilot vs mode chi-squared test: p = %.4g (df = %d, n_pilot = %d) ---\n\n",
        engine_label,
        pilot_chisq$p_value, pilot_chisq$df, pilot_chisq$n_pilot
      ))
    }

    fixef_init <- .two_block_fixef_colmeans(
      pilot_raw$fixef_draws, re_names, fixef_mode
    )

    tau2_start_main <- .two_block_tau2_start_from_dispersion_draws(
      pilot_raw$dispersion_fixef_draws, re_names
    )

    if (run_ub) {
      pilot_ub <- .two_block_pilot_ub_from_coefficients(
        pilot_coefficients = pilot_raw$coefficients,
        n_pilot            = n_pilot,
        re_names           = re_names,
        group_levels       = group_levels,
        group_name         = group_name,
        x                  = inp$x,
        block              = block,
        x_hyper            = inp$x_hyper,
        prior_list         = prior_list_block1_rate,
        pfamily_list       = pfamily_list,
        family             = family,
        tv_tol             = tv_tol,
        dispersion         = rate_inputs$dispersion_scalar,
        weights            = rate_inputs$weights
      )
      if (pilot_ub$m_min_upper > m_convergence_used) {
        m_convergence_used <- pilot_ub$m_min_upper
      }
      convergence_info$lambda_star_upper <- pilot_ub$rate_upper$lambda_star
      convergence_info$eigenvalues_upper <- pilot_ub$max_eigenvalues
      convergence_info$m_min_upper       <- pilot_ub$m_min_upper
      convergence_info$i_max_rate        <- pilot_ub$i_max_rate
      convergence_info$lambda_star_vec   <- pilot_ub$lambda_star_vec
      convergence_info$m_convergence     <- m_convergence_used
    }

    if (isTRUE(stage_verbose) && run_ub) {
      .two_block_print_pilot_stage_diagnostics(
        n_pilot            = n_pilot,
        n_main             = n,
        pilot_ub           = pilot_ub,
        rate_calibration   = rate_calibration,
        m_convergence_used = m_convergence_used
      )
    } else if (isTRUE(verbose)) {
      cat(sprintf(
        "--- %s [sweep-outer]: pilot complete; main stage (%d chains; m_convergence = %d) ---\n\n",
        engine_label, n, m_convergence_used
      ))
    }

    pilot_res <- .rLMM_format_ing_sweep_out(
      sweep_out    = pilot_raw,
      n            = n_pilot,
      re_names     = re_names,
      group_levels = group_levels,
      fixef_mode   = fixef_mode,
      fixef_init   = fixef_mode
    )
  } else if (isTRUE(verbose)) {
    cat(sprintf(
      "--- %s [sweep-outer]: main stage (%d chains; m_convergence = %d) ---\n\n",
      engine_label, n, m_convergence_used
    ))
  }

  main_raw <- do.call(
    .rGLMM_sweep_ing_block1,
    c(
      list(
        n_chains     = n,
        start_fixef  = fixef_init,
        inner_sweeps = m_convergence_used,
        stage_label  = "main",
        diag_sweeps  = isTRUE(diag_sweeps),
        tau2_start   = tau2_start_main
      ),
      sweep_common
    )
  )

  draw_engine_args <- c(
    list(
      n_chains     = n,
      start_fixef  = fixef_init,
      inner_sweeps = m_convergence_used,
      stage_label  = "main",
      diag_sweeps  = isTRUE(diag_sweeps),
      tau2_start   = tau2_start_main
    ),
    sweep_common
  )

  main_res <- .rLMM_format_ing_sweep_out(
    sweep_out    = main_raw,
    n            = n,
    re_names     = re_names,
    group_levels = group_levels,
    fixef_mode   = fixef_mode,
    fixef_init   = fixef_init
  )

  main_res$call                <- cl
  main_res$n_pilot             <- n_pilot
  main_res$gap_tol             <- gap_tol
  main_res$m_convergence       <- m_convergence_used
  main_res$m_convergence_pilot <- if (run_pilot) m_convergence_pilot else NULL
  main_res$convergence_info    <- convergence_info
  main_res$draw_engine         <- "rGLMM_sweep_ing_block1_ind"
  main_res$draw_engine_call    <- quote(.rGLMM_sweep_ing_block1)
  main_res$draw_engine_args    <- draw_engine_args
  main_res$pfamily_list        <- pfamily_list
  main_res$family              <- family
  main_res$prior_list          <- ing_prior_list
  main_res$dispersion_fix      <- dispersion_fix
  main_res$ing_prior_list      <- ing_prior_list
  main_res$ranef.mode          <- ranef_mode
  main_res$icm_info            <- icm_info
  main_res$ptypes              <- pf_summary$ptypes
  main_res$any_non_normal      <- any_non_normal

  if (run_pilot) {
    main_res$pilot       <- pilot_res
    main_res$pilot_chisq <- pilot_chisq
  }
  if (run_ub) {
    main_res$pilot_ub <- pilot_ub
    main_res$tv_tol   <- tv_tol
  }

  class(main_res) <- c(result_class, "rLMMindepNormalGamma_reg", "list")
  main_res
}

#' Shared sampling pipeline for \code{rLMMNormal_reg_*_vcov} routes
#' @noRd
.rLMMNormal_reg_run <- function(
    inp,
    block,
    P,
    dispersion,
    pfamily_list,
    pf_summary,
    icm_tol,
    icm_maxit,
    progbar,
    verbose,
    engine_label,
    any_non_normal,
    draw_engine,
    result_class,
    cl,
    fixef_start = NULL
) {
  prior_list_block1 <- list(
    P          = P,
    dispersion = dispersion,
    ddef       = FALSE
  )
  .two_block_validate_block1_prior(prior_list_block1, family = gaussian())

  icm_info   <- NULL
  ranef_mode <- NULL

  if (is.null(fixef_start)) {
    icm <- .rLMM_icm_at_start(
      y                 = inp$y,
      x                 = inp$x,
      block             = block,
      x_hyper           = inp$x_hyper,
      prior_list_block1 = prior_list_block1,
      pfamily_list      = pfamily_list,
      re_names          = inp$re_names,
      group_levels      = inp$group_levels,
      group_name        = inp$group_name,
      icm_tol           = icm_tol,
      icm_maxit         = icm_maxit,
      verbose           = verbose,
      engine_label      = engine_label
    )
    fixef_mode <- icm$start
    ranef_mode <- icm$b_start
    icm_info   <- icm$icm
  } else {
    if (!is.list(fixef_start) || is.null(names(fixef_start))) {
      stop("'fixef_start' must be a named list.", call. = FALSE)
    }
    if (!setequal(names(fixef_start), inp$re_names)) {
      stop("names(fixef_start) must match re_coef_names.", call. = FALSE)
    }
    fixef_mode <- fixef_start[inp$re_names]
  }

  calib <- .rLMM_calibrate_m_convergence(
    x                 = inp$x,
    block             = block,
    x_hyper           = inp$x_hyper,
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
    group_levels      = inp$group_levels,
    tv_tol            = inp$tv_tol,
    any_non_normal    = any_non_normal,
    engine_label      = engine_label,
    verbose           = verbose
  )
  m_convergence    <- calib$m_convergence
  convergence_info <- calib$convergence_info
  convergence_info$draw_engine <- draw_engine

  out <- two_block_rNormal_reg(
    n                 = inp$n,
    y                 = inp$y,
    x                 = inp$x,
    block             = block,
    x_hyper           = inp$x_hyper,
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
    fixef_start       = fixef_mode,
    re_coef_names     = inp$re_names,
    group_levels      = inp$group_levels,
    group_name        = inp$group_name,
    family            = gaussian(),
    m_convergence     = m_convergence,
    progbar           = progbar
  )

  staged <- .rLMM_format_v2_out(
    v2_out       = out,
    n            = inp$n,
    re_names     = inp$re_names,
    group_levels = inp$group_levels,
    fixef_mode   = fixef_mode,
    fixef_init   = fixef_mode
  )

  staged$call             <- cl
  staged$m_convergence    <- m_convergence
  staged$convergence_info <- convergence_info
  staged$draw_engine      <- draw_engine
  staged$pfamily_list     <- pfamily_list
  staged$prior_list       <- prior_list_block1
  staged$family           <- gaussian()
  staged$ranef.mode       <- ranef_mode
  staged$icm_info         <- icm_info
  staged$ptypes           <- pf_summary$ptypes
  staged$any_non_normal   <- any_non_normal

  class(staged) <- c(result_class, "rLMMNormal_reg", "list")
  staged
}

#' @describeIn rLMM_reg Dispatcher for fixed \eqn{\sigma^2}: routes to
#'   \code{\link{rLMMNormal_reg_known_vcov}} or
#'   \code{\link{rLMMNormal_reg_estimated_vcov}} by Block~2 \code{pfamily_list}.
#' @export
rLMMNormal_reg <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    P,
    prior_list,
    pfamily_list,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    re_coef_names   = colnames(x),
    group_levels    = levels(block),
    group_name      = NULL,
    progbar         = TRUE,
    verbose         = FALSE
) {
  cl <- match.call()

  inp <- .rLMM_validate_matrix_inputs(
    n, y, x, x_hyper, tv_tol,
    re_coef_names, group_levels, group_name, block
  )
  P <- .rLMM_validate_P(P, length(inp$re_names))
  dispersion <- .rLMM_validate_fixed_dispersion_prior_list(prior_list)
  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)

  route_fn <- if (pf_summary$all_dNormal) {
    rLMMNormal_reg_known_vcov
  } else {
    rLMMNormal_reg_estimated_vcov
  }
  mc <- match.call(expand.dots = FALSE)
  mc[[1L]] <- route_fn
  out <- eval(mc, parent.frame())
  out$call <- cl
  out
}

#' @describeIn rLMM_reg Fixed \eqn{\sigma^2}; all Block~2 \code{dNormal} (known
#'   \eqn{\tau^2_k}). Exact Theorem~3 rate calibration.
#' @export
rLMMNormal_reg_known_vcov <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    P,
    prior_list,
    pfamily_list,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    re_coef_names   = colnames(x),
    group_levels    = levels(block),
    group_name      = NULL,
    progbar         = TRUE,
    verbose         = FALSE
) {
  cl <- match.call()
  fn_name <- "rLMMNormal_reg_known_vcov"

  inp <- .rLMM_validate_matrix_inputs(
    n, y, x, x_hyper, tv_tol,
    re_coef_names, group_levels, group_name, block
  )
  P <- .rLMM_validate_P(P, length(inp$re_names), fn_name = fn_name)
  dispersion <- .rLMM_validate_fixed_dispersion_prior_list(
    prior_list, fn_name = fn_name
  )
  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (!pf_summary$all_dNormal) {
    stop(
      fn_name, "(): all Block~2 components must be dNormal(); ",
      "use rLMMNormal_reg_estimated_vcov() or rLMMNormal_reg().",
      call. = FALSE
    )
  }

  .rLMMNormal_reg_run(
    inp              = inp,
    block            = block,
    P                = P,
    dispersion       = dispersion,
    pfamily_list     = pfamily_list,
    pf_summary       = pf_summary,
    icm_tol          = icm_tol,
    icm_maxit        = icm_maxit,
    progbar          = progbar,
    verbose          = verbose,
    engine_label     = fn_name,
    any_non_normal   = FALSE,
    draw_engine      = "two_block_rNormal_reg_known_vcov",
    result_class     = "rLMMNormal_reg_known_vcov",
    cl               = cl
  )
}

#' @describeIn rLMM_reg Fixed \eqn{\sigma^2}; ING Block~2 (estimated \eqn{\tau^2_k}).
#'   Optional pilot stage; conservative \code{disp_lower} rate bound.
#' @export
rLMMNormal_reg_estimated_vcov <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    P,
    prior_list,
    pfamily_list,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    re_coef_names   = colnames(x),
    group_levels    = levels(block),
    group_name      = NULL,
    progbar         = TRUE,
    verbose         = FALSE,
    gap_tol         = 0.0196,
    mode_gap_max    = 1.0,
    diag_sweeps     = FALSE,
    stage_verbose   = FALSE
) {
  cl <- match.call()
  fn_name <- "rLMMNormal_reg_estimated_vcov"

  inp <- .rLMM_validate_matrix_inputs(
    n, y, x, x_hyper, tv_tol,
    re_coef_names, group_levels, group_name, block
  )
  P <- .rLMM_validate_P(P, length(inp$re_names), fn_name = fn_name)
  dispersion <- .rLMM_validate_fixed_dispersion_prior_list(
    prior_list, fn_name = fn_name
  )
  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (pf_summary$all_dNormal) {
    stop(
      fn_name, "(): at least one Block~2 component must not be dNormal(); ",
      "use rLMMNormal_reg_known_vcov() or rLMMNormal_reg().",
      call. = FALSE
    )
  }

  .rLMMNormal_reg_run_with_pilot(
    inp            = inp,
    block          = block,
    P              = P,
    dispersion     = dispersion,
    pfamily_list   = pfamily_list,
    pf_summary     = pf_summary,
    icm_tol        = icm_tol,
    icm_maxit      = icm_maxit,
    progbar        = progbar,
    verbose        = verbose,
    stage_verbose  = stage_verbose,
    gap_tol        = gap_tol,
    mode_gap_max   = mode_gap_max,
    diag_sweeps    = diag_sweeps,
    engine_label   = fn_name,
    result_class   = "rLMMNormal_reg_estimated_vcov",
    cl             = cl
  )
}

#' @describeIn rLMM_reg Legacy outer-loop engine: draws \eqn{\sigma^2} via
#'   \code{\link{rGamma_reg}}, then calls \code{\link{rLMMNormal_reg}} each
#'   replicate. Prefer \code{\link{rLMMindepNormalGamma_reg_known_vcov}} or
#'   \code{\link{rLMMindepNormalGamma_reg_estimated_vcov}} for the ING Block~1
#'   sweep engine (formula-level entry is in **lmebayes**).
#' @export
rLMMindepNormalGamma_reg <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    P,
    prior_list,
    pfamily_list,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    re_coef_names   = colnames(x),
    group_levels    = levels(block),
    group_name      = NULL,
    progbar         = TRUE,
    verbose         = FALSE
) {
  cl <- match.call()
  inp <- .rLMM_validate_matrix_inputs(
    n, y, x, x_hyper, tv_tol,
    re_coef_names, group_levels, group_name, block
  )
  P <- .rLMM_validate_P(P, length(inp$re_names))
  prior_list <- .rLMM_validate_dGamma_dispersion_prior_list(prior_list)
  dispersion_fix <- .rLMM_dispersion_fix_from_prior_list(
    prior_list, fn_name = "rLMMindepNormalGamma_reg"
  )

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)

  prior_list_block1_cal <- list(
    P          = P,
    dispersion = dispersion_fix,
    ddef       = FALSE
  )
  .two_block_validate_block1_prior(prior_list_block1_cal, family = gaussian())

  icm <- .rLMM_icm_at_start(
    y                     = inp$y,
    x                     = inp$x,
    block                 = block,
    x_hyper               = inp$x_hyper,
    prior_list_block1     = prior_list_block1_cal,
    pfamily_list          = pfamily_list,
    re_names              = inp$re_names,
    group_levels          = inp$group_levels,
    group_name            = inp$group_name,
    icm_tol               = icm_tol,
    icm_maxit             = icm_maxit,
    verbose               = verbose,
    engine_label          = "rLMMindepNormalGamma_reg"
  )
  fixef_start <- icm$start
  ranef_mode  <- icm$b_start
  icm_info    <- icm$icm

  calib <- .rLMM_calibrate_m_convergence(
    x                 = inp$x,
    block             = block,
    x_hyper           = inp$x_hyper,
    prior_list_block1 = prior_list_block1_cal,
    pfamily_list      = pfamily_list,
    group_levels      = inp$group_levels,
    tv_tol            = inp$tv_tol,
    any_non_normal    = pf_summary$any_non_normal,
    engine_label      = "rLMMindepNormalGamma_reg",
    verbose           = verbose
  )
  m_convergence    <- calib$m_convergence
  convergence_info <- calib$convergence_info
  convergence_info$draw_engine <- "rLMMindepNormalGamma_reg_outer"

  if (is.null(ranef_mode)) {
    ranef_mode <- matrix(
      0,
      nrow = length(inp$group_levels),
      ncol = length(inp$re_names),
      dimnames = list(inp$group_levels, inp$re_names)
    )
  }
  b_mat <- ranef_mode

  n_obs <- length(inp$y)
  x_disp <- matrix(1, nrow = n_obs, ncol = 1L)
  wt     <- rep(1, n_obs)

  fixef_cur <- fixef_start
  re_names  <- inp$re_names
  n         <- inp$n

  fixef_draws <- stats::setNames(
    lapply(re_names, function(k) {
      q_k <- length(fixef_cur[[k]])
      matrix(NA_real_, nrow = n, ncol = q_k,
             dimnames = list(NULL, names(fixef_cur[[k]])))
    }),
    re_names
  )
  dispersion_ranef_draws <- numeric(n)
  coef_parts             <- vector("list", n)
  disp_fixef_draws       <- NULL
  iters_fixef_draws      <- NULL

  if (isTRUE(progbar) && n > 1L) {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_len(n)) {
    mu_hat <- .rLMM_observation_mu(
      inp$x, block, b_mat, inp$group_levels
    )
    gamma_out <- rGamma_reg(
      n           = 1L,
      y           = inp$y,
      x           = x_disp,
      prior_list  = prior_list,
      offset      = mu_hat,
      weights     = wt,
      family      = gaussian(),
      progbar     = FALSE,
      verbose     = FALSE
    )
    sigma2_i <- as.numeric(gamma_out$dispersion[1L])
    dispersion_ranef_draws[i] <- sigma2_i

    draw_engine <- if (pf_summary$all_dNormal) {
      "two_block_rNormal_reg_known_vcov"
    } else {
      "two_block_rNormal_reg_estimated_vcov"
    }
    inp_i <- inp
    inp_i$n <- 1L
    lmm_i <- .rLMMNormal_reg_run(
      inp            = inp_i,
      block          = block,
      P              = P,
      dispersion     = sigma2_i,
      pfamily_list   = pfamily_list,
      pf_summary     = pf_summary,
      icm_tol        = icm_tol,
      icm_maxit      = icm_maxit,
      progbar        = FALSE,
      verbose        = FALSE,
      engine_label   = "rLMMindepNormalGamma_reg",
      any_non_normal = pf_summary$any_non_normal,
      draw_engine    = draw_engine,
      result_class   = "rLMMNormal_reg",
      cl             = NULL,
      fixef_start    = fixef_cur
    )

    for (k in re_names) {
      fixef_draws[[k]][i, ] <- lmm_i$fixef[[k]][1L, ]
    }
    coef_parts[[i]] <- lmm_i$coefficients
    fixef_cur <- lapply(lmm_i$fixef, function(m) {
      stats::setNames(m[1L, ], colnames(m))
    })
    b_mat <- .rLMM_b_matrix_from_coefficients(
      lmm_i$coefficients,
      re_names,
      inp$group_levels,
      inp$group_name
    )
    if (i == 1L) {
      disp_fixef_draws  <- lmm_i$fixef.dispersion
      iters_fixef_draws <- lmm_i$fixef.iters
    } else {
      disp_fixef_draws  <- rbind(disp_fixef_draws, lmm_i$fixef.dispersion)
      iters_fixef_draws <- rbind(iters_fixef_draws, lmm_i$fixef.iters)
    }

    if (isTRUE(progbar) && n > 1L) {
      utils::setTxtProgressBar(pb, i / n)
    }
  }

  coefficients <- do.call(rbind, coef_parts)
  rownames(coefficients) <- NULL

  v2_out <- list(
    fixef_draws            = fixef_draws,
    coefficients           = coefficients,
    dispersion_fixef_draws = disp_fixef_draws,
    iters_fixef_draws      = iters_fixef_draws,
    mu_all_last            = build_mu_all(
      list(
        X_hyper       = inp$x_hyper,
        re_coef_names = re_names,
        groups        = factor(block, levels = inp$group_levels)
      ),
      fixef_cur,
      group_levels = inp$group_levels
    )$mu_all
  )

  staged <- .rLMM_format_v2_out(
    v2_out       = v2_out,
    n            = n,
    re_names     = re_names,
    group_levels = inp$group_levels,
    fixef_mode   = fixef_start,
    fixef_init   = fixef_start
  )

  staged$call                  <- cl
  staged$m_convergence         <- m_convergence
  staged$convergence_info      <- convergence_info
  staged$draw_engine           <- "rLMMindepNormalGamma_reg_outer"
  staged$pfamily_list          <- pfamily_list
  staged$prior_list            <- prior_list
  staged$dispersion_prior_list <- prior_list
  staged$dispersion_ranef      <- dispersion_ranef_draws
  staged$dispersion_ranef.mean <- mean(dispersion_ranef_draws)
  staged$family                <- gaussian()
  staged$ranef.mode            <- ranef_mode
  staged$icm_info              <- icm_info
  staged$P                     <- P
  staged$dispersion_fix          <- dispersion_fix
  staged$ptypes         <- pf_summary$ptypes
  staged$any_non_normal <- pf_summary$any_non_normal

  class(staged) <- c("rLMMindepNormalGamma_reg", "list")
  staged
}

#' @describeIn rLMM_reg Random \eqn{\sigma^2} (per-group ING Block~1); all
#'   Block~2 \code{dNormal}. Used when \code{dispersion_ranef} is dGamma and
#'   Block~2 is all \code{dNormal}.
#' @export
rLMMindepNormalGamma_reg_known_vcov <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    P,
    prior_list,
    pfamily_list,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    re_coef_names   = colnames(x),
    group_levels    = levels(block),
    group_name      = NULL,
    progbar         = TRUE,
    verbose         = FALSE
) {
  cl <- match.call()
  fn_name <- "rLMMindepNormalGamma_reg_known_vcov"

  inp <- .rLMM_validate_matrix_inputs(
    n, y, x, x_hyper, tv_tol,
    re_coef_names, group_levels, group_name, block
  )
  P <- .rLMM_validate_P(P, length(inp$re_names), fn_name = fn_name)
  ing_prior_list <- .rLMM_validate_ing_measurement_prior_list(
    prior_list, length(inp$re_names), fn_name = fn_name,
    group_levels = inp$group_levels
  )

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (!pf_summary$all_dNormal) {
    stop(
      fn_name, "(): all Block~2 components must be dNormal(); ",
      "use rLMMindepNormalGamma_reg_estimated_vcov().",
      call. = FALSE
    )
  }

  .rLMMIngNormal_reg_run_with_pilot(
    inp                = inp,
    block              = block,
    P                  = P,
    ing_prior_list     = ing_prior_list,
    pfamily_list       = pfamily_list,
    pf_summary         = pf_summary,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    progbar            = progbar,
    verbose            = verbose,
    any_non_normal     = FALSE,
    random_measurement = TRUE,
    engine_label       = fn_name,
    result_class       = "rLMMindepNormalGamma_reg_known_vcov",
    cl                 = cl
  )
}

#' @describeIn rLMM_reg Random \eqn{\sigma^2} (per-group ING Block~1); ING Block~2.
#'   Pilot/UB calibration via sweep-outer engine. Default path for dGamma
#'   \code{dispersion_ranef} with ING Block~2.
#' @export
rLMMindepNormalGamma_reg_estimated_vcov <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    P,
    prior_list,
    pfamily_list,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    re_coef_names   = colnames(x),
    group_levels    = levels(block),
    group_name      = NULL,
    progbar         = TRUE,
    verbose         = FALSE,
    gap_tol         = 0.0196,
    mode_gap_max    = 1.0,
    diag_sweeps     = FALSE,
    stage_verbose   = FALSE
) {
  cl <- match.call()
  fn_name <- "rLMMindepNormalGamma_reg_estimated_vcov"

  inp <- .rLMM_validate_matrix_inputs(
    n, y, x, x_hyper, tv_tol,
    re_coef_names, group_levels, group_name, block
  )
  P <- .rLMM_validate_P(P, length(inp$re_names), fn_name = fn_name)
  ing_prior_list <- .rLMM_validate_ing_measurement_prior_list(
    prior_list, length(inp$re_names), fn_name = fn_name,
    group_levels = inp$group_levels
  )

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (pf_summary$all_dNormal) {
    stop(
      fn_name, "(): at least one Block~2 component must not be dNormal(); ",
      "use rLMMindepNormalGamma_reg_known_vcov().",
      call. = FALSE
    )
  }

  .rLMMIngNormal_reg_run_with_pilot(
    inp                = inp,
    block              = block,
    P                  = P,
    ing_prior_list     = ing_prior_list,
    pfamily_list       = pfamily_list,
    pf_summary         = pf_summary,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    progbar            = progbar,
    verbose            = verbose,
    stage_verbose      = stage_verbose,
    gap_tol            = gap_tol,
    mode_gap_max       = mode_gap_max,
    diag_sweeps        = diag_sweeps,
    any_non_normal     = TRUE,
    random_measurement = TRUE,
    engine_label       = fn_name,
    result_class       = "rLMMindepNormalGamma_reg_estimated_vcov",
    cl                 = cl
  )
}
#' Format v2 batch output for staged \code{fixef.*} naming
#' @noRd
.rLMM_format_v2_out <- function(
    v2_out,
    n,
    re_names,
    group_levels,
    fixef_mode,
    fixef_init
) {
  x <- list(
    fixef_draws            = v2_out$fixef_draws,
    coefficients           = v2_out$coefficients,
    dispersion_fixef_draws = v2_out$dispersion_fixef_draws,
    iters_fixef_draws      = v2_out$iters_fixef_draws,
    iters_ranef_draws      = v2_out$iters_ranef_draws,
    mu_all_last            = v2_out$mu_all_last,
    re_coef_names          = re_names,
    group_levels           = group_levels,
    n                      = n
  )
  .two_block_as_staged_names(
    x,
    fixef_mode = fixef_mode,
    fixef_init = fixef_init
  )
}
