#' Matrix-level replicate-chain Gibbs engines for Bayesian GLMMs
#'
#' @description
#' Non-Gaussian (and optionally Gaussian) two-block GLMM samplers at matrix
#' level. Each stored draw runs \code{m_convergence} inner Gibbs sweeps via
#' \code{\link{rGLMM_sweep}}. Formula-level fitting lives in **lmebayes**
#' (\code{glmerb}/\code{rglmerb}); Gaussian models with observation dispersion
#' typically use the \code{\link{rLMM_reg}} routes.
#'
#' @section Two route engines:
#' Both routes run a \strong{pilot stage for non-Gaussian} families (local-Gaussian
#' rate calibration and chain-mean initialisation; skip only with
#' \code{n_pilot = 0L}). The route split is \strong{not} whether a pilot runs,
#' but how \strong{eigenvalue bounds} are set for Theorem~3 and post-pilot
#' upper-bound calibration:
#' \describe{
#'   \item{\code{rGLMM_reg_known_vcov}}{
#'     All Block~2 \code{dNormal} (known \eqn{\tau^2_k}): standard fixed-dispersion
#'     rate at the mode; post-pilot eigenvalue upper bound from pilot
#'     \code{coefficients} without ING \code{disp_lower} conservatism.}
#'   \item{\code{rGLMM_reg_estimated_vcov}}{
#'     At least one ING Block~2 component: conservative \code{disp_lower} plug-in
#'     in \code{\link{two_block_rate_from_pfamily_list}} (upper bound over truncated
#'     \eqn{\tau^2}); pilot updates \eqn{\tau^2} starts from dispersion draws.}
#' }
#'
#' @section Dispatcher:
#' \code{\link{rGLMM_reg}} validates Block~2 \code{pfamily_list} and delegates
#' to the appropriate route.
#'
#' @param n Number of stored main-stage draws. If \code{length(n) > 1}, the
#'   length is used.
#' @param y Response vector of length \code{nrow(x)}.
#' @param x Level-1 design matrix \code{Z} (\code{l2 x p_re}).
#' @param block Grouping factor or block partition of length \code{l2}.
#' @param x_hyper Named list of group-level design matrices (\code{J x q_k}),
#'   one per column of \code{x}.
#' @param prior_list Prior for Block~1: \code{P} or \code{Sigma},
#'   \code{dispersion} (required for \code{gaussian()}), optional \code{ddef}.
#' @param pfamily_list Named list of \code{pfamily} objects for Block~2.
#' @param icm_tol,icm_maxit ICM convergence controls for the internal Block~2 start.
#' @param offset,weights,family Passed to Block~1 (length \code{l2} or recycled).
#' @param re_coef_names Character vector naming columns of \code{x}.
#' @param group_levels Character vector defining row order of Block~1 draws.
#' @param group_name Name for the grouping column in \code{coefficients}.
#' @param gap_tol Legacy mode--mean gap for pilot chain count when \code{tv_tol}
#'   is \code{NULL}.
#' @param tv_tol Total-variation tolerance for Theorem~3 calibration.
#'   \code{n_pilot}, \code{m_convergence_pilot}, and main-stage
#'   \code{m_convergence} are derived internally.
#' @param mode_gap_max Pilot inner-sweep calibration when a pilot stage runs.
#' @param Gridtype,n_envopt,use_parallel,use_opencl Reserved (not yet forwarded).
#' @param verbose Print stage headers and convergence calibration lines.
#' @param progbar Progress bars during sampling.
#' @param stage_verbose Print pilot chi-squared and post-pilot UB diagnostics.
#' @param rate_calibration Optional rate object for \code{stage_verbose}.
#' @param collect_block1 Collect Block~1 \code{coefficients} from main chains.
#' @family simfuncs
#' @seealso \code{\link{rGLMM_sweep}}, \code{\link{rLMM_reg}}
#' @name rGLMM_reg
NULL

#' Shared matrix-level validation for GLMM replicate-chain engines
#' @noRd
.rGLMM_validate_matrix_inputs <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    tv_tol,
    re_coef_names,
    group_levels,
    group_name,
    family,
    mode_gap_max,
    gap_tol,
    prior_list,
    pfamily_list
) {
  family <- .two_block_normalize_family(family)
  is_gaussian <- identical(family$family, "gaussian")

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  gap_tol <- .two_block_validate_gap_tol(gap_tol)

  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be a single positive finite number.",
           call. = FALSE)
    }
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
      stop("'x' must have at least 1 column.", call. = FALSE)
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

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, re_names, J = length(group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)

  if (!is.null(tv_tol)) {
    if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
        !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
      stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
    }
  }

  .two_block_validate_block1_prior(prior_list, family = family)

  list(
    n              = n,
    y              = y,
    x              = x,
    block          = block,
    x_hyper        = x_hyper,
    re_names       = re_names,
    group_levels   = group_levels,
    group_name     = group_name,
    family         = family,
    is_gaussian    = is_gaussian,
    prior_list     = prior_list,
    pfamily_list   = pfamily_list,
    pf_summary     = pf_summary,
    ptypes         = pf_summary$ptypes,
    any_non_normal = pf_summary$any_non_normal,
    tv_tol         = tv_tol,
    gap_tol        = gap_tol
  )
}

#' Main GLMM sampling pipeline entry (ICM, optional pilot, main sweep)
#' @noRd
.rGLMM_reg_run <- function(
    inp,
    icm_tol,
    icm_maxit,
    mode_gap_max,
    verbose,
    progbar,
    stage_verbose,
    rate_calibration,
    collect_block1,
    engine_label,
    result_class,
    cl
) {
  .rGLMM_reg_run_pipeline(
    inp                = inp,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    mode_gap_max       = mode_gap_max,
    verbose            = verbose,
    progbar            = progbar,
    stage_verbose      = stage_verbose,
    rate_calibration   = rate_calibration,
    collect_block1     = collect_block1,
    engine_label       = engine_label,
    result_class       = result_class,
    cl                 = cl
  )
}

#' Full GLMM sampling pipeline (same as \code{.rGLMM_reg_run}; reserved for
#' future pilot-policy split)
#' @noRd
.rGLMM_reg_run_with_pilot <- function(
    inp,
    icm_tol,
    icm_maxit,
    mode_gap_max,
    verbose,
    progbar,
    stage_verbose,
    rate_calibration,
    collect_block1,
    engine_label,
    result_class,
    cl
) {
  .rGLMM_reg_run_pipeline(
    inp                = inp,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    mode_gap_max       = mode_gap_max,
    verbose            = verbose,
    progbar            = progbar,
    stage_verbose      = stage_verbose,
    rate_calibration   = rate_calibration,
    collect_block1     = collect_block1,
    engine_label       = engine_label,
    result_class       = result_class,
    cl                 = cl
  )
}

#' Shared GLMM replicate-chain pipeline (sweep-outer driver)
#' @noRd
.rGLMM_reg_run_pipeline <- function(
    inp,
    icm_tol,
    icm_maxit,
    mode_gap_max,
    verbose,
    progbar,
    stage_verbose,
    rate_calibration,
    collect_block1,
    engine_label,
    result_class,
    cl
) {
  n              <- inp$n
  y              <- inp$y
  x              <- inp$x
  block          <- inp$block
  x_hyper        <- inp$x_hyper
  re_names       <- inp$re_names
  group_levels   <- inp$group_levels
  group_name     <- inp$group_name
  family         <- inp$family
  is_gaussian    <- inp$is_gaussian
  prior_list     <- inp$prior_list
  pfamily_list   <- inp$pfamily_list
  pf_summary     <- inp$pf_summary
  ptypes         <- inp$ptypes
  any_non_normal <- inp$any_non_normal
  tv_tol         <- inp$tv_tol
  gap_tol        <- inp$gap_tol
  n_pilot_arg         <- NULL
  m_convergence_pilot <- NULL

  will_pilot <- .two_block_pilot_will_run(
    is_gaussian, n_pilot_arg, gap_tol, tv_tol,
    any_non_normal = any_non_normal
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

  icm_info <- NULL
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
    prior_list   = prior_list,
    pfamily_list = pfamily_list,
    re_names     = re_names,
    family       = family,
    tol          = icm_tol,
    maxit        = icm_maxit
  )
  fixef_mode <- icm$start
  b_start    <- icm$b_start
  ranef_mode <- b_start
  icm_info   <- icm$icm
  if (isTRUE(verbose)) {
    if (isTRUE(any_non_normal)) {
      icm_what <- "Block 2 start at lmer tau^2 plug-in"
    } else if (is_gaussian) {
      icm_what <- "ICM posterior mean"
    } else {
      icm_what <- "ICM posterior mode"
    }
    cat(sprintf(
      "  %s: %s (converged: %s, %d iter, delta = %.2e)\n\n",
      engine_label,
      icm_what,
      icm_info$converged,
      icm_info$iterations,
      icm_info$delta
    ))
  }

  design <- list(
    y             = y,
    Z             = x,
    groups        = factor(block, levels = group_levels),
    X_hyper       = x_hyper,
    re_coef_names = re_names,
    group_name    = group_name
  )

  fixef_mode_ref <- fixef_mode
  b_mode_ref     <- b_start
  progbar_use    <- isTRUE(progbar) || isTRUE(verbose) || isTRUE(stage_verbose)

  rate <- .rGLMM_rate_at_mode(
    design       = design,
    prior_list   = prior_list,
    pfamily_list = pfamily_list,
    family       = family,
    b_mode       = b_start,
    group_levels = group_levels,
    is_gaussian  = is_gaussian
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
    any_non_normal      = any_non_normal
  )
  n_pilot          <- pilot_plan$n_pilot
  m_convergence    <- pilot_plan$m_convergence
  pilot_cost_opt   <- pilot_plan$pilot_cost_opt
  run_pilot        <- n_pilot > 0L
  run_ub           <- run_pilot && !is.null(tv_tol)

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

  calib_label <- if (is_gaussian) {
    "exact (Gaussian posterior)"
  } else {
    sprintf("approximate (local-Gaussian at mode, %s)", family$family)
  }
  if (isTRUE(any_non_normal)) {
    calib_label <- paste0(
      calib_label,
      "; conservative: non-dNormal RE dispersion (disp_lower plug-in)"
    )
  }

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

  method_label <- if (is_gaussian) "exact" else "local_gaussian_mode"
  if (isTRUE(any_non_normal)) {
    method_label <- paste0(method_label, "+disp_lower_bound")
  }

  convergence_info <- list(
    method              = method_label,
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
    draw_engine         = "rGLMM_sweep"
  )

  m_convergence_used <- m_convergence
  fixef_init         <- fixef_mode
  pilot_res          <- NULL
  pilot_chisq        <- NULL
  pilot_ub           <- NULL
  tau2_start_main    <- .two_block_tau2_start_from_pfamily(pfamily_list, re_names)

  if (run_pilot) {
    if (isTRUE(verbose)) {
      cat(sprintf(
        "--- %s [sweep-outer]: pilot stage (%d chains; m_convergence_pilot = %d) ---\n\n",
        engine_label, n_pilot, m_convergence_pilot
      ))
    }

    pilot_raw <- rGLMM_sweep(
      n_chains       = n_pilot,
      start_fixef    = fixef_mode,
      inner_sweeps   = m_convergence_pilot,
      design         = design,
      block1_prior   = prior_list,
      pfamily_list   = pfamily_list,
      family         = family,
      re_names       = re_names,
      group_levels   = group_levels,
      collect_block1 = collect_block1,
      progbar        = progbar_use,
      stage_label    = "pilot",
      fixef_mode     = fixef_mode_ref,
      b_mode         = b_mode_ref,
      b_start        = b_mode_ref,
      ptypes         = ptypes
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

    if (isTRUE(any_non_normal)) {
      tau2_start_main <- .two_block_tau2_start_from_dispersion_draws(
        pilot_raw$dispersion_fixef_draws, re_names
      )
    }

    if (run_ub) {
      pilot_ub <- .two_block_pilot_ub_from_coefficients(
        pilot_coefficients = pilot_raw$coefficients,
        n_pilot            = n_pilot,
        re_names           = re_names,
        group_levels       = group_levels,
        group_name         = group_name,
        x                  = x,
        block              = block,
        x_hyper            = x_hyper,
        prior_list         = prior_list,
        pfamily_list       = pfamily_list,
        family             = family,
        tv_tol             = tv_tol
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

    pilot_res <- .rGLMM_format_v6_out(
      v6_out       = pilot_raw,
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

  main_raw <- rGLMM_sweep(
    n_chains       = n,
    start_fixef    = fixef_init,
    inner_sweeps   = m_convergence_used,
    design         = design,
    block1_prior   = prior_list,
    pfamily_list   = pfamily_list,
    family         = family,
    re_names       = re_names,
    group_levels   = group_levels,
    collect_block1 = collect_block1,
    progbar        = progbar_use,
    stage_label    = "main",
    fixef_mode     = fixef_mode_ref,
    b_mode         = b_mode_ref,
    b_start        = b_mode_ref,
    ptypes         = ptypes,
    tau2_start     = tau2_start_main
  )

  draw_engine_args <- list(
    n_chains       = n,
    start_fixef    = fixef_init,
    inner_sweeps   = m_convergence_used,
    design         = design,
    block1_prior   = prior_list,
    pfamily_list   = pfamily_list,
    family         = family,
    re_names       = re_names,
    group_levels   = group_levels,
    collect_block1 = collect_block1,
    progbar        = progbar_use,
    stage_label    = "main",
    fixef_mode     = fixef_mode_ref,
    b_mode         = b_mode_ref,
    b_start        = b_mode_ref,
    ptypes         = ptypes,
    tau2_start     = tau2_start_main
  )

  main_res <- .rGLMM_format_v6_out(
    v6_out       = main_raw,
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
  main_res$draw_engine         <- "rGLMM_sweep"
  main_res$draw_engine_call    <- quote(rGLMM_sweep)
  main_res$draw_engine_args    <- draw_engine_args
  main_res$pfamily_list        <- pfamily_list
  main_res$family              <- family
  main_res$prior_list          <- prior_list
  main_res$ranef.mode          <- ranef_mode
  main_res$icm_info            <- icm_info
  main_res$ptypes              <- pf_summary$ptypes
  main_res$any_non_normal      <- pf_summary$any_non_normal

  if (run_pilot) {
    main_res$pilot       <- pilot_res
    main_res$pilot_chisq <- pilot_chisq
  }
  if (run_ub) {
    main_res$pilot_ub <- pilot_ub
    main_res$tv_tol   <- tv_tol
  }

  class(main_res) <- c(result_class, "list")
  main_res
}

#' Local-Gaussian rate at the ICM mode
#' @noRd
.rGLMM_rate_at_mode <- function(
    design,
    prior_list,
    pfamily_list,
    family,
    b_mode,
    group_levels,
    is_gaussian
) {
  if (is_gaussian) {
    two_block_rate_from_pfamily_list(
      x                 = design$Z,
      block             = design$groups,
      x_hyper           = design$X_hyper,
      prior_list_block1 = prior_list,
      pfamily_list      = pfamily_list,
      family            = gaussian(),
      group_levels      = group_levels
    )
  } else {
    mode_w <- two_block_mode_weights(
      x            = design$Z,
      block        = design$groups,
      b_mode       = b_mode,
      family       = family,
      group_levels = group_levels
    )
    two_block_rate_from_pfamily_list(
      x                 = design$Z,
      block             = design$groups,
      x_hyper           = design$X_hyper,
      prior_list_block1 = prior_list,
      pfamily_list      = pfamily_list,
      weights           = mode_w$weights,
      family            = family,
      group_levels      = group_levels
    )
  }
}

#' Format v6 batch output for staged \code{fixef.*} naming
#' @noRd
.rGLMM_format_v6_out <- function(
    v6_out,
    n,
    re_names,
    group_levels,
    fixef_mode,
    fixef_init
) {
  x <- list(
    fixef_draws            = v6_out$fixef_draws,
    coefficients           = v6_out$coefficients,
    dispersion_fixef_draws = v6_out$dispersion_fixef_draws,
    iters_fixef_draws      = v6_out$iters_fixef_draws,
    iters_ranef_draws      = v6_out$iters_ranef_draws,
    mu_all_last            = v6_out$mu_all_last,
    sweep_history          = v6_out$sweep_history,
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

#' @describeIn rGLMM_reg All Block~2 \code{dNormal} (known \eqn{\tau^2_k}).
#'   Non-Gaussian: pilot always (unless \code{n_pilot = 0L}); standard
#'   fixed-dispersion eigenvalue rate bounds (no ING \code{disp_lower} path).
#' @export
rGLMM_reg_known_vcov <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    prior_list,
    pfamily_list,
    icm_tol             = 1e-10,
    icm_maxit           = 200L,
    offset              = NULL,
    weights             = 1,
    family              = gaussian(),
    re_coef_names       = colnames(x),
    group_levels        = levels(block),
    group_name          = NULL,
    gap_tol             = 0.0196,
    tv_tol              = 0.01,
    mode_gap_max        = 1.0,
    Gridtype            = 2,
    n_envopt            = NULL,
    use_parallel        = TRUE,
    use_opencl          = FALSE,
    verbose             = FALSE,
    progbar             = FALSE,
    stage_verbose       = FALSE,
    rate_calibration    = NULL,
    collect_block1      = TRUE
) {
  cl <- match.call()
  fn_name <- "rGLMM_reg_known_vcov"

  inp <- .rGLMM_validate_matrix_inputs(
    n, y, x, block, x_hyper, tv_tol,
    re_coef_names, group_levels, group_name, family, mode_gap_max,
    gap_tol, prior_list, pfamily_list
  )
  if (!inp$pf_summary$all_dNormal) {
    stop(
      fn_name, "(): all Block~2 components must be dNormal(); ",
      "use rGLMM_reg_estimated_vcov() or rGLMM_reg().",
      call. = FALSE
    )
  }

  .rGLMM_reg_run(
    inp                = inp,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    mode_gap_max       = mode_gap_max,
    verbose            = verbose,
    progbar            = progbar,
    stage_verbose      = stage_verbose,
    rate_calibration   = rate_calibration,
    collect_block1     = collect_block1,
    engine_label       = fn_name,
    result_class       = "rGLMM_reg_known_vcov",
    cl                 = cl
  )
}

#' @describeIn rGLMM_reg ING Block~2 (estimated \eqn{\tau^2_k}).
#'   Non-Gaussian: pilot always (unless \code{n_pilot = 0L}); conservative
#'   \code{disp_lower} eigenvalue bounds and pilot-updated \eqn{\tau^2} starts.
#' @export
rGLMM_reg_estimated_vcov <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    prior_list,
    pfamily_list,
    icm_tol             = 1e-10,
    icm_maxit           = 200L,
    offset              = NULL,
    weights             = 1,
    family              = gaussian(),
    re_coef_names       = colnames(x),
    group_levels        = levels(block),
    group_name          = NULL,
    gap_tol             = 0.0196,
    tv_tol              = 0.01,
    mode_gap_max        = 1.0,
    Gridtype            = 2,
    n_envopt            = NULL,
    use_parallel        = TRUE,
    use_opencl          = FALSE,
    verbose             = FALSE,
    progbar             = FALSE,
    stage_verbose       = FALSE,
    rate_calibration    = NULL,
    collect_block1      = TRUE
) {
  cl <- match.call()
  fn_name <- "rGLMM_reg_estimated_vcov"

  inp <- .rGLMM_validate_matrix_inputs(
    n, y, x, block, x_hyper, tv_tol,
    re_coef_names, group_levels, group_name, family, mode_gap_max,
    gap_tol, prior_list, pfamily_list
  )
  if (inp$pf_summary$all_dNormal) {
    stop(
      fn_name, "(): at least one Block~2 component must not be dNormal(); ",
      "use rGLMM_reg_known_vcov() or rGLMM_reg().",
      call. = FALSE
    )
  }

  .rGLMM_reg_run_with_pilot(
    inp                = inp,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    mode_gap_max       = mode_gap_max,
    verbose            = verbose,
    progbar            = progbar,
    stage_verbose      = stage_verbose,
    rate_calibration   = rate_calibration,
    collect_block1     = collect_block1,
    engine_label       = fn_name,
    result_class       = "rGLMM_reg_estimated_vcov",
    cl                 = cl
  )
}

#' @describeIn rGLMM_reg Route by Block~2 \code{pfamily_list} to known or
#'   estimated \eqn{\tau^2} engines.
#' @export
rGLMM_reg <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    prior_list,
    pfamily_list,
    icm_tol             = 1e-10,
    icm_maxit           = 200L,
    offset              = NULL,
    weights             = 1,
    family              = gaussian(),
    re_coef_names       = colnames(x),
    group_levels        = levels(block),
    group_name          = NULL,
    gap_tol             = 0.0196,
    tv_tol              = 0.01,
    mode_gap_max        = 1.0,
    Gridtype            = 2,
    n_envopt            = NULL,
    use_parallel        = TRUE,
    use_opencl          = FALSE,
    verbose             = FALSE,
    progbar             = FALSE,
    stage_verbose       = FALSE,
    rate_calibration    = NULL,
    collect_block1      = TRUE
) {
  cl <- match.call()

  inp <- .rGLMM_validate_matrix_inputs(
    n, y, x, block, x_hyper, tv_tol,
    re_coef_names, group_levels, group_name, family, mode_gap_max,
    gap_tol, prior_list, pfamily_list
  )

  route_fn <- if (inp$pf_summary$all_dNormal) {
    rGLMM_reg_known_vcov
  } else {
    rGLMM_reg_estimated_vcov
  }
  mc <- match.call(expand.dots = FALSE)
  mc[[1L]] <- route_fn
  out <- eval(mc, parent.frame())
  out$call <- cl
  out
}
