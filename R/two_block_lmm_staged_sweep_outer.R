#' Format v6 sweep-outer batch output for LMM staged \code{fixef.*} naming
#' @noRd
.rLMM_format_sweep_out <- function(
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
    sweep_history          = v6_out$sweep_history,
    n                      = n
  )
  .two_block_as_staged_names(
    x,
    fixef_mode = fixef_mode,
    fixef_init = fixef_init
  )
}

#' ING LMM replicate chains: pilot then main via sweep-outer v6 (mirrors rGLMM_reg)
#' @noRd
.rLMMNormal_reg_run_with_pilot <- function(
    inp,
    group,
    P,
    dispersion,
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
    engine_label  = "rLMMNormal_reg_estimated_vcov",
    result_class  = "rLMMNormal_reg_estimated_vcov",
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
  any_non_normal   <- TRUE
  ptypes           <- pf_summary$ptypes

  prior_list_block1 <- list(
    P          = P,
    dispersion = dispersion,
    ddef       = FALSE
  )
  .two_block_validate_block1_prior(prior_list_block1, family = family)

  gap_tol <- .two_block_validate_gap_tol(gap_tol)

  will_pilot <- .two_block_pilot_will_run(
    is_gaussian,
    n_pilot_arg,
    gap_tol,
    tv_tol,
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

  icm <- .rLMM_icm_at_start(
    y                 = inp$y,
    D                 = inp$D,
    group             = group,
    W                 = inp$W,
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
    re_names          = re_names,
    group_levels      = group_levels,
    group_name        = group_name,
    icm_tol           = icm_tol,
    icm_maxit         = icm_maxit,
    verbose           = verbose,
    engine_label      = engine_label
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
    D             = inp$D,
    group         = factor(group, levels = group_levels),
    W             = inp$W,
    groupef.names = re_names,
    group_name    = group_name
  )

  fixef_mode_ref <- fixef_mode
  b_mode_ref     <- b_start
  progbar_use    <- isTRUE(progbar) || isTRUE(verbose) || isTRUE(stage_verbose)

  rate <- two_block_rate_from_pfamily_list(
    x                 = inp$D,
    group             = group,
    x_hyper           = inp$W,
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
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
    any_non_normal      = any_non_normal
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

  calib_label <- paste0(
    "exact (Gaussian posterior)",
    "; conservative: non-dNormal RE dispersion (disp_lower plug-in)"
  )

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
    method              = "exact+disp_lower_bound",
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
  tau2_start_main    <- if (!is.null(icm) && !is.null(icm$tau2_start)) {
    icm$tau2_start
  } else {
    .two_block_tau2_start_from_pfamily(pfamily_list, re_names)
  }

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
      block1_prior   = prior_list_block1,
      pfamily_list   = pfamily_list,
      family         = family,
      re_names       = re_names,
      group_levels   = group_levels,
      collect_block1 = collect_block1,
      progbar        = progbar_use,
      stage_label    = "pilot",
      diag_sweeps    = isTRUE(diag_sweeps),
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
        x                  = inp$D,
        group              = group,
        x_hyper            = inp$W,
        prior_list         = prior_list_block1,
        pfamily_list       = pfamily_list,
        family             = family,
        tv_tol             = tv_tol,
        dispersion         = dispersion
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

    pilot_res <- .rLMM_format_sweep_out(
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
    block1_prior   = prior_list_block1,
    pfamily_list   = pfamily_list,
    family         = family,
    re_names       = re_names,
    group_levels   = group_levels,
    collect_block1 = collect_block1,
    progbar        = progbar_use,
    stage_label    = "main",
    diag_sweeps    = isTRUE(diag_sweeps),
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
    block1_prior   = prior_list_block1,
    pfamily_list   = pfamily_list,
    family         = family,
    re_names       = re_names,
    group_levels   = group_levels,
    collect_block1 = collect_block1,
    progbar        = progbar_use,
    stage_label    = "main",
    diag_sweeps    = isTRUE(diag_sweeps),
    fixef_mode     = fixef_mode_ref,
    b_mode         = b_mode_ref,
    b_start        = b_mode_ref,
    ptypes         = ptypes,
    tau2_start     = tau2_start_main
  )

  main_res <- .rLMM_format_sweep_out(
    v6_out       = main_raw,
    n            = n,
    re_names     = re_names,
    group_levels = group_levels,
    fixef_mode   = fixef_mode,
    fixef_init   = fixef_init
  )

  .lmebayes_assemble_reg_result(
    staged              = main_res,
    call                = cl,
    m_convergence       = m_convergence_used,
    convergence_info    = convergence_info,
    pfamily_list        = pfamily_list,
    dispprior_list      = prior_list_block1,
    family              = family,
    groupef.mode        = ranef_mode,
    any_non_normal      = TRUE,
    design              = design,
    result_class        = result_class,
    parent_class        = "rLMMNormal_reg",
    draw_engine         = "rGLMM_sweep",
    sim_method_used     = "TWO_BLOCK_GIBBS",
    icm_info            = icm_info,
    pilot_draws         = if (run_pilot) pilot_res else NULL,
    n_pilot             = if (run_pilot) n_pilot else NULL,
    m_convergence_pilot = if (run_pilot) m_convergence_pilot else NULL,
    pilot_chisq         = if (run_pilot) pilot_chisq else NULL,
    pilot_ub            = if (run_ub) pilot_ub else NULL,
    tv_tol              = if (run_ub) tv_tol else NULL,
    offset              = inp$offset,
    weights             = if (is.null(inp$weights)) 1 else inp$weights
  )
}
