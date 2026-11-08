#' Two-block Gibbs sweep for replicate-chain sampling (\code{rGLMM_reg} engine)
#'
#' Executes \code{n_chains} independent two-block Gibbs chains, each
#' initialised at \code{start_fixef} and run for \code{inner_sweeps} sweeps.
#' Loop order is sweep-outer: for each inner sweep, all chains run Block~1,
#' then all chains run Block~2. Each chain stores exactly one draw
#' (\code{n = 1L} per chain), so the result is a collection of
#' \code{n_chains} approximately independent draws from the target posterior
#' (for large enough \code{inner_sweeps}).
#'
#' @param n_chains Integer. Number of independent chains to run.
#' @param start_fixef Named list of starting hyper-parameter vectors, one
#'   named numeric vector per RE component (matches \code{re_names}).
#' @param inner_sweeps Integer. Number of inner Gibbs sweeps per chain.
#' @param design List with components \code{y}, \code{D}, \code{group},
#'   \code{W}, \code{groupef.names}, and optional \code{group_name}.
#' @param block1_prior Block~1 prior list (\code{P} or \code{Sigma},
#'   \code{dispersion}, optional \code{ddef}).
#' @param pfamily_list Named list of \code{\link[glmbayesCore]{pfamily}} objects (one per
#'   RE component).
#' @param family A \code{\link[stats]{family}} object for the response model.
#' @param re_names Character vector of RE coefficient names.
#' @param group_levels Character vector of group levels.
#' @param collect_block1 Logical. If \code{TRUE}, collect and rbind Block~1
#'   (\code{coefficients}) draws from every chain. Default \code{TRUE}.
#' @param progbar Logical. When \code{TRUE}, show text progress bars over chains
#'   during each inner sweep. Default \code{FALSE}.
#' @param stage_label Character label stored on \code{$sweep_history} (e.g.
#'   \code{"pilot"} or \code{"main"}).
#' @param diag_sweeps When \code{TRUE}, print one combined Block~2 fixef
#'   chain-mean table for the stage when the inner-sweep loop finishes (same
#'   layout as \code{print()} on \code{$sweep_history} with all sweeps).
#'   Sweep history is always collected; use \code{print()} on
#'   \code{$sweep_history} to display tables later.
#' @param fixef_mode ICM mode reference stored on \code{$sweep_history}.
#' @param b_mode ICM mode reference for random-effect diagnostics.
#' @param b_start Initial random-effect matrix for all chains (\code{J x p_re}).
#'   Defaults to \code{b_mode} when \code{NULL}.
#' @param ptypes Per-component pfamily names (optional; derived from
#'   \code{pfamily_list} when \code{NULL}).
#' @param tau2_start Optional named numeric vector of plug-in
#'   \eqn{\tau^2_k} values for chain initialisation (one per \code{re_names}).
#'   When \code{NULL}, derived from \code{pfamily_list} prior fields
#'   (\code{dNormal} dispersion or ING \code{rate/shape} = \eqn{1/E[1/\tau^2]}).
#' @param use_cpp_block2 Must be \code{TRUE}. Block~2 uses
#'   \code{\link{two_block_block2_one_chain_cpp}} (native C++ align + \code{rglmb}).
#' @param use_cpp_tau2_row Accepted and ignored. Block~1 indexes \code{tau2} in R.
#' @param use_cpp_b_slice Accepted and ignored. Block~1 stores \code{b} in R.
#' @param use_cpp_iters_ranef_add Accepted and ignored. Block~1 accumulates
#'   envelope iteration counts in R.
#' @param rwt Ridge weight in \code{(0, 1]}, default \code{0.01}. Forwarded to
#'   Block~1 and not used by the current draw or sweep count. The default
#'   leaves existing callers unchanged.
#' @return A list with components \code{fixef_draws}, \code{dispersion_fixef_draws},
#'   \code{iters_fixef_draws}, \code{iters_ranef_draws}, \code{coefficients},
#'   \code{mu_all_last}, and \code{sweep_history} (class
#'   \code{"two_block_sweep_history"}).
#' @family simfuncs
#' @seealso \code{\link{two_block_block2_one_chain_cpp}}
#' @keywords internal
rGLMM_sweep <- function(
    n_chains,
    start_fixef,
    inner_sweeps,
    design,
    block1_prior,
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
    use_cpp_tau2_row = TRUE,
    use_cpp_b_slice = TRUE,
    use_cpp_iters_ranef_add = TRUE,
    use_cpp_block2 = TRUE,
    rwt = 0.01
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
  sweep_cov   <- vector("list", inner_sweeps)
  sweep_disp  <- vector("list", inner_sweeps)

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
    # .two_block_print_sweep_boundary(
    #   stage_label  = stage_label,
    #   sweep        = m,
    #   inner_sweeps = inner_sweeps,
    #   phase        = "Block1",
    #   boundary     = "enter"
    # )
    b1 <- .two_block_block1_all_chains(
      n                       = batch$n,
      fixef                   = batch$fixef,
      tau2                    = batch$tau2,
      b                       = batch$b,
      iters_ranef             = batch$iters_ranef,
      re_names                = batch$re_names,
      group_levels            = batch$group_levels,
      design                  = design,
      block1_prior            = block1_prior,
      family                  = family,
      ptypes                  = ptypes,
      progbar                 = progbar_use,
      progbar_prefix          = prefix_b1,
      progbar_finish_newline  = FALSE,
      use_cpp_tau2_row        = use_cpp_tau2_row,
      use_cpp_b_slice         = use_cpp_b_slice,
      use_cpp_iters_ranef_add = use_cpp_iters_ranef_add,
      rwt                     = rwt
    )

    batch$b           <- b1$b
    batch$iters_ranef <- b1$iters_ranef

      b2 <- .two_block_block2_all_chains(
      n                       = batch$n,
      b                       = batch$b,
      fixef                   = batch$fixef,
      tau2                    = batch$tau2,
      iters                   = batch$iters,
      re_names                = batch$re_names,
      group_levels            = batch$group_levels,
      design                  = design,
      pfamily_list            = pfamily_list,
      ptypes                  = ptypes,
      use_cpp_block2          = use_cpp_block2,
      progbar                 = progbar_use,
      progbar_prefix          = prefix_b2,
      progbar_finish_newline  = (m == inner_sweeps)
    )
    batch$fixef <- b2$fixef
    batch$tau2  <- b2$tau2
    batch$iters <- b2$iters
    # if (m <= 2L) {
    #   .two_block_print_sweep_early_diagnostics(
    #     sweep          = m,
    #     stage_label    = stage_label,
    #     batch          = batch,
    #     design         = design,
    #     re_names       = re_names,
    #     group_levels   = group_levels,
    #     use_cpp_mu_all = FALSE
    #   )
    # }
    # .two_block_print_sweep_boundary(
    #   stage_label  = stage_label,
    #   sweep        = m,
    #   inner_sweeps = inner_sweeps,
    #   phase        = "Block2",
    #   boundary     = "exit"
    # )
    sweep_stats[[m]] <- .two_block_snapshot_fixef_stats(
      fixef    = batch$fixef,
      re_names = re_names
    )
    sweep_cov[[m]] <- .two_block_snapshot_fixef_cov(
      fixef    = batch$fixef,
      re_names = re_names
    )
    sweep_disp[[m]] <- .two_block_snapshot_disp_stats(
      tau2     = batch$tau2,
      ptypes   = ptypes,
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

  out <- .rGLMM_sweep_save(
    n              = batch$n,
    fixef          = batch$fixef,
    tau2           = batch$tau2,
    b              = batch$b,
    iters          = batch$iters,
    iters_ranef    = batch$iters_ranef,
    re_names       = batch$re_names,
    group_levels   = batch$group_levels,
    design         = design,
    collect_block1 = collect_block1
  )
  out$sweep_history <- .two_block_build_sweep_history(
    stage_label = stage_label,
    sweep_stats = sweep_stats,
    fixef_mode  = fixef_mode,
    re_names    = re_names,
    sweep_cov   = sweep_cov,
    disp_stats  = sweep_disp
  )
  if (isTRUE(diag_sweeps)) {
    print(out$sweep_history)
  }
  invisible(out)
}
