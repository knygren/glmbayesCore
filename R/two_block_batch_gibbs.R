# Pure-R batch driver for two-block Gibbs sampling (sweep-outer layout).
#
# Two-block Gibbs alternates:
#   Block 1 — given hyperparameters (fixef gamma, tau2), sample random effects b
#             per neighborhood via rNormalGLM_reg_group / rNormal_reg_group.
#   Block 2 — given b, treat each RE column as a Gaussian pseudo-response and
#             update hyperparameters via rglmb (one call per RE component).
#
# Sweep-outer: ALL chains run Block 1, then ALL chains run Block 2, per inner
# sweep. Block 1 is split into prep (mu_all + prior_list per chain) and draw
# (rNormalGLM_reg_group / rNormal_reg_group per chain); both phases are embarrassingly
# parallel over chains (optional n_cores on Unix/macOS).


#' Build a progress-bar prefix for sweep-outer Gibbs stages
#' @noRd
.two_block_progbar_prefix <- function(stage_label, sweep, inner_sweeps, phase) {
  phase_label <- if (identical(phase, "Block1")) {
    "RE"
  } else {
    "fixef"
  }
  stage_label <- as.character(stage_label)[1L]
  if (nzchar(stage_label)) {
    sprintf("[%s] sweep %d/%d %s: ", stage_label, sweep, inner_sweeps, phase_label)
  } else {
    sprintf("sweep %d/%d %s: ", sweep, inner_sweeps, phase_label)
  }
}

#' Text progress bar matching glmbayesCore C++ style
#' @param current Completed step (1-based, up to \code{total}).
#' @param total Total number of steps.
#' @param prefix Optional label printed before the bar (not cleared by \code{\r}).
#' @noRd
.two_block_progress_bar <- function(current, total, prefix = "") {
  total <- as.integer(total[1L])
  current <- as.integer(current[1L])
  if (!is.finite(total) || total <= 0L) {
    return(invisible())
  }
  if (!is.finite(current)) {
    current <- 0L
  }
  current <- max(0L, min(current, total))
  totaldotz <- 40L
  fraction  <- current / total
  if (!is.finite(fraction)) {
    return(invisible())
  }
  dotz      <- round(fraction * totaldotz)
  cat("\r", strrep(" ", 100L), "\r", sep = "")
  if (nzchar(prefix)) {
    cat(prefix)
  }
  cat(sprintf("%3.0f%% [", fraction * 100), sep = "")
  cat(paste0(rep("=", dotz), collapse = ""))
  cat(paste0(rep(" ", totaldotz - dotz), collapse = ""))
  cat("]", sep = "")
  utils::flush.console()
}

#' Finish a progress bar started by \code{.two_block_progress_bar}
#' @param newline If \code{FALSE}, leave the completed bar on the current
#'   line so the next bar can overwrite it with \code{\r}.
#' @noRd
.two_block_progress_bar_finish <- function(newline = TRUE) {
  if (isTRUE(newline)) {
    cat("\n")
  }
}


#' Lookup ICM mode for one fixef (re_name, covariate) pair
#' @noRd
.two_block_fixef_mode_at <- function(fixef_mode, re_name, cov_name, col_names) {
  if (is.null(fixef_mode) || is.null(fixef_mode[[re_name]])) {
    return(NA_real_)
  }
  mode_k <- fixef_mode[[re_name]]
  if (!is.null(names(mode_k)) && cov_name %in% names(mode_k)) {
    return(unname(mode_k[[cov_name]]))
  }
  j <- match(cov_name, col_names)
  if (is.finite(j) && length(mode_k) >= j) {
    return(unname(mode_k[j]))
  }
  NA_real_
}


#' Snapshot chain colMeans and SDs of Block 2 fixef after one sweep
#' @param fixef Named list of fixed-effect matrices (\code{n x p_k} per RE block).
#' @param re_names Character vector of random-effect block names.
#' @return Named list of \code{list(mean = ..., sd = ...)} per \code{re_names} entry.
#' @noRd
.two_block_snapshot_fixef_stats <- function(fixef, re_names) {
  out <- list()
  for (k in re_names) {
    mat <- fixef[[k]]
    cn  <- colnames(mat)
    if (is.null(cn)) {
      cn <- paste0("V", seq_len(ncol(mat)))
    }
    out[[k]] <- list(
      mean = stats::setNames(colMeans(mat), cn),
      sd   = stats::setNames(apply(mat, 2, stats::sd), cn)
    )
  }
  out
}

#' Snapshot chain mean/SD of RE precision (\code{1/tau^2}) after one sweep,
#' ING components only
#'
#' Companion to \code{\link{.two_block_snapshot_fixef_stats}} for the
#' estimated-vcov case: \code{Sigma_ranef = diag(tau2_k)}, and only
#' \code{dIndependent_Normal_Gamma} (ING) components actually sample
#' \code{tau2_k} each sweep (\code{dNormal} components keep a fixed
#' dispersion, so there is nothing to track for them). Reports the
#' precision \eqn{1/\tau^2_k}, not \eqn{\tau^2_k} itself, because
#' \eqn{E[1/\tau^2_k]} (a Gamma mean) is always well-defined while
#' \eqn{E[\tau^2_k]} (an inverse-Gamma mean) is not for weak/low-shape
#' priors (see \code{\link{.two_block_tau2_start_from_dispersion_draws}}'s
#' use of the same \code{1/mean(1/x)} identity). Cross-chain *covariance*
#' (between precision components, or between precision and fixef) is not
#' captured here -- see \code{\link{plot_var_convergence}}'s
#' \code{component} docs for why \code{whitened = TRUE} does not (yet)
#' support \code{component = "precision"}.
#' @param tau2 \code{n_chains x p_re} matrix (\code{batch$tau2}).
#' @param ptypes Named character vector of Block~2 pfamily names, one per
#'   \code{re_names} entry (\code{"dNormal"} or
#'   \code{"dIndependent_Normal_Gamma"}).
#' @param re_names Character vector of random-effect block names (defines
#'   \code{colnames(tau2)} order).
#' @return Named list (one entry per ING component in \code{re_names}, in
#'   \code{re_names} order; empty list if none are ING) of
#'   \code{list(mean = <scalar>, sd = <scalar>)} for that component's
#'   cross-chain precision.
#' @noRd
.two_block_snapshot_disp_stats <- function(tau2, ptypes, re_names) {
  ing <- re_names[vapply(re_names, function(k) {
    identical(ptypes[[k]], "dIndependent_Normal_Gamma")
  }, logical(1))]
  out <- list()
  for (k in ing) {
    prec <- 1 / tau2[, k]
    out[[k]] <- list(mean = mean(prec), sd = stats::sd(prec))
  }
  out
}

#' Snapshot the full cross-chain covariance of stacked Block 2 fixef after one sweep
#'
#' Stacks \code{fixef[[k]]} (\code{n_chains x p_k}, one per RE component) into
#' a single \code{n_chains x P_total} matrix (same column order as
#' \code{.lmerb_posterior_normal_system()}'s \code{idx}: loop over
#' \code{re_names}, then \code{colnames(fixef[[k]])}) and returns its
#' cross-chain covariance -- the empirical analogue of Claim 3's
#' \eqn{\Sigma^{(l)}_{11}} (see \code{\link{plot_var_convergence}}).
#' Companion to \code{\link{.two_block_snapshot_fixef_stats}}, which only
#' keeps the marginal (diagonal) mean/sd.
#' @param fixef Named list of fixed-effect matrices (\code{n_chains x p_k} per
#'   RE block).
#' @param re_names Character vector of RE component names (defines stacking
#'   order).
#' @return List with \code{cov} (\code{P_total x P_total} matrix, dimnamed
#'   \code{"re_component | covariate"}) and \code{coef_index} (a data frame
#'   with one row per column of \code{cov}, giving \code{re_component} and
#'   \code{covariate} in stacking order).
#' @noRd
.two_block_snapshot_fixef_cov <- function(fixef, re_names) {
  mats <- lapply(re_names, function(k) {
    mat <- fixef[[k]]
    cn  <- colnames(mat)
    if (is.null(cn)) {
      cn <- paste0("V", seq_len(ncol(mat)))
    }
    colnames(mat) <- cn
    mat
  })
  coef_index <- do.call(rbind, lapply(seq_along(re_names), function(i) {
    data.frame(
      re_component = re_names[[i]],
      covariate    = colnames(mats[[i]]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(coef_index) <- NULL
  stacked <- do.call(cbind, mats)
  lbl <- paste(coef_index$re_component, coef_index$covariate, sep = " | ")
  cov_mat <- stats::cov(stacked)
  dimnames(cov_mat) <- list(lbl, lbl)
  list(cov = cov_mat, coef_index = coef_index)
}


#' Initialize sweep state for \code{\link{rGLMM_sweep}}
#' @noRd
.rGLMM_sweep_initialize <- function(
    n_chains,
    start_fixef,
    b_start,
    tau2_start,
    re_names,
    group_levels
) {
  p_re <- length(re_names)
  J    <- length(group_levels)

  fixef <- lapply(start_fixef, function(beta0) {
    matrix(
      beta0,
      nrow = n_chains,
      ncol = length(beta0),
      byrow = TRUE,
      dimnames = list(NULL, names(beta0))
    )
  })
  names(fixef) <- re_names

  tau2 <- matrix(
    tau2_start,
    nrow = n_chains,
    ncol = p_re,
    byrow = TRUE,
    dimnames = list(NULL, re_names)
  )

  b <- array(b_start, dim = c(J, p_re, n_chains),
             dimnames = list(group_levels, re_names, NULL))

  iters <- matrix(0, nrow = n_chains, ncol = p_re)
  colnames(iters) <- re_names

  iters_ranef <- numeric(n_chains)

  list(
    n     = n_chains,
    fixef = fixef,
    tau2  = tau2,
    b     = b,
    iters = iters,
    iters_ranef = iters_ranef,
    re_names     = re_names,
    group_levels = group_levels
  )
}

#' Restore \code{dim} and \code{dimnames} on batch \code{b} (Block 1 downstream)
#' @noRd
.two_block_ensure_batch_b_dimnames <- function(b, group_levels, re_names, n) {
  J <- length(group_levels)
  p_re <- length(re_names)
  expected_len <- J * p_re * n
  if (length(b) != expected_len) {
    stop(
      "batch$b length (", length(b), ") != J * p_re * n (",
      J, " * ", p_re, " * ", n, ").",
      call. = FALSE
    )
  }
  dim(b) <- c(J, p_re, n)
  dimnames(b) <- list(group_levels, re_names, NULL)
  b
}


#' Extract chain-i fixef list from batch state
#' @noRd
.two_block_batch_fixef_chain <- function(batch, i) {
  lapply(batch$fixef, function(mat) {
    v <- mat[i, , drop = TRUE]
    stats::setNames(as.numeric(v), colnames(mat))
  })
}

#' Refresh Block 1 prior precision for ING components (R reference; v5 any_ing)
#' @noRd
.two_block_block1_prior_with_tau2_r <- function(
    base_prior,
    tau2_vec,
    ptypes,
    re_names,
    mu_all
) {
  out <- list(
    mu         = mu_all,
    dispersion = base_prior$dispersion,
    ddef       = base_prior$ddef
  )

  if (!any(ptypes == "dIndependent_Normal_Gamma")) {
    out$P <- base_prior$P
    return(out)
  }

  P1 <- base_prior$P
  for (k in seq_along(re_names)) {
    if (ptypes[[k]] != "dIndependent_Normal_Gamma") next
    P1[k, ] <- 0
    P1[, k] <- 0
    P1[k, k] <- 1 / tau2_vec[k]
  }
  out$P <- P1
  out
}

#' Refresh Block 1 prior precision for ING components (mirrors C++ twoBlockGibbs)
#' @param use_cpp lmerb passes \code{FALSE}. \code{TRUE} is not implemented.
#' @noRd
.two_block_block1_prior_with_tau2 <- function(
    base_prior,
    tau2_vec,
    ptypes,
    re_names,
    mu_all,
    use_cpp = TRUE
) {
  if (isTRUE(use_cpp)) {
    stop("'use_cpp' must be FALSE.", call. = FALSE)
  }
  .two_block_block1_prior_with_tau2_r(
    base_prior, tau2_vec, ptypes, re_names, mu_all
  )
}

#' One-chain Block 1 prep: fixef -> mu_all -> prior_list (no sampling)
#' @noRd
.two_block_block1_prep_one_chain <- function(
    batch,
    i,
    design,
    block1_prior,
    ptypes,
    use_cpp_mu_all = TRUE,
    use_cpp_prior_tau2 = TRUE
) {
  fixef_i <- .two_block_batch_fixef_chain(batch, i)
  mu_all  <- as.matrix(build_mu_all(
    design, fixef_i, batch$group_levels, use_cpp = use_cpp_mu_all
  )$mu_all)
  tau2_i  <- batch$tau2[i, ]
  prior_list <- .two_block_block1_prior_with_tau2(
    block1_prior, tau2_i, ptypes, batch$re_names, mu_all,
    use_cpp = use_cpp_prior_tau2
  )
  list(mu_all = mu_all, prior_list = prior_list)
}

#' All-chain Block 1 prep: mu_all and prior_list for every chain
#' @noRd
.two_block_block1_prep_all_chains <- function(
    batch,
    design,
    block1_prior,
    ptypes,
    n_cores = NULL,
    progbar = FALSE,
    progbar_finish_newline = TRUE,
    use_cpp_mu_all = TRUE,
    use_cpp_prior_tau2 = TRUE
) {
  n <- batch$n
  show_bar <- isTRUE(progbar) && n > 1L &&
    (is.null(n_cores) || as.integer(n_cores[1L]) < 2L)

  prep_i <- function(i) {
    if (show_bar) .two_block_progress_bar(i, n)
    .two_block_block1_prep_one_chain(
      batch              = batch,
      i                  = i,
      design             = design,
      block1_prior       = block1_prior,
      ptypes             = ptypes,
      use_cpp_mu_all     = use_cpp_mu_all,
      use_cpp_prior_tau2 = use_cpp_prior_tau2
    )
  }

  prep_list <- .two_block_lapply_chains(n, prep_i, n_cores = n_cores)
  if (show_bar) .two_block_progress_bar_finish(newline = progbar_finish_newline)

  structure(
    list(
      mu_all      = lapply(prep_list, `[[`, "mu_all"),
      prior_lists = lapply(prep_list, `[[`, "prior_list")
    ),
    class = "two_block_block1_prep"
  )
}

#' Mean envelope candidates per group from a Block~1 draw (R reference)
#' @noRd
.two_block_block1_iters_mean_r <- function(block_out) {
  br <- block_out$group_results
  if (is.null(br) || !length(br)) {
    return(1)
  }
  vals <- vapply(br, function(b) {
    if (is.null(b$iters)) {
      return(1)
    }
    it <- b$iters
    if (is.matrix(it)) {
      return(as.numeric(it[1, 1]))
    }
    as.numeric(it[1])
  }, numeric(1))
  mean(vals)
}

#' Mean envelope candidates per group from a Block~1 draw
#' @param block_out Output from \code{rNormalGLM_reg_group} or \code{rNormal_reg_group}.
#' @param use_cpp lmerb passes \code{FALSE}. \code{TRUE} is not implemented.
#' @noRd
.two_block_block1_iters_mean <- function(block_out, use_cpp = TRUE) {
  if (isTRUE(use_cpp)) {
    stop("'use_cpp' must be FALSE.", call. = FALSE)
  }
  .two_block_block1_iters_mean_r(block_out)
}

#' Reorder Block~1 \code{b} rows to \code{group_levels} (R reference)
#' @noRd
.two_block_block1_reorder_b_r <- function(b_draw, group_levels) {
  rn <- rownames(b_draw)
  if (is.null(rn)) {
    return(b_draw)
  }
  ord <- match(group_levels, rn)
  if (any(is.na(ord))) {
    stop("Block 1 group ids do not match group_levels.", call. = FALSE)
  }
  b_draw[ord, , drop = FALSE]
}

#' Reorder Block~1 \code{b} rows to \code{group_levels}
#' @param b_draw Coefficient matrix from a Block~1 draw.
#' @param group_levels Target row order.
#' @param block_ids Row ids for \code{b_draw}; defaults to \code{rownames(b_draw)}.
#' @param use_cpp lmerb passes \code{FALSE}. \code{TRUE} is not implemented.
#' @noRd
.two_block_block1_reorder_b <- function(
    b_draw,
    group_levels,
    block_ids = rownames(b_draw),
    use_cpp = TRUE
) {
  if (isTRUE(use_cpp)) {
    stop("'use_cpp' must be FALSE.", call. = FALSE)
  }
  .two_block_block1_reorder_b_r(b_draw, group_levels)
}

#' Block~1 envelope draw for one chain (\code{rNormal_reg_group}/\code{rNormalGLM_reg_group})
#' @noRd
.two_block_block1_draw_block <- function(
    prior_list,
    design,
    family,
    is_gaussian
) {
  if (is_gaussian) {
    return(rNormal_reg_group(
      n          = 1L,
      y          = design$y,
      x          = design$D,
      group      = design$group,
      prior_list = prior_list
    ))
  }
  rNormalGLM_reg_group(
    n            = 1L,
    y            = design$y,
    x            = design$D,
    group        = design$group,
    prior_list   = prior_list,
    family       = family,
    use_parallel = FALSE,
    verbose      = FALSE,
    progbar      = FALSE
  )
}

#' Reorder Block~1 draw coefficients and summarize envelope iterations
#' @noRd
.two_block_block1_draw_reorder <- function(
    block_out,
    group_levels,
    use_cpp_reorder = TRUE,
    use_cpp_iters = TRUE
) {
  b_draw <- block_out$coefficients
  b_draw <- .two_block_block1_reorder_b(
    b_draw       = b_draw,
    group_levels = group_levels,
    use_cpp      = use_cpp_reorder
  )
  list(
    b          = b_draw,
    iters_mean = .two_block_block1_iters_mean(block_out, use_cpp = use_cpp_iters)
  )
}

#' Block~1 random-effect redraw for replicate chains
#'
#' Runs one \code{\link{rNormalGLM_reg_group}} or \code{\link{rNormal_reg_group}}
#' envelope draw per chain from prepared Block~1 \code{prior_list} objects.
#' This is the Block~1 draw step inside \code{\link{rGLMM_sweep}} after
#' per-chain prep (\code{fixef} \eqn{\to} \code{mu_all}, \eqn{\tau^2} refresh
#' of \code{P}).
#'
#' Each draw uses \code{n = 1L} per chain. Outputs are raw block-sampler
#' results; reorder rows of \code{$coefficients} to \code{group_levels} and
#' summarize envelope iterations separately (as \code{rGLMM_sweep} does after
#' this step).
#'
#' @param n Integer. Number of independent chains; must equal
#'   \code{length(prior_lists)}.
#' @param prior_lists List of length \code{n}. Each element is a Block~1
#'   \code{prior_list} with \code{mu} (\code{p_re x J} matrix, rows =
#'   \code{re_names}, columns = group levels in block order),
#'   \code{P} (\code{p_re x p_re} precision matrix), \code{dispersion}, and
#'   optional \code{ddef}. Typically built by the prep step that pairs
#'   \code{fixef[[k]][i, ]} with \code{tau2[i, k]} for chain \code{i}.
#' @param design List with at least \code{y} (response), \code{D} (level-1
#'   design, \code{l2 x p_re}), and \code{group} (block partition of length
#'   \code{l2}). Same layout as the \code{design} argument to
#'   \code{\link{rGLMM_sweep}}.
#' @param family A \code{\link[stats]{family}} object. When \code{gaussian()},
#'   uses \code{\link{rNormal_reg_group}}; otherwise
#'   \code{\link{rNormalGLM_reg_group}} with default offset zero and unit weights.
#' @param progbar Logical. When \code{TRUE} and \code{n > 1}, show a text
#'   progress bar over chains. Default \code{FALSE}.
#' @param progbar_prefix Character prefix for progress-bar lines. Default
#'   \code{""}.
#' @param progbar_finish_newline Logical. Print a newline after the progress
#'   bar finishes. Default \code{TRUE}.
#' @return List of length \code{n}. Each element is a \code{rNormalGLM_reg_group}
#'   or \code{rNormal_reg_group} object (\code{$coefficients}, \code{$coef.mode},
#'   \code{$group_info}, \code{$group_results}, \ldots).
#' @family simfuncs
#' @seealso \code{\link{rGLMM_sweep}}, \code{\link{rNormalGLM_reg_group}},
#'   \code{\link{rNormal_reg_group}}
#' @keywords internal
rGLMM_Re_Draw <- function(
    n,
    prior_lists,
    design,
    family,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE
) {
  if (length(prior_lists) != n) {
    stop(
      "length(prior_lists) (", length(prior_lists), ") must equal n (", n, ").",
      call. = FALSE
    )
  }
  is_gaussian <- identical(family$family, "gaussian")
  show_bar <- isTRUE(progbar) && n > 1L
  block_outs <- vector("list", n)
  for (i in seq_len(n)) {
    if (show_bar) .two_block_progress_bar(i, n, prefix = progbar_prefix)
    block_outs[[i]] <- .two_block_block1_draw_block(
      prior_list  = prior_lists[[i]],
      design      = design,
      family      = family,
      is_gaussian = is_gaussian
    )
  }
  if (show_bar) {
    .two_block_progress_bar_finish(newline = progbar_finish_newline)
  }
  block_outs
}


#' All-chain Block 1 draw from prepared prior_lists
#' @noRd
.two_block_block1_draw_all_chains <- function(
    batch,
    prep,
    design,
    family,
    n_cores = NULL,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE,
    use_cpp_reorder = TRUE,
    use_cpp_iters = TRUE
) {
  n <- batch$n
  show_bar <- isTRUE(progbar) && n > 1L &&
    (is.null(n_cores) || as.integer(n_cores[1L]) < 2L)
  prior_lists <- prep$prior_lists
  if (length(prior_lists) != n) {
    stop("length(prep$prior_lists) must equal batch$n.", call. = FALSE)
  }

  block_outs <- rGLMM_Re_Draw(
    n                      = n,
    prior_lists            = prior_lists,
    design                 = design,
    family                 = family,
    progbar                = show_bar,
    progbar_prefix         = progbar_prefix,
    progbar_finish_newline = progbar_finish_newline
  )

  reorder_i <- function(i) {
    .two_block_block1_draw_reorder(
      block_out       = block_outs[[i]],
      group_levels    = batch$group_levels,
      use_cpp_reorder = use_cpp_reorder,
      use_cpp_iters   = use_cpp_iters
    )
  }

  b_draws <- .two_block_lapply_chains(n, reorder_i, n_cores = NULL)

  for (i in seq_len(n)) {
    batch$b[, , i] <- b_draws[[i]]$b
    batch$iters_ranef[i] <- batch$iters_ranef[i] + b_draws[[i]]$iters_mean
  }
  batch
}

#' Apply FUN to each chain index, optionally in parallel (Unix/macOS only)
#' @noRd
.two_block_lapply_chains <- function(n, FUN, n_cores = NULL) {
  idx <- seq_len(n)
  if (is.null(n_cores)) {
    return(lapply(idx, FUN))
  }
  n_cores <- as.integer(n_cores[1L])
  if (!is.finite(n_cores) || n_cores < 2L) {
    return(lapply(idx, FUN))
  }
  n_cores <- min(n_cores, n)
  if (.Platform$OS.type == "windows") {
    warning(
      "Chain-parallel Block 1 (n_cores > 1) is not supported on Windows; ",
      "using sequential lapply.",
      call. = FALSE
    )
    return(lapply(idx, FUN))
  }
  parallel::mclapply(idx, FUN, mc.cores = n_cores)
}


#' Block 2 one-chain update via C++ (native align + \code{rglmb})
#'
#' \code{b} is aligned to \code{X_hyper} rows in C++ before each \code{rglmb()}
#' call. This is the population update used by \code{\link{rGLMM_sweep}}.
#'
#' @param i Chain index (\code{1..n}).
#' @param b \code{J x p_re x n} array of random effects (read-only).
#' @param fixef Named list of fixed-effect matrices (\code{n x p_k} per RE block).
#' @param tau2 \code{n x p_re} matrix of random-effect variances.
#' @param iters \code{n x p_re} matrix of Block~2 iteration counts.
#' @param re_names Character vector of random-effect block names.
#' @param group_levels Character vector of group level labels.
#' @param design Model design list (\code{W}, etc.).
#' @param pfamily_list Named list of Block~2 \code{pfamily} objects.
#' @param ptypes Named character vector of \code{pfamily} types.
#' @return List with \code{fixef}, \code{tau2}, and \code{iters}.
#' @keywords internal
two_block_block2_one_chain_cpp <- function(
    i,
    b,
    fixef,
    tau2,
    iters,
    re_names,
    group_levels,
    design,
    pfamily_list,
    ptypes
) {
  b_i <- b[, , i, drop = FALSE]
  b_i <- matrix(
    b_i, nrow = nrow(b_i), ncol = ncol(b_i),
    dimnames = dimnames(b_i)[1:2]
  )
  fixef_rows <- lapply(re_names, function(k) fixef[[k]][i, ])
  names(fixef_rows) <- re_names
  x_hyper <- lapply(design$W, as.matrix)
  out <- .two_block_block2_one_chain_cpp(
    b_i            = b_i,
    fixef_rows     = fixef_rows,
    tau2_i         = tau2[i, ],
    iters_i        = iters[i, ],
    x_hyper        = x_hyper,
    group_levels   = group_levels,
    pfamily_list   = pfamily_list,
    ptypes         = ptypes,
    re_names       = re_names
  )

  fixef_out <- fixef
  for (k in re_names) {
    fixef_out[[k]][i, ] <- out$fixef[[k]]
  }
  tau2_out <- tau2
  tau2_out[i, ] <- out$tau2
  iters_out <- iters
  iters_out[i, ] <- out$iters

  list(
    fixef = fixef_out,
    tau2  = tau2_out,
    iters = iters_out
  )
}


#' Assign one chain slice of \code{batch$b} (R reference; all-chains step C)
#' @noRd
.two_block_batch_b_assign_slice_r <- function(b_store, chain_i, b_draw) {
  b_store[, , chain_i] <- b_draw
  b_store
}

#' Assign \code{batch$b[, , chain_i] <- b_draw} (all-chains step C)
#'
#' Default \code{TRUE} on \code{\link{rGLMM_sweep}}; pass
#' \code{use_cpp_b_slice = FALSE} to revert to R \code{batch$b[, , i] <- b_draw}.
#'
#' @param b_store 3-D array \code{batch$b} (\code{dim = c(J, p_re, n)}).
#' @param chain_i Chain index (\code{1..n}).
#' @param b_draw \code{J x p_re} matrix (\code{out$b} from one-chain export).
#' @param use_cpp_b_slice When \code{TRUE}, use C++ slice export; \code{FALSE}
#'   uses the R reference subassignment.
#' @noRd
.two_block_batch_b_assign_slice <- function(
    b_store,
    chain_i,
    b_draw,
    use_cpp_b_slice = FALSE
) {
  if (isTRUE(use_cpp_b_slice)) {
    stop("'use_cpp_b_slice' must be FALSE.", call. = FALSE)
  }
  .two_block_batch_b_assign_slice_r(b_store, chain_i, b_draw)
}

#' Add one-chain Block~1 envelope iters to \code{batch$iters_ranef} (R reference; step D)
#' @noRd
.two_block_batch_iters_ranef_add_r <- function(iters_ranef, chain_i, iters_mean) {
  stopifnot(length(iters_ranef) >= chain_i)
  iters_ranef[chain_i] <- iters_ranef[chain_i] + as.numeric(iters_mean)
  iters_ranef
}

#' \code{batch$iters_ranef[chain_i] <- batch$iters_ranef[chain_i] + iters_mean}
#' (all-chains step D)
#'
#' Default \code{TRUE} on \code{\link{rGLMM_sweep}}; pass
#' \code{use_cpp_iters_ranef_add = FALSE} to revert to R accumulation.
#'
#' @param iters_ranef Length-\code{n} vector (\code{batch$iters_ranef}).
#' @param chain_i Chain index (\code{1..n}).
#' @param iters_mean Scalar from \code{out$iters_mean} (one-chain export).
#' @param use_cpp_iters_ranef_add When \code{TRUE}, use C++ add export; \code{FALSE}
#'   uses the R reference.
#' @noRd
.two_block_batch_iters_ranef_add <- function(
    iters_ranef,
    chain_i,
    iters_mean,
    use_cpp_iters_ranef_add = FALSE
) {
  if (isTRUE(use_cpp_iters_ranef_add)) {
    stop("'use_cpp_iters_ranef_add' must be FALSE.", call. = FALSE)
  }
  .two_block_batch_iters_ranef_add_r(iters_ranef, chain_i, iters_mean)
}


#' Block 1 batch: update random effects for all chains
#'
#' Pre-\code{73097c2} orchestration: \code{prep_all_chains} then
#' \code{draw_all_chains}, each draw via \code{draw_one_chain} \textrightarrow{}
#' \code{rNormalGLM_reg_group}/\code{rNormal_reg_group} (no per-chain C++ wrapper).
#' Defaults use the R reference for prep/reorder/iters; set \code{use_cpp_*} flags
#' to \code{TRUE} to re-enable piecewise C++ helpers.
#'
#' @param use_cpp_reorder Passed to \code{draw_all_chains} / \code{draw_one_chain}.
#' @param use_cpp_iters Passed to \code{draw_all_chains} / \code{draw_one_chain}.
#' @param use_cpp_mu_all Passed to \code{prep_all_chains}.
#' @param use_cpp_prior_tau2 Passed to \code{prep_all_chains}.
#' @param use_cpp_tau2_row Retained for \code{rGLMM_sweep} API compatibility (unused).
#' @param use_cpp_b_slice Retained for \code{rGLMM_sweep} API compatibility (unused).
#' @param use_cpp_iters_ranef_add Retained for API compatibility (unused).
#' @return List with \code{b} and \code{iters_ranef}.
#' @noRd
.two_block_block1_all_chains <- function(
    n,
    fixef,
    tau2,
    b,
    iters_ranef,
    re_names,
    group_levels,
    design,
    block1_prior,
    family,
    ptypes,
    n_cores = NULL,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE,
    use_cpp_tau2_row = TRUE,
    use_cpp_b_slice = TRUE,
    use_cpp_iters_ranef_add = TRUE,
    use_cpp_reorder = FALSE,
    use_cpp_iters = FALSE,
    use_cpp_mu_all = FALSE,
    use_cpp_prior_tau2 = FALSE
) {
  if (!is.null(n_cores) && as.integer(n_cores[1L]) >= 2L) {
    warning(
      "Chain-parallel Block 1 (n_cores > 1) is not supported; ",
      "running sequential loop.",
      call. = FALSE
    )
  }
  show_bar <- isTRUE(progbar) && n > 1L &&
    (is.null(n_cores) || as.integer(n_cores[1L]) < 2L)

  batch <- list(
    n            = n,
    fixef        = fixef,
    tau2         = tau2,
    re_names     = re_names,
    group_levels = group_levels,
    b            = .two_block_ensure_batch_b_dimnames(b, group_levels, re_names, n),
    iters_ranef  = iters_ranef + 0
  )

  prep <- .two_block_block1_prep_all_chains(
    batch                = batch,
    design               = design,
    block1_prior         = block1_prior,
    ptypes               = ptypes,
    n_cores              = n_cores,
    progbar              = FALSE,
    use_cpp_mu_all       = use_cpp_mu_all,
    use_cpp_prior_tau2   = use_cpp_prior_tau2
  )

  batch <- .two_block_block1_draw_all_chains(
    batch                  = batch,
    prep                   = prep,
    design                 = design,
    family                 = family,
    n_cores                = n_cores,
    progbar                = show_bar,
    progbar_prefix         = progbar_prefix,
    progbar_finish_newline = progbar_finish_newline,
    use_cpp_reorder        = use_cpp_reorder,
    use_cpp_iters          = use_cpp_iters
  )

  list(
    b           = .two_block_ensure_batch_b_dimnames(
      batch$b, group_levels, re_names, n
    ),
    iters_ranef = batch$iters_ranef
  )
}


#' Block 2 batch: update fixed effects / tau2 for all chains
#'
#' Reads batch slice inputs, returns only Block~2 outputs (\code{fixef},
#' \code{tau2}, \code{iters}). Caller assigns \code{fixef}, \code{tau2}, and \code{iters}
#'   into batch state.
#'
#' @param n Number of replicate chains.
#' @param b \code{J x p_re x n} array of random effects (read-only).
#' @param fixef Named list of fixed-effect matrices.
#' @param tau2 \code{n x p_re} matrix of random-effect variances.
#' @param iters \code{n x p_re} matrix of Block~2 iteration counts.
#' @param re_names Character vector of random-effect block names.
#' @param group_levels Character vector of group level labels.
#' @param design Model design list.
#' @param pfamily_list Named list of Block~2 \code{pfamily} objects.
#' @param ptypes Named character vector of \code{pfamily} types.
#' @param use_cpp_block2 lmerb passes \code{TRUE}
#'   (\code{\link{two_block_block2_one_chain_cpp}}). \code{FALSE} is not
#'   implemented.
#' @return List with \code{fixef}, \code{tau2}, and \code{iters}.
#' @noRd
.two_block_block2_all_chains <- function(
    n,
    b,
    fixef,
    tau2,
    iters,
    re_names,
    group_levels,
    design,
    pfamily_list,
    ptypes,
    use_cpp_block2 = TRUE,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE
) {
  show_bar <- isTRUE(progbar) && n > 1L
  if (!isTRUE(use_cpp_block2)) {
    stop("'use_cpp_block2' must be TRUE.", call. = FALSE)
  }
  block2_fn <- two_block_block2_one_chain_cpp

  fixef_out <- fixef
  tau2_out <- tau2
  iters_out <- iters

  for (i in seq_len(n)) {
    if (show_bar) .two_block_progress_bar(i, n, prefix = progbar_prefix)
    chain_out <- block2_fn(
      i            = i,
      b            = b,
      fixef        = fixef_out,
      tau2         = tau2_out,
      iters        = iters_out,
      re_names     = re_names,
      group_levels = group_levels,
      design       = design,
      pfamily_list = pfamily_list,
      ptypes       = ptypes
    )
    fixef_out <- chain_out$fixef
    tau2_out <- chain_out$tau2
    iters_out <- chain_out$iters
  }
  if (show_bar) .two_block_progress_bar_finish(newline = progbar_finish_newline)

  list(
    fixef = fixef_out,
    tau2  = tau2_out,
    iters = iters_out
  )
}

#' Starting tau2 vector from pfamily prior fields (plug-in dispersions)
#' @noRd
.two_block_tau2_start_from_pfamily <- function(pfamily_list, re_names) {
  vapply(re_names, function(k) {
    .two_block_tau2_ref_from_pfamily(pfamily_list[[k]])
  }, numeric(1))
}

#' Pack \code{rGLMM_sweep} outputs (draws, coefficients, \code{mu_all_last})
#' @param n Number of replicate chains.
#' @param fixef Named list of fixed-effect matrices.
#' @param tau2 \code{n x p_re} matrix of random-effect variances.
#' @param b \code{J x p_re x n} array of random effects.
#' @param iters \code{n x p_re} matrix of Block~2 iteration counts.
#' @param iters_ranef Length-\code{n} vector of random-effect iteration counts.
#' @param re_names Character vector of random-effect block names.
#' @param group_levels Character vector of group level labels.
#' @param design Model design list.
#' @param collect_block1 When \code{TRUE}, build long \code{coefficients} table from \code{b}.
#' @return List with \code{fixef_draws}, \code{dispersion_fixef_draws},
#'   \code{iters_fixef_draws}, \code{iters_ranef_draws}, \code{coefficients},
#'   and \code{mu_all_last}.
#' @noRd
.rGLMM_sweep_save <- function(
    n,
    fixef,
    tau2,
    b,
    iters,
    iters_ranef,
    re_names,
    group_levels,
    design,
    collect_block1 = TRUE
) {
  J <- length(group_levels)

  fixef_draws <- lapply(fixef, function(mat) {
    out <- mat
    rownames(out) <- NULL
    out
  })
  names(fixef_draws) <- re_names

  if (isTRUE(collect_block1)) {
    grp_col <- design$group_name
    if (is.null(grp_col) || !nzchar(grp_col)) grp_col <- "group"
    coef_rows <- vector("list", n)
    for (i in seq_len(n)) {
      draw_df <- data.frame(
        draw = rep(i, J),
        stringsAsFactors = FALSE
      )
      draw_df[[grp_col]] <- group_levels
      for (k in re_names) {
        draw_df[[k]] <- b[, k, i]
      }
      coef_rows[[i]] <- draw_df
    }
    coefficients <- do.call(rbind, coef_rows)
    rownames(coefficients) <- NULL
  } else {
    coefficients <- NULL
  }

  fixef_mean <- lapply(fixef, colMeans)
  mu_all_last <- as.matrix(build_mu_all(
    design, fixef_mean, group_levels
  )$mu_all)

  list(
    fixef_draws            = fixef_draws,
    dispersion_fixef_draws = tau2,
    iters_fixef_draws      = iters,
    iters_ranef_draws      = iters_ranef,
    coefficients           = coefficients,
    mu_all_last            = mu_all_last
  )
}
