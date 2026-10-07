# Return-list assembly for the rLMM* engines.
# .lmebayes_groupef_rank_from_Z is the per-group rank flag the ING
# sweep uses when a group's design is not full rank.

#' Assemble a public \code{rLMM_reg}/\code{rGLMM_reg} return list
#'
#' Enforces group/pop slot names, nests pilot and convergence metadata, and
#' returns components in a fixed glm/glmbayes-style order.
#' @noRd
.lmebayes_assemble_reg_result <- function(
    staged,
    call,
    m_convergence,
    convergence_info,
    pfamily_list,
    dispprior_list,
    family,
    groupef.mode,
    any_non_normal,
    design,
    result_class,
    parent_class,
    draw_engine       = NULL,
    sim_method_used   = NULL,
    icm_info          = NULL,
    pilot_draws       = NULL,
    n_pilot           = NULL,
    m_convergence_pilot = NULL,
    pilot_chisq       = NULL,
    pilot_ub          = NULL,
    tv_tol            = NULL,
    offset            = NULL,
    weights           = 1
) {
  if (!is.list(convergence_info)) {
    convergence_info <- list()
  }
  if (!is.null(draw_engine)) {
    convergence_info$draw_engine <- draw_engine
  }
  if (!is.null(sim_method_used)) {
    convergence_info$sim_method_used <- sim_method_used
  }
  if (!is.null(icm_info)) {
    convergence_info$icm_info <- icm_info
  }
  if (!is.null(pilot_ub)) {
    convergence_info$pilot_ub <- pilot_ub
  }
  if (!is.null(tv_tol)) {
    convergence_info$tv_tol <- tv_tol
  }

  run_pilot <- !is.null(pilot_draws) ||
    (!is.null(n_pilot) && is.finite(n_pilot) && as.integer(n_pilot) > 0L)
  pilot <- if (isTRUE(run_pilot)) {
    list(
      n              = if (!is.null(n_pilot)) as.integer(n_pilot) else NULL,
      m_convergence  = m_convergence_pilot,
      chisq          = pilot_chisq,
      draws          = pilot_draws
    )
  } else {
    NULL
  }

  ## Echo glmbayes-style prior.weights / offset / offset2 (normalized).
  ## Not yet consumed by the mixed-model sampling paths.
  n_obs <- length(design$y)
  prior.weights <- .lmebayes_normalize_weights(weights, n_obs)
  offset_v <- .lmebayes_normalize_offset(offset, n_obs)
  design$weights <- prior.weights
  design$offset <- offset_v

  ## Fixed public order (omit absent optional slots).
  out <- list()
  out$groupef <- staged$groupef
  out$groupef.mode <- groupef.mode
  if (!is.null(staged$groupef.iters)) {
    out$groupef.iters <- staged$groupef.iters
  }
  out$popef <- staged$popef
  out$popef.mode <- staged$popef.mode
  out$popef.init <- staged$popef.init
  if (!is.null(staged$popef.dispersion)) {
    out$popef.dispersion <- staged$popef.dispersion
  }
  if (!is.null(staged$popef.iters)) {
    out$popef.iters <- staged$popef.iters
  }
  if (!is.null(staged[["group.dispersion"]])) {
    out[["group.dispersion"]] <- staged[["group.dispersion"]]
    out[["group.dispersion.mean"]] <- staged[["group.dispersion.mean"]]
    if (is.null(out[["group.dispersion.mean"]])) {
      dr <- out[["group.dispersion"]]
      out[["group.dispersion.mean"]] <- if (is.matrix(dr)) {
        colMeans(dr)
      } else if (length(dr) > 1L) {
        mean(dr)
      } else {
        dr
      }
    }
  } else if (identical(family$family, "gaussian") &&
             !is.null(dispprior_list$dispersion)) {
    ## Fixed observation dispersion: expose on the return even when the
    ## sampler did not draw it.
    out[["group.dispersion"]] <- dispprior_list$dispersion
    out[["group.dispersion.mean"]] <- dispprior_list$dispersion
  }
  if (!is.null(staged[["group.dispersion.iters"]])) {
    out[["group.dispersion.iters"]] <- staged[["group.dispersion.iters"]]
  }
  out$pfamily_list <- pfamily_list
  out$dispprior_list <- dispprior_list
  out$prior.weights <- prior.weights
  out$offset <- offset_v
  out$offset2 <- offset_v
  out$any_non_normal <- any_non_normal
  out$family <- family
  out$design <- design
  out$n <- staged$n
  out$call <- call
  out$m_convergence <- m_convergence
  out$convergence_info <- convergence_info
  if (!is.null(pilot)) {
    out$pilot <- pilot
  }
  if (!is.null(staged$sweep_history)) {
    out$sweep_history <- staged$sweep_history
  }

  class(out) <- unique(c(result_class, parent_class, "list"))
  out
}

#' Per-group full column-rank flag for Block~1 \code{Z_j} (same rule as
#' \code{model_setup()$groupef.rank}).
#' @noRd
.lmebayes_groupef_rank_from_Z <- function(Z, groups, group_levels = NULL) {
  Z <- as.matrix(Z)
  g_chr <- as.character(groups)
  levs <- if (is.null(group_levels)) {
    unique(g_chr)
  } else {
    as.character(group_levels)
  }
  p_re <- ncol(Z)
  stats::setNames(
    vapply(
      levs,
      function(lev) {
        rows <- which(g_chr == lev)
        Z_j  <- Z[rows, , drop = FALSE]
        nrow(Z_j) >= p_re &&
          Matrix::rankMatrix(Z_j, method = "qr")[1L] == p_re
      },
      logical(1L)
    ),
    levs
  )
}
