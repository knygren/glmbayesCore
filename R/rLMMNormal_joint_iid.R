#' Exact iid draws from the two-block Gaussian posterior (no Gibbs sweeps)
#'
#' @description
#' Exact iid simulation for the known-variance Gaussian LMM, restricted to
#' \code{family = gaussian()} with fixed (known) observation dispersion and every population
#' component \code{\link[glmbayesCore]{dNormal}} (known/fixed variance
#' components) -- the "known \eqn{\tau^2_k}" route documented in
#' \code{inst/README_KNOWN_VCOV_GAUSSIAN.md}. Under these conditions the
#' joint posterior over \eqn{(\gamma, b_1, ..., b_J)} is \emph{exactly}
#' multivariate normal (see \code{\link{lmerb_posterior_mean}} Details), so
#' every stored draw is generated directly from that Gaussian -- no
#' block-conditional Gibbs sweeps, no burn-in, and no residual
#' autocorrelation between draws.
#' Not exported: call via \code{glmbayesCore:::rLMMNormal_joint_iid()}.
#' Production callers should use
#' \code{\link{rLMMNormal_reg_known_vcov_iid}} (or
#' \code{\link{rLMMNormal_reg_known_vcov}} with \code{sim_method = "DEFAULT"}).
#'
#' @details
#' \code{.lmerb_posterior_normal_system()} builds the posterior
#' precision \code{M} of the stacked Block~2 hyperparameter vector
#' \code{gamma_full} and, per group \eqn{j}, the conditional precision
#' \code{post_P_j} of \eqn{b_j \mid \gamma} (both independent of
#' \code{gamma}). \code{.lmerb_posterior_system_cholesky()}
#' Cholesky-factors \code{M} and every \code{post_P_j} \emph{once}
#' (verifying \code{M} is numerically symmetric first -- see that
#' function's Details for when this can fail); each of the \code{n}
#' replicates then costs one small triangular solve for \code{gamma_full}
#' plus one per group for \code{b_j} -- no iteration, tolerance, or
#' non-convergence is possible.
#'
#' @param n Number of iid draws (integer, at least 1).
#' @param y,x,group,x_hyper,pfamily_list
#'   Response, level-1 design, grouping factor, level-2 designs, and
#'   population priors. Every \code{pfamily_list} component must be
#'   \code{dNormal()}. \code{x} must have unique, non-empty
#'   \code{colnames(x)} and \code{group} must be a factor (there are no
#'   separate \code{groupef.names}/\code{group_levels} arguments). The
#'   grouping-column name (\code{group_name}) is resolved from
#'   \code{attr(group, "group_name")} if set, otherwise from \code{group}'s
#'   own variable name via \code{substitute()}.
#' @param prior_list_block1 Block~1 prior: a fixed \code{dispersion} (a
#'   single positive scalar, or a numeric vector of length
#'   \code{length(group_levels)} giving one fixed, known dispersion per
#'   group). Random (\code{dIndependent_Normal_Gamma}) measurement
#'   dispersion is not supported. The Block~1 random-effect prior precision
#'   (formerly a separate \code{P}/\code{Sigma} field) is always derived
#'   internally from \code{pfamily_list}; \code{prior_list_block1} must not
#'   contain \code{P} or \code{Sigma}.
#' @param progbar Show a text progress bar while drawing.
#' @param verbose Currently unused; accepted so callers can pass the same
#'   flag they pass to \code{\link{rLMMNormal_reg_known_vcov}}.
#' @return Object of class \code{"rLMMNormal_joint_iid"}, a list with
#'   components \code{fixef_draws}, \code{coefficients}, \code{fixef_last},
#'   \code{b_last}, \code{mu_all_last}, \code{dispersion_fixef_draws},
#'   \code{iters_fixef_draws}, \code{pfamily_list}, \code{family}, \code{n},
#'   \code{m_convergence} (always \code{1L}), \code{sampling},
#'   \code{fixef_start}, \code{groupef.names}, \code{group_levels},
#'   \code{group_name}, \code{call}, plus \code{fixef_mean} -- the exact
#'   posterior mean that every draw is centered on (what
#'   \code{\link{lmerb_posterior_mean}} would return as \code{fixef}) -- and
#'   \code{b_mean} (\code{J x p_re}, dimnamed by \code{group_levels}/
#'   \code{re_names}) -- the exact posterior mean of \eqn{\beta_j}, the
#'   group-\eqn{b_j} analogue of \code{fixef_mean}. Unlike \code{b_last}
#'   (the random effects from the final stored draw only), \code{b_mean} is
#'   a stable point estimate: it is computed once, directly from
#'   \code{fixef_mean}, with no Monte Carlo noise.
#' @seealso \code{\link{lmerb_posterior_mean}},
#'   \code{\link{rLMMNormal_reg_known_vcov}}
#' @family simfuncs
#' @keywords internal
rLMMNormal_joint_iid <- function(
    n,
    y,
    x,
    group,
    x_hyper,
    prior_list_block1,
    pfamily_list,
    progbar       = TRUE,
    verbose       = FALSE
) {
  cl <- match.call()

  group_name <- .lmebayes_resolve_group_name(
    group, substitute(group), fn_name = "rLMMNormal_joint_iid"
  )

  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }

  y <- as.vector(y)
  x <- as.matrix(x)
  l2 <- nrow(x)
  if (length(y) != l2) {
    stop("length(y) must equal nrow(x).", call. = FALSE)
  }

  re_names <- colnames(x)
  if (is.null(re_names) || length(re_names) != ncol(x) || anyNA(re_names) ||
      any(!nzchar(re_names)) || anyDuplicated(re_names)) {
    stop(
      "'x' must have unique, non-empty column names (colnames(x)); ",
      "there is no 'groupef.names' argument to override this.",
      call. = FALSE
    )
  }
  p_re <- length(re_names)

  if (!is.factor(group)) {
    stop(
      "'group' must be a factor (wrap with factor(group, levels = ...) ",
      "to control level order or supply a fixed superset of levels); ",
      "there is no 'group_levels' argument to override this.",
      call. = FALSE
    )
  }
  group_levels <- levels(group)
  if (length(group_levels) < 1L) {
    stop("'group' must have at least one level.", call. = FALSE)
  }
  J <- length(group_levels)

  if (!is.list(x_hyper) || is.data.frame(x_hyper)) {
    stop("'x_hyper' must be a list of design matrices.", call. = FALSE)
  }
  if (length(x_hyper) != length(re_names)) {
    stop(
      "length(x_hyper) must equal ncol(x) = ", length(re_names), ".",
      call. = FALSE
    )
  }
  if (!setequal(names(x_hyper), re_names)) {
    stop(
      "names(x_hyper) must match colnames(x): ",
      paste(re_names, collapse = ", "), ".", call. = FALSE
    )
  }
  x_hyper <- x_hyper[re_names]
  x_hyper <- lapply(x_hyper, as.matrix)

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, re_names,
    J = J
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (!pf_summary$all_dNormal) {
    stop(
      "rLMMNormal_joint_iid(): all Block~2 components must be dNormal() ",
      "(known/fixed variance components). Use ",
      "rLMMNormal_reg_estimated_vcov() when a component is ",
      "dIndependent_Normal_Gamma().",
      call. = FALSE
    )
  }

  if (!is.null(prior_list_block1$P) || !is.null(prior_list_block1$Sigma)) {
    stop(
      "'prior_list_block1' must not contain 'P'/'Sigma'; the Block~1 ",
      "random-effect prior precision is derived internally from ",
      "'pfamily_list'.",
      call. = FALSE
    )
  }
  prior_list_block1$P <- .rLMM_P_from_pfamily_list(pfamily_list, re_names)

  .two_block_validate_block1_prior(prior_list_block1, family = gaussian())

  design <- list(
    y             = y,
    D             = x,
    group         = factor(group, levels = group_levels),
    W             = x_hyper,
    groupef.names = re_names,
    group_name    = group_name
  )

  measurement_prior_list <- .two_block_measurement_prior_list(
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
    re_names          = re_names,
    x_hyper           = x_hyper,
    family            = gaussian()
  )

  system   <- .lmerb_posterior_normal_system(design, measurement_prior_list)
  chol_sys <- .lmerb_posterior_system_cholesky(system)

  gamma_mean <- as.vector(solve(system$M, system$v))
  fixef_mean <- stats::setNames(
    lapply(re_names, function(k) {
      g <- gamma_mean[system$idx[[k]]]
      names(g) <- colnames(x_hyper[[k]])
      g
    }),
    re_names
  )

  ## Exact posterior mean of beta_j (b_mean), the group-b_j analogue of
  ## fixef_mean/gamma_mean above: since (gamma, b_1, ..., b_J) is jointly
  ## Gaussian, E[b_j | data] = E[b_j | gamma = gamma_mean, data] (the
  ## conditional-mean map is linear in gamma), so this reuses the same
  ## per-group conjugate-mean formula as the sampling loop below, evaluated
  ## at fixef_mean instead of a random gamma draw -- no Monte Carlo needed.
  mu_all_mean <- as.matrix(
    build_mu_all(design, fixef_mean, group_levels = group_levels)$mu_all
  )
  b_mean <- matrix(
    NA_real_, nrow = J, ncol = p_re,
    dimnames = list(group_levels, re_names)
  )
  for (jj in seq_len(J)) {
    lev          <- group_levels[jj]
    mu_j         <- mu_all_mean[, jj]
    post_v_j     <- system$Zty_scaled[[lev]] + system$P_b %*% mu_j
    b_mean[jj, ] <- as.numeric(solve(system$post_P_j_list[[lev]], post_v_j))
  }

  tau2 <- vapply(re_names, function(k) {
    as.numeric(measurement_prior_list$pop.prior_list[[k]]$dispersion)
  }, numeric(1L))

  fixef_draws <- stats::setNames(
    lapply(re_names, function(k) {
      matrix(
        NA_real_, nrow = n, ncol = length(fixef_mean[[k]]),
        dimnames = list(NULL, names(fixef_mean[[k]]))
      )
    }),
    re_names
  )
  b_arr <- array(NA_real_, dim = c(J, p_re, n))

  if (isTRUE(progbar) && n > 1L) {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  fixef_i  <- fixef_mean
  mu_all_i <- NULL
  for (i in seq_len(n)) {
    z_gamma      <- stats::rnorm(system$P_total)
    gamma_full_i <- gamma_mean + backsolve(chol_sys$R_M, z_gamma)

    fixef_i <- stats::setNames(
      lapply(re_names, function(k) {
        g <- gamma_full_i[system$idx[[k]]]
        names(g) <- colnames(x_hyper[[k]])
        g
      }),
      re_names
    )
    for (k in re_names) {
      fixef_draws[[k]][i, ] <- fixef_i[[k]]
    }

    mu_all_i <- as.matrix(
      build_mu_all(design, fixef_i, group_levels = group_levels)$mu_all
    )
    for (jj in seq_len(J)) {
      lev      <- group_levels[jj]
      mu_j     <- mu_all_i[, jj]
      post_v_j <- system$Zty_scaled[[lev]] + system$P_b %*% mu_j
      mean_j   <- solve(system$post_P_j_list[[lev]], post_v_j)
      z_j      <- stats::rnorm(p_re)
      b_arr[jj, , i] <- as.numeric(mean_j) +
        backsolve(chol_sys$R_j_list[[lev]], z_j)
    }

    if (isTRUE(progbar) && n > 1L) utils::setTxtProgressBar(pb, i)
  }

  fixef_last <- fixef_i
  b_last <- matrix(
    b_arr[, , n], nrow = J, ncol = p_re,
    dimnames = list(group_levels, re_names)
  )
  mu_all_last <- mu_all_i
  dimnames(mu_all_last) <- list(re_names, group_levels)

  dispersion_fixef_draws <- matrix(
    tau2, nrow = n, ncol = p_re, byrow = TRUE,
    dimnames = list(NULL, re_names)
  )
  iters_fixef_draws <- matrix(
    1L, nrow = n, ncol = p_re,
    dimnames = list(NULL, re_names)
  )

  coef_cols <- c("draw", group_name, re_names)
  draw_rows <- vector("list", n)
  for (i in seq_len(n)) {
    draw_df <- data.frame(draw = rep(i, J), stringsAsFactors = FALSE)
    draw_df[[group_name]] <- group_levels
    for (jj in seq_len(p_re)) {
      draw_df[[re_names[jj]]] <- b_arr[, jj, i]
    }
    draw_rows[[i]] <- draw_df
  }
  coefficients <- do.call(rbind, draw_rows)
  rownames(coefficients) <- NULL
  coefficients <- coefficients[, coef_cols, drop = FALSE]

  structure(
    list(
      fixef_draws            = fixef_draws,
      coefficients           = coefficients,
      fixef_last             = fixef_last,
      fixef_mean             = fixef_mean,
      b_last                 = b_last,
      b_mean                 = b_mean,
      mu_all_last            = mu_all_last,
      dispersion_fixef_draws = dispersion_fixef_draws,
      iters_fixef_draws      = iters_fixef_draws,
      pfamily_list           = pfamily_list,
      family                 = gaussian(),
      n                      = n,
      m_convergence          = 1L,
      sampling               = "replicate",
      fixef_start            = fixef_mean,
      groupef.names          = re_names,
      group_levels           = group_levels,
      group_name             = group_name,
      call                   = cl
    ),
    class = c("rLMMNormal_joint_iid", "list")
  )
}
