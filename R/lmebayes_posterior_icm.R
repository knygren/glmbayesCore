#' Exact joint posterior mean of the Gaussian two-block model
#'
#' @description
#' \code{lmerb_posterior_mean()} solves for the exact joint Gaussian
#' posterior mean of the Block~2 hyperparameters and the Block~1 random
#' effects, in closed form.
#'
#' @details
#' \strong{Shared Block~2 update.}
#' For each RE component \eqn{k}:
#' \deqn{
#'   E[\gamma_k \mid b_k] =
#'   \bigl(X_k^\top X_k / \tau^2_k + P_{\gamma_k}\bigr)^{-1}
#'   \bigl(X_k^\top b_k / \tau^2_k + P_{\gamma_k} \mu_{\gamma_k}\bigr)
#' }
#' where \eqn{b_k} is the \eqn{k}-th column of the current Block~1 matrix,
#' \eqn{X_k =} \code{design$W[[k]]},
#' \eqn{\tau^2_k =} \code{dispersion}, and
#' \eqn{P_{\gamma_k} = \Sigma_{\gamma_k}^{-1}} from \code{prior_list}.
#'
#' \strong{Block~1 update.} See \code{\link{lmerb_posterior_mean}}: an exact
#' closed-form Gaussian solve.
#'
#' @param design Design list with \code{y}, \code{D}, \code{group},
#'   \code{W}, and \code{groupef.names}.
#' @param measurement_prior_list List with \code{group.Sigma},
#'   \code{pop.prior_list}, and \code{group.dispersion} (\eqn{\sigma^2}).
#'   Each \code{pop.prior_list[[k]]} must contain \code{mu},
#'   \code{Sigma}, and \code{dispersion}.
#' @param tol Accepted and unused. The mean is an exact solve, so the
#'   return value always has \code{converged = TRUE}, \code{iterations = 1L},
#'   and \code{delta = 0}. Default \code{1e-10}.
#' @param maxit Accepted and unused. Default \code{200L}. See \code{tol}.
#' @return A list with components \code{fixef}, \code{b_mean}, \code{converged},
#'   \code{iterations}, and \code{delta} (the Mahalanobis-distance stopping
#'   statistic described under \code{tol}; always \code{0} for
#'   \code{lmerb_posterior_mean()}).
#' @seealso \code{\link{build_mu_all}}
#' @name lmebayes_posterior_icm
#' @aliases lmerb_posterior_mean
#' @keywords internal
NULL

#' @describeIn lmebayes_posterior_icm Joint posterior \emph{mean} of the
#'   two-block Gaussian model (= joint mode when variance components are
#'   fixed), computed \emph{exactly} in closed form -- no ICM iteration.
#'   Block~1's conditional mean per group \eqn{j},
#'   \deqn{
#'     E[b_j \mid \gamma] =
#'     \bigl(Z_j^\top Z_j / \sigma^2 + P_b\bigr)^{-1}
#'     \bigl(Z_j^\top y_j / \sigma^2 + P_b \,\mu_j(\gamma)\bigr),
#'   }
#'   is affine in \eqn{\gamma} (via \eqn{\mu_j(\gamma)} from
#'   \code{\link{build_mu_all}}), and \eqn{b_j} couples only to the shared
#'   \eqn{\gamma} -- never to another group's \eqn{b_{j'}}.  Substituting this
#'   affine relationship into the Block~2 update (see Details) eliminates all
#'   \eqn{b_j} algebraically, leaving one linear system in \eqn{\gamma} alone
#'   (dimension = total hyperparameter count, independent of the number of
#'   groups \eqn{J}); solving it once gives the exact joint mean, and one
#'   back-substitution pass gives \code{b_mean}.  Forming this system costs
#'   \eqn{O(J)} (one \eqn{p_{re} \times p_{re}} solve per group plus a small
#'   accumulation), never a \eqn{J \times J} or \eqn{J p_{re}}-dimensional
#'   matrix, so it scales to large numbers of groups.
#' @keywords internal
lmerb_posterior_mean <- function(design,
                                 measurement_prior_list,
                                 tol   = 1e-10,
                                 maxit = 200L) {

  .lmerb_validate_design(design)
  .lmerb_validate_measurement_prior_list(measurement_prior_list)

  if (is.null(design$y) || is.null(design$D)) {
    stop("'design' must contain 'y' and 'D'.", call. = FALSE)
  }

  system <- .lmerb_posterior_normal_system(design, measurement_prior_list)
  gamma_full <- as.vector(solve(system$M, system$v))

  fixef <- stats::setNames(
    lapply(system$re_names, function(k) {
      g <- gamma_full[system$idx[[k]]]
      names(g) <- colnames(design$W[[k]])
      g
    }),
    system$re_names
  )

  b <- .lmerb_posterior_b_given_gamma(system, design, fixef)

  list(
    fixef      = fixef,
    b_mean     = b$b,
    converged  = TRUE,
    iterations = 1L,
    delta      = 0
  )
}

#' Exact posterior covariance of the stacked Block~2 hyperparameter vector
#'
#' @description
#' Companion to \code{\link{lmerb_posterior_mean}}: instead of solving for the
#' posterior mean of the stacked Block~2 hyperparameter vector
#' \code{gamma_full}, returns its exact posterior \emph{covariance}
#' \eqn{\Sigma_{11} = M^{-1}}, where \code{M} (the Schur-complement system's
#' posterior precision of \code{gamma_full}) is built by the package-internal
#' \code{.lmerb_posterior_normal_system()}. This is the \eqn{\Sigma_{11}}
#' of Claim~3 in the two-block Gibbs ergodicity reference (see
#' \code{\link{plot_var_convergence}}): the exact target covariance
#' that a two-block Gibbs sampler's cross-chain covariance converges to as the
#' number of inner sweeps grows.
#'
#' Only meaningful when \code{measurement_prior_list$group.dispersion} and
#' \code{measurement_prior_list$group.Sigma} are both \emph{fixed} (not
#' sampled) -- i.e. the same restriction as
#' \code{.lmerb_posterior_normal_system()} itself. For models with
#' estimated dispersion or estimated random-effect variance components, no
#' single exact \eqn{\Sigma_{11}} exists (it would vary by posterior draw);
#' use the empirical cross-chain covariance instead (see
#' \code{\link{plot_var_convergence}}'s fallback).
#'
#' @inheritParams lmerb_posterior_mean
#' @return A \code{P_total x P_total} covariance matrix (\code{P_total} = the
#'   total Block~2 hyperparameter count, summed over RE components), dimnamed
#'   \code{"re_component | covariate"} in the same stacking order as the
#'   package-internal \code{.two_block_snapshot_fixef_cov()}'s \code{coef_index}.
#' @seealso \code{\link{lmerb_posterior_mean}}
#' @keywords internal
lmerb_posterior_covariance <- function(design, measurement_prior_list) {
  .lmerb_validate_design(design)
  .lmerb_validate_measurement_prior_list(measurement_prior_list)

  system <- .lmerb_posterior_normal_system(design, measurement_prior_list)
  Sigma  <- solve(system$M)

  lbl <- unlist(lapply(system$re_names, function(k) {
    paste(k, colnames(design$W[[k]]), sep = " | ")
  }), use.names = FALSE)
  dimnames(Sigma) <- list(lbl, lbl)
  Sigma
}

#' Shared Schur-complement linear system for the exact two-block Gaussian posterior
#'
#' Builds the pieces of the exact joint Gaussian posterior over
#' \code{(gamma, b_1, ..., b_J)} that do not depend on any particular value
#' of \code{gamma}: the \eqn{P_{total} \times P_{total}} linear system
#' \code{M \%*\% gamma_full = v} (\code{M} is the posterior \emph{precision}
#' of the stacked Block~2 hyperparameter vector \code{gamma_full}, \code{v}
#' its precision-weighted mean contribution) and, per group \eqn{j}, the
#' Block~1 conditional precision \code{post_P_j} and \code{Zty_scaled}
#' needed to recover \eqn{E[b_j \mid \gamma]} for any \code{gamma} via
#' \code{.lmerb_posterior_b_given_gamma()}. Derived by eliminating
#' \eqn{b_1, ..., b_J} algebraically from the joint posterior (see
#' \code{\link{lmerb_posterior_mean}} Details) -- cost \code{O(J)}, no
#' \eqn{J \times J} or \eqn{J p_{re}}-dimensional matrix is ever formed.
#'
#' Shared by \code{lmerb_posterior_mean()} (solves \code{M \%*\% gamma = v}
#' once for the exact posterior mean) and
#' \code{\link{rLMMNormal_joint_iid}} (Cholesky-factors \code{M} and each
#' \code{post_P_j} once, then draws \code{n} iid samples from the same
#' system -- no Gibbs sweeps, no burn-in, since the target is exactly
#' Gaussian).
#' @noRd
.lmerb_posterior_normal_system <- function(design, measurement_prior_list) {
  re_names     <- design$groupef.names
  group_levels <- levels(design$group)
  J            <- length(group_levels)
  p_re         <- length(re_names)
  g_chr        <- as.character(design$group)

  sigma2 <- measurement_prior_list$group.dispersion
  if (is.null(sigma2)) {
    stop(
      "'measurement_prior_list' must contain 'group.dispersion' ",
      "for the exact two-block Gaussian posterior.",
      call. = FALSE
    )
  }
  sigma2 <- as.numeric(sigma2)
  if (!(length(sigma2) %in% c(1L, J))) {
    stop(
      "'measurement_prior_list$group.dispersion' must have length 1 or the ",
      "number of groups (", J, ").",
      call. = FALSE
    )
  }
  P_b    <- solve(measurement_prior_list$group.Sigma)

  P_gamma  <- stats::setNames(
    lapply(re_names, function(k) {
      solve(measurement_prior_list$pop.prior_list[[k]]$Sigma)
    }),
    re_names
  )
  mu_gamma <- stats::setNames(
    lapply(re_names, function(k) measurement_prior_list$pop.prior_list[[k]]$mu),
    re_names
  )
  tau2 <- stats::setNames(
    lapply(re_names, function(k) measurement_prior_list$pop.prior_list[[k]]$dispersion),
    re_names
  )

  ZtZ_scaled <- vector("list", J)
  Zty_scaled <- vector("list", J)
  names(ZtZ_scaled) <- names(Zty_scaled) <- group_levels

  for (jj in seq_len(J)) {
    lev  <- group_levels[jj]
    rows <- which(g_chr == lev)
    Z_j  <- design$D[rows, , drop = FALSE]
    y_j  <- design$y[rows]
    ## Per-group fixed dispersion vector (mode = "fixed_vector") is already
    ## ordered to match group_levels by the resolver in mixed_rmerb_helpers.R;
    ## the scalar "fixed" mode broadcasts the same value to every group.
    sigma2_j <- if (length(sigma2) > 1L) sigma2[[jj]] else sigma2
    ZtZ_scaled[[lev]] <- crossprod(Z_j) / sigma2_j
    Zty_scaled[[lev]] <- crossprod(Z_j, y_j) / sigma2_j
  }

  ## --- Exact joint posterior via block (Schur-complement) elimination.
  ##
  ## The target is exactly jointly Gaussian, so there is no need to alternate
  ## Block~1/Block~2 conditional means to a tolerance: eliminate Block~1
  ## algebraically and solve one small linear system for gamma_full instead.
  ##
  ## Block~1 (per group j): b_j = a_j + D_j %*% mu_j(gamma), with
  ##   post_P_j = ZtZ_scaled[[j]] + P_b  (independent of gamma),
  ##   a_j      = solve(post_P_j, Zty_scaled[[j]]),
  ##   D_j      = solve(post_P_j, P_b).
  ## mu_j(gamma) is linear in the stacked hyperparameter vector gamma_full,
  ## so b_j is affine in gamma_full.  Substituting into the Block~2 update
  ## for each RE component k (b_j only couples to the *shared* gamma_full,
  ## never to another group's b_j', so this elimination is exact and costs
  ## O(J) -- no J x J or J*p_re-dimensional matrix is ever formed):
  ##
  ##   post_P_gamma_k %*% gamma_k - sum_k2 (T_{k,k2} / tau2_k) %*% gamma_k2
  ##     = s_k / tau2_k + P_gamma_k %*% mu_gamma_k
  ##
  ## with post_P_gamma_k = crossprod(X_k) / tau2_k + P_gamma_k,
  ##      s_k            = crossprod(X_k, a_k)          (a_k = a_j[k] over j),
  ##      T_{k,k2}       = crossprod(X_k, d_{k,k2} * X_k2), d_{k,k2}[j] = D_j[k, k2].
  ##
  ## Stacking over k gives one P_total x P_total system (P_total = sum of
  ## hyperparameter counts, independent of J); solving it once yields the
  ## *exact* gamma_full, then one back-substitution pass gives b_mean.  No
  ## iteration, tolerance, or non-convergence is possible for this model.
  ## M is exactly the posterior *precision* of gamma_full, and each
  ## post_P_j the posterior precision of b_j | gamma -- this is what makes
  ## the system directly re-usable for iid sampling, not just the mean.
  post_P_j_list <- stats::setNames(vector("list", J), group_levels)
  a_list        <- stats::setNames(vector("list", J), group_levels)
  D_list        <- stats::setNames(vector("list", J), group_levels)

  for (lev in group_levels) {
    post_P_j <- ZtZ_scaled[[lev]] + P_b
    post_P_j_list[[lev]] <- post_P_j
    a_list[[lev]] <- as.vector(solve(post_P_j, Zty_scaled[[lev]]))
    D_list[[lev]] <- solve(post_P_j, P_b)
  }

  p_k     <- vapply(re_names, function(k) ncol(design$W[[k]]), integer(1L))
  P_total <- sum(p_k)
  idx     <- stats::setNames(vector("list", p_re), re_names)
  cursor  <- 0L
  for (k in re_names) {
    idx[[k]] <- seq.int(cursor + 1L, cursor + p_k[[k]])
    cursor   <- cursor + p_k[[k]]
  }

  post_P_gamma_list <- stats::setNames(
    lapply(re_names, function(k) {
      crossprod(design$W[[k]]) / tau2[[k]] + P_gamma[[k]]
    }),
    re_names
  )

  a_mat <- vapply(group_levels, function(lev) a_list[[lev]], numeric(p_re))
  if (p_re == 1L) a_mat <- matrix(a_mat, nrow = 1L)

  M <- matrix(0.0, P_total, P_total)
  v <- numeric(P_total)

  for (ik in seq_len(p_re)) {
    k    <- re_names[[ik]]
    X_k  <- design$W[[k]]
    a_k  <- a_mat[ik, ]
    s_k  <- as.vector(crossprod(X_k, a_k))

    v[idx[[k]]] <- s_k / tau2[[k]] + as.vector(P_gamma[[k]] %*% mu_gamma[[k]])
    M[idx[[k]], idx[[k]]] <- M[idx[[k]], idx[[k]]] + post_P_gamma_list[[k]]

    for (ik2 in seq_len(p_re)) {
      k2    <- re_names[[ik2]]
      X_k2  <- design$W[[k2]]
      d_vec <- vapply(group_levels, function(lev) D_list[[lev]][ik, ik2], numeric(1L))
      T_kk2 <- crossprod(X_k, X_k2 * d_vec)
      M[idx[[k]], idx[[k2]]] <- M[idx[[k]], idx[[k2]]] - T_kk2 / tau2[[k]]
    }
  }

  list(
    M             = M,
    v             = v,
    idx           = idx,
    p_k           = p_k,
    P_total       = P_total,
    post_P_j_list = post_P_j_list,
    Zty_scaled    = Zty_scaled,
    P_b           = P_b,
    re_names      = re_names,
    group_levels  = group_levels,
    J             = J,
    p_re          = p_re
  )
}

#' Block~1 conditional mean given a Block~2 hyperparameter vector
#'
#' Back-substitution step of the exact two-block Gaussian posterior: given
#' any \code{fixef} (a named list of Block~2 hyperparameter vectors, shaped
#' like \code{lmerb_posterior_mean()}'s own \code{fixef} return value --
#' either the exact posterior mean, or one iid draw from
#' \code{\link{rLMMNormal_joint_iid}}), returns \eqn{E[b_j \mid \gamma]} for
#' every group \eqn{j}, using the \code{post_P_j} / \code{Zty_scaled} /
#' \code{P_b} already built by \code{.lmerb_posterior_normal_system()}.
#' @noRd
.lmerb_posterior_b_given_gamma <- function(system, design, fixef) {
  group_levels <- system$group_levels
  J    <- system$J
  p_re <- system$p_re

  mu_all <- as.matrix(
    build_mu_all(design, fixef, group_levels = group_levels)$mu_all
  )
  b <- matrix(
    0.0, nrow = J, ncol = p_re,
    dimnames = list(group_levels, system$re_names)
  )
  for (jj in seq_len(J)) {
    lev      <- group_levels[jj]
    mu_j     <- mu_all[, jj]
    post_v_j <- system$Zty_scaled[[lev]] + system$P_b %*% mu_j
    b[jj, ]  <- solve(system$post_P_j_list[[lev]], post_v_j)
  }

  list(b = b, mu_all = mu_all)
}

#' Cholesky factors of the exact two-block Gaussian posterior, for iid sampling
#'
#' \code{M} (posterior precision of \code{gamma_full}) and each
#' \code{post_P_j} (posterior precision of \eqn{b_j \mid \gamma}) are always
#' \emph{mathematically} symmetric -- but \code{.lmerb_posterior_normal_system()}
#' assembles \code{M} one Block~2 component-row at a time using each
#' component's own \eqn{\tau^2_k} (see that function's Details), which is
#' only guaranteed to reproduce a symmetric \code{M} when
#' \code{measurement_prior_list$group.Sigma} is diagonal \emph{and} its
#' \eqn{k}-th diagonal entry equals component \eqn{k}'s \code{dispersion}
#' (\eqn{\tau^2_k}) -- i.e. \code{group.Sigma = diag(tau2_k)}, exactly how
#' \code{lmerb()}/\code{glmerb()} always construct it (see
#' \code{lmebayes::lmerb} Details), but not a precondition enforced at this
#' matrix level. \code{solve(M, v)} (the posterior \emph{mean}, used by
#' \code{\link{lmerb_posterior_mean}}) is unaffected by this -- it is a
#' correct linear-system solve regardless of \code{M}'s symmetry. Drawing
#' \emph{samples} via \code{chol(M)}, however, is only valid when \code{M}
#' actually is (the precision of) a symmetric covariance, so this check is
#' required before any iid sampling (\code{\link{rLMMNormal_joint_iid}}).
#'
#' @param tol Maximum tolerated relative asymmetry
#'   \code{max(abs(M - t(M))) / max(abs(M))} before erroring. Default
#'   \code{1e-6} (well above floating-point roundoff, which is \code{~1e-16},
#'   but far below the asymmetry produced by a genuinely inconsistent
#'   \code{group.Sigma}/\code{pfamily_list} pairing).
#' @return List with \code{R_M} (\code{chol(M)}, upper triangular) and
#'   \code{R_j_list} (named list of \code{chol(post_P_j)} per group level).
#' @noRd
.lmerb_posterior_system_cholesky <- function(system, tol = 1e-6) {
  M <- system$M
  asym <- max(abs(M - t(M)))
  scale <- max(abs(M))
  if (scale > 0 && asym / scale > tol) {
    stop(
      "rLMMNormal_joint_iid(): the Block~2 posterior precision is not ",
      "symmetric (relative asymmetry ", signif(asym / scale, 3), "). ",
      "Exact iid sampling requires 'group.Sigma' (from 'prior_list_block1') ",
      "to be diagonal with its k-th entry exactly equal to component k's ",
      "dNormal() dispersion (tau2_k) in 'pfamily_list' -- i.e. ",
      "group.Sigma = diag(tau2_k), as lmerb()/glmerb() always construct it. ",
      "Use sim_method = \"TWO_BLOCK_GIBBS\" (or call ",
      "rLMMNormal_reg_known_vcov_two_bg()) ",
      "for a 'group.Sigma' that is correlated or inconsistent with 'pfamily_list'.",
      call. = FALSE
    )
  }
  M_sym <- 0.5 * (M + t(M))

  R_M <- chol(M_sym)
  R_j_list <- stats::setNames(
    lapply(system$group_levels, function(lev) chol(system$post_P_j_list[[lev]])),
    system$group_levels
  )

  list(R_M = R_M, R_j_list = R_j_list)
}

#' @noRd
.lmerb_validate_measurement_prior_list <- function(mpl) {
  if (!is.list(mpl)) {
    stop("'measurement_prior_list' must be a list.", call. = FALSE)
  }
  for (nm in c("group.Sigma", "pop.prior_list")) {
    if (is.null(mpl[[nm]])) {
      stop(
        "'measurement_prior_list' must contain '", nm, "'.",
        call. = FALSE
      )
    }
  }
  invisible(mpl)
}
