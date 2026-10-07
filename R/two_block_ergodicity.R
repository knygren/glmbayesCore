## Convergence rate (Remark 8 eigenvalues) for the two-block Gibbs sampler.
##
## Notation follows Nygren (2020), "On the total variation distance between
## multivariate normal densities with applications to two-block Gibbs
## samplers", mapped onto two_block_rNormal_reg(): the paper's x1 (block
## updated SECOND in each sweep) is our Block 2 hyper vector gamma (dim q);
## the paper's x2 (updated FIRST) is our Block 1 random-effect stack b
## (dim J * p_re).  The joint posterior precision blocks are
##
##   P11 = sum_j H_j' P_b H_j + blockdiag_k(V_k^-1)          (q x q)
##   P22 = blockdiag_j( Z_j' W_j Z_j + P_b )                 (J p_re x J p_re)
##   P12 = [ -H_1' P_b | ... | -H_J' P_b ]                   (q x J p_re)
##
## where H_j (p_re x q) maps gamma to the Block 1 prior mean of group j
## (row k of H_j is X_k[j, ] in gamma_k's columns).  The Remark 8 rate is
## the maximal eigenvalue of A = P11^{-1/2} P12 P22^{-1} P21 P11^{-1/2},
## a q x q matrix.  Nothing larger than q x q is ever decomposed here:
## P12 P22^{-1} P21 = sum_j H_j' [P_b B_j^{-1} P_b] H_j is accumulated per
## group with p_re x p_re solves, exploiting the one-nonzero-segment-per-
## component structure of H_j.

#' Shared setup for the Remark-8 two-block rate diagnostics
#' @noRd
.two_block_rate_inputs <- function(x,
                                   group,
                                   x_hyper,
                                   prior_list_block1,
                                   prior_list_block2,
                                   weights = NULL,
                                   family = gaussian(),
                                   group_levels = levels(group)) {
  x <- as.matrix(x)
  l2 <- nrow(x)
  p_re <- ncol(x)
  if (p_re < 1L) {
    stop("'x' must have at least one column.", call. = FALSE)
  }
  re_names <- colnames(x)
  if (is.null(re_names) || length(re_names) != p_re) {
    re_names <- paste0("RE", seq_len(p_re))
    colnames(x) <- re_names
  }

  family <- .two_block_normalize_family(family)

  ## --- group partition (same path as the sampler) ---------------------------
  group_info <- normalize_group(group, l2)
  if (is.null(group_levels)) {
    group_levels <- group_info$ids
  }
  group_levels <- as.character(group_levels)
  idx_map <- match(group_levels, as.character(group_info$ids))
  if (anyNA(idx_map)) {
    stop(
      "group_levels not found in block ids: ",
      paste(group_levels[is.na(idx_map)], collapse = ", "),
      call. = FALSE
    )
  }
  row_idx <- group_info$rows[idx_map]
  J <- length(group_levels)

  ## --- Block 1 prior precision P_b ------------------------------------------
  if (!is.null(prior_list_block1$P)) {
    P_b <- as.matrix(prior_list_block1$P)
  } else if (!is.null(prior_list_block1$Sigma)) {
    S_b <- as.matrix(prior_list_block1$Sigma)
    P_b <- chol2inv(chol(S_b))
  } else {
    stop("prior_list_block1 must contain 'P' or 'Sigma'.", call. = FALSE)
  }
  if (nrow(P_b) != p_re || ncol(P_b) != p_re) {
    stop("dim(prior_list_block1 P/Sigma) must be ", p_re, " x ", p_re, ".",
         call. = FALSE)
  }
  P_b <- 0.5 * (P_b + t(P_b))

  ## --- per-observation working weights ---------------------------------------
  if (is.null(weights)) {
    if (!identical(family$family, "gaussian")) {
      stop(
        "For family \"", family$family, "\" supply explicit 'weights' ",
        "(e.g. IRLS weights at the posterior mode); the Gaussian theory ",
        "then applies as a local approximation only.",
        call. = FALSE
      )
    }
    disp <- prior_list_block1$dispersion
    if (is.null(disp)) {
      stop(
        "prior_list_block1 must contain 'dispersion' to derive weights ",
        "for gaussian().",
        call. = FALSE
      )
    }
    disp <- as.numeric(disp)
    if (length(disp) == 1L) {
      w <- rep(1 / disp, l2)
    } else if (length(disp) == J) {
      ## Fixed per-group dispersion vector, already ordered to match
      ## `group_levels` (see .lmebayes_resolve_dispersion_ranef_fixed_vector()
      ## in lmebayesCore's mixed_rmerb_helpers.R): expand to one weight per
      ## observation via the same group->row mapping used everywhere else in
      ## this function (row_idx, aligned to group_levels).
      w <- numeric(l2)
      for (g in seq_len(J)) {
        w[row_idx[[g]]] <- 1 / disp[g]
      }
    } else {
      stop(
        "length(prior_list_block1$dispersion) must be 1 or the number of ",
        "groups (", J, ").",
        call. = FALSE
      )
    }
    weights_source <- "dispersion"
  } else {
    w <- as.numeric(weights)
    if (length(w) == 1L) w <- rep(w, l2)
    if (length(w) != l2) {
      stop("length(weights) must be 1 or nrow(x).", call. = FALSE)
    }
    if (any(!is.finite(w)) || any(w < 0)) {
      stop("'weights' must be finite and non-negative.", call. = FALSE)
    }
    weights_source <- "user"
  }

  ## --- level-2 designs aligned to group_levels ------------------------------
  if (!is.list(x_hyper) || length(x_hyper) != p_re) {
    stop("'x_hyper' must be a list of length ncol(x) = ", p_re, ".",
         call. = FALSE)
  }
  if (!is.null(names(x_hyper)) && setequal(names(x_hyper), re_names)) {
    x_hyper <- x_hyper[re_names]
  }
  X_hyper <- vector("list", p_re)
  for (k in seq_len(p_re)) {
    X_k <- as.matrix(x_hyper[[k]])
    rn <- rownames(X_k)
    if (is.null(rn)) {
      if (nrow(X_k) != J) {
        stop("nrow(x_hyper[[", k, "]]) must equal length(group_levels).",
             call. = FALSE)
      }
    } else {
      if (!all(group_levels %in% rn)) {
        stop("x_hyper[[", k, "]] is missing rows for some group levels.",
             call. = FALSE)
      }
      X_k <- X_k[group_levels, , drop = FALSE]
    }
    if (is.null(colnames(X_k))) {
      colnames(X_k) <- paste0("X", seq_len(ncol(X_k)))
    }
    X_hyper[[k]] <- X_k
  }
  names(X_hyper) <- re_names
  q_k <- vapply(X_hyper, ncol, integer(1L))
  q <- sum(q_k)
  gamma_cols <- split(seq_len(q), rep(seq_len(p_re), q_k))
  gamma_names <- unlist(lapply(seq_len(p_re), function(k) {
    paste0(re_names[k], "::", colnames(X_hyper[[k]]))
  }), use.names = FALSE)

  ## --- Block 2 prior precisions V_k^-1 ---------------------------------------
  if (!is.list(prior_list_block2) || length(prior_list_block2) != p_re) {
    stop("'prior_list_block2' must be a list of length ncol(x) = ", p_re, ".",
         call. = FALSE)
  }
  if (!is.null(names(prior_list_block2)) &&
      setequal(names(prior_list_block2), re_names)) {
    prior_list_block2 <- prior_list_block2[re_names]
  }
  V_inv <- vector("list", p_re)
  for (k in seq_len(p_re)) {
    pl <- prior_list_block2[[k]]
    if (!is.null(pl$P)) {
      Vk <- as.matrix(pl$P)
    } else if (!is.null(pl$Sigma)) {
      Vk <- chol2inv(chol(as.matrix(pl$Sigma)))
    } else {
      stop("prior_list_block2[[", k, "]] must contain 'P' or 'Sigma'.",
           call. = FALSE)
    }
    if (nrow(Vk) != q_k[k] || ncol(Vk) != q_k[k]) {
      stop("prior_list_block2[[", k, "]] P/Sigma must be ",
           q_k[k], " x ", q_k[k], ".", call. = FALSE)
    }
    V_inv[[k]] <- 0.5 * (Vk + t(Vk))
  }
  names(V_inv) <- re_names

  list(
    x = x,
    re_names = re_names,
    group_levels = group_levels,
    row_idx = row_idx,
    P_b = P_b,
    w = w,
    weights_source = weights_source,
    X_hyper = X_hyper,
    q_k = q_k,
    q = q,
    gamma_cols = gamma_cols,
    gamma_names = gamma_names,
    V_inv = V_inv,
    family = family$family,
    dims = list(J = J, p_re = p_re, q = q, l2 = l2)
  )
}

#' Build Remark-8 joint precision blocks and \eqn{S = P_{12} P_{22}^{-1} P_{21}}
#' @keywords internal
.two_block_S_P11 <- function(inp) {
  p_re <- inp$dims$p_re
  J <- inp$dims$J
  q <- inp$q
  P_b <- inp$P_b
  cols <- inp$gamma_cols
  b_names <- unlist(
    lapply(inp$group_levels, function(g) paste(g, inp$re_names, sep = "::")),
    use.names = FALSE
  )

  S <- matrix(0, q, q)
  G <- matrix(0, q, q)
  P22 <- matrix(0, J * p_re, J * p_re)
  P12 <- matrix(0, q, J * p_re)

  for (j in seq_len(J)) {
    rows <- inp$row_idx[[j]]
    Z_j <- inp$x[rows, , drop = FALSE]
    w_j <- inp$w[rows]
    B_j <- crossprod(Z_j, Z_j * w_j) + P_b
    B_j <- 0.5 * (B_j + t(B_j))
    C_j <- P_b %*% solve(B_j, P_b)
    C_j <- 0.5 * (C_j + t(C_j))

    bc <- (j - 1L) * p_re + seq_len(p_re)
    P22[bc, bc] <- B_j

    x_j <- lapply(seq_len(p_re), function(k) inp$X_hyper[[k]][j, , drop = TRUE])
    H_j <- matrix(0, p_re, q)
    for (k in seq_len(p_re)) {
      H_j[k, cols[[k]]] <- x_j[[k]]
    }
    P12[, bc] <- -t(H_j) %*% P_b

    for (i in seq_len(p_re)) {
      for (k in i:p_re) {
        out_ik <- outer(x_j[[i]], x_j[[k]])
        S[cols[[i]], cols[[k]]] <- S[cols[[i]], cols[[k]]] + C_j[i, k] * out_ik
        G[cols[[i]], cols[[k]]] <- G[cols[[i]], cols[[k]]] + P_b[i, k] * out_ik
        if (k > i) {
          S[cols[[k]], cols[[i]]] <- t(S[cols[[i]], cols[[k]]])
          G[cols[[k]], cols[[i]]] <- t(G[cols[[i]], cols[[k]]])
        }
      }
    }
  }

  P11 <- G
  for (k in seq_len(p_re)) {
    P11[cols[[k]], cols[[k]]] <- P11[cols[[k]], cols[[k]]] + inp$V_inv[[k]]
  }
  S <- 0.5 * (S + t(S))
  P11 <- 0.5 * (P11 + t(P11))
  P22 <- 0.5 * (P22 + t(P22))
  P21 <- t(P12)

  dimnames(S) <- list(inp$gamma_names, inp$gamma_names)
  dimnames(P11) <- list(inp$gamma_names, inp$gamma_names)
  dimnames(P12) <- list(inp$gamma_names, b_names)
  dimnames(P21) <- list(b_names, inp$gamma_names)
  dimnames(P22) <- list(b_names, b_names)

  list(S = S, P11 = P11, P12 = P12, P21 = P21, P22 = P22, P_b = P_b)
}

#' Print joint precision blocks when \code{lambda*} reaches the numerical ceiling
#' @keywords internal
.two_block_print_remark8_precision_stop <- function(
    ev_raw,
    blocks,
    inp,
    S,
    tol = .two_block_lambda_star_one_tol()
) {
  lam_max <- max(ev_raw)
  schur_ev <- eigen(
    0.5 * ((S - blocks$P11) + t(S - blocks$P11)),
    symmetric = TRUE,
    only.values = TRUE
  )$values

  cat("\n=== two_block_rate: lambda* at ceiling (Remark 8 eigenvalue >= 1 - tol) ===\n\n")
  cat(sprintf(
    paste0(
      "  stop tolerance: 1 - tol with tol = %s\n",
      "  raw max eigenvalue of P11^{-1/2} S P11^{-1/2} = %s\n",
      "  full raw spectrum (q = %d) = %s\n",
      "  max eigenvalue of (S - P11) = %s  ",
      "(should be <= 0 for a valid joint precision)\n\n"
    ),
    format(tol, digits = 3),
    format(lam_max, digits = 17),
    inp$q,
    paste(format(sort(ev_raw, decreasing = TRUE), digits = 10), collapse = ", "),
    format(max(schur_ev), digits = 10)
  ))

  cat("Rate-calibration context (working weights embed 1/dispersion):\n")
  cat(sprintf(
    "  family = %s, weights_source = %s, J = %d, p_re = %d, q = %d, n_obs = %d\n",
    inp$family, inp$weights_source, inp$dims$J, inp$dims$p_re, inp$q, inp$dims$l2
  ))
  cat(sprintf(
    "  working weights w: min = %s, max = %s, mean = %s\n",
    format(min(inp$w), digits = 10),
    format(max(inp$w), digits = 10),
    format(mean(inp$w), digits = 10)
  ))
  if (identical(inp$weights_source, "user") && inp$dims$J >= 1L) {
    w_by_group <- vapply(inp$row_idx, function(rows) mean(inp$w[rows]), numeric(1L))
    names(w_by_group) <- inp$group_levels
    cat("  mean working weight 1/sigma^2 plug-in by group:\n")
    print(w_by_group)
  }
  cat("  P_b (random-effect prior precision, p_re x p_re):\n")
  print(blocks$P_b)

  cat("\nJoint precision blocks (Nygren 2020; gamma = x_1 updated second, b = x_2):\n\n")

  cat("P11  (q x q; gamma-gamma block):\n")
  print(blocks$P11)

  cat("\nP12  (q x J*p_re; gamma-to-b cross block; one p_re-wide strip per group):\n")
  print(blocks$P12)

  cat("\nP21  (= t(P12); J*p_re x q):\n")
  print(blocks$P21)

  cat("\nP22  (J*p_re x J*p_re; block-diagonal with blocks B_j = Z_j' W_j Z_j + P_b):\n")
  print(blocks$P22)

  cat("\nS = P12 %*% solve(P22, P21)  (q x q; should match fast-path accumulation):\n")
  print(S)

  cat("\nPer-group B_j diagonal blocks on P22 (each embeds 1/dispersion via W_j):\n")
  p_re <- inp$dims$p_re
  for (j in seq_len(inp$dims$J)) {
    bc <- (j - 1L) * p_re + seq_len(p_re)
    cat(sprintf("  [%s]\n", inp$group_levels[[j]]))
    print(blocks$P22[bc, bc, drop = FALSE])
  }

  invisible(NULL)
}

#' Practical tolerance for treating \code{lambda*} as numerically 1
#'
#' Values within this distance of 1 round to \code{1.0000} in calibration
#' output and imply astronomically large TV sweep counts, but are far above
#' \code{.Machine$double.eps} (e.g. \code{0.99999994} on the 5-school dGamma
#' fixture). \code{two_block_rate()} stops when the raw maximum reaches this
#' tolerance.
#' @noRd
.two_block_lambda_star_one_tol <- function() 1e-7

#' Threshold above which cross-block coupling is treated as very strong.
#' @noRd
.two_block_lambda_star_slow_threshold <- function() 0.95

#' Warn when Remark-8 \code{lambda*} indicates slow two-block mixing.
#' @noRd
.two_block_warn_lambda_star_slow <- function(lambda_star, warn = TRUE, call = NULL) {
  if (!isTRUE(warn)) {
    return(invisible(NULL))
  }
  thr <- .two_block_lambda_star_slow_threshold()
  if (is.finite(lambda_star) && lambda_star > thr) {
    warning(
      "lambda* = ", signif(lambda_star, 4),
      " exceeds ", thr,
      ": correlation between the two Gibbs blocks is very high ",
      "and convergence may be slow.",
      call. = call
    )
  }
  invisible(NULL)
}

#' Generalized eigenvalues of the Remark-8 rate matrix
#'
#' @param S,P11 Remark-8 blocks (see \code{\link{.two_block_S_P11}}).
#' @param blocks,inp Optional; only used to print diagnostics if the
#'   \code{strict} ceiling check fails.
#' @param strict If \code{TRUE} (default), treat \code{lam_max >= 1 - tol} as
#'   an implementation error and \code{stop()} (correct for
#'   \code{\link{two_block_rate}}'s exact base system, where Remark~8
#'   guarantees \code{lambda_star < 1}). If \code{FALSE}, report the raw
#'   spectrum instead. A local, state-dependent Hessian at a single reference
#'   state (for example residuals from a reference fit) is not guaranteed to
#'   satisfy that invariant -- \code{lambda_star >= 1} there is a legitimate
#'   diagnostic finding, not a computation failure. Nonzero residuals only
#'   raise the local rate; they cannot lower it below the \code{u = 0} /
#'   \code{e = 0} case.
#' @return Numeric vector of eigenvalues, descending.
#' @noRd
.two_block_gen_eigen <- function(S, P11, blocks = NULL, inp = NULL, strict = TRUE) {
  q <- nrow(P11)
  R <- chol(P11)
  Rinv <- backsolve(R, diag(q))
  M <- t(Rinv) %*% S %*% Rinv
  M <- 0.5 * (M + t(M))
  ev_raw <- eigen(M, symmetric = TRUE, only.values = TRUE)$values

  lam_max   <- max(ev_raw)
  one_tol   <- .two_block_lambda_star_one_tol()
  lam_stop  <- 1 - one_tol
  lam_clamp <- 1 - .Machine$double.eps
  if (lam_max >= lam_stop) {
    if (!strict) {
      return(sort(pmax(ev_raw, 0), decreasing = TRUE))
    }
    if (!is.null(blocks) && !is.null(inp)) {
      .two_block_print_remark8_precision_stop(
        ev_raw, blocks, inp, S, tol = one_tol
      )
    }
    stop(
      "lambda_star equals 1 (raw Remark-8 eigenvalue ",
      format(lam_max, digits = 17),
      " >= 1 - ", format(one_tol, digits = 3),
      "); see printed P11/P12/P21/P22 above.",
      call. = FALSE
    )
  }

  ev <- pmin(pmax(ev_raw, 0), lam_clamp)
  sort(ev, decreasing = TRUE)
}

#' Two-block Gibbs sampler convergence rate (Remark 8 eigenvalues)
#'
#' Computes the eigenvalues of
#' \eqn{A = P_{11}^{-1/2} P_{12} P_{22}^{-1} P_{21} P_{11}^{-1/2}} for the
#' joint Gaussian posterior of the LMM engines
#' with \code{family = gaussian()} and fixed variance components.  The
#' maximal eigenvalue \eqn{\lambda^*} is the geometric contraction rate of
#' the two-block Gibbs sampler (\insertCite{Nygren2020}{glmbayesCore},
#' Claim 2 / Remark 8 / Corollary 1): the total variation distance between
#' the \eqn{l}-step kernel and the target decays like \eqn{(\lambda^*)^l}.
#'
#' Block convention: the paper's \eqn{x_1} (the block updated second in
#' each sweep) is the Block 2 hyper vector \eqn{\gamma} of dimension
#' \eqn{q = \sum_k q_k}; the paper's \eqn{x_2} is the Block 1
#' random-effect stack.  \eqn{A} is therefore \eqn{q \times q} and is
#' computed without ever forming the \eqn{Jp_{re} \times Jp_{re}}
#' Block 1 precision: per group only \eqn{p_{re} \times p_{re}} solves
#' are needed, followed by a single \eqn{q \times q} symmetric
#' eigendecomposition.
#'
#' For non-Gaussian families the joint posterior is not normal; supplying
#' IRLS-style \code{weights} evaluated at the posterior mode (via the internal
#' \code{two_block_mode_weights()} helper used by \code{rGLMM_reg()} in
#' \pkg{lmebayesCore}) yields a
#' local-Gaussian heuristic rate (no theorem applies).
#'
#' @section From lambda-star to sweeps:
#' \eqn{\lambda^*} alone gives only the geometric proxy
#' \eqn{(\lambda^*)^m \le} \code{tol}, which is loose because it waits for the
#' slowest eigendirection while ignoring that the others converged sweeps
#' earlier.  The sampler instead inverts the exact Theorem 3 bound, which sums
#' closed-form per-eigendirection terms \eqn{d_i^{(l)}} over the full spectrum
#' returned in \code{eigenvalues}:
#' \deqn{\mathrm{TV}\big(K^l(x^{(0)}, \cdot),\, \pi\big) \;\le\;
#'       \sum_{i=1}^{q} d_i^{(l)} \;+\; (\lambda^*)^l \sqrt{D_0 / (2\pi)},}
#' capped at 1, with the eigenvalues taken in ascending order (Lemma 2) and
#' \eqn{D_0} the squared standardized distance of the chain start from the
#' posterior mean.  Because \code{rlmerb} starts every replicate chain at the
#' exact ICM posterior mean, \eqn{D_0 = 0} and the mean term vanishes
#' identically; only the variance-convergence sum has to be certified.
#'
#' The pipeline is therefore
#' \code{two_block_rate()} \eqn{\rightarrow}
#' \code{\link{two_block_tv_bound}()} \eqn{\rightarrow}
#' \code{\link{two_block_l_for_tv}()}, and the sampler adds \strong{one} sweep
#' to the result: the bound applies to the block updated second in each sweep
#' (the hyper vector \eqn{\gamma}), and the stored Block 1 draw lags by a
#' half-step.  Calling \code{print()} on the returned object shows all three
#' answers --- geometric proxy, Corollary 1 envelope, and Theorem 3 --- side
#' by side for a range of tolerances.
#'
#' Cost grows only logarithmically in \eqn{1/}\code{tol}, so certifying an
#' entire \code{n}-draw sample at level \eqn{\alpha} (pass
#' \code{tol = alpha / n}) usually costs a handful of extra sweeps.  When
#' chains cannot start at the posterior mean --- non-Gaussian families, or
#' sampled variance components --- a pilot stage supplies a start whose
#' \eqn{D_0} is certified from the pilot size; see
#' \code{\link{two_block_optimize_pilot_cost}} and
#' \code{vignette("Chapter-B03", package = "lmebayesCore")}.
#'
#' @param x Level-1 RE design matrix \code{Z} (\code{l2 x p_re}), as passed
#'   to \code{\link{rGLMM_sweep}}.
#' @param group Grouping factor or row-group partition of length \code{l2}.
#' @param x_hyper Named list of group-level design matrices \code{X_k}
#'   (\code{J x q_k}), one per column of \code{x}.
#' @param prior_list_block1 Block 1 prior: \code{P} or \code{Sigma}
#'   (\code{p_re x p_re}); \code{dispersion} required for
#'   \code{gaussian()} when \code{weights} is \code{NULL}.
#' @param prior_list_block2 Named list of Block 2 prior lists (each with
#'   \code{P} or \code{Sigma}, \code{q_k x q_k}).
#' @param weights Optional per-observation working weights (length 1 or
#'   \code{l2}).  Default for \code{gaussian()}:
#'   \code{1 / prior_list_block1$dispersion}.  Required for non-Gaussian
#'   families.
#' @param family Response family (default \code{gaussian()}).
#' @param group_levels Character vector defining group order (default
#'   \code{levels(group)}); must match the row order of \code{x_hyper}
#'   when rownames are absent.
#' @param warn_slow If \code{TRUE} (default), issue a \code{warning()} when
#'   \code{lambda_star} exceeds \code{0.95} (strong cross-block coupling and
#'   potentially slow convergence). Set \code{FALSE} for repeated internal
#'   rate evaluations (e.g. pilot upper-bound scans).
#' @return For \code{two_block_rate()} / \code{two_block_rate_from_pfamily_list()},
#'   an object of class \code{"two_block_rate"}: a list with
#'   \code{lambda_star} (the Remark 8 rate), \code{eigenvalues} (full
#'   spectrum \eqn{a_1 \ge \dots \ge a_q} in \eqn{[0,1)}),
#'   \code{m_for_tol} (function: iterations needed so that
#'   \eqn{(\lambda^*)^m \le} \code{tol}), \code{S}
#'   (\eqn{P_{12}P_{22}^{-1}P_{21}}), \code{P11}, \code{dims},
#'   \code{re_names}, \code{gamma_names}, \code{group_levels},
#'   \code{family}, \code{weights_source}, and \code{call}.
#' @references
#'   \insertAllCited{}
#' @family simfuncs
#' @seealso \code{\link{two_block_tv_bound}}, \code{\link{rGLMM_sweep}},
#'   \code{\link{two_block_optimize_pilot_cost}};
#'   \code{vignette("Chapter-B02", package = "lmebayesCore")} for the
#'   derivation and \code{vignette("Chapter-B03", package = "lmebayesCore")}
#'   for the sweep calibration.
#' @example inst/examples/Ex_two_block_rate.R
#' @name two_block_rate
#' @aliases two_block_rate_from_pfamily_list print.two_block_rate
#' @export
two_block_rate <- function(x,
                           group,
                           x_hyper,
                           prior_list_block1,
                           prior_list_block2,
                           weights = NULL,
                           family = gaussian(),
                           group_levels = levels(group),
                           warn_slow = TRUE) {
  cl <- match.call()
  inp <- .two_block_rate_inputs(
    x = x, group = group, x_hyper = x_hyper,
    prior_list_block1 = prior_list_block1,
    prior_list_block2 = prior_list_block2,
    weights = weights, family = family, group_levels = group_levels
  )
  sp <- .two_block_S_P11(inp)
  ev <- .two_block_gen_eigen(sp$S, sp$P11, blocks = sp, inp = inp)
  lambda_star <- ev[1L]
  .two_block_warn_lambda_star_slow(lambda_star, warn = warn_slow, call = cl)

  m_for_tol <- function(tol) {
    if (!is.numeric(tol) || length(tol) != 1L || tol <= 0 || tol >= 1) {
      stop("'tol' must be a single value in (0, 1).", call. = FALSE)
    }
    if (lambda_star <= 0) return(1L)
    as.integer(ceiling(log(tol) / log(lambda_star)))
  }

  structure(
    list(
      lambda_star = lambda_star,
      eigenvalues = ev,
      m_for_tol = m_for_tol,
      S = sp$S,
      P11 = sp$P11,
      P12 = sp$P12,
      P21 = sp$P21,
      P22 = sp$P22,
      dims = inp$dims,
      re_names = inp$re_names,
      gamma_names = inp$gamma_names,
      group_levels = inp$group_levels,
      family = inp$family,
      weights_source = inp$weights_source,
      call = cl
    ),
    class = "two_block_rate"
  )
}

#' @describeIn two_block_rate Print a \code{"two_block_rate"} object
#'   (eigenvalues, \eqn{\lambda^*}, and implied sweeps for common TV
#'   tolerances).
#' @param tols For \code{print.two_block_rate}: numeric vector of TV
#'   tolerances for the implied \code{m_convergence} table.
#' @param ... Ignored by the print methods.
#' @method print two_block_rate
#' @export
print.two_block_rate <- function(x, tols = c(1e-2, 1e-3, 1e-6), ...) {
  d <- x$dims
  cat("Two-block Gibbs convergence rate (Nygren 2020, Remark 8)\n")
  cat(sprintf("  groups J = %d, p_re = %d, q = %d, n_obs = %d  [family: %s]\n",
              d$J, d$p_re, d$q, d$l2, x$family))
  if (identical(x$weights_source, "user")) {
    cat("  (user-supplied weights: local-Gaussian heuristic, not a theorem)\n")
  }
  cat("\nEigenvalues of A (q x q):\n")
  print(signif(x$eigenvalues, 6))
  cat(sprintf("\nlambda* = %.6g\n", x$lambda_star))
  cat("\nImplied sweeps per tolerance:\n")
  cat(sprintf("  %-10s  %18s  %18s  %18s\n",
              "tol", "(lambda*)^m proxy", "Theorem 3 bound", "Corollary 1 bound"))
  for (tol in tols) {
    m_t3 <- tryCatch(two_block_l_for_tv(x, tol, method = "theorem3"),
                     error = function(e) NA_integer_)
    m_c1 <- tryCatch(two_block_l_for_tv(x, tol, method = "corollary1"),
                     error = function(e) NA_integer_)
    cat(sprintf("  %-10g  %18d  %18s  %18s\n",
                tol, x$m_for_tol(tol),
                ifelse(is.na(m_t3), "-", format(m_t3)),
                ifelse(is.na(m_c1), "-", format(m_c1))))
  }
  cat("  (TV bounds assume start at the posterior mean: D0 = 0)\n")
  invisible(x)
}

#' @describeIn two_block_rate Convergence rate from a \code{pfamily_list}
#'   Block~2 spec: thin wrapper that accepts Block~2 priors as \code{pfamily}
#'   objects.  For \code{dNormal} components the fixed \code{dispersion} is
#'   used; for \code{dIndependent_Normal_Gamma} components the conservative
#'   \code{disp_lower} plug-in is used (the rate is then an upper bound over
#'   the truncated tau^2 range).
#' @param pfamily_list Named list of \code{pfamily} objects (one per RE
#'   component), as in \code{\link{rGLMM_sweep}}.
#' @export
two_block_rate_from_pfamily_list <- function(x,
                                              group,
                                              x_hyper,
                                              prior_list_block1,
                                              pfamily_list,
                                              weights = NULL,
                                              family = gaussian(),
                                              group_levels = levels(group),
                                              warn_slow = TRUE) {
  re_names <- names(x_hyper)
  pfamily_list <- .two_block_validate_pfamily_list(pfamily_list, re_names)
  prior_list_block2 <- lapply(pfamily_list, function(pf) {
    pl <- pf$prior_list
    list(
      mu = pl$mu,
      Sigma = pl$Sigma,
      dispersion = if (identical(pf$pfamily, "dNormal")) {
        pl$dispersion
      } else {
        pl$disp_lower
      }
    )
  })
  two_block_rate(
    x = x, group = group, x_hyper = x_hyper,
    prior_list_block1 = prior_list_block1,
    prior_list_block2 = prior_list_block2,
    weights = weights, family = family, group_levels = group_levels,
    warn_slow = warn_slow
  )
}
## Likelihood precision (IRLS/Fisher weights) at the random-effects posterior
## mode, for the local-Gaussian heuristic extension of two_block_rate() to
## non-Gaussian families.
##
## The per-observation weight is the (expected) negative second derivative of
## that observation's log-likelihood with respect to its own linear predictor:
##
##   w_i = wt_i * mu'(eta_i)^2 / (V(mu_i) * dispersion)
##
## evaluated at eta_i = offset_i + Z_i b_mode[j(i), ].  For canonical links
## (gaussian-identity, poisson-log, binomial-logit) this is the exact observed
## Hessian; otherwise it is the expected (Fisher) information.  The weights
## are computed generically from the stats family object (mu.eta, variance),
## which reproduces the family-specific glmbfamfunc()$f7 curvature weights on
## the branches where f7 is correct and avoids its copy-pasted logistic
## weights in the probit / cloglog / Gamma-log branches.  The matrix-level
## precision with respect to b_j is the assembly B_j^lik = Z_j' W_j Z_j.

#' @noRd
two_block_mode_weights <- function(x,
                                   group,
                                   b_mode,
                                   family = gaussian(),
                                   wt = 1,
                                   offset = 0,
                                   dispersion = NULL,
                                   group_levels = levels(group)) {
  cl <- match.call()

  x <- as.matrix(x)
  l2 <- nrow(x)
  p_re <- ncol(x)
  re_names <- colnames(x)
  if (is.null(re_names) || length(re_names) != p_re) {
    re_names <- paste0("RE", seq_len(p_re))
    colnames(x) <- re_names
  }

  family <- .two_block_normalize_family(family)

  ## --- dispersion ------------------------------------------------------------
  needs_disp <- family$family %in% c("gaussian", "Gamma")
  if (is.null(dispersion)) {
    if (needs_disp) {
      stop(
        "'dispersion' is required for family \"", family$family, "\".",
        call. = FALSE
      )
    }
    dispersion <- 1
  }
  dispersion <- as.numeric(dispersion)[1L]
  if (!is.finite(dispersion) || dispersion <= 0) {
    stop("'dispersion' must be a single positive number.", call. = FALSE)
  }

  ## --- group partition (same path as two_block_rate) --------------------------
  group_info <- normalize_group(group, l2)
  if (is.null(group_levels)) {
    group_levels <- group_info$ids
  }
  group_levels <- as.character(group_levels)
  idx_map <- match(group_levels, as.character(group_info$ids))
  if (anyNA(idx_map)) {
    stop(
      "group_levels not found in block ids: ",
      paste(group_levels[is.na(idx_map)], collapse = ", "),
      call. = FALSE
    )
  }
  row_idx <- group_info$rows[idx_map]
  J <- length(group_levels)

  ## --- b_mode aligned to group_levels ----------------------------------------
  b_mode <- as.matrix(b_mode)
  if (ncol(b_mode) != p_re) {
    stop("ncol(b_mode) must equal ncol(x) = ", p_re, ".", call. = FALSE)
  }
  rn <- rownames(b_mode)
  if (!is.null(rn)) {
    if (!all(group_levels %in% rn)) {
      stop("b_mode is missing rows for some group levels.", call. = FALSE)
    }
    b_mode <- b_mode[group_levels, , drop = FALSE]
  } else if (nrow(b_mode) != J) {
    stop("nrow(b_mode) must equal length(group_levels) = ", J, ".",
         call. = FALSE)
  }

  ## --- wt / offset -------------------------------------------------------------
  wt <- as.numeric(wt)
  if (length(wt) == 1L) wt <- rep(wt, l2)
  if (length(wt) != l2) {
    stop("length(wt) must be 1 or nrow(x).", call. = FALSE)
  }
  if (any(!is.finite(wt)) || any(wt < 0)) {
    stop("'wt' must be finite and non-negative.", call. = FALSE)
  }
  offset <- as.numeric(offset)
  if (length(offset) == 1L) offset <- rep(offset, l2)
  if (length(offset) != l2) {
    stop("length(offset) must be 1 or nrow(x).", call. = FALSE)
  }

  ## --- curvature weights at b_mode --------------------------------------------
  ## eta_i = offset_i + z_i' b_mode[j(i), ], evaluated group by group.
  eta <- numeric(l2)
  for (j in seq_len(J)) {
    rows <- row_idx[[j]]
    eta[rows] <- offset[rows] +
      as.vector(x[rows, , drop = FALSE] %*% b_mode[j, ])
  }
  mu <- family$linkinv(eta)
  w <- wt * family$mu.eta(eta)^2 / (family$variance(mu) * dispersion)
  if (any(!is.finite(w))) {
    stop(
      "non-finite curvature weights at b_mode (fitted means on the ",
      "boundary of the parameter space?).",
      call. = FALSE
    )
  }

  ## --- per-group and total likelihood precision blocks ------------------------
  B_lik <- vector("list", J)
  info_total <- matrix(0, p_re, p_re, dimnames = list(re_names, re_names))
  for (j in seq_len(J)) {
    rows <- row_idx[[j]]
    Z_j <- x[rows, , drop = FALSE]
    B_j <- crossprod(Z_j, Z_j * w[rows])
    B_j <- 0.5 * (B_j + t(B_j))
    dimnames(B_j) <- list(re_names, re_names)
    B_lik[[j]] <- B_j
    info_total <- info_total + B_j
  }
  names(B_lik) <- group_levels

  structure(
    list(
      weights = w,
      eta = eta,
      mu = mu,
      B_lik = B_lik,
      info_total = info_total,
      family = family$family,
      link = family$link,
      dispersion = dispersion,
      group_levels = group_levels,
      b_mode = b_mode,
      call = cl
    ),
    class = "two_block_mode_weights"
  )
}

#' @noRd
#' @method print two_block_mode_weights
print.two_block_mode_weights <- function(x, ...) {
  w <- x$weights
  cat("Likelihood precision weights at posterior mode\n")
  cat(sprintf("  family: %s(link = %s),  dispersion = %g\n",
              x$family, x$link, x$dispersion))
  cat(sprintf("  n_obs = %d,  groups J = %d,  p_re = %d\n",
              length(w), length(x$group_levels), ncol(x$b_mode)))
  cat(sprintf("  weights: min = %.4g, median = %.4g, max = %.4g\n",
              min(w), stats::median(w), max(w)))
  if (!identical(x$family, "gaussian")) {
    cat("  (non-Gaussian: local-Gaussian heuristic input for two_block_rate)\n")
  }
  invisible(x)
}
## Total-variation bounds for the two-block Gibbs sampler (Nygren 2020,
## Theorem 3 and Corollary 1), evaluated from the Remark 8 eigenvalue
## spectrum computed by two_block_rate().
##
## With the chain started at the exact posterior mean (as lmerb does via
## lmerb_posterior_mean), the mean term of both bounds is identically zero
## and only the variance-convergence sum remains:
##
##   ||Q^(l) - pi||_TV  <=  sum_{i=1}^n d_i^(l),
##   r_i^(l) = (1 - a_{i-1}^{2l}) / (1 - a_i^{2l}),  a_0 = 0,
##
## with the eigenvalues a_i in ASCENDING order (so that k_i^(l) =
## 1/(1 - a_i^{2l}) is nondecreasing, as Lemma 2 requires).
##
## The generalized error function has the closed form (paper Remark 3)
## erf_n(x) = P(chi_n <= sqrt(2) x) = pchisq(2 x^2, df = n), which makes the
## exact Theorem 3 terms elementary.  The Corollary 1 envelope replaces each
## d_i^(l) with the relaxation of Remarks 5 and 17:
##
##   d_i^(l) <= [a_i^{2l} / (sqrt(1-a_i^2) sqrt(1-a_{i-1}^2))]
##              * sqrt((n+1-i)/2) * c_{n-i},
##   c_m = 2 e^{-m/2} (m/2)^{m/2} / Gamma((m+1)/2),
##
## (the sqrt((n+1-i)/2) factor follows from the Remark 5 derivation; the
## Corollary 1 display omits it, but it is required for the chain
## d_i <= [x1-x2] c_{n-i} with x1-x2 <= [sqrt(r)-sqrt(1/r)] sqrt((n+1-i)/2)).
## The optional mean term uses D0 = (x0-mu)' Sigma11^{-1} (x0-mu):
## exact erf_1 for "theorem3", the linear envelope lambda*^l sqrt(D0/(2 pi))
## for "corollary1".

#' Generalized n-dimensional error function
#'
#' \code{erf_n(x) = pchisq(2 x^2, df = n)} (Nygren 2020, Remark 3; Brown
#' 1963).  \code{n = 1} reduces to the classical error function.
#'
#' @param x Non-negative numeric vector.
#' @param n Dimension (positive integer).
#' @return Numeric vector of probabilities.
#' @keywords internal
.two_block_erfn <- function(x, n) {
  stats::pchisq(2 * x^2, df = n)
}

#' Single-\code{l} TV-bound term (Theorem~3 or Corollary~1)
#' @noRd
.two_block_tv_bound_one <- function(a_asc, l, method, D0, lambda_star) {
  n <- length(a_asc)
  ## a_i^{2l} computed in log space; a = 0 -> 0.
  u <- ifelse(a_asc > 0, exp(2 * l * log(a_asc)), 0)
  u_prev <- c(0, u[-n])

  if (identical(method, "theorem3")) {
    d <- numeric(n)
    for (i in seq_len(n)) {
      m_i <- n + 1L - i
      num <- u[i] - u_prev[i]
      den <- 1 - u[i]
      if (num <= 0 || den <= 0) next
      rm1 <- num / den                  # r - 1
      lr <- log1p(rm1)                  # log(r)
      if (!is.finite(rm1) || rm1 <= 0) next
      t_hi <- m_i * lr * (1 + rm1) / rm1   # m * ln(r) * r/(r-1)
      t_lo <- m_i * lr / rm1               # m * ln(r) / (r-1)
      d[i] <- stats::pchisq(t_hi, df = m_i) - stats::pchisq(t_lo, df = m_i)
    }
    var_term <- sum(d)
    mean_term <- if (D0 > 0) {
      .two_block_erfn(0.5 * lambda_star^l * sqrt(D0) / sqrt(2), 1L)
    } else 0
  } else {
    ## Corollary 1 envelope (Remarks 5 + 17)
    s_i <- sqrt(1 - a_asc^2)
    s_prev <- c(1, s_i[-n])
    d <- numeric(n)
    for (i in seq_len(n)) {
      if (u[i] <= 0) next
      m_i <- n + 1L - i
      m_e <- m_i - 1L                    # erf order exponent (n - i)
      c_m <- if (m_e == 0L) {
        2 / gamma(0.5)                   # = 2/sqrt(pi)
      } else {
        2 * exp(-m_e / 2) * (m_e / 2)^(m_e / 2) / gamma((m_e + 1) / 2)
      }
      d[i] <- u[i] / (s_i[i] * s_prev[i]) * sqrt(m_i / 2) * c_m
    }
    var_term <- sum(d)
    mean_term <- if (D0 > 0) {
      lambda_star^l * sqrt(D0 / (2 * pi))
    } else 0
  }

  min(var_term + mean_term, 1)
}

#' Total-variation bound for the two-block Gibbs sampler
#'
#' Evaluates the bound on the total variation distance between the
#' \eqn{l}-step kernel of the two-block Gibbs sampler and its target
#' (Nygren 2020), from the eigenvalue spectrum computed by
#' \code{\link{two_block_rate}}.
#'
#' \code{method = "theorem3"} evaluates the exact Theorem 3 terms
#' \eqn{d_i^{(l)}} using the closed form
#' \eqn{\mathrm{erf}_n(x) = P(\chi^2_n \le 2x^2)} with
#' \eqn{r_i^{(l)} = (1 - a_{i-1}^{2l})/(1 - a_i^{2l})}.
#' \code{method = "corollary1"} evaluates the looser geometric envelope of
#' Corollary 1 (via Remarks 5 and 17), which decays like
#' \eqn{a_i^{2l}} with explicit constants.
#'
#' When the chain is started at the exact posterior mean (as
#' \code{lmerb} does), \code{D0 = 0} and the mean term of both bounds
#' vanishes identically; only the variance-convergence sum remains.  The
#' returned bound is capped at 1.
#'
#' Note the bound applies to the block updated \emph{second} in each sweep
#' (the Block 2 hyper vector \eqn{\gamma}); the stored Block 1 draw lags by
#' a half-step, so evaluate at \code{l - 1} when calibrating
#' \code{m_convergence} for the random-effect draws.
#'
#' @param rate Object from \code{\link{two_block_rate}}.
#' @param l Integer vector of sweep counts (each \code{>= 1}).
#' @param method \code{"theorem3"} (exact terms) or \code{"corollary1"}
#'   (geometric envelope).
#' @param D0 Optional squared standardized distance of the starting point
#'   from the posterior mean,
#'   \eqn{(x^{(0)}-\mu)^\top \Sigma_{11}^{-1} (x^{(0)}-\mu)}.  Default 0
#'   (start at the posterior mean).
#' @return Numeric vector of TV bounds, one per element of \code{l}, capped
#'   at 1.
#' @references Nygren, K. (2020). \emph{On the total variation distance
#'   between multivariate normal densities with applications to two-block
#'   Gibbs samplers.} Unpublished manuscript.
#' @family simfuncs
#' @seealso \code{\link{two_block_rate}}
#' @name two_block_tv_bound
#' @aliases two_block_l_for_tv
#' @keywords internal
two_block_tv_bound <- function(rate,
                               l,
                               method = c("theorem3", "corollary1"),
                               D0 = 0) {
  if (!inherits(rate, "two_block_rate")) {
    stop("'rate' must be a two_block_rate object.", call. = FALSE)
  }
  method <- match.arg(method)
  l <- as.integer(l)
  if (length(l) < 1L || any(!is.finite(l)) || any(l < 1L)) {
    stop("'l' must contain integers >= 1.", call. = FALSE)
  }
  if (!is.numeric(D0) || length(D0) != 1L || D0 < 0) {
    stop("'D0' must be a single non-negative number.", call. = FALSE)
  }

  a_asc <- sort(rate$eigenvalues)   # Lemma 2 requires ascending order
  vapply(
    l,
    function(li) {
      .two_block_tv_bound_one(a_asc, li, method, D0, rate$lambda_star)
    },
    numeric(1L)
  )
}

#' @describeIn two_block_tv_bound Smallest \code{l} such that
#'   \code{two_block_tv_bound(rate, l, method, D0) <= tol}.  The bound is
#'   decreasing in \code{l}, so a doubling search followed by bisection is
#'   exact.  If the bound does not reach \code{tol} within \code{l_max}
#'   sweeps (near-degenerate \code{rate$lambda_star}, close to 1), a
#'   \code{warning()} is issued (when \code{warn = TRUE}) and the
#'   uncertified geometric proxy \code{l} from \code{rate$m_for_tol(tol)} is
#'   returned (not \code{l_max}).  If \code{lambda_star >= 1 -
#'   .two_block_lambda_star_one_tol()}, \code{\link{two_block_rate}} stops
#'   before sweep calibration proceeds.
#' @param tol Target total-variation tolerance in (0, 1).
#' @param l_max Search cap for the doubling/bisection search only; not used
#'   as the returned sweep count when certification fails.
#' @param warn If \code{TRUE} (default), emit a single warning when the
#'   fallback path is taken.  Set \code{FALSE} for repeated internal calls
#'   (e.g. pilot cost optimization) to avoid warning spam.
#' @return Integer: the required number of sweeps, or an uncertified
#'   fallback (see Details) with a warning when certification fails.
#' @keywords internal
two_block_l_for_tv <- function(rate,
                               tol,
                               method = c("theorem3", "corollary1"),
                               D0 = 0,
                               l_max = 1000000L,
                               warn = TRUE) {
  if (!is.numeric(tol) || length(tol) != 1L || tol <= 0 || tol >= 1) {
    stop("'tol' must be a single value in (0, 1).", call. = FALSE)
  }
  method <- match.arg(method)
  bnd <- function(li) two_block_tv_bound(rate, li, method = method, D0 = D0)

  if (bnd(1L) <= tol) return(1L)
  lo <- 1L
  hi <- 2L
  while (bnd(hi) > tol) {
    lo <- hi
    if (hi >= l_max) {
      fallback <- .two_block_uncertified_l_fallback(rate, tol)
      if (isTRUE(warn)) {
        warning(
          "TV bound does not reach tol = ", tol, " by ", l_max,
          " sweeps (lambda* = ", signif(rate$lambda_star, 4),
          "); using uncertified fallback l = ", fallback, ".",
          call. = FALSE
        )
      }
      return(fallback)
    }
    hi <- min(2L * hi, as.integer(l_max))
  }
  while (hi - lo > 1L) {
    mid <- lo + (hi - lo) %/% 2L
    if (bnd(mid) <= tol) hi <- mid else lo <- mid
  }
  hi
}

#' Coerce an inner-sweep count to a valid positive integer (no upper cap).
#' @noRd
.two_block_cap_inner_sweeps <- function(m) {
  m <- as.integer(m[1L])
  if (length(m) < 1L || !is.finite(m) || m < 1L) return(10L)
  m
}

#' Fallback sweep count when the TV bound cannot certify \code{tol}
#' @noRd
.two_block_uncertified_l_fallback <- function(rate, tol) {
  lam <- rate$lambda_star
  l_proxy <- 10L
  if (is.finite(lam) && lam > 0 && lam < 1) {
    m_proxy <- tryCatch(
      as.integer(rate$m_for_tol(tol)),
      error = function(e) NA_integer_
    )
    if (is.finite(m_proxy) && m_proxy >= 1L) {
      l_proxy <- max(1L, m_proxy - 1L)
    }
  }
  as.integer(max(1L, l_proxy))
}

#' Pilot-stage inner sweeps from mode-gap calibration
#' @noRd
.two_block_m_pilot_from_gap <- function(rate, D_max, c_tol, m_min) {
  m_min <- .two_block_cap_inner_sweeps(m_min)
  if (D_max <= c_tol || !is.finite(rate$lambda_star) || rate$lambda_star <= 0) {
    return(m_min)
  }
  denom <- log(1 / rate$lambda_star)
  if (!is.finite(denom) || denom <= 0) {
    return(m_min)
  }
  gap_m <- as.integer(ceiling(log(D_max / c_tol) / denom))
  if (length(gap_m) < 1L || !is.finite(gap_m) || gap_m < 1L) {
    return(m_min)
  }
  .two_block_cap_inner_sweeps(max(m_min, gap_m))
}
