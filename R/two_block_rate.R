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

#' @keywords internal
.two_block_rate_inputs <- function(x,
                                   block,
                                   x_hyper,
                                   prior_list_block1,
                                   prior_list_block2,
                                   weights = NULL,
                                   family = gaussian(),
                                   group_levels = levels(block)) {
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
  block_info <- normalize_block(block, l2)
  if (is.null(group_levels)) {
    group_levels <- block_info$ids
  }
  group_levels <- as.character(group_levels)
  idx_map <- match(group_levels, as.character(block_info$ids))
  if (anyNA(idx_map)) {
    stop(
      "group_levels not found in block ids: ",
      paste(group_levels[is.na(idx_map)], collapse = ", "),
      call. = FALSE
    )
  }
  row_idx <- block_info$rows[idx_map]
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
    w <- rep(1 / as.numeric(disp)[1L], l2)
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

#' @keywords internal
.two_block_S_P11 <- function(inp) {
  p_re <- inp$dims$p_re
  J <- inp$dims$J
  q <- inp$q
  P_b <- inp$P_b
  cols <- inp$gamma_cols

  S <- matrix(0, q, q)
  G <- matrix(0, q, q)

  for (j in seq_len(J)) {
    rows <- inp$row_idx[[j]]
    Z_j <- inp$x[rows, , drop = FALSE]
    w_j <- inp$w[rows]
    B_j <- crossprod(Z_j, Z_j * w_j) + P_b
    C_j <- P_b %*% solve(B_j, P_b)
    C_j <- 0.5 * (C_j + t(C_j))

    x_j <- lapply(seq_len(p_re), function(k) inp$X_hyper[[k]][j, ])

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
  dimnames(S) <- list(inp$gamma_names, inp$gamma_names)
  dimnames(P11) <- list(inp$gamma_names, inp$gamma_names)

  list(S = S, P11 = P11)
}

#' @keywords internal
.two_block_gen_eigen <- function(S, P11) {
  q <- nrow(P11)
  R <- chol(P11)
  Rinv <- backsolve(R, diag(q))
  M <- t(Rinv) %*% S %*% Rinv
  M <- 0.5 * (M + t(M))
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values

  if (any(ev >= 1)) {
    warning(
      "eigenvalue(s) >= 1 (max = ", format(max(ev), digits = 6),
      "); clamping to < 1. The joint precision may not be positive definite."
    )
  }
  ev <- pmin(pmax(ev, 0), 1 - .Machine$double.eps)
  sort(ev, decreasing = TRUE)
}

#' Two-block Gibbs sampler convergence rate (Remark 8 eigenvalues)
#'
#' Computes the eigenvalues of
#' \eqn{A = P_{11}^{-1/2} P_{12} P_{22}^{-1} P_{21} P_{11}^{-1/2}} for the
#' joint Gaussian posterior targeted by \code{\link{two_block_rNormal_reg}}
#' with \code{family = gaussian()} and fixed variance components.  The
#' maximal eigenvalue \eqn{\lambda^*} is the geometric contraction rate of
#' the two-block Gibbs sampler (Nygren 2020, Claim 2 / Remark 8 /
#' Corollary 1): the total variation distance between the \eqn{l}-step
#' kernel and the target decays like \eqn{(\lambda^*)^l}.
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
#' IRLS-style \code{weights} evaluated at the posterior mode yields a
#' local-Gaussian heuristic rate (no theorem applies).
#'
#' @param x Level-1 RE design matrix \code{Z} (\code{l2 x p_re}), as passed
#'   to \code{\link{two_block_rNormal_reg}}.
#' @param block Grouping factor or block partition of length \code{l2}.
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
#'   \code{levels(block)}); must match the row order of \code{x_hyper}
#'   when rownames are absent.
#' @return Object of class \code{"two_block_rate"}: a list with
#'   \code{lambda_star} (the Remark 8 rate), \code{eigenvalues} (full
#'   spectrum \eqn{a_1 \ge \dots \ge a_q} in \eqn{[0,1)}),
#'   \code{m_for_tol} (function: iterations needed so that
#'   \eqn{(\lambda^*)^m \le} \code{tol}), \code{S}
#'   (\eqn{P_{12}P_{22}^{-1}P_{21}}), \code{P11}, \code{dims},
#'   \code{re_names}, \code{gamma_names}, \code{group_levels},
#'   \code{family}, \code{weights_source}, and \code{call}.
#' @references Nygren, K. (2020). \emph{On the total variation distance
#'   between multivariate normal densities with applications to two-block
#'   Gibbs samplers.} Unpublished manuscript.
#' @family simfuncs
#' @seealso \code{\link{two_block_rNormal_reg}}
#' @export
two_block_rate <- function(x,
                           block,
                           x_hyper,
                           prior_list_block1,
                           prior_list_block2,
                           weights = NULL,
                           family = gaussian(),
                           group_levels = levels(block)) {
  cl <- match.call()
  inp <- .two_block_rate_inputs(
    x = x, block = block, x_hyper = x_hyper,
    prior_list_block1 = prior_list_block1,
    prior_list_block2 = prior_list_block2,
    weights = weights, family = family, group_levels = group_levels
  )
  sp <- .two_block_S_P11(inp)
  ev <- .two_block_gen_eigen(sp$S, sp$P11)
  lambda_star <- ev[1L]

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

#' Print method for two_block_rate objects
#'
#' @param x Object of class \code{"two_block_rate"}.
#' @param tols Numeric vector of TV tolerances for the implied
#'   \code{m_convergence} table.
#' @param ... Ignored.
#' @return \code{x} invisibly.
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
