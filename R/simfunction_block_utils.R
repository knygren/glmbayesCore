# Internal utilities for conditionally independent block simulation.
# See inst/DESIGN_RGLM_BLOCKS.md.

#' Normalize a row-block partition for BY-style fits
#'
#' @param block Block partition: \code{factor} or integer vector of length
#'   \code{l2}, \code{l2_blocks} counts summing to \code{l2}, or a list of
#'   disjoint row-index vectors covering \code{1:l2}.
#' @param l2 Number of observations (rows) after \code{model.frame}.
#' @return List with \code{k}, \code{ids}, \code{l2_blocks}, \code{starts},
#'   and \code{rows} (per-block row indices).
#' @example inst/examples/Ex_normalize_block.R
#' @export
normalize_block <- function(block, l2) {
  l2 <- as.integer(l2)
  if (length(l2) != 1L || l2 < 1L) {
    stop("'l2' must be a positive integer (length of y).", call. = FALSE)
  }

  if (is.list(block)) {
    if (length(block) < 1L) {
      stop("'block' list must have at least one element.", call. = FALSE)
    }
    rows <- lapply(block, function(idx) {
      idx <- as.integer(idx)
      if (anyNA(idx) || any(idx < 1L) || any(idx > l2)) {
        stop("Row indices in 'block' must be integers in 1:l2.", call. = FALSE)
      }
      unique(idx)
    })
    all_idx <- unlist(rows, use.names = FALSE)
    if (any(duplicated(all_idx))) {
      stop("Row indices in 'block' list must be disjoint.", call. = FALSE)
    }
    if (!identical(sort(all_idx), seq_len(l2))) {
      stop("Row indices in 'block' list must cover exactly 1:l2.", call. = FALSE)
    }
    k <- length(rows)
    ids <- names(block)
    if (is.null(ids) || any(ids == "")) {
      ids <- paste0("block", seq_len(k))
    }
    l2_blocks <- vapply(rows, length, integer(1L))
    starts <- c(1L, cumsum(l2_blocks)[-k] + 1L)
    return(list(
      k = k,
      ids = ids,
      l2_blocks = l2_blocks,
      starts = starts,
      rows = rows
    ))
  }

  # If block is already a factor, preserve its level order before any coercion.
  # as.vector() on a factor returns integer indices, and re-calling factor() on
  # those integers sorts the labels lexicographically, silently destroying the
  # caller-supplied level order.
  if (is.factor(block) && length(block) == l2) {
    blk <- block
    k <- nlevels(blk)
    rows <- split(seq_len(l2), blk)
    ids <- levels(blk)
    l2_blocks <- vapply(rows, length, integer(1L))
    starts <- c(1L, cumsum(l2_blocks)[-k] + 1L)
    return(list(
      k = k,
      ids = ids,
      l2_blocks = l2_blocks,
      starts = starts,
      rows = rows
    ))
  }

  block <- as.vector(block)

  if (length(block) == l2) {
    blk <- factor(block)
    k <- nlevels(blk)
    rows <- split(seq_len(l2), blk)
    ids <- levels(blk)
    l2_blocks <- vapply(rows, length, integer(1L))
    starts <- c(1L, cumsum(l2_blocks)[-k] + 1L)
    return(list(
      k = k,
      ids = ids,
      l2_blocks = l2_blocks,
      starts = starts,
      rows = rows
    ))
  }

  if (length(block) >= 1L && length(block) < l2 && all(block >= 1L)) {
    l2_blocks <- as.integer(block)
    if (sum(l2_blocks) != l2) {
      stop(
        "When 'block' has length k < l2, it is treated as l2_blocks; ",
        "sum(block) must equal length(y) (", l2, ").",
        call. = FALSE
      )
    }
    k <- length(l2_blocks)
    ends <- cumsum(l2_blocks)
    starts <- c(1L, ends[-k] + 1L)
    rows <- lapply(seq_len(k), function(j) seq.int(starts[j], ends[j]))
    ids <- paste0("block", seq_len(k))
    return(list(
      k = k,
      ids = ids,
      l2_blocks = l2_blocks,
      starts = starts,
      rows = rows
    ))
  }

  stop(
    "'block' must be a factor or integer vector of length l2, ",
    "a list of row indices, or an integer vector of l2_blocks counts.",
    call. = FALSE
  )
}

#' @keywords internal
.prior_list_to_P_Sigma <- function(pl) {
  mu <- pl$mu
  if (is.null(mu)) stop("prior_list must contain 'mu'.", call. = FALSE)
  mu <- as.numeric(mu)
  if (!is.null(pl$P)) {
    P <- as.matrix(pl$P)
    if (!isSymmetric(P)) {
      stop("prior precision matrix P must be symmetric.", call. = FALSE)
    }
    Sigma <- tryCatch(solve(P), error = function(e) {
      stop("Could not invert prior precision P: ", conditionMessage(e), call. = FALSE)
    })
    return(list(mu = mu, P = P, Sigma = Sigma))
  }
  if (!is.null(pl$Sigma)) {
    Sigma <- as.matrix(pl$Sigma)
    if (!isSymmetric(Sigma)) {
      stop("prior covariance Sigma must be symmetric.", call. = FALSE)
    }
    R <- chol(Sigma)
    Pinv <- chol2inv(R)
    P <- 0.5 * (Pinv + t(Pinv))
    return(list(mu = mu, P = P, Sigma = Sigma))
  }
  stop("prior_list must contain 'P' or 'Sigma'.", call. = FALSE)
}

#' @keywords internal
.check_P_pd <- function(P, label = "P") {
  tol <- 1e-6
  es <- eigen(P, symmetric = TRUE, only.values = TRUE)
  ev <- es$values
  if (!all(ev >= -tol * abs(ev[1L]))) {
    stop("'", label, "' is not positive definite.", call. = FALSE)
  }
}

#' @keywords internal
normalize_prior_for_blocks <- function(prior_list,
                                       prior_lists,
                                       block_info,
                                       l1) {
  k <- block_info$k
  ids <- block_info$ids
  l1 <- as.integer(l1)

  base_pl <- function(pl) {
    ps <- .prior_list_to_P_Sigma(pl)
    .check_P_pd(ps$P)
    if (length(ps$mu) != l1) {
      stop("length(mu) must equal ncol(x) (", l1, ").", call. = FALSE)
    }
    if (nrow(ps$P) != l1 || ncol(ps$P) != l1) {
      stop("dim(P) or dim(Sigma) must be ", l1, " x ", l1, ".", call. = FALSE)
    }
    out <- list(
      mu = ps$mu,
      Sigma = ps$Sigma,
      P = ps$P
    )
    if (!is.null(pl$dispersion)) out$dispersion <- pl$dispersion
    if ("ddef" %in% names(pl)) out$ddef <- pl$ddef
    out
  }

  if (!is.null(prior_lists)) {
    if (!is.list(prior_lists)) {
      stop("'prior_lists' must be a list.", call. = FALSE)
    }
    if (length(prior_lists) == 1L) {
      one <- base_pl(prior_lists[[1L]])
      return(rep(list(one), k))
    }
    if (length(prior_lists) != k) {
      stop("'prior_lists' must have length 1 or k = ", k, ".", call. = FALSE)
    }
    return(lapply(prior_lists, base_pl))
  }

  if (missing(prior_list) || is.null(prior_list)) {
    stop("Provide 'prior_list' or 'prior_lists'.", call. = FALSE)
  }

  if (!is.null(prior_list$blocks)) {
    bl <- prior_list$blocks
    if (!is.list(bl)) stop("'prior_list$blocks' must be a list.", call. = FALSE)
    if (length(bl) != k) {
      stop("'prior_list$blocks' must have length k = ", k, ".", call. = FALSE)
    }
    if (!is.null(names(bl)) && all(names(bl) != "")) {
      if (!is.null(ids) && all(ids %in% names(bl))) {
        return(lapply(ids, function(nm) base_pl(bl[[nm]])))
      }
    }
    return(lapply(seq_len(k), function(j) base_pl(bl[[j]])))
  }

  mu <- prior_list$mu
  if (is.matrix(mu)) {
    if (nrow(mu) != l1) {
      stop("nrow(prior_list$mu) must equal ncol(x) (", l1, ").", call. = FALSE)
    }
    if (ncol(mu) == 1L) {
      one <- base_pl(c(prior_list, list(mu = as.numeric(mu[, 1L]))))
      return(rep(list(one), k))
    }
    if (ncol(mu) != k) {
      stop("ncol(prior_list$mu) must equal number of blocks k = ", k, ".", call. = FALSE)
    }
    P_list <- NULL
    Sigma_list <- NULL
    if (!is.null(prior_list$P)) {
      if (is.list(prior_list$P) && length(prior_list$P) %in% c(1L, k)) {
        P_list <- prior_list$P
      }
    }
    if (!is.null(prior_list$Sigma)) {
      if (is.list(prior_list$Sigma) && length(prior_list$Sigma) %in% c(1L, k)) {
        Sigma_list <- prior_list$Sigma
      }
    }
    if (is.null(P_list) && is.null(Sigma_list)) {
      if (!is.null(prior_list$P) && is.matrix(prior_list$P)) {
        P_list <- list(prior_list$P)
      } else if (!is.null(prior_list$Sigma) && is.matrix(prior_list$Sigma)) {
        Sigma_list <- list(prior_list$Sigma)
      }
    }
    out <- vector("list", k)
    for (j in seq_len(k)) {
      pl_j <- list(mu = mu[, j])
      if (!is.null(prior_list$dispersion)) {
        disp_j <- prior_list$dispersion
        if (length(disp_j) == k) pl_j$dispersion <- disp_j[j]
        else pl_j$dispersion <- disp_j
      }
      if ("ddef" %in% names(prior_list)) pl_j$ddef <- prior_list$ddef
      if (!is.null(P_list)) {
        pl_j$P <- P_list[[min(j, length(P_list))]]
      } else {
        pl_j$Sigma <- Sigma_list[[min(j, length(Sigma_list))]]
      }
      out[[j]] <- base_pl(pl_j)
    }
    return(out)
  }

  one <- base_pl(prior_list)
  rep(list(one), k)
}

#' @noRd
.prior_payload_for_rNormalGLMBlocks_cpp <- function(prior_block, l1, k) {
  pb1 <- prior_block[[1L]]
  disp_v <- vapply(
    prior_block,
    function(pb) {
      if (is.null(pb$dispersion)) {
        1
      } else {
        as.numeric(pb$dispersion)[1L]
      }
    },
    numeric(1)
  )
  differs <- function(a, b) {
    !isTRUE(all.equal(a$mu, b$mu)) ||
      !isTRUE(all.equal(a$P, b$P)) ||
      !isTRUE(all.equal(a$dispersion, b$dispersion, check.attributes = FALSE))
  }
  prior_by_block <- any(vapply(
    prior_block[-1L],
    function(pb) differs(pb, pb1),
    logical(1L)
  ))

  if (!prior_by_block) {
    return(list(
      mu = matrix(pb1$mu, nrow = l1, ncol = 1L),
      P_blocks = list(pb1$P),
      dispersion = disp_v[1L],
      prior_by_block = FALSE
    ))
  }

  list(
    mu = do.call(cbind, lapply(prior_block, `[[`, "mu")),
    P_blocks = lapply(prior_block, `[[`, "P"),
    dispersion = disp_v,
    prior_by_block = TRUE
  )
}
