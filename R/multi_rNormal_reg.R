#' Multi-response Normal regression simulation
#'
#' Runs \code{\link{rNormal_reg}} once per column of \code{y}.  Argument
#' \code{x} may be a single shared design matrix or a \strong{list of matrices}
#' with one entry per response column.  The list path allows each column of
#' \code{y} to have a different number of predictors.
#'
#' @inheritParams rNormal_reg
#' @param n Number of draws per column of \code{y}.
#' @param y Numeric matrix with one column per random-effect component.
#' @param x Shared design matrix, or list of per-column design matrices.
#' @param prior_list List of per-column prior lists (\code{mu}, \code{Sigma}
#'   or \code{P}, optional \code{dispersion}).
#' @return For shared \code{x}, an object of class \code{"mrglmb"}.  For list
#'   \code{x}, a plain named list of \code{rNormal_reg} results.
#' @family simfuncs
#' @seealso \code{\link{rNormal_reg}}
#' @export
multi_rNormal_reg <- function(n,
                              y,
                              x,
                              prior_list,
                              offset = NULL,
                              weights = 1,
                              family = gaussian(),
                              Gridtype = 2,
                              n_envopt = NULL,
                              use_parallel = TRUE,
                              use_opencl = FALSE,
                              verbose = FALSE,
                              progbar = TRUE) {
  call  <- match.call()
  n_draw <- .mrglmb_n_draw(n)

  x_is_list <- is.list(x) && !is.data.frame(x)

  if (x_is_list) {
    y_mat <- as.matrix(y)
    l1    <- ncol(y_mat)
    n_obs <- nrow(y_mat)

    if (l1 < 1L) stop("y must have at least one column.", call. = FALSE)

    coef_names <- colnames(y_mat)
    if (is.null(coef_names) || length(coef_names) != l1) {
      coef_names <- paste0("Y", seq_len(l1))
    }

    if (length(x) != l1) {
      stop(
        "When x is a list, length(x) must equal ncol(y) = ", l1, ".",
        call. = FALSE
      )
    }
    x_list <- lapply(x, as.matrix)
    for (j in seq_len(l1)) {
      if (nrow(x_list[[j]]) != n_obs) {
        stop(
          "nrow(x[[", j, "]]) (", nrow(x_list[[j]]),
          ") must equal nrow(y) (", n_obs, ").",
          call. = FALSE
        )
      }
    }
    p_vec <- vapply(x_list, ncol, integer(1L))

    if (!is.list(prior_list)) {
      stop(
        "prior_list must be a list of length ncol(y) = ", l1, ".",
        call. = FALSE
      )
    }
    if (!is.null(prior_list$mu) || !is.null(prior_list$Sigma)) {
      stop(
        "prior_list must be a list of prior_list objects (one per column ",
        "of y), not a single prior_list with components mu and Sigma.",
        call. = FALSE
      )
    }
    if (length(prior_list) != l1) {
      stop(
        "length(prior_list) must equal ncol(y) = ", l1, ".",
        call. = FALSE
      )
    }
    prior_lists <- lapply(seq_len(l1), function(j) {
      .validate_normal_prior_list(prior_list[[j]], j = j, p = p_vec[j])
    })

    block_results <- vector("list", l1)
    names(block_results) <- coef_names
    for (j in seq_len(l1)) {
      block_results[[j]] <- rNormal_reg(
        n            = n_draw,
        y            = y_mat[, j],
        x            = x_list[[j]],
        prior_list   = prior_lists[[j]],
        offset       = offset,
        weights      = weights,
        family       = family,
        Gridtype     = Gridtype,
        n_envopt     = n_envopt,
        use_parallel = use_parallel,
        use_opencl   = use_opencl,
        verbose      = verbose,
        progbar      = progbar && (j == 1L)
      )
    }
    block_results

  } else {
    inp         <- .mrglmb_check_inputs(y, x, prior_list)
    prior_lists <- .mrglmb_normalize_prior_lists(
      prior_list, inp$l1, inp$p, .validate_normal_prior_list
    )

    block_results <- vector("list", inp$l1)
    for (j in seq_len(inp$l1)) {
      block_results[[j]] <- rNormal_reg(
        n            = n_draw,
        y            = inp$y_mat[, j],
        x            = inp$x,
        prior_list   = prior_lists[[j]],
        offset       = offset,
        weights      = weights,
        family       = family,
        Gridtype     = Gridtype,
        n_envopt     = n_envopt,
        use_parallel = use_parallel,
        use_opencl   = use_opencl,
        verbose      = verbose,
        progbar      = progbar && (j == 1L)
      )
    }
    .mrglmb_assemble(
      block_results,
      inp$coef_names,
      call,
      inp$y_mat,
      inp$x,
      inp$l1,
      inp$p,
      prior_lists,
      inp$pred_names
    )
  }
}

#' @keywords internal
.mrglmb_check_inputs <- function(y, x, spec_list, spec_name = "prior_list") {
  if (missing(spec_list)) {
    stop("'", spec_name, "' is required.", call. = FALSE)
  }
  y_mat <- as.matrix(y)
  x <- as.matrix(x)
  l1 <- ncol(y_mat)
  if (l1 < 1L) {
    stop("y must have at least one column.", call. = FALSE)
  }
  p <- ncol(x)
  if (p < 1L) {
    stop("x must have at least one column.", call. = FALSE)
  }
  if (nrow(x) != nrow(y_mat)) {
    stop("nrow(x) must equal nrow(y).", call. = FALSE)
  }
  coef_names <- colnames(y_mat)
  if (is.null(coef_names) || length(coef_names) != l1) {
    coef_names <- paste0("Y", seq_len(l1))
  }
  pred_names <- colnames(x)
  if (is.null(pred_names) || length(pred_names) != p) {
    pred_names <- paste0("X", seq_len(p))
  }
  list(
    y_mat = y_mat,
    x = x,
    l1 = l1,
    p = p,
    coef_names = coef_names,
    pred_names = pred_names
  )
}

#' @keywords internal
.mrglmb_n_draw <- function(n) {
  n_draw <- if (length(n) > 1L) length(n) else as.integer(n)
  if (!is.finite(n_draw) || n_draw < 1L) {
    stop(
      "'n' must be a positive scalar or a vector whose length defines the number of draws.",
      call. = FALSE
    )
  }
  n_draw
}

#' @keywords internal
.mrglmb_normalize_prior_lists <- function(prior_list, l1, p, validate_fn) {
  if (!is.list(prior_list)) {
    stop(
      "prior_list must be a list of length ncol(y) of per-column prior lists.",
      call. = FALSE
    )
  }
  if (!is.null(prior_list$mu) || !is.null(prior_list$Sigma)) {
    stop(
      "prior_list must be a list of prior_list objects (one per column of y), ",
      "not a single prior_list with components mu and Sigma.",
      call. = FALSE
    )
  }
  if (length(prior_list) != l1) {
    stop("length(prior_list) must equal ncol(y) = ", l1, ".", call. = FALSE)
  }
  lapply(seq_len(l1), function(j) {
    validate_fn(prior_list[[j]], j = j, p = p)
  })
}

#' @keywords internal
.mrglmb_assemble <- function(block_results,
                             coef_names,
                             call,
                             y_mat,
                             x,
                             l1,
                             p,
                             prior_lists,
                             pred_names,
                             pfamily_lists = NULL) {
  outlist <- setNames(block_results, coef_names)
  attr(outlist, "call")       <- call
  attr(outlist, "y")          <- y_mat
  attr(outlist, "x")          <- x
  attr(outlist, "l1")         <- l1
  attr(outlist, "p")          <- p
  attr(outlist, "coef_names") <- coef_names
  attr(outlist, "pred_names") <- pred_names
  if (!is.null(prior_lists)) {
    attr(outlist, "prior_lists") <- prior_lists
  }
  if (!is.null(pfamily_lists)) {
    attr(outlist, "pfamily_lists") <- pfamily_lists
  }
  class(outlist) <- "mrglmb"
  outlist
}

#' @keywords internal
.validate_normal_prior_list <- function(pl, j, p) {
  if (!is.list(pl)) {
    stop("prior_list[[", j, "]] must be a list.", call. = FALSE)
  }
  if (is.null(pl$mu)) {
    stop("prior_list[[", j, "]] must contain 'mu'.", call. = FALSE)
  }
  if (is.null(pl$Sigma) && is.null(pl$P)) {
    stop("prior_list[[", j, "]] must contain 'Sigma' or 'P'.", call. = FALSE)
  }

  mu <- as.numeric(pl$mu)
  if (length(mu) != p) {
    stop(
      "prior_list[[", j, "]]$mu must have length ncol(x) = ", p, ".",
      call. = FALSE
    )
  }

  out <- list(mu = mu)
  if (!is.null(pl$Sigma)) {
    S <- as.matrix(pl$Sigma)
    if (nrow(S) != p || ncol(S) != p) {
      stop(
        "prior_list[[", j, "]]$Sigma must be ", p, " x ", p, ".",
        call. = FALSE
      )
    }
    .check_symmetric_pd(S, label = paste0("prior_list[[", j, "]]$Sigma"))
    out$Sigma <- S
  }
  if (!is.null(pl$P)) {
    P <- as.matrix(pl$P)
    if (nrow(P) != p || ncol(P) != p) {
      stop(
        "prior_list[[", j, "]]$P must be ", p, " x ", p, ".",
        call. = FALSE
      )
    }
    .check_symmetric_pd(P, label = paste0("prior_list[[", j, "]]$P"))
    out$P <- P
  }
  if (!is.null(pl$dispersion)) {
    out$dispersion <- pl$dispersion
  }
  if (!is.null(pl$shape)) {
    out$shape <- pl$shape
  }
  if (!is.null(pl$rate)) {
    out$rate <- pl$rate
  }
  if (!is.null(pl$ddef)) {
    out$ddef <- pl$ddef
  }
  out
}

#' @keywords internal
.check_symmetric_pd <- function(M, label) {
  if (!isSymmetric(M)) {
    stop(label, " must be symmetric.", call. = FALSE)
  }
  tol <- 1e-6
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if (!all(ev >= -tol * abs(ev[1L]))) {
    stop(label, " is not positive definite.", call. = FALSE)
  }
  invisible(TRUE)
}
