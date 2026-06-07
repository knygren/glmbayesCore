#' One Gibbs block update via \code{block_rNormalReg}
#'
#' Draw a single blockwise Gaussian posterior sample (\code{n = 1}) and
#' return coefficient draws for two-block Gibbs samplers where each group
#' is its own block (e.g. school-level random effects \eqn{b_j}).
#'
#' This is a thin wrapper around \code{\link{block_rNormalReg}} that always
#' uses \code{n = 1}. When \code{prior_lists} and \code{prior_list} are both
#' omitted, per-block \code{dNormal} priors are built from \code{mu_all},
#' \code{P}, and \code{dispersion}.
#'
#' @param mu_all Numeric matrix \code{l1 x k} or vector of length \code{l1}
#'   (recycled to all blocks): prior means per block (e.g. \eqn{X_\text{hyper}
#'   \gamma}). Required unless \code{prior_lists} or \code{prior_list} is
#'   supplied.
#' @param P Prior precision matrix \code{l1 x l1} (inverse of
#'   \eqn{\Sigma_{\text{ranef}}}). Shared across all blocks when building
#'   per-block priors from \code{mu_all}.
#' @param dispersion Residual variance \eqn{\sigma^2}. Shared across all blocks
#'   when building per-block priors from \code{mu_all}.
#' @param y Response vector of length \code{nrow(x)}.
#' @param x Design matrix \code{nrow(x)} by \code{ncol(x)}.
#' @param block Block partition passed to \code{\link{normalize_block}}.
#' @param prior_list Single prior specification (with \code{mu}, \code{P}/\code{Sigma},
#'   \code{dispersion}). Used when not building priors from \code{mu_all}.
#' @param prior_lists List of length \code{k} (or 1) of per-block prior specifications.
#' @param offset Optional numeric vector (length \code{1} or \code{length(y)}).
#' @param weights Optional weights vector.
#' @param Gridtype Passed to each block's sampler (Armadillo Gridtype).
#' @param seed Optional; passed to \code{\link{set.seed}} before sampling.
#' @param coef_cols Column indices of \code{coefficients} to return as
#'   \code{b_draws} (default \code{NULL} returns all columns).
#' @return A list with:
#'   \describe{
#'     \item{b_draws}{Matrix \code{k x length(coef_cols)} of coefficient draws.}
#'     \item{coefficients,coef.mode}{Matrices from \code{block_rNormalReg}.}
#'     \item{block_rNormalReg}{Full block sampler output.}
#'   }
#' @seealso \code{\link{block_rNormalReg}}, \code{\link{normalize_block}}
#' @example inst/examples/Ex_block_rNormalReg_update.R
#' @rdname block_reg_simfuncs
#' @export
block_rNormalReg_update <- function(mu_all,
                                    P          = NULL,
                                    dispersion = NULL,
                                    y,
                                    x,
                                    block,
                                    prior_list  = NULL,
                                    prior_lists = NULL,
                                    offset  = NULL,
                                    weights = 1,
                                    Gridtype = 2L,
                                    seed = NULL,
                                    coef_cols = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (is.null(prior_lists) && is.null(prior_list)) {
    if (missing(mu_all)) {
      stop(
        "Provide 'mu_all' (with 'P' and 'dispersion'), 'prior_lists', or 'prior_list'.",
        call. = FALSE
      )
    }
    if (is.null(P)) {
      stop("'P' (prior precision) is required when building per-block priors from 'mu_all'.",
           call. = FALSE)
    }
    if (is.null(dispersion)) {
      stop("'dispersion' (residual variance) is required when building per-block priors from 'mu_all'.",
           call. = FALSE)
    }
    mu_mat <- as.matrix(mu_all)
    if (ncol(mu_mat) == 1L || is.null(dim(mu_all))) {
      mu_vec <- as.numeric(mu_all)
      prior_list <- list(mu = mu_vec, P = P, dispersion = dispersion, ddef = FALSE)
    } else {
      k <- ncol(mu_mat)
      prior_lists <- lapply(seq_len(k), function(j) {
        list(mu = mu_mat[, j], P = P, dispersion = dispersion, ddef = FALSE)
      })
    }
  }

  out <- block_rNormalReg(
    n           = 1L,
    y           = y,
    x           = x,
    block       = block,
    prior_list  = prior_list,
    prior_lists = prior_lists,
    offset      = offset,
    weights     = weights,
    Gridtype    = as.integer(Gridtype)
  )

  cols <- if (is.null(coef_cols)) seq_len(ncol(out$coefficients)) else coef_cols
  b_draws <- out$coefficients[, cols, drop = FALSE]

  list(
    b_draws          = b_draws,
    coefficients     = out$coefficients,
    coef.mode        = out$coef.mode,
    block_rNormalReg = out
  )
}
