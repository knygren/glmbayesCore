#' Conditionally independent block Gaussian regression simulation (Gibbs)
#'
#' @description
#' Draw from blockwise full conditionals for **Gaussian** responses when the
#' posterior factorizes across observation blocks, via
#' \code{.rNormalRegBlocks_cpp()} (each block calls \code{rNormalReg}).
#' This is the Gaussian counterpart of \code{\link{block_rNormalGLM}}.
#'
#' Typical use is **block Gibbs** (\code{n = 1} per outer iteration, updating
#' the random-effects block).  \code{n > 1} yields iid draws from the product
#' conditional.
#'
#' @details
#' **Output layout:** \code{coefficients} and \code{coef.mode} are matrices
#' with **rows = blocks** and **columns = predictors** (same column set for
#' every block).
#'
#' **Per-block prior mean:** pass \code{prior_list$mu} as an \code{l1 x k}
#' matrix (one column per block) to supply a distinct prior mean per block.
#' This is the standard form for Block 1 of the lmebayes Gibbs sampler, where
#' \eqn{\mu_j = X_{\text{hyper}} \gamma} depends on the current fixed-effects
#' draw.
#'
#' **Dispersion:** must be supplied explicitly via
#' \code{prior_list$dispersion} (the residual variance \eqn{\sigma^2}).
#' Unlike GLM families, \code{gaussian()} has no implicit unit dispersion.
#'
#' @param n Number of iid draws per block (\code{n = 1} typical for Gibbs).
#' @param y Response vector of length \code{nrow(x)}.
#' @param x Design matrix \code{nrow(x)} by \code{ncol(x)}.
#' @param block Block partition: \code{factor}/integer vector of length
#'   \code{l2}, integer counts summing to \code{l2}, or list of row-index
#'   vectors. See \code{\link{normalize_block}}.
#' @param prior_list Prior specification. Must contain \code{mu}, one of
#'   \code{P} or \code{Sigma}, and \code{dispersion}. \code{mu} may be an
#'   \code{l1 x k} matrix for block-specific prior means.
#' @param prior_lists Optional list of length \code{k} (or \code{1}) of
#'   per-block \code{prior_list} objects.
#' @param offset Optional numeric vector (length \code{1} or \code{length(y)});
#'   partitioned across blocks like \code{y}.
#' @param weights Optional weights; same recycling as \code{offset}.
#' @param Gridtype Passed to each block's sampler (Armadillo Gridtype).
#' @return A list with class \code{"block_rNormalReg"} including:
#'   \describe{
#'     \item{coefficients}{Matrix \code{k x p}; row \code{b} is the draw for block \code{b}.}
#'     \item{coef.mode}{Matrix \code{k x p}; posterior mode per block.}
#'     \item{dispersion}{Numeric vector of length \code{k}; residual dispersion per block.}
#'     \item{block_info}{Block partition metadata from \code{\link{normalize_block}}.}
#'     \item{block_results}{List of length \code{k} with each block's sampler output.}
#'   }
#' @seealso \code{\link{block_rNormalGLM}}, \code{\link{rNormal_reg}},
#'   \code{\link{normalize_block}}
#' @example inst/examples/Ex_block_rNormalReg.R
#' @name block_reg_simfuncs
#' @aliases block_rNormalReg
#' @family block_simfuncs
NULL

#' @rdname block_reg_simfuncs
#' @export
block_rNormalReg <- function(n,
                              y,
                              x,
                              block,
                              prior_list  = NULL,
                              prior_lists = NULL,
                              offset  = NULL,
                              weights = 1,
                              Gridtype = 2L) {
  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  y <- as.numeric(y)
  x <- as.matrix(x)
  l2 <- length(y)
  l1 <- ncol(x)
  if (nrow(x) != l2) stop("nrow(x) must equal length(y).", call. = FALSE)

  offset2 <- offset
  wt <- weights
  if (is.null(offset2)) {
    offset2 <- rep(0, l2)
  } else {
    offset2 <- as.numeric(offset2)
    if (length(offset2) == 1L) offset2 <- rep(offset2, l2)
    if (length(offset2) != l2) stop("length(offset) must be 1 or length(y).", call. = FALSE)
  }
  if (length(wt) == 1L) wt <- rep(wt, l2)
  if (length(wt) != l2) stop("length(weights) must be 1 or length(y).", call. = FALSE)

  block_info <- normalize_block(block, l2)
  k <- block_info$k
  prior_block <- normalize_prior_for_blocks(
    prior_list  = prior_list,
    prior_lists = prior_lists,
    block_info  = block_info,
    l1          = l1
  )

  for (j in seq_len(k)) {
    if (is.null(prior_block[[j]]$dispersion)) {
      stop(
        "prior_list must contain 'dispersion' (residual variance sigma^2) for ",
        "block_rNormalReg. Block ", j, " has no dispersion.",
        call. = FALSE
      )
    }
  }

  famfunc   <- glmbfamfunc(gaussian())
  prior_cpp <- .prior_payload_for_rNormalGLMBlocks_cpp(prior_block, l1, k)

  cpp_out <- .rNormalRegBlocks_cpp(
    n             = n,
    y             = y,
    x             = x,
    offset        = offset2,
    wt            = wt,
    dispersion    = prior_cpp$dispersion,
    mu            = prior_cpp$mu,
    P_blocks      = prior_cpp$P_blocks,
    prior_by_block = prior_cpp$prior_by_block,
    row_blocks    = block_info$rows,
    f2            = famfunc$f2,
    f3            = famfunc$f3,
    Gridtype      = as.integer(Gridtype)
  )

  coef_draw     <- cpp_out$coefficients
  coef_mode     <- cpp_out$coef.mode
  disp_block    <- as.numeric(cpp_out$dispersion)
  block_results <- cpp_out$block_results

  cn <- colnames(x)
  if (!is.null(cn)) {
    colnames(coef_draw) <- cn
    colnames(coef_mode) <- cn
  }
  rn <- block_info$ids
  if (!is.null(rn)) {
    rownames(coef_draw) <- rn
    rownames(coef_mode) <- rn
  }

  outlist <- list(
    coefficients  = coef_draw,
    coef.mode     = coef_mode,
    dispersion    = disp_block,
    n             = n,
    k             = k,
    l1            = l1,
    l2            = l2,
    block_info    = block_info,
    block_results = block_results,
    y             = y,
    x             = x,
    offset        = offset2,
    prior.weights = wt,
    prior_lists   = prior_block,
    call          = match.call()
  )
  class(outlist) <- c("block_rNormalReg", "list")
  outlist
}

