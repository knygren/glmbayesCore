#' Conditionally independent block simulation (GLM envelope path)
#'
#' @description
#' Draw from blockwise full conditionals when the posterior factorizes across
#' observation blocks. Each block is sampled with a single draw (\code{n = 1})
#' via [rNormal_reg()] (binomial, Poisson, Gamma, and quasi variants through
#' the \code{rNormalGLM} path). Intended for block Gibbs and similar workflows;
#' convergence of the outer chain is the user's responsibility.
#'
#' @details
#' **Output layout:** \code{coefficients} and \code{coef.mode} are matrices
#' with **rows = blocks** and **columns = regression dimensions** (\code{ncol(x)}).
#'
#' Not wired into [rglmb()] or [rlmb()]; see \code{inst/DESIGN_RGLM_BLOCKS.md}.
#'
#' @param n Number of draws per block; **must be \code{1}** in this implementation.
#' @param y Response vector of length \code{nrow(x)}.
#' @param x Design matrix \code{nrow(x)} by \code{ncol(x)}; same \code{ncol} in every block.
#' @param block Block partition: \code{factor}/integer length \code{l2}, \code{l2_blocks}
#'   counts summing to \code{l2}, or list of row index vectors.
#' @param prior_list Single prior specification recycled to all blocks, or with
#'   \code{mu} as \code{l1} by \code{k} matrix or \code{blocks} sublist.
#' @param prior_lists Optional list of length \code{k} (or \code{1}) of per-block
#'   \code{prior_list} objects.
#' @param offset Optional numeric vector (length \code{1} or \code{length(y)});
#'   partitioned across blocks like \code{y}.
#' @param weights Optional weights; same recycling and blocking as \code{offset}.
#' @param family GLM \code{\link{family}} (not \code{gaussian()}); passed to [rNormal_reg()].
#' @param Gridtype,use_parallel,use_opencl,verbose,progbar Passed to each block's [rNormal_reg()].
#' @param n_envopt Passed to each block; defaults to \code{1} when \code{NULL}.
#' @return A list with class \code{"rNormalGLM_reg_block"} including:
#'   \describe{
#'     \item{coefficients}{Matrix \code{k * p}; row \code{b} is the draw for block \code{b}.}
#'     \item{coef.mode}{Matrix \code{k * p}; posterior mode per block.}
#'     \item{block_info}{Block partition metadata.}
#'     \item{block_results}{List of length \code{k} with each block's [rNormal_reg()] output.}
#'   }
#' @seealso [rNormal_reg], [simfuncs], \code{inst/DESIGN_RGLM_BLOCKS.md}
#' @name simfuncs_block
#' @family simfuncs_block
NULL

#' @rdname simfuncs_block
#' @export
rNormalGLM_reg_block <- function(n,
                                 y,
                                 x,
                                 block,
                                 prior_list = NULL,
                                 prior_lists = NULL,
                                 offset = NULL,
                                 weights = 1,
                                 family = gaussian(),
                                 Gridtype = 2L,
                                 n_envopt = NULL,
                                 use_parallel = TRUE,
                                 use_opencl = FALSE,
                                 verbose = FALSE,
                                 progbar = FALSE) {
  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n != 1L) {
    stop(
      "rNormalGLM_reg_block currently requires n = 1 (one Gibbs draw per block).",
      call. = FALSE
    )
  }

  y <- as.numeric(y)
  x <- as.matrix(x)
  l2 <- length(y)
  l1 <- ncol(x)
  if (nrow(x) != l2) {
    stop("nrow(x) must equal length(y).", call. = FALSE)
  }

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) family <- family()
  if (is.null(family$family)) stop("'family' not recognized.", call. = FALSE)

  if (family$family == "gaussian") {
    stop(
      "rNormalGLM_reg_block is for the GLM envelope path only; ",
      "use block loops with rNormal_reg() for gaussian() or add a Gaussian block helper later.",
      call. = FALSE
    )
  }

  okfamilies <- c(
    "poisson", "binomial", "quasipoisson", "quasibinomial", "Gamma"
  )
  if (!family$family %in% okfamilies) {
    stop(
      "family \"", family$family, "\" is not supported by rNormalGLM_reg_block.",
      call. = FALSE
    )
  }

  offset2 <- offset
  wt <- weights
  if (is.null(offset2)) {
    offset2 <- rep(0, l2)
  } else {
    offset2 <- as.numeric(offset2)
    if (length(offset2) == 1L) offset2 <- rep(offset2, l2)
    if (length(offset2) != l2) {
      stop("length(offset) must be 1 or length(y).", call. = FALSE)
    }
  }
  if (length(wt) == 1L) wt <- rep(wt, l2)
  if (length(wt) != l2) {
    stop("length(weights) must be 1 or length(y).", call. = FALSE)
  }

  block_info <- normalize_block(block, l2)
  k <- block_info$k
  prior_block <- normalize_prior_for_blocks(
    prior_list = prior_list,
    prior_lists = prior_lists,
    block_info = block_info,
    l1 = l1
  )

  coef_draw <- matrix(NA_real_, nrow = k, ncol = l1)
  coef_mode <- matrix(NA_real_, nrow = k, ncol = l1)
  block_results <- vector("list", k)
  dispersion_block <- numeric(k)

  for (b in seq_len(k)) {
    rows_b <- block_info$rows[[b]]
    out_b <- rNormal_reg(
      n = 1L,
      y = y[rows_b],
      x = x[rows_b, , drop = FALSE],
      prior_list = prior_block[[b]],
      offset = offset2[rows_b],
      weights = wt[rows_b],
      family = family,
      Gridtype = Gridtype,
      n_envopt = if (is.null(n_envopt)) 1L else n_envopt,
      use_parallel = use_parallel,
      use_opencl = use_opencl,
      verbose = verbose,
      progbar = progbar
    )
    block_results[[b]] <- out_b
    cb <- out_b$coefficients
    if (is.matrix(cb)) {
      cb <- cb[1L, , drop = TRUE]
    } else {
      cb <- as.numeric(cb)
    }
    if (length(cb) != l1) {
      stop("Block ", b, ": expected ", l1, " coefficients, got ", length(cb), ".",
           call. = FALSE)
    }
    coef_draw[b, ] <- cb
    cm <- out_b$coef.mode
    if (is.matrix(cm)) {
      cm <- as.vector(cm)
    } else {
      cm <- as.numeric(cm)
    }
    if (length(cm) != l1) {
      stop("Block ", b, ": coef.mode length mismatch.", call. = FALSE)
    }
    coef_mode[b, ] <- cm
    dispersion_block[b] <- as.numeric(out_b$dispersion)[1L]
  }

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
    coefficients = coef_draw,
    coef.mode = coef_mode,
    dispersion = dispersion_block,
    n = n,
    k = k,
    l1 = l1,
    l2 = l2,
    block_info = block_info,
    block_results = block_results,
    y = y,
    x = x,
    offset = offset2,
    prior.weights = wt,
    family = family,
    prior_lists = prior_block,
    call = match.call()
  )
  class(outlist) <- c("rNormalGLM_reg_block", "list")
  outlist
}
