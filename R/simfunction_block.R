#' Conditionally independent block simulation (Gibbs / product likelihood)
#'
#' @description
#' Draw from blockwise full conditionals when the posterior factorizes across
#' observation blocks.  \code{\link{block_rNormalReg}} uses
#' \code{.rNormalRegBlocks_cpp()} (Gaussian; each block calls \code{rNormalReg}).
#' \code{\link{block_rNormalGLM}} uses \code{.rNormalGLMBlocks_cpp()} (GLM envelope;
#' each block calls \code{rNormalGLM}).  Typical use is **block Gibbs**
#' (\code{n = 1} per outer step); \code{n > 1} gives iid draws from the product
#' conditional.
#'
#' @details
#' **Output layout:** \code{coefficients} and \code{coef.mode} are matrices with
#' **rows = blocks** and **columns = predictors**.
#'
#' See \code{inst/DESIGN_RGLM_BLOCKS.md}.
#'
#' @param n Number of iid draws per block (\code{n = 1} typical for Gibbs).
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
#' @param family GLM \code{\link{family}} (not \code{gaussian()}).
#' @param Gridtype Passed to each block's sampler (Armadillo Gridtype).
#' @param use_parallel,use_opencl,verbose,progbar Passed to each block's GLM sampler.
#' @param n_envopt Passed to each block; defaults to \code{1} when \code{NULL}.
#' @return A list with class \code{"block_rNormalGLM"} including:
#'   \describe{
#'     \item{coefficients}{Matrix \code{k * p}; row \code{b} is the draw for block \code{b}.}
#'     \item{coef.mode}{Matrix \code{k * p}; posterior mode per block.}
#'     \item{block_info}{Block partition metadata.}
#'     \item{block_results}{List of length \code{k} with each block's sampler output.}
#'   }
#' @seealso \code{\link{rNormal_reg}}, \code{\link{simfunction}},
#'   \code{\link{normalize_block}}, \code{inst/DESIGN_RGLM_BLOCKS.md}
#' @example inst/examples/Ex_block_rNormalGLM.R
#' @name block_simfuncs
#' @aliases block_rNormalGLM block_rNormalReg
#' @family block_simfuncs
NULL

#' @rdname block_simfuncs
#' @export
block_rNormalGLM <- function(n,
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
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
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
      "block_rNormalGLM is for the GLM envelope path only; ",
      "use block loops with rNormal_reg() for gaussian() or add a Gaussian block helper later.",
      call. = FALSE
    )
  }

  okfamilies <- c(
    "poisson", "binomial", "quasipoisson", "quasibinomial", "Gamma"
  )
  if (!family$family %in% okfamilies) {
    stop(
      "family \"", family$family, "\" is not supported by block_rNormalGLM.",
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

  oklinks <- switch(
    family$family,
    poisson = "log",
    quasipoisson = "log",
    binomial = c("logit", "probit", "cloglog"),
    quasibinomial = c("logit", "probit", "cloglog"),
    Gamma = "log",
    character(0)
  )
  if (!family$link %in% oklinks) {
    stop(
      "link \"", family$link, "\" not available for family \"",
      family$family, "\".",
      call. = FALSE
    )
  }

  famfunc <- glmbfamfunc(family)
  prior_cpp <- .prior_payload_for_rNormalGLMBlocks_cpp(prior_block, l1, k)
  n_envopt_use <- if (is.null(n_envopt)) 1L else as.integer(n_envopt)

  cpp_out <- .rNormalGLMBlocks_cpp(
    n = n,
    y = y,
    x = x,
    offset = offset2,
    wt = wt,
    dispersion = prior_cpp$dispersion,
    mu = prior_cpp$mu,
    P_blocks = prior_cpp$P_blocks,
    prior_by_block = prior_cpp$prior_by_block,
    row_blocks = block_info$rows,
    f2 = famfunc$f2,
    f3 = famfunc$f3,
    family = family$family,
    link = family$link,
    Gridtype = as.integer(Gridtype),
    n_envopt = n_envopt_use,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = verbose
  )

  coef_draw <- cpp_out$coefficients
  coef_mode <- cpp_out$coef.mode
  dispersion_block <- as.numeric(cpp_out$dispersion)
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
  class(outlist) <- c("block_rNormalGLM", "list")
  outlist
}

#' @describeIn block_simfuncs Gaussian blockwise full conditionals via
#'   \code{.block_rNormalReg_cpp()} (C++ partition, prior payload, and
#'   \code{rNormalReg()} per block).
#'   This is the Gaussian counterpart of \code{\link{block_rNormalGLM}}.
#'
#' @details
#' **Per-block prior mean:** pass \code{prior_list$mu} as an \code{l1 x k}
#' matrix (one column per block).  Standard for Block~1 of the lmebayes Gibbs
#' sampler, where \eqn{\mu_j = X_{\text{hyper}} \gamma} depends on the current
#' fixed-effects draw.
#'
#' **Dispersion:** must be supplied via \code{prior_list$dispersion} (residual
#' variance \eqn{\sigma^2}).
#'
#' @return A list with class \code{"block_rNormalReg"} including
#'   \code{coefficients}, \code{coef.mode}, \code{dispersion}, and
#'   \code{block_info}.
#' @example inst/examples/Ex_block_rNormalReg.R
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

  famfunc <- glmbfamfunc(gaussian())

  cpp_out <- .block_rNormalReg_cpp(
    n           = n,
    y           = y,
    x           = x,
    block       = block,
    prior_list  = prior_list,
    prior_lists = prior_lists,
    offset      = offset2,
    wt          = wt,
    f2          = famfunc$f2,
    f3          = famfunc$f3,
    Gridtype    = as.integer(Gridtype)
  )

  block_info    <- cpp_out$block_info
  coef_draw     <- cpp_out$coefficients
  coef_mode     <- cpp_out$coef.mode
  disp_block    <- as.numeric(cpp_out$dispersion)
  block_results <- cpp_out$block_results
  prior_block   <- cpp_out$prior_lists
  k             <- cpp_out$k
  l1            <- cpp_out$l1

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

