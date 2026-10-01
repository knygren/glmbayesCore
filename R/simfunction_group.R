#' Conditionally independent row-group simulation (Gibbs / product likelihood)
#'
#' @description
#' Draw from groupwise full conditionals when the posterior factorizes across
#' observation groups.  \code{\link{rNormal_reg_group}} uses the Gaussian
#' \code{rNormalReg()} path per group; \code{\link{rNormalGLM_reg_group}} uses
#' the GLM envelope path per group. Typical use is block Gibbs
#' (\code{n = 1} per outer step); \code{n > 1} gives iid draws from the product
#' conditional.
#'
#' @details
#' **Output layout:** \code{coefficients} and \code{coef.mode} are matrices with
#' **rows = groups** and **columns = predictors**.
#'
#' @param n Number of iid draws per group (\code{n = 1} typical for Gibbs).
#' @param y Response vector of length \code{nrow(x)}.
#' @param x Design matrix \code{nrow(x)} by \code{ncol(x)}; same \code{ncol} in every group.
#' @param group Group partition: \code{factor}/integer length \code{l2}, \code{l2_blocks}
#'   counts summing to \code{l2}, or list of row index vectors.
#' @param prior_list Single prior specification recycled to all groups, or with
#'   \code{mu} as \code{l1} by \code{k} matrix or \code{blocks} sublist.
#' @param prior_lists Optional list of length \code{k} (or \code{1}) of per-group
#'   \code{prior_list} objects.
#' @param offset Optional numeric vector (length \code{1} or \code{length(y)});
#'   partitioned across groups like \code{y}.
#' @param weights Optional weights; same recycling and blocking as \code{offset}.
#' @param family GLM \code{\link{family}} (not \code{gaussian()}).
#' @param Gridtype Passed to each group's sampler (Armadillo Gridtype).
#' @param use_parallel,use_opencl,verbose,progbar Passed to each group's GLM sampler.
#' @param n_envopt Passed to each group; defaults to \code{1} when \code{NULL}.
#' @return For \code{rNormalGLM_reg_group}, a list with class
#'   \code{"rNormalGLM_reg_group"} including \code{coefficients},
#'   \code{coef.mode}, \code{group_info}, and \code{group_results}.
#' @seealso \code{\link{rNormal_reg}}, \code{\link{simfunction}},
#'   \code{\link{normalize_group}}, \code{\link{Prior_SetupGroup}},
#'   \code{\link{pfamily_list}}
#' @example inst/examples/Ex_rNormalGLM_reg_group.R
#' @name simfuncs_group
#' @aliases rNormalGLM_reg_group rNormal_reg_group
NULL

#' @rdname simfuncs_group
#' @export
rNormalGLM_reg_group <- function(n,
                                 y,
                                 x,
                                 group,
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
      "rNormalGLM_reg_group is for the GLM envelope path only; ",
      "use rNormal_reg_group() for gaussian().",
      call. = FALSE
    )
  }

  okfamilies <- c(
    "poisson", "binomial", "quasipoisson", "quasibinomial", "Gamma"
  )
  if (!family$family %in% okfamilies) {
    stop(
      "family \"", family$family, "\" is not supported by rNormalGLM_reg_group.",
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
  n_envopt_use <- if (is.null(n_envopt)) 1L else as.integer(n_envopt)

  cpp_out <- .group_rNormalGLM_cpp(
    n            = n,
    y            = y,
    x            = x,
    group        = group,
    prior_list   = prior_list,
    prior_lists  = prior_lists,
    offset       = offset2,
    wt           = wt,
    f2           = famfunc$f2,
    f3           = famfunc$f3,
    family       = family$family,
    link         = family$link,
    Gridtype     = as.integer(Gridtype),
    n_envopt     = n_envopt_use,
    use_parallel = use_parallel,
    use_opencl   = use_opencl,
    verbose      = verbose
  )

  group_info <- cpp_out$group_info
  k <- cpp_out$k
  prior_block <- cpp_out$prior_lists
  coef_draw <- cpp_out$coefficients
  coef_mode <- cpp_out$coef.mode
  dispersion_block <- as.numeric(cpp_out$dispersion)
  group_results <- cpp_out$group_results

  cn <- colnames(x)
  if (!is.null(cn)) {
    colnames(coef_draw) <- cn
    colnames(coef_mode) <- cn
  }
  rn <- group_info$ids
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
    group_info = group_info,
    group_results = group_results,
    y = y,
    x = x,
    offset = offset2,
    prior.weights = wt,
    family = family,
    prior_lists = prior_block,
    call = match.call()
  )
  class(outlist) <- c("rNormalGLM_reg_group", "list")
  outlist
}

#' @describeIn simfuncs_group Gaussian groupwise full conditionals via
#'   \code{rNormalReg()} per group.
#'
#' @details
#' **Per-group prior mean:** pass \code{prior_list$mu} as an \code{l1 x k}
#' matrix (one column per group).
#'
#' **Dispersion:** must be supplied via \code{prior_list$dispersion} (residual
#' variance \eqn{\sigma^2}).
#'
#' @return A list with class \code{"rNormal_reg_group"} including
#'   \code{coefficients}, \code{coef.mode}, \code{dispersion}, and
#'   \code{group_info}.
#' @example inst/examples/Ex_rNormal_reg_group.R
#' @export
rNormal_reg_group <- function(n,
                              y,
                              x,
                              group,
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

  cpp_out <- .group_rNormalReg_cpp(
    n           = n,
    y           = y,
    x           = x,
    group       = group,
    prior_list  = prior_list,
    prior_lists = prior_lists,
    offset      = offset2,
    wt          = wt,
    f2          = famfunc$f2,
    f3          = famfunc$f3,
    Gridtype    = as.integer(Gridtype)
  )

  group_info    <- cpp_out$group_info
  coef_draw     <- cpp_out$coefficients
  coef_mode     <- cpp_out$coef.mode
  disp_block    <- as.numeric(cpp_out$dispersion)
  group_results <- cpp_out$group_results
  prior_block   <- cpp_out$prior_lists
  k             <- cpp_out$k
  l1            <- cpp_out$l1

  cn <- colnames(x)
  if (!is.null(cn)) {
    colnames(coef_draw) <- cn
    colnames(coef_mode) <- cn
  }
  rn <- group_info$ids
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
    group_info    = group_info,
    group_results = group_results,
    y             = y,
    x             = x,
    offset        = offset2,
    prior.weights = wt,
    prior_lists   = prior_block,
    call          = match.call()
  )
  class(outlist) <- c("rNormal_reg_group", "list")
  outlist
}
