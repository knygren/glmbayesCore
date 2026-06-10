#' Two-block Gibbs sampler for hierarchical regression
#'
#' Runs the coupled Block~1 / Block~2 Gibbs sampler for two-block mixed models.
#' Block~1 draws group-level random effects \eqn{b_j}; Block~2 updates hyper
#' means \eqn{\gamma_k} via \code{\link{multi_rNormal_reg}} (always Gaussian).
#'
#' Block~1 follows the same path as \code{\link{rNormal_reg}}:
#' \code{\link{block_rNormalReg}} when \code{family = gaussian()}, otherwise
#' \code{\link{block_rNormalGLM}} for the GLM envelope path.
#'
#' @param n Number of stored draws.
#' @param y Response vector of length \code{nrow(x)}.
#' @param x Level-1 design matrix \code{Z} (\code{l2 x p_re}).
#' @param block Grouping factor or block partition of length \code{l2}.
#' @param x_hyper Named list of group-level design matrices \code{X_k}
#'   (\code{J x q_k}), one per column of \code{x}.
#' @param prior_list_block1 Prior for Block~1: \code{P} or \code{Sigma},
#'   \code{dispersion} (required for \code{gaussian()}), optional \code{ddef}.
#'   \code{mu} is updated internally.
#' @param prior_list_block2 Named list of Block~2 prior lists passed to
#'   \code{\link{multi_rNormal_reg}}.
#' @param fixef_start Named list of hyper-parameter vectors at which each inner
#'   chain is initialised.
#' @param re_coef_names Character vector naming columns of \code{x}.
#' @param group_levels Character vector defining row order of Block~1 draws.
#' @param group_name Name for the grouping column in \code{coefficients}.
#' @param m_convergence Number of inner Gibbs steps per stored draw.
#' @param sampling Sampling scheme; only \code{"replicate"} is implemented.
#' @param family Response \code{\link[stats]{family}} for Block~1 (default
#'   \code{gaussian()}). Block~2 always uses \code{gaussian()}.
#' @param offset,weights Passed to Block~1 (length \code{l2} or recycled).
#' @param Gridtype,n_envopt,use_parallel,use_opencl,verbose Passed to Block~1
#'   when \code{family} is not Gaussian.
#' @param seed Optional RNG seed.
#' @param progbar Logical; show a text progress bar.
#' @return Object of class \code{"two_block_rNormal_reg"}.
#' @family simfuncs
#' @seealso \code{\link{build_mu_all}}, \code{\link{two_block_rNormal_reg_v2}},
#'   \code{\link{block_rNormalReg}}, \code{\link{block_rNormalGLM}},
#'   \code{\link{multi_rNormal_reg}}, \code{\link{rNormal_reg}}
#' @export
two_block_rNormal_reg <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    prior_list_block1,
    prior_list_block2,
    fixef_start,
    re_coef_names = colnames(x),
    group_levels = levels(block),
    group_name = NULL,
    m_convergence = 10L,
    sampling = c("replicate", "chain"),
    family = gaussian(),
    offset = NULL,
    weights = 1,
    Gridtype = 2L,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE,
    seed = NULL,
    progbar = TRUE) {

  cl <- match.call()
  sampling <- match.arg(sampling)
  if (!identical(sampling, "replicate")) {
    stop("Only sampling = \"replicate\" is implemented.", call. = FALSE)
  }

  family <- .two_block_normalize_family(family)
  is_gaussian <- identical(family$family, "gaussian")

  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }
  m_convergence <- as.integer(m_convergence[1L])
  if (m_convergence < 1L) {
    stop("'m_convergence' must be at least 1.", call. = FALSE)
  }

  y <- as.vector(y)
  x <- as.matrix(x)
  l2 <- nrow(x)
  if (length(y) != l2) {
    stop("length(y) must equal nrow(x).", call. = FALSE)
  }

  if (is.null(re_coef_names) || length(re_coef_names) != ncol(x)) {
    re_coef_names <- if (ncol(x) >= 1L) {
      cn <- colnames(x)
      if (is.null(cn) || length(cn) != ncol(x)) paste0("RE", seq_len(ncol(x))) else cn
    } else {
      stop("'x' must have at least one column.", call. = FALSE)
    }
  }
  colnames(x) <- re_coef_names
  re_names <- re_coef_names

  group_levels <- as.character(group_levels)
  if (length(group_levels) < 1L) {
    stop("'group_levels' must contain at least one level.", call. = FALSE)
  }

  if (is.null(group_name) || !nzchar(group_name)) {
    group_name <- tryCatch(
      deparse(substitute(block))[1L],
      error = function(e) "group"
    )
    if (!nzchar(group_name)) group_name <- "group"
  }

  if (!is.list(x_hyper) || is.data.frame(x_hyper)) {
    stop("'x_hyper' must be a list of design matrices.", call. = FALSE)
  }
  if (length(x_hyper) != length(re_names)) {
    stop(
      "length(x_hyper) must equal ncol(x) = ", length(re_names), ".",
      call. = FALSE
    )
  }
  if (!setequal(names(x_hyper), re_names)) {
    x_hyper <- x_hyper[re_names]
  }

  if (!is.list(prior_list_block2)) {
    stop("'prior_list_block2' must be a named list.", call. = FALSE)
  }
  if (!setequal(names(prior_list_block2), re_names)) {
    stop(
      "names(prior_list_block2) must match re_coef_names.",
      call. = FALSE
    )
  }

  if (!is.list(fixef_start) || is.null(names(fixef_start))) {
    stop("'fixef_start' must be a named list.", call. = FALSE)
  }
  if (!setequal(names(fixef_start), re_names)) {
    stop("names(fixef_start) must match re_coef_names.", call. = FALSE)
  }
  fixef_start <- fixef_start[re_names]

  block1_prior_meta <- .two_block_validate_block1_prior(
    prior_list_block1,
    family = family
  )

  ## Observation-level vectors normalized once (same rules as the block
  ## wrappers, which used to do this per inner Gibbs step).
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

  ## Family closures resolved once: Block 1 uses the response family, Block 2
  ## is always Gaussian (multi_rNormal_reg/rNormal_reg used to resolve these
  ## per call).
  famfunc_block1 <- glmbfamfunc(if (is_gaussian) gaussian() else family)
  famfunc_gauss <- glmbfamfunc(gaussian())
  n_envopt_use <- if (is.null(n_envopt)) 1L else as.integer(n_envopt)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  x_hyper_mats <- lapply(x_hyper, as.matrix)

  ## The n x m_convergence Gibbs loop runs in C++
  ## (.two_block_rNormal_reg_cpp -> two_block_rNormal_reg_cpp_export):
  ## mu_all -> Block 1 (block_rNormalReg/block_rNormalGLM C++ exports) ->
  ## Block 2 (rNormalReg per RE component), exactly as the former R loop.
  cpp_out <- .two_block_rNormal_reg_cpp(
    n                 = n,
    m_convergence     = m_convergence,
    y                 = y,
    x                 = x,
    block             = block,
    x_hyper           = x_hyper_mats,
    prior_list_block1 = prior_list_block1,
    dispersion_block1 = block1_prior_meta$dispersion,
    ddef_block1       = block1_prior_meta$ddef,
    prior_list_block2 = prior_list_block2,
    fixef_start       = fixef_start,
    group_levels      = group_levels,
    family            = family$family,
    link              = family$link,
    f2                = famfunc_block1$f2,
    f3                = famfunc_block1$f3,
    f2_gauss          = famfunc_gauss$f2,
    f3_gauss          = famfunc_gauss$f3,
    offset            = offset2,
    wt                = wt,
    Gridtype          = as.integer(Gridtype),
    n_envopt          = n_envopt_use,
    use_parallel      = use_parallel,
    use_opencl        = use_opencl,
    verbose           = verbose,
    progbar           = isTRUE(progbar)
  )

  J <- length(group_levels)
  p_re <- length(re_names)
  group_ids <- as.character(cpp_out$group_ids)

  fixef_draws <- stats::setNames(cpp_out$fixef_draws, re_names)
  for (k in re_names) {
    dimnames(fixef_draws[[k]]) <- list(NULL, names(fixef_start[[k]]))
  }

  fixef <- stats::setNames(cpp_out$fixef_last, re_names)

  b_arr <- array(as.numeric(cpp_out$b_draws), dim = c(J, p_re, n))
  b_i <- matrix(b_arr[, , n], nrow = J, ncol = p_re,
                dimnames = list(group_ids, re_names))

  mu_all <- cpp_out$mu_all_last
  dimnames(mu_all) <- list(re_names, group_levels)

  coef_cols <- c("draw", group_name, re_names)
  draw_rows <- vector("list", n)
  for (i in seq_len(n)) {
    draw_df <- data.frame(
      draw = rep(i, J),
      stringsAsFactors = FALSE
    )
    draw_df[[group_name]] <- group_ids
    for (jj in seq_len(p_re)) {
      draw_df[[re_names[jj]]] <- b_arr[, jj, i]
    }
    draw_rows[[i]] <- draw_df
  }

  coefficients <- do.call(rbind, draw_rows)
  rownames(coefficients) <- NULL
  coefficients <- coefficients[, coef_cols, drop = FALSE]

  structure(
    list(
      fixef_draws   = fixef_draws,
      coefficients  = coefficients,
      fixef_last    = fixef,
      b_last        = b_i,
      mu_all_last   = mu_all,
      family        = family,
      n             = n,
      m_convergence = m_convergence,
      sampling      = sampling,
      fixef_start   = fixef_start,
      re_coef_names = re_names,
      group_levels  = group_levels,
      group_name    = group_name,
      call          = cl
    ),
    class = "two_block_rNormal_reg"
  )
}

#' @noRd
.two_block_normalize_family <- function(family) {
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    stop("'family' not recognized.", call. = FALSE)
  }

  okfamilies <- c(
    "gaussian", "poisson", "binomial",
    "quasipoisson", "quasibinomial", "Gamma"
  )
  if (!family$family %in% okfamilies) {
    stop(
      "family \"", family$family, "\" is not supported by two_block_rNormal_reg.",
      call. = FALSE
    )
  }

  oklinks <- switch(
    family$family,
    gaussian = "identity",
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

  family
}

#' @noRd
.two_block_validate_block1_prior <- function(prior_list_block1, family) {
  if (!is.list(prior_list_block1)) {
    stop("'prior_list_block1' must be a list.", call. = FALSE)
  }
  if (is.null(prior_list_block1$P) && is.null(prior_list_block1$Sigma)) {
    stop("prior_list_block1 must contain 'P' or 'Sigma'.", call. = FALSE)
  }

  ddef <- if ("ddef" %in% names(prior_list_block1)) {
    prior_list_block1$ddef
  } else {
    is.null(prior_list_block1$dispersion)
  }

  dispersion <- prior_list_block1$dispersion

  if (identical(family$family, "gaussian")) {
    if (is.null(dispersion)) {
      stop(
        "prior_list_block1 must contain 'dispersion' for gaussian() Block~1.",
        call. = FALSE
      )
    }
    if (isTRUE(ddef)) {
      stop(
        "For gaussian() Block~1, dNormal() requires an explicit dispersion.",
        call. = FALSE
      )
    }
  }

  if (family$family %in% c("gaussian", "Gamma") && isTRUE(ddef)) {
    stop(
      "For gaussian() and Gamma() models, dNormal() requires an explicit dispersion.",
      call. = FALSE
    )
  }

  list(ddef = ddef, dispersion = dispersion)
}

#' @noRd
.two_block_block1_prior_list <- function(prior_list_block1, mu_all, meta) {
  out <- list(
    mu         = mu_all,
    dispersion = meta$dispersion,
    ddef       = meta$ddef
  )
  if (!is.null(prior_list_block1$P)) {
    out$P <- prior_list_block1$P
  }
  if (!is.null(prior_list_block1$Sigma)) {
    out$Sigma <- prior_list_block1$Sigma
  }
  if (is.null(out$P) && is.null(out$Sigma)) {
    stop("prior_list_block1 must contain 'P' or 'Sigma'.", call. = FALSE)
  }
  out
}

#' @noRd
.two_block_mu_all <- function(fixef, x_hyper, re_names, group_levels) {
  p_re <- length(re_names)
  J    <- length(group_levels)
  mu_all <- matrix(NA_real_, nrow = p_re, ncol = J,
                   dimnames = list(re_names, group_levels))
  for (i in seq_len(p_re)) {
    k       <- re_names[i]
    gamma_k <- fixef[[k]]
    X_k     <- as.matrix(x_hyper[[k]])
    rn      <- rownames(X_k)
    if (is.null(rn)) {
      if (nrow(X_k) != J) {
        stop(
          "nrow(x_hyper[[", k, "]]) must equal length(group_levels).",
          call. = FALSE
        )
      }
      for (j in seq_len(J)) {
        mu_all[i, j] <- sum(X_k[j, , drop = TRUE] * gamma_k)
      }
    } else {
      for (j in seq_len(J)) {
        mu_all[i, j] <- sum(X_k[group_levels[j], , drop = TRUE] * gamma_k)
      }
    }
  }
  mu_all
}
