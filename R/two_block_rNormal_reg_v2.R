#' Two-block Gibbs sampler with pfamily Block 2 priors (development v2)
#'
#' Development version of \code{\link{two_block_rNormal_reg}} where the
#' Block~2 priors are supplied as \code{pfamily} objects instead of bare
#' prior lists.  Each component of \code{pfamily_list} may be a
#' \code{\link{dNormal}} prior (conjugate gamma_k draw at fixed dispersion,
#' identical to v1) or a \code{\link{dIndependent_Normal_Gamma}} prior, in
#' which case Block~2 makes a joint (gamma_k, tau^2_k) draw via the
#' likelihood-subgradient envelope sampler (the same path as \code{rglmb}
#' with an ING pfamily) and the sampled tau^2_k is fed back into the
#' Block~1 prior precision on the next inner step.
#'
#' With \code{dNormal} priors throughout, this function produces draws that
#' are identical to \code{\link{two_block_rNormal_reg}} under the same seed;
#' that equivalence is the regression gate for the v2 development track.
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
#' @param pfamily_list Named list of \code{pfamily} objects, one per column
#'   of \code{x}: \code{\link{dNormal}} or
#'   \code{\link{dIndependent_Normal_Gamma}}.
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
#' @return Object of class \code{c("two_block_rNormal_reg_v2",
#'   "two_block_rNormal_reg")}.  Same fields as
#'   \code{\link{two_block_rNormal_reg}}, plus
#'   \code{dispersion_fixef_draws}: an \code{n x p_re} matrix of the Block~2
#'   dispersion (tau^2_k) at each stored draw (constant columns for
#'   \code{dNormal} components).
#' @family simfuncs
#' @seealso \code{\link{two_block_rNormal_reg}}, \code{\link{dNormal}},
#'   \code{\link{dIndependent_Normal_Gamma}}
#' @export
two_block_rNormal_reg_v2 <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    prior_list_block1,
    pfamily_list,
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

  pfamily_list <- .two_block_validate_pfamily_list(pfamily_list, re_names)

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

  famfunc_block1 <- glmbfamfunc(if (is_gaussian) gaussian() else family)
  famfunc_gauss <- glmbfamfunc(gaussian())
  n_envopt_use <- if (is.null(n_envopt)) 1L else as.integer(n_envopt)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  x_hyper_mats <- lapply(x_hyper, as.matrix)

  cpp_out <- .two_block_rNormal_reg_v2_cpp(
    n                 = n,
    m_convergence     = m_convergence,
    y                 = y,
    x                 = x,
    block             = block,
    x_hyper           = x_hyper_mats,
    prior_list_block1 = prior_list_block1,
    dispersion_block1 = block1_prior_meta$dispersion,
    ddef_block1       = block1_prior_meta$ddef,
    pfamily_list      = pfamily_list,
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

  dispersion_fixef_draws <- cpp_out$dispersion_fixef_draws
  dimnames(dispersion_fixef_draws) <- list(NULL, re_names)

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
      fixef_draws            = fixef_draws,
      coefficients           = coefficients,
      fixef_last             = fixef,
      b_last                 = b_i,
      mu_all_last            = mu_all,
      dispersion_fixef_draws = dispersion_fixef_draws,
      pfamily_list           = pfamily_list,
      family                 = family,
      n                      = n,
      m_convergence          = m_convergence,
      sampling               = sampling,
      fixef_start            = fixef_start,
      re_coef_names          = re_names,
      group_levels           = group_levels,
      group_name             = group_name,
      call                   = cl
    ),
    class = c("two_block_rNormal_reg_v2", "two_block_rNormal_reg")
  )
}

#' Validate a Block 2 pfamily list (v2 contract)
#'
#' Checks that each component is a supported \code{pfamily} object and that
#' the type-specific hyperparameters needed by the v2 driver are present.
#'
#' @param pfamily_list Named list of \code{pfamily} objects.
#' @param re_names Character vector of RE component names (defines order).
#' @return The validated list, reordered to match \code{re_names}.
#' @noRd
.two_block_validate_pfamily_list <- function(pfamily_list, re_names) {
  if (!is.list(pfamily_list)) {
    stop("'pfamily_list' must be a named list of pfamily objects.",
         call. = FALSE)
  }
  if (!setequal(names(pfamily_list), re_names)) {
    stop("names(pfamily_list) must match re_coef_names.", call. = FALSE)
  }
  pfamily_list <- pfamily_list[re_names]

  supported <- c("dNormal", "dIndependent_Normal_Gamma")
  for (k in re_names) {
    pf <- pfamily_list[[k]]
    if (!inherits(pf, "pfamily") ||
        is.null(pf$pfamily) || is.null(pf$prior_list)) {
      stop(
        "pfamily_list[[\"", k, "\"]] must be a pfamily object ",
        "(e.g. dNormal() or dIndependent_Normal_Gamma()).",
        call. = FALSE
      )
    }
    if (!pf$pfamily %in% supported) {
      stop(
        "pfamily_list[[\"", k, "\"]]: unsupported pfamily \"", pf$pfamily,
        "\" (allowed: ", paste(supported, collapse = ", "), ").",
        call. = FALSE
      )
    }
    pl <- pf$prior_list
    if (identical(pf$pfamily, "dNormal")) {
      if (is.null(pl$dispersion) || isTRUE(pl$ddef)) {
        stop(
          "pfamily_list[[\"", k, "\"]]: dNormal() requires an explicit ",
          "dispersion for the Block 2 (gaussian) update.",
          call. = FALSE
        )
      }
    } else {
      if (is.null(pl$shape) || is.null(pl$rate) ||
          !is.numeric(pl$shape) || !is.numeric(pl$rate) ||
          pl$shape[1L] <= 0 || pl$rate[1L] <= 0) {
        stop(
          "pfamily_list[[\"", k, "\"]]: dIndependent_Normal_Gamma() requires ",
          "positive 'shape' and 'rate'.",
          call. = FALSE
        )
      }
      if (is.null(pl$disp_lower) || !is.numeric(pl$disp_lower) ||
          pl$disp_lower[1L] <= 0) {
        stop(
          "pfamily_list[[\"", k, "\"]]: dIndependent_Normal_Gamma() requires ",
          "a positive 'disp_lower' (initial/calibration tau^2_k plug-in).",
          call. = FALSE
        )
      }
    }
  }
  pfamily_list
}

#' Convergence rate for the v2 (pfamily) two-block sampler
#'
#' Thin wrapper around \code{\link{two_block_rate}} that accepts the Block~2
#' priors as \code{pfamily} objects.  For \code{dNormal} components the fixed
#' \code{dispersion} is used; for \code{dIndependent_Normal_Gamma} components
#' the conservative \code{disp_lower} plug-in is used (the rate is then an
#' upper bound over the truncated tau^2 range).
#'
#' @param x,block,x_hyper,prior_list_block1,weights,family,group_levels As in
#'   \code{\link{two_block_rate}}.
#' @param pfamily_list Named list of \code{pfamily} objects (one per RE
#'   component), as in \code{\link{two_block_rNormal_reg_v2}}.
#' @return Object of class \code{"two_block_rate"}.
#' @family simfuncs
#' @seealso \code{\link{two_block_rate}},
#'   \code{\link{two_block_rNormal_reg_v2}}
#' @export
two_block_rate_v2 <- function(x,
                              block,
                              x_hyper,
                              prior_list_block1,
                              pfamily_list,
                              weights = NULL,
                              family = gaussian(),
                              group_levels = levels(block)) {
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
    x = x, block = block, x_hyper = x_hyper,
    prior_list_block1 = prior_list_block1,
    prior_list_block2 = prior_list_block2,
    weights = weights, family = family, group_levels = group_levels
  )
}
