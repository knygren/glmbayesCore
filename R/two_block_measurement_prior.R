#' Build measurement-prior list from matrix-level two-block inputs
#'
#' Converts Block~1 \code{prior_list} and Block~2 \code{pfamily_list} into the
#' \code{measurement_prior_list} structure required by
#' \code{\link{lmerb_posterior_mean}}.
#'
#' @param prior_list_block1 Block~1 prior with \code{P} or \code{Sigma} and,
#'   for Gaussian models, \code{dispersion}.
#' @param pfamily_list Named list of \code{pfamily} objects (validated).
#' @param re_names Random-effect component names.
#' @param x_hyper Named list of hyper design matrices (column names for
#'   Block~2 parameters).
#' @param family Response \code{family} object.
#' @return List with \code{group.dispersion}, \code{group.Sigma}, and
#'   \code{pop.prior_list}.
#' @noRd
.two_block_measurement_prior_list <- function(
    prior_list_block1,
    pfamily_list,
    re_names,
    x_hyper,
    family
) {
  family <- .two_block_normalize_family(family)
  is_gaussian <- identical(family$family, "gaussian")

  if (!is.null(prior_list_block1$Sigma)) {
    group.Sigma <- as.matrix(prior_list_block1$Sigma)
  } else if (!is.null(prior_list_block1$P)) {
    group.Sigma <- solve(as.matrix(prior_list_block1$P))
  } else {
    stop(
      "prior_list must contain 'P' or 'Sigma' for ICM.",
      call. = FALSE
    )
  }

  group.dispersion <- if (is_gaussian) {
    d <- prior_list_block1$dispersion
    if (is.null(d)) {
      stop(
        "prior_list must contain 'dispersion' for gaussian() ICM.",
        call. = FALSE
      )
    }
    as.numeric(d)
  } else {
    NULL
  }

  prior_list <- stats::setNames(
    vector("list", length(re_names)),
    re_names
  )

  for (k in re_names) {
    pf <- pfamily_list[[k]]
    pl <- pf$prior_list
    par_names <- colnames(x_hyper[[k]])
    if (is.null(par_names) || !length(par_names)) {
      stop(
        "x_hyper[[\"", k, "\"]] must have colnames for ICM.",
        call. = FALSE
      )
    }

    mu_k <- as.numeric(pl$mu)
    if (length(mu_k) != length(par_names)) {
      stop(
        "pfamily_list[[\"", k, "\"]]$prior_list$mu length (", length(mu_k),
        ") must match ncol(x_hyper[[\"", k, "\"]]) (", length(par_names), ").",
        call. = FALSE
      )
    }
    mu_nms <- rownames(pl$mu)
    if (!is.null(mu_nms) && all(nzchar(mu_nms))) {
      if (!setequal(mu_nms, par_names)) {
        stop(
          "Parameter names of pfamily_list[[\"", k,
          "\"]] do not match x_hyper[[\"", k, "\"]] columns.",
          call. = FALSE
        )
      }
      ord <- match(par_names, mu_nms)
      mu_k <- mu_k[ord]
      Sigma_k <- as.matrix(pl$Sigma)[ord, ord, drop = FALSE]
    } else {
      Sigma_k <- as.matrix(pl$Sigma)
    }
    names(mu_k) <- par_names
    dimnames(Sigma_k) <- list(par_names, par_names)

    d_k <- if (identical(pf$pfamily, "dNormal")) {
      pl$dispersion
    } else {
      .two_block_tau2_plug_in_from_pfamily(pf)
    }
    if (is.null(d_k) || !is.numeric(d_k) || length(d_k) != 1L ||
        !is.finite(d_k) || d_k <= 0) {
      stop(
        "pfamily_list[[\"", k, "\"]] must supply a positive plug-in tau^2 ",
        "(dNormal dispersion or ING rate/shape = 1/E[1/tau^2] from the prior spec).",
        call. = FALSE
      )
    }

    prior_list[[k]] <- list(
      mu         = mu_k,
      Sigma      = Sigma_k,
      dispersion = as.numeric(d_k)
    )
  }

  list(
    group.dispersion = group.dispersion,
    group.Sigma      = group.Sigma,
    pop.prior_list   = prior_list
  )
}

#' Iterated conditional means/modes at matrix-level two-block inputs
#'
#' @param design Design list for ICM (\code{y}, \code{D}, \code{group},
#'   \code{W}, \code{groupef.names}).
#' @param prior_list Block~1 prior list.
#' @param pfamily_list Validated Block~2 \code{pfamily_list}.
#' @param re_names Random-effect component names.
#' @param family Response family (default \code{gaussian()}).
#' @param tol,maxit ICM convergence controls.
#' @return List with \code{start}, \code{b_start}, and \code{icm}.
#' @noRd
.two_block_icm_at_start <- function(
    design,
    prior_list,
    pfamily_list,
    re_names,
    family = gaussian(),
    tol   = 1e-10,
    maxit = 200L
) {
  family <- .two_block_normalize_family(family)
  mpl <- .two_block_measurement_prior_list(
    prior_list_block1 = prior_list,
    pfamily_list      = pfamily_list,
    re_names          = re_names,
    x_hyper           = design$W,
    family            = family
  )

  if (!identical(family$family, "gaussian")) {
    stop("'family' must be gaussian.", call. = FALSE)
  }
  pm <- lmerb_posterior_mean(
    design,
    mpl,
    tol   = tol,
    maxit = maxit
  )

  list(
    start   = pm$fixef,
    b_start = pm$b_mean,
    icm     = list(
      converged  = pm$converged,
      iterations = pm$iterations,
      delta      = pm$delta
    )
  )
}

#' Validate gap tolerance for pilot chain count derivation
#' @noRd
.two_block_validate_gap_tol <- function(gap_tol) {
  if (is.null(gap_tol)) {
    return(invisible(NULL))
  }
  if (!is.numeric(gap_tol) || length(gap_tol) != 1L ||
      !is.finite(gap_tol) || gap_tol <= 0 || gap_tol >= 1) {
    stop("'gap_tol' must be NULL or a single value in (0, 1).", call. = FALSE)
  }
  invisible(as.numeric(gap_tol))
}

