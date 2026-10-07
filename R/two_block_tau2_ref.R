#' ICM / joint-mode tau^2 plug-in from a Block~2 \code{pfamily} (prior spec)
#'
#' For \code{dNormal}, returns fixed \code{dispersion}.  For ING, returns the
#' prior mean \eqn{E[\tau^2_k] = rate/(shape - 1)} from the calibrated
#' \code{pfamily_list(Prior_Setup_GLMM(...))} object.
#' @noRd
.two_block_tau2_plug_in_from_pfamily <- function(pf) {
  pl <- pf$prior_list
  if (identical(pf$pfamily, "dNormal")) {
    d <- as.numeric(pl$dispersion)
    if (!is.finite(d) || d <= 0) {
      stop(
        "dNormal pfamily requires a positive finite 'dispersion' plug-in.",
        call. = FALSE
      )
    }
    return(d)
  }
  if (identical(pf$pfamily, "dIndependent_Normal_Gamma")) {
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (!is.finite(shape) || shape <= 1 || !is.finite(rate) || rate <= 0) {
      stop(
        "ING pfamily requires shape > 1 and positive rate for tau^2 plug-in ",
        "rate/(shape - 1).",
        call. = FALSE
      )
    }
    return(rate / (shape - 1))
  }
  stop("Unsupported pfamily: ", pf$pfamily, call. = FALSE)
}

#' Named tau^2 plug-in vector for ICM / joint mode from \code{pfamily_list}
#' @noRd
.two_block_tau2_plug_in_vector <- function(pfamily_list, re_names) {
  stats::setNames(
    vapply(re_names, function(k) {
      .two_block_tau2_plug_in_from_pfamily(pfamily_list[[k]])
    }, numeric(1)),
    re_names
  )
}

#' Random-effect prior precision matrix derived from a Block~2 pfamily_list
#'
#' \code{P = solve(diag(tau2_k))}, where each \eqn{\tau^2_k} is the same
#' ICM/joint-mode plug-in \code{.two_block_tau2_plug_in_vector()} returns
#' (fixed \code{dispersion} for \code{dNormal}, prior mean
#' \eqn{rate/(shape - 1)} for \code{dIndependent_Normal_Gamma}). Used by the
#' \code{rLMMNormal_reg}/\code{rLMMindepNormalGamma_reg} exports so \code{P}
#' is always consistent with \code{pfamily_list} instead of being supplied
#' (and never cross-checked) separately.
#' @noRd
.rLMM_P_from_pfamily_list <- function(pfamily_list, re_names) {
  p_re <- length(re_names)
  tau2 <- .two_block_tau2_plug_in_vector(pfamily_list, re_names)
  Sigma_ranef <- diag(tau2, nrow = p_re, ncol = p_re)
  dimnames(Sigma_ranef) <- list(re_names, re_names)
  P <- solve(Sigma_ranef)
  dimnames(P) <- list(re_names, re_names)
  P
}

#' Rate-calibration plug-in tau^2_k from a Block~2 pfamily (not ICM starts)
#'
#' For \code{dNormal}, returns fixed \code{dispersion}.  For ING, returns
#' \code{rate / shape} = \eqn{1/E[1/\tau^2]} (precision-mean plug-in) used
#' only for \eqn{\lambda^*} / conservative rate bounds — not for ICM starts.
#' @noRd
.two_block_tau2_ref_from_pfamily <- function(pf) {
  pl <- pf$prior_list
  if (identical(pf$pfamily, "dNormal")) {
    return(as.numeric(pl$dispersion))
  }
  if (identical(pf$pfamily, "dIndependent_Normal_Gamma")) {
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (is.finite(shape) && shape > 0 && is.finite(rate) && rate > 0) {
      return(rate / shape)
    }
    stop(
      "pfamily dIndependent_Normal_Gamma requires positive shape and rate ",
      "so the precision-mean plug-in tau^2 = rate/shape is defined.",
      call. = FALSE
    )
  }
  stop("Unsupported pfamily: ", pf$pfamily, call. = FALSE)
}


#' Plug-in tau^2 vector from pilot dispersion draws (precision-mean plug-in)
#'
#' Returns \code{1 / colMeans(1 / dispersion_fixef_draws)} per RE component —
#' the sample analogue of \eqn{1/E[1/\tau^2]} for main-stage chain starts.
#' @noRd
.two_block_tau2_start_from_dispersion_draws <- function(
    dispersion_fixef_draws,
    re_names = colnames(dispersion_fixef_draws)
) {
  if (is.null(dispersion_fixef_draws) || !is.matrix(dispersion_fixef_draws)) {
    stop(
      "'dispersion_fixef_draws' must be a matrix of positive tau^2 values.",
      call. = FALSE
    )
  }
  if (is.null(re_names) || length(re_names) != ncol(dispersion_fixef_draws)) {
    stop(
      "'re_names' must match ncol(dispersion_fixef_draws).",
      call. = FALSE
    )
  }
  tau2 <- vapply(seq_along(re_names), function(j) {
    x <- dispersion_fixef_draws[, j]
    if (!all(is.finite(x)) || any(x <= 0)) {
      stop(
        "dispersion_fixef_draws[*, ", re_names[j],
        "] must be finite and strictly positive.",
        call. = FALSE
      )
    }
    1 / mean(1 / x)
  }, numeric(1))
  stats::setNames(tau2, re_names)
}
