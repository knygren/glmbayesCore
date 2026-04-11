#' Compute Gaussian Normal-Gamma Calibration Components
#'
#' Gaussian Normal-Gamma calibration function used by \code{\link{Prior_Setup}}.
#' Given Gaussian-model sufficient inputs and a dispersion-independent
#' coefficient prior covariance \eqn{\Sigma_0} (Chapter 11 framing), this
#' function computes calibrated Gaussian prior quantities:
#' \code{dispersion}, \code{shape}, \code{rate}, and \code{Sigma}, and returns
#' the input \code{Sigma_0} as \code{Sigma_0}.
#'
#' The function is structured as a step-wise pipeline:
#' \enumerate{
#'   \item Validate dimensions and numeric inputs.
#'   \item Compute weighted residual sum of squares at \code{bhat}.
#'   \item Build \eqn{X'WX}, invert it, and construct \eqn{S_{marg}}.
#'   \item Map \code{n_prior} to Gamma shape and calibrate rate.
#'   \item Calibrate \code{dispersion}, \code{rate}, and coefficient covariance.
#'   \item Return calibrated outputs.
#' }
#'
#' @details
#' Chapter 11 emphasizes a parameterization where \eqn{\Sigma_0} is constant
#' with respect to dispersion. In that framing, \eqn{\Sigma_0} is the prior
#' covariance on the precision-weighted coefficient scale, while the returned
#' \code{Sigma} is the covariance on the coefficient scale after calibration.
#'
#' A common choice is to set \eqn{\Sigma_0} proportional to the inverse weighted
#' Gram matrix, \eqn{(X^T W X)^{-1}}. In the scalar-\code{pwt} (single prior
#' weight shared by all coefficients) Zellner-style setup used in
#' \code{\link{Prior_Setup}}, the implied form is:
#' \deqn{
#' \Sigma_0 = \frac{1-\mathrm{pwt}}{\mathrm{pwt}} (X^T W X)^{-1}.
#' }
#' More generally, users may pass any positive-definite \eqn{\Sigma_0} encoding
#' alternative prior structure on coefficients.
#'
#' The function computes:
#' \itemize{
#'   \item \eqn{S_{marg} = RSS_w + (\hat\beta-\mu)^T(\Sigma_0 + (X^T W X)^{-1})^{-1}(\hat\beta-\mu)}
#'   \item Prior Gamma \code{shape} is \eqn{(n_{\mathrm{prior}}+1)/2}{}. Posterior Gamma shape is
#'   \eqn{a_n = (n_{\mathrm{prior}}+n_w+1)/2}{} with \eqn{n_w}{} = \code{n_effective}.
#'   \item Calibrated \code{dispersion} equals \eqn{S_{marg}/(n_w-p)}{} with \eqn{p=\texttt{ncol}(X)}{},
#'   i.e.\ the usual weighted residual-df denominator.
#'   \item Prior rate \eqn{b_0 = \frac{1}{2} S_{marg}(n_{\mathrm{prior}}+p-1)/(n_w-p)}{} so that
#'   \eqn{b_n=b_0+S_{marg}/2}{} gives \eqn{E[\sigma^2\mid y]=b_n/(a_n-1)=S_{marg}/(n_w-p)}{}.
#'   The \eqn{+1}{} in the prior shape and \eqn{p-1}{} in the rate numerator pair so
#'   \eqn{a_n-1}{} and \eqn{b_n}{} share a factor \eqn{n_{\mathrm{prior}}+n_w-1}{}; this
#'   preserves the posterior mean dispersion while improving prior propriety when
#'   \eqn{p=1}{} and \eqn{n_{\mathrm{prior}}}{} is small.
#'   \item \eqn{\Sigma = (n_w/n_{\mathrm{prior}}) E[\sigma^2 \mid y] (X^T W X)^{-1}}{}
#' }
#'
#' Limiting behavior:
#' \itemize{
#'   \item As \eqn{n_{prior} \to \infty}, prior information dominates and the
#'   Gamma prior on precision becomes increasingly concentrated.
#'   \item The expected coefficient location remains centered at \eqn{\mu}; this
#'   function does not change the mean vector, only scale components.
#'   \item As \eqn{n_{prior}\to 0^+}, the returned expected dispersion tends to
#'   \eqn{S_{marg}/(n_{\mathrm{effective}}-p)} when \eqn{n_{\mathrm{effective}}>p}.
#'   With \eqn{S_{marg}\to \mathrm{RSS}_w}, this matches the usual weighted
#'   Gaussian unbiased dispersion \eqn{\mathrm{RSS}_w/(n_{\mathrm{effective}}-p)}.
#'   \item Strict positivity of \eqn{b_0}{} requires \eqn{n_{\mathrm{prior}}+p>1}{}.
#'   Prior \code{shape} is at least \eqn{1/2}{} because \eqn{\texttt{shape}=(n_{\mathrm{prior}}+1)/2}{}.
#'   As \eqn{n_{prior}\to 0^+}, terms can become very small; this
#'   corresponds to a very diffuse precision prior and may be numerically
#'   delicate near the boundary.
#'   \item In the Chapter 11 scalar-\code{pwt} path
#'   (\eqn{pwt=n_{prior}/(n_{prior}+n_w)} with fixed \eqn{n_w>p}), one has
#'   \eqn{pwt\to 0}, \eqn{S_{marg}\to \mathrm{RSS}_w},
#'   \eqn{E[\sigma^2\mid y]\to \mathrm{RSS}_w/(n_w-p)}, and \eqn{V_n\to (X^TWX)^{-1}}
#'   under the conjugate Zellner setup, so marginal coefficient covariance
#'   aligns with \eqn{\mathrm{RSS}_w/(n_w-p)\,(X^TWX)^{-1}} in the same limit.
#' }
#'
#' @param X Numeric model matrix with \code{nrow(X) == length(Y)}.
#' @param Y Numeric response vector.
#' @param weights Numeric case weights vector of length \code{nrow(X)}.
#' @param offset Numeric offset vector of length \code{nrow(X)}.
#' @param dispersion Optional scalar dispersion input (default \code{NULL}).
#' Must be \code{NULL} or a single positive finite numeric value.
#' @param n_effective Positive scalar effective sample size, typically
#'   \code{sum(weights)} for Gaussian models in this package.
#' @param bhat Numeric coefficient vector (typically full-model MLE), length
#'   \code{ncol(X)}.
#' @param mu Numeric prior mean vector (or one-column matrix coercible to vector)
#'   of length \code{ncol(X)}.
#' @param Sigma_0 Dispersion-independent prior covariance matrix on coefficients,
#'   dimension \code{[p x p]} where \code{p = ncol(X)}.
#' @param Sigma Optional coefficient-scale covariance matrix from upstream
#'   \code{Prior_Setup()} plumbing (default \code{NULL}). When provided, the
#'   returned \code{Sigma} is set to this matrix and the returned
#'   \code{Sigma_0} is set to \code{Sigma / dispersion} using the returned
#'   \code{dispersion}. If \code{Sigma} is \code{NULL} but \code{dispersion} is
#'   provided, the returned \code{dispersion} is set to that input and the
#'   returned \code{Sigma} is set to \code{Sigma_0 * dispersion}. When both are
#'   \code{NULL}, the existing calibrated \code{Sigma} / input \code{Sigma_0}
#'   path is used.
#' @param n_prior Positive scalar effective prior sample size.
#' @param shape_df Defunct. If non-\code{NULL}, a warning is issued; the argument is ignored.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{dispersion}: calibrated Gaussian dispersion scalar.
#'   \item \code{shape}: Gamma shape for residual precision.
#'   \item \code{rate}: Gamma rate for residual precision.
#'   \item \code{Sigma}: calibrated coefficient prior covariance matrix.
#'   \item \code{Sigma_0}: the dispersion-independent prior covariance matrix passed in
#'     via argument \code{Sigma_0} (same matrix, with \code{dimnames} taken from
#'     \code{colnames(X)} when available).
#' }
#'
#' @export
compute_gaussian_prior <- function(
    X,
    Y,
    weights,
    offset,
    dispersion = NULL,
    n_effective,
    bhat,
    mu,
    Sigma_0,
    Sigma = NULL,
    n_prior,
    shape_df = NULL
) {
  if (!is.null(shape_df)) {
    warning(
      "compute_gaussian_prior(): argument shape_df is defunct and ignored; ",
      "shape and rate use n_prior only.",
      call. = FALSE
    )
  }
  ## ---------------------------------------------------------------------------
  ## Gaussian calibration pipeline:
  ## Step A: Validate inputs and dimensions.
  ## Step B: Compute weighted RSS from (Y, X, bhat, offset, weights).
  ## Step C: Build Gram terms and S_marg using Sigma_0 + (X'WX)^{-1}.
  ## Step D: Gamma shape = (n_prior + 1) / 2; rate from n_prior (not shape_df).
  ## Step E: Calibrate dispersion/rate and map to coefficient Sigma.
  ## Step F: Return calibrated terms.
  ## ---------------------------------------------------------------------------
  if (!is.null(dispersion)) {
    if (!is.numeric(dispersion) || length(dispersion) != 1L ||
        !is.finite(dispersion) || dispersion <= 0) {
      stop("compute_gaussian_prior: dispersion must be NULL or a single positive finite numeric value.", call. = FALSE)
    }
  }
  dispersion_input <- dispersion
  Sigma_input <- Sigma

  ## Step A: validate all required Gaussian inputs.
  n_obs <- NROW(Y)
  if (!is.matrix(X) || NROW(X) != n_obs) {
    stop("compute_gaussian_prior: X must be a matrix with nrow(X) == length(Y).", call. = FALSE)
  }
  if (!is.numeric(Y) || length(Y) != n_obs) {
    stop("compute_gaussian_prior: Y must be a numeric vector with length equal to nrow(X).", call. = FALSE)
  }
  if (!is.numeric(weights) || length(weights) != n_obs) {
    stop("compute_gaussian_prior: weights must be a numeric vector with length equal to nrow(X).", call. = FALSE)
  }
  if (!is.numeric(offset) || length(offset) != n_obs) {
    stop("compute_gaussian_prior: offset must be a numeric vector with length equal to nrow(X).", call. = FALSE)
  }
  p <- NCOL(X)
  if (!is.numeric(bhat) || length(bhat) != p || any(!is.finite(bhat))) {
    stop("compute_gaussian_prior: bhat must be a finite numeric vector with length ncol(X).", call. = FALSE)
  }
  mu_num <- as.numeric(mu)
  if (length(mu_num) != p || any(!is.finite(mu_num))) {
    stop("compute_gaussian_prior: mu must be a finite numeric vector with length ncol(X).", call. = FALSE)
  }
  if (!is.matrix(Sigma_0) || nrow(Sigma_0) != p || ncol(Sigma_0) != p || anyNA(Sigma_0)) {
    stop("compute_gaussian_prior: Sigma_0 must be a numeric [p x p] matrix with no missing values.", call. = FALSE)
  }
  if (!is.numeric(n_prior) || length(n_prior) != 1L || !is.finite(n_prior) || n_prior <= 0) {
    stop("compute_gaussian_prior: n_prior must be a single positive finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(n_effective) || length(n_effective) != 1L || !is.finite(n_effective) || n_effective <= 0) {
    stop("compute_gaussian_prior: n_effective must be a single positive finite numeric value.", call. = FALSE)
  }

  ## Step B: weighted residual sum of squares at bhat.
  res <- as.numeric(Y) - as.numeric(X %*% bhat) - as.numeric(offset)
  rss_weighted <- sum(as.numeric(weights) * res^2)
  if (!is.finite(rss_weighted) || rss_weighted <= 0) {
    stop("compute_gaussian_prior: weighted RSS must be strictly positive.", call. = FALSE)
  }
  if (!is.finite(p) || p < 1L) {
    stop("compute_gaussian_prior: require ncol(X) >= 1.", call. = FALSE)
  }
  if (n_effective <= p) {
    stop(
      "compute_gaussian_prior: require n_effective > p (number of coefficients) for Gaussian dispersion (denominator n_effective - p). ",
      "Got n_effective = ", n_effective, " and p = ", p, ".",
      call. = FALSE
    )
  }

  ## Step C: weighted Gram inverse and S_marg quadratic augmentation.
  XtW <- sweep(X, 1, as.numeric(weights), `*`)
  Gm <- crossprod(XtW, X)
  Ginv <- tryCatch(
    solve(Gm),
    error = function(e) {
      stop("compute_gaussian_prior: cannot invert weighted Gram matrix X'WX. ", conditionMessage(e), call. = FALSE)
    }
  )
  dlt <- matrix(bhat, ncol = 1L) - matrix(mu_num, ncol = 1L)
  M <- Sigma_0 + Ginv
  Mi <- tryCatch(
    solve(M),
    error = function(e) {
      stop("compute_gaussian_prior: cannot invert Sigma_0 + (X'WX)^{-1}. ", conditionMessage(e), call. = FALSE)
    }
  )
  quad <- as.numeric(crossprod(dlt, Mi %*% dlt))
  if (!is.finite(quad) || quad < 0) {
    stop("compute_gaussian_prior: S_marg quadratic form is not finite or nonnegative.", call. = FALSE)
  }
  S_marg <- rss_weighted + quad

  ## Step D: Gamma shape and rate from n_prior only.
  shape <- (n_prior + 1L) / 2
  if (!is.finite(shape) || shape <= 0) {
    stop("compute_gaussian_prior: computed shape must be strictly positive.", call. = FALSE)
  }

  ## Step E: calibrate Gaussian dispersion/rate and implied Sigma.
  ## E[sigma^2|y] = b_n/(a_n-1) with a_n = (n_prior + n_effective + 1)/2,
  ## b_n = b_0 + S_marg/2, dispersion = S_marg/(n_effective - p),
  ## b_0 = (S_marg/2) * (n_prior + p - 1) / (n_effective - p).
  den_resid_df <- n_effective - p
  dispersion_cal <- S_marg / den_resid_df
  b_0_S_marg_formula <- 0.5 * S_marg * (n_prior + p - 1L) / den_resid_df

  if (!is.finite(dispersion_cal) || dispersion_cal <= 0) {
    stop("compute_gaussian_prior: calibrated dispersion (S_marg/(n_effective-p)) is missing or not positive.", call. = FALSE)
  }
  if (!is.finite(b_0_S_marg_formula) || b_0_S_marg_formula <= 0) {
    stop(
      "compute_gaussian_prior: prior rate b_0 is missing or not positive. ",
      "Require n_prior + p > 1. Got n_prior = ", n_prior, ", p = ", p, ".",
      call. = FALSE
    )
  }
  rate <- b_0_S_marg_formula
  Sigma_calibrated <- (n_effective / n_prior) * dispersion_cal * Ginv
  dimnames(Sigma_calibrated) <- list(colnames(X), colnames(X))

  dispersion <- dispersion_cal
  Sigma <- Sigma_calibrated
  Sigma_0_out <- Sigma_0
  dimnames(Sigma_0_out) <- list(colnames(X), colnames(X))

  if (!is.null(Sigma_input)) {
    if (!is.matrix(Sigma_input) || nrow(Sigma_input) != p || ncol(Sigma_input) != p || anyNA(Sigma_input)) {
      stop("compute_gaussian_prior: Sigma must be NULL or a numeric [p x p] matrix with no missing values.", call. = FALSE)
    }
    Sigma <- Sigma_input
    if (!is.null(dispersion_input)) {
      dispersion <- dispersion_input
    }
    dimnames(Sigma) <- list(colnames(X), colnames(X))
    Sigma_0_out <- Sigma / dispersion
    dimnames(Sigma_0_out) <- list(colnames(X), colnames(X))
  } else if (!is.null(dispersion_input)) {
    dispersion <- dispersion_input
    Sigma <- Sigma_0_out * dispersion
    dimnames(Sigma) <- list(colnames(X), colnames(X))
  }

  ## Step F: return calibrated outputs.
  list(
    dispersion = dispersion,
    shape = shape,
    rate = rate,
    Sigma = Sigma,
    Sigma_0 = Sigma_0_out
  )
}
