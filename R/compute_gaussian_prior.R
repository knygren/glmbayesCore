#' Compute Gaussian Normal-Gamma Calibration Components
#'
#' Gaussian Normal-Gamma calibration function used by \code{\link{Prior_Setup}}.
#' Given Gaussian-model sufficient inputs and a dispersion-independent
#' coefficient prior covariance \eqn{\Sigma_0} (Chapter 11 framing), this
#' function computes calibrated Gaussian prior quantities:
#' \code{dispersion}, \code{shape}, \code{rate}, and \code{Sigma}.
#'
#' The function is structured as a step-wise pipeline:
#' \enumerate{
#'   \item Validate dimensions and numeric inputs.
#'   \item Compute weighted residual sum of squares at \code{bhat}.
#'   \item Build \eqn{X'WX}, invert it, and construct \eqn{S_{marg}}.
#'   \item Map \code{shape_df} and \code{n_prior} to Gamma shape.
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
#'   \item \eqn{b_0 = \frac{1}{2}\frac{n_{prior}}{n_{effective}}S_{marg}}
#'   \item \eqn{E[\sigma^2 \mid y] = S_{marg}\frac{n_{effective}+n_{prior}}{n_{effective}(n_{effective}+n_{prior}-2)}}
#'   \item \eqn{\Sigma = \frac{n_{effective}}{n_{prior}} E[\sigma^2 \mid y] (X^T W X)^{-1}}
#' }
#' with Gamma shape controlled by \code{shape_df}.
#'
#' Limiting behavior:
#' \itemize{
#'   \item As \eqn{n_{prior} \to \infty}, prior information dominates and the
#'   Gamma prior on precision becomes increasingly concentrated.
#'   \item The expected coefficient location remains centered at \eqn{\mu}; this
#'   function does not change the mean vector, only scale components.
#'   \item As \eqn{n_{prior}\to 0^+}, the returned expected dispersion tends to
#'   \eqn{S_{marg}/(n_{effective}-2)} from the formula
#'   \eqn{E[\sigma^2|y]=S_{marg}(n_{effective}+n_{prior})/[n_{effective}(n_{effective}+n_{prior}-2)]},
#'   i.e. a finite data-dominated limit when \eqn{n_{effective}>2}.
#'   \item As \eqn{n_{prior}\to 0^+}, \eqn{b_0\to 0} and shape terms can become
#'   very small depending on \code{shape_df}; this corresponds to a very diffuse
#'   precision prior and may be numerically delicate near the boundary.
#'   \item In the Chapter 11 scalar-\code{pwt} special-case path
#'   (\eqn{pwt=n_{prior}/(n_{prior}+n_w)} with fixed \eqn{n_w>2}), one has
#'   \eqn{pwt\to 0}, \eqn{S_{marg}\to RSS}, \eqn{b_n\to RSS/2},
#'   \eqn{a_n-1\to (n_w-2)/2}, and \eqn{V_n\to (X^TWX)^{-1}. Therefore}
#'   \deqn{
#'   \mathrm{Cov}(\beta\mid y)\to \frac{RSS}{n_w-2}(X^TWX)^{-1},
#'   }
#'   matching the classical weighted least-squares covariance under the
#'   Gaussian \code{glm.fit} dispersion convention used in \code{Prior_Setup}.
#'   \item The \code{shape_df} choice adjusts how coefficient dimension \eqn{p}
#'   enters the precision shape (\code{"n_prior"}, \code{"n_prior+p"}, \code{"n_prior-p"}).
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
#' @param n_prior Positive scalar effective prior sample size.
#' @param shape_df Character string controlling the Gamma shape numerator:
#'   \code{"n_prior"}, \code{"n_prior+p"}, or \code{"n_prior-p"}.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{dispersion}: calibrated Gaussian dispersion scalar.
#'   \item \code{shape}: Gamma shape for residual precision.
#'   \item \code{rate}: Gamma rate for residual precision.
#'   \item \code{Sigma}: calibrated coefficient prior covariance matrix.
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
    n_prior,
    shape_df = c("n_prior", "n_prior+p", "n_prior-p")
) {
  ## ---------------------------------------------------------------------------
  ## Gaussian calibration pipeline:
  ## Step A: Validate inputs and dimensions.
  ## Step B: Compute weighted RSS from (Y, X, bhat, offset, weights).
  ## Step C: Build Gram terms and S_marg using Sigma_0 + (X'WX)^{-1}.
  ## Step D: Build Gamma shape from n_prior + shape_df rule.
  ## Step E: Calibrate dispersion/rate and map to coefficient Sigma.
  ## Step F: Return calibrated terms.
  ## ---------------------------------------------------------------------------
  shape_df <- match.arg(shape_df)
  if (!is.null(dispersion)) {
    if (!is.numeric(dispersion) || length(dispersion) != 1L ||
        !is.finite(dispersion) || dispersion <= 0) {
      stop("compute_gaussian_prior: dispersion must be NULL or a single positive finite numeric value.", call. = FALSE)
    }
  }
  dispersion_input <- dispersion

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
  if (n_effective <= 2) {
    stop("compute_gaussian_prior: require n_effective > 2 for Gaussian dispersion (denominator n_effective - 2).", call. = FALSE)
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

  ## Step D: Gamma shape via shape_df rule.
  n_shape_num <- switch(
    shape_df,
    "n_prior"   = n_prior,
    "n_prior+p" = n_prior + p,
    "n_prior-p" = {
      if (!is.finite(n_prior) || !is.finite(p) || n_prior <= p) {
        stop(
          "compute_gaussian_prior: shape_df = \"n_prior-p\" requires n_prior > p (number of coefficients). ",
          "Got n_prior = ", n_prior, " and p = ", p, ".",
          call. = FALSE
        )
      }
      n_prior - p
    }
  )
  shape <- n_shape_num / 2
  if (!is.finite(shape) || shape <= 0) {
    stop("compute_gaussian_prior: computed shape must be strictly positive.", call. = FALSE)
  }

  ## Step E: calibrate Gaussian dispersion/rate and implied Sigma.
  b_0_S_marg_formula <- 0.5 * (n_prior / n_effective) * S_marg
  den_phi <- n_prior + n_effective - 2
  E_phi_sigma2_special <- NA_real_
  if (is.finite(den_phi) && den_phi > 0) {
    E_phi_sigma2_special <- S_marg * (n_effective + n_prior) / n_effective / den_phi
  }

  if (!is.finite(E_phi_sigma2_special) || E_phi_sigma2_special <= 0) {
    stop("compute_gaussian_prior: E[sigma^2|y] (special) is missing or not positive.", call. = FALSE)
  }
  if (!is.finite(b_0_S_marg_formula) || b_0_S_marg_formula <= 0) {
    stop("compute_gaussian_prior: b_0 (special) is missing or not positive.", call. = FALSE)
  }
  dispersion <- E_phi_sigma2_special
  rate <- b_0_S_marg_formula
  Sigma <- (n_effective / n_prior) * dispersion * Ginv
  dimnames(Sigma) <- list(colnames(X), colnames(X))

  ## Step F: return calibrated outputs.
  list(
    dispersion = dispersion,
    shape = shape,
    rate = rate,
    Sigma = Sigma
  )
}
