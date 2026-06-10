## Total-variation bounds for the two-block Gibbs sampler (Nygren 2020,
## Theorem 3 and Corollary 1), evaluated from the Remark 8 eigenvalue
## spectrum computed by two_block_rate().
##
## With the chain started at the exact posterior mean (as lmerb does via
## lmerb_posterior_mean), the mean term of both bounds is identically zero
## and only the variance-convergence sum remains:
##
##   ||Q^(l) - pi||_TV  <=  sum_{i=1}^n d_i^(l),
##   r_i^(l) = (1 - a_{i-1}^{2l}) / (1 - a_i^{2l}),  a_0 = 0,
##
## with the eigenvalues a_i in ASCENDING order (so that k_i^(l) =
## 1/(1 - a_i^{2l}) is nondecreasing, as Lemma 2 requires).
##
## The generalized error function has the closed form (paper Remark 3)
## erf_n(x) = P(chi_n <= sqrt(2) x) = pchisq(2 x^2, df = n), which makes the
## exact Theorem 3 terms elementary.  The Corollary 1 envelope replaces each
## d_i^(l) with the relaxation of Remarks 5 and 17:
##
##   d_i^(l) <= [a_i^{2l} / (sqrt(1-a_i^2) sqrt(1-a_{i-1}^2))]
##              * sqrt((n+1-i)/2) * c_{n-i},
##   c_m = 2 e^{-m/2} (m/2)^{m/2} / Gamma((m+1)/2),
##
## (the sqrt((n+1-i)/2) factor follows from the Remark 5 derivation; the
## Corollary 1 display omits it, but it is required for the chain
## d_i <= [x1-x2] c_{n-i} with x1-x2 <= [sqrt(r)-sqrt(1/r)] sqrt((n+1-i)/2)).
## The optional mean term uses D0 = (x0-mu)' Sigma11^{-1} (x0-mu):
## exact erf_1 for "theorem3", the linear envelope lambda*^l sqrt(D0/(2 pi))
## for "corollary1".

#' Generalized n-dimensional error function
#'
#' \code{erf_n(x) = pchisq(2 x^2, df = n)} (Nygren 2020, Remark 3; Brown
#' 1963).  \code{n = 1} reduces to the classical error function.
#'
#' @param x Non-negative numeric vector.
#' @param n Dimension (positive integer).
#' @return Numeric vector of probabilities.
#' @keywords internal
.two_block_erfn <- function(x, n) {
  stats::pchisq(2 * x^2, df = n)
}

#' @keywords internal
.two_block_tv_bound_one <- function(a_asc, l, method, D0, lambda_star) {
  n <- length(a_asc)
  ## a_i^{2l} computed in log space; a = 0 -> 0.
  u <- ifelse(a_asc > 0, exp(2 * l * log(a_asc)), 0)
  u_prev <- c(0, u[-n])

  if (identical(method, "theorem3")) {
    d <- numeric(n)
    for (i in seq_len(n)) {
      m_i <- n + 1L - i
      num <- u[i] - u_prev[i]
      den <- 1 - u[i]
      if (num <= 0 || den <= 0) next
      rm1 <- num / den                  # r - 1
      lr <- log1p(rm1)                  # log(r)
      if (!is.finite(rm1) || rm1 <= 0) next
      t_hi <- m_i * lr * (1 + rm1) / rm1   # m * ln(r) * r/(r-1)
      t_lo <- m_i * lr / rm1               # m * ln(r) / (r-1)
      d[i] <- stats::pchisq(t_hi, df = m_i) - stats::pchisq(t_lo, df = m_i)
    }
    var_term <- sum(d)
    mean_term <- if (D0 > 0) {
      .two_block_erfn(0.5 * lambda_star^l * sqrt(D0) / sqrt(2), 1L)
    } else 0
  } else {
    ## Corollary 1 envelope (Remarks 5 + 17)
    s_i <- sqrt(1 - a_asc^2)
    s_prev <- c(1, s_i[-n])
    d <- numeric(n)
    for (i in seq_len(n)) {
      if (u[i] <= 0) next
      m_i <- n + 1L - i
      m_e <- m_i - 1L                    # erf order exponent (n - i)
      c_m <- if (m_e == 0L) {
        2 / gamma(0.5)                   # = 2/sqrt(pi)
      } else {
        2 * exp(-m_e / 2) * (m_e / 2)^(m_e / 2) / gamma((m_e + 1) / 2)
      }
      d[i] <- u[i] / (s_i[i] * s_prev[i]) * sqrt(m_i / 2) * c_m
    }
    var_term <- sum(d)
    mean_term <- if (D0 > 0) {
      lambda_star^l * sqrt(D0 / (2 * pi))
    } else 0
  }

  min(var_term + mean_term, 1)
}

#' Total-variation bound for the two-block Gibbs sampler
#'
#' Evaluates the bound on the total variation distance between the
#' \eqn{l}-step kernel of the two-block Gibbs sampler and its target
#' (Nygren 2020), from the eigenvalue spectrum computed by
#' \code{\link{two_block_rate}}.
#'
#' \code{method = "theorem3"} evaluates the exact Theorem 3 terms
#' \eqn{d_i^{(l)}} using the closed form
#' \eqn{\mathrm{erf}_n(x) = P(\chi^2_n \le 2x^2)} with
#' \eqn{r_i^{(l)} = (1 - a_{i-1}^{2l})/(1 - a_i^{2l})}.
#' \code{method = "corollary1"} evaluates the looser geometric envelope of
#' Corollary 1 (via Remarks 5 and 17), which decays like
#' \eqn{a_i^{2l}} with explicit constants.
#'
#' When the chain is started at the exact posterior mean (as
#' \code{lmerb} does), \code{D0 = 0} and the mean term of both bounds
#' vanishes identically; only the variance-convergence sum remains.  The
#' returned bound is capped at 1.
#'
#' Note the bound applies to the block updated \emph{second} in each sweep
#' (the Block 2 hyper vector \eqn{\gamma}); the stored Block 1 draw lags by
#' a half-step, so evaluate at \code{l - 1} when calibrating
#' \code{m_convergence} for the random-effect draws.
#'
#' @param rate Object from \code{\link{two_block_rate}}.
#' @param l Integer vector of sweep counts (each \code{>= 1}).
#' @param method \code{"theorem3"} (exact terms) or \code{"corollary1"}
#'   (geometric envelope).
#' @param D0 Optional squared standardized distance of the starting point
#'   from the posterior mean,
#'   \eqn{(x^{(0)}-\mu)^\top \Sigma_{11}^{-1} (x^{(0)}-\mu)}.  Default 0
#'   (start at the posterior mean).
#' @return Numeric vector of TV bounds, one per element of \code{l}, capped
#'   at 1.
#' @references Nygren, K. (2020). \emph{On the total variation distance
#'   between multivariate normal densities with applications to two-block
#'   Gibbs samplers.} Unpublished manuscript.
#' @family simfuncs
#' @seealso \code{\link{two_block_rate}}, \code{\link{two_block_l_for_tv}}
#' @export
two_block_tv_bound <- function(rate,
                               l,
                               method = c("theorem3", "corollary1"),
                               D0 = 0) {
  if (!inherits(rate, "two_block_rate")) {
    stop("'rate' must be a two_block_rate object.", call. = FALSE)
  }
  method <- match.arg(method)
  l <- as.integer(l)
  if (length(l) < 1L || any(!is.finite(l)) || any(l < 1L)) {
    stop("'l' must contain integers >= 1.", call. = FALSE)
  }
  if (!is.numeric(D0) || length(D0) != 1L || D0 < 0) {
    stop("'D0' must be a single non-negative number.", call. = FALSE)
  }

  a_asc <- sort(rate$eigenvalues)   # Lemma 2 requires ascending order
  vapply(
    l,
    function(li) {
      .two_block_tv_bound_one(a_asc, li, method, D0, rate$lambda_star)
    },
    numeric(1L)
  )
}

#' Sweeps required to reach a TV tolerance
#'
#' Smallest \code{l} such that
#' \code{two_block_tv_bound(rate, l, method, D0) <= tol}.  The bound is
#' decreasing in \code{l}, so a doubling search followed by bisection is
#' exact.
#'
#' @inheritParams two_block_tv_bound
#' @param tol Target total-variation tolerance in (0, 1).
#' @param l_max Search cap (error if the bound stays above \code{tol}).
#' @return Integer: the required number of sweeps.
#' @family simfuncs
#' @seealso \code{\link{two_block_tv_bound}}
#' @export
two_block_l_for_tv <- function(rate,
                               tol,
                               method = c("theorem3", "corollary1"),
                               D0 = 0,
                               l_max = 1000000L) {
  if (!is.numeric(tol) || length(tol) != 1L || tol <= 0 || tol >= 1) {
    stop("'tol' must be a single value in (0, 1).", call. = FALSE)
  }
  method <- match.arg(method)
  bnd <- function(li) two_block_tv_bound(rate, li, method = method, D0 = D0)

  if (bnd(1L) <= tol) return(1L)
  lo <- 1L
  hi <- 2L
  while (bnd(hi) > tol) {
    lo <- hi
    if (hi >= l_max) {
      stop("bound does not reach tol = ", tol, " within l_max = ", l_max,
           " sweeps.", call. = FALSE)
    }
    hi <- min(2L * hi, as.integer(l_max))
  }
  while (hi - lo > 1L) {
    mid <- lo + (hi - lo) %/% 2L
    if (bnd(mid) <= tol) hi <- mid else lo <- mid
  }
  hi
}
