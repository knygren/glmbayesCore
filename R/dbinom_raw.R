#' Binomial density for non-integer counts
#'
#' Evaluates the binomial density through R's C-level \code{dbinom_raw()}
#' (Rmath), which, unlike \code{\link[stats]{dbinom}}, does not require
#' \code{x} and \code{n} to be integers. For integer \code{x} and \code{n} the
#' result matches \code{stats::dbinom(x, n, p, log)}.
#'
#' @param x Vector of success counts; may be non-integer.
#' @param n Vector of trial counts (\code{n >= 0}); may be non-integer.
#' @param p Vector of success probabilities in \eqn{[0, 1]}.
#' @param q Vector of failure probabilities, \code{1 - p} by default. Supplying
#'   \code{q} computed directly (for example \code{plogis(-eta)}) is more
#'   accurate when \code{p} is close to 1. Must satisfy
#'   \eqn{|p + q - 1| \le 10^{-10}}.
#' @param log Logical; if \code{TRUE}, return the log density.
#'
#' @details
#' Arguments are recycled to the length of the longest. For non-integer
#' \code{x} and \code{n} the density is the continuous extension of the
#' binomial probability mass function used internally by R (via Stirling-error
#' and deviance terms), as in \code{\link[stats]{dbeta}}.
#'
#' Conventions follow \code{stats::dbinom()}: \code{x < 0} or \code{x > n}
#' gives density 0 (\code{-Inf} on the log scale); \code{NA}/\code{NaN} inputs
#' propagate; invalid parameters (\code{n < 0}, \code{p} or \code{q} outside
#' \eqn{[0, 1]}, or \code{p + q} not equal to 1) give \code{NaN} with a
#' warning.
#'
#' @return Numeric vector of densities (log densities if \code{log = TRUE}),
#'   with length equal to the longest of \code{x}, \code{n}, \code{p} and
#'   \code{q}.
#'
#' @seealso \code{\link[stats]{dbinom}}
#'
#' @examples
#' dbinom_raw(3, 10, 0.4)
#' stats::dbinom(3, 10, 0.4)
#'
#' ## Non-integer counts, e.g. tempered weights (1 - rwt) * wt
#' dbinom_raw(0.75 * 3, 0.75 * 10, 0.4, log = TRUE)
#'
#' @export
dbinom_raw <- function(x, n, p, q = 1 - p, log = FALSE) {
  len <- max(length(x), length(n), length(p), length(q))
  if (len == 0L || min(length(x), length(n), length(p), length(q)) == 0L) {
    return(numeric(0))
  }
  x <- rep_len(as.numeric(x), len)
  n <- rep_len(as.numeric(n), len)
  p <- rep_len(as.numeric(p), len)
  q <- rep_len(as.numeric(q), len)
  out <- dbinom_raw_cpp_export(x, n, p, q, isTRUE(log))
  if (any(is.nan(out) & !is.na(x + n + p + q))) {
    warning("NaNs produced")
  }
  out
}
