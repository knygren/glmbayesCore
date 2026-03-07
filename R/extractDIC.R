#' Extract DIC from a Fitted Bayesian Model
#'
#' Computes the Deviance Information Criterion (DIC) and the effective
#' number of parameters for Bayesian generalized linear models fitted
#' via \code{glmb()} or \code{rglmb()}. The DIC, introduced by
#' \insertCite{Spiegelhalter2002}{glmbayes}, provides a Bayesian analog to the
#' AIC \insertCite{Akaike1974}{glmbayes}.
#'
#' @param fit A fitted model of class \code{"glmb"} or \code{"rglmb"}.
#' @param ... Additional arguments passed to or from methods.
#'
#' @return A named numeric vector with components:
#'   \describe{
#'     \item{pD}{Estimated effective number of parameters}
#'     \item{DIC}{Deviance Information Criterion}
#'   }
#'
#' @seealso \code{\link{summary.glmb}}, \code{\link{glmb}}, \code{\link{lmb}}
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_extractAIC.glmb.R
#'
#' @rdname extractDIC
#' @export
#' @method extractAIC glmb
extractAIC.glmb <- function(fit, ...) {
  c(pD = fit$pD, DIC = fit$DIC)
}

#' @rdname extractDIC
#' @export
#' @method extractAIC rglmb
extractAIC.rglmb <- function(fit, ...) {
  fit2 <- summary(fit)
  extractAIC(fit2, ...)
}

#' @rdname extractDIC
#' @export
extractDIC <- function(fit, ...) UseMethod(extractAIC)