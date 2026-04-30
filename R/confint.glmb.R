#' Credible Intervals for Model Parameters
#'
#' Computes Bayesian credible intervals for model parameters from posterior draws.
#' Intervals are constructed as quantiles of the posterior distribution, following
#' the standard approach for posterior summaries \insertCite{Gelman2013}{glmbayes}.
#'
#' @param object a fitted model object of class \code{"glmb"}. Typically the result of a call to \link{glmb}.
#' @param parm a specification (not yet implemented) of which parameters are to be given credible sets,
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the credible interval required.
#' @param \ldots additional argument(s) for methods.
#' @return A matrix (or vector) with columns giving lower and
#' upper credible limits for each parameter. These will be labeled
#' (1-level)/2 and 1-(1-level)/2 in \% (by default 2.5\% and 97.5\%).
#' @seealso \code{\link{summary.glmb}}, \code{\link{vcov.glmb}}, \code{\link{glmb}},
#'   \code{\link{glmbayes-package}}; \code{\link{rglmb}}, \code{\link{rlmb}}, \code{\link{lmb}};
#'   \code{\link[stats]{confint}} for classical confidence intervals
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_confint.glmb.R
#' @export 
#' @method confint glmb


## This method follows stats::confint() method conventions while reporting
## posterior credible intervals. See inst/COPYRIGHTS.
confint.glmb<-function(object,parm,level=0.95,...)
{
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  ci <- t(apply(object$coefficients, 2, FUN = quantile, probs = a))
  return(ci)
}
