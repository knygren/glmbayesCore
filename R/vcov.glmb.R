#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' Returns the posterior variance-covariance matrix of the regression coefficients
#' from a fitted Bayesian GLM object \insertCite{Gelman2013}{glmbayes}.
#'
#' @param object fitted model object, typically the result of a call to \code{\link{glmb}}.
#' @param \ldots additional arguments for method functions.
#' @return A matrix of estimated covariances between the parameter estimates
#' in the linear or non-linear predictor of the model. This should have
#' row and column names corresponding to the parameter names given by the
#' \code{\link{coef}} method.
#' @seealso \code{\link{confint.glmb}}, \code{\link{summary.glmb}}, \code{\link{glmb}},
#' \code{\link[stats]{vcov}}
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_vcov.glmb.R
#' @export
#' @method vcov glmb

vcov.glmb<-function(object,...)
{
  return(cov(object$coefficients))
  
}
