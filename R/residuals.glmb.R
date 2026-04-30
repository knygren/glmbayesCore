#' Accessing Bayesian Generalized Linear Model Fits
#'
#' Extract deviance residuals from fitted Bayesian GLM objects. The residuals
#' use the family's deviance residuals function as in \code{\link[stats]{residuals.glm}}
#' \insertCite{McCullagh1989}{glmbayes}.
#'
#' These functions are all \link{methods} for class \code{glmb}, \code{lmb}, or \code{summary.glmb} objects.
#' @param object an object of class \code{glmb}, typically the result of a call to \link{glmb}
#' @param ysim Optional simulated data for the data y.
#' @param \ldots further arguments to or from other methods
#' @return A matrix \code{DevRes} of dimension \code{n} times \code{p} containing
#' the Deviance residuals for each draw. If ysim is provided, the residuals are based
#' on a comparison to the simulated data instead. The credible intervals
#' for residuals based on simulated data should be a more appropriate measure of
#' whether individual residuals represent outliers or not.
#' @seealso \code{\link{predict.glmb}}, \code{\link{summary.glmb}}, \code{\link{glmb}},
#'   \code{\link{glmbayes-package}}; \code{\link{rglmb}}, \code{\link{rlmb}}, \code{\link{lmb}};
#'   \code{\link[stats]{residuals.glm}}
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_residuals.glmb.R
#' @export 
#' @method residuals glmb

## This method follows stats::residuals.glm() conventions while computing
## residuals across posterior draws. See inst/COPYRIGHTS.
residuals.glmb<-function(object,ysim=NULL,...)
{
  y<-object$y	
  n<-length(object$coefficients[,1])
  
  ## Updated to use prior.weights - likely matters for binomial data
  ## Need to verify this performs as expected
  wts <- object$prior.weights
  
  
  fitted.values<-object$fitted.values
  dev.residuals<-object$family$dev.resids
  DevRes<-matrix(0,nrow=n,ncol=length(y))
  
  for(i in 1:n)
  {
    if(is.null(ysim))    DevRes[i,]<-sign(y-fitted.values[i,])*sqrt(dev.residuals(y,fitted.values[i,],wts))
    else(DevRes[i,]<-sign(ysim[i,]-fitted.values[i,])*sqrt(dev.residuals(ysim[i,],fitted.values[i,],wts)))
  }
  
  colnames(DevRes)<-names(y)
  DevRes
}


#' @rdname residuals.glmb
#' @export 
#' @method residuals rglmb

residuals.rglmb <- function(object, ysim = NULL, ...) {
  y   <- object$y
  # n simulations (or rows in the posterior coef matrix)
  n   <- nrow(object$coefficients)
  wts <- object$prior.weights
  
  # 1) build a matrix (n * length(y)) of linear predictors and fitted values
  lp_mat <- t(object$x %*% t(object$coefficients))
  fv_mat <- object$family$linkinv(lp_mat)
  
  # 2) grab the family's deviance-resids function
  devfun <- object$family$dev.resids
  
  # 3) allocate
  DevRes <- matrix(0, nrow = n, ncol = length(y))
  
  # 4) fill it
  for (i in seq_len(n)) {
    if (is.null(ysim)) {
      mu_vec <- fv_mat[i, ]
    } else {
      mu_vec <- ysim[i, ]
    }
    
    # call the C-level deviance-resids with exactly (y, mu, wts)
    DevRes[i, ] <- sign(y - mu_vec) * sqrt(devfun(y, mu_vec, wts))
  }
  
  colnames(DevRes) <- names(y)
  DevRes
}


#' @rdname residuals.glmb
#' @export 
#' @method residuals lmb

residuals.lmb<-function(object,ysim=NULL,...)
{
  return(residuals.lm(object,ysim,...))
  }

