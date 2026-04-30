#' Case and Variable Names of Fitted Models
#'
#' Simple utilities returning (non-missing) case names, and (non-eliminated) variable names.
#' @param object a fitted model object, typically of class \code{"glmb"}. 
#' @param full logical; if \code{TRUE}, all names (including zero weights, ...) are returned.
#' @param \ldots further arguments passed to or from other methods.
#' @return A character vector
#' @seealso \code{\link{glmb}}, \code{\link{glmbayes-package}}; \code{\link{rglmb}}, \code{\link{rlmb}},
#'   \code{\link{lmb}}; \code{\link[stats]{case.names}}.
#' @example inst/examples/Ex_case.names.glmb.R
#' @rdname case.names.glmb
#' @export
#' @method case.names glmb


## These methods follow stats model-object naming conventions. See
## inst/COPYRIGHTS for derived-code notices.
case.names.glmb<-function(object,full=FALSE,...){
  w=weights(object)
  dn=rownames(object$x)
  if(full || is.null(w)) dn else dn[w!=0]
}

#' @rdname case.names.glmb
#' @export
#' @method variable.names glmb

variable.names.glmb <- function(object, full = FALSE, ...)
{
  # Currently just returns results based on glm object
  return(colnames(object$coefficients))
}