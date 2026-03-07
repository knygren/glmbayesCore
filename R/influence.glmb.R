#' Bayesian Regression Diagnostics
#'
#' This function provides the basic quantities which are used in forming a wide variety of diagnostics for checking
#' the quality of Bayesian regression fits. These methods delegate to \code{\link[stats]{influence}},
#' \code{\link[stats]{cooks.distance}}, \code{\link[stats]{dfbetas}}, and related functions in the \pkg{stats}
#' package, applied to the fitted GLM component stored in the model object.
#'
#' @inheritParams stats::lm.influence
#' @return a \code{\link{list}} with components as returned by \code{\link[stats]{influence}}.
#' @details
#' Cook's distance was introduced by \insertCite{Cook1977}{glmbayes}. The dfbetas, dffits, and covratio
#' diagnostics follow the framework of \insertCite{BelsleyKuhWelsch1980}{glmbayes}. Because \code{glmb}
#' and \code{lmb} store coefficient draws rather than a single mode, these methods use the fitted
#' \code{fit} component (from the underlying \code{glm}/\code{lm} fit at the posterior mode) for
#' influence calculations.
#' @seealso \code{\link[stats]{influence}}, \code{\link[stats]{influence.measures}},
#' \code{\link[stats]{cooks.distance}}, \code{\link[stats]{dfbetas}}, \code{\link{summary.glmb}}
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_glmb.wfit.R
#' @method influence glmb
#' @export 


influence.glmb<-function(model,...){
  
  # Just tell function to use fit component
  # necessary because coefficients are draws and not modes
  # so not all items returned by fitting function can be include in glmb returned list
  
  return(influence(model$fit,...)) 
  
}


#' Bayesian Regression Diagnostics
#'
#' These functions compute Cook's distance, dfbetas, dffits, covratio, and standardized/studentized
#' residuals for fitted Bayesian GLM objects. They delegate to the corresponding \pkg{stats} functions
#' applied to the fitted GLM component.
#' @param infl influence structure as returned by \code{influence.glmb}
#' (the latter only for the glm method of \code{rstudent} and \code{cooks.distance}).
#' @inheritParams stats::influence.measures
#' @return a \code{\link{list}} with components as returned by the underlying \pkg{stats} functions.
#' @example inst/examples/Ex_glmb.wfit.R
#' @export 


glmb.influence.measures<-function (model, infl = influence(model)) 
{

  if(is.null(infl)) infl=influence(model)
  
  return(influence.measures(model$fit,infl))
}

#' @export 
#' @method rstandard glmb
#' @rdname glmb.influence.measures 

rstandard.glmb<-function(model,...,infl=influence(model)){

  if(is.null(infl)) infl=influence(model)
  
  return(rstandard(model$fit,infl))

}
  
#' @export 
#' @method rstudent glmb
#' @rdname glmb.influence.measures 

rstudent.glmb<-function(model,...,infl=influence(model)){
  
  if(is.null(infl)) infl=influence(model)
  
  return(rstudent(model$fit,infl))
  
}


#' @export 
#' @rdname glmb.influence.measures 

glmb.dffits<-function(model,infl=influence(model)){
  
  if(is.null(infl)) infl=influence(model)
  
  return(dffits(model$fit,infl))
  
}


#' @export 
#' @method dfbetas glmb
#' @rdname glmb.influence.measures 

dfbetas.glmb<-function(model,...,infl=influence(model)){
  
  if(is.null(infl)) infl=influence(model)
  
  return(dfbetas(model$fit,infl))
  
}

#' @export 
#' @rdname glmb.influence.measures 

glmb.covratio<-function(model,infl=influence(model)){
  
  if(is.null(infl)) infl=influence(model)
  
  return(covratio(model$fit,infl))
  
}


#' @export 
#' @method cooks.distance glmb
#' @rdname glmb.influence.measures 

cooks.distance.glmb<-function(model,...,infl=influence(model)){
  
  if(is.null(infl)) infl=influence(model)
  
  return(cooks.distance(model$fit,infl))
  
}

