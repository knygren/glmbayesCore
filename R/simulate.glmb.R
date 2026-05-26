#' Simulate Responses
#'
#' Simulate responses from the posterior predictive distribution corresponding to a
#' fitted \code{glmb} object \insertCite{Gelman2013}{glmbayes}.
#' @param object An object of class \code{glmb}, typically the result of a call to the 
#' function \code{glmb}.
#' @param nsim Defunct (see below).
#' @param seed an object specifying if and how the random number generator should be 
#' initialized (seeded).
#' @param \ldots Additional arguments passed to the function. Will 
#' frequently include a matrix pred of simulated predictions from
#' the predict function, the family (e.g., binomial) and an optional
#' vector of weights specifying prior.weights for the simulated values (default is 1)
#' @return Simulated values for data corresponding to simulated model predictions that correspond either
#' to the original data or to a \code{newdata} data frame provided to the predict function.
#' @seealso \code{\link{predict.glmb}}, \code{\link{glmb}}, \code{\link{glmbayes-package}};
#'   \code{\link{rglmb}}, \code{\link{rlmb}}, \code{\link{lmb}};
#'   \code{\link[stats]{simulate}} (e.g. \code{simulate.glm}, \code{simulate.lm} for classical fits);
#'   see \insertCite{glmbayesChapter06}{glmbayes} for model statistics.
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_residuals.glmb.R
#' @export
#' @method  simulate glmb

## This method follows stats::simulate() method conventions while using
## posterior predictive draws from glmb objects. See inst/COPYRIGHTS.
simulate.glmb<-function(object,nsim=1,seed=NULL,...){
  if (!is.null(seed)) set.seed(seed)

  family=object$family$family
  method_args=list(...)
  
  if(!is.null(method_args[['pred']])) { pred=method_args[['pred']]}
  else{pred=predict(object)}
  if(!is.null(method_args[['prior.weights']])) { wt=method_args[['prior.weights']]}
  else{ wt <- object$prior.weights }
  
  nvars=ncol(pred)
  nsims=nrow(pred)
  y_temp<-matrix(0,nrow=nrow(pred),ncol=ncol(pred))
  
  for(i in 1:nsims){
    
    if(family=="poisson") y_temp[i,1:nvars]=rpois(n=nvars,pred[i,1:nvars])
    if(family=="quasipoisson") y_temp[i,1:nvars]=rpois(n=nvars,pred[i,1:nvars])             
    ### Verify this part - rather complicated
    if(family=="Gamma") y_temp[i,1:nvars]=rgamma(n=nvars,shape=wt,(1/wt)*pred[i,1:nvars])              

    ## Simulate from binomial and then divide by the weight!    
    if(family=="binomial") y_temp[i,1:nvars]=rbinom(n=nvars,size=round(wt),prob=pred[i,1:nvars])/round(wt)              
    if(family=="quasibinomial") y_temp[i,1:nvars]=rbinom(n=nvars,size=round(wt),prob=pred[i,1:nvars])/round(wt)                            
    
    ## Verify this part - rather complicated    
    if(family=="gaussian") y_temp[i,1:nvars]=rnorm(n=nvars,mean=pred[i,1:nvars],sd=sqrt(1/wt))              
    
  }
  
  return(y_temp)
}
