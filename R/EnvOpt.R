#' Optimizes Envelope function for simulation
#'
#' Optimizes the size of the grid to try to limit the combined time of
#' the envelope construction and the simulation phase.
#' @param a1 Diagonal elements of data precision matrix for a model in standard form
#' @param n  Number of draws to generate
#' @param core_cnt Integer; number of OpenCL cores or parallel workers available.
#'   Defaults to 1 for backward compatibility. When >1, the function treats
#'   envelope build cost as reduced by a factor of core_cnt because construction
#'   can be parallelized.
#' @details This function attempts to find a computationally optimal 
#' gridsize by using information on the strength of the prior and the 
#' number of iterations desired. Generally, more data (i.e., larger values 
#' for the diagonal elements of the precision matrix) will require a larger grid.
#' The same also holds when the number of desired draws is higher 
#' (as the setup costs associated with the larger grid is offset by the 
#' savings in the number of candidates per sample).
#' @return A vector containing information on how many component each 
#' dimension should be split into.
#' @seealso \code{\link{rglmb}}, \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeSort}} 
#' @example inst/examples/Ex_EnvelopeOpt.R
#' @export
#' @keywords internal


EnvelopeOpt<-function(a1,n,core_cnt=1L){
  
  core_cnt <- as.integer(core_cnt)
  if (is.na(core_cnt) || core_cnt < 1L) core_cnt <- 1L
  
  
  a1rank<-rank(1/(1+a1))
  l1<-length(a1)
  
  dimcount<-matrix(0,(l1+1),l1)
  scaleest<-matrix(0,(l1+1),l1)
  intest<-c(1:(l1+1))
  slopeest<-c(1:(l1+1))
  
  dimcount[1,]<-diag(diag(l1))
  scaleest[1,]<-sqrt(1+a1)
  slopeest[1]<-prod(scaleest[1,])
  
  for(i in 2:(l1+1)){
    dimcount[i,]<-dimcount[i-1,]
    scaleest[i,]<-scaleest[i-1,]
    for(j in 1:l1){
      if(a1rank[j]==i-1){ 
        dimcount[i,j]<-3
        scaleest[i,j]<-2/sqrt(pi) 
      }
    }
##    intest[i]<-3^(i-1)
    intest[i]<-(3^(i-1))
    slopeest[i]<-prod(scaleest[i,])
  }
  evalest<-(intest/core_cnt)+n*slopeest
  minindex<-0
  for(j in 1:(l1+1)){if(evalest[j]==min(evalest)){minindex<-j}}
  
  message("Estimated draws per Acceptance: ", slopeest[minindex])
  
  dimcount[minindex,]
  
}