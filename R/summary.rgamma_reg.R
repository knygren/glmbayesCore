#' Summarizing Bayesian gamma_reg Distribution Functions
#'
#' These functions are all \code{\link{methods}} for class \code{rgamna_reg} or \code{summary.rgamma_reg} objects.
#' 
#' @aliases 
#' summary.rgamma_reg
#' print.summary.rgamma_reg
#' @param object an object of class \code{"rgamma_reg"} for which a 
#' summary is desired.
#' @param x an object of class \code{"summary.rgamma_reg"} for which a printed output is desired.
#' @param digits the number of significant digits to use when printing.
#' @param \ldots Additional optional arguments
#' @export
#' @rdname summary.rgamma_reg
#' @order 1
#' @method summary rGamma_reg


summary.rGamma_reg<-function(object,...){
  
  n<-length(object$dispersion)  
  percentiles<-matrix(0,nrow=1,ncol=7)
  me=mean(object$dispersion)
  se<-sqrt(var(object$dispersion))
  mc<-se/n
  Priorwt<-(se/(sqrt(object$Prior$shape)/object$Prior$rate))^2
  percentiles[1,]<-quantile(object$dispersion,probs=c(0.01,0.025,0.05,0.5,0.95,0.975,0.99))
  test<-append(object$dispersion,object$Prior$shape/object$Prior$rate)
  test2<-rank(test)
  priorrank<-test2[n+1]
  pval1<-priorrank/(n+1)
  pval2<-min(pval1,1-pval1)
  
  shape <- object$Prior$shape
  rate  <- object$Prior$rate
  
  # Prior mean of dispersion
  if (shape > 1) {
    prior_mean_disp <- rate / (shape - 1)
  } else {
    prior_mean_disp <- NA
  }
  
  # Prior SD of dispersion
  if (shape > 2) {
    prior_sd_disp <- rate / ((shape - 1) * sqrt(shape - 2))
  } else {
    prior_sd_disp <- NA
  }
  
  Tab1 <- cbind(
    "Prior.Mean" = prior_mean_disp,
    "Prior.Sd"   = prior_sd_disp,
    "Approx.Prior.wt" = Priorwt
  )
  
  rownames(Tab1) <- "dispersion"
  colnames(Tab1) <- c("Prior.Mean", "Prior.Sd", "Approx.Prior.wt")
  
  
  
  
  #  Tab1<-cbind("Prior.Mean"=object$Prior$shape/object$Prior$rate,"Prior.Sd"=sqrt(object$Prior$shape)/object$Prior$rate
  #              ,"Approx.Prior.wt"=Priorwt  )
  
  TAB<-cbind(
    #"Post.Mode"=as.numeric(object$PostMode),
    "Post.Mean"=me,
    "Post.Sd"=se,
    "MC Error"=as.numeric(mc)
    ,"Pr(tail)"=as.numeric(pval2)
  )
  TAB2<-cbind("1.0%"=percentiles[,1],"2.5%"=percentiles[,2],"5.0%"=percentiles[,3],Median=as.numeric(percentiles[,4]),"95.0%"=percentiles[,5],"97.5%"=as.numeric(percentiles[,6]),"99.0%"=as.numeric(percentiles[,7]))
  
  rownames(TAB)=c("dispersion")
  rownames(Tab1)=c("dispersion")
  rownames(TAB2)=c("dispersion")
  
  colnames(Tab1) <- c("Prior.Mean", "Prior.Sd", "Approx.Prior.wt")
  colnames(TAB) <- c("Post.Mean", "Post.Sd", "MC.Error", "Pr(tail)")
  
  res<-list(call=object$call,
            n=n,
            coefficients1=Tab1,
            coefficients=TAB,
            Percentiles=TAB2
  )
  
  # Reuse summary.rglmb class
  
  class(res) <- "summary.rGamma_reg"  
  res
  
}


#' @export
#' @rdname summary.rgamma_reg
#' @order 2
#' @method print summary.rGamma_reg


print.summary.rGamma_reg <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  
  ## --- Call ---
  cat("Call\n")
  print(x$call)
  cat("\n")
  
  ## --- Prior Estimates ---
  cat("Prior Estimates with Standard Deviations\n\n")
  print(round(x$coefficients1, digits))
  cat("\n")
  
  ## --- Posterior Estimates ---
  cat("Bayesian Estimates Based on", x$n, "iid draws\n\n")
  
  # Extract posterior table
  TAB <- round(x$coefficients, digits)
  
  # Compute significance stars
  pvals <- x$coefficients[, "Pr(tail)"]
  stars <- ifelse(pvals < 0.001, "***",
                  ifelse(pvals < 0.01,  "**",
                         ifelse(pvals < 0.05,  "*",
                                ifelse(pvals < 0.1,   ".", " "))))
  
  # Build final table with stars
  TAB2 <- cbind(TAB, Signif = stars)
  
  print(TAB2)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  ## --- Percentiles ---
  cat("Distribution Percentiles\n\n")
  print(round(x$Percentiles, digits))
  cat("\n")
  
  ## --- Dispersion summary ---
  disp.mean <- x$coefficients["dispersion", "Post.Mean"]
  cat("Expected Mean dispersion:", round(disp.mean, digits), "\n")
  cat("Sq.root of Expected Mean dispersion:", round(sqrt(disp.mean), digits), "\n\n")
  
  invisible(x)
}

