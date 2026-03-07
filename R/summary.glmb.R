#' Summarizing Bayesian Generalized Linear Model Fits
#'
#' These functions are all \code{\link{methods}} for class \code{glmb} or \code{summary.glmb} objects.
#' 
#' @aliases 
#' summary.glmb
#' print.summary.glmb
#' @param object an object of class \code{"glmb"} for which a summary is desired.
#' @param x an object of class \code{"summary.glmb"} for which a printed output is desired.
#' @param digits the number of significant digits to use when printing.
#' @param \ldots Additional optional arguments
#' @return \code{summary.glmb} returns a object of class \code{"summary.glmb"}, a 
#' list with components: 
#' \item{call}{the component from \code{object}}
#' \item{n}{number of draws generated}
#' \item{residuals}{vector of mean deviance residuals}
#' \item{coefficients1}{Matrix with the prior mean and maximum likelihood coefficients with associated standard deviations}
#' \item{coefficients}{Matrix with columns for the posterior mode, posterior mean, posterior standard
#' deviation, monte carlo error, and tail probabilities (posterior probability of observing a 
#' value for the coefficient as extreme as the prior mean)}
#' \item{dir_tail}{List containing information related to the directional tail relative to the Prior}
#' \item{dir_tail_null}{List containing information related to the directional tail relative to the Null Model}
#' \item{Percentiles}{Matrix with estimated percentiles associated with the posterior density}
#' \item{pD}{Estimated effective number of parameters}
#' \item{deviance}{Vector with draws for the deviance}
#' \item{DIC}{Estimated DIC statistic}
#' \item{iters}{Average number of candidates per generated draws}
#' @details The \code{summary.glmb} function summarizes the output from the \code{glmb} function.
#' Key output includes mean residuals, information related to the prior, mean coefficients
#' with associated stats, percentiles for the coefficients, as well as the effective number of
#' parameters and the DIC statistic. The \code{dir_tail} component reports the directional tail
#' diagnostic; see \code{\link{directional_tail}} and \insertCite{glmbayesChapterA04}{glmbayes}
#' for interpretation.
#'
#' @seealso \code{\link{directional_tail}}, \code{\link{lmb}}, \code{\link{glmb}}, \code{\link{summary}},
#' \code{[stats]\link{summary.lm}},\code{[stats]\link{summary.glm}} 
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_summary.glmb.R
#' @export
#' @method summary glmb


summary.glmb<-function(object,...){
  
    
  #est.disp <- FALSE
  
  #res<-residuals(object)
  #nobs=length(object$y)
  #n=nrow(object$coefficients)
  #df.r <- length(object$y)-object$pD
  
  
  ##  
  
  dispersion=get_dispersion(object)
  
  dir_tail <- directional_tail(object)
  
  

  #if(!is.null(object$DIC)) DIC=object$DIC
  #else DIC=NA
  #print("Mean Dispersion")
  mres<-colMeans(residuals(object))
  

  # Ensure null_est is a full-length vector matching Prior$mean
  coef_names <- names(object$Prior$mean)
  null_est_full <- numeric(length(coef_names))
  names(null_est_full) <- coef_names
  
  

  
  l1<-length(object$coef.means)
  n<-length(object$coefficients[,1])
  percentiles<-matrix(0,nrow=l1,ncol=7)

  se<-sqrt(diag(var(object$coefficients)))

  #  if(object$family$family=="quasipoisson") se=se*sqrt(mean(dispersion))
  #  mc<-se/sqrt(object$n)
  mc<-se/sqrt(n)
  priorrank<-matrix(0,nrow=l1,ncol=1)
  pval1<-matrix(0,nrow=l1,ncol=1)
  pval2<-matrix(0,nrow=l1,ncol=1)
  
  pval_null_dir <- matrix(0, nrow = l1, ncol = 1)  # Directional tail probability (left tail) relative to null
  pval_null_tail <- matrix(0, nrow = l1, ncol = 1) # Two-sided tail probability relative to null
  
  ## Note: This restricts function to classes glmb and lmb
  ## Could break calls to this function from other classes

  

  if(!is.null(object$glm)){
    
    glmsummary<-summary(object$glm)
    se1<-sqrt(diag(glmsummary$cov.scaled))
    ml=as.numeric(object$glm$coefficients)


    lhs <- formula(object)[[2]]
    null_formula <- as.formula(paste(deparse(lhs), "~ 1"), env = environment(object$glm))
 
    lhs <- formula(object$glm)[[2]]
    null_formula <- as.formula(paste(deparse(lhs), "~ 1"), env = environment(object$glm))


    # Refit intercept-only GLM to get null estimates
    null_glm <- glm(
      formula = null_formula,
      data = object$data,
      family = object$glm$family
    )
    null_est <- coef(null_glm)
    

  } 


    
  ##
  

  if (!is.null(object$lm)) {
    lm_summary <- try(summary.lm(object$lm), silent = TRUE)
    
    if (inherits(lm_summary, "try-error")) {
      warning("summary.lm failed due to singular covariance; coefficients set to NA")
      
      se1 <- rep(NA_real_, length(coef(object$lm)))
      ml  <- coef(object$lm)
      
      # still build a null model if you want
      response_var <- all.vars(eval(object$lm$call$formula))[1]
      null_formula <- reformulate("1", response = response_var)
      null_lm <- lm(formula = null_formula, data = object$lm$model)
      null_est <- coef(null_lm)
      
    } else {
      se1 <- lm_summary$coefficients[, 2]
      ml  <- lm_summary$coefficients[, 1]
      
      response_var <- all.vars(eval(object$lm$call$formula))[1]
      null_formula <- reformulate("1", response = response_var)
      null_lm <- lm(formula = null_formula, data = object$lm$model)
      null_est <- coef(null_lm)
    }
  }


  # Fill intercept from null model, others set to 0
  intercept_name <- "(Intercept)"
  if (intercept_name %in% names(null_est)) {
    null_est_full[intercept_name] <- null_est[intercept_name]
  }  
  
  dir_tail_null <- directional_tail(object,mu0=null_est_full)
  
  

    
  for(i in 1:l1){
    percentiles[i,]<-quantile(object$coefficients[,i],probs=c(0.01,0.025,0.05,0.5,0.95,0.975,0.99))
    ##tail_center <- if (use_mode) object$coef.mode[i] else object$coef.means[i]
    #test <- append(object$coefficients[,i], tail_center)
    test<-append(object$coefficients[,i],object$Prior$mean[i])
    test2<-rank(test)
    priorrank[i,1]<-test2[n+1]
    priorrank[i,1]<-test2[n+1]
    #pval1[i,1]<-priorrank[i,1]/(n+1)
    #pval2[i,1]<-min(pval1[i,1],1-pval1[i,1])
    
    # Directional tail probability (left tail) 
    pval1[i,1] <- mean(object$coefficients[,i] < object$Prior$mean[i])
    pval2[i,1] <- min(pval1[i,1], 1 - pval1[i,1])
    
    # Directional tail probability relative to null mode
    pval_null_dir[i, 1] <- mean(object$coefficients[, i] < null_est_full[i])
    
    # Two-sided tail probability relative to null mode
    pval_null_tail[i, 1] <- min(pval_null_dir[i, 1], 1 - pval_null_dir[i, 1])
    
  }
  
  se_pval2 <- sqrt(pval2 * (1 - pval2) / n)
  se_pvalnull <- sqrt(pval_null_tail * (1 - pval_null_tail) / n)
  
  
  
    Tab1 <- cbind(
      "Null Mode"= as.numeric(null_est_full),
      "Prior Mean"=as.numeric(object$Prior$mean),
      "Prior.sd"= as.numeric(sqrt(diag(object$Prior$Variance))),
      "Max Like."= as.numeric(ml),
      "Like.sd"= as.numeric(se1)
    )
    rownames(Tab1) <- names(object$Prior$mean)
  TAB<-cbind("Post.Mode"=as.numeric(object$coef.mode),
             "Post.Mean"=as.numeric(object$coef.means),
             "Post.Sd"=as.numeric(se),
             "MC Error"=as.numeric(mc),             
##             "SE(Null_tail)"  = as.numeric(se_pvalnull),
             "Pr(Null_tail)"=as.numeric(pval_null_tail),  
             "SE(tail)"  = as.numeric(se_pval2),
             "Pr(Prior_tail)"=as.numeric(pval2)  
  )
  rownames(TAB) <- names(object$Prior$mean)
  
  TAB2<-cbind("1.0%"=percentiles[,1],"2.5%"=percentiles[,2],"5.0%"=percentiles[,3],Median=as.numeric(percentiles[,4]),"95.0%"=percentiles[,5],"97.5%"=as.numeric(percentiles[,6]),"99.0%"=as.numeric(percentiles[,7]))
  
  rownames(TAB2)<-rownames(TAB)
  
  res<-list(
    call=object$call,
    n=n,
    residuals=mres,
    coefficients1=Tab1,
    coefficients=TAB,
    dir_tail=dir_tail,
    dir_tail_null=dir_tail_null,
    Percentiles=TAB2,
    pD=object$pD,
    deviance=object$deviance,
    DIC=object$DIC,
    #          DIC=DIC,
    dispersion=mean(object$dispersion),
    #          dispersion=mean(dispersion),
    iters=mean(object$iters)
  )
  
  class(res)<-"summary.glmb"
  
  res
  
}

#' @rdname summary.glmb
#' @export
#' @method print summary.glmb

print.summary.glmb<-function(x,digits = max(3, getOption("digits") - 3),...){
  cat("Call\n")
  print(x$call)
  primary_class <- class(x)[1]
  
  if (primary_class == "lmb") {
    cat("\nExpected Residuals:\n")
  } else if (primary_class == "glmb") {
    fam <- tryCatch(x$family$family, error = function(e) NA)
    if (!is.na(fam) && fam == "gaussian") {
      cat("\nExpected Residuals:\n")
    } else {
      cat("\nExpected Deviance Residuals:\n")
    }
  } else {
    cat("\nExpected Residuals:\n")  # fallback for unknown class
  }
#  cat("\nExpected Deviance Residuals:\n")
  if (length(x$residuals) > 5) {
    fn   <- fivenum(x$residuals)
    names(fn) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(fn)
  } else {
    print(x$residuals)
  }
  cat("\nPrior and Maximum Likelihood Estimates with Standard Deviations\n\n")
  printCoefmat(x$coefficients1,digits=digits)
  cat("\nBayesian Estimates Based on",x$n,"iid draws\n\n")
  printCoefmat(x$coefficients,digits=digits,P.values=TRUE,has.Pvalue=TRUE)
  cat("  Directional Tail Summaries:\n\n")
  dir_table <- data.frame(
    Metric = c("Mahalanobis Distance", "Tail Probability"),
    `vs Null` = c(
      formatC(x$dir_tail_null$mahalanobis_shift, digits = 4, format = "f"),
      formatC(x$dir_tail_null$p_directional, digits = 4, format = "f")
    ),
    `vs Prior` = c(
      formatC(x$dir_tail$mahalanobis_shift, digits = 4, format = "f"),
      formatC(x$dir_tail$p_directional, digits = 4, format = "f")
    ),
    check.names = FALSE
  )
  
  print(dir_table, row.names = FALSE)
  cat("  [Tail probabilities are P(delta^T * Z <= 0) in whitened space]\n\n")
  cat("\nDistribution Percentiles\n\n")
  printCoefmat(x$Percentiles,digits=digits)
  cat("\nEffective Number of Parameters:",x$pD,"\n")
  cat("Expected Residual Deviance:",mean(x$deviance),"\n")
  cat("DIC:",x$DIC,"\n\n")
  cat("Expected Mean dispersion:",x$dispersion,"\n")
  cat("Sq.root of Expected Mean dispersion:",sqrt(x$dispersion),"\n\n")
  cat("Mean Likelihood Subgradient Candidates Per iid sample:",x$iters,"\n\n")
  
  
}

# Helpers --------------------------------------------------------------------


get_dispersion<-function(object){
  
  df.r <- length(object$y)-object$pD
  
  
  if(!is.null(object$dispersion)) object$dispersion
  else{
    
    if(object$family$family=="quasipoisson"){
      n=nrow(object$coefficients)
      disp_temp=rep(0,n)
      m=length(object$y)  
      k=ncol(object$x)    
      res_temp=matrix(0,nrow=n,ncol=m)
      fit_temp=object$x%*%t(object$coefficients)
      for(l in 1:n){
        if(is.null(object$fit$offset))         fit_temp[1:m,l]=exp(fit_temp[1:m,l])
        else fit_temp[1:m,l]=exp(object$fit$offset+fit_temp[1:m,l])
        
        res_temp[l,1:m]=(object$y-fit_temp[1:m,l])
        disp_temp[l]=(1/(m-k))*sum(res_temp[l,1:m]^2*object$prior.weights/fit_temp[1:m,l])
  
      }
    
     return(disp_temp) 
    }
    
  }
  
}

