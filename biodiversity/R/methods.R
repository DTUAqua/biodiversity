##' Log likelihood of biodiv object 
##' @method logLik biodiv 
##' @param  object biodiv fitted object (result from biodiv)
##' @param  ... extra arguments
##' @details ...
##' @export
logLik.biodiv<-function(object, ...){
  ret<- -object$opt$objective
  attr(ret,"df")<-length(object$opt$par)
  class(ret)<-"logLik"
  ret
}

##' Print biodiv object 
##' @method print biodiv 
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
print.biodiv<-function(x, ...){
  cat("Biodiversity model: log likelihood is", logLik.biodiv(x,...),"Convergence", ifelse(0==x$opt$convergence, "OK\n", "failed\n"))
}

##' Extract fixed coefficients of biodiversity model object 
##' @method coef biodiv 
##' @param  object biodiv fitted object (result from biodiv)
##' @param  ... extra arguments
##' @details ...
##' @importFrom stats coef
##' @export
coef.biodiv <- function(object, ...){
  ret <- object$sdrep$par.fixed
  attr(ret,"cov") <- object$sdrep$cov.fixed
  attr(ret,"sd") <- sqrt(diag(object$sdrep$cov.fixed))
  class(ret)<-"biodivcoef"
  ret
}

##' Extract residuals from biodiv object 
##' @method residuals biodiv 
##' @param object biodiv fitted object (result from biodiv)
##' @param ... extra arguments for TMB's oneStepPredict
##' @importFrom stats residuals
##' @importFrom TMB oneStepPredict
##' @details ...
##' @export
residuals.biodiv<-function(object, ...){
  cat("One-observation-ahead residuals. Total number of observations: ", nobs(object), "\n")
  ## ooa-residuals
  res <- oneStepPredict(object$obj,
                        method = "oneStepGeneric",
                        discrete = TRUE, 
                        observation.name = "nsp", 
                        data.term.indicator = "keep", discreteSupport = 2*max(object$data$nsp),...)
  cat("One-observation-ahead residuals. Done\n")  
  ret <- cbind(object$data, res)
  class(ret)<-"biodivres"
  ret
}

