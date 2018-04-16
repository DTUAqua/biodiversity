##' Fit biodiversity model 
##' @param data A data set in same format as the species data set supplied with this package 
##' @param conf A code identifying the model configuration (-2=constant model, -1=saturated model, 1=Neutral, 2=Best, 3=Latitude, and 4=Metabolic)
##' @param fixK optionally a value used to fix the k parameter.
##' @param run if FALSE return AD object without running the optimization
##' @param ... extra arguments to MakeADFun
##' @return an object of class \code{biodiv}
##' @details The model ...
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats nlminb optimHess
##' @useDynLib biodiversity
##' @export
##' @examples
##' data(species)
##' fit <- biodiv(species, conf=1)
biodiv <- function(data, conf, fixK=NULL, run=TRUE, ...){
  origdata<-data
  #dat<-data[(data$density*data$asurv)>1.0e-6,]
  dat<-data  
  data <- list()
  data$code <- conf
  data$lat <- dat$lat
  data$lon <- dat$lon+100  
  data$nsp <- dat$nsp
  data$abund <- dat$density*dat$asurv
  data$asampl <- dat$aswept
  data$temp <- (dat$sst+273.15)*8.62/100000
  data$cat <- as.integer(dat$mlgr)-1
  data$npp <- dat$npp  
  data$mesh <- dat$mesh
  data$siz <- dat$siz
  data$depth <- dat$depth
  data$density <- dat$density
  param <- list()
  if(data$code==(-2)){
    param$loga7 <- 14
    param$logk <- 3
  }
  if(data$code==(-1)){   
    param$logb5 <- rep(-1,length(data$nsp))
    param$logk <- 3
  }
  if(data$code==1){
    param$b5 <- rep(-.5,length(unique(data$cat)))
    param$b0 <- 20
    param$b1 <- -0.1
    param$b2 <- 0.4
    param$logb3 <- 1
    param$b4 <- 0.2
    param$logk <- 5
    param$dummy <- 0
  }
  if(data$code==2){
    param$loga7 <- 20
    param$b8 <- rep(0,length(unique(data$cat)))
    param$logb1 <- -1
    param$logb2 <- -5
    param$logb3 <- -6
    param$logb5 <- -6
    param$logb6 <- -6
    param$logb7 <- 1
    param$logk <-1
  }
  if(data$code==3){
    param$loga7 <- -2
    param$b5 <- rep(0,length(unique(data$cat)))
    param$b0 <- 0
    param$b1 <- 0
    param$b2 <- 0
    param$b3 <- 0.1
    param$b4 <- 0.1
    param$logk <- 0
  }
  if(data$code==4){
    param$loga7 <- 16
    param$b5 <- rep(0,length(unique(data$cat)))
    param$logb1 <- -1
    param$logb2 <- -1
    param$logb3 <- -1
    param$logb4 <- 0.1
    param$logk <- 1
  }

  mymap=list()
  if(!missing(fixK)){
    mymap$logk=factor(NA)
    param$logk=log(fixK)
  }
  
  obj <- MakeADFun(data, param, DLL="biodiversity", map=mymap, silent=TRUE)
  
  low <- rep(-Inf,length(obj$par))
  hig <- rep(Inf,length(obj$par))
  if(data$code==(-2)){
    low <- c(-100, -Inf)
    hig <- c(100, Inf)
  }  
  if(data$code==(-1)){
  }
  if(data$code==1){
    low <- c(-10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -Inf, -Inf, -Inf, -100, 0, -100)
    hig <- c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, Inf, Inf, Inf, 100, Inf, 100)
  }
  if(data$code==2){
    low <- c(-Inf, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -Inf, -Inf, -Inf, -Inf, -Inf, -10, -100)
    hig <- c(Inf, 10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10, Inf, Inf,  Inf,  Inf,  Inf, 10,  100)
  }
  if(data$code==3){
    low <- c(-Inf,-10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,-Inf, -Inf, -Inf, -1, -100)
    hig <- c(Inf, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, Inf, Inf, Inf, 1, 100)
  }
  if(data$code==4){
    low <- c(0, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,-10, -Inf, -10, -100)
    hig <- c(Inf, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 2, 10, 10, 10, 100)
  }
  opt <- nlminb(obj$par, obj$fn, obj$gr, lower=low, upper=hig, control=list(eval.max=10000, iter.max=10000))
  if(opt$convergence!=0){ # try again 
    opt <- nlminb(opt$par, obj$fn, obj$gr, lower=low, upper=hig, control=list(eval.max=10000, iter.max=10000))
  }
  rep <- obj$report()
  sdrep <- sdreport(obj,opt$par)
  ret <- list(opt=opt, obj=obj, data=data, rep=rep, sdrep=sdrep)
  attr(ret,"od")<-origdata
  class(ret)<-"biodiv"
  return(ret)
}

##' R2 of bioviv output
##' @param fit as returned from the biodiv function 
##' @return R2
##' @details ...
##' @export
R2 <- function(fit, ...){
  1-sum((fit$data$nsp-fit$rep$mu)^2)/sum((fit$data$nsp-mean(fit$data$nsp))^2)
}


##' Deviance fraction of bioviv output
##' @param fit as returned from the biodiv function.
##' @return 
##' @details ...
##' @export
devi <- function(fit, ...){
  fit0<-biodiv(attr(fit,"od"), conf=-2, fixK=exp(fit$opt$par["logk"]))
  fits<-biodiv(attr(fit,"od"), conf=-1, fixK=exp(fit$opt$par["logk"]))
  D<-(1-(logLik(fit)-logLik(fits))/(logLik(fit0)-logLik(fits)))
  attributes(D)<-NULL
  D
}



