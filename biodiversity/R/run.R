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
biodiv <- function(data, conf, fixK=NULL, useTotCatch=TRUE, useAswept=TRUE, run=TRUE, ...){
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
  #data$temp <- (dat$sst+273.15)*8.62/100000
  data$temp <- (dat$ult+273.15)*8.62/100000
  data$cat <- as.integer(dat$mlgr)-1
  data$npp <- dat$npp  
  data$mesh <- dat$mesh
  data$siz <- dat$lml
  data$depth <- dat$depth
  data$density <- dat$density
  data$totCatch <- dat$catch
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
    param$b3 <- 0.5
    param$b4 <- 0.2
    param$b8 <- 0.1    
    param$logk <- 3
    #param$dummy <- 0
  }
  if(data$code==2){
    param$loga7 <- 20
    param$b8 <- rep(0,length(unique(data$cat)))
    param$b1 <- exp(-1)
    param$b2 <- exp(-5)
    param$b3 <- -exp(-6)
    param$b4 <- 0.1
    param$b5 <- exp(-6)
    param$b6 <- exp(-6)
    param$b7 <- -exp(1)
    param$logk <- 4
  }
  if(data$code==3){
    param$loga7 <- -2
    param$b5 <- rep(0,length(unique(data$cat)))
    param$b0 <- 0
    param$b1 <- 0
    param$b2 <- 0
    param$b3 <- -0.1
    param$b4 <- 0.1
    param$b8 <- 0.1    
    param$logk <- 0
  }
  if(data$code==4){
    param$loga7 <- 16
    param$b5 <- rep(0,length(unique(data$cat)))
    param$b1 <- exp(-1)
    param$b2 <- exp(-1)
    param$b3 <- exp(-1)
    param$b4 <- exp(0.1)
    param$b8 <- 0.1        
    param$logk <- 1
  }

  mymap=list()
  if(!missing(fixK)){
    mymap$logk=factor(NA)
    param$logk=log(fixK)
  }
  if(!useTotCatch){
    if(data$code==2){
        mymap$b4=factor(NA)
        param$b4 <- 0
    }
    if(data$code%in%c(1,3,4)){
        mymap$b8=factor(NA)
        param$b8 <- 0
    }
    if(data$code%in%c(-2,-1)){
      warning("option not relevant for this configuration")
    }
  }
  if(!useAswept){
    if(data$code==2){
        mymap$b5=factor(NA)
        param$b5 <- 0
    }
    if(data$code%in%c(1,3,4)){
        mymap$b4=factor(NA)
        param$b4 <- 0
    }
    if(data$code%in%c(-2,-1)){
      warning("option not relevant for this configuration")
    }
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
    low <- c(-10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -Inf, -Inf, 0, 1.0e-6, if(useAswept){0}else{NULL}, if(useTotCatch){1}else{NULL}, -100)
    hig <- c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, Inf, Inf, Inf, 1-1.0e-6, if(useAswept){Inf}else{NULL}, if(useTotCatch){1}else{NULL}, 100)
  }
  if(data$code==2){
    low <- c(-Inf, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, 0, 0, -Inf, if(useTotCatch){0}else{NULL}, if(useAswept){0}else{NULL}, 0, -Inf, -100)
    hig <- c(Inf, 10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10, Inf, Inf,  0, if(useTotCatch){1}else{NULL}, if(useAswept){Inf}else{NULL},  Inf, 0,  100)
  }
  if(data$code==3){
    low <- c(-Inf,-10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,-Inf, -Inf, -Inf, if(useAswept){0}else{NULL}, if(useTotCatch){0}else{NULL}, -100)
    hig <- c(Inf, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, Inf, Inf, Inf, if(useAswept){Inf}else{NULL}, if(useTotCatch){1}else{NULL}, 100)
  }
  if(data$code==4){
    low <- c(0, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -Inf,-Inf, -Inf, if(useAswept){0}else{NULL}, if(useTotCatch){0}else{NULL}, -100)
    hig <- c(Inf, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, Inf, Inf, Inf, if(useAswept){Inf}else{NULL}, if(useTotCatch){1}else{NULL}, 100)
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
##' @details Calculates R-squared coefficient of determination as: 1-SS(model residuals)/SS(total residuals)
##' @export
R2 <- function(fit, ...){
  1-sum((fit$data$nsp-fit$rep$mu)^2)/sum((fit$data$nsp-mean(fit$data$nsp))^2)
}

##' Deviance fraction of bioviv output
##' @param fit as returned from the biodiv function.
##' @return deviance fraction (see details)
##' @details Deviance fraction calculated by comparing the deviance between the estimated model (M) and the saturated model (S) to the deviance between the constant model(C) and the saturated model (S). The constant and saturated models are optimized with same overdispersion parameter as the estimated model. The deviance fraction is defined as: 1-(logLik(M)-logLik(S))/(logLik(C)-logLik(S))
##' @export
devi <- function(fit, ...){
  fit0<-biodiv(attr(fit,"od"), conf=-2, fixK=exp(fit$opt$par["logk"]))
  fits<-biodiv(attr(fit,"od"), conf=-1, fixK=exp(fit$opt$par["logk"]))
  D<-(1-(logLik(fit)-logLik(fits))/(logLik(fit0)-logLik(fits)))
  attributes(D)<-NULL
  D
}

##' Re-run a model without a given covariate  
##' @param fit as returned from the biodiv function.
##' @param what quoted name of covariate to exclude (done by assigning it to its average value) 
##' @return The fraction that the average variance is reduced by including the named covariate 
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats nlminb
##' @useDynLib biodiversity
runwithout<-function(fit, what){
  data2<-fit$data
  param2<- as.list(fit$sdrep, "Est")
  data2[[what]]<-rep(mean(data2[[what]]), length(data2[[what]]))
  if(what=="lml"){
    data2[["cat"]][] <- 0
  }
  mymap <- fit$obj$env$map
  paramlist <- unlist(param2[!names(param2)%in%names(mymap)])
  npar <- length(paramlist)
  mylow <- rep(-Inf,npar)
  if(fit$data$code==2){
    mylow[which(names(paramlist)=="b4")] <- 0
  }else{
    mylow[which(names(paramlist)=="b8")] <- 0
  }
  myhig <- rep(Inf,npar)
  obj2 <- MakeADFun(data2, param2, DLL = "biodiversity", silent = TRUE, map=mymap)
  opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, control=list(eval.max=10000, iter.max=10000), lower=mylow, upper=myhig)
  frac = (mean(obj2$report()$var) - mean(fit$obj$report()$var)) / mean(obj2$report()$var)
  frac
}

##' Variance reduction by each term 
##' @param fit as returned from the biodiv function.
##' @param lookat a vector of quoted covariate names to be looked at 
##' @return A vector of variance reduction fractions 
##' @export 
##' @examples 
##' order<-c(2,4,1,3)
##' fits <- lapply(order, function(i)biodiv(species,i))
##' cols <- lapply(fits, function(fit){vr<-varianceReduction(fit); d<-devi(fit); c(vr/sum(vr,na.rm=TRUE)*d,unexplained=1-d)})
##' bars <- do.call(cbind, cols)
##' colnames(bars) <- c("Neutral", "Best", "Latitude", "Metabolic")[order]
##' par(mar=c(3,4,1,9)+0.5,oma=c(0,0,0,3), ask=FALSE)
##' colvec <- rgb(rbind(c(105,0,0),c(204,0,0),c(204,98,0),c(204,196,0),
##'                     c(0,204,130),c(140,204,204),c(0,131,204),c(65,0,204),c(163,0,204),
##'                     c(120,120,120)),alpha=204,max=204)
##' op <- palette(colvec)
##' barplot(bars,col=c(op,"white"),names.arg=colnames(bars),
##'       legend.text=rownames(bars),ylab="Proportion of deviance",
##'         args.legend=list(x=7.5,y=1.05,bty="n"))
##' 
varianceReduction<-function(fit, lookat=c("lml", "temp", "asampl", "mesh", "density", "npp", "depth", "abund", "lat", "lon")){
  vr <- sapply(lookat, function(x)runwithout(fit,x))
  vr[vr<1.0e-6] <- 0
  vr
}
