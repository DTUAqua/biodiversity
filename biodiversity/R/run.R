##' Fit biodiversity model 
##' @param data A data set in same format as the species data set supplied with this package 
##' @param conf A code identifying the model configuration (now only one is possible) 
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
biodiv <- function(data, conf, run=TRUE, ...){
  dat<-data[(data$density*data$asurv)>1.0e-6,]
  data <- list()
  data$code <- conf
  if(data$code==1){
    data$lat <- dat$lat
    data$nsp <- dat$nsp
    data$abund <- dat$density*dat$asurv
    data$logsiz <- dat$siz
    data$asampl <- dat$aswept
    data$mesh <- dat$mesh
    data$temp <- (dat$sst+273.15)*8.62/100000
    data$npp <- dat$npp
    data$cat <- as.integer(dat$mlgr)-1
  }
  if(data$code==2){
    data$nsp <- dat$nsp
  }
  if(data$code==3){
    data$lat <- dat$lat
    data$lon <- dat$lon+100
    data$nsp <- dat$nsp
    data$abund <- dat$density*dat$asurv
    data$siz <- dat$siz
    data$asampl <- dat$aswept
    data$mesh <- dat$mesh
    data$depth <- dat$depth
    data$npp <- dat$npp
    data$cat <- as.integer(dat$mlgr)-1
  }
  if(data$code==4){
    data$nsp <- dat$nsp
    data$density <- dat$density
    data$siz <- dat$siz
    data$asampl <- dat$aswept
    data$mesh <- dat$mesh
    data$temp <- (dat$sst+273.15)*8.62/100000
    data$npp <- dat$npp
    data$cat <- as.integer(dat$mlgr)-1
  }

  param <- list()
  if(data$code==1){
    param$b5<- rep(-.5,length(unique(data$cat)))
    param$b0 <- 20
    param$b1 <- -0.1
    param$b2 <- 0.4
    param$logb3 <- 1
    param$b4 <- 0.2
    param$logk<-3
    param$dummy<-0
  }
  if(data$code==2){
    param$logb5<- rep(-1,length(data$nsp))
  }
  if(data$code==3){
    param$loga7 <- -2
    param$b5<- rep(0,length(unique(data$cat)))
    param$b0 <- 0
    param$b1 <- 0
    param$b2 <- 0
    param$b4 <- 0.1
    param$logk<-0
  }
  if(data$code==4){
    param$loga7 <- 16
    param$b5<- rep(0,length(unique(data$cat)))
    param$logb1 <- -1
    param$logb2 <- -1
    param$logb3 <- -1
    param$logb4 <- 0.1
    param$logk<-1
  }
  
  obj <- MakeADFun(data, param, DLL="biodiversity", silent=TRUE)
  
  low <- rep(-Inf,length(obj$par))
  hig <- rep(Inf,length(obj$par))
  if(data$code==1){
    low <- c(-10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -Inf, -Inf, -Inf, -100, 0, -100)
    hig <- c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, Inf, Inf, Inf, 100, Inf, 100)
  }
  if(data$code==2){
  }
  if(data$code==3){
    low <- c(-Inf,-10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,-Inf, -Inf, -1, -100)
    hig <- c(Inf, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, Inf, Inf, 1, 100)
  }
  if(data$code==4){
    low <- c(0, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,-10, -Inf, -10, -100)
    hig <- c(Inf, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 2, 10, 10, 10, 100)
  }
  opt <- nlminb(obj$par, obj$fn, obj$gr, lower=low, upper=hig, control=list(eval.max=10000, iter.max=10000))
  rep <- obj$report()
  sdrep <- sdreport(obj,opt$par)
  ret <- list(opt=opt, obj=obj, data=data, rep=rep, sdrep=sdrep)
  class(ret)<-"biodiv"
  return(ret)
}

