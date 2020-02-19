library(biodiversity)

m1<-biodiv(species, conf=1, useTotCatch=FALSE)
m2<-biodiv(species, conf=2, useTotCatch=FALSE)
m3<-biodiv(species, conf=3, useAswept=FALSE)
m4<-biodiv(species, conf=4, useTotCatch=FALSE)

fits<-list(m1=m1, m2=m2, m3=m3, m4=m4)

f<-function(fit){
    options(digits=3)
    print(coef(fit))
    options(digits=7)
    print(paste("Devi:", signif(devi(fit),3)))
    print(paste("R2:", signif(R2(fit),3)))
    print(paste("AIC:", signif(AIC(fit),4)))
}

sink(file="res.out")
  dummy<-lapply(fits, f)
sink()
