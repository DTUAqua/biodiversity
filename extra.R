rm(list=ls())
par(ask=FALSE)

library("gtable")
library("foreign")
library("MASS")
library("lmtest")
library("mgcv")
library("car")
library("lattice")
library("modEvA")
library("ggplot2")
library("plyr")
library("car")
library("marmap")
library("grid")

# Read basic input data

trends0<-read.table("https://raw.githubusercontent.com/DTUAqua/biodiversity/master/biodiversity/data/species.tab?token=AFH0d1WQ-rvtPg9Qh-oiRe3fGxVZrPS0ks5bF8L6wA", header=TRUE,row.names = NULL) 
names(trends0)
trends0$lnaswept<-log(trends0$aswept)
trends0$lnnhauls<-log(trends0$nhauls)
trends0$lnasurv<-log(trends0$asurv)
trends0$lnmesh<-log(trends0$mesh)
trends0$emesh<-exp(trends0$mesh)
trends0$lnhorop<-log(trends0$horop)
trends0$lnvertop<-log(trends0$vertop)
trends0$lntowsp<-log(trends0$towsp)
trends0$lnlat<-log(trends0$lat)


#Sea surface temperatures
trends0$sstemp<-trends0$sst
trends0$sst<-(trends0$sst+273.15)
# temperatures can be converted to eV by multiplication with 8.62/100000
trends0$sstdif<-(trends0$sstdif)
trends0$invsst<-1/(trends0$sst)
trends0$einvsst<-exp(trends0$invsst)
trends0$invsstdif<-1/(-0.5*trends0$sstdif+trends0$sst)-1/(0.5*trends0$sstdif+trends0$sst)
trends0$invsstdif2<-trends0$invsstdif*trends0$invsstdif


#bottom temperatures
trends0$sbt<-(trends0$sbt+273.15)
trends0$sbtdif<-(trends0$sbtdif)
trends0$invsbt<-1/(trends0$sbt)
trends0$einvsbt<-exp(trends0$invsbt)
trends0$invsbtdif<-1/(-0.5*trends0$sbtdif+trends0$sbt)-1/(0.5*trends0$sbtdif+trends0$sbt)
trends0$invsbtdif2<-trends0$invsbtdif*trends0$invsbtdif

#Upper 200 m average temperaturew
trends0$ult<-(trends0$ult+273.15)
trends0$ultdif<-(trends0$ultdif)
trends0$invult<-1/(trends0$ult)
trends0$einvult<-exp(trends0$invult)
trends0$invultdif<-1/(-0.5*trends0$ultdif+trends0$ult)-1/(0.5*trends0$ultdif+trends0$ult)
trends0$invultdif2<-trends0$invultdif*trends0$invultdif

trends0$lnnpp<-log(trends0$npp)
trends0$lndepth<-log(trends0$depth)
trends0$lnsiz<-log(trends0$siz)
trends0$siz2<-trends0$siz*trends0$siz
trends0$siz3<-trends0$siz2*trends0$siz

#conversion of zero richness and density to 1/nhauls to avoid having to take the log of zero
trends0$lnnsp<-ifelse(trends0$nsp<1e-10,log(trends0$nsp+1/(trends0$nhauls)),log(trends0$nsp))
trends0$lndensity<-ifelse(trends0$density<1e-10,log(trends0$density+1/(trends0$nhauls)),log(trends0$density))
trends0$lnabundance<-(trends0$lnasurv+trends0$lndensity)


names(trends0)
summary(trends0)
par(mfrow=c(1,1))


# Colors
colvec <- rgb(rbind(c(160,0,0),c(204,0,0),c(204,98,0),c(204,196,0),c(134,204,0),c(16,204,0),
                 c(0,204,130),c(0,131,204),c(65,0,204),c(163,0,204),c(120,120,120)),
              alpha=204,max=204)
palette(colvec)

# Figure 1 Map of survey areas with pies showing relative number of species in each survey

getNOAA.bathy(lon1=-80,lon2=40,lat1=5,lat2=85,resolution=4,keep=TRUE)->nea
plot.new()
old.par <- par(no.readonly=TRUE)
grey<-colorRampPalette(c("grey100","grey100"))
blues<-colorRampPalette(c("steelblue4","powderblue"))
pdf("Fig1_Survey_map.pdf",width=7,height=5.712 )
plot.bathy(nea,image=TRUE,land=TRUE,n=0,bpal=list(c(0,max(nea),grey(10)),c(min(nea),0,blues(100))),asp=1)
neadat<-read.table("https://raw.githubusercontent.com/DTUAqua/biodiversity/master/biodiversity/data/pies.tab?token=AFH0d4c4nknpN_SvR9UvhJyznxPapzPKks5bF9A5wA",header=TRUE)
points(neadat[,2:3],pch=19,col="blue", cex=0.5)
space.pies(neadat$lon, neadat$lat,
           pie.slices=neadat[,10:14], pie.colors=neadat[,16:20], pie.radius=2, pie.space=0.01)
col1 <- c("#2DFF09","#E8DD08","#FFA804","#E84108","#A818FF")
namvec <- c(" <16cm ","16-43cm  ","43-116cm ","116-314cm "," >314cm")
par(plt=c(0.72,0.82,0.22,0.32),new=TRUE)
pie(rep(1,5),labels=as.character(namvec),col=col1,cex=0.7,radius=1)
par(new=F)
par(old.par)
dev.off()


#Dataset with species lengths collapsed

tr1<-read.table("https://raw.githubusercontent.com/DTUAqua/biodiversity/master/biodiversity/data/areas.tab?token=AFH0d-Aa_LraFP5qcpHjBHO2Z1K54kRBks5bF8jWwA",header=TRUE,row.names=NULL)
names(tr1)
summary(tr1)

#Correlation between independent variables

tr1$lnnsp<-log(tr1$nsp)
tr1$lndensity<-log(tr1$density)
tr1$sst<-(273.15+tr1$sst)
tr1$ult<-(273.15+tr1$ult)
tr1$sbt<-(273.15+tr1$sbt)

#latitude
cor.test(tr1$lat,tr1$sst)
cor.test(tr1$lat,tr1$sstdif)
cor.test(tr1$lat,tr1$ult)
cor.test(tr1$lat,tr1$ultdif)
cor.test(tr1$lat,tr1$sbt)
cor.test(tr1$lat,tr1$sbtdif)
cor.test(tr1$lat,tr1$npp)
cor.test(tr1$lat,tr1$lndensity)
cor.test(tr1$lat,tr1$depth)
cor.test(tr1$lat,tr1$asurv)
cor.test(tr1$lat,tr1$aswept)
cor.test(tr1$lat,tr1$mesh)
cor.test(tr1$lat,tr1$vertop)
#longitude
cor.test(tr1$lon,tr1$sst)
cor.test(tr1$lon,tr1$sstdif)
cor.test(tr1$lon,tr1$ult)
cor.test(tr1$lon,tr1$ultdif)
cor.test(tr1$lon,tr1$sbt)
cor.test(tr1$lon,tr1$sbtdif)
cor.test(tr1$lon,tr1$npp)
cor.test(tr1$lon,tr1$lndensity)
cor.test(tr1$lon,tr1$depth)
cor.test(tr1$lon,tr1$asurv)
cor.test(tr1$lon,tr1$aswept)
cor.test(tr1$lon,tr1$mesh)
cor.test(tr1$lon,tr1$vertop)
#sst
cor.test(tr1$sst,tr1$sstdif)
cor.test(tr1$sst,tr1$ult)
cor.test(tr1$sst,tr1$ultdif)
cor.test(tr1$sst,tr1$sbt)
cor.test(tr1$sst,tr1$sbtdif)
cor.test(tr1$sst,tr1$npp)
cor.test(tr1$sst,tr1$density)
cor.test(tr1$sst,tr1$depth)
cor.test(tr1$sst,tr1$asurv)
cor.test(tr1$sst,tr1$aswept)
cor.test(tr1$sst,tr1$mesh)
cor.test(tr1$sst,tr1$vertop)
#sstdif
cor.test(tr1$sstdif,tr1$ult)
cor.test(tr1$sstdif,tr1$ultdif)
cor.test(tr1$sstdif,tr1$sbt)
cor.test(tr1$sstdif,tr1$sbtdif)
cor.test(tr1$sstdif,tr1$npp)
cor.test(tr1$sstdif,tr1$density)
cor.test(tr1$sstdif,tr1$depth)
cor.test(tr1$sstdif,tr1$asurv)
cor.test(tr1$sstdif,tr1$aswept)
cor.test(tr1$sstdif,tr1$mesh)
cor.test(tr1$sstdif,tr1$vertop)
#ult
cor.test(tr1$ult,tr1$ultdif)
cor.test(tr1$ult,tr1$sbt)
cor.test(tr1$ult,tr1$sbtdif)
cor.test(tr1$ult,tr1$npp)
cor.test(tr1$ult,tr1$density)
cor.test(tr1$ult,tr1$depth)
cor.test(tr1$ult,tr1$asurv)
cor.test(tr1$ult,tr1$aswept)
cor.test(tr1$ult,tr1$mesh)
cor.test(tr1$ult,tr1$vertop)
#ultdif
cor.test(tr1$ultdif,tr1$sbt)
cor.test(tr1$ultdif,tr1$sbtdif)
cor.test(tr1$ultdif,tr1$npp)
cor.test(tr1$ultdif,tr1$density)
cor.test(tr1$ultdif,tr1$depth)
cor.test(tr1$ultdif,tr1$asurv)
cor.test(tr1$ultdif,tr1$aswept)
cor.test(tr1$ultdif,tr1$mesh)
cor.test(tr1$ultdif,tr1$vertop)
#sbt
cor.test(tr1$sbt,tr1$sbtdif)
cor.test(tr1$sbt,tr1$npp)
cor.test(tr1$sbt,tr1$density)
cor.test(tr1$sbt,tr1$depth)
cor.test(tr1$sbt,tr1$asurv)
cor.test(tr1$sbt,tr1$aswept)
cor.test(tr1$sbt,tr1$mesh)
cor.test(tr1$sbt,tr1$vertop)
#sbtdif
cor.test(tr1$sbtdif,tr1$npp)
cor.test(tr1$sbtdif,tr1$density)
cor.test(tr1$sbtdif,tr1$depth)
cor.test(tr1$sbtdif,tr1$asurv)
cor.test(tr1$sbtdif,tr1$aswept)
cor.test(tr1$sbtdif,tr1$mesh)
#npp
cor.test(tr1$npp,tr1$lndensity)
cor.test(tr1$npp,tr1$depth)
cor.test(tr1$npp,tr1$asurv)
cor.test(tr1$npp,tr1$aswept)
cor.test(tr1$npp,tr1$mesh)
cor.test(tr1$npp,tr1$vertop)
#density
cor.test(tr1$density,tr1$depth)
cor.test(tr1$density,tr1$asurv)
cor.test(tr1$density,tr1$aswept)
cor.test(tr1$density,tr1$mesh)
cor.test(tr1$density,tr1$vertop)
#depth
cor.test(tr1$depth,tr1$asurv)
cor.test(tr1$depth,tr1$aswept)
cor.test(tr1$depth,tr1$density)
cor.test(tr1$depth,tr1$mesh)
cor.test(tr1$depth,tr1$vertop)
#asurv
cor.test(tr1$asurv,tr1$aswept)
cor.test(tr1$asurv,tr1$density)
cor.test(tr1$asurv,tr1$mesh)
cor.test(tr1$asurv,tr1$vertop)
#aswept
cor.test(tr1$aswept,tr1$density)
cor.test(tr1$aswept,tr1$mesh)
cor.test(tr1$aswept,tr1$vertop)
#mesh
cor.test(tr1$mesh,tr1$aswept)
cor.test(tr1$mesh,tr1$vertop)
#vertop
cor.test(tr1$vertop,tr1$aswept)

# Figure S1 Log number of species versus independent variables 

MyV<-c("lat","lon","sst","sstdif","ult","ultdif","sbt","sbtdif","npp","depth","asurv","aswept","vertop","mesh","lndensity")
Myxyplot <- function(Z, MyV, NameY1) {
  AllX  <- as.vector(as.matrix(Z[,MyV]))
  AllY  <- rep(Z[,NameY1] , length(MyV))
  AllID <- rep(MyV, each = nrow(Z))
  
  
  
  P <- xyplot(AllY ~ AllX|factor(AllID), col = 1,
              xlab = "Independent variables",
              ylab = "Ln No of Species",
              strip = function(bg='white', ...)
                strip.default(bg='white', ...),
              scales = list(alternating = T,
                            x = list(relation = "free"),
                            y = list(relation = "same")),
              panel=function(x, y){
                panel.grid(h=-1, v= 2)
                panel.points(x, y, col = "black", pch=1, cex=0.5)
                panel.loess(x, y, span = 0.7,col = "black", lwd = 2)})
  
  print(P)
}
Myxyplot(tr1,MyV,"lnnsp")
pdf("FigS1_Independent_variables.pdf",width=7,height=5 )
Myxyplot(tr1,MyV,"lnnsp")
dev.off()


#Figure S2 Correlations

pair<-read.table("https://raw.githubusercontent.com/DTUAqua/biodiversity/master/biodiversity/data/pairs.tab?token=AFH0d2ZpwdIB32WbrDPH_aa--ckqU8URks5bF8mlwA",header=TRUE,row.names=NULL)
names(pair)
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  cex.axis=0.5
  if(missing(cex.cor)) cex <- 0.70/strwidth(txt)
  #  text(0.5, 0.5, txt, cex = cex * sqrt(abs(r)))
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", " ")) 
  text(0.5, 0.45, txt, cex = 1.0*cex*0.8 
       #        * (abs(r))^(1/2)
  ) 
  text(.7, .8, Signif, cex=cex, col="red")
}
panel.smooth<-function (x, y, 
                        cex = 0.7, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}
pairs(pair,cex.labels=1, upper.panel=panel.smooth, lower.panel=panel.cor, gap=0.2)
pdf("FigS2_Correlations.pdf",width=9,height=7 )
pairs(pair,cex.labels=1,upper.panel=panel.smooth, lower.panel=panel.cor, gap=0.2)
dev.off()




# Figure 2 Plot of density and richness

#Summary function
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {

  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars,.drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

par(mfrow=c(1,1))
tr<-trends0
tr$temperature<-tr$band
tr
any(is.na(tr))
tgc<-summarySE(tr,measurevar="nsp",groupvars=c("temperature","mlgr"),na.rm=FALSE)
any(is.na(tgc))
which(is.na(tgc) == TRUE, arr.ind=TRUE)
tgc
pd<-position_dodge(0.1)
tgc$temperature<-factor(tgc$temperature,levels=c(">20.5","14-20.5","7.5-14","<7.5"))
p1<-ggplot(tgc, aes(x=mlgr, y=nsp,group=temperature,colour=temperature)) +
  geom_errorbar(aes(ymin=nsp-ci, ymax=nsp+ci), width=.5, position=pd) +
  geom_line(position=pd,size=1.5) +
  geom_point(position=pd,size=4)+
  scale_colour_manual(values=c(colvec[2],colvec[10],colvec[8],colvec[11]))+
  theme_bw()+
  theme(legend.title = element_text(size=16))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.justification=c(1,1), legend.position=c(0.95,0.95))+
  ylab("No. of species")+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())+
  theme(legend.key=element_rect(fill='NA'))+
  scale_colour_manual(name="Temperature",values=c(colvec[2],colvec[10],colvec[8],colvec[11]))+
  theme(axis.title.y = element_text(size=20),
      axis.text.y  = element_text(size=16))
tgc<-summarySE(tr,measurevar="lndensity",groupvars=c("temperature","mlgr"))
pd<-position_dodge(0.1)
tgc$temperature<-factor(tgc$temperature,levels=c(">20.5","14-20.5","7.5-14","<7.5"))
p2<-ggplot(tgc, aes(x=mlgr, y=lndensity,group=temperature,colour=temperature)) +
  geom_errorbar(aes(ymin=lndensity-ci, ymax=lndensity+ci), width=.5, position=pd) +
  geom_line(position=pd,size=1.5) +
  geom_point(position=pd,size=4)+
  scale_x_discrete(labels=c("3-6","6-9","9-16","16-26","26-43","43-70","70-116","116-191","191-314","314-518","518-854"))+
  theme_bw()+
  theme(legend.position="none")+
  scale_colour_manual(name="Temperature",values=c(colvec[2],colvec[10],colvec[8],colvec[11]))+
  ylab("Log density")+xlab("Maximum length")+
  theme(axis.title.x = element_text(size=20),
        axis.text.x  = element_text(angle=45, vjust=0.5, size=16))+
  theme(axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=16))


g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
maxwidths <- grid::unit.pmax(g1$widths[2:5], g2$widths[2:5])
g1$widths[2:5] <- as.list(maxwidths)
g2$widths[2:5] <- as.list(maxwidths)
g <- gtable_matrix(name = "demo",
                   grobs = matrix(list(g1, g2), nrow = 2), 
                   widths = unit(13, "cm"),
                   heights = unit(c(8, 7), "cm"))
grid.newpage()

grid.draw(g) 

pdf("Fig2_Density_and_richness_N.pdf",width=6,height=6)
grid.draw(g) 
dev.off()

#Plotting concurvity
vis.concurvity <- function(b, type="estimate"){
  # arguments:
  #   b     -- a fitted gam
  #   type  -- concurvity measure to plot, see ?concurvity
  cc <- concurvity(b, full=FALSE)[[type]]
  
  diag(cc) <- NA
  cc[lower.tri(cc)]<-NA
  
  layout(matrix(1:2, ncol=2), widths=c(5,1))
  opar <- par(mar=c(5, 6, 5, 0) + 0.1)
  # main plot
  image(z=cc, x=1:ncol(cc), y=1:nrow(cc), ylab="", xlab="",
        axes=FALSE, asp=1, zlim=c(0,1))
  axis(1, at=1:ncol(cc), labels = colnames(cc), las=2)
  axis(2, at=1:nrow(cc), labels = rownames(cc), las=2)
  # legend
  opar <- par(mar=c(5, 0, 4, 3) + 0.1)
  image(t(matrix(rep(seq(0, 1, len=100), 2), ncol=2)),
        x=1:3, y=1:101, zlim=c(0,1), axes=FALSE, xlab="", ylab="")
  axis(4, at=seq(1,101,len=5), labels = round(seq(0,1,len=5),1), las=2)
  par(opar)
}




#Analysis of SPECIES RICHNESS

#What does a GAM tell us?
#Selecting the best fit by backwards elimination of insignificant terms and terms with an estimated concurvity>0.85
  
#Sea surface temperature  
  m01a<-gam(nsp~s(sst,k=4)+s(sstdif,k=4)+s(asurv,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(vertop,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01a)  
  concurvity(m01a)
  vis.concurvity(m01a)
  m01a<-gam(nsp~s(sst,k=4)+s(asurv,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(vertop,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01a)  
  concurvity(m01a)
  vis.concurvity(m01a)
  m01a<-gam(nsp~s(sst,k=4)+s(asurv,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01a)  
  concurvity(m01a)
  vis.concurvity(m01a)
  m01a<-gam(nsp~s(sst,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01a)  
  concurvity(m01a)
  vis.concurvity(m01a)
  par(mfrow=c(3,3))
  plot(m01a)
  anova.gam(m01a)
  
 #Bottom temperature
  m01b<-gam(nsp~s(sbt,k=4)+s(sbtdif,k=4)+s(asurv,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(vertop,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01b)
  concurvity(m01b)
  vis.concurvity(m01b)
  m01b<-gam(nsp~s(sbt,k=4)+s(sbtdif,k=4)+s(asurv,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01b)
  concurvity(m01b)
  vis.concurvity(m01b)
  m01b<-gam(nsp~s(sbt,k=4)+s(asurv,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01b)
  concurvity(m01b)
  vis.concurvity(m01b)
  m01b<-gam(nsp~s(sbt,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01b)
  concurvity(m01b)
  vis.concurvity(m01b) 
  
 #Temperature in the upper 200 m of the watercolumn 
  m01c<-gam(nsp~s(ult,k=4)+s(ultdif,k=4)+s(asurv,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(vertop,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01c)
  concurvity(m01c)
  vis.concurvity(m01c)
  m01c<-gam(nsp~s(ult,k=4)+s(ultdif,k=4)+s(asurv,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01c)
  concurvity(m01c)
  vis.concurvity(m01c)
  m01c<-gam(nsp~s(ult,k=4)+s(ultdif,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01c)
  concurvity(m01c)
  vis.concurvity(m01c)
  m01c<-gam(nsp~s(ult,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01c)
  concurvity(m01c)
  vis.concurvity(m01c)
  
  #Latitude and longitude
  m01d<-gam(nsp~s(lat,k=4)+s(lon,k=4)+s(asurv,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(vertop,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01d)
  concurvity(m01d)
  vis.concurvity(m01d)
  m01d<-gam(nsp~s(lat,k=4)+s(lon,k=4)+s(asurv,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01d)
  concurvity(m01d)
  vis.concurvity(m01d)
  m01d<-gam(nsp~s(lat,k=4)+s(lon,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
            +s(npp,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
  summary(m01d)
  concurvity(m01d)
  vis.concurvity(m01d)
  
AIC(m01a,m01b,m01c,m01d)
#Model a (sst ) is the best 
sum(resid(m01a,type="pearson")^2)/m01a$df.res
#and the variation of the residuals are almost exactly as expected in a negative binomial.
#The model with Sea surface temperature is the best

#Checking whether a Poisson model would provide a better fit
m01ap<-gam(nsp~s(sst,k=4)+s(aswept,k=4)+s(density,k=4)+s(depth,k=4)+s(siz,k=4)
           +s(npp,k=4)+s(mesh,k=4,by=mlgr),family=poisson,method="REML",select="TRUE",data=trends0)
summary(m01ap)
concurvity(m01ap)
vis.concurvity(m01ap)
par(mfrow=c(4,3))
plot(m01ap)
sum(resid(m01ap,type="pearson")^2)/m01ap$df.res
sum(resid(m01a,type="pearson")^2)/m01a$df.res
AIC(m01a,m01ap)
#The negative binomial model has the lowest AIC and describes the variance better than the Poisson model

#Figure 3 Estimated smoothing curves obtained by a GAM model using sea surface temperature
#par(mfrow=c(3,3))
#plot(m01a,residuals=FALSE,rug=TRUE,se=TRUE,shade=TRUE,shade.col="gray80",seWithMean=TRUE)
pdf("Fig3_GAM_results.pdf",width=7,height=8)
par(mfrow=c(3,3))
plot(m01a,residuals=FALSE,rug=TRUE,se=TRUE,shade=TRUE,shade.col="gray80",seWithMean=TRUE)
dev.off()

#Figure S4 Survey residualsfrom GAM

par(mfrow=c(1,1))
old.par <- par(no.readonly=T)
par(las=3)
par(mar=c(8,8,1,1),oma=c(3,0,0,0))
plot(x=trends0$survey,y=resid(m01a),ylab="Residuals",cex.axis=0.7,pch=19,col=1)
pdf("FigS4_Log_survey_residuals_from_GAM.pdf",width=7,height=7)
par(mfrow=c(1,1))
par(las=3)
par(mar=c(8,8,1,1),oma=c(3,0,0,0))
plot(x=trends0$survey,y=resid(m01a),ylab="Residuals",cex.axis=0.7,pch=19,col=1)
dev.off()
par<-old.par 
par(las=1)



#Figure S5 QQ-plot

par(mfrow=c(2,2))
gam.check(m01a, col=trends0$mlgr)
pdf("FigS5_QQplot_from_GAM.pdf")
par(mfrow=c(2,2))
gam.check(m01a, col=trends0$mlgr)
dev.off()

# Figure S6 Log residuals from GAM model
par(mfrow=c(3,3))
plot(x=trends0$siz,y=resid(m01a),xlab="ln(max length)",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=fitted(m01a),y=resid(m01a),xlab="Fitted values",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$sst,y=resid(m01a),xlab="Sea surface temperature",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$aswept,y=resid(m01a),xlab="Area swept",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$density,y=resid(m01a),xlab="Density",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$depth,y=resid(m01a),xlab="Depth",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$npp,y=resid(m01a),xlab="Net primary production",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$mesh,y=resid(m01a),xlab="Mesh size",ylab="Residuals",pch=19,col=trends0$mlgr)


pdf("FigS6_Log_residuals_from_GAM.pdf",width=7,height=8)
par(mfrow=c(3,3))
plot(x=trends0$siz,y=resid(m01a),xlab="ln(max length)",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=fitted(m01a),y=resid(m01a),xlab="Fitted values",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$sst,y=resid(m01a),xlab="Sea surface temperature",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$aswept,y=resid(m01a),xlab="Area swept",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$density,y=resid(m01a),xlab="Density",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$depth,y=resid(m01a),xlab="Depth",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$npp,y=resid(m01a),xlab="Net primary production",ylab="Residuals",pch=19,col=trends0$mlgr)
plot(x=trends0$mesh,y=resid(m01a),xlab="Mesh size",ylab="Residuals",pch=19,col=trends0$mlgr)
dev.off()

#Linearising the response of the sst model
m01al<-gam(nsp~s(sst,k=4)+s(lnaswept,k=4)+s(lndensity,k=4)+s(lndepth,k=4)+s(siz,k=4)+s(siz2,k=4)
          +s(lnnpp,k=4)+s(mesh,k=4,by=mlgr),family=nb(),method="REML",select="TRUE",data=trends0)
AIC(m01a,m01al)
summary(m01al)
concurvity(m01al)
sum(resid(m01al,type="pearson")^2)/m01al$df.res
par(mfrow=c(3,3))
plot(m01al,residuals=T,pch=1,scale=0)

#How welll does a linear version of the model work
m01all<-glm.nb(nsp~invsst+lnaswept+lndensity+lndepth+siz+siz2+lnnpp+mesh:mlgr,link="log",data=trends0)
summary(m01all)
sum(resid(m01all,type="pearson")^2)/m01all$df.res

vif(m01all)
par(mfrow=c(2,2))
plot(m01all,residuals=T,pch=1,scale=0)

quit()


