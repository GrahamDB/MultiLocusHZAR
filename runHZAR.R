## Make locations to store output

if(!dir.exists("cline-test"))
  dir.create("cline-test")
if(!dir.exists("cline-models"))
  dir.create("cline-models")
if(!dir.exists("cline-analysis"))
  dir.create("cline-analysis")


## Load general parameter about zone,
## such as traits and loci to analyze,
## range of distances for sampling sites
source("params.R")

## Generic helper code to use subplex package to find MLE
## as well as select models (currently set to all models
## within 2 AIC of model with lowest AIC for a given trait
## or locus)
source("mleHelper.R")

## Specialized helper code construct models (including
## estimation of ML parameters) and save for deeper analysis
source("modelLoader.R")


if(FALSE){
  ## Try one (to make sure code works before entering parallel code)
  res <- cline.loadLocusModels(loci[1])
  print(tmp <- do.call(data.frame,c(as.list(res),row.names=loci[1])))
  tmp<- rbind(tmp,
              foreach(id=loci[-1],.combine=rbind) %dopar% {
                res <-  cline.loadLocusModels(id);
                do.call(data.frame,c(as.list(res),row.names=id))
              })
  tmp$dAIC=tmp$null.model- apply(tmp[,1:15],1,min)
  tmp$isClinal=tmp$dAIC>2
  tmp$best.model <- apply(tmp[,1:15],1,function(x,lN) lN[which.min(x)],lN=colnames(tmp))
  tmp$center <- foreach(m=tmp$best.model,id=rownames(tmp),.combine=c) %do%
    cline.getModels(id)$MLEs[[m]]$param.all$center
  tmp$width <- foreach(m=tmp$best.model,id=rownames(tmp),.combine=c) %do%
    cline.getModels(id)$MLEs[[m]]$param.all$width
  summary(tmp)
  densityplot(~center+width,data=tmp)
}
if(FALSE){
  # Q.A.QT<-do.call(data.frame,c(as.list(cline.loadPhenModels(traits[1])),row.names=traits[1]))
  
  ## Find MLE for QT traits (make sure params.R and modelLoader.R are correct and sourced!)
  ## (make sure mleHelper.R is sourced as well)
  Q.A.QT<- #rbind(Q.A.QT,
    foreach(id=traits,.combine=rbind) %dopar% {
      res <-  cline.loadPhenModels(id);
      do.call(data.frame,c(as.list(res),row.names=id))
    }#)
} else {
  
  Q.A.QT<- foreach(id=traits,.combine=rbind) %do% {
    res<-cline.getModels(id)$AIC
    do.call(data.frame,c(as.list(res),row.names=id))
  }
  
}

postProcQT <- function(Q.A.QT){
  Q.A.QT$best.model <- apply(Q.A.QT[,1:10],1,function(x,lN) lN[which.min(x)],lN=colnames(Q.A.QT))
  Q.A.QT$center <- foreach(m=Q.A.QT$best.model,id=rownames(Q.A.QT),.combine=c) %do%
    cline.getModels(id)$MLEs[[m]]$param.all$center
  Q.A.QT$width <- foreach(m=Q.A.QT$best.model,id=rownames(Q.A.QT),.combine=c) %do%
    cline.getModels(id)$MLEs[[m]]$param.all$width
  Q.A.QT$logLike <- foreach(m=Q.A.QT$best.model,id=rownames(Q.A.QT),.combine=c) %do%
    cline.getModels(id)$MLEs[[m]]$logLike
  print(Q.A.QT)
  Q.A.QT
}

Q.A.QT <- postProcQT(Q.A.QT)

## Plotting Trait Clines
## Add space for right axis
par(mar=c(5, 4, 4, 3) + 0.1)
hzar.plot.obsData(cline.getModels(traits[3])$obsData,type="n"  ,ylim=c(-2,2))
foreach(id=traits[-5],m=Q.A.QT$best.model[-5],c.col=c("red","orange", "blue","purple")) %do% hzar.plot.cline(cline.getModels(id)$MLEs[[m]],add=TRUE,col=c.col)
old.xlim=par("usr")[1:2]
par(new=TRUE,usr=c(old.xlim,525,650));
hzar.plot.cline(cline.getModels(traits[5])$MLEs[["free.none"]],ylim=c(437.225,801.875),add=TRUE,lty="dashed",xlim=c(0,1187.49746))
axis(side=4)


#Updating clines to be characteristic ascending
par(mar=c(5, 4, 4, 3) + 0.1)
hzar.plot.obsData(cline.getModels(traits[3])$obsData,type="n"  ,ylim=c(0,1), ylab="Q")
plotMLEs = foreach(id=traits,m=Q.A.QT$best.model,.combine=c) %do% { res=list(cline.getModels(id)$MLEs[[m]]); names(res)=id; res;}
#$Q.A
#$Q.diagnostic
plotMLEs$Comp.2$clineFunc = function(x) (1/(1 + exp(-4 * (x - 285.005002253929)/6.89580973160708)))
plotMLEs$Comp.3$clineFunc = function (x)    (1/(1 + exp(-4 * (x - 595.495798640652)/803.989052915233)))
plotMLEs$WHLEN$clineFunc = function (x) (1/(1 + exp(-4 * (x - 283.162409728904)/158.553100921555)))
qt.col=c("red","orange", "blue","purple","darkgreen")
qt.lty=c("solid","solid","solid","solid","solid" )
## Modify to change legend labels (do not change order here!)
qt.traits=c("Q.A", "Q.diagnostic", "Comp.2", "Comp.3", "WHLEN") 
## example:
## qt.traits=c("Q SNP", "Q Diagnostic", "Song PC2", "Song PC3", "Whistle Length") 
## Modify to change legend order
qtlegend.order=c(1, 2, 3, 4, 5)
legend("bottomright",legend=qt.traits[qtlegend.order],col=qt.col[qtlegend.order],lty=qt.lty[qtlegend.order])

## Modify to change appearance order
qtshow.order=c(1, 2, 3, 4, 5)
## example, centers left to right
## qtshow.order=c(2, 1, 4, 5, 3)

foreach(id=traits[qtshow.order],c.col=qt.col[qtshow.order], c.lty=qt.lty[qtshow.order]) %do% {
  print(hzar.plot.cline(plotMLEs[[id]],add=TRUE,col=c.col,lty=c.lty))
  dev.print(pdf,file=paste("iterativeMLEclines20160120legend_",id,".pdf",sep=""),width=7,height=7,pointsize=12);
}




## ## Try one (to make sure code works before entering parallel code)
## cline.loadLocusModels(loci[1])
## ##cline.loadPhenModels(traits[1])

## ##tmp<- foreach(id=traits[-1]) %dopar% { cline.loadPhenModels(id); }
## tmp<- foreach(id=loci[-1]) %dopar% { cline.loadLocusModels(id); }


## Helper code to run traces and estimate ranges of center and width
## for each selected model.
source("mcmcHelper.R")
if(FALSE){
  all.selected <- getSelectedModels(traits)
  unix.time({
    runChainsAndAnalyze(clines=traits,
                        build.trace=TRUE,
                        build.analysis=TRUE,
                        doPar=TRUE)->res
  })
  print(res)
  
  hzar.plot.cline(Q.A.analysis <- cline.loadAnalysis("Q.A"))
  
  hzar.plot.cline(Q.A.analysis, ylab="Q value",xlab="Distance from San Francisco (km)",pch=16)
  site.map=loc
  rownames(site.map) <- as.character(site.map$Site)
  site.map$top.label.y=1.055; site.map[c("D","E","M","J"),"top.label.y"]=1.005;
  obs.dummy<-Q.A.analysis$obsData
  
  ##wcs.plot.sites <- function(obs.d){
  hzar.plot.obsData(Q.A.analysis,pch=16,type="n",
                    ylim=extendrange(c(0,1)),
                    yaxp=c(0,1,5),mgp=c(2,0.7,0),
                    ylab="Q value",xlab="Transect Distance (km)")
  tck.sz=0.8
  key=c(2,4,6,8,9,11,12,14:17)
  axis(at= site.map[key,]$distance,side=3,mgp=c(3,0.7,0),
       adj=c(0),padj=0,las=0,cex.axis=tck.sz,col="black",
       labels=as.numeric(site.map[key,]$Site)) #} 
  axis(at= site.map[-key,]$distance,side=3,mgp=c(3,1.4,0),
       adj=c(0),padj=0,las=0,cex.axis=tck.sz,col="black",
       labels=as.numeric(site.map[-key,]$Site)) 
  axis(at= site.map[4,]$distance,side=3,mgp=c(3,0.7,0),
       adj=c(0),padj=0,las=0,cex.axis=tck.sz,col="black",
       labels=as.numeric(site.map[4,]$Site))
  axis(at= site.map[5,]$distance,side=3,mgp=c(3,1.4,0),
       adj=c(0),padj=0,las=0,cex.axis=tck.sz,col="black",
       labels=as.numeric(site.map[5,]$Site))#} 
  #hzar.plot.fzCline(Q.A.analysis$data.groups$fixed.none,add=TRUE,pch=16,fzCol = hsv(0,0,0.7,1))
  hzar.plot.fzCline(Q.A.analysis$data.groups$fixed.right,add=TRUE,pch=16,fzCol = hsv(0,0,0.7,1))
  dev.print(pdf,file="AncestralProportionClinePlot20160108.pdf",height=8,width=8,pointsize=12)
  # }
  #   wcs.plot.sites(obs.dummy)
  # hzar.plot.cline(Q.A.analysis$data.groups$fixed.none$ML.cline, add=TRUE)
}

## MCMC analysis of loci (a bit complicated for odd reasons)
if(FALSE){
  n.selected <- sapply(getSelectedModels(loci),length)
  loci.single <- names( which(n.selected==1))
  
  ## Try one (to make sure code works before entering parallel code)
  all.selected <- getSelectedModels(loci)
  runChainsAndAnalyze(clines=loci.single[1],
                      build.trace=TRUE,
                      build.analysis=TRUE,
                      doPar=FALSE) -> res
  
  unix.time({
    rbind(res,
          runChainsAndAnalyze(clines=loci.single[2:10],
                              build.trace=TRUE,
                              build.analysis=TRUE,
                              doPar=TRUE)) -> res })
  res
  unix.time({
    rbind(res,
          runChainsAndAnalyze(clines=loci.single[-(1:10)],
                              build.trace=TRUE,
                              build.analysis=TRUE,
                              doPar=TRUE)) -> res })
  summary(res)
}


