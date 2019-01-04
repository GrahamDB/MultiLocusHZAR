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

## Test Loci against Null Model (Code to drop non clinal loci)
rebuild.tests <- FALSE
rebuild.clinal.loci <- rebuild.tests
if(rebuild.clinal.loci){
  ## Try one (to make sure code works before entering parallel code)
  ## Uncomment the following lines 
  
  ## res <- cline.testLocusModels2(loci[1])
  ## print(tmp <- do.call(data.frame,c(as.list(res),row.names=loci[1])))
  ##cline.loadPhenModels(traits[1])
  
  ##tmp<- foreach(id=traits[-1]) %dopar% { cline.loadPhenModels(id); }
  ## tmp<- rbind(tmp,
  foreach(id=loci,.combine=rbind) %dopar% {
    res <-  cline.fitFastLocusModels(id, rebuild=rebuild.tests);
    do.call(data.frame,c(as.list(res),row.names=id))
  }  ->cline.test1
  foreach(id=loci,.combine=rbind) %dopar% {
    res <-  cline.testLocusModels2(id, rebuild=rebuild.tests);
    do.call(data.frame,c(as.list(res),row.names=id))
  }  ->cline.test2
  if(FALSE){
    ## reload from saved files
    foreach(id=loci,.combine=rbind) %dopar% {
      res <-  cline.getFastTest(id)$AIC;
      do.call(data.frame,c(as.list(res),row.names=id))
    }  ->cline.test1
    foreach(id=loci,.combine=rbind) %dopar% {
      res <-  cline.getFastTest2(id)$AIC;
      do.call(data.frame,c(as.list(res),row.names=id))
    }  ->cline.test2
  }
  cline.test1$dAIC=cline.test1$null.model-cline.test1$fixed.none
  cline.test1$isClinal=cline.test1$dAIC>2
  cline.test1$center <- cline.getParamDist(rownames(cline.test1),param="center")
  cline.test1$width <- cline.getParamDist(rownames(cline.test1),param="width")
  summary(cline.test1)
  loci.test1<-rownames(cline.test1)[which(cline.test1$isClinal)]
  summary(cline.test1[loci.test1,])
  
  write.csv(cline.test1,file="fast_cline_selection.csv")
  write.table(rownames(cline.test1)[which(cline.test1$isClinal)],
              file="clinal_loci_names.txt",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE)
  
  cline.test2$dAIC=cline.test2$null.model-cline.test2$free.none
  cline.test2$isClinal=cline.test2$dAIC>2
  cline.test2$center <- cline.getParamDist2(rownames(cline.test2),param="center")
  cline.test2$width <- cline.getParamDist2(rownames(cline.test2),param="width")
  summary(cline.test2)
  loci.test2<-rownames(cline.test2)[which(cline.test2$isClinal)]
  summary(cline.test2[loci.test2,])
  
  write.csv(cline.test2,file="fast_cline_selection2.csv")
  write.table(rownames(cline.test2)[which(cline.test2$isClinal)],
              file="clinal_loci_names2.txt",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE)
}

##Collect and update loci to run from null model test results
#if(FALSE){
loci.test1 <- read.table("clinal_loci_names.txt",stringsAsFactors=FALSE)[[1]]
loci.test2 <- read.table("clinal_loci_names2.txt",stringsAsFactors=FALSE)[[1]]
loci.all <- loci ## read.table("significate_loci_names.txt",stringsAsFactors=FALSE)[[1]]

summary(data.frame(T1Good    = loci.test1 %in% loci.all,
                   isT1subT2 = loci.test1 %in% loci.test2))
summary(data.frame(T2Good    = loci.test2 %in% loci.all,
                   isT2subT1 = loci.test2 %in% loci.test1))

loci <- sort(unique(c( loci.test1,loci.test2)))
if(FALSE){
  write.table(loci,
              file="clinal_loci.txt",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE)
}
}
loci <- read.table("clinal_loci.txt",stringsAsFactors=FALSE)[[1]]

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

##Ignore this code
if(FALSE){
  rbind.grow.helper <- function(aL,lenA=length(aL),fill.value=NA) {
    ## cat("H");
    if(lenA>10)
      return(rbind.grow(rbind.grow.helper(aL[1:floor(lenA/2)],
                                          lenA=floor(lenA/2),
                                          fill.value=fill.value),
                        rbind.grow.helper(aL[(floor(lenA/2)+1):lenA],
                                          lenA=lenA-floor(lenA/2),
                                          fill.value=fill.value) ));
    ## cat("L");
    if(lenA==1) return(aL[[1]]);
    nL <- unique(do.call(c,lapply(aL,names)));
    for(argI in 1:lenA) {
      ## cat(argI);
      if(!inherits(aL[[argI]],"data.frame" ))
        aL[[argI]] <- as.data.frame(aL[[argI]]);
      for(iter in nL[!nL%in%names(aL[[argI]])] ){
        ## cat(iter);
        aL[[argI]][[iter]] <- fill.value;
      }
      if(argI==1) {
        res <-  aL[[argI]]
      } else {
        res <- rbind(res,aL[[argI]][,names(res)])
      };
    } ;
    ## cat("|");
    res;
  };
  rbind.grow <- function(...,fill.value=NA){
    ## cat("G");
    aL <- list(...);lenA=length(aL);
    if(any(sapply(lapply(aL,names),is.null)))
      stop("argument with no names encountered.");
    if(lenA>10)
      return(rbind.grow(rbind.grow.helper(aL[1:floor(lenA/2)],
                                          lenA=floor(lenA/2),
                                          fill.value=fill.value),
                        rbind.grow.helper(aL[(floor(lenA/2)+1):lenA],
                                          lenA=lenA-floor(lenA/2),
                                          fill.value=fill.value) ));
    if(lenA==1) return(aL[[1]]);
    ## cat("P");
    nL <- unique(do.call(c,lapply(aL,names)));
    ## print(nL)
    for(argI in 1:lenA) {
      ## cat(argI);
      if(!inherits(aL[[argI]],"data.frame" ))
        aL[[argI]] <- as.data.frame(aL[[argI]]);
      ## print(names(aL[[argI]]))
      for(iter in nL[!nL%in%names(aL[[argI]])] ){
        ## cat(iter);
        aL[[argI]][,iter] <- fill.value;
      }
      ## print(names(aL[[argI]]))
      if(argI==1) {
        res <-  aL[[argI]]
      } else {
        res <- rbind(res,aL[[argI]][,names(res)])
      };
    } ;
    ## cat(">\n");
    res;
  }
  
  
  if(FALSE){
    cline.singleModelSummary <- function(id) {
      m <-  cline.getModels(id)$MLEs[[all.selected[[id]]]];
      data.frame(do.call(data.frame,m$param.all),
                 pLeft=m$clineFunc(0),
                 pRight=m$clineFunc(1200),
                 row.names=id);}
    single.res.tmp <- lapply(loci.single,cline.singleModelSummary)
    summary(single.res <- do.call(rbind.grow,single.res.tmp ))
    single.res$reversed <- (single.res$pLeft > single.res$pRight)
    cline.single.forward <- sapply(rownames(single.res)[!single.res$reversed],
                                   function(id) cline.getModels(id)$MLEs[[all.selected[[id]]]],
                                   simplify=FALSE)
    cline.single.reverse <- sapply(rownames(single.res)[single.res$reversed],
                                   function(id) cline.getModels(id)$MLEs[[all.selected[[id]]]],
                                   simplify=FALSE)
    cline.single <- cline.single.forward
    for( iter in names(cline.single.reverse) ) {
      cline.single[[iter]] <- cline.single.reverse[[iter]];
      body(cline.single[[iter]]$clineFunc ) <-
        bquote(1-.(body(cline.single.reverse[[iter]]$clineFunc)))
    }
    
    
    obs.dummy <- cline.getModels(loci[1])$obsData
    
    hzar.plot.obsData(obs.dummy,type="n",ylim=extendrange(c(0,1)))
    foo <- lapply(cline.single,hzar.plot.cline,add=TRUE,col=hsv(0,0,0,0.5))
    
    hzar.plot.obsData(obs.dummy,type="n",ylim=extendrange(c(0,1)))
    foo <- lapply(cline.single,hzar.plot.cline,add=TRUE,col=hsv(0,0,0,0.1))
    
    hzar.plot.obsData(obs.dummy,type="n",ylim=extendrange(c(0,1)),xlim=c(-600,1800))
    foo <- lapply(cline.single,hzar.plot.cline,add=TRUE,col=hsv(0,0,0,0.1))
    
  }
  
}

## Finding and Plotting Diagnostic Loci
if(FALSE) {
  ## Build set of top models:
  getLociMLE.best <- function(clines=all.clines,
                              restrict=c(0.2,0.8),
                              old.sel=FALSE){
    tmp<- foreach(id=clines,.combine=c) %dopar% {
      junk<-list() ;
      m <- cline.getModels(id)
      m.sel <-  names(which.min(m$AIC))
      if(length(m.sel>1))
        m.sel<-sample(m.sel,1)
      m <- m$MLEs[[m.sel]];
      attr(m,"pLeft")=m$clineFunc(0)
      attr(m,"pRight")=m$clineFunc(1200)
      if(attr(m,"pLeft") > attr(m,"pRight") ){ 
        body(m$clineFunc ) <-
          bquote(1-.(body(m$clineFunc)))
        attr(m,"cline.rev" ) <- TRUE
        if(old.sel)
          attr(m,"cline.plot") <- is.null(restrict)||
          attr(m,"pLeft" ) >= restrict[2] &&
          attr(m,"pRight" ) <= restrict[1]
      }else{
        attr(m,"cline.rev" ) <-  FALSE
        if(old.sel)
          attr(m,"cline.plot") <- is.null(restrict)||
            attr(m,"pLeft" ) <= restrict[1] &&
            attr(m,"pRight" ) >= restrict[2]
      }
      if(!old.sel)
        attr(m,"cline.plot") <- is.null(restrict)||
        m$param.all$pMin <= restrict[1] &&
        m$param.all$pMax >= restrict[2]
      if(is.na(attr(m,"cline.plot")))
        attr(m,"cline.plot") <- FALSE
      if(inherits(m,"hzar.cline")){
        junk[[id]] <- m
      } else {
        res=list()
        attr(res,"cline.plot") <- FALSE
        junk[[id]] <- res
      }
      junk
    }
  }
  
  getLociMLE.selected <- function(clines=all.clines,
                                  max.dAIC=2,restrict=c(0.2,0.8),old.sel=FALSE){
    foreach(id=clines,.combine=c) %dopar% {
      m <- cline.getModels(id)
      m.sel <-  names(which(m$AIC<=min(m$AIC+max.dAIC)))
      foreach(m.id=m.sel) %do% {
        mIter=m$MLEs[[m.id]]
        if( ( attr(mIter,"pLeft" ) <-mIter$clineFunc(0) )
            > ( attr(mIter,"pRight" ) <-mIter$clineFunc(1200)) ) {
          body(mIter$clineFunc ) <-
            bquote(1-.(body(mIter$clineFunc)))
          attr(mIter,"cline.rev" ) <- TRUE
          if(old.sel)
            attr(mIter,"cline.plot") <- is.null(restrict)||
            attr(mIter,"pLeft" ) >= restrict[2] &&
            attr(mIter,"pRight" ) <= restrict[1]
        } else{
          attr(mIter,"cline.rev" ) <-  FALSE
          if(old.sel)
            attr(mIter,"cline.plot") <- is.null(restrict)||
              attr(mIter,"pLeft" ) <= restrict[1] &&
              attr(mIter,"pRight" ) >= restrict[2]
        }
        if(!old.sel)
          attr(mIter,"cline.plot") <- is.null(restrict)||
            mIter$param.all$pMin <= restrict[1] &&
            mIter$param.all$pMax >= restrict[2]
        if(is.na(attr(mIter,"cline.plot")))
          attr(mIter,"cline.plot") <- FALSE
        attr(mIter,"n.eq") <- length(m.sel)
        mIter } -> junk
      names(junk) <- paste(id,m.sel,sep=".")
      junk
    }
  }
  groupLociMLE.sets <- function(clineMLEs,
                                plot=TRUE,
                                alpha.max=0.1,
                                col.v=0){
    print(n.eq.lvl <- sort(unique( sapply(clineMLEs,attr,which="n.eq"))))
    cline.sets <- foreach(key=n.eq.lvl,.combine=c ) %:%
      foreach(cline.model=clineMLEs,cline.id=names(clineMLEs),
              .combine=c,
              .final=function(myL) {
                res <- list(clineMLEs[myL]);
                names(res) <- paste("G",key,sep="");
                res; }) %:%
      when(attr(cline.model,"n.eq")==key) %do% {
        if(plot && (is.null(attr(cline.model,"cline.plot")) ||
                    attr(cline.model,"cline.plot")) )
          hzar.plot.cline(cline.model,
                          add=TRUE,
                          col=hsv(0,0,col.v,alpha.max/key));
        cline.id
      }
    if(length(cline.sets)<10)
      print(summary(cline.sets))
    invisible(cline.sets)
  }
  
  
  
  loci.clines.all <- getLociMLE.best(loci,restrict=c(0.2,0.8))
  loci.clines <- loci.clines.all[sapply(loci.clines.all,attr,which="cline.plot")]
  clines.selected <- getLociMLE.selected (loci,restrict=c(0.2,0.8))
  loci.models.sel <- names(clines.selected)[sapply(clines.selected,
                                                   attr,
                                                   which="cline.plot")]
  cat(paste("# Clinal loci:",length(loci),
            "# Diagnostic loci:",length(loci.clines),
            "# Candidate loci:",length(unique(gsub("\\..*$","",loci.models.sel))),
            "\n",sep="\t"))
  if(FALSE){
    write.table(names(loci.clines),
                file="diagnostic_loci.txt",
                col.names=FALSE,
                row.names=FALSE,
                quote=FALSE)
  }
  GenFreq[,c(2,3,6)] -> site.map
  rownames(site.map) <- as.character(site.map$Site)
  site.map$top.label.y=1.025; site.map[c("D","E","M","J"),"top.label.y"]=1.005;
  
  wcs.plot.sites <- function(){
    hzar.plot.obsData(obs.dummy,type="n",
                      ylim=extendrange(c(0,1)),xlim=c(-200,1400))
    segments(x0=site.map$distance,
             y0=1.001,y1=site.map$top.label.y,
             col="gray")
    for(iter in 1:2){
      text(x= site.map$distance,y=site.map$top.label.y,
           adj=c(0,NA),srt=45,cex=0.5,col="black",
           vfont=c("sans serif","bold"),
           labels=gsub(pattern="_",
                       replace=" ",
                       as.character(site.map$SiteName))) } }
  
  wcs.plot.sites()
  segments(x0= site.map$distance,y0=-0.001,y1=-0.025,col="gray")
  foo <- lapply(loci.clines.all,hzar.plot.cline,add=TRUE,col=hsv(0,0,0,0.05))
  
  dev.print(pdf,file="BestClineModelsFor2763Loci20160121.pdf",width=8,height=8)
  
  
  wcs.plot.sites()
  segments(x0= site.map$distance,y0=1,y1=-0.025,col="gray",lty="dotted")
  foo <- lapply(loci.clines,hzar.plot.cline,add=TRUE,col=hsv(0,0,0,0.1))
  
  dev.print(pdf,file="BestClineModelsFor241DiagnosticLoci20160121.pdf",width=8,height=8)
  
  
  wcs.plot.sites()
  segments(x0= site.map$distance,y0=1,y1=-0.025,col="gray",lty="dotted")
  loci.cline.sets <-groupLociMLE.sets(clines.selected,col.v=0,alpha.max=0.1)
  
  dev.print(pdf,file="SelectedClineModelsFor936CandidateLoci20160121.pdf",
            width=8,height=8)
  
  
}







if(FALSE){
  ##length(test3.selected <- cline.getDiagnostic(loci.all,restrict=c(0.3,0.7)))
  
  ## length(test3.selected <- cline.getDiagnostic(loci.all,restrict=c(0.3,0.7),r.range=c(-200,1400)))
  ## loci.test3 <- gsub(pattern="\\..*",replacement="",names(test3.selected[sapply(test3.selected,attr,which="cline.plot")]))
  ## names( loci.test3 ) <- loci.test3
  length(test3.selected <- cline.getDiagnostic(loci.all,restrict=c(0.3,0.7),r.range=c(-200,1400)));
  loci.test3 <- gsub(pattern="\\..*",replacement="",names(test3.selected[sapply(test3.selected,attr,which="cline.plot")]));
  names( loci.test3 ) <- loci.test3;
  summary(loci.test3 %in% loci)
  
  write.table(loci.test3,
              file="diagnostic_loci.txt",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE)
  
  summary(cline.MLE.tmp <- foreach(id=loci.test3[!loci.test3 %in% loci],.combine=rbind) %dopar% {res <-  cline.loadLocusModels(id);
  do.call(data.frame,c(as.list(res),row.names=id))
  });cline.MLE.tmp$dAIC=cline.MLE.tmp$null.model- apply(cline.MLE.tmp[,1:15],1,min);cline.MLE.tmp$isClinal=cline.MLE.tmp$dAIC>2;summary(cline.MLE.tmp )
  
  
  loci.test3 <- unique(loci.test3)
  loci.94clines <- getLociMLE.best(loci.test3,restrict=c(0.3,0.7))
  loci.94clines.diag <- loci.94clines[sapply(loci.94clines,attr,which="cline.plot")]
  loci.94clines.sig <- loci.94clines[cline.MLE.tmp [loci.test3,"isClinal"] ]
  loci.94clines.sig.diag <- loci.94clines.sig[sapply(loci.94clines.sig,attr,which="cline.plot")]
  length(loci.94clines.diag)
  length(loci.94clines.sig.diag)
  hzar.plot.obsData(obs.dummy,type="n",ylim=extendrange(c(0,1)),xlim=c(-600,1800))
  points(x=obs.dummy$frame$dist,y=y.jitter,
         pch=rownames(obs.dummy$frame),col="black",cex=0.5)
  points(x=obs.dummy$frame$dist,y=1.092+y.jitter,
         pch=rownames(obs.dummy$frame),col="black",cex=0.5)
  foo <- lapply(loci.94clines,hzar.plot.cline,add=TRUE,col=hsv(0,0,0.5,1))
  dev.print(pdf,file="BestClineModelsFor94Loci4.pdf",width=8,height=8)
  hzar.plot.obsData(obs.dummy,type="n",ylim=extendrange(c(0,1)),xlim=c(-600,1800))
  points(x=obs.dummy$frame$dist,y=y.jitter,
         pch=rownames(obs.dummy$frame),col="black",cex=0.5)
  points(x=obs.dummy$frame$dist,y=1.092+y.jitter,
         pch=rownames(obs.dummy$frame),col="black",cex=0.5)
  foo <- lapply(loci.94clines.diag,hzar.plot.cline,add=TRUE,col=hsv(0,0,0.5,1))
  dev.print(pdf,file="BestClineModelsFor7Loci4diag.pdf",width=8,height=8)
  hzar.plot.obsData(obs.dummy,type="n",ylim=extendrange(c(0,1)),xlim=c(-600,1800))
  points(x=obs.dummy$frame$dist,y=y.jitter,
         pch=rownames(obs.dummy$frame),col="black",cex=0.5)
  points(x=obs.dummy$frame$dist,y=1.092+y.jitter,
         pch=rownames(obs.dummy$frame),col="black",cex=0.5)
  foo <- lapply(loci.94clines.sig,hzar.plot.cline,add=TRUE,col=hsv(0,0,0.5,1))
  dev.print(pdf,file="BestClineModelsFor69Loci4sig.pdf",width=8,height=8)
  
  clines.selected <- getLociMLE.selected (loci,restrict=c(0.3,0.7))
  hzar.plot.obsData(obs.dummy,type="n",ylim=extendrange(c(0,1)),xlim=c(-600,1800))
  ## abline(v=obs.dummy$frame$dist,lty="dotted",col="red",lwd=1)
  points(x=obs.dummy$frame$dist,y=y.jitter,
         pch=rownames(obs.dummy$frame),col="black",cex=0.5)
  points(x=obs.dummy$frame$dist,y=1.092+y.jitter,
         pch=rownames(obs.dummy$frame),col="black",cex=0.5)
  loci.cline.sets <-groupLociMLE.sets(clines.selected,col.v=0.5,alpha.max=1)
  
  dev.print(pdf,file="SelectedClineModelsForXXLoci3.pdf",width=8,height=8)
  
}
