
##source("params.R")




temp.getMLE.old <- function(fitR){
  mP <- fitR$modelParam
  res <- subplex(mP$init,
                 function(theta) -fitR$llFunc(theta),
                 control=list(parscale=sapply(names(mP$init),
                                              function(x) (mP$upper[[x]]-mP$lower[[x]])/50)))
  hzar.gen.cline(res$par, fitR)
}

temp.getMLE <- function(fitR, nFits=30,doPar=FALSE){
    mP <- fitR$modelParam
    resF <- subplex(mP$init,
                 function(theta) -fitR$llFunc(theta),
                 control=list(parscale=sapply(names(mP$init),
                                              function(x) (mP$upper[[x]]-mP$lower[[x]])/50)))
    res <- list(hzar.gen.cline(resF$par, fitR))
    if(doPar ){
        res <- c(res,times(nFits-1) %dopar% {
            tmp <- subplex(temp.safe.init(fitR),
                           function(theta) -fitR$llFunc(theta),
                           control=list(parscale=sapply(names(mP$init),
                                            function(x) (mP$upper[[x]]-mP$lower[[x]])/50)))
            list(hzar.gen.cline(tmp$par, fitR)) })
    } else {
        res <- c(res,times(nFits-1) %do% {
            tmp <- subplex(temp.safe.init(fitR),
                           function(theta) -fitR$llFunc(theta),
                           control=list(parscale=sapply(names(mP$init),
                                            function(x) (mP$upper[[x]]-mP$lower[[x]])/50)))
            list(hzar.gen.cline(tmp$par, fitR) )})
    }
    res[[which.max(sapply(res,function(x) x$logLike))[1]]]
}
   

temp.getNULL <- function(obsData){
    pExpected <- sum(obsData$frame$n * obsData$frame$obsFreq)/sum(obsData$frame$n)
    free.parameters=list(pVal=pExpected);
    cline.func <- eval(substitute(function(x) rep(p, length(x)), list(p = pExpected)))
    hzar.make.cline(free.parameters,  parameters=free.parameters,
                    func=cline.func, LL = obsData$model.LL(cline.func) )
}
                    


cline.getModels<-function(id){
    load(file=file.path("cline-models",paste(id,"MLE.RData",sep="")))
    list(obsData=obsData,models=clineModels,init=clineInit,MLEs=clineMLE,AIC=clineAIC)
}

cline.getFastTest <- function(id){
    load(file=file.path("cline-test",paste(id,"NULLtest.RData",sep="")))
    list(obsData=obsData,models=clineModels,init=clineInit,MLEs=clineMLE,AIC=clineAIC)
}

cline.getFastTest2 <- function(id){
    load(file=file.path("cline-test",paste(id,"NULLtest2.RData",sep="")))
    list(obsData=obsData,models=clineModels,init=clineInit,MLEs=clineMLE,AIC=clineAIC)
}

cline.getParamDist <- function(idL, param="center"){
    sapply(idL,function(id) cline.getFastTest(id)$MLEs$fixed.none$param.all[[param]])
}                                                    

cline.getParamDist2 <- function(idL, param="center"){
    sapply(idL,function(id) cline.getFastTest2(id)$MLEs$free.none$param.all[[param]])
}  
getSelectedModels <- function(clines=all.clines){
    tmp<- foreach(id=clines,.combine=c) %dopar% {
        junk<-list() ;
        junk[[id]]<-cline.getModels(id)$AIC;
        junk }
    sapply(tmp, function(x) names(which(x<min(x+2))),simplify=FALSE)
}



cline.getDiagnostic <- function(clines=all.clines,restrict=c(0.3,0.7),r.range=c(0,900)){
    foreach(id=clines,.combine=c) %dopar% {
        m <- cline.getFastTest2(id)
        m.ex <-  cline.getFastTest(id)
        m$MLEs <- c(m$MLEs,m.ex$MLEs["fixed.none"])
            m.sel <-  names(m$MLEs)
            noRestrict <-  is.null(restrict)
            if(is.null(r.range))
                r.range=c(0,1200)
            foreach(m.id=m.sel) %do% {
                mIter=m$MLEs[[m.id]]
                if( ( attr(mIter,"pLeft" ) <-mIter$clineFunc(min(r.range)) )
                   > ( attr(mIter,"pRight" ) <-mIter$clineFunc(max(r.range))) ) {
                    body(mIter$clineFunc ) <-
                        bquote(1-.(body(mIter$clineFunc)))
                    attr(mIter,"cline.rev" ) <- TRUE
                    attr(mIter,"cline.plot") <- noRestrict ||
                        attr(mIter,"pLeft" ) >= restrict[2] &&
                            attr(mIter,"pRight" ) <= restrict[1]
                } else{
                    attr(mIter,"cline.rev" ) <-  FALSE
                    attr(mIter,"cline.plot") <- noRestrict ||
                        attr(mIter,"pLeft" ) <= restrict[1] &&
                            attr(mIter,"pRight" ) >= restrict[2]
                }
                attr(mIter,"n.eq") <- length(m.sel)
                mIter } -> junk
            names(junk) <- paste(id,m.sel,sep=".")
            junk
        }
}


temp.random.init <-
    function(modelParam) {
        sapply(names(modelParam$init),
               function(p) {
                   if(p %in% c("varR","varL") )
                       return(modelParam$init[[p]])
                   runif(1,modelParam$lower[[p]],modelParam$upper[[p]] )
               },simplify=FALSE)
    }
temp.safe.init <- function(fitR,init=temp.random.init(fitR$modelParam),nAttempts=10){
    mLL=-1e8
    try(mLL <- fitR$llFunc(init) )
    while(mLL <= -1e8 && nAttempts >0){
        nAttempts=nAttempts-1
        init <- temp.random.init(fitR$modelParam)
        try(mLL <- fitR$llFunc(init) )
    }
    init
}
    

cline.MLE.support <- function(fitR, nFits=30,doPar=TRUE,compact=FALSE){
    mP <- fitR$modelParam
    res <- NULL
    if(doPar ){
        res <- times(nFits) %dopar% {
            tmp <- subplex(temp.safe.init(fitR),
                           function(theta) -fitR$llFunc(theta),
                           control=list(parscale=sapply(names(mP$init),
                                            function(x) (mP$upper[[x]]-mP$lower[[x]])/50)))
            list(hzar.gen.cline(tmp$par, fitR)) }
    } else {
        res <- times(nFits) %do% {
            tmp <- subplex(temp.safe.init(fitR),
                           function(theta) -fitR$llFunc(theta),
                           control=list(parscale=sapply(names(mP$init),
                                            function(x) (mP$upper[[x]]-mP$lower[[x]])/50)))
            list(hzar.gen.cline(tmp$par, fitR) )}
    }
    if(!compact)
        return(res)
    resC <- foreach(MLE=res,.combine=rbind) %do% {
        data.frame(as.list(MLE$param.free), MLE["logLike"]) }
    return(resC )
}
    
cline.MLE.check <- function(id,m.name,...){
    m <- cline.getModels(id)
    if(any(m.name %in% "all"))
        m.name <- names( m$init)
    if(length(m.name ) ==1){
        mF <- m$init[[m.name]]
        print(summary(res <- cline.MLE.support(mF,compact=TRUE,...)))
        print(summary(res$logLike - m$MLEs[[m.name]]$logLike))
        return(invisible(res ))
    }
    res <- foreach(m.id=m.name,.combine=c)%do%{
        mF <- m$init[[m.id]]
        tmp <- list()
        tmp[[m.id]] <- cline.MLE.support(mF,compact=TRUE,...)
        tmp
    }
    resS <- foreach(m.id=m.name,.combine=cbind)%do%{
        tmp <- list()
        tmp[[m.id]] <- res[[m.id]]$logLike - m$MLEs[[m.id]]$logLike
        as.data.frame(tmp )}
    print(summary(resS ))
    return(invisible(res ))
}
