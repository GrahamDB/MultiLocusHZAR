
##source("mleHelper.R")

cline.getTrace<- function(id,
                          model.selected=all.selected[[id]],
                          doPar=FALSE,tuneVarEdge=1,tuneCW=1){
  models=cline.getModels(id)
  cat("A")
  obsData=models$obsData
  for( m.id in model.selected){
      if("varL" %in% names(models$init[[m.id]]$modelParam))
      models$init[[m.id]]$modelParam$tune$varL   <-
          models$init[[m.id]]$modelParam$tune$varL * tuneVarEdge
      if("varR" %in% names(models$init[[m.id]]$modelParam))
      models$init[[m.id]]$modelParam$tune$varR   <-
          models$init[[m.id]]$modelParam$tune$varR * tuneVarEdge
      if("center" %in% names(models$init[[m.id]]$modelParam))
      models$init[[m.id]]$modelParam$tune$center   <-
          models$init[[m.id]]$modelParam$tune$center * tuneCW
      if("width" %in% names(models$init[[m.id]]$modelParam))
      models$init[[m.id]]$modelParam$tune$width   <-
          models$init[[m.id]]$modelParam$tune$width * tuneCW
    
  }
  init=hzar.multiFitRequest(models$init[model.selected],
                            skip=which(id==all.clines)*100);
  print(init)
  runInit=hzar.doFit.multi(init,doPar=doPar)
  cline.id=id;
  save(cline.id,model.selected,obsData,init,runInit,
       file=file.path("cline-models",paste(id,"Init.RData",sep="")))
  cat("B")
  if(doPar){
    chainInit=mclapply(runInit,
                       hzar.next.fitRequest,mc.cores=3)
  } else {
    chainInit=lapply(runInit,
                     hzar.next.fitRequest) 
  }
  chainInit=hzar.multiFitRequest(chainInit,
                                 each=3,
                                 baseSeed=NULL)
  chainRuns=hzar.doChain.multi(chainInit,doPar=doPar)
  
  save(cline.id,model.selected,obsData,chainInit,chainRuns,file=file.path("cline-models",paste(id,"Chains.RData",sep="")))
  cat("C")
}



cline.loadInit<- function(id) {
  local({
    load(file=file.path("cline-models",paste(id,"Init.RData",sep="")))
    list(cline.name=cline.id, m.names=model.selected,obsData=obsData,
         init=init,
         runs=runInit)
  
  })
  #load(file=file.path("cline-models",paste(id,"Chains.RData",sep="")))
}

cline.loadTrace<- function(id) {
  local({
    load(file=file.path("cline-models",paste(id,"Init.RData",sep="")))
    init
  }) -> models
  load(file=file.path("cline-models",paste(id,"Chains.RData",sep="")))
  list(cline.name=cline.id, m.names=model.selected,obsData=obsData,
       initDG=sapply(models,hzar.dataGroup.add,simplify=FALSE),
       runs=chainRuns)
}
cline.buildAnalysis <- function(id,doPar = TRUE) {
  local({
    cline.loadTrace(id)-> trace
      chainDG<-lapply(trace$runs,
                      function(x) hzar.dataGroup.add(x[[3]],doPar = doPar))
    oDG<-hzar.make.obsDataGroup(chainDG)
    oDG<-hzar.copyModelLabels(trace$initDG,oDG)
    oDG
  })->analysis.oDG
  
  save(analysis.oDG,file=file.path("cline-analysis",paste(id,"model-data.RData",sep="")))
}

cline.loadAnalysis <- function(id){
  load(file=file.path("cline-analysis",paste(id,"model-data.RData",sep="")))
  analysis.oDG
}

runChainsAndAnalyze <- function(clines=all.clines,
                                    build.trace=FALSE,
                                build.analysis=build.trace,
                                doPar=FALSE){
    if(build.trace){
        if(!doPar){
            tmp<-foreach(id=clines) %do% { cline.getTrace(id, doPar=FALSE); NULL}
        }else if(length(clines)<3){
            tmp<-foreach(id=clines) %do% { cline.getTrace(id, doPar=TRUE); NULL}
        }else {
            tmp<-foreach(id=clines) %dopar% { cline.getTrace(id, doPar=FALSE); NULL}
        }
    }
            
    if(build.analysis)
        tmp<-foreach(id=clines) %do% { cline.buildAnalysis(id,doPar = doPar); NULL}
    
    foreach(id=clines,.combine=rbind) %do% {
        load(file=file.path("cline-analysis",paste(id,"model-data.RData",sep="")))
        res <- hzar.getLLCutParam (analysis.oDG$data.groups,c("center", "width"))
        rownames(res) <- paste(id,rownames(res),sep=".")
        res
    }
}
