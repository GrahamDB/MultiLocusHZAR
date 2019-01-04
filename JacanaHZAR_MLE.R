## Estimate the Maximum Likelihood parameters for a given compiled model
cline.getMLE <- function(fitR, nFits=30,doPar=FALSE){
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
