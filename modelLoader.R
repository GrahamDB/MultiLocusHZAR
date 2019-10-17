
summary(Phen <- read.csv("jac_K2_Qst.csv", row.names=1))

## Sites A, B form the "left" side of the cline, so pull the
## initial values from there.
obsData.getMeanLeft <- function(trait,
                                obsData=cline.phenObsData(trait))
  mean(na.omit(Phen[[trait]][Phen$SiteID %in% c("A","B")]))
  
obsData.getVarLeft <- function(trait,
                               obsData=cline.phenObsData(trait))
  var(na.omit(Phen[[trait]][Phen$SiteID %in% c("A","B")]))

## Sites ZC, ZG form the "right" side of the cline, so pull the
## initial values from there.
obsData.getMeanRight <- function(trait,
                                obsData=cline.phenObsData(trait))
      mean(na.omit(Phen[[trait]][Phen$SiteID %in% c("ZC","ZG")]))
obsData.getVarRight <- function(trait,
                               obsData=cline.phenObsData(trait))
    var(na.omit(Phen[[trait]][Phen$SiteID %in% c("ZC","ZG")]))




## Genotypic (aka frequency) models
GenFreq<-read.csv("table_hwe_cull_076750_jac_SNP_Freq.csv",row.names=1)

#Generate loci names
#colnames(GenFreq)



locus_key <- data.frame(    #can also be generated from colnames(GenFreq)
  loci_SNP = c("S1_593993909_G","S1_308056105_T" , "S1_308117620_A" , "S1_200011937_G" , "S1_200030362_G" , "S1_200030392_G",
               "S1_360310864_A" , "S1_378861502_C" , "S1_378885504_C" , "S1_378927984_C" ,  "S1_379048952_G" , "S1_1001972416_C", 
               "S1_15585748_C"  , "S1_15586630_C" , "S1_219827999_G" , "S1_219862918_T" , "S1_219862925_C" , "S1_219863010_G" , 
               "S1_219863019_T" , "S1_219863059_G" , "S1_219863081_C" ,  "S1_219863089_G" ,  "S1_219863121_C" ,
               "S1_219863124_T" , "S1_709471828_C" , "S1_709471832_T" ,"S1_500438201_G" , "S1_500440953_T" , "S1_500440954_C" ,
               "S1_500440986_C" , "S1_500458867_A" , "S1_500458877_G" , "S1_500477580_G" , "S1_1002907772_A", "S1_178131011_C" , 
               "S1_239440237_T" , "S1_239451071_A" , "S1_484964082_T" , "S1_938036466_G" , "S1_272303766_C" , "S1_272320592_C" ,
               "S1_442143779_G" , "S1_442184589_A" ,  "S1_442184604_A" , "S1_442207270_G" , "S1_442255483_G" , "S1_442255566_G" ,
               "S1_442255568_A" , "S1_442255583_T" , "S1_442270424_A" , "S1_744502908_G" , "S1_825984260_A" , "S1_825987966_G" , 
               "S1_825987969_T" ,  "S1_825997510_T" , "S1_825997532_C" , "S1_825997542_C",  "S1_825997612_G" , "S1_825997691_A" , 
               "S1_825997696_A" , "S1_825997728_T" , "S1_825997745_A" , "S1_825997752_C",  "S1_825997779_T" , "S1_825998212_G" , 
               "S1_825998215_G" , "S1_825998216_T" , "S1_825998237_G" ,"S1_26838967_C"),
 
  loci_CNT = c("S1_593993909_N" , "S1_308056105_N" , "S1_308117620_N" , "S1_200011937_N" , "S1_200030362_N" , "S1_200030392_N" , 
               "S1_360310864_N" , "S1_378861502_N" , "S1_378885504_N" , "S1_378927984_N" , "S1_379048952_N", "S1_1001972416_N" , 
               "S1_15585748_N" , "S1_15586630_N" , "S1_219827999_N" ,  "S1_219862918_N" , "S1_219862925_N" ,"S1_219863010_N" , 
               "S1_219863019_N" , "S1_219863059_N" ,"S1_219863081_N" ,  "S1_219863089_N" , "S1_219863121_N" ,
               "S1_219863124_N" , "S1_709471828_N" , "S1_709471832_N" , "S1_500438201_N" , "S1_500440953_N" , "S1_500440954_N" ,
               "S1_500440986_N", "S1_500458867_N" , "S1_500458877_N" , "S1_500477580_N" , "S1_1002907772_N",  "S1_178131011_N" ,
               "S1_239440237_N" , "S1_239451071_N" , "S1_484964082_N" , "S1_938036466_N" , "S1_272303766_N" , "S1_272320592_N" ,
               "S1_442143779_N" , "S1_442184589_N" , "S1_442184604_N" , "S1_442207270_N" , "S1_442255483_N" , "S1_442255566_N" ,
               "S1_442255568_N" , "S1_442255583_N" , "S1_442270424_N" , "S1_744502908_N" , "S1_825984260_N" , "S1_825987966_N" ,
               "S1_825987969_N" , "S1_825997510_N" , "S1_825997532_N" ,"S1_825997542_N", "S1_825997612_N" ,  "S1_825997691_N" ,
               "S1_825997696_N" , "S1_825997728_N" , "S1_825997745_N", "S1_825997752_N", "S1_825997779_N" ,  "S1_825998212_N" ,
               "S1_825998215_N" , "S1_825998216_N" , "S1_825998237_N", "S1_26838967_N" ),
  row.names=loci)

cline.buildObsFreq<-function(locus) {
  hzar.doMolecularData1DPops(GenFreq$distance,
                             GenFreq[[locus_key[locus,"loci_SNP"]]],
                             GenFreq[[locus_key[locus,"loci_CNT"]]],
                             siteID=row.names(GenFreq))
    }


cline.buildObsNormal <- function(trait) {
  hzar.doNormalData1DRaw(hzar.mapSiteDist(siteID = rownames(GenFreq), 
                                          distance = GenFreq$distance), # Use the unique site ids and distances from the molecular table
                         traitSite = Phen$SiteID,  # Column of each individual's site. Will only use sites from the map 
                         traitValue = Phen[[trait]])
}

#Load data from already existing file
#Picking an allele for a locus

#Load Data - first create file



cline.fitFastLocusModels <- function(locus,obsData=cline.buildObsFreq(locus),rebuild=FALSE) {
  if(!rebuild&& isTRUE(0==as.numeric(file.access( file.path("cline-test",paste(locus,"NULLtest.RData",sep=""))))) ){
  
    cat("Found Saved File\n")
    load(file=file.path("cline-test",paste(locus,"NULLtest.RData",sep="")))
    
    print(do.call(data.frame,c(as.list(clineAIC),row.names=locus)))
    return(clineAIC)
  }
  
  clineModels<-list()
  cline.loadModel <- function(scaling,tails,
                                  id=paste(scaling,tails,sep="."))
    clineModels[[id]] <<- hzar.makeCline1DFreq(obsData, scaling, tails)
  cline.loadModel("fixed","none")
  ##cline.loadModel("free", "none")
  
  ## Modify all models to focus on the region where the observed
  ## data were collected.
  clineModels <- sapply(clineModels,
                                        hzar.model.addBoxReq,
                                        far.left , far.right,
                                        simplify=FALSE)
#This is where we're building the MCMC functions to the model
  ## Compile each of the models to prepare for fitting
  ## cat("Compile.\n")
  clineInit <- sapply(clineModels,
           hzar.first.fitRequest.old.ML,
           obsData=obsData,
           verbose=FALSE,
           simplify=FALSE)
  clineMLE <- sapply(clineInit, temp.getMLE, simplify=FALSE)
  cat("\n")
  ## cat("Null model.\n")
  clineMLE <- c(clineMLE ,list(null.model=temp.getNULL(obsData)))
  ##str(clineMLE)
  ## cat("AIC score:\n")
  ##print(
  clineAIC<- sapply( clineMLE, hzar.AIC.hzar.cline)
  ##)
  print(do.call(data.frame,c(as.list(clineAIC),row.names=locus)))
  save(obsData,clineModels,clineInit,clineMLE,clineAIC,file=file.path("cline-test",paste(locus,"NULLtest.RData",sep="")))
  clineAIC
}
## cline.fitFastLocusModels("S2741018")

cline.testLocusModels2 <- function(locus,obsData=cline.buildObsFreq(locus),rebuild=FALSE) {
  if(!rebuild&& isTRUE(0==as.numeric(file.access( file.path("cline-test",paste(locus,"NULLtest2.RData",sep=""))))) ){
  
    cat("Found Saved File\n")
    load(file=file.path("cline-test",paste(locus,"NULLtest2.RData",sep="")))
    
    print(do.call(data.frame,c(as.list(clineAIC),row.names=locus)))
    return(clineAIC)
  }
  
  clineModels<-list()
  cline.loadModel <- function(scaling,tails,
                                  id=paste(scaling,tails,sep="."))
    clineModels[[id]] <<- hzar.makeCline1DFreq(obsData, scaling, tails)
  cline.loadModel("free","none")
  ##cline.loadModel("free", "none")
  
  ## Modify all models to focus on the region where the observed
  ## data were collected.
  clineModels <- sapply(clineModels,
                                        hzar.model.addBoxReq,
                                        far.left , far.right,
                                        simplify=FALSE)
#This is where we're building the MCMC functions to the model
  ## Compile each of the models to prepare for fitting
  ## cat("Compile.\n")
  clineInit <- sapply(clineModels,
           hzar.first.fitRequest.old.ML,
           obsData=obsData,
           verbose=FALSE,
           simplify=FALSE)
  clineMLE <- sapply(clineInit, temp.getMLE, simplify=FALSE)
  cat("\n")
  ## cat("Null model.\n")
  clineMLE <- c(clineMLE ,list(null.model=temp.getNULL(obsData)))
  ##str(clineMLE)
  ## cat("AIC score:\n")
  ##print(
  clineAIC<- sapply( clineMLE, hzar.AIC.hzar.cline)
  ##)
  print(do.call(data.frame,c(as.list(clineAIC),row.names=locus)))
  save(obsData,clineModels,clineInit,clineMLE,clineAIC,file=file.path("cline-test",paste(locus,"NULLtest2.RData",sep="")))
  clineAIC
}

cline.loadLocusModels <- function(locus,obsData=cline.buildObsFreq(locus),rebuild=FALSE) {
  if(!rebuild&& isTRUE(0==as.numeric(file.access( file.path("cline-models",paste(locus,"MLE.RData",sep="")) ))) ){
  
    cat("Found Saved File\n")
    load(file=file.path("cline-models",paste(locus,"MLE.RData",sep="")))
    
    print(do.call(data.frame,c(as.list(clineAIC),row.names=locus)))
    return(clineAIC)
  }
  clineModels<-list()
  cline.loadModel <- function(scaling,tails,
                                  id=paste(scaling,tails,sep="."))
    clineModels[[id]] <<- hzar.makeCline1DFreq(obsData, scaling, tails)
  cline.loadModel("none", "none")
  cline.loadModel("fixed","none")
  cline.loadModel("free", "none")
  cline.loadModel("none", "mirror")
  cline.loadModel("fixed","mirror")
  cline.loadModel("free", "mirror")
  cline.loadModel("none", "right")
  cline.loadModel("fixed","right")
  cline.loadModel("free", "right")
  cline.loadModel("none", "left")
  cline.loadModel("fixed","left")
  cline.loadModel("free", "left")
  cline.loadModel("none", "both")
  cline.loadModel("fixed","both")
  cline.loadModel("free", "both")
  
  ## Modify all models to focus on the region where the observed
  ## data were collected.
  clineModels <- sapply(clineModels,
                                        hzar.model.addBoxReq,
                                        far.left , far.right,
                                        simplify=FALSE)
  #This is where we're building the MCMC functions to the model
  ## Compile each of the models to prepare for fitting
  clineInit <- sapply(clineModels,
           hzar.first.fitRequest.old.ML,
           obsData=obsData,
           verbose=FALSE,
           simplify=FALSE)
  clineMLE <- sapply(clineInit, temp.getMLE, simplify=FALSE)
  cat("\n")
  clineMLE <- c(clineMLE ,list(null.model=temp.getNULL(obsData)))
  clineAIC<- sapply( clineMLE, hzar.AIC.hzar.cline)
  print(do.call(data.frame,c(as.list(clineAIC),row.names=locus)))
  save(obsData,clineModels,clineInit,clineMLE,clineAIC,file=file.path("cline-models",paste(locus,"MLE.RData",sep="")))
  clineAIC
}

