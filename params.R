

library(hzar)
library(subplex)

if(require(doMC)){
  ## If you have doMC, use foreach in parallel mode
  ## to speed up computation.
  registerDoMC(3)
} else {
  ## Use foreach in sequential mode
  registerDoSEQ();
}


## loci<-read.table("significate_loci_names.txt",stringsAsFactors=FALSE)[[1]]
#loci<-read.table("jacanaMolecular.csv",stringsAsFactors=FALSE)[[1]]
loci<-("mtDNA")
traits<-c("GenPCF","MassF", "AvgSpurF")

all.clines <- c(loci,traits)


## cline range (assuming observations between 0 and 1200 km)

far.left  = -200
far.right = 900
