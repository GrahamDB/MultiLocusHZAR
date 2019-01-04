## Jacana HZAR

# Transform Lat Long into Distance
Locations <- read.csv("~/Dropbox/Jacanas/Jacana_HZAR/jacanaLatLong.csv")
lat<-Locations$Lat
lon<-Locations$Long
earthDist <- function (lon1, lat1, lon2, lat2){
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- lon1 * rad
  b1 <- lat2 * rad
  b2 <- lon2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}

earthDist(lon[1], lat[1], lon, lat)
write.table(earthDist(lon[1], lat[1], lon, lat), file="~/Dropbox/Jacanas/Jacana_HZAR/jacanaDistance.csv", sep=",")

## Load the package
# install.packages("hzar")
library(hzar)

## Write out internal example trait data 
jacanaMorphological <- read.csv("~/Dropbox/Jacanas/Jacana_HZAR/jacanaMorphological.csv")
data(jacanaMorphological) 
  # Warning message:In data(jacanaMorphological) : data set ‘jacanaMorphological’ not found
print(summary(jacanaMorphological))

write.table(jacanaMorphological, # The data we just loaded
            file="jacExTrait.txt", # The file to overwrite
            col.names=TRUE, # The columns are named
            row.names=FALSE, # The rows are not named
            sep="\t", # The file will be tab-delimited
            quote=TRUE) # Use quotes as needed.

## As we no longer need the in-memory copy, drop the local reference
jacanaMorphological <- NULL

## Write out internal example site data
jacanaLocations <- read.csv("~/Dropbox/Jacanas/Jacana_HZAR/jacanaLocations.csv")
data(jacanaLocations) # Warning message: In data(jacanaLocations) : data set ‘jacanaLocations’ not found
print(jacanaLocations)

write.table(jacanaLocations, # The data we just loaded
            file="jacExSite.txt", # The file to overwrite
            col.names=TRUE, # The columns are named
            row.names=FALSE, # The rows are not named
            sep="\t", # The file will be tab-delimited
            quote=TRUE) # Use quotes as needed.

## As we no longer need the in-memory copy, drop the local reference
jacanaLocations<- NULL

## Save all plots in a series of png files
png(width=900, height=900, res=200, family="Arial", filename="jExQTPlot%03d.png",pointsize=8)

## A typical chain length
chainLength=1e5; 

## Make each model run off a separate seed
# install.packages("doMC")
mainSeed=list(A=c(978,544,99,596,528,124),B=c(544,99,596,528,124,978),C=c(99,596,528,124,978,544))

if(require(doMC)){  ## If you have doMC, use foreach in parallel mode to speed up computation.
    registerDoMC()
  } else {  ## Use foreach in sequential mode
    registerDoSEQ();
  }

## Trait Analysis
  ## Load example Quantitative trait data from the package
  ## Note that this is the individual data, with labels
  ## identifying the locality of each sample.
jacanaMorphological <- read.table("jacExTrait.txt",header=TRUE)
## Print a summary of the trait data
print(summary(jacanaMorphological))

## Load example locality data, matching each localility to a site ID and a transect distance.
jacanaLocations <- read.table("jacExSite.txt",header=TRUE)

## Print the locality information
print(jacanaLocations)

## Blank out space in memory to hold morphological analysis
if(length(apropos("^jac$",ignore.case=FALSE)) == 0 || !is.list(jac) ) jac <- list()
## We are doing just the one quantitative trait, but it is good to stay organized.
jac$Mass <- list();
## Space to hold the observed data
jac$Mass$obs <- list();
## Space to hold the models to fit
jac$Mass$models <- list();
## Space to hold the compiled fit requests
jac$Mass$fitRs <- list();
## Space to hold the output data chains
jac$Mass$runs <- list();
## Space to hold the analysed data
jac$Mass$analysis <- list();

## Beard Length Trait from Brumfield et al 2001 - change to Jacana Mass
jac$Mass$obs <- hzar.doNormalData1DRaw(hzar.mapSiteDist(jacanaLocations$LocalityID,
  jacanaLocations$distance), jacanaMorphological$Locality, jacanaMorphological$Mass)

## Look at a graph of the observed data
hzar.plot.obsData(jac$Mass$obs)
  ## This didn't work
