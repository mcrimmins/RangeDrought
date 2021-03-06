# assess climate/drought index relationships with rpms/smNDVI
# by NRCS LRUs, CRAs
# MAC 04/17/2020
# ordinal regression to predict drought categories
# with higher order models and model validation

library(raster)
library(rgdal)
library(reshape2)
library(ggplot2)
library(SPEI)
#library(leaps)
#library(caret)
#library(tidyverse)
library(gtools)
library(MASS)
library(rms)
library(pracma)

# set rasteroptions
rasterOptions(progress = 'text')

# custom formulas
# replace empty with NA in list
nullToNA <- function(x) {
  x[sapply(x, is.null)] <- NA
  return(x)
}

#####


# load and process grids ----
  # production
  rpmsSW<-stack("/scratch/crimmins/USDA_NPP/v2019/AZNM_RPMS_8419_WGS84.grd")
  # monthly gridmet, orig res
  precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_precip_1979_2019.grd")
  pet<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_pet_1979_2019.grd")

  # # try out nClimGrid
  # nClimPrecip<-stack("/scratch/crimmins/climgrid/processed/WESTmonthly.prcp.conus.pnt_1895_2019.grd")
  # nClimPrecip<-nClimPrecip[[((nlayers(nClimPrecip)-nlayers(precip))+1):nlayers(nClimPrecip)]]
  # nprecip<-crop(nClimPrecip, extent(precip))  

  # try out smoothed NDVI
    # smNDVI<-stack("/scratch/crimmins/vhi/processed/ANNUAL_MAX_smNDVI_complete_1982-2019_SWUS.grd") 
    # smNDVI<-resample(smNDVI, precip[[1]])
     #smNDVI<-smNDVI[[3:38]]
    # rpmsSW<-smNDVI      

  # create RPMS NA mask
  maskRPMS<-calc(rpmsSW, sum)
  maskRPMS[maskRPMS>=0] <- 1
  maskRPMSgm <- resample(maskRPMS,precip[[1]],method='ngb')
  #precip<-crop(precip, maskRPMS)
#  precip<-raster::mask(precip, maskRPMSgm)  
#  pet<-raster::mask(pet, maskRPMSgm)
  # nprecip mask only
  #maskRPMS<-calc(rpmsSW, sum)
  #  maskRPMS[maskRPMS>=0] <- 1
  #maskRPMSn <- resample(maskRPMS,nprecip[[1]],method='ngb')
  #nprecip<-raster::mask(nprecip, maskRPMSn)

# County as LMUs ----
#  us<-getData('GADM', country='USA', level=2)
#  counties<-subset(us, NAME_1=="Arizona" | NAME_1=="New Mexico")
  #county<-subset(counties, NAME_2=="Santa Cruz")
  
# ---- load NRCS CRAs for SW
AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
# combine
SW_CRA<- rbind(AZ_CRA,NM_CRA)

# AZ LRU from Emilio
AZ_LRU<-readOGR(dsn = "./shapes/AZ_LRU", layer = "mlra_a_az")
AZ_LRU<-spTransform(AZ_LRU, crs(AZ_CRA))

# poly grids
#halfDeg<-readOGR(dsn = "./shapes", layer = "polyGrid_eighthDeg")

# switch to land management unit
LMU<-AZ_LRU
#LMU<-counties
# add id code to data
LMU$ID<-seq(1,nrow(LMU),1)

# put results in dataframe
resultsFrame = data.frame(ID=rep(0, nrow(LMU)),
                          rpmsCoverage=rep(0,nrow(LMU)),
                          r2_spi=rep(0,nrow(LMU)),
                          form_spi=rep(0,nrow(LMU)),
                          spi_pVal=rep(0,nrow(LMU)),
                          r2_spei=rep(0,nrow(LMU)),
                          form_spei=rep(0,nrow(LMU)),
                          spei_pVal=rep(0,nrow(LMU)))
                          
# save monthly SPI and SPEI in lists
spiData<-list()
speiData<-list()
rpmsData<-list()
rpmsData_detrend<-list()
bestSPImodels<-list()
bestSPEImodels<-list()

# switch to forecast type, needs to match dataset
quantCuts<-c(0,0.33,0.66,1); quantNames<-c("Below","Normal","Above") # for 3 cat
#quantCuts<-c(0,0.33,1); quantNames<-c("Below","Normal-Above") # for 2 cat
#quantCuts<-c(0,0.03,0.06,0.11,0.21,0.30,1); quantNames<- c("D4","D3","D2","D1","D0","No Drought") # for 3 cat

# earliest year for RS data?
yearRS<-1984


# Start the clock!
ptm <- proc.time()

# loop through shapefile to get regression results with SPI and SPEI
for(i in 1:nrow(LMU)){
  poly<-LMU[i,]
  # summary stat ts for LRU mean or median
  beginCluster(n=6)
    rpmsCoverage<-raster::extract(maskRPMS, poly, df=TRUE)
    rpmsTS<-raster::extract(rpmsSW, poly, fun=median, df=TRUE, na.rm=TRUE)
    precipTS<-raster::extract(precip, poly, fun=median, df=TRUE, na.rm=TRUE)
    #nprecipTS<-extract(nprecip, LRU, fun=median, df=TRUE, na.rm=TRUE)
    petTS<-raster::extract(pet, poly, fun=median, df=TRUE, na.rm=TRUE)
  endCluster()
  
  # percent coverage of data in polygon - for rpms
  rpmsCoverage<-(sum(rpmsCoverage$layer, na.rm = TRUE)/nrow(rpmsCoverage))*100
  
  # format time series data frames
  rpmsTS<-as.data.frame(t(rpmsTS))
  rpmsTS = as.data.frame(rpmsTS[-1,]) # ONLY IF MULTIPLE POLYGONS
  rpmsTS$date<-seq(as.Date(paste0(yearRS,"-12-01"), "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="years") # 1982 or 1984
    # weighted mean of polygons by area
    #wgts<-poly@data[["AREA"]]/sum(poly@data[["AREA"]])
    #rpmsTS$wgtMean<-apply(as.matrix(rpmsTS[,1:(ncol(rpmsTS)-1)]), 1, function(x) weighted.mean(x, wgts))
  colnames(rpmsTS)[1]<-"wgtMean"
  
  # format climate dataframes
  precipTS<-as.data.frame(t(precipTS))
  precipTS = as.data.frame(precipTS[-1,])
  precipTS$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
    # weighted mean of polygons by area, comment out for CRAs
    #wgts<-poly@data[["AREA"]]/sum(poly@data[["AREA"]])
    #precipTS$wgtMean<-apply(as.matrix(precipTS[,1:(ncol(precipTS)-1)]), 1, function(x) weighted.mean(x, wgts))
    # add in for CRAs/comment out for LRU
    colnames(precipTS)[1]<-"wgtMean"
  # pet  
  petTS<-as.data.frame(t(petTS))
  petTS = as.data.frame(petTS[-1,])
  petTS$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
    # weighted mean of polygons by area
    #wgts<-poly@data[["AREA"]]/sum(poly@data[["AREA"]])
    #petTS$wgtMean<-apply(as.matrix(petTS[,1:(ncol(petTS)-1)]), 1, function(x) weighted.mean(x, wgts))
    colnames(petTS)[1]<-"wgtMean"
  # nClimGrid ONLY
  # nprecipTS<-as.data.frame(t(nprecipTS))
  # nprecipTS = as.data.frame(nprecipTS[-1,])
  # nprecipTS$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
  #   # weighted mean of polygons by area
  #   wgts<-poly@data[["AREA"]]/sum(poly@data[["AREA"]])
  #   nprecipTS$wgtMean<-apply(as.matrix(nprecipTS[,1:(ncol(nprecipTS)-1)]), 1, function(x) weighted.mean(x, wgts))
    
  # weighted mean climate dataframe - try different windows
    moClimData<-cbind.data.frame(precipTS$date, precipTS$wgtMean, petTS$wgtMean)
    moClimData$precip3mo<-movingFun(moClimData$`precipTS$wgtMean`, 3,type = "to", sum)
    moClimData$pet3mo<-movingFun(moClimData$`petTS$wgtMean`, 3,type = "to", sum)
    moClimData$SPI3<-spi(moClimData$`precipTS$wgtMean`, 3)$fitted  
    moClimData$SPEI3<-spei(moClimData$`precipTS$wgtMean`- moClimData$`petTS$wgtMean`, 3)$fitted
    moClimData$month<-as.numeric(format(moClimData$`precipTS$date`, format="%m"))
    moClimData$year<-as.numeric(format(moClimData$`precipTS$date`, format="%Y"))
    
  # grab SPI,SPEI,rpms data
    spiData[[i]]<-moClimData$SPI3
    speiData[[i]]<-moClimData$SPEI3
    rpmsData[[i]]<-rpmsTS$wgtMean
  # detrend remote sensing data
    rpmsData_detrend[[i]]<-detrend(rpmsTS$wgtMean, 'linear')
    rpmsTS$wgtMean<-detrend(rpmsTS$wgtMean, 'linear')
    
    # all subsets to build models with all seasons
    # SPI = 6, 
    allSeas<-moClimData[,c(8,9,6)]
    allSeas<-dcast(allSeas, year~month)
    colnames(allSeas)<-c("year",paste0(month.abb[11],"_",month.abb[1]),paste0(month.abb[12],"_",month.abb[2]),
                         paste0(month.abb[seq(1,12-2,1)],"_",month.abb[seq(1+2,12,1)]))
      # all subsets regressions 
      trimSeas<-subset(allSeas, year>=yearRS) # change to 1982 for smNDVI or 1984 for RPMS
      trimSeas$rpms<-rpmsTS$wgtMean
      trimSeas<-trimSeas[,-1]
      # create factors
      #trimSeas$rpms<- quantcut(trimSeas$rpms, c(0,0.03,0.06,0.11,0.21,0.30,1)) # USDM percentiles
       trimSeas$rpms<- quantcut(trimSeas$rpms, quantCuts)
      # rename levels
      #levels(trimSeas$rpms) <- c("D4","D3","D2","D1","D0","No Drought")
       levels(trimSeas$rpms) <- quantNames
      
      # find best subset of 1 and 2 model combinations
      regMat <- expand.grid(c(TRUE,FALSE), c(TRUE,FALSE),
                            c(TRUE,FALSE), c(TRUE,FALSE),
                            c(TRUE,FALSE), c(TRUE,FALSE),
                            c(TRUE,FALSE), c(TRUE,FALSE),
                            c(TRUE,FALSE), c(TRUE,FALSE),
                            c(TRUE,FALSE), c(TRUE,FALSE))
                            
      # find one and two factor models
      regMat <- regMat[which(rowSums(regMat == "TRUE")<=3),]
      regMat <- regMat[-(dim(regMat)[1]),]
      names(regMat) <- colnames(trimSeas[1:ncol(trimSeas)-1])
      regressors <- colnames(trimSeas[1:ncol(trimSeas)-1])
      # all model list
      allModelsList <- apply(regMat, 1, function(x) as.formula(
        paste(c("rpms ~", regressors[x]),
              collapse=" + ")) )
      # execute all models with lrm
      allM <- lapply(allModelsList,
                     function(x) lrm(x, data=trimSeas, x=TRUE, y=TRUE, tol = 1e-6)) # adjusted tol from 1e-7
      # extract best model
      stats<-lapply(X = allM, FUN = `[[`, "stats")
        r2<-lapply(X=stats, FUN = `[[`, 10)
      chiP<-lapply(X = allM, FUN = `[[`, "stats")
        chiP<-lapply(X=chiP, FUN = `[[`, 5)  
      model<-lapply(X = allM, FUN = `[[`, "terms")
        model<-lapply(X = model, FUN = `[[`, 3)
      dev<-lapply(X = allM, FUN = `[[`, c("deviance"))
        dev<-lapply(X = dev, FUN = `[[`, 2)
      r2<-as.data.frame(unlist(nullToNA(r2)))
        colnames(r2)<-"r2value"
      r2$model<-unlist(nullToNA(model))
       r2$dev<-unlist(nullToNA(dev))
       r2$chiP<-unlist(nullToNA(chiP))
      # save model in list
      bestmodels<-allM[[which(r2$r2value==max(r2$r2value, na.rm = TRUE))]] 
       # stepwise reduced model
        tempModel<-fastbw(bestmodels, rule="p", sls = 0.05) # or AIC rule="p", sls = 0.05
        if(length(tempModel$names.kept)==0){
          tempModel<-bestmodels
        }else{
        tempModel<-lrm(as.formula(paste(c("rpms ~", tempModel$names.kept ),
                collapse=" + ")), data=trimSeas, x=TRUE, y=TRUE, tol = 1e-6)
        }
      # vars for table 
      r2_spi<-tempModel$stats[["R2"]]
      form_spi<-tempModel$terms[[3]]
      spi_pVal<-tempModel$stats[["P"]]
      # save model in list
      bestSPImodels[[i]]<-tempModel
      # clean up
      rm(allM, allSeas, trimSeas, regMat)
        
      # all subsets to build models with all seasons
      # SPEI = 7, 
      allSeas<-moClimData[,c(8,9,7)]
      allSeas<-dcast(allSeas, year~month)
      colnames(allSeas)<-c("year",paste0(month.abb[11],"_",month.abb[1]),paste0(month.abb[12],"_",month.abb[2]),
                           paste0(month.abb[seq(1,12-2,1)],"_",month.abb[seq(1+2,12,1)]))
      # all subsets regressions 
      trimSeas<-subset(allSeas, year>=yearRS) # change to 1982 for smNDVI or 1984 for RPMS
      trimSeas$rpms<-rpmsTS$wgtMean
      trimSeas<-trimSeas[,-1]
      # create factors
      #trimSeas$rpms<- quantcut(trimSeas$rpms, c(0,0.03,0.06,0.11,0.21,0.30,1)) # USDM percentiles
       trimSeas$rpms<- quantcut(trimSeas$rpms, quantCuts)
      # rename levels
      #levels(trimSeas$rpms) <- c("D4","D3","D2","D1","D0","No Drought")
       levels(trimSeas$rpms) <- quantNames
      
      # find best subset of 1 and 2 model combinations
      regMat <- expand.grid(c(TRUE,FALSE), c(TRUE,FALSE),
                            c(TRUE,FALSE), c(TRUE,FALSE),
                            c(TRUE,FALSE), c(TRUE,FALSE),
                            c(TRUE,FALSE), c(TRUE,FALSE),
                            c(TRUE,FALSE), c(TRUE,FALSE),
                            c(TRUE,FALSE), c(TRUE,FALSE))
                        
      # find one and two factor models
      regMat <- regMat[which(rowSums(regMat == "TRUE")<=3),]
      regMat <- regMat[-(dim(regMat)[1]),]
      names(regMat) <- colnames(trimSeas[1:ncol(trimSeas)-1])
      regressors <- colnames(trimSeas[1:ncol(trimSeas)-1])
      # all model list
      allModelsList <- apply(regMat, 1, function(x) as.formula(
        paste(c("rpms ~", regressors[x]),
              collapse=" + ")) )
      # execute all models with lrm
      allM <- lapply(allModelsList,
                     function(x) lrm(x, data=trimSeas, x=TRUE, y=TRUE,tol = 1e-6))
      # extract best model
      stats<-lapply(X = allM, FUN = `[[`, "stats")
      r2<-lapply(X=stats, FUN = `[[`, 10)
      chiP<-lapply(X = allM, FUN = `[[`, "stats")
      chiP<-lapply(X=chiP, FUN = `[[`, 5)  
      model<-lapply(X = allM, FUN = `[[`, "terms")
      model<-lapply(X = model, FUN = `[[`, 3)
      dev<-lapply(X = allM, FUN = `[[`, c("deviance"))
      dev<-lapply(X = dev, FUN = `[[`, 2)
      r2<-as.data.frame(unlist(nullToNA(r2)))
      colnames(r2)<-"r2value"
      r2$model<-unlist(nullToNA(model))
      r2$dev<-unlist(nullToNA(dev))
      r2$chiP<-unlist(nullToNA(chiP))
      # save model in list
      bestmodels<-allM[[which(r2$r2value==max(r2$r2value, na.rm = TRUE))]] 
      # stepwise reduced model
      tempModel<-fastbw(bestmodels, rule="p", sls = 0.05) # or AIC rule="p", sls = 0.05 / rule="aic"
      if(length(tempModel$names.kept)==0){
        tempModel<-bestmodels
      }else{
        tempModel<-lrm(as.formula(paste(c("rpms ~", tempModel$names.kept ),
                                        collapse=" + ")), data=trimSeas, x=TRUE, y=TRUE, tol = 1e-6)
      }
      # vars for table 
      r2_spei<-tempModel$stats[["R2"]]
      form_spei<-tempModel$terms[[3]]
      spei_pVal<-tempModel$stats[["P"]]
      # save model in list
      bestSPEImodels[[i]]<-tempModel
      # clean up
      rm(allM, allSeas, trimSeas, regMat)
      
    # put results in dataframe
    resultsFrame[i, ] = c(i,rpmsCoverage,r2_spi, paste0(form_spi, collapse=', ' ), spi_pVal,
                          r2_spei, paste0(form_spei, collapse=', ' ), spei_pVal)
    
    print(i)
    
}

# Stop the clock
proc.time() - ptm

# save a snapshot before breaking stuff
save.image("~/RProjects/RangeDrought/tempWorkspace.RData")
load("~/RProjects/RangeDrought/tempWorkspace.RData")

# create data frames from saved lists ##### FIX
spiDataFrame =as.data.frame(do.call(cbind, spiData))
  colnames(spiDataFrame)<-paste0("LRU-",seq(1,ncol(spiDataFrame),1))
    spiDataFrame$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
speiDataFrame =as.data.frame(do.call(cbind, speiData))
  colnames(speiDataFrame)<-paste0("LRU-",seq(1,ncol(speiDataFrame),1))
    speiDataFrame$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
rpmsDataFrame = as.data.frame(do.call(cbind, rpmsData))
  colnames(rpmsDataFrame)<-paste0("LRU-",seq(1,ncol(rpmsDataFrame),1))
    rpmsDataFrame$date<-seq(as.Date(paste0(yearRS,"-12-01"), "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="years")
  
# merge with spatial dataframe
resultsFrame[c(1,2,3,5,6,8)] <- sapply(resultsFrame[c(1,2,3,5,6,8)],as.numeric)
resultsFrame[c(4,7)] <- lapply(resultsFrame[c(4,7)],factor)
LMU@data<-merge(LMU@data, resultsFrame, by.x="ID", by.y="ID")
# save into file for later use
save(LMU,rpmsData_detrend,spiDataFrame,speiDataFrame,rpmsDataFrame,bestSPEImodels,
     bestSPImodels,file="./results/AZ_LRU_smNDVI_detrend_stepwise_pval_gridmet_SPI3_SPEI3_3cats_results.Rdata")

# plot of results by polygon
load(file="./results/AZ_LMU_rpms_gridmet_SPI3_SPEI3_USDMcats_results.Rdata")
library(rasterVis)
library(RColorBrewer)

# plot individual polygons
plot(smNDVI[[1]]); 
plot(LMU, add=TRUE); plot(LMU[47,], add=TRUE, border='red')

# look at distributions of values
plot(density(LMU@data$r2_spi))
lines(density(LMU@data$r2_spei))
# Shapiro-Wilk normality test for the differences
  shapiro.test(LMU@data$r2_2spi-LMU@data$r2_2spei) # => p-value >0.05, not different from normal
  t.test(LMU@data$r2_2spi, LMU@data$r2_2spei, paired = TRUE, alternative = "two.sided")
  summary(LMU@data$r2_2spei)
  summary(LMU@data$r2_2spi)
  
# plot r2 values
  
pal <- RColorBrewer::brewer.pal(n=9, 'OrRd')
spplot(LMU, c("r2_spi","r2_spei"), at=seq(0.1, 1, 0.1),
       main=list(label="R2- Best Regression to predict detrended MAX-smNDVI based on 3mo-SPI/SPEI",cex=1),
       layout=c(2, 1),scales=list(draw = TRUE),  col.regions = pal)+
  layer(sp.polygons(AZ_LRU, line="grey"))

# plot formulas
library(RColorBrewer)
# random colors -- https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#
spplot(LMU, c("form_spi"), col.regions=sample(col_vector, length(levels(LMU$form_spi))),
       main=list(label="Best Regression to predict RPMS based on 3mo-SPI",cex=1),scales=list(draw = TRUE))
spplot(LMU, c("form2_spi"), col.regions=sample(col_vector, length(levels(LMU$form2_spi))),
       main=list(label="2nd Best Regression to predict RPMS based on 3mo-SPI",cex=1),scales=list(draw = TRUE))

spplot(LMU, c("form_spei"), col.regions=sample(col_vector, length(levels(LMU$form_spei))),
       main=list(label="Best Regression to predict RPMS based on 3mo-SPEI",cex=1),scales=list(draw = TRUE))
spplot(LMU, c("form2_spei"), col.regions=sample(col_vector, length(levels(LMU$form2_spei))),
       main=list(label="2nd Best Regression to predict RPMS based on 3mo-SPEI",cex=1),scales=list(draw = TRUE))
