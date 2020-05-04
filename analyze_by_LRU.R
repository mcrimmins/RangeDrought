# assess climate/drought index relationships with rpms/smNDVI
# by NRCS LRUs, CRAs
# MAC 04/17/2020

library(raster)
library(rgdal)
library(reshape2)
library(ggplot2)
library(SPEI)
library(leaps)
library(caret)
library(tidyverse)

# set rasteroptions
rasterOptions(progress = 'text')

# custom formulas
get_model_formula <- function(id, object, outcome){
  # get models data
  models <- summary(object)$which[id,-1]
  # Get outcome variable
  #form <- as.formula(object$call[[2]])
  #outcome <- all.vars(form)[1]
  # Get model predictors
  predictors <- names(which(models == TRUE))
  predictors <- paste(predictors, collapse = "+")
  # Build model formula
  as.formula(paste0(outcome, "~", predictors))
}
get_cv_error <- function(model.formula, data){
  set.seed(1)
  train.control <- trainControl(method = "cv", number = 5)
  cv <- train(model.formula, data = data, method = "lm",
              trControl = train.control)
  cv$results$RMSE
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
   # smNDVI<-smNDVI[[3:38]]
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
  us<-getData('GADM', country='USA', level=2)
  counties<-subset(us, NAME_1=="Arizona" | NAME_1=="New Mexico")
  #county<-subset(counties, NAME_2=="Santa Cruz")
  
# ---- load NRCS CRAs for SW
AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
# combine
SW_CRA<- rbind(AZ_CRA,NM_CRA)

# AZ LRU from Emilio
AZ_LRU<-readOGR(dsn = "./shapes/AZ_LRU", layer = "mlra_a_az")
AZ_LRU<-spTransform(AZ_LRU, crs(AZ_CRA))

# switch to land management unit
#LMU<-SW_CRA
LMU<-counties
# add id code to data
LMU$ID<-seq(1,nrow(LMU),1)

# loop through by polygon
resultsFrame = data.frame(ID=rep(0, nrow(LMU)),
                          rpmsCoverage=rep(0,nrow(LMU)),
                          form1_spi=rep(0,nrow(LMU)),
                          r2_1spi=rep(0,nrow(LMU)),
                          form2_spi=rep(0,nrow(LMU)),
                          r2_2spi=rep(0,nrow(LMU)),
                          form1_spei=rep(0,nrow(LMU)),
                          r2_1spei=rep(0,nrow(LMU)),
                          form2_spei=rep(0,nrow(LMU)),
                          r2_2spei=rep(0,nrow(LMU)),
                          spi_var1=rep(0,nrow(LMU)),
                          spi_var2=rep(0,nrow(LMU)),
                          spi_fPval=rep(0,nrow(LMU)),
                          spei_var1=rep(0,nrow(LMU)),
                          spei_var2=rep(0,nrow(LMU)),
                          spei_fPval=rep(0,nrow(LMU)))

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
  rpmsTS$date<-seq(as.Date("1984-12-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="years")
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
    
    # all subsets to build models with all seasons
    # SPI = 6, 
    allSeas<-moClimData[,c(8,9,6)]
    allSeas<-dcast(allSeas, year~month)
    colnames(allSeas)<-c("year",paste0(month.abb[11],"_",month.abb[1]),paste0(month.abb[12],"_",month.abb[2]),
                         paste0(month.abb[seq(1,12-2,1)],"_",month.abb[seq(1+2,12,1)]))
      # all subsets regressions 
      trimSeas<-subset(allSeas, year>=1984)
      trimSeas$rpms<-rpmsTS$wgtMean
    # find best subset
    models <- regsubsets(rpms~., data = trimSeas, nvmax = 2)
    res.sum <- summary(models)
      r2_1spi<-res.sum$adjr2[1]
      r2_2spi<-res.sum$adjr2[2]
      form1_spi<-as.character(Reduce(paste, deparse(get_model_formula(1, models, "rpms"))))
      form2_spi<-as.character(Reduce(paste, deparse(get_model_formula(2, models, "rpms"))))
      # Compute cross-validation error
      model.ids <- 1:2
      cv.errors <-  map(model.ids, get_model_formula, models, "rpms") %>%
        map(get_cv_error, data = trimSeas) %>%
        unlist()
      tempModel<-lm(get_model_formula(which.min(cv.errors), models, "rpms"), data=trimSeas)
        spi_var1<-coef(tempModel)[2]
        spi_var2<-coef(tempModel)[3]
          fstat <- summary(tempModel)$fstatistic
        spi_fPval<- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
    rm(models, allSeas, trimSeas, tempModel) # clean out vars
    
    # all subsets to build models with all seasons
    # SPEI = 7
    allSeas<-moClimData[,c(8,9,7)]
    allSeas<-dcast(allSeas, year~month)
    colnames(allSeas)<-c("year",paste0(month.abb[11],"_",month.abb[1]),paste0(month.abb[12],"_",month.abb[2]),
                         paste0(month.abb[seq(1,12-2,1)],"_",month.abb[seq(1+2,12,1)]))
      # all subsets regressions 
      trimSeas<-subset(allSeas, year>=1984)
      trimSeas$rpms<-rpmsTS$wgtMean
    # find best subset
    models <- regsubsets(rpms~., data = trimSeas, nvmax = 2)
    res.sum <- summary(models)
      r2_1spei<-res.sum$adjr2[1]
      r2_2spei<-res.sum$adjr2[2]
      form1_spei<-as.character(Reduce(paste, deparse(get_model_formula(1, models, "rpms"))))
      form2_spei<-as.character(Reduce(paste, deparse(get_model_formula(2, models, "rpms"))))
      # Compute cross-validation error
      model.ids <- 1:2
      cv.errors <-  map(model.ids, get_model_formula, models, "rpms") %>%
        map(get_cv_error, data = trimSeas) %>%
        unlist()
      tempModel<-lm(get_model_formula(which.min(cv.errors), models, "rpms"), data=trimSeas)
      spei_var1<-coef(tempModel)[2]
      spei_var2<-coef(tempModel)[3]
      fstat <- summary(tempModel)$fstatistic
      spei_fPval<- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
    rm(models, allSeas, trimSeas, tempModel) # clean out vars
      
    # ADD IN CROSS VALIDATION FOR MODELS
    
    # put results in dataframe
    resultsFrame[i, ] = c(i,rpmsCoverage,form1_spi,r2_1spi,form2_spi,r2_2spi,form1_spei,r2_1spei,form2_spei,r2_2spei,
                          spi_var1,spi_var2,spi_fPval,spei_var1,spei_var2,spei_fPval)
    
    print(i)
    
}

# merge with spatial dataframe
resultsFrame[c(1,2,4,6,8,10)] <- sapply(resultsFrame[c(1,2,4,6,8,10)],as.numeric)
resultsFrame[c(3,5,7,9)] <- lapply(resultsFrame[c(3,5,7,9)],factor)
LMU@data<-merge(LMU@data, resultsFrame, by.x="ID", by.y="ID")
# save into file for later use
save(LMU, file="./results/AZNM_counties_rpms_gridmet_SPI_SPEI_results.Rdata")

# plot of results by polygon
load(file="./results/SW_CRA_rpms_gridmet_SPI_SPEI_results.Rdata")
library(rasterVis)
library(RColorBrewer)

# plot individual polygons
plot(LMU);plot(LMU[55,], add=TRUE, border='red')

# look at distributions of values
plot(density(LMU@data$r2_1spi))
lines(density(LMU@data$r2_1spei))
# Shapiro-Wilk normality test for the differences
  shapiro.test(LMU@data$r2_2spi-LMU@data$r2_2spei) # => p-value >0.05, not different from normal
  t.test(LMU@data$r2_2spi, LMU@data$r2_2spei, paired = TRUE, alternative = "two.sided")
  summary(LMU@data$r2_2spei)
  summary(LMU@data$r2_2spi)
  
# plot r2 values
pal <- RColorBrewer::brewer.pal(n=9, 'OrRd')
spplot(LMU, c("r2_1spi","r2_2spi","r2_1spei","r2_2spei"), at=seq(0.1, 1, 0.1),
       main=list(label="R2- Best Regression to predict RPMS based on 3mo-SPI/SPEI",cex=1),
       layout=c(2, 2),scales=list(draw = TRUE),  col.regions = pal)

# plot formulas
# random colors -- https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#
spplot(LMU, c("form1_spi"), col.regions=sample(col_vector, length(levels(LMU$form1_spi))),
       main=list(label="Best Regression to predict RPMS based on 3mo-SPI",cex=1),scales=list(draw = TRUE))
spplot(LMU, c("form2_spi"), col.regions=sample(col_vector, length(levels(LMU$form2_spi))),
       main=list(label="2nd Best Regression to predict RPMS based on 3mo-SPI",cex=1),scales=list(draw = TRUE))

spplot(LMU, c("form1_spei"), col.regions=sample(col_vector, length(levels(LMU$form1_spei))),
       main=list(label="Best Regression to predict RPMS based on 3mo-SPEI",cex=1),scales=list(draw = TRUE))
spplot(LMU, c("form2_spei"), col.regions=sample(col_vector, length(levels(LMU$form2_spei))),
       main=list(label="2nd Best Regression to predict RPMS based on 3mo-SPEI",cex=1),scales=list(draw = TRUE))
