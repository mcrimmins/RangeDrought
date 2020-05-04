# extract RPMS amd climate time series by land management unit
# MAC 04/09/2020

library(raster)
library(rgdal)
library(reshape2)
library(ggplot2)
library(SPEI)

# set rasteroptions
rasterOptions(progress = 'text')

# ---- load NRCS CRAs for SW
AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")

# AZ LRU from Emilio
AZ_LRU<-readOGR(dsn = "./shapes/AZ_LRU", layer = "mlra_a_az")
AZ_LRU<-spTransform(AZ_LRU, crs(AZ_CRA))

# load reprojected RPMS data and gridmet
#pet<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_monthly_sum_pet_1979_2019.grd")
#precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_monthly_sum_precip_1979_2019.grd")
rpmsSW<-stack("/scratch/crimmins/USDA_NPP/v2019/AZNM_RPMS_8419_WGS84.grd")
# monthly gridmet, orig res
precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_precip_1979_2019.grd")
pet<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_pet_1979_2019.grd")

# try out nClimGrid
nClimPrecip<-stack("/scratch/crimmins/climgrid/processed/WESTmonthly.prcp.conus.pnt_1895_2019.grd")
  nClimPrecip<-nClimPrecip[[((nlayers(nClimPrecip)-nlayers(precip))+1):nlayers(nClimPrecip)]]
  nprecip<-crop(nClimPrecip, extent(precip))  

# try out smoothed NDVI
smNDVI<-stack("/scratch/crimmins/vhi/processed/ANNUAL_MAX_smNDVI_complete_1982-2019_SWUS.grd") 
  smNDVI<-resample(smNDVI, precip[[1]])
  smNDVI<-smNDVI[[3:38]]
rpmsSW<-smNDVI      
  
# create RPMS NA mask
maskRPMS<-calc(rpmsSW, sum)
  maskRPMS[maskRPMS>=0] <- 1
  maskRPMSgm <- resample(maskRPMS,precip[[1]],method='ngb')
  #precip<-crop(precip, maskRPMS)
  precip<-raster::mask(precip, maskRPMSgm)  
  pet<-raster::mask(pet, maskRPMSgm)
# nprecip mask only
#maskRPMS<-calc(rpmsSW, sum)
#  maskRPMS[maskRPMS>=0] <- 1
  maskRPMSn <- resample(maskRPMS,nprecip[[1]],method='ngb')
  nprecip<-raster::mask(nprecip, maskRPMSn)
# subset LRU
LRU<-subset(AZ_LRU, LRU=="41-3AZ")
  #LRU<-LRU[5,]
# plot individual units
# plot(LRU[1,])
  plot(AZ_LRU)
  plot(LRU, add=TRUE, border='red')
  
  
# # rpms-climate pixel wise results
#   pixelResults<-stack("./results/rpms_AZNM_precip3mo_pearson_1984_2019.grd")
#   pixelResults<-extract(pixelResults, LRU, df=TRUE)
#     boxplot(as.matrix(pixelResults[,2:13]), use.cols=TRUE)
  
# extract all pixels
# beginCluster(n=6)
#   ts<-extract(rpmsSW, LRU, df=TRUE)
# endCluster()
# boxplot(test[,2:37], use.cols = TRUE)

# summary stat ts for LRU mean or median
beginCluster(n=6)
  rpmsTS<-extract(rpmsSW, LRU, fun=median, df=TRUE, na.rm=TRUE)
  precipTS<-extract(precip, LRU, fun=median, df=TRUE, na.rm=TRUE)
  nprecipTS<-extract(nprecip, LRU, fun=median, df=TRUE, na.rm=TRUE)
  petTS<-extract(pet, LRU, fun=median, df=TRUE, na.rm=TRUE)
endCluster()

# format rpms dataframe
rpmsTS<-as.data.frame(t(rpmsTS))
  rpmsTS = as.data.frame(rpmsTS[-1,]) # ONLY IF MULTIPLE POLYGONS
  rpmsTS$date<-seq(as.Date("1984-12-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="years")
# corr of polygons if multiple  
  cor(rpmsTS[,1:(ncol(rpmsTS)-1)], method = "pearson", use = "complete.obs")
# weighted mean of polygons by area
  wgts<-LRU@data[["AREA"]]/sum(LRU@data[["AREA"]])
  rpmsTS$wgtMean<-apply(as.matrix(rpmsTS[,1:(ncol(rpmsTS)-1)]), 1, function(x) weighted.mean(x, wgts))
  # test plot
  #test<- melt(rpmsTS, id.vars = "date")
  #ggplot(test, aes(date,value,color=variable))+
  #  geom_line()+
  #  ggtitle("Median RPMS by Year LRU 41-3AZ ")
    
# format climate dataframes
  precipTS<-as.data.frame(t(precipTS))
  precipTS = as.data.frame(precipTS[-1,])
  precipTS$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
  # corr of polygons if multiple  
  cor(precipTS[,1:(ncol(precipTS)-1)], method = "pearson", use = "complete.obs")
  # weighted mean of polygons by area
  wgts<-LRU@data[["AREA"]]/sum(LRU@data[["AREA"]])
  precipTS$wgtMean<-apply(as.matrix(precipTS[,1:(ncol(precipTS)-1)]), 1, function(x) weighted.mean(x, wgts))
# pet  
  petTS<-as.data.frame(t(petTS))
  petTS = as.data.frame(petTS[-1,])
  petTS$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
  # corr of polygons if multiple  
  cor(petTS[,1:(ncol(petTS)-1)], method = "pearson", use = "complete.obs")
  # weighted mean of polygons by area
  wgts<-LRU@data[["AREA"]]/sum(LRU@data[["AREA"]])
  petTS$wgtMean<-apply(as.matrix(petTS[,1:(ncol(petTS)-1)]), 1, function(x) weighted.mean(x, wgts))
# nClimGrid ONLY
  nprecipTS<-as.data.frame(t(nprecipTS))
  nprecipTS = as.data.frame(nprecipTS[-1,])
  nprecipTS$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
  # corr of polygons if multiple  
  cor(nprecipTS[,1:(ncol(nprecipTS)-1)], method = "pearson", use = "complete.obs")
  # weighted mean of polygons by area
  wgts<-LRU@data[["AREA"]]/sum(LRU@data[["AREA"]])
  nprecipTS$wgtMean<-apply(as.matrix(nprecipTS[,1:(ncol(nprecipTS)-1)]), 1, function(x) weighted.mean(x, wgts))
  
  # precip compare
  # precipData<-cbind.data.frame(precipTS$date, precipTS$wgtMean, nprecipTS$wgtMean)
  #   precipData$precip3mo<-movingFun(precipData$`precipTS$wgtMean`,3,type="to",sum)
  #   precipData$nprecip3mo<-movingFun(precipData$`nprecipTS$wgtMean`,3,type="to",sum)
  #   precipData$month<-as.numeric(format(precipData$`precipTS$date`, format="%m"))
  # precipData$diff<-precipData$precip3mo-precipData$nprecip3mo
  # qplot(precipData$`precipTS$date`,precipData$diff, geom = "line")
  # qplot(precipData$precip3mo,precipData$nprecip3mo)
  #   temp<-subset(precipData, month==9)
  #   qplot(temp$`precipTS$date`,temp$diff, geom="line")
  # # raster plots
  #   nprecip<-crop(nprecip, extent(precip))
  #   plot(sum(nprecip[[379:381]]), zlim=c(0,600))
  #   plot(sum(precip[[379:381]]), zlim=c(0,600))
    
# weighted mean climate dataframe
moClimData<-cbind.data.frame(precipTS$date, precipTS$wgtMean, petTS$wgtMean)
  moClimData$precip3mo<-movingFun(moClimData$`precipTS$wgtMean`, 3,type = "to", sum)
  moClimData$pet3mo<-movingFun(moClimData$`petTS$wgtMean`, 3,type = "to", sum)
  moClimData$SPI3<-spi(moClimData$`precipTS$wgtMean`, 3)$fitted  
  moClimData$SPEI3<-spei(moClimData$`precipTS$wgtMean`- moClimData$`petTS$wgtMean`, 3)$fitted
  moClimData$month<-as.numeric(format(moClimData$`precipTS$date`, format="%m"))
  moClimData$year<-as.numeric(format(moClimData$`precipTS$date`, format="%Y"))
  
# assess relationships between rpms and climate ts
  library(Hmisc)
trimClim<-subset(moClimData, `precipTS$date`>="1984-01-01")
  resultsList<-list()
  
  for(i in 1:12){
    temp<-subset(trimClim, month==i)
    temp$rpms<-rpmsTS$wgtMean
    res<-rcorr(as.matrix(temp[,2:9]))
    resultsList$corrList[[i]]<-res$r
    resultsList$pList[[i]]<-res$P
    # precip and pet model
    model<-lm(rpms~precip3mo+pet3mo, data=temp) # these are likely autocorrelated
      resultsList$coeffList[[i]]<-summary(model)$coefficients
      resultsList$r2List[i]<-summary(model)$adj.r.squared
      fstat <- summary(model)$fstatistic
      resultsList$fList[i]<- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
    # spi model
    model<-lm(rpms~SPI3, data=temp)
      resultsList$SPIr2List[i]<-summary(model)$adj.r.squared
    # spei model
    model<-lm(rpms~SPEI3, data=temp)
      resultsList$SPEIr2List[i]<-summary(model)$adj.r.squared
  }
    
# format output
# lapply(modelset, function (x) x[c('likelihood', 'fixef')])  
# resultsList[["corrList"]][[1]][,8]
rm(r2results)
r2results<-cbind.data.frame(seq(1,12,1),resultsList[["SPIr2List"]], resultsList[["SPEIr2List"]])

# regression diagnostics
library(performance)
library(lindia)
  temp<-subset(trimClim, month==8)
    temp$rpms<-rpmsTS$wgtMean
    model<-lm(rpms~SPI3, data=temp)
    summary(model)
    check_model(model) # performance
    gg_diagnose(model) # lindia
    model_performance(model) # performance
  # plot time series plot
    # add 'fit', 'lwr', and 'upr' columns to dataframe (generated by predict)
    temp.predict <- cbind(temp, predict(model, interval = 'confidence'))
    
    # plot the points (actual observations), regression line, and confidence interval
    p <- ggplot(temp.predict, aes(SPI3,rpms))
    p <- p + geom_point()
    p <- p + geom_line(aes(SPI3, rpms))
     p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3)+
      ggtitle("41-3AZ Annual Production (RPMS)v JAS-SPI3, nClimGrid")
    # time series
    p<-ggplot(temp.predict, aes(as.Date(temp.predict$`precipTS$date`),rpms))+
       geom_line()
    p+geom_line(aes(as.Date(temp.predict$`precipTS$date`),fit), color='red')+
      ggtitle("41-3AZ RPMS Production and JAS-SPI3 nClimGrid (obs=black, pred=red)")
    ggplot(temp.predict, aes(as.Date(temp.predict$`precipTS$date`),rpms-fit))+
      geom_line()+
      ggtitle("41-3AZ RPMS Production by JAS-SPI3 nClimGrid model  (obs-pred)")+
      geom_hline(yintercept=0)
    
 # model cross validation 
 #  http://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/#k-fold-cross-validation  
    library(caret)
    # Define training control
    set.seed(123) 
    train.control <- trainControl(method = "repeatedcv", number = 5, repeats = 3)
    #train.control <- trainControl(method = "LOOCV")
    # Train the model
    model2 <- train(rpms~SPI3, data=temp, method = "lm",
                   trControl = train.control)
    # Summarize the results
    print(model2)
  summary(model2)

# all subsets to build models with all seasons
  allSeas<-moClimData[,c(8,9,7)]
  allSeas<-dcast(allSeas, year~month)
  colnames(allSeas)<-c("year",paste0(month.abb[11],"_",month.abb[1]),paste0(month.abb[12],"_",month.abb[2]),
                               paste0(month.abb[seq(1,12-2,1)],"_",month.abb[seq(1+2,12,1)]))
# precip and pet
  precip3mo<-moClimData[,c(8,9,4)]
    precip3mo<-dcast(precip3mo, year~month)
    colnames(precip3mo)<-c("year",paste0("p_",month.abb[11],"_",month.abb[1]),paste0("p_",month.abb[12],"_",month.abb[2]),
                         paste0("p_",month.abb[seq(1,12-2,1)],"_",month.abb[seq(1+2,12,1)]))
  pet3mo<-moClimData[,c(8,9,5)]
    pet3mo<-dcast(pet3mo, year~month)
    colnames(pet3mo)<-c("year",paste0("pet_",month.abb[11],"_",month.abb[1]),paste0("pet_",month.abb[12],"_",month.abb[2]),
                           paste0("pet_",month.abb[seq(1,12-2,1)],"_",month.abb[seq(1+2,12,1)]))
  # or just precip and pet alone
  allSeas<-cbind.data.frame(precip3mo,pet3mo[,2:13])
  # water balance
  allSeas<-cbind.data.frame(precip3mo[,1],precip3mo[,2:13]-pet3mo[,2:13])
    colnames(allSeas)[1]<-"year"
  
  # regressions
  # http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/155-best-subsets-regression-essentials-in-r/
   trimSeas<-subset(allSeas, year>=1984)
    trimSeas$rpms<-rpmsTS$wgtMean
 rm(models, model)
 library(leaps)
    models <- regsubsets(rpms~., data = trimSeas, nvmax = 12)
  summary(models)  
  res.sum <- summary(models)
  data.frame(
    Adj.R2 = which.max(res.sum$adjr2),
    CP = which.min(res.sum$cp),
    BIC = which.min(res.sum$bic)
  )
  
  model<-lm(rpms~Jun_Aug, data=trimSeas)
  summary(model)
    
  # best regression results
  # id: model id
  # object: regsubsets object
  # data: data used to fit regsubsets
  # outcome: outcome variable
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
  
  get_model_formula(1, models, "rpms")
  
  get_cv_error <- function(model.formula, data){
    set.seed(1)
    train.control <- trainControl(method = "cv", number = 5)
    cv <- train(model.formula, data = data, method = "lm",
                trControl = train.control)
    cv$results$RMSE
  }
  
  # Compute cross-validation error
  library(tidyverse)
  model.ids <- 1:12
  cv.errors <-  map(model.ids, get_model_formula, models, "rpms") %>%
    map(get_cv_error, data = trimSeas) %>%
    unlist()
  cv.errors

  # Select the model that minimize the CV error
  which.min(cv.errors)
    
  coef(models, which.min(cv.errors))
  
# # stuff to try
# x assess how much of LRU is actually covered by RPMS data...if not the corrs/models will be for different areas
# x mask LRU for RPMS coverage first? Then extract clim vars so only contributing areas are included
# - assess variability within and between LRUs with climate, rpms, VHI...are these good units to summarize drought conditions?
# - look at random forest or structure to assess multiple seasons contributions to variance in rpms?
# - is there a multivariate combo of timescales and vars that outperforms simple one timescale SPI or SPEI?
# - all subsets regression, but check to see if time series are parametric
# - try nClimDiv as precip...is it better than daily PRISM? 
# - quantreg with all cell values from RPMS and GridMet
# - develop drought tracker for each LRU - historical drought summaries, near-real time drought summaries?
# - use existing county extract range res code for LRUs