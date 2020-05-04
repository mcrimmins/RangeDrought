# use NOAA NESDIS NDVI data
# update data files with getVHI.R
# process VHI/NDVI into raster files with processVHI.R


library(raster)
library(tidyr)

# set rasteroptions
rasterOptions(progress = 'text')

smNDVI<-stack('/scratch/crimmins/vhi/processed/smNDVI_complete_1982-2019_SWUS.grd')

dates<-as.data.frame(names(smNDVI))
colnames(dates)<-"dateCode"
  dates<-separate(dates,dateCode,c("X","date"), sep ="X")
  dates<-separate(dates,date,c("year","month","day"))
  dates$fullDate<-as.Date(paste0(dates$year,"-",dates$month,"-",dates$day), format="%Y-%m-%d")
  dates$week<-as.numeric(strftime(dates$fullDate, format = "%V"))
  dates$doy<-as.numeric(format(dates$fullDate,"%j"))
  # max NDVI
  maxNDVI<- stackApply(smNDVI, (as.integer(dates$year)-1981), fun = max, na.rm = TRUE)
  names(maxNDVI)<-seq(1982,2019,1)
  # write Raster to file
  writeRaster(maxNDVI,filename='/scratch/crimmins/vhi/processed/ANNUAL_MAX_smNDVI_complete_1982-2019_SWUS.grd', overwrite=TRUE)
  
  # which max for each year - timing of max NDVI
  which.max2 <- function(x, ...) ifelse( length(x) ==sum(is.na(x) ), NA, which.max(x))          # helper function to absorb "na.rm"
  whenMaxNDVI <- stackApply(smNDVI, (as.integer(dates$year)-1981), fun=which.max2, na.rm=NULL)
  # write Raster to file
  writeRaster(whenMaxNDVI,filename='/scratch/crimmins/vhi/processed/ANNUAL_MAX_timing_smNDVI_complete_1982-2019_SWUS.grd', overwrite=TRUE)
  

# sum NDVI - like iNDVI? problem with year with missing periods
  sumNDVI<- stackApply(smNDVI, (as.integer(dates$year)-1981), fun = sum, na.rm = TRUE)
  names(sumNDVI)<-seq(1982,2019,1)
  
# write Raster to file
  writeRaster(sumNDVI,filename='/scratch/crimmins/vhi/processed/ANNUAL_SUM_smNDVI_complete_1982-2019_SWUS.grd', overwrite=TRUE)
  
# maps of max NDVI timing
library(rasterVis)
library(rgdal)
library(dplyr)
  # ---- load NRCS CRAs for SW
  AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
  NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
  # combine
  SW_CRA<- rbind(AZ_CRA,NM_CRA)
  
  # map of timing for diff years
  names(whenMaxNDVI)<-seq(1982,2019,1)
  levelplot(whenMaxNDVI[[19:21]], margin=FALSE)+  
    layer(sp.polygons(SW_CRA))
  
  # average max timing
  avgMaxTiming<-calc(whenMaxNDVI, mean)

  # dates from above, get average doy for week
  avgDOY <- dates %>% group_by(week) %>% 
                      summarise(doy=mean(doy))
  avgDOY$doy[1]<-2
  avgDOY$doy<-round(avgDOY$doy,0)
  # round values to closest week integer
  avgMaxTiming<-round(avgMaxTiming,0)
  
  # reclassify to day of year
  # cheap way to multiply by 7, or reclassify with avgDOY
  levelplot(avgMaxTiming*7, margin=FALSE, main="Mean Day of Year of maxNDVI")+  
    layer(sp.polygons(SW_CRA, col="grey"))
  
  # max day of year by LMU
  
  
      
    
    
  