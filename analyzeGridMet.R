# analyze climate patterns with GridMet data
# MAC 04/04/20

library(raster)
library(rasterVis)
library(rgdal)

# set rasteroptions
rasterOptions(progress = 'text')

# load states
states<-getData('GADM', country='USA', level=1)  

# # ---- load NRCS CRAs for SW
AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
#CA_CRA <- readOGR(dsn = "./shapes/CA_CRA/soils", layer = "cra_a_ca")
CO_CRA <- readOGR(dsn = "./shapes/CO_CRA/soils", layer = "cra_a_co")
NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
#NV_CRA <- readOGR(dsn = "./shapes/NV_CRA/soils", layer = "cra_a_nv")
UT_CRA <- readOGR(dsn = "./shapes/UT_CRA/soils", layer = "cra_a_ut")
# combine
SW_CRA<- rbind(AZ_CRA,CO_CRA,NM_CRA,UT_CRA)
#SW_CRA<- rbind(AZ_CRA,NM_CRA)

# load broader SW clipped data
precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_precip_1979_2019.grd")
pet<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_pet_1979_2019.grd")
# switch var
climVar<-pet

# get average annnual total precip
dates<-seq(as.Date(paste0(1979,"-01-01"),format="%Y-%m-%d"),as.Date(paste0(2019,"-12-01"),format="%Y-%m-%d"),by="month")
  month<-as.numeric(format(dates,"%m"))
  year<-as.numeric(format(dates,"%Y"))-1978
  annPrecip<- (stackApply(climVar, year, fun = sum, na.rm = TRUE))
  annAvgPrecip<- calc(annPrecip, mean)

# create seasonal means
library(cluster)
  beginCluster(n=6)
    seasPrecip <- clusterR(climVar,fun=calc,filename="./temp/precip_moving_sum.grd",
                args=list(fun = function(x) movingFun(x, 3,type = "to", sum)))
  endCluster()

# get long-term seasonal means
  seasPrecipMean<- (stackApply(seasPrecip, month, fun = mean, na.rm = TRUE))
  
# 3-month labels
  dateSeq<-c(paste0(month.abb[11],"-",month.abb[1]),paste0(month.abb[12],"-",month.abb[2]),
             paste0(month.abb[seq(1,12-2,1)],"-",month.abb[seq(1+2,12,1)]))
  temp<-(seasPrecipMean/annAvgPrecip)*100
  names(temp)<-dateSeq
  p0<-levelplot(temp, margin=FALSE, par.settings = viridisTheme, at=seq(1, 70, 5), main="Percent of total ann PET by 3-mo season (1984-2019)")+
    layer(sp.polygons(states))
  # plot to png
  png("/home/crimmins/RProjects/RangeDrought/figs/GridMet_seasonality_perc_of_annual_PET.png", width = 10, height = 8, units = "in", res = 300L)
  print(p0, newpage = FALSE)
  dev.off()
