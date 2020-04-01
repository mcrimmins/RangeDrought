# summarize gridmet data for monthly calculations
# MAC 03/24/2020

library(raster)
library(rgdal)

# set rasteroptions
rasterOptions(progress = 'text')

# load states
states<-getData('GADM', country='USA', level=1)  

# ---- load NRCS CRAs for SW
AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
#CA_CRA <- readOGR(dsn = "./shapes/CA_CRA/soils", layer = "cra_a_ca")
#CO_CRA <- readOGR(dsn = "./shapes/CO_CRA/soils", layer = "cra_a_co")
NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
#NV_CRA <- readOGR(dsn = "./shapes/NV_CRA/soils", layer = "cra_a_nv")
#UT_CRA <- readOGR(dsn = "./shapes/UT_CRA/soils", layer = "cra_a_ut")
# combine
#SW_CRA<- rbind(AZ_CRA,CA_CRA,CO_CRA,NM_CRA,NV_CRA,UT_CRA)
SW_CRA<- rbind(AZ_CRA,NM_CRA)


# create monthly sums for region
monthStack=list()
i<-1

for(yr in 1979:2019){
  gridMet<-stack(paste0("/scratch/crimmins/gridmet/update_Aug2019/pr_",yr,".nc"))
  # crop to region
  gridMet <- crop(gridMet, extent(SW_CRA))
  dates<-seq(as.Date(paste0(yr,"-01-01"),format="%Y-%m-%d"),as.Date(paste0(yr,"-12-31"),format="%Y-%m-%d"),1)
  month<-as.numeric(format(dates,"%m"))
  monthStack[[i]]<- stackApply(gridMet, month, fun = sum, na.rm = TRUE)
  print(yr)
  i=i+1
}

# combine into stack
monthSums = stack(monthStack)

# save cropped raster
writeRaster(monthSums, filename = "/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_precip_1979_2019.grd",
            overwrite=TRUE)

