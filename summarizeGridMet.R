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
prefix<-"pet"

for(yr in 1979:2019){
  gridMet<-stack(paste0("/scratch/crimmins/gridmet/update_Aug2019/",prefix,"_",yr,".nc"))
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
writeRaster(monthSums, filename = "/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_pet_1979_2019.grd",
            overwrite=TRUE)

# create SPI/SPEI grids
library(SPEI)
precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_precip_1979_2019.grd")
pet<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_pet_1979_2019.grd")
# monthly water bal
wtrbal<-precip-pet

# dates
dates=seq(as.Date("1979-01-01"), as.Date("2019-12-01"), by="month")
names(precip)<-dates

# calculate SPI ----
funSPI <- function(x, scale=6, na.rm=TRUE,...) as.numeric((spi(x, scale=scale, na.rm=na.rm, ...))$fitted)
# parallell calc
ptm <- proc.time()
  beginCluster(6)
    spi<- clusterR(precip, calc, args=list(fun=funSPI))
  endCluster()
proc.time() - ptm

# save raster
names(spi)<-dates
writeRaster(spi, filename = "/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_6mo_SPI_1979_2019.grd",overwrite=TRUE)
# ----

# calculate SPEI ----
funSPEI <- function(x, scale=6, na.rm=TRUE,...) as.numeric((spei(x, scale=scale, na.rm=na.rm, ...))$fitted)
# parallell calc
ptm <- proc.time()
  beginCluster(6)
    spei <- clusterR(wtrbal, calc, args=list(fun=funSPEI))
  endCluster()
proc.time() - ptm

# save raster
names(spei)<-dates
writeRaster(spei, filename = "/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_6mo_SPEI_1979_2019.grd",overwrite=TRUE)
# ----




# interactive plots of data
# precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_precip_1979_2019.grd")
# library(leaflet)
# library(leafem)
# library(htmlwidgets)
# grid<-sum(precip[[1:12]])
# pal=colorNumeric(c("blue", "yellow", "orange"), values(grid),
#                  na.color = "transparent")
# leaflet() %>% addProviderTiles(providers$OpenTopoMap)%>%
#   addRasterImage(grid, colors = pal, maxBytes=20000000,opacity = 0.8, layerId = "precip") %>%
#   addMouseCoordinates() %>%
#   addImageQuery(grid, type="mousemove", layerId = "precip", digits = 2, prefix = "") %>%
#   addLegend(pal = pal, values = values(grid),
#             title="precip")

