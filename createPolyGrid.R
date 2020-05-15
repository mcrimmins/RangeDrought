# create polygon grid for regression modeling 
# MAC 05/06/2020

library(raster)
library(rgdal)

# set rasteroptions
rasterOptions(progress = 'text')

# load states
states<-getData('GADM', country='USA', level=1)
state<-subset(states, NAME_1=="Arizona" )

# get gridMet elevation data
#download.file("https://climate.northwestknowledge.net/METDATA/data/metdata_elevationdata.nc", destfile = "gridMet_elev.nc")

pet<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_pet_1979_2019.grd")
gridMet<-pet[[1]]

# grab gridmet 
#gridMet<-raster("gridMet_elev.nc")
  grid<- gridMet[[1]]
  grid<- crop(grid, extent(state))
  grid[grid==0]<-NA
  grid<-aggregate(grid,3, function(i, ...) sum(!is.na(i))) # 12 is 0.5 degree grid, 6 ~ 0.25 deg grid
  grid[grid < 9] <- NA # 144 for factor 12, 36 for factor 6, 9 for 3
  polyGrid<-rasterToPolygons(grid, na.rm = TRUE)
  
  library(rgdal)
  writeOGR(obj=polyGrid, dsn="shapes", layer="polyGrid_eighthDeg", driver="ESRI Shapefile", overwrite_layer = TRUE) # this is in geographical projection
  