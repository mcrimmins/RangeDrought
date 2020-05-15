# read in USGS Pheno RS data
# https://phenology.cr.usgs.gov/get_data_1km.php
# MAC 05/06/2020

library(raster)

# set rasteroptions
rasterOptions(progress = 'text')

# get info from VHI/NDVI data
vhi<-stack('/scratch/crimmins/vhi/processed/VHI_complete_1982-2019_SWUS.grd')

# loop through to get stack of ann max NDVI
stackGrid<-stack()
for (i in seq(2012,2012, by=1)) {
  tempRast<-raster(paste0("/scratch/crimmins/USGS_Pheno/annMax/av_MAXN",i,"v4.bsq"), 
                   crs="+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
    tempRast<-projectRaster(tempRast, crs = crs(vhi[[1]]))
  e <- extent(-116.2267,-99.21573,28.9193,42.76954) # WESTUS extent(-125, -100, 25, 49)
  tempRast <- crop(tempRast, e)
  tempRast[tempRast>=254]<-NA # water bodies
  #tempRast[tempRast==0]<-NA # No data pixels are 0, set to NA if necessary
  stackGrid <- stack(stackGrid, tempRast)
  print(names(tempRast))
}

# To revert to unscaled NDVI, a scale factor ((MAXN-100)*0.01 )
# file:///C:/Users/Crimmins/AppData/Local/Temp/Temp1_av_MAXN2000v4.zip/av_maxn2000v4_metadata.htm


# write Raster to file
writeRaster(stackGrid,filename='/scratch/crimmins/USGS_Pheno/processed/av_MAXN_1989-2014_SWUS.grd', overwrite=TRUE)
