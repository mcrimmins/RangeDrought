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
  
  # AZ LRU from Emilio
  AZ_LRU<-readOGR(dsn = "./shapes/AZ_LRU", layer = "mlra_a_az")
  AZ_LRU<-spTransform(AZ_LRU, crs(AZ_CRA))
  
  # load stack
  whenMaxNDVI<-stack('/scratch/crimmins/vhi/processed/ANNUAL_MAX_timing_smNDVI_complete_1982-2019_SWUS.grd')
   
  # map of timing for diff years
  names(whenMaxNDVI)<-seq(1982,2019,1)
  levelplot(whenMaxNDVI[[19:21]], margin=FALSE)+  
    layer(sp.polygons(SW_CRA))
  
  # average max timing
  avgMaxTiming<-calc(whenMaxNDVI, mean)

  # dates from above, get average doy for week
  #avgDOY <- dates %>% group_by(week) %>% 
  #                    summarise(doy=mean(doy))
  #avgDOY$doy[1]<-2
  #avgDOY$doy<-round(avgDOY$doy,0)
  # round values to closest week integer
  avgMaxTiming<-round(avgMaxTiming,0)
  
  # reclassify to day of year
  # cheap way to multiply by 7, or reclassify with avgDOY
  levelplot(avgMaxTiming*7, margin=FALSE, main="Mean Day of Year of maxNDVI")+  
    layer(sp.polygons(SW_CRA, col="grey"))
  
  # max day of year by LMU
    # switch to land management unit
    LMU<-AZ_LRU
    #LMU<-counties
    # add id code to data
    LMU$ID<-seq(1,nrow(LMU),1)
  
    resultsFrame = data.frame(ID=rep(0, nrow(LMU)),
                              doy=rep(0,nrow(LMU)))
    
    for(i in 1:nrow(LMU)){
      poly<-LMU[i,]
      # summary stat ts for LRU mean or median
      beginCluster(n=6)
        timingTS<-raster::extract(avgMaxTiming, poly, df=TRUE)
      endCluster()
      meanDOY<-round(mean(timingTS[,2], na.rm=TRUE)*7,0)
      resultsFrame[i, ] = c(i,meanDOY)
      print(i)
    }
    
    # merge with spatial dataframe
    LMU@data<-merge(LMU@data, resultsFrame, by.x="ID", by.y="ID")
    
    save(LMU,file="./results/AZ_LRU_smNDVI_maxDOY_map.Rdata")
    
    pal <- RColorBrewer::brewer.pal(n=10, 'YlGnBu')
    spplot(LMU, c("doy"), at=seq(90, 270, 20),
           main=list(label="Avg Max DOY (smNDVI)",cex=1),
           scales=list(draw = TRUE),  col.regions = pal)
    
  # cluster on maxNDVI timing
    library(raster)
    library(RStoolbox)
    library(rasterVis)
    library(rgdal)
    
    # set rasteroptions
    rasterOptions(progress = 'text')
    
    # AZ LRU from Emilio
    # ---- load NRCS CRAs for SW
    AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
    AZ_LRU<-readOGR(dsn = "./shapes/AZ_LRU", layer = "mlra_a_az")
    AZ_LRU<-spTransform(AZ_LRU, crs(AZ_CRA))
    
    # load stack
    whenMaxNDVI<-stack('/scratch/crimmins/vhi/processed/ANNUAL_MAX_timing_smNDVI_complete_1982-2019_SWUS.grd')
    whenMaxNDVI<-crop(whenMaxNDVI, extent(AZ_LRU))
    # average max timing
    avgMaxTiming<-calc(whenMaxNDVI, mean)
    
    
    
    # lat lon grids
    # lonGrid <- init(whenMaxNDVI, 'x')
    # latGrid <- init(whenMaxNDVI, 'y')
    # whenMaxNDVI<-stack(whenMaxNDVI,lonGrid,latGrid)
     
    # rsToolbox options
    rsOpts(verbose=TRUE)
    clusterN<-12
    #unC <- unsuperClass(avgMaxTiming, nSamples = 1000, nClasses = clusterN, nStarts = 25, nIter = 1000, norm = FALSE) # or allGDD
    unC <- unsuperClass(whenMaxNDVI, nSamples = 1000, nClasses = clusterN, nStarts = 25, nIter = 1000, norm = FALSE) # or allGDD
    
    # more colors
    # https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    darkcols<-sample(col_vector, clusterN)
    
    # smooth or clump?
    #test <- focal(classMap, w=matrix(1, 7, 7), fun=min)
    # test<-clump(classMap)
    
    #cluster12 <- readShapePoly("./mapFiles/cluster12poly")
    
    classMap<-as.factor(unC$map)
    rat <- levels(classMap)[[1]]
    # cluster names
    #rat[["cluster"]]<-c("N Rockies","Northeast","S Plains","S Rockies","N Plains",
    #
    #                    "Midwest","Southeast","Mid Atlantic","Upper Midwest","Southwest",
    #                    "Gulf Coast","Pacific NW")
    rat[["cluster"]] <- as.character(seq(1, clusterN, by=1))
    levels(classMap) <- rat 
    
    # create polygon shapefile
    #library(rgdal)
    #library(rgeos)
    #library(maptools)
    # clusterpoly<-rasterToPolygons(classMap, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
    # proj4string(clusterpoly) <- "+proj=longlat +datum=WGS84"
    # writeOGR(clusterpoly, ".", "./mapFiles/cluster12final", driver="ESRI Shapefile")
    
    # plot classified map
    levelplot(classMap, col.regions=darkcols, par.settings=list(panel.background=list(col="white")),
              margin=FALSE, main="smNDVI Max NDVI timing - kmeans clusters")+
      layer(sp.polygons(AZ_LRU))
      
    
    
  