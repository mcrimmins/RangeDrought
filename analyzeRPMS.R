# analyze RPMS with climate data
# relies on data files created in processRPMS.R and summarizeGrid
# MAC 03/25/2020

library(raster)
library(rgdal)
library(rasterVis)

# set rasteroptions
rasterOptions(progress = 'text')

# load states
states<-getData('GADM', country='USA', level=1)  

# ---- load NRCS CRAs for SW
AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
# combine
SW_CRA<- rbind(AZ_CRA,NM_CRA)

# # ---- load NRCS CRAs for SW
# AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
# #CA_CRA <- readOGR(dsn = "./shapes/CA_CRA/soils", layer = "cra_a_ca")
# #CO_CRA <- readOGR(dsn = "./shapes/CO_CRA/soils", layer = "cra_a_co")
# NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
# #NV_CRA <- readOGR(dsn = "./shapes/NV_CRA/soils", layer = "cra_a_nv")
# #UT_CRA <- readOGR(dsn = "./shapes/UT_CRA/soils", layer = "cra_a_ut")
# # combine
# #SW_CRA<- rbind(AZ_CRA,CA_CRA,CO_CRA,NM_CRA,NV_CRA,UT_CRA)
# SW_CRA<- rbind(AZ_CRA,NM_CRA)

# REPROJECT RPMS data ----
# load data
# precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_precip_1979_2019.grd")
# rpmsSW<-stack("./temp/SWUS_RPMS_8419.grd")
# 
# # clusterR parallel
#   library(cluster)
#   beginCluster(type="SOCK",n=3)
#   ## 6 out of 8 cores
#     rpmsSW_WGS84 <-projectRaster(rpmsSW,crs=crs(precip[[1]]), file="./temp/SWUS_RPMS_8419_WGS84.grd")
#   ## Using cluster with 3 nodes
#   endCluster() 
# moved _WGS84 version back to scratch drive
# ----  
  
# # load reprojected RPMS data and gridmet ----
# rpmsSW<-stack("/scratch/crimmins/USDA_NPP/v2019/SWUS_RPMS_8419_WGS84.grd")
# precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/SWUS_gridmet_monthly_sum_precip_1979_2019.grd")
# 
 # levelplot((precip[[1]]), margin=FALSE)+
 #   layer(sp.polygons(states))
# 
# # crop to common AZ/NM bounding box
# # -116 to -102 lon, 37.5 to 30.5 lat
#   e <- extent(-116,-102,30.5,37.5)
#   rpmsSW <- crop(rpmsSW, e)
#   precip<- crop(precip, e)
# # resmple GridMet to RPMS res
#   # # clusterR parallel
#   library(cluster)
#      beginCluster(type="SOCK",n=4)
#      ## 3 out of 8 cores
#       precip <- resample(precip,rpmsSW[[1]],method='bilinear')
#         ## Using cluster with 3 nodes
#      endCluster()
#      # save cropped raster
#      writeRaster(precip, filename = "/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_monthly_sum_precip_1979_2019.grd", overwrite=TRUE)
  # rpmsSW<-stack("/scratch/crimmins/USDA_NPP/v2019/SWUS_RPMS_8419_WGS84.grd")
  #  e <- extent(-116,-102,30.5,37.5)
  # rpmsSW <- crop(rpmsSW, e)     
  # writeRaster(rpmsSW, 
  #             filename = "/scratch/crimmins/USDA_NPP/v2019/AZNM_RPMS_8419_WGS84.grd", overwrite=TRUE)
# ----

# summary stats/diagnostics of rpms
rpmsSW<-stack("/scratch/crimmins/USDA_NPP/v2019/AZNM_RPMS_8419_WGS84.grd")

minYear<-(which.min(rpmsSW))+1983
 levels(minYear)=data.frame(ID=1984:2019, code=1984:2019)
levelplot(minYear, par.settings = PuOrTheme, main="Year of Min RPMS") +
  layer(sp.polygons(states))   
#----
# leaflet
library(rgdal)
library(leaflet)
library(leafem)
library(htmlwidgets)
pal=colorNumeric(c("blue", "yellow", "orange"), values(minYear),
                 na.color = "transparent")
leafMap<-leaflet() %>% addTiles() %>%
  addRasterImage(minYear, colors = pal, maxBytes=20000000,opacity = 0.8, layerId = "minYear") %>%
  addMouseCoordinates() %>%
  addImageQuery(minYear, type="mousemove", layerId = "minYear", digits = 2, prefix = "") %>%
  addLegend(pal = pal, values = values(minYear),
            title="Min Year - RPMS")%>%
  addPolygons(data=SW_CRA, color = "#444444", weight = 0.5, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.1,
              highlightOptions = highlightOptions(color = "red", weight = 2,
                                                  bringToFront = TRUE),
              label=as.character(SW_CRA$CRA),
              labelOptions = labelOptions(direction = "auto"))
saveWidget(leafMap, file="/home/crimmins/RProjects/RangeDrought/maps/min_RPMS_year.html", selfcontained = TRUE)

maxYear<-(which.max(rpmsSW))+1983
levels(maxYear)=data.frame(ID=1984:2019, code=1984:2019)
levelplot(maxYear, par.settings = PuOrTheme, main="Year of Max RPMS") +
  layer(sp.polygons(states))  
 
# load reprojected RPMS data and gridmet
precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_monthly_sum_precip_1979_2019.grd")
rpmsSW<-stack("/scratch/crimmins/USDA_NPP/v2019/AZNM_RPMS_8419_WGS84.grd")

    
# summarize precip data into 
# dates<-seq(as.Date("1979-01-01",format="%Y-%m-%d"),as.Date("2019-12-01",format="%Y-%m-%d"),
#            by="month")
# month<-as.data.frame(as.numeric(format(dates, "%m")))
#   colnames(month)<-"month"
#   month$year<-as.numeric(format(dates, "%Y"))
# library(dplyr)
# season<-month %>%
#   mutate(
#     season = case_when(
#       month %in% 10:12 ~ 4,
#       month %in%  1:3  ~ 1,
#       month %in%  4:6  ~ 2,
#       TRUE ~ 3))
# seasCode<-paste0(season$season,"-",month$year)

  #seasSums<- stackApply(precip, rep(1:(nlayers(precip)/3), each=3), fun = sum, na.rm = TRUE)
  #writeRaster(seasSums, filename = "/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_seasonal_sum_precip_1979_2019.grd", overwrite=TRUE)
  
  # correlate summer precip with annual production
  seasSums<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_seasonal_sum_precip_1979_2019.grd")
  # extract summer precip for rpms por
  names(seasSums)<-seq(as.Date("1979-01-01","%Y-%m-%d"),as.Date("2019-12-01","%Y-%m-%d"),by="quarter")
  seasPrecip<-seasSums[[seq(23,163,4)]] # summer
  #seasPrecip<-seasSums[[seq(22,162,4)]] # summer
  
  
  # https://gis.stackexchange.com/questions/279535/pixel-correlations-from-two-raster-datasets-in-r
  # combine your two raster stacks. Now layers 1:(n/2) will be 
  # correlated with layers ((n/2)+1):n. Make sure both stacks have the same 
  # number of layers, i.e. nlayers(l_all) is even!
  
  #mskSeasPrecip = mask(seasPrecip, rpmsSW)
  #varStack<-stack( mskSeasPrecip, rpmsSW)
  varStack<-stack(seasPrecip,rpmsSW)
  
  # write a little function to do the correlation. The input is the vector
  # of values at cell n, which has length() == nlayers(l_all)
  corvec <- function(vec = NULL) {
    cor(
      # 'top' half of stack
      x      = vec[1:(length(vec)/2)],
      # 'bottom' half of stack
      y      = vec[((length(vec)/2) + 1):length(vec)],
      use    = 'na.or.complete',
      method = 'pearson'
    )
  }
  
  corlyrs <- calc(
    varStack,
    fun = function(x) {
      # handle areas where all cell vals are NA, e.g. ocean
      if (all(is.na(x))) {
        NA_real_
      } else {
        corvec(vec = x)
      }
    }
  )
  
  library(cluster)
  beginCluster(n = 4)
  corlyrs <- clusterR(x         = varStack,
                      fun       = calc,
                      args      = list(fun = function(cell) {
                        if (all(is.na(cell))) {
                          NA_real_
                        } else {
                          corvec(cell)
                        }
                      }),
                      export    = 'corvec'
  )
  endCluster()
  
   levelplot(corlyrs, margin=FALSE, par.settings = BuRdTheme, at=seq(-1, 1, 0.1))+
     layer(sp.polygons(states))
  
# linear model of precip/rpms
  # slope  
  #fun_slope<-function(x) {lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)])$coefficients[2];}
  fun_slope<-function(x) { m<-lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)]);summary(m)$r.squared}
  beginCluster(n = 5)
  slope <- clusterR(x         = varStack,
                    fun       = calc,
                    args      = list(fun = function(x) {
                      if (sum(is.na(x))>=(length(x)/2)) {
                        NA
                      } else {
                        fun_slope(x)
                      }
                    }),
                    export    = 'fun_slope'
  )
  endCluster()

  levelplot(slope, margin=FALSE, par.settings = BuRdTheme, at=seq(0, 1, 0.05))+
    layer(sp.polygons(states))   

  # leaflet
  library(rgdal)
  library(leaflet)
  library(leafem)
  library(htmlwidgets)
  
  # ---- load NRCS CRAs for SW
  AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
  NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
  # combine
  SW_CRA<- rbind(AZ_CRA,NM_CRA)
  pal=colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(slope),
               na.color = "transparent")
  leafMap<-leaflet() %>% addTiles() %>%
    addRasterImage(slope, colors = pal, maxBytes=20000000,opacity = 0.8, layerId = "R2") %>%
    addMouseCoordinates() %>%
    addImageQuery(slope, type="mousemove", layerId = "R2", digits = 2, prefix = "") %>%
    addLegend(pal = pal, values = values(slope),
              title="JAS Precip vs RPMS - R2")%>%
    addPolygons(data=SW_CRA, color = "#444444", weight = 0.5, smoothFactor = 0.5,
                opacity = 1.0, fillOpacity = 0.1,
                highlightOptions = highlightOptions(color = "red", weight = 2,
                                                    bringToFront = TRUE),
                label=as.character(SW_CRA$CRA),
                labelOptions = labelOptions(direction = "auto"))
  saveWidget(leafMap, file="./rpms_precip_test.html", selfcontained = FALSE)
      
   # ## p-value
   # fun_pvalue <- function(y) { 
   #   if(all(is.na(y))) {
   #     NA
   #   } else {
   #     m = lm(y ~ Date); summary(m)$coefficients[8] 
   #   }
   # }
   # 
   # # intercept
   # fun=function(x) { if (is.na(x[1])){ NA } else { lm(x[1:120] ~ x[121:240])$coefficients[1] }}
   # intercept <- calc(s, fun)
   # # r2
   # fun=function(x) { if (is.na(x[1])){ NA } else { m <- lm(x[1:120] ~ x[121:240]);summary(m)$r.squared }}
   # r.squared <- calc(s, fun)
   # 
   
  #showTmpFiles()
  #removeTmpFiles(h=0) # in hours
  