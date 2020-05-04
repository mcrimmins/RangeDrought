# map CRAs across SW US
# MAC 03/24/2020

library(rgdal)
library(leaflet)
library(htmltools)
library(htmlwidgets)
library(raster)

# ---- load NRCS CRAs for SW
AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
CA_CRA <- readOGR(dsn = "./shapes/CA_CRA/soils", layer = "cra_a_ca")
CO_CRA <- readOGR(dsn = "./shapes/CO_CRA/soils", layer = "cra_a_co")
NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
NV_CRA <- readOGR(dsn = "./shapes/NV_CRA/soils", layer = "cra_a_nv")
UT_CRA <- readOGR(dsn = "./shapes/UT_CRA/soils", layer = "cra_a_ut")
# combine
SW_CRA<- rbind(AZ_CRA,CA_CRA,CO_CRA,NM_CRA,NV_CRA,UT_CRA)

# AZ LRU from Emilio
AZ_LRU<-readOGR(dsn = "./shapes/AZ_LRU", layer = "mlra_a_az")
  AZ_LRU<-spTransform(AZ_LRU, crs(SW_CRA))

# AZ and NM MLRA from NRCS gateway
AZ_MLRA<-readOGR(dsn = "./shapes/AZ_MLRA/soils", layer = "mlra_a_az")
NM_MLRA<-readOGR(dsn = "./shapes/NM_MLRA/soils", layer = "mlra_a_nm")
  
# plot maps
leafmap2<-leaflet() %>% addProviderTiles(providers$OpenTopoMap) %>%
  addPolygons(data=SW_CRA, color = "#444444", weight = 0.5, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.1,
              highlightOptions = highlightOptions(color = "red", weight = 2,
                                                  bringToFront = TRUE),
              label=as.character(SW_CRA$CRA),
              labelOptions = labelOptions(direction = "auto"))

saveWidget(leafmap, file="./SW_NRCS_CRAs.html", selfcontained = FALSE)

leafmap2<-leaflet() %>% addProviderTiles(providers$OpenTopoMap) %>%
  addPolygons(data=AZ_LRU, color = "#444444", weight = 0.5, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.1,
              highlightOptions = highlightOptions(color = "red", weight = 2,
                                                  bringToFront = TRUE),
              label=as.character(AZ_LRU$LRU),
              labelOptions = labelOptions(direction = "auto"))

saveWidget(leafmap2, file="/home/crimmins/RProjects/RangeDrought/maps/AZ_NRCS_LRUs.html", selfcontained = TRUE)
