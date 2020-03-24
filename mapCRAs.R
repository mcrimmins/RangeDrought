# map CRAs across SW US
# MAC 03/24/2020

library(rgdal)
library(leaflet)
library(htmltools)
library(htmlwidgets)

# ---- load NRCS CRAs for SW
AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
CA_CRA <- readOGR(dsn = "./shapes/CA_CRA/soils", layer = "cra_a_ca")
CO_CRA <- readOGR(dsn = "./shapes/CO_CRA/soils", layer = "cra_a_co")
NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
NV_CRA <- readOGR(dsn = "./shapes/NV_CRA/soils", layer = "cra_a_nv")
UT_CRA <- readOGR(dsn = "./shapes/UT_CRA/soils", layer = "cra_a_ut")
# combine
SW_CRA<- rbind(AZ_CRA,CA_CRA,CO_CRA,NM_CRA,NV_CRA,UT_CRA)

# plot maps
leafmap<-leaflet() %>% addProviderTiles(providers$OpenTopoMap) %>%
  addPolygons(data=SW_CRA, color = "#444444", weight = 0.5, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.1,
              highlightOptions = highlightOptions(color = "red", weight = 2,
                                                  bringToFront = TRUE),
              label=as.character(SW_CRA$CRA),
              labelOptions = labelOptions(direction = "auto"))

saveWidget(leafmap, file="./SW_NRCS_CRAs.html", selfcontained = FALSE)


