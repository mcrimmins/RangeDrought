# Correlate different windows of climate variables with RPMS
# MAC 04/01/20
# builds on code from test_analyzeRPMS.R

library(raster)
library(rgdal)
library(rasterVis)
library(cluster)

# set rasteroptions
rasterOptions(progress = 'text')

# load states
states<-getData('GADM', country='USA', level=1)  
#us<-getData('GADM', country='USA', level=2)
#counties<-subset(us, NAME_1=="Arizona")
#county<-subset(counties, NAME_2=="Santa Cruz")

# load reprojected RPMS data and gridmet
pet<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_monthly_sum_pet_1979_2019.grd")
precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_monthly_sum_precip_1979_2019.grd")
rpmsSW<-stack("/scratch/crimmins/USDA_NPP/v2019/AZNM_RPMS_8419_WGS84.grd")

# evaluating time periods ----
# subset to RPMS data period
#dates<-seq(as.Date("1979-01-01","%Y-%m-%d"),as.Date("2019-12-01","%Y-%m-%d"),by="month")
#precip<-precip[[61:492]]
pet<-pet[[61:492]]
# crop to county
#precip<-crop(precip, county)
#rpmsSW<-crop(rpmsSW, county)
# sum over monthly periods
library(cluster)
  beginCluster(n=5)
    pet <- clusterR(pet,fun=calc,filename="./temp/pet_moving_sum.grd",
                       args=list(fun = function(x) movingFun(x, 3,type = "to", sum)))
  endCluster()
  writeRaster(pet, 
              filename = "/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_3monthly_sum_pet_1984_2019.grd", overwrite=TRUE)
  
# or calc SPI/SPEI
# water bal for SPEI  
  precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_3monthly_sum_precip_1984_2019.grd")
  pet<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_3monthly_sum_pet_1984_2019.grd")
  wtrBal<- overlay(precip,
                   pet,
                   fun=function(r1, r2){return(r1-r2)})
  writeRaster(wtrBal, 
              filename = "/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_3monthly_precip_minus_pet_1984_2019.grd", overwrite=TRUE)
  
   
# function
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
library(cluster)
# calc with custom function
resultsStack=list()
beginCluster(n = 6)
for(i in 7:12){
  layers<-seq(i,nlayers(wtrBal)-(12-i),12)
  varStack<-stack(wtrBal[[layers]],rpmsSW)
  corlyrs <- clusterR(x         = varStack,
                      fun       = calc,
                      args      = list(fun = function(x) {
                        if (sum(is.na(x))>=(length(x)/2)) {
                          NA_real_
                        } else {
                          corvec(x)
                        }
                      }),
                      export    = 'corvec'
  )
  resultsStack[[i]]<-corlyrs
  print(i)
}
endCluster()
# combine into stack
resultsStack = stack(resultsStack)
# 3-month labels
dateSeq<-c(paste0(month.abb[11],"-",month.abb[1]),paste0(month.abb[12],"-",month.abb[2]),
          paste0(month.abb[seq(1,12-2,1)],"-",month.abb[seq(1+2,12,1)]))
# 6-month labels
# dateSeq<-c(paste0(month.abb[8:12],"-",month.abb[1:5]),
#           paste0(month.abb[seq(1,12-5,1)],"-",month.abb[seq(1+5,12,1)]))
 names(resultsStack)<-dateSeq
writeRaster(resultsStack, 
            filename = "/home/crimmins/RProjects/RangeDrought/results/rpms_AZNM_precip_minus_pet_3mo_pearson_1984_2019.grd", overwrite=TRUE)

# plot  
resultsStack<-stack("/home/crimmins/RProjects/RangeDrought/results/rpms_AZNM_precip6mo_pearson_1984_2019.grd")
  bwplot(resultsStack, main="RPMS/3mo Precip-PET Pearson r (1984-2019)")
  plot(cellStats(resultsStack, stat=mean), type="b")
p0<-levelplot(resultsStack, margin=FALSE, par.settings = BuRdTheme, at=seq(-1, 1, 0.1), main="RPMS/3mo Precip-PET Pearson r (1984-2019)")+
  layer(sp.polygons(states))
# plot to png
 png("/home/crimmins/RProjects/RangeDrought/figs/RPMS_3MoPrecip_minus_PET_Corrs.png", width = 10, height = 8, units = "in", res = 300L)
 print(p0, newpage = FALSE)
 dev.off()

p0<-levelplot(which.max(resultsStack), margin=FALSE, main="3mo period with highest r (RPMS-3mo P-PET corr, 1984-2019")+
  layer(sp.polygons(states))
# plot to png
png("/home/crimmins/RProjects/RangeDrought/figs/RPMS_3MoPrecip_minus_PET_Highest_Corrs.png", width = 10, height = 8, units = "in", res = 300L)
print(p0, newpage = FALSE)
dev.off()

# compare climate corrs
# 3-month labels
dateSeq<-c(paste0(month.abb[11],"-",month.abb[1]),paste0(month.abb[12],"-",month.abb[2]),
           paste0(month.abb[seq(1,12-2,1)],"-",month.abb[seq(1+2,12,1)]))

precipCorr<-stack("/home/crimmins/RProjects/RangeDrought/results/rpms_AZNM_precip3mo_pearson_1984_2019.grd")
petCorr<-stack("/home/crimmins/RProjects/RangeDrought/results/rpms_AZNM_precip_minus_pet_3mo_pearson_1984_2019.grd")
#  bwplot(precipCorr-petCorr, main="RPMS/3mo Precip-PET Pearson r (1984-2019)")
diff<-precipCorr-petCorr
names(diff)<-dateSeq

p0<-levelplot(diff, margin=FALSE, par.settings = BuRdTheme, at=seq(-0.25, 0.25, 0.05), main="Precip-(Precip-PET) Pearson r (1984-2019)")+
  layer(sp.polygons(states))
# plot to png
png("/home/crimmins/RProjects/RangeDrought/figs/RPMS_3Mo_diff_P_PPET_Corrs.png", width = 10, height = 8, units = "in", res = 300L)
print(p0, newpage = FALSE)
dev.off()  

# compare vals
temp<-(stack(max(precipCorr),max(petCorr)))
  names(temp)<-c("Precip","P-PET")
densityplot(temp, main="Max 3mo Corrs")

# compare 3 to 6 month corr density
temp3<-stack("/home/crimmins/RProjects/RangeDrought/results/rpms_AZNM_precip3mo_pearson_1984_2019.grd")
temp6<-stack("/home/crimmins/RProjects/RangeDrought/results/rpms_AZNM_precip6mo_pearson_1984_2019.grd")
  p0<-levelplot(temp3[[c(12,3,6,9)]], layout=c(4, 1), margin=FALSE, par.settings = BuRdTheme, at=seq(-1, 1, 0.05), main="3-mo Precip Pearson r (1984-2019)")+
    layer(sp.polygons(states))
  p1<-levelplot(temp6[[c(3,9)]], layout=c(2, 1), margin=FALSE, par.settings = BuRdTheme, at=seq(-1, 1, 0.05), main="6-mo Precip Pearson r (1984-2019)")+
    layer(sp.polygons(states))
  library(gridExtra)
  grid.arrange(p0, p1, nrow=2)
temp<-stack(temp3[[c(12,3,6,9)]],temp6[[c(3,9)]])
  densityplot(temp, main="Density of 3 and 6mo Precip/RPMS Pearson r")
  # just AZ
  densityplot(crop(temp,extent(subset(states, NAME_1=="Arizona"))),main="Density of 3 and 6mo Precip/RPMS Pearson r - AZ")
  # all months
  densityplot(stack(temp3,temp6),main="Density of 3 and 6mo Precip/RPMS Pearson r")

# leaflet ----
#maxCorrMo<-which.max(resultsStack)
maxCorrMo<-max(petCorr)

library(rgdal)
library(leaflet)
library(leafem)
library(htmlwidgets)

# ---- load NRCS CRAs for SW
AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
# combine
SW_CRA<- rbind(AZ_CRA,NM_CRA)

pal=colorNumeric(c("blue","grey", "red"), values(maxCorrMo),
                 na.color = "transparent")
leafMap<-leaflet() %>% addTiles() %>%
  addRasterImage(maxCorrMo, colors = pal, maxBytes=20000000,opacity = 0.8, layerId = "maxCorrMo") %>%
  addMouseCoordinates() %>%
  addImageQuery(maxCorrMo, type="mousemove", layerId = "maxCorrMo", digits = 2, prefix = "") %>%
  addLegend(pal = pal, values = values(maxCorrMo),
            title="Highest Corr")%>%
  addPolygons(data=SW_CRA, color = "#444444", weight = 0.5, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.1,
              highlightOptions = highlightOptions(color = "red", weight = 2,
                                                  bringToFront = TRUE),
              label=as.character(SW_CRA$CRA),
              labelOptions = labelOptions(direction = "auto"))
saveWidget(leafMap, file="/home/crimmins/RProjects/RangeDrought/maps/highest_corr_rpms_wtrbal.html", selfcontained = TRUE)
# ----