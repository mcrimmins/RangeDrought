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

# or calc SPI/SPEI
# water bal for SPEI  
precip<-  stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_3monthly_sum_precip_1984_2019.grd")
pet<-     stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_3monthly_sum_pet_1984_2019.grd")
wtrBal<-  stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_3monthly_precip_minus_pet_1984_2019.grd") 
rpmsSW<-stack("/scratch/crimmins/USDA_NPP/v2019/AZNM_RPMS_8419_WGS84.grd")

  
# lm function
fun_slope<-function(x) {lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)])$coefficients[2];}
#fun_slope<-function(x) { m<-lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)]);summary(m)$r.squared}

# calc with custom function
resultsStack=list()
beginCluster(n = 6)
for(i in 10:12){
  layers<-seq(i,nlayers(precip)-(12-i),12)
  varStack<-stack(precip[[layers]],rpmsSW)
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
  resultsStack[[i]]<-slope
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
            filename = "/home/crimmins/RProjects/RangeDrought/results/rpms_AZNM_precip3mo_lm_slope_coeff_1984_2019.grd", overwrite=TRUE)

# plot  
#resultsStack<-stack("/home/crimmins/RProjects/RangeDrought/results/rpms_AZNM_precip6mo_pearson_1984_2019.grd")
bwplot(resultsStack, main="RPMS/3mo Precip Slope Coeff (1984-2019)")
plot(cellStats(resultsStack, stat=mean), type="b", main="Mean slope coeff PRMS~3moPrecip")
p0<-levelplot(resultsStack, margin=FALSE, par.settings = BuRdTheme, at=seq(-0.5, 0.5, 0.01), main="Mean slope coeff PRMS~3moPrecip (1984-2019)")+
  layer(sp.polygons(states))
# plot to png
png("/home/crimmins/RProjects/RangeDrought/figs/RPMS_3MoPrecip_slope_coeff.png", width = 10, height = 8, units = "in", res = 300L)
print(p0, newpage = FALSE)
dev.off()



