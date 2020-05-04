# Linear modeling of different windows of climate variables with RPMS
# MAC 04/01/20
# builds on code from test_analyzeRPMS.R
# CYVERSE VERSION

library(raster)

# set rasteroptions
rasterOptions(progress = 'text')

# or calc SPI/SPEI
# water bal for SPEI  
#precip<-  stack("/scratch/rpmsRes_AZNM_gridmet_3monthly_sum_precip_1984_2019.grd")
#pet<-     stack("/scratch/rpmsRes_AZNM_gridmet_3monthly_sum_pet_1984_2019.grd")
wtrBal<-  stack("/scratch/rpmsRes_AZNM_gridmet_3monthly_precip_minus_pet_1984_2019.grd") 
rpmsSW<-stack("/scratch/AZNM_RPMS_8419_WGS84.grd")

#                           r2<-summary(m)$r.squared;
#                           slope<-summary(m)$coefficients[2];
#                           pval<-summary(m)$coefficients[8];

# lm function
#fun_slope<-function(x) {lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)])$coefficients[2];} # slope coeff
fun_slope<-function(x) { m<-lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)]);summary(m)$r.squared} # r2

# calc with custom function
#resultsStack=list()
beginCluster(n = 7)
for(i in 1:12){
  layers<-seq(i,nlayers(wtrBal)-(12-i),12)
  varStack<-stack(wtrBal[[layers]],rpmsSW)
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
  writeRaster(slope, 
              filename = paste0("./tmpFiles/layer",i,".grd"), overwrite=TRUE)
  print(i)
}
endCluster()
# combine into stack
#temp = stack(resultsStack)
# 3-month labels
dateSeq<-c(paste0(month.abb[11],"-",month.abb[1]),paste0(month.abb[12],"-",month.abb[2]),
           paste0(month.abb[seq(1,12-2,1)],"-",month.abb[seq(1+2,12,1)]))



#writeRaster(resultsStack, 
#            filename = "./rpms/rpms_AZNM_precip3mo_lm_r2_1984_2019.grd", overwrite=TRUE)