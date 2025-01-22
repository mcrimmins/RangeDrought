# test code to sort out linear modeling with RPMS/climate raster stacks
# mac 03/30/20
# builds on code from analyzeRPMS.R

library(raster)
library(rgdal)
library(rasterVis)

# set rasteroptions
rasterOptions(progress = 'text')

# load states
states<-getData('GADM', country='USA', level=1)  
us<-getData('GADM', country='USA', level=2)
  counties<-subset(us, NAME_1=="Arizona")
  county<-subset(counties, NAME_2=="Santa Cruz")

# load reprojected RPMS data and gridmet
precip<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_monthly_sum_precip_1979_2019.grd")
rpmsSW<-stack("/scratch/crimmins/USDA_NPP/v2019/AZNM_RPMS_8419_WGS84.grd")

# evaluating time periods ----
# subset
  #dates<-seq(as.Date("1979-01-01","%Y-%m-%d"),as.Date("2019-12-01","%Y-%m-%d"),by="month")
  precip<-precip[[61:492]]
# crop to county
precip<-crop(precip, county)
rpmsSW<-crop(rpmsSW, county)
# sum over monthly periods
#precip <- calc(precip, function(x) movingFun(x, 3,type = "to", sum))
library(cluster)
beginCluster(n=5)
  precip <- clusterR(precip,fun=calc,args=list(fun = function(x) movingFun(x, 3,type = "to", sum)))
endCluster()
  
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
  beginCluster(n = 4)
  for(i in 1:12){
    layers<-seq(i,nlayers(precip)-(12-i),12)
    varStack<-stack(precip[[layers]],rpmsSW)
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
  levelplot(resultsStack, margin=FALSE, par.settings = BuRdTheme, at=seq(-1, 1, 0.1))+
    layer(sp.polygons(county))
  levelplot(which.max(resultsStack), margin=FALSE)+
    layer(sp.polygons(county))
  
# end of   
  
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

# ----


# correlate summer precip with annual production
seasSums<-stack("/scratch/crimmins/gridmet/update_Aug2019/processed/rpmsRes_AZNM_gridmet_seasonal_sum_precip_1979_2019.grd")
# extract summer precip for rpms por
names(seasSums)<-seq(as.Date("1979-01-01","%Y-%m-%d"),as.Date("2019-12-01","%Y-%m-%d"),by="quarter")
#seasPrecip<-seasSums[[seq(23,163,4)]] # summer
seasPrecip<-seasSums[[seq(22,162,4)]] # summer

#mskSeasPrecip = mask(seasPrecip, rpmsSW)
#varStack<-stack( mskSeasPrecip, rpmsSW)
varStack<-stack(seasPrecip,rpmsSW)

# crop to Mohave
varStack<-crop(varStack, county)

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


# linear model of precip/rpms
# scratch
fun_slope=function(x) {lm(x[1:36] ~ x[37:72])$coefficients[2]}
slope <- calc(varStack, fun= function(x) {
  if (sum(is.na(x))>=(length(x)/2)) {
    NA
  } else {
    fun_slope(x)
  }
})

fun_slope<-function(x) {lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)])$coefficients[2];}
fun_slope<-function(x) {m<-lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)]);summary(m)$r.squared}


# fun_slope<-function(x) {lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)])$coefficients[2];
#                         lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)])$coefficients[8];
#                         m<-lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)]);summary(m)$r.squared
#                         }

# fun_slope<-function(x) {m<-lm(x[1:(length(x)/2)] ~ x[((length(x)/2) + 1):length(x)]);
#                           r2<-summary(m)$r.squared;
#                           slope<-summary(m)$coefficients[2];
#                           pval<-summary(m)$coefficients[8];
#                           return(list(r2,slope,pval))
# }

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

