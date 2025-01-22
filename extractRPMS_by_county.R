# extract RPMS data for AZ counties
# MAC 11/10/2020

library(raster)

# set rasteroptions
rasterOptions(progress = 'text')

# load rpms data
rpmsSW<-stack("/scratch/crimmins/USDA_NPP/v2019/AZNM_RPMS_8419_WGS84.grd")
names(rpmsSW)<-seq(1984,2019,1)

# load states
us<-getData('GADM', country='USA', level=2)
#county<-subset(us,NAME_2==countyName)
county<-subset(us, NAME_1=="Arizona")

# county boundaries
countyRas<- rasterize(county, rpmsSW[[1]])

# rpms zonal stats
meanRPMS<-t(zonal(rpmsSW, countyRas, 'mean'))
sdRPMS<-t(zonal(rpmsSW, countyRas, 'sd'))

meanRPMS<-meanRPMS[-1,]
sdRPMS<-sdRPMS[-1,]

colnames(meanRPMS)<-county$NAME_2
colnames(sdRPMS)<-county$NAME_2