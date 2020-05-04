# process VHI into rasters
# MAC 02/26/19

library(raster)
library(stringr)
library(lubridate)
library(plyr)
library(grid)
library(cowplot)

# set rasteroptions
rasterOptions(progress = 'text')

# process with list of filenames
fileNames <- dir("/scratch/crimmins/vhi/smNDVI", "*.tif", full.names = TRUE)
#fileNames <- dir("/scratch/crimmins/vhi/VHItiff", "*.tif", full.names = TRUE)

# full dates list for missing layer check ----
AllfileNames<-list()
allYears<-list()
allWeeks<-list()
l<-1

# create full file list
for (j in 1982:2019) {
  for (i in 1:52) {
    AllfileNames[l]<-paste0(j,str_pad(i, 3, pad = "0"))
    allYears[l]<-j
    allWeeks[l]<-i
    l<-l+1
  }
}
# full dates
fullDates<-as.data.frame(cbind(unlist(allWeeks),unlist(allYears)))
fullDates$date<-as.Date(paste(fullDates$V2, fullDates$V1, 1, sep="-"), "%Y-%U-%u")

# insert empty layers for missing periods...addLayer?  
fileString<-as.data.frame(t(as.data.frame(strsplit(substr(fileNames,30,63),"[.]"))))
rownames(fileString) <- c()
colnames(fileString)<-c("prefix","res","sat","comPeriod","yearWeek","prod","type","fileExt")
fileString$year<-as.numeric(substr(fileString$yearWeek, 2,5))   
fileString$week<-as.numeric(substr(fileString$yearWeek, 6,8))
fileString$date<-as.Date(paste(fileString$year, fileString$week, 1, sep="-"), "%Y-%U-%u")
# merge with full dates
mergedList<-merge(fileString,fullDates,by="date", all.y=TRUE)
# find missing days
missingWeeks<-mergedList[which(is.na(mergedList$prefix)==TRUE),c("V1","V2")]
indexNA<-which(is.na(mergedList$prefix)==TRUE)
indexNotNA<-which(is.na(mergedList$prefix)==FALSE)
  orderLyrs<-as.integer(c(indexNotNA,indexNA))
# ----

# for Range Drought project use SW_CRA extent
#> extent(SW_CRA)
#class      : Extent 
#xmin       : -116.2267 
#xmax       : -99.21573 
#ymin       : 28.9193 
#ymax       : 42.76954 
# extent(-116.2267,-99.21573,28.9193,42.76954)

# read layers into stack seq(1, length(fileNames), by=1)
l<-1
  smNDVI <- stack()
  for (i in seq(1, length(fileNames), by=1)) {
    tempRast<-raster(fileNames[i])
    e <- extent(-116.2267,-99.21573,28.9193,42.76954) # WESTUS extent(-125, -100, 25, 49)
    tempRast <- crop(tempRast, e)
    tempRast[tempRast<=-999]<-NA # or set to 0 for NDVI?
    smNDVI <- stack(smNDVI, tempRast)
    print(names(tempRast))
    l<-l+1
  }
# fix values
#smNDVI[smNDVI < 0] <- NA

# add NA layers to make complete dates raster stack
  layerNA<-smNDVI[[1]]
    layerNA[layerNA>=0]<-NA
  smNDVI<-raster::subset(stack(smNDVI, stack(replicate(length(indexNA), layerNA))), order(orderLyrs))
  names(smNDVI)<-mergedList$date
# write Raster to file
writeRaster(smNDVI,filename='/scratch/crimmins/vhi/processed/VHI_complete_1982-2019_SWUS.grd', overwrite=TRUE)




# ---- EXTRACT time series  
# things to try: mask out forest areas, evaluate min/max, stdev of vals in county, do compare 1/2/4 weeks per month

# CHANGE NDVI PRODUCT TYPE
#smNDVI<-stack('/scratch/crimmins/vhi/processed/smNDVI_1982-2018_WESTUS.grd')
smNDVI<-stack('/scratch/crimmins/vhi/processed/VHI_1982-2018_WESTUS.grd')

# counties
c("San Luis Obispo","Santa Barbara","Humboldt","Rich","Garfield","Kane County","Weld","Logan",
  "Morgan","Luna","Hidalgo","Sierra","Grant","Gila","Cherry","Grant")



# county time series
# labels 
countyName<-"Garfield"
# get boundary
us<-getData('GADM', country='USA', level=2)
#county<-subset(us,NAME_2==countyName)
county<-subset(us,NAME_2==countyName & NAME_1=="Utah")

# extract time series
#test <- cellStats(mask(smNDVI,county,inverse=FALSE),stat = "mean" ,na.rm=TRUE)
tsNDVI <- t(extract(smNDVI, county, fun=mean, df=TRUE, na.rm=TRUE))
# add file info, dates from file name
fileString<-as.data.frame(t(as.data.frame(strsplit((rownames(tsNDVI)),"[.]"))))
  rownames(fileString) <- c()
  colnames(fileString)<-c("prefix","res","sat","comPeriod","yearWeek","prod","type")
  fileString$year<-as.numeric(substr(fileString$yearWeek, 2,5))   
  fileString$week<-as.numeric(substr(fileString$yearWeek, 6,8))
  fileString$date<-as.Date(paste(fileString$year, fileString$week, 1, sep="-"), "%Y-%U-%u")
# combine
tsNDVI<-cbind(tsNDVI,fileString)  
  tsNDVI<-tsNDVI[-1,]
  rownames(tsNDVI) <- c()
tsNDVI$month<-as.numeric(format(tsNDVI$date, "%m"))
tsNDVI$weekYearCode<-paste0(tsNDVI$week,"-",tsNDVI$year)
colnames(tsNDVI)[1]<-"NDVI"
#plot(tsNDVI$date,tsNDVI$NDVI, type="l")

# get yearly ranks
temp<-tsNDVI[,c(1,9)]
yearlyNDVI <- ddply(temp, .(year), summarize, yearAvgNDVI = mean(NDVI, na.rm=TRUE))

# ---- NDVI Time series ----

# plot NDVI time series
# pTSndvi<-ggplot(tsNDVI, aes(date, NDVI))+
#   geom_line()+
#   scale_x_date(date_labels = "%b %Y", date_breaks="2 years")+
#   #theme_bw()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   geom_hline(yintercept=mean(tsNDVI$NDVI, na.rm=TRUE), linetype="dashed", color = "gray")+
#   #geom_hline(yintercept=40, color='burlywood4', size=0.5)+
#   #geom_hline(yintercept=60, color='darkgreen', size=0.5)+
#   #ylim(0,1)+
#   ggtitle(paste0("Monthly Normalized Difference Vegetation Index (NDVI,1981-2018): ",countyName," County"))+
#   labs(x = "Date",y="NDVI")

# plot VHI time series
pTSvhi<-ggplot(tsNDVI, aes(date, NDVI))+
  geom_line()+
  scale_x_date(date_labels = "%b %Y", date_breaks="2 years")+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept=mean(tsNDVI$NDVI, na.rm=TRUE), linetype="dashed", color = "gray")+
  #geom_hline(yintercept=40, color='burlywood4', size=0.5)+
  #geom_hline(yintercept=60, color='darkgreen', size=0.5)+
  #ylim(0,100)+
  ggtitle(paste0("Monthly Vegetation Health Index (VHI,1981-2018): ",countyName," County"))+
  labs(x = "Date",y="VHI")



### --- SEASONAL PLOTS ---

# thinned out
thinNDVI<-tsNDVI[,c(12,9,1)]

# get monthly percentiles
NDVIQuant<- ddply(thinNDVI,.(month),summarise,
                 q05 = quantile(NDVI,0.05,na.rm='TRUE'),
                 q25 = quantile(NDVI,0.25,na.rm='TRUE'),
                 q10 = quantile(NDVI,0.10,na.rm='TRUE'),
                 q50 = quantile(NDVI,0.50,na.rm='TRUE'),
                 q90 = quantile(NDVI,0.90,na.rm='TRUE'),
                 q75 = quantile(NDVI,0.75,na.rm='TRUE'),
                 q95 = quantile(NDVI,0.95,na.rm='TRUE'),
                 min = min(NDVI,na.rm='TRUE'),
                 max = max(NDVI,na.rm='TRUE'),
                 avg = mean(NDVI,na.rm='TRUE'))
NDVIQuant$date <- as.Date(strptime(paste("2012", NDVIQuant$month, "1"), format="%Y %m %d"))

# plot monthly NDVI - FULL PERCENTILES
p1<-ggplot(NDVIQuant, aes(x=date, y=avg)) +
  geom_line(size=.5) + 
  geom_ribbon(aes(ymin=min, ymax=max, x=date, fill = "min/max"), alpha = 0.7) +
  geom_ribbon(aes(ymin=q05, ymax=q95, x=date, fill = "5th-95th"), alpha = 0.7) +
  geom_ribbon(aes(ymin=q25, ymax=q75, x=date, fill = "25th-75th"), alpha = 0.7) +
  geom_ribbon(aes(ymin=q50, ymax=q50, x=date, fill = "Median"), alpha = 0.7) +
  scale_fill_manual(name='Percentiles', values=c("chartreuse","chartreuse2","black","chartreuse4"))+
  geom_line(size=.5)+
  scale_x_date(date_labels = "%B")+
  ggtitle(paste0("Monthly NDVI (1981-2018): ",countyName," County"))+
  labs(x = "Month",y="NDVI-Greenness")

# MIN NDVI months
#tempMin<-tsNDVI[which(tsNDVI$year==yearlyNDVI[which.min(yearlyNDVI$yearAvgNDVI),"year"]),] # auto select
tempMin<-tsNDVI[which(tsNDVI$year==2002),] # manual select
tempMin$date <- as.Date(strptime(paste("2012", tempMin$month, "1"), format="%Y %m %d"))
ggplot(tempMin, aes(x=date, y=NDVI)) +
         geom_line(size=.5)
# plot min and median stats
ggplot()+
  geom_line(data=NDVIQuant, aes(x=date, y=avg,color="darkgrey"))+
  geom_ribbon(data=NDVIQuant, aes(ymin=q05, ymax=q95, x=date, fill = "5th-95th"), alpha = 0.3) +
  scale_fill_manual(name='Percentiles', values=c("grey"))+
  geom_line(data=tempMin, aes(x=date, y=NDVI, color="darkgreen"))+
  scale_x_date(date_labels = "%B")+
  ggtitle(paste0("Monthly NDVI (1981-2018): ",countyName," County"))+
  labs(x = "Month",y="NDVI-Greenness")+
  #scale_color_manual(name = "NDVI", labels = c(paste0(yearlyNDVI[which.min(yearlyNDVI$yearAvgNDVI),"year"]),"average"),
  #                   values=c("darkgreen","darkgrey"))
  scale_color_manual(name = "NDVI", labels = c("2002","average"),
                   values=c("darkgreen","darkgrey"))


# MAX NDVI months
tempMax<-tsNDVI[which(tsNDVI$year==yearlyNDVI[which.max(yearlyNDVI$yearAvgNDVI),"year"]),]
tempMax$date <- as.Date(strptime(paste("2012", tempMax$month, "1"), format="%Y %m %d"))
ggplot(tempMax, aes(x=date, y=NDVI)) +
  geom_line(size=.5)
# plot min and median stats
ggplot()+
  geom_line(data=NDVIQuant, aes(x=date, y=avg,color="darkgrey"))+
  geom_ribbon(data=NDVIQuant, aes(ymin=q05, ymax=q95, x=date, fill = "5th-95th"), alpha = 0.3) +
  scale_fill_manual(name='Percentiles', values=c("grey"))+
  geom_line(data=tempMax, aes(x=date, y=NDVI, color="darkgreen"))+
  scale_x_date(date_labels = "%B")+
  ggtitle(paste0("Monthly NDVI (1981-2018): ",countyName," County"))+
  labs(x = "Month",y="NDVI-Greenness")+
  scale_color_manual(name = "NDVI", labels = c(paste0(yearlyNDVI[which.max(yearlyNDVI$yearAvgNDVI),"year"]),"average"),
                     values=c("darkgreen","darkgrey"))

# do for Precip
# load raw gridded climate data
prec<-stack("/scratch/crimmins/livneh/processed/WESTmonthlyLivneh_prec_1915_2015.grd")
# dates 1895-2017 PRISM data
dates=seq(as.Date("1915-01-01"), as.Date("2015-12-01"), by="month")
prec<-prec[[which(dates=="1981-01-01"):which(dates=="2015-12-01")]] 
# extract prec
tsPREC <- t(extract(prec, county, fun=mean, df=TRUE, na.rm=TRUE))
  tsPREC<-as.data.frame(tsPREC[c(-1),])
dates=as.data.frame(seq(as.Date("1981-01-01"), as.Date("2015-12-01"), by="month"))
tsPREC<-cbind(dates,tsPREC)
  rownames(tsPREC)<-c()
  colnames(tsPREC)<-c("date","prec")
  tsPREC$month<-as.numeric(format(tsPREC$date, "%m"))
  tsPREC$year<-as.numeric(format(tsPREC$date, "%Y"))
  
# get yearly ranks
  temp<-tsPREC[,c(2,4)]
  yearlyPrec <- ddply(temp, .(year), summarize, yearSumPrecip = sum(prec, na.rm=TRUE))  
  
  # get monthly percentiles
precQuant<- ddply(tsPREC,.(month),summarise,
                    q05 = quantile(prec,0.05,na.rm='TRUE'),
                    q25 = quantile(prec,0.25,na.rm='TRUE'),
                    q10 = quantile(prec,0.10,na.rm='TRUE'),
                    q50 = quantile(prec,0.50,na.rm='TRUE'),
                    q90 = quantile(prec,0.90,na.rm='TRUE'),
                    q75 = quantile(prec,0.75,na.rm='TRUE'),
                    q95 = quantile(prec,0.95,na.rm='TRUE'),
                    min = min(prec,na.rm='TRUE'),
                    max = max(prec,na.rm='TRUE'),
                    avg = mean(prec,na.rm='TRUE'))
  precQuant$date <- as.Date(strptime(paste("2012", precQuant$month, "1"), format="%Y %m %d"))

p2<-ggplot(precQuant, aes(x=date, y=avg)) +
    geom_line(size=.5) + 
    geom_ribbon(aes(ymin=min, ymax=max, x=date, fill = "min/max"), alpha = 0.7) +
    geom_ribbon(aes(ymin=q05, ymax=q95, x=date, fill = "5th-95th"), alpha = 0.7) +
    geom_ribbon(aes(ymin=q25, ymax=q75, x=date, fill = "25th-75th"), alpha = 0.7) +
    geom_ribbon(aes(ymin=q50, ymax=q50, x=date, fill = "Median"), alpha = 0.7) +
    scale_fill_manual(name='Percentiles', values=c("darkslategray1","darkslategray2","black","darkslategray3"))+
    geom_line(size=.5)+
    scale_x_date(date_labels = "%B")+
    ggtitle(paste0("Monthly Precipitation (1981-2015): ",countyName," County"))+
    labs(x = "Month",y="Precip (inches)")

# MIN PRECIP months
tempMin<-tsPREC[which(tsPREC$year==yearlyPrec[which.min(yearlyPrec$yearSumPrecip),"year"]),]
tempMin$date <- as.Date(strptime(paste("2012", tempMin$month, "1"), format="%Y %m %d"))
ggplot(tempMin, aes(x=date, y=prec)) +
  geom_line(size=.5)
# plot min and median stats
ggplot()+
  geom_line(data=precQuant, aes(x=date, y=avg,color="darkgrey"))+
  geom_ribbon(data=precQuant, aes(ymin=q05, ymax=q95, x=date, fill = "5th-95th"), alpha = 0.3) +
  scale_fill_manual(name='Percentiles', values=c("grey"))+
  geom_line(data=tempMin, aes(x=date, y=prec, color="darkgreen"))+
  scale_x_date(date_labels = "%B")+
  ggtitle(paste0("Monthly Precip (1981-2015): ",countyName," County"))+
  labs(x = "Month",y="Precip (in)")+
  scale_color_manual(name = "Precip", labels = c(paste0(yearlyPrec[which.min(yearlyPrec$yearSumPrecip),"year"]),"average"),
                     values=c("darkblue","darkgrey"))

# MAX PRECIP months
tempMax<-tsPREC[which(tsPREC$year==yearlyPrec[which.max(yearlyPrec$yearSumPrecip),"year"]),]
tempMax$date <- as.Date(strptime(paste("2012", tempMax$month, "1"), format="%Y %m %d"))
ggplot(tempMax, aes(x=date, y=prec)) +
  geom_line(size=.5)
# plot min and median stats
ggplot()+
  geom_line(data=precQuant, aes(x=date, y=avg,color="darkgrey"))+
  geom_ribbon(data=precQuant, aes(ymin=q05, ymax=q95, x=date, fill = "5th-95th"), alpha = 0.3) +
  scale_fill_manual(name='Percentiles', values=c("grey"))+
  geom_line(data=tempMax, aes(x=date, y=prec, color="darkgreen"))+
  scale_x_date(date_labels = "%B")+
  ggtitle(paste0("Monthly Precip (1981-2015): ",countyName," County"))+
  labs(x = "Month",y="Precip (in)")+
  scale_color_manual(name = "Precip", labels = c(paste0(yearlyPrec[which.max(yearlyPrec$yearSumPrecip),"year"]),"average"),
                     values=c("darkblue","darkgrey"))

# precip monthly time series
ggplot(tsPREC, aes(date, prec))+
  geom_bar(stat = "identity")+
  scale_x_date(date_labels = "%b %Y", date_breaks="2 years")+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept=mean(tsPREC$prec, na.rm=TRUE), linetype="dashed", color = "gray")+
  #geom_hline(yintercept=40, color='burlywood4', size=0.5)+
  #geom_hline(yintercept=60, color='darkgreen', size=0.5)+
  #ylim(0,1)+
  ggtitle(paste0("Monthly Precipitation (1981-2015): ",countyName," County"))+
  labs(x = "Date",y="precip (in)")

# plot yearly precip
pPrecip<-ggplot(yearlyPrec, aes(year, yearSumPrecip))+
  geom_bar(stat = "identity", fill="darkblue")+
  geom_hline(yintercept=mean(yearlyPrec$yearSumPrecip, na.rm=TRUE), linetype="dashed", color = "gray")+
  ggtitle(paste0("Total Annual Precipitation (1981-2015): ",countyName," County"))+
  labs(x = "Date",y="precip (in)")


# ---- PRINT OUT SEASONAL MULTIPLOTS ----
   
pAll<-plot_grid(p1, p2, nrow = 2, align = "v")

#png(paste0("./vhi/seasFigs/",countyName,"_SeasonalClim.png"), width = 7, height = 5, units = "in", res = 300L)
png(paste0("./vhi/seasFigs/SeasonalClim.png"), width = 7, height = 5, units = "in", res = 300L)
#grid.newpage()
print(pAll, newpage = FALSE)
# add footer
## now add the text 
grid.text(paste0("NDVI: https://www.star.nesdis.noaa.gov/smcd/emb/vci/VH/vh_ftp.php 
                 Precip: http://ciresgroups.colorado.edu/livneh/data"),
          x=.85, y=.03, 
          gp = gpar(col=1, 
                    fontfamily="Arial", cex=0.4)) 
dev.off()

# ----

# ---- PRINT OUT TIME SERIES MULTIPLOTS ----

pAllts<-plot_grid(pPrecip, pTSvhi, nrow = 2, align = "v")

#png(paste0("./vhi/seasFigs/",countyName,"_ClimTS.png"), width = 7, height = 5, units = "in", res = 300L)
png(paste0("./vhi/seasFigs/ClimTS.png"), width = 7, height = 5, units = "in", res = 300L)
#grid.newpage()
print(pAllts, newpage = FALSE)
# add footer
## now add the text 
grid.text(paste0("VHI: https://www.star.nesdis.noaa.gov/smcd/emb/vci/VH/vh_ftp.php 
                 Precip: http://ciresgroups.colorado.edu/livneh/data"),
          x=.85, y=.03, 
          gp = gpar(col=1, 
                    fontfamily="Arial", cex=0.4)) 
dev.off()

# ----
  
# TRYING TO CREATE FULL WEEK LIST TO ID MISSING WEEKS
# # create full week list
# l<-1
# allYears<-list()
# allWeeks<-list()
# for (j in 1982:2018) {
#   for (i in 1:52) {
#     allYears[l]<-j
#     allWeeks[l]<-i
#     l<-l+1
#   }
# }
# # full dates
# fullDates<-as.data.frame(cbind(unlist(allWeeks),unlist(allYears)))
# fullDates$date<-as.Date(paste(fullDates$V2, fullDates$V1, 1, sep="-"), "%Y-%U-%u")
#   fullDates$weekYearCode<-paste0(fullDates$V1,"-",fullDates$V2)
#   #fullDates<-as.data.frame(fullDates[seq(1, nrow(fullDates), by=4),])
# 
# # # get full date sequence  
# # fullDates<-as.data.frame(seq(min(tsNDVI$date),max(tsNDVI$date),by="week"))
# #   colnames(fullDates)<-"date"
# # # match with sequence in stack, every 2nd or 4th
# # fullDates<-as.data.frame(fullDates[seq(1, length(fileNames), by=4),])
# #   colnames(fullDates)<-"date"
# # merge with full dates
# tsNDVI<-merge(tsNDVI,fullDates,by="weekYearCode", all.y=TRUE)
