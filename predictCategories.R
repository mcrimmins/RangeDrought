# predict rpms/ndvi drought/climate categories using gridmet SPI/SPEI
# need output from analyze_by_LRU_ordReg.R
# MAC 05/01/2020

library(rms)
library(MASS)
library(reshape2)
library(gtools)
library(stringr)
library(tidyr)
library(verification)


# results from analyze_by_LRU_ordReg.R
#load("~/RProjects/RangeDrought/results/AZ_LMU_smNDVI_gridmet_SPI3_SPEI3_USDMcats_results.Rdata")
load("~/RProjects/RangeDrought/results/AZ_LMU_smNDVI_gridmet_SPI3_SPEI3_AbvBloCats_results.Rdata")
#load("~/RProjects/RangeDrought/results/AZ_LMU_rpms_gridmet_SPI3_SPEI3_USDMcats_results.Rdata")

# add month/year to climate dataframes
spiDataFrame$month<-as.numeric(format(spiDataFrame$date, format="%m"))
spiDataFrame$year<-as.numeric(format(spiDataFrame$date, format="%Y"))
speiDataFrame$month<-as.numeric(format(speiDataFrame$date, format="%m"))
speiDataFrame$year<-as.numeric(format(speiDataFrame$date, format="%Y"))

# empty list
droughtCat<-list()
droughtProb<-list()
allDroughtProb<-list()

# swith to target climate dataset
climData<-spiDataFrame # CHANGE TO SPI OR SPEI
# switch to forecast type, needs to match dataset
quantCuts<-c(0,0.33,0.66,1); quantNames<-c("Below","Normal","Above") # for 3 cat
#quantCuts<-c(0,0.03,0.06,0.11,0.21,0.30,1); quantNames<- c("D4","D3","D2","D1","D0","No Drought") # for 3 cat


# loop through all LMU
for(i in 1:nrow(LMU)){
  # wrangle data
  allSeas<-climData[,c(paste0("LRU-",i),"month","year")] 
  allSeas<-dcast(allSeas, year~month, value.var = paste0("LRU-",i))
    colnames(allSeas)<-c("year",paste0(month.abb[11],"_",month.abb[1]),paste0(month.abb[12],"_",month.abb[2]),
                       paste0(month.abb[seq(1,12-2,1)],"_",month.abb[seq(1+2,12,1)]))
  allSeas<-subset(allSeas, year>=1984)
  allSeas$rpms<-rpmsDataFrame[,i]
  # create factors
  #  allSeas$rpms<- quantcut(allSeas$rpms, c(0,0.03,0.06,0.11,0.21,0.30,1)) # USDM percentiles
   allSeas$rpms<- quantcut(allSeas$rpms, quantCuts) 
  # rename levels
  #  levels(allSeas$rpms) <- c("D4","D3","D2","D1","D0","No Drought")
    levels(allSeas$rpms) <- quantNames
  
  # develop prediction model  
  periods<-unlist(str_split(str_remove_all(as.character(LMU@data$form_spi[i]),"[+,]")," "))
    form<-as.formula(paste(c("rpms ~", periods[2:3]),collapse=" + "))
    
    # error catch for model failures ----
    tryCatch({
      # this is the test (just run the call in question)
      catModel<-polr(form, data = allSeas)
    }, warning = function(w) {
      # this is how tryCatch handles a stop event
      message("Model doesn't fit") # message to console
      modelFail <<- TRUE # set a variable to test                 
    }, error = function(e) {
      # this is tryCatch handing an error event
      message("error!")
    }, finally = {
      # this is what to do in the event of an.. event
      # I'm not using this approach
    })
    
    if(exists("modelFail")){
        cat("DON'T PANIC! The model fit failed, moving on...")
        rm(modelFail) # I clear the flag to reuse in case of another failure in loop
        # store results in list for each LMU
        droughtCat[[i]]<-NA
        droughtProb[[i]]<-NA
        #return()
    } else {
    # now it is safe to run the fit as we have handled an error above                  
      catModel<-polr(form, data = allSeas) # build model
    # create custom monthly time series
      newClim<-climData[,c(paste0("LRU-",i),"month","year")] 
      # conditional switch if year is in formula
      if(which(colnames(allSeas)==periods[2])==1){
        # only var2 with year already in frame
        newClim$var1<-climData[,1]
        newClim$var1[which(newClim$month!=(which(colnames(allSeas)==periods[3])-1))]<-NA
        newClim<-fill(newClim, var1)
        colnames(newClim)[4]<-periods[3]
      }else{
        # var1
        newClim$var1<-climData[,1]
        newClim$var1[which(newClim$month!=(which(colnames(allSeas)==periods[2])-1))]<-NA
        newClim<-fill(newClim, var1)
        colnames(newClim)[4]<-periods[2]
        # var2
        newClim$var1<-climData[,1]
        newClim$var1[which(newClim$month!=(which(colnames(allSeas)==periods[3])-1))]<-NA
        newClim<-fill(newClim, var1)
        colnames(newClim)[5]<-periods[3]
      }
    
      # make predictions  
      predCats<-predict(catModel,newClim,type="probs")
          z<-apply(as.data.frame(predCats),1,which.max)
          z<-max.col(predCats)
          newClim$catNames<-colnames(predCats)[z]
          newClim$probs<-apply(predCats,1,max)
          
      # store results in list for each LMU
          droughtCat[[i]]<-newClim$catNames
          droughtProb[[i]]<-newClim$probs
          allDroughtProb[[i]]<-predCats
    }
    print(i)
}

# turn lists into data frames
droughtCatdf =as.data.frame(do.call(cbind, droughtCat))
  colnames(droughtCatdf)<-paste0("LRU-",seq(1,ncol(droughtCatdf),1))
  droughtCatdf$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")

droughtProbdf =as.data.frame(do.call(cbind, droughtProb))
  colnames(droughtProbdf)<-paste0("LRU-",seq(1,ncol(droughtProbdf),1))
  droughtProbdf$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
  
# check drought category frequencies
table(unlist(droughtCatdf[,1:91]))
  
# verification stats
  qfunc<-function(x)  { qcut<-quantcut(x, quantCuts)
                        levels(qcut) <- seq(1,length(quantCuts)-1,1)
                        qcut<-as.numeric(qcut)
                        return(qcut)
  }
  
rpmsCats<- apply((rpmsDataFrame[,1:ncol(rpmsDataFrame)-1]), 2, FUN = qfunc)

verifList<-list()
for(i in 1:ncol(rpmsCats)){
  temp<- as.data.frame(allDroughtProb[[i]])
  if(nrow(temp)!=0){  
    temp$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
      temp<- temp[seq(0, nrow(temp), 12), ]
      temp<-temp[((nrow(temp)-nrow(rpmsCats))+1):nrow(temp),]
      verifList[[i]]<-rps( obs = rpmsCats[,i], pred = base::as.matrix(temp[,1:(length(quantCuts)-1)]))
  }else{
    verifList[[i]]<-verifList[[i-1]]
    verifList[[i]][1:3]<-NA
  }
}

verifList <-as.data.frame(matrix(unlist(verifList), nrow=length(verifList), ncol=3, byrow=TRUE))
  colnames(verifList)<-c("rps","rpss","rps.clim")
  mean(verifList$rps, na.rm = TRUE); median(verifList$rps, na.rm = TRUE)
  mean(verifList$rpss, na.rm = TRUE); median(verifList$rpss, na.rm = TRUE)
  hist(verifList$rps)
# map rps/rpss
  verifLMU<-LMU
  verifLMU@data<-cbind(verifLMU@data,verifList)
  
  #pal <- RColorBrewer::brewer.pal(n=9, 'OrRd')
  #library(rasterVis)
  pal <- RColorBrewer::brewer.pal(n=9, 'PRGn')
  spplot(verifLMU, c("rps"), 
         main=list(label=paste0("Rank Probability Score - NDVI/SPI/3 cat (",round(mean(verifList$rps, na.rm = TRUE),2),")"),
                   cex=1),scales=list(draw = TRUE),
         at=seq(0, 0.9, 0.1),col.regions = pal)
 
  spplot(verifLMU, c("rpss"), 
         main=list(label=paste0("Rank Probability Skill Score - NDVI/SPI/3 Cat (",round(mean(verifList$rpss, na.rm = TRUE),2),")"),
                   cex=1),scales=list(draw = TRUE),
         at=seq(-1, 1, 0.2),col.regions = rev(brewer.pal(n = 10, 'RdBu')))
  
   
  
    
# map category results
temp<-as.data.frame(t(droughtCatdf[,1:91]))
temp[] <- lapply(temp, factor, 
               levels=c("Above","Below","Normal"), 
               labels = c("Above","Below","Normal"))

mapCats<-LMU
#colnames(temp)<-paste0("X",droughtCatdf[,92])
  mapCats@data<-cbind.data.frame(LMU@data,temp)
  dateNames<-format(seq(as.Date("1979.01.01", "%Y.%m.%d"),
                        as.Date("2019.12.01", "%Y.%m.%d"), by="months"),"%m.%d.%Y")
  colnames(mapCats@data)[(ncol(LMU@data)+1):(ncol(temp)+(ncol(LMU@data)))]<-paste0("X",dateNames)
  
  # plot drought cats
  dates<-seq(as.Date("1979.01.01", "%Y.%m.%d"),
           as.Date("2019.12.01", "%Y.%m.%d"), by="months")
  dateIdx<-which(dates>=as.Date("2000-01-01","%Y-%m-%d") & dates<=as.Date("2000-12-01","%Y-%m-%d"))+(ncol(LMU@data)+1)
  library(RColorBrewer)
  my.palette <- brewer.pal(n = 3, name = "RdYlGn")
  my.palette<-c("#91CF60","#FC8D59","#FFFFBF")
    spplot(mapCats, dateIdx, 
         main=list(label="Drought Categories",cex=1),scales=list(draw = TRUE), col.regions=my.palette)
  