# predict rpms/ndvi drought/climate categories using gridmet SPI/SPEI
# need output from analyze_by_LRU_ordReg.R
# MAC 05/01/2020
# use with analyze_by_LRU_ordReg2.R, deals with variable number of coeffs, no year as predictor

library(rms)
library(MASS)
library(reshape2)
library(gtools)
library(stringr)
library(tidyr)
library(verification)
library(rgdal)
library(pracma)
library(RColorBrewer)
library(rasterVis)

# replace empty with NA in list
nullToNA <- function(x) {
  x[sapply(x, is.null)] <- NA
  return(x)
}
# get lm pval
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
####

# ---- load NRCS CRAs for SW ----
#library(rgdal)
#library(raster)
#AZ_CRA <- readOGR(dsn = "./shapes/AZ_CRA/soils", layer = "cra_a_az")
#NM_CRA <- readOGR(dsn = "./shapes/NM_CRA/soils", layer = "cra_a_nm")
# combine
#SW_CRA<- rbind(AZ_CRA,NM_CRA)

# AZ LRU from Emilio
#AZ_LRU<-readOGR(dsn = "./shapes/AZ_LRU", layer = "mlra_a_az")
#AZ_LRU<-spTransform(AZ_LRU, crs(AZ_CRA))
# ----

# results from analyze_by_LRU_ordReg2.R
#load("~/RProjects/RangeDrought/results/AZ_LRU_smNDVI_detrend_stepwiseAIC_gridmet_SPI3_SPEI3_3cats_results.Rdata")
#load("~/RProjects/RangeDrought/results/AZ_LRU_smNDVI_detrend_stepwise_pval_gridmet_SPI3_SPEI3_3cats_results.Rdata")
load("~/RProjects/RangeDrought/results/AZ_LRU_rpms_detrend_stepwise_pval_gridmet_SPI3_SPEI3_3cats_results.Rdata")

# add month/year to climate dataframes
spiDataFrame$month<-as.numeric(format(spiDataFrame$date, format="%m"))
spiDataFrame$year<-as.numeric(format(spiDataFrame$date, format="%Y"))
speiDataFrame$month<-as.numeric(format(speiDataFrame$date, format="%m"))
speiDataFrame$year<-as.numeric(format(speiDataFrame$date, format="%Y"))

# empty list
droughtCat<-list()
droughtProb<-list()
allDroughtProb<-list()
lowProb<-list()
trends<-list()
coefs<-list()
lowProbHitMiss<-list()

# swith to target climate dataset
climData<-spiDataFrame # CHANGE TO SPI OR SPEI
# switch to forecast type, needs to match dataset
quantCuts<-c(0,0.33,0.66,1); quantNames<-c("Below","Normal","Above") # for 3 cat
#quantCuts<-c(0,0.03,0.06,0.11,0.21,0.30,1); quantNames<- c("D4","D3","D2","D1","D0","No Drought") # for 3 cat

# earliest year for RS data?
yearRS<-1984

# loop through all LMU
for(i in 1:nrow(LMU)){
  # wrangle data
  allSeas<-climData[,c(paste0("LRU-",i),"month","year")] 
  allSeas<-dcast(allSeas, year~month, value.var = paste0("LRU-",i))
    colnames(allSeas)<-c("year",paste0(month.abb[11],"_",month.abb[1]),paste0(month.abb[12],"_",month.abb[2]),
                       paste0(month.abb[seq(1,12-2,1)],"_",month.abb[seq(1+2,12,1)]))
  allSeas<-subset(allSeas, year>=yearRS)
  allSeas$rpms<-rpmsDataFrame[,i]
  # detrend rs data
   tr<-lm(allSeas$rpms~seq(1,nrow(allSeas),1))
      trends[[i]]<-c(coef(tr)[2],lmp(tr))
   allSeas$rpms<-detrend(allSeas$rpms, 'linear')
  # create factors
  #  allSeas$rpms<- quantcut(allSeas$rpms, c(0,0.03,0.06,0.11,0.21,0.30,1)) # USDM percentiles
   allSeas$rpms<- quantcut(allSeas$rpms, quantCuts) 
  # rename levels
  #  levels(allSeas$rpms) <- c("D4","D3","D2","D1","D0","No Drought")
    levels(allSeas$rpms) <- quantNames
  
    # diag plots - not run ----
    #par(mfrow=c(3,1))
    #  plot(allSeas$year, allSeas$rpms)
    #  plot(allSeas$year, allSeas$Feb_Apr, type='b', ylim=c(-2,2))
    #  plot(allSeas$year, allSeas$Mar_May, type='b', ylim=c(-2,2))
    # ----  
      
  # develop prediction model  
  periods<-unlist(str_split(str_remove_all(as.character(LMU@data$form_spi[i]),"[+,]")," "))
  periods<-periods[which(lapply(periods, nchar)!=0)]
    form<-as.formula(paste(c("rpms ~", periods),collapse=" + "))
    
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
        allDroughtProb[[i]]<-matrix(NA, nrow = nrow(climData), ncol = length(quantNames))
        lowProb[[i]]<-NA
        coefs[[i]]<-NA
        lowProbHitMiss[[i]]<-NA
        #return()
    } else {
    # now it is safe to run the fit as we have handled an error above                  
      catModel<-polr(form, data = allSeas) # build model
      coefs[[i]]<-coef(catModel)
    # create custom monthly time series
      newClim<-climData[,c(paste0("LRU-",i),"month","year")] 
    # create time series for each var
      for(j in 1:length(periods)){
        newClim$var1<-newClim[,1]
        newClim$var1[which(newClim$month!=(which(colnames(allSeas)==periods[j])-1))]<-NA
        newClim<-fill(newClim, var1)
        colnames(newClim)[j+3]<-periods[j]
      }
        
      # make predictions  
      predCats<-predict(catModel,newClim,type="probs")
          z<-apply(as.data.frame(predCats),1,which.max)
          z<-max.col(predCats)
          newClim$catNames<-colnames(predCats)[z]
          newClim$probs<-apply(predCats,1,max)
          
      # add in observed rs category for each year
        newClim<-merge(newClim, allSeas[,c("year","rpms")], by="year", all.x=TRUE)
      # find low category hits misses
        newClim$bloHit<-ifelse(newClim$rpms==quantNames[1],
                               ifelse(newClim$rpms==quantNames[1] & newClim$catNames==quantNames[1],"HIT","MISS"),"No Drought")
              
      # store results in list for each LMU
          droughtCat[[i]]<-newClim$catNames
          droughtProb[[i]]<-newClim$probs
          allDroughtProb[[i]]<-predCats
          lowProb[[i]]<-predCats[,1]
          lowProbHitMiss[[i]]<-newClim$bloHit
          
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
# get only below probs
lowProbdf =as.data.frame(do.call(cbind, lowProb))
  colnames(lowProbdf)<-paste0("LRU-",seq(1,ncol(lowProbdf),1))
  lowProbdf$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
# get only below probs
lowProbHitMissdf =as.data.frame(do.call(cbind, lowProbHitMiss))
  colnames(lowProbHitMissdf)<-paste0("LRU-",seq(1,ncol(lowProbHitMissdf),1))
  lowProbHitMissdf$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")  
  
# check drought category frequencies
#table(unlist(droughtCatdf[,1:91]))
  
# verification stats
  qfunc<-function(x)  { qcut<-quantcut(x, quantCuts)
                        levels(qcut) <- seq(1,length(quantCuts)-1,1)
                        qcut<-as.numeric(qcut)
                        return(qcut)
  }
  
rpmsCats<- apply((rpmsDataFrame[,1:ncol(rpmsDataFrame)-1]), 2, FUN = qfunc)

verifList<-list()
contList<-list()
percCorrect<-list()
for(i in 1:ncol(rpmsCats)){
  temp<- as.data.frame(allDroughtProb[[i]])
  if(nrow(temp)!=0 & sum(is.na(temp))!=(nrow(climData)*length(quantNames))){  
    temp$date<-seq(as.Date("1979-01-01", "%Y-%m-%d"), as.Date("2019-12-01", "%Y-%m-%d"), by="months")
      temp<- temp[seq(0, nrow(temp), 12), ]
      temp<-temp[((nrow(temp)-nrow(rpmsCats))+1):nrow(temp),]
      compare<-cbind(rpmsCats[,i],temp[,1:(length(quantCuts)-1)])
      verifList[[i]]<-rps( obs = rpmsCats[,i], pred = base::as.matrix(temp[,1:(length(quantCuts)-1)]))
      # percent correct
      percCorrect[[i]]<- sum(rpmsCats[,i]==apply((temp[,1:(length(quantCuts)-1)]),1,which.max))/length(rpmsCats[,i])*100
      # verif stats from contingency table, needs to be complete cases
      temp<-cbind.data.frame(apply((temp[,1:(length(quantCuts)-1)]),1,which.max),rpmsCats[,i])
        colnames(temp)<-c("forecast","obs")
          diffs<-setdiff(unique(temp$obs),unique(temp$forecast))  
          add<-data.frame("forecast"=diffs, "obs"=rep(NA,length(diffs))) 
          temp<-rbind.data.frame(temp,add)
          contList[[i]]<- multi.cont(as.matrix(table(temp$forecast,temp$obs)))
          
  }else{
    verifList[[i]]<-verifList[[i-1]]
      verifList[[i]][1:3]<-NA
    percCorrect[[i]]<-NA
    contList[[i]]<-contList[[i-1]]
      contList[[i]]<-rapply(contList[[i]],function(x) NA, how = "replace")
  }
}

# get forecast verifications 
hss<-unlist(nullToNA(lapply(X = contList, FUN = `[[`, "hss")))
pss<-unlist(nullToNA(lapply(X = contList, FUN = `[[`, "pss")))
gs<-unlist(nullToNA(lapply(X = contList, FUN = `[[`, "gs")))

# coefficients to dataframe
n.obs <- sapply(coefs, length)
seq.max <- seq_len(max(n.obs))
coefsDF <- as.data.frame(t(sapply(coefs, "[", i = seq.max)))
colnames(coefsDF)<-c("coef1","coef2","coef3")

percCorrect =as.data.frame(do.call(rbind, percCorrect)); mean(percCorrect$V1, na.rm=TRUE)
verifList <-as.data.frame(matrix(unlist(verifList), nrow=length(verifList), ncol=3, byrow=TRUE))
  colnames(verifList)<-c("rps","rpss","rps.clim")
  mean(verifList$rps, na.rm = TRUE); median(verifList$rps, na.rm = TRUE)
  mean(verifList$rpss, na.rm = TRUE); median(verifList$rpss, na.rm = TRUE)
  hist(verifList$rps); hist(hss); hist(pss); hist(gs); hist(verifList$rpss); 
# map rps/rpss
  verifLMU<-LMU
  verifLMU@data<-cbind(verifLMU@data,verifList, percCorrect$V1,hss,pss,gs,coefsDF)
  
# plot LMUs if needed
  plot(LMU)
  plot(LMU[29,], add=TRUE, col='red')
  
  #pal <- RColorBrewer::brewer.pal(n=9, 'OrRd')
  library(rasterVis)
  pal <- RColorBrewer::brewer.pal(n=9, 'OrRd')
  spplot(verifLMU, c("rps"), 
         main=list(label=paste0("Rank Probability Score - NDVI/SPI/3 cat (",round(mean(verifList$rps, na.rm = TRUE),2),")"),
                   cex=1),scales=list(draw = TRUE),
         at=seq(0, 0.5, 0.1),col.regions = pal)
  
  spplot(verifLMU, c("r2_spei"), 
         main=list(label=paste0("R2 - NDVI/SPI/3 cat (",round(mean(LMU@data$r2_spei, na.rm = TRUE),2),")"),
                   cex=1),scales=list(draw = TRUE),
         at=seq(0, 0.9, 0.1),col.regions = pal)
 
  spplot(verifLMU, c("rpss"), 
         main=list(label=paste0("Rank Probability Skill Score - NDVI/SPEI/3 Cat (",round(mean(verifList$rpss, na.rm = TRUE),2),")"),
                   cex=1),scales=list(draw = TRUE),
         at=seq(-1, 1, 0.2),col.regions = rev(brewer.pal(n = 10, 'RdBu')))+
    layer(sp.polygons(AZ_LRU, line="black"))
  
  pal <- RColorBrewer::brewer.pal(n=9, 'OrRd')
  spplot(verifLMU, c("percCorrect$V1"), 
         main=list(label=paste0("Percent Correct - RPMS/SPEI/USDM cat (",round(mean(percCorrect$V1, na.rm=TRUE),2),")"),
                   cex=1),scales=list(draw = TRUE),
         at=seq(0, 90, 10),col.regions = pal)
  
  spplot(verifLMU, c("hss"), 
         main=list(label=paste0("Heidke Skill Score - smNDVI/SPEI/3 Cat (",round(mean(hss, na.rm = TRUE),2),")"),
                   cex=1),scales=list(draw = TRUE),
         at=seq(-1, 1, 0.2),col.regions = rev(brewer.pal(n = 10, 'RdBu')))

  speirpss<-verifLMU@data$rpss
  verifLMU@data$rpssdiff<-verifLMU@data$rpss-speirpss
  spplot(verifLMU, c("rpssdiff"), 
         main=list(label=paste0("SPI-SPEI RPSS - smNDVI/3 Cat (",round(mean( verifLMU@data$rpssdiff, na.rm = TRUE),2),")"),
                   cex=1),scales=list(draw = TRUE),
         at=seq(-0.2, 0.2, 0.05),col.regions = rev(brewer.pal(n = 10, 'RdBu')))
  
  # coefs
  spplot(verifLMU, c("coef1"), 
         main=list(label=paste0("Coef1 - detrended pVal-Step smNDVI/3 Cat"),
                   cex=1),scales=list(draw = TRUE),
         at=seq(-3.5, 3.5, 0.75),col.regions = rev(brewer.pal(n = 10, 'RdBu')))
  
  
# map trends in rs
  slope<-unlist(nullToNA(lapply(X = trends, FUN = `[[`, 1)))
  pval<-unlist(nullToNA(lapply(X = trends, FUN = `[[`, 2)))
  slope[which(pval>0.05)]<-NA  
  verifLMU@data<-cbind(verifLMU@data, slope)
  spplot(verifLMU, c("slope"), 
         main=list(label="Significant slopes in smNDVI data",
                   cex=1),scales=list(draw = TRUE),
         at=seq(-0.005, 0.005, 0.001),col.regions = rev(brewer.pal(n = 10, 'RdBu')))

# different models for SPI/SPEI?
verifLMU@data$diffModels<-as.factor(mapply(identical, as.character(verifLMU@data$form_spi), as.character(verifLMU@data$form_spei)))
spplot(verifLMU, c("diffModels"), 
       main=list(label="Are SPI and SPEI models different?",
                 cex=1),scales=list(draw = TRUE))
diffModels<-as.data.frame(table(verifLMU@data$form_spi,verifLMU@data$form_spei))
         diffModels<-diffModels[which(diffModels$Freq!=0),]
         diffModels<-diffModels[which(mapply(identical, as.character(diffModels$Var1), as.character(diffModels$Var2))==FALSE),]
         
         
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
  
 
# MAP only lowest category probabilities
    temp<-as.data.frame(t(lowProbdf[,1:91]))
    mapLo<-LMU
    #colnames(temp)<-paste0("X",droughtCatdf[,92])
    mapLo@data<-cbind.data.frame(LMU@data,temp)
    dateNames<-format(seq(as.Date("1979.01.01", "%Y.%m.%d"),
                          as.Date("2019.12.01", "%Y.%m.%d"), by="months"),"%m.%d.%Y")
    colnames(mapLo@data)[(ncol(LMU@data)+1):(ncol(temp)+(ncol(LMU@data)))]<-paste0("X",dateNames)
    
    # plot drought cats
    dates<-seq(as.Date("1979.01.01", "%Y.%m.%d"),
               as.Date("2019.12.01", "%Y.%m.%d"), by="months")
    dateIdx<-which(dates>=as.Date("1998-01-01","%Y-%m-%d") & dates<=as.Date("2003-12-01","%Y-%m-%d"))+(ncol(LMU@data))
    library(RColorBrewer)
    library(rasterVis)
    p<-spplot(mapLo, dateIdx, par.settings = list(strip.background=list(col="lightgrey")),
           main=list(label="Probability of being in lowest tercile - 3mo SPI predicts max-NDVI",cex=1),scales=list(draw = FALSE),
           at=seq(0.3, 1, 0.1),col.regions = (brewer.pal(n = 9, 'Oranges')),as.table=TRUE, layout=c(12,6))+
      layer(sp.polygons(mapLo[which(mapLo@data$spi_pVal>=0.05),], fill='grey'))
    
    # plot to png
     png("/home/crimmins/RProjects/RangeDrought/figs/prob_lowest_tercile_SPI3mo_maxNDVI.png", width = 21, height = 15, units = "in", res = 300L)
     #grid.newpage()
     print(p, newpage = FALSE)
     dev.off()
     
# hit or miss of monthly lowest tercile probability
     # MAP only lowest category probabilities
     temp<-as.data.frame(t(lowProbHitMissdf[,1:91]))
     temp[] <- lapply(temp, factor, 
                      levels=c("HIT","MISS","No Drought"), 
                      labels = c("HIT","MISS","No Drought"))
     mapLo<-LMU
     #colnames(temp)<-paste0("X",droughtCatdf[,92])
     mapLo@data<-cbind.data.frame(LMU@data,temp)
     dateNames<-format(seq(as.Date("1979.01.01", "%Y.%m.%d"),
                           as.Date("2019.12.01", "%Y.%m.%d"), by="months"),"%m.%d.%Y")
     colnames(mapLo@data)[(ncol(LMU@data)+1):(ncol(temp)+(ncol(LMU@data)))]<-paste0("X",dateNames)
     
     # plot drought cats
     dates<-seq(as.Date("1979.01.01", "%Y.%m.%d"),
                as.Date("2019.12.01", "%Y.%m.%d"), by="months")
     dateIdx<-which(dates>=as.Date("2015-01-01","%Y-%m-%d") & dates<=as.Date("2019-12-01","%Y-%m-%d"))+(ncol(LMU@data))

     p<-spplot(mapLo, dateIdx, par.settings = list(strip.background=list(col="lightgrey")),
               main=list(label="Hit/Miss Below Outlook - 3mo SPI predicts max-NDVI",cex=1),scales=list(draw = FALSE),
               col.regions = c("green","red","white"),as.table=TRUE, layout=c(12,5))+
       layer(sp.polygons(mapLo[which(mapLo@data$spi_pVal>=0.05),], fill='grey'))
     
     # plot to png
     png("/home/crimmins/RProjects/RangeDrought/figs/hit_miss_lowest_tercile_SPI3mo_maxNDVI_2015_2019.png", width = 21, height = 15, units = "in", res = 300L)
     #grid.newpage()
     print(p, newpage = FALSE)
     dev.off()    
     
   