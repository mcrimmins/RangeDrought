# scratch code for all subsets polr
# https://ryouready.wordpress.com/2009/02/06/r-calculating-all-possible-linear-regression-models-for-a-given-set-of-predictors/

regressors <- c("y1", "y2", "y3", "y4")

vec <- c(T, F, T, F)
paste(regressors[vec])

paste(regressors[vec], collapse=" + ")

paste(c("y ~ 1", regressors[vec]), collapse=" + ")

as.formula(paste(c("y ~ 1", regressors[vec]), collapse=" + "))

regMat <- expand.grid(c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE))



regMat <- regMat[-(dim(regMat)[1]),]

# let's name the columns
names(regMat) <- paste("y", 1:4, sep="")

# x1 will be dependent
regressors <- c("y1", "y2", "y3", "y4")

allModelsList <- apply(regMat, 1, function(x) as.formula(
  paste(c("x1 ~ 1", regressors[x]),
        collapse=" + ")) )

# for Range Drought
regMat <- expand.grid(c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE))

# find one and two factor models
regMat <- regMat[which(rowSums(regMat == "TRUE")<=2),]
regMat <- regMat[-(dim(regMat)[1]),]
names(regMat) <- colnames(trimSeas[1:ncol(trimSeas)-1])

# rename
regressors <- colnames(trimSeas[1:ncol(trimSeas)-1])

allModelsList <- apply(regMat, 1, function(x) as.formula(
  paste(c("rpms ~", regressors[x]),
        collapse=" + ")) )


# stepwise model selection
m <- polr(rpms~., data = trimSeas, Hess=TRUE)
stepM<-stepAIC(m, trace=TRUE, direction="both")
stepM$anova

new_data <- data.frame("Dec_Feb"= -1,"Feb_Apr"=-1,"Mar_May"=0,"Jul_Sep"=0,"Sep_Nov"=0)
round(predict(stepM,new_data,type = "p"), 3)

#cor(trimSeas[,2:13])

# rms package
library(rms)
(lrmFit <- lrm(rpms~., data=trimSeas))
test<-fastbw(lrmFit, rule = "aic")

lrm(rpms~Jul_Sep, data = trimSeas)

# bestglm package
library(bestglm)



