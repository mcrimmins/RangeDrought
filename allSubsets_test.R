# all subsets test
# http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/155-best-subsets-regression-essentials-in-r/

library(tidyverse)
library(caret)
library(leaps)

# Load the data
data("swiss")
# Inspect the data
sample_n(swiss, 3)

models <- regsubsets(Fertility~., data = swiss, nvmax = 5)
summary(models)

res.sum <- summary(models)
data.frame(
  Adj.R2 = which.max(res.sum$adjr2),
  CP = which.min(res.sum$cp),
  BIC = which.min(res.sum$bic)
)

# id: model id
# object: regsubsets object
# data: data used to fit regsubsets
# outcome: outcome variable
get_model_formula <- function(id, object, outcome){
  # get models data
  models <- summary(object)$which[id,-1]
  # Get outcome variable
  #form <- as.formula(object$call[[2]])
  #outcome <- all.vars(form)[1]
  # Get model predictors
  predictors <- names(which(models == TRUE))
  predictors <- paste(predictors, collapse = "+")
  # Build model formula
  as.formula(paste0(outcome, "~", predictors))
}

get_model_formula(3, models, "Fertility")

get_cv_error <- function(model.formula, data){
  set.seed(1)
  train.control <- trainControl(method = "cv", number = 5)
  cv <- train(model.formula, data = data, method = "lm",
              trControl = train.control)
  cv$results$RMSE
}
# Compute cross-validation error
model.ids <- 1:5
cv.errors <-  map(model.ids, get_model_formula, models, "Fertility") %>%
  map(get_cv_error, data = swiss) %>%
  unlist()
cv.errors

# Select the model that minimize the CV error
which.min(cv.errors)

coef(models, 4)

