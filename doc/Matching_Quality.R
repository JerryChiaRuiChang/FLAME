## ----setup, include=FALSE, cache=TRUE------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
 set.seed(500)
library(FLAME)

## ---- echo=FALSE---------------------------------------------------------
data <- FLAME::Data_Generation(num_control = 1000, num_treated = 1000,
                               num_cov_dense = 10, num_cov_unimportant = 5, U = 5)
holdout <- data 
head(data)

## ------------------------------------------------------------------------
linear <- FLAME::FLAME_bit(data = data, holdout = holdout, model = "Linear")

## ------------------------------------------------------------------------
ridge <- FLAME::FLAME_bit(data = data, holdout = holdout, model = "Ridge", ridge_reg = 0.1)

## ------------------------------------------------------------------------
lasso <- FLAME::FLAME_bit(data = data, holdout = holdout, model = "Lasso", lasso_reg = 0.1)

## ------------------------------------------------------------------------
tree <- FLAME::FLAME_bit(data = data, holdout = holdout, model = "DecisionTree", tree_depth = 8)

## ------------------------------------------------------------------------
SVM_PE <- function(outcome_treated, outcome_control, covs_treated, covs_control) {
  
  library(e1071) #load e1071 library
  
  # MSE for treated
  model_svm <- svm(outcome_treated ~ covs_treated) # fit the data to SVM model
  pred_treated <- predict(model_svm, covs_treated) #get predicted values
  MSE_treated <- sum((outcome_treated - pred_treated)^2)/length(outcome_treated) # compute mean squared error
  
  # MSE for control
  model_svm <- svm(outcome_control ~ covs_control) # fit the data to SVM model
  pred_control <- predict(model_svm, covs_control) #get predicted values
  MSE_control <- sum((outcome_control - pred_control)^2)/length(outcome_control) # compute mean squared error
  
  return(MSE_treated + MSE_control)
}

## ------------------------------------------------------------------------
SVM <- FLAME::FLAME_bit(data = data, holdout = holdout, PE_function = SVM_PE)

