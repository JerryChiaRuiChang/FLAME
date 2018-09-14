## ----setup, include = FALSE, cache=TRUE----------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
 set.seed(500)
library(FLAME)

## ------------------------------------------------------------------------
data <- FLAME::Data_Generation(num_control = 1000, num_treated = 1000,
                               num_cov_dense = 10, num_cov_unimportant = 5, U = 5)
holdout <- data #Assume holdout training data is the same as input data
result_bit <- FLAME::FLAME_bit(data = data, holdout = holdout, compute_var = TRUE)

## ---- echo=FALSE---------------------------------------------------------
result_bit[[1]][[8]]

## ------------------------------------------------------------------------
head(FLAME::CATE(FLAME_object = result_bit, num_covs  = 8))

## ---- echo=FALSE---------------------------------------------------------
cov_df <- as.matrix(result_bit[[2]][[8]][1,])
cov_name = colnames(cov_df)[1:8]
cov_val = as.vector(cov_df[1,1:8])
cov_df[,1:8]

## ------------------------------------------------------------------------
#covariate names
cov_name 

#covariate values in character R data type
cov_val 

FLAME::CATE(FLAME_object = result_bit, num_covs  = 8, cov_name = cov_name, cov_val = cov_val)

## ------------------------------------------------------------------------
FLAME::MATCH(FLAME_object = result_bit, cov_name = cov_name, cov_val = cov_val)

## ------------------------------------------------------------------------
CATE_object <- CATE(FLAME_object = result_bit, 
                    cov_name = c("x1", "x3", "x5", "x7", "x9"), 
                    cov_val = c("1", "0", "1", "0", "1"))
CATE_object

## ------------------------------------------------------------------------
# Estimated treatment effects 
AVG_EFFECT(CATE_object) 

CATE_plot(CATE_object)

## ------------------------------------------------------------------------
FLAME::ATE(result_bit)

## ------------------------------------------------------------------------
FLAME::summary_plot(result_bit)

