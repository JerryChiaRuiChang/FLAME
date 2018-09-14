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
holdout <- data 

head(data)

## ---- warning=FALSE------------------------------------------------------
result_bit <- FLAME::FLAME_bit(data = data, holdout = holdout, compute_var = TRUE)

## ------------------------------------------------------------------------
# Connect to PostgreSQL
drv <- dbDriver('PostgreSQL')

# Name the connection as db
db <- dbConnect(drv, user="postgres", dbname="FLAME", host='localhost',
             port=5432, password = 'new_password')

# Run FLAME_PostgreSQL 
result_PostgreSQL <- FLAME::FLAME_PostgreSQL(db = db, data = data, holdout = holdout, compute_var = TRUE)

# Disconnect from db
dbDisconnect(db)

## ------------------------------------------------------------------------
#Name the connection as conn
db <- dbConnect(SQLite(),"tempdb") 

#Run FLAME_SQLite
result_SQLite <- FLAME::FLAME_SQLite(db = db, data = data, holdout = holdout, compute_var = TRUE)

#Disconnect from db
dbDisconnect(db)

## ------------------------------------------------------------------------
FLAME_summary(result_bit)
FLAME_summary(result_SQLite)
FLAME_summary(result_PostgreSQL)

## ------------------------------------------------------------------------
result_bit[[1]][[10]] # bit vectors
result_PostgreSQL[[1]][[10]] #PostgreSQL
result_SQLite[[1]][[10]] #SQLite

## ------------------------------------------------------------------------
head(result_bit[[2]][[10]]) # bit vectors
head(result_PostgreSQL[[2]][[10]]) #PostgreSQL
head(result_SQLite[[2]][[10]]) #SQLite

## ------------------------------------------------------------------------
result_bit[[3]] # bit vectors
result_PostgreSQL[[3]] #PostgreSQL
result_SQLite[[3]] #SQLite

