#' Generate Data
#'
#' @param num_control Number of controlled samples
#' @param num_treated Number of treated samples
#' @param num_cov_dense Number of covariates
#' @param num_covs_unimportant Number of unimportant covariates
#' @return Data Frame
#' @import reticulate
#' @export
#' @examples
#' Data_Generation(10,10,10,0)
#' Data_Generation(20,20,10,10)

#Data Generation Function
Data_Generation <- function(num_control, num_treated, num_cov_dense, num_covs_unimportant) {
  source_python(system.file("python","Data_Generation.py",package = "FLAME"))
  #Convert Data Type from Numeric to Integer
  num_control <- as.integer(num_control)
  num_treated <- as.integer(num_treated)
  num_cov_dense <- as.integer(num_cov_dense)
  num_covs_unimportant <- as.integer(num_covs_unimportant)
  #Pass Arguement to data_generation_dense_2
  dat <- data_generation_dense_2(num_control, num_treated, num_cov_dense, num_covs_unimportant)
  return(dat)
}
