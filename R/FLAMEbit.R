#' bit vectors implementation
#'
#' \code{FLAME_bit} applies FLAME matching algorithm based on bit vectors
#' implementation.
#'
#' @param data input data
#' @param holdout holdout training data
#' @param num_covs number of covariates
#' @param num_treated number of units in treated group
#' @param num_control number of units in control group
#' @param covs_max_list list indicates each covariate is binary/ternary/...
#' @param tradeoff tradeoff parameter to compute Matching Quality
#' @return (1) list of covariates used for matching at each iteration (2) list
#'   of dataframe showing all matched units, size of each matched group, and its
#'   conditional average treatment effect (CATE).
#' @import reticulate
#' @export

FLAME_bit <- function(data, holdout, num_covs, num_treated, num_control, covs_max_list, tradeoff) {

  #Type Conversion and Preprocessing

  data <- data.frame(data) # Convert input data to data.frame if not already converted
  holdout <- data.frame(holdout) #Convert holdout data to data.frame if not already converted
  column <- colnames(data)

  data$matched <- 0 #add column matched to input data
  holdout$matched <- 0 #add column matched to holdout data

  # Convert each covariate and treated into type integer
  data[,c(1:num_covs,num_covs+2)] <- sapply(data[,c(1:num_covs,num_covs+2)],as.integer)
  holdout[,c(1:num_covs,num_covs+2)] <- sapply(holdout[,c(1:num_covs,num_covs+2)],as.integer)

  covs <- as.integer(seq(0,num_covs-1)) #Create list of covariate
  covs_max_list <- as.integer(covs_max_list) #Convert covs_max_list to integer array
  num_treated <- as.integer(num_treated) #Convert numeric to integer
  num_control <- as.integer(num_control) #Convert numeric to integer

  #Call Python File
  source_python(system.file("run_bit.py",package = "FLAME"))

  #Run_bit
  result <- run_bit(r_to_py(data), r_to_py(holdout), covs, covs_max_list,
                    num_treated, num_control, tradeoff)

  #Convert column name to match with covariate returned in each iteration
  result_cov <- NULL
  result_df <- result[[2]]
  for (i in 1:length(result_df)) {
    result_df[[i]] <- data.frame(result_df[[i]])
    colnames(result_df[[i]]) <- c(column[(result[[1]][[i]] + 1)],"effect","size")
    result_cov[[i]] <- column[(result[[1]][[i]] + 1)]
  }

  return_list <- NULL
  return_list[[1]] <- result_cov
  return_list[[2]] <- result_df

  return(return_list)
}



