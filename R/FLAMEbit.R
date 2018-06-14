#' Bit Vectors Implementation
#'
#' \code{FLAME_bit} applies FLAME matching algorithm based on bit vectors implementation.
#'
#' @param data Input data
#' @param holdout Holdout training data
#' @param num_covs Number of covariates
#' @param num_treated Number of units in treated group
#' @param num_control Number of units in control group
#' @param covs_max_list List indicates each covariate is binary/ternary/...
#' @param tradeoff Tradeoff parameter to compute Match Quality
#' @return (1) List of covariates matched at each iteration (2) list of data frame
#' showing the size of matched group and its conditional average treatment effect (CATE)
#' @import reticulate
#' @export

FLAME_bit <- function(data, holdout, num_covs, num_treated, num_control, covs_max_list, tradeoff) {

  #Type Conversion and Preprocessing

  data <- data.frame(data) #Convert input data to data.frame if not already converted
  holdout <- data.frame(holdout) #Convert holdout data to data.frame if not already converted
  column <- colnames(data)

  data$matched <- 0 #add column matched to input data
  holdout$matched <- 0 #add column matched to holdout data

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



