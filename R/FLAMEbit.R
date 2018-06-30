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

  data <- data.frame(data) #Convert input data to data.frame if not already converted
  holdout <- data.frame(holdout) #Convert holdout data to data.frame if not already converted

  data$matched <- 0 #add column matched to input data
  column <- colnames(data)

  # Convert each covariate and treated into type integer
  data[,c(1:num_covs, num_covs+2, num_covs+3)] <- sapply(data[,c(1:num_covs, num_covs+2, num_covs+3)],as.integer)
  holdout[,c(1:num_covs,num_covs+2)] <- sapply(holdout[,c(1:num_covs,num_covs+2)],as.integer)

  # Convert outcome variable to numeric
  data[,num_covs + 1] <- as.numeric(data[,num_covs + 1])
  holdout[,num_covs + 1] <- as.numeric(holdout[,num_covs + 1])

  covs <- as.integer(seq(0,num_covs-1)) #Create list of covariate
  covs_max_list <- as.integer(covs_max_list) #Convert covs_max_list to integer array
  num_treated <- as.integer(num_treated) #Convert numeric to integer
  num_control <- as.integer(num_control) #Convert numeric to integer

  #Call Python File
  source_python(system.file("run_bit.py",package = "FLAME"))

  #source_python("run_bit.py")

  #Run_bit
  result <- run_bit(r_to_py(data), r_to_py(holdout), covs, covs_max_list,
                    num_treated, num_control, tradeoff)


  #Convert column name to match with covariate returned in each iteration
  result_cov <- NULL
  CATE_df <- result[[2]]
  for (i in 1:length(CATE_df)) {
    CATE_df[[i]] <- data.frame(CATE_df[[i]])
    colnames(CATE_df[[i]]) <- c(column[(result[[1]][[i]] + 1)],"effect","size")
    result_cov[[i]] <- column[(result[[1]][[i]] + 1)]
  }

  # Convert return dataframe column name
  result_df <- as.data.frame(result[[3]])
  colnames(result_df) <- column

  # Convert from double ==> integer
  result_df[,c(1:num_covs, num_covs+2, num_covs+3)] <- sapply(result_df[,c(1:num_covs, num_covs+2, num_covs+3)],as.integer)

  return_list <- NULL
  return_list[[1]] <- result_cov
  return_list[[2]] <- CATE_df
  return_list[[3]] <- result_df

  return(return_list)
}



