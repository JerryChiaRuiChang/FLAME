#' Bit Vector Algorithm
#'
#' @param df Data Frame
#' @param holdout Holdout Training Data
#' @param covs List of covariates
#' @param covs_max_list List indicates each covariate is binary/ternary/...
#' @param num_treated Number of units in treated group
#' @param num_control Number of units in control group
#' @return Data Frame
#' @import reticulate
#' @export

FLAME_bit <- function(df, holdout, covs, covs_max_list,num_treated,num_control) {
  df <- data.frame(df)
  covs <- as.integer(covs)
  holdout <- data.frame(holdout)
  covs_max_list <- as.integer(covs_max_list)
  source_python(system.file("run_bit.py",package = "FLAME"))
  return(run_bit(r_to_py(df),r_to_py(holdout),covs,covs_max_list,as.integer(num_treated),
                 as.integer(num_control)))
}
