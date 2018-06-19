#' detailed information of each matched group
#'
#' Given number of covariates used for matching, \code{CATE} returns all matched
#' units and their conditional average treatment effect. If no
#' units are matched, \code{CATE} will return nothing.
#'
#' @param FLAME_object object returned by applying the FLAME algorithm
#'   (\code{\link{FLAME_bit}} or \code{\link{FLAME_PostgreSQL}} or
#'   \code{\link{FLAME_SQLite}})
#' @param num_covs number of covariates used for matching
#' @return data frame with all matched units, the size and CATE of each matched group
#' @export

CATE <- function(FLAME_object,num_covs) {

  len <- lengths(FLAME_object[[1]])
  index <- which(len == num_covs)

  #If there are no matches for the number of covariates, then return strings
  if (length(index) == 0) {
    return(paste("There are no matched units for", toString(num_covs), "covariate(s)",sep = " "))
  }

  # Get data frame from FLAME object
  CATE_df <- FLAME_object[[2]][[which(len == num_covs)]]

  #If there are no matches for the number of covariates, then return strings
  if (nrow(CATE_df) == 0) {
    return(paste("There are no matched units for", toString(num_covs), "covariate(s)",sep = " "))
  }

  return(CATE_df)
}
