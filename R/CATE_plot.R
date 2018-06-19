#' summarize CATE of all matched groups by boxplot
#'
#' Given number of covariates used for matching, \code{CATE_plot} visualizes
#' CATE of all matched groups by boxplot
#'
#' @param FLAME_object object returned by applying the FLAME algorithm
#'   (\code{\link{FLAME_bit}} or \code{\link{FLAME_PostgreSQL}} or
#'   \code{\link{FLAME_SQLite}})
#' @param num_covs number of covariates used for matching
#' @return boxplot
#' @export

CATE_plot <- function(FLAME_object,num_covs) {

  # get a list of integer from 1:num_covs
  len <- lengths(FLAME_object[[1]])
  index <- which(len == num_covs)

  #If there are no matches for the number of covariates, then return strings
  if (length(index) == 0) {
    return(paste("There are no matches for", toString(num_covs), "covariate(s)",sep = " "))
  }

  #If there are no matches for the number of covariates, then return strings
  CATE_df <- FLAME_object[[2]][[which(len == num_covs)]]
  if (nrow(CATE_df) == 0) {
    return(paste("There are no matches for", toString(num_covs), "covariate(s)",sep = " "))
  }

  effect <- CATE_df[,which(colnames(CATE_df) == "effect")]

  boxplot(effect, main = paste("CATE summary of", toString(num_covs), "covariate(s)",sep = " "),
          ylab = "CATE", xlab = paste(toString(num_covs), "covariate(s)", sep = " ") )
}

