#' Compute Conditional Average Treatment Effect
#'
#' \code{CATE} computes conditional average treatment effect given number of
#' covariates.
#'
#' @param FLAME_object object returned by applying FLAME matching algorithm
#'   (\code{\link{FLAME_bit}} or \code{\link{FLAME_PostgreSQL}} or
#'   \code{\link{FLAME_SQLite}})
#' @param num_covs number of covariates
#' @return Conditional Average Treatment Effect (CATE)
#' @export

CATE <- function(FLAME_object,num_covs) {

  len <- lengths(FLAME_object[[1]])
  index <- which(len == num_covs)

  #If there are no matches for the number of covariates, then return strings
  if (length(index) == 0) {
    return(paste("There are no matches for", toString(num_covs), "covariate",sep = " "))
  }

  #If there are no matches for the number of covariates, then return strings
  CATE_df <- FLAME_object[[2]][[which(len == num_covs)]]
  if (nrow(CATE_df) == 0) {
    return(paste("There are no matches for", toString(num_covs), "covariate",sep = " "))
  }
  effect <- CATE_df[,which(colnames(CATE_df) == "effect")]
  size <- CATE_df[,which(colnames(CATE_df) == "size")]
  return(sum(effect * size)/sum(size))
}

simp <- function(x) {
  return(x[,which(colnames(x) == "effect"):which(colnames(x) == "size")])}

#'Compute Average Treatment Effect
#'
#' \code{ATE} computes average treatment effect for all matched groups.
#'
#' @param FLAME_object object returned by applying FLAME matching algorithm
#' (\code{\link{FLAME_bit}} or \code{\link{FLAME_PostgreSQL}} or
#' \code{\link{FLAME_SQLite}})
#' @return Average Treatment Effect (ATE)
#' @export

ATE <- function(FLAME_object) {
  CATE_df <- do.call(rbind,lapply(FLAME_object[[2]],simp))
  effect <- CATE_df[,which(colnames(CATE_df) == "effect")]
  size <- CATE_df[,which(colnames(CATE_df) == "size")]
  return(sum(effect * size)/sum(size))
}
