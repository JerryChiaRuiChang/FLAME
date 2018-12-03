# Simp function extract "effect" and "size columns from data frame

simp <- function(x) {
  return(x[,which(colnames(x) == "effect"):which(colnames(x) == "size")])}

#' Compute Average Treatment Effect
#'
#' \code{ATE} computes average treatment effect of the matched subsamples by
#' a weighted average of the estimated treatment effects in each matched group.
#' The weight is the number of matched units.
#'
#' @param FLAME_object object returned by applying the FLAME algorithm
#'   (\code{\link{FLAME_bit}}, \code{\link{FLAME_PostgreSQL}}, or
#'   \code{\link{FLAME_SQLite}})
#' @examples
#' data(toy_data)
#' result <- FLAME::FLAME_bit(data = toy_data, holdout = toy_data)
#' FLAME::ATE(result)
#' @return average treatment effect (ATE) of the matched subsamples
#' @export

ATE <- function(FLAME_object) {

  # Get summary data frame with effects and size from all matched units
  CATE_df <- do.call(rbind,lapply(FLAME_object[[2]],simp))

  effect <- CATE_df[,which(colnames(CATE_df) == "effect")]
  size <- CATE_df[,which(colnames(CATE_df) == "size")]

  return(sum(effect * size)/sum(size))
}
