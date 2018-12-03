#' Compute Estimated Treatment Effects
#'
#' \code{AVG_EFFECT} computes estimated treatment effects. Estimated treatment
#' effect is the weighted average of CATEs, with weight being the number of units
#' in each matched group.
#' @examples
#' data(toy_data)
#' result <- FLAME::FLAME_bit(data = toy_data, holdout = toy_data)
#' CATE_object <- FLAME::CATE(FLAME_object = result, cov_name = c("X1", "X2"), cov_val = c("2", "2"))
#' FLAME::AVG_EFFECT(CATE_object)
#' @param CATE_object object returned by applying \code{\link{CATE}} function
#' @return estimated treatment effects
#' @export

AVG_EFFECT <- function(CATE_object) {

  if(is.data.frame(CATE_object)) {
    effect <- CATE_object[,which(colnames(CATE_object) == "effect")]
    size <- CATE_object[,which(colnames(CATE_object) == "size")]
  }

  else {
    effect <- unlist(sapply(CATE_object, function(x) x[,which(colnames(x) == "effect")]))
    size <- unlist(sapply(CATE_object, function(x) x[,which(colnames(x) == "size")]))
  }

  return(sum(effect * size)/sum(size))
}



