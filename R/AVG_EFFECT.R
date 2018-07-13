#' Compute Estimated Treatment Effects
#'
#' \code{AVG_EFFECT} computes estimated treatment effects. Estimated treatment
#' effect is the weighted average of CATEs, with weight being the number of units
#' in each matched group.
#'
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



