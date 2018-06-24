#' compute average CATE
#'
#' \code{CATE_AVG} computes (weighted) average CATE given CATE_object
#'
#' @param CATE_object object returned by applying \code{\link{CATE}} function
#' @return (weighted) average CATE
#' @export

CATE_AVG <- function(FLAME_object,num_covs) {

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



