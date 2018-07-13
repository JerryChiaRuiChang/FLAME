#' Summarize CATEs of All Matched Groups by Boxplot
#'
#' Given CATE_object, \code{CATE_plot} visualizes CATEs of all matched groups by
#' boxplot
#'
#' @param CATE_object object returned by applying \code{\link{CATE}} function
#' @return boxplot
#' @export

CATE_plot <- function(CATE_object) {

  if(is.data.frame(CATE_object)) {
    effect <- CATE_object[,which(colnames(CATE_object) == "effect")]
    size <- CATE_object[,which(colnames(CATE_object) == "size")]
  }

  else {
    effect <- unlist(sapply(CATE_object, function(x) x[,which(colnames(x) == "effect")]))
    size <- unlist(sapply(CATE_object, function(x) x[,which(colnames(x) == "size")]))
  }

  boxplot(effect, main = paste("CATE summary with", toString(sum(size)), "matched units",sep = " "),
          ylab = "CATE" )
}
