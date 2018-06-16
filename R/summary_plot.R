summary_data <- function(num_covs,FLAME_object) {
  len <- lengths(FLAME_object[[1]])
  index <- which(len == num_covs)

  CATE_df <- FLAME_object[[2]][[which(len == num_covs)]]
  if (nrow(CATE_df) == 0) {
    return(c(NA,NA))
  }

  effect <- CATE_df[,which(colnames(CATE_df) == "effect")]
  size <- CATE_df[,which(colnames(CATE_df) == "size")]

  return(c(sum(size),sum(effect * size)/sum(size)))
}

#' Summary Plot
#'
#' \code{summary_plot} visualizes the covariate dropped, number of matched units, and its
#' conditional average treatment effect at each iteration.
#'
#' @param FLAME_object object returned by applying FLAME matching algorithm
#'   (\code{\link{FLAME_bit}} or \code{\link{FLAME_PostgreSQL}} or
#'   \code{\link{FLAME_SQLite}})
#' @return summary plot of FLAME algorithm
#' @import latticeExtra
#' @export

summary_plot <- function(FLAME_object) {
  len <- lengths(FLAME_object[[1]])
  covs_drop = c("NA")
  for(i in 1:(length(FLAME_object[[1]])-1)) {
    covs_drop <- c(covs_drop,setdiff(FLAME_object[[1]][[i]],FLAME_object[[1]][[i+1]]))
  }
  summary <- data.frame(covs_drop,t(sapply(len,summary_data,FLAME_object)))
  colnames(summary) <- c("covs_dropped","size","CATE")
  ate <- FLAME::ATE(FLAME_object)


  # --> construct separate plots for each series
  obj1 <- barchart(size ~ covs_dropped, data = summary, xlab="covariate dropped", ylab = "number of matched units", main = "Summary Plot")
  obj2 <- xyplot(CATE ~ covs_dropped, summary, type = "p", pch = 20, lwd=5, xlab="covariate dropped", main = "Summary Plot")

  # --> Make the plot with second y axis:
  doubleYScale(obj1, obj2, add.ylab2 = TRUE)
}
