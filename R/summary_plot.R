
summary_data <- function(num_covs,FLAME_object) {

  # get a list of integer from 1:num_covs in order to know hte index
  len <- lengths(FLAME_object[[1]])
  index <- which(len == num_covs)

  # retrieve more detailed information from FLAME_object[[2]]
  CATE_df <- FLAME_object[[2]][[which(len == num_covs)]]

  # if the data frame is empty, then return NA
  if (nrow(CATE_df) == 0) {
    return(c(NA,NA))
  }

  # get the CATE and its size
  effect <- CATE_df[,which(colnames(CATE_df) == "effect")]
  size <- CATE_df[,which(colnames(CATE_df) == "size")]

  return(c(sum(size),sum(effect * size)/sum(size)))
}

#' visualize matching process
#'
#' \code{summary_plot} visualizes matching process, summarizing
#' each covariate dropped, number of matched units, and average
#' CATE at each iteration.
#'
#' @param FLAME_object object returned by applying the FLAME algorithm
#'   (\code{\link{FLAME_bit}}, \code{\link{FLAME_PostgreSQL}}, or
#'   \code{\link{FLAME_SQLite}})
#' @examples
#' data(toy_data)
#' result <- FLAME::FLAME_bit(data = toy_data, holdout = toy_data)
#' FLAME::summary_plot(result)
#' @return summary plot of the FLAME algorithm
#' @import latticeExtra
#' @importFrom lattice xyplot barchart
#' @export

summary_plot <- function(FLAME_object) {

  # get a list of integer from 1:num_covs
  len <- lengths(FLAME_object[[1]])

  # since nothing is dropped in the beginning, covs_drop = NA
  covs_drop = c("NA")

  # get the covariate dropped at each level
  for(i in 1:(length(FLAME_object[[1]])-1)) {
    cov = setdiff(FLAME_object[[1]][[i]],FLAME_object[[1]][[i+1]])
    if (length(cov) > 1) {
      cov <- paste(cov[1], "+", sep = "")
    }
    covs_drop <- c(covs_drop,cov)
  }

  # get average CATE at each level and combine it with covs_drop into a dataframe
  summary <- data.frame(covs_drop,t(sapply(len,summary_data,FLAME_object)))
  colnames(summary) <- c("covs_dropped","size","CATE")

  # coerece level order for plotting purpose
  summary$covs_dropped <- factor(summary$covs_dropped, levels=unique( as.character(summary$covs_dropped)))

  # construct separate plots for each series
  obj1 <- barchart(size ~ covs_dropped, data = summary, xlab="covariate(s) dropped", ylab = "number of units matched", main = "Summary Plot")
  obj2 <- xyplot(CATE ~ covs_dropped, summary, type = "p", pch = 20, lwd=5, xlab="covariate dropped", ylab = "estimated treatment effects", main = "Summary Plot")

  # make the plot with second y axis
  doubleYScale(obj1, obj2, add.ylab2 = TRUE)
}
