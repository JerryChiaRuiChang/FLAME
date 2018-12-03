#' Summary
#'
#' \code{FLAME_summary} provides brief summary of FLAME implementation.
#'
#' @param FLAME_object object returned by applying the FLAME algorithm
#'   (\code{\link{FLAME_bit}}, \code{\link{FLAME_PostgreSQL}}, or
#'   \code{\link{FLAME_SQLite}})
#' @return (1) Number of units matched (2) Average treatment effect (ATE)
#' @examples
#' data(toy_data)
#' result <- FLAME::FLAME_bit(data = toy_data, holdout = toy_data)
#' FLAME::summary(result)
#' @export

summary <- function(FLAME_object) {

  #number of matched units
  print(paste("Number of units matched = ", sum(FLAME_object[[4]]['matched'] >= 1)))

  # ATE
  df <- do.call(rbind,lapply(FLAME_object[[2]],simp))

  effect <- df[,which(colnames(df) == "effect")]
  size <- df[,which(colnames(df) == "size")]

  print(paste("Average treatment effect = ", sum(effect * size)/sum(size)))
}

