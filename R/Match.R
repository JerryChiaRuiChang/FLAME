#' get all matched units given certain covariate combination
#'
#' @param FLAME_object object returned by applying the FLAME algorithm
#'   (\code{\link{FLAME_bit}} or \code{\link{FLAME_PostgreSQL}} or
#'   \code{\link{FLAME_SQLite}})
#' @param cov_name a vector of covariate names
#' @param cov_val  a vector of covariate values, where the value position should match cov_name position
#' @return data frame with all matched units
#' @export

MATCH <- function(FLAME_object, cov_name, cov_val) {

  len <- lengths(FLAME_object[[1]])
  index <- which(len == length(cov_name))

  # If there are no matches for the number of covariates, then return strings
  if (length(index) == 0) {
    return("There are no matched units for such combination")
  }

  cov_matched <- FLAME_object[[1]][[index]] # which covariates were used for matching

  # If covariate names provided are not used for matching, then return strings
  if(!identical(sort(cov_matched), sort(cov_name))) {
    return("There are no matched units for such combination")
  }

  if (length(cov_val) != length(cov_name)) {
    return("Length of covariate names does not match with length of covariate values")
  }

  # Sort the covariate values by column order
  cov_val_sort <- as.integer(cov_val[match(cov_matched,cov_name)])
  cov_val_sort <- data.frame(t(unlist(cov_val_sort)))
  colnames(cov_val_sort) <- cov_matched

  df <- FLAME_object[[3]]
  df <- df[df$matched == length(cov_name),]

  L <- sapply(1:nrow(df), function(x) identical(unlist(df[x,cov_matched]),unlist(cov_val_sort)))

  if (nrow(df[L,]) == 0) {
    return("no CATE for such combination is available using this procedure")
  }
  else {
    return_df <- df[L,-ncol(df)]
    row.names(return_df) <- NULL
    return(return_df)
  }
}


#cov_name <- c("x1", "x2", "x3", "x4", "x5", "x9")
#cov_val <- c(0,0,0,1,1,0)
#Match(result_bit, cov_name, cov_val)
#Match(result_PostgreSQL, cov_name, cov_val)
#Match(result_SQLite, cov_name, cov_val)
