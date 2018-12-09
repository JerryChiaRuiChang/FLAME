
partial_CATE <- function(FLAME_object, num_covs = NULL, cov_name = NULL, cov_val = NULL) {

  len <- lengths(FLAME_object[[1]])
  index <- which(len == num_covs)

  #If there are no matches for the number of covariates, then return strings
  if (length(index) == 0) {
    return(paste("There are no matched units for", toString(num_covs), "covariate(s)",sep = " "))
  }

  # Get data frame from FLAME object
  CATE_df <- FLAME_object[[2]][[which(len == num_covs)]]

  #If there are no matches for the number of covariates, then return strings
  if (nrow(CATE_df) == 0) {
    return(paste("There are no matched units for", toString(num_covs), "covariate(s)",sep = " "))
  }

  # If the user does not provide cov_name and cov_val, then return data frame
  if (is.null(cov_name) & is.null(cov_val)) {
    return(CATE_df)
  }

  # If the user provides cov_name and cov_val, then return what's being matched with that combination
  if (!is.null(cov_name) & !is.null(cov_val)) {

    # Consider all edge cases

    if (length(cov_name) != num_covs) {
      return("Number of covariates provided does not match with length of covariate names")
    }

    if (length(cov_val) != num_covs) {
      return("Number of covariates provided does not match with length of covariate values")
    }

    if (length(cov_val) != length(cov_name)) {
      return("Length of covariate names does not match with length of covariate values")
    }

    if (!SameElements(colnames(CATE_df)[1:num_covs],cov_name)) {
      return("no CATE for such combination is available using this procedure")
    }

    # Sort cov_val based on the order of colnames in CATE_df
    cov_val_sort <- cov_val[match(colnames(CATE_df)[1:num_covs],cov_name)]
    cov_val_sort <- as.data.frame(t(cov_val_sort), stringsAsFactors=FALSE)
    colnames(cov_val_sort) <- colnames(CATE_df)[1:num_covs]

    L <- sapply(1:nrow(CATE_df), function(x)
      isTRUE(all.equal(cov_val_sort,CATE_df[x,1:num_covs],check.attributes = FALSE)))
    if (nrow(CATE_df[L,]) == 0) {
      return("no CATE for such combination is available using this procedure")
    }
    else {
      rownames(CATE_df) <- NULL
      return(CATE_df[L,])
    }
  }
}

SameElements <- function(a, b) return(identical(sort(a), sort(b)))



all_CATE <- function(FLAME_object, cov_name, cov_val) {

  # If user does not provide cov_name and cov_val, then return String
  if (is.null(cov_name) | is.null(cov_val)) {
    return("Please provide correct input.")
  }

  # If the lengths of cov_name and cov_val are different, then return String
  if (length(cov_val) != length(cov_name)) {
    return("Length of covariate names does not match with length of covariate values")
  }

  # Get all levels which use cov_name for matching
  index <- sapply(FLAME_object[[1]], function(x) return(Reduce("&",cov_name %in% x)))
  all_df <- FLAME_object[[2]][index]
  all_df <- all_df[sapply(all_df,function(x) return(nrow(x) != 0))] #if no matches, then rm it

  #If there are no matches for such procedure, then return String
  if (length(all_df) == 0) {
    return("no CATE for such combination is available using this procedure")
  }

  return_df <- sapply(all_df, find_match, cov_name, cov_val)
  return_df <- return_df[sapply(return_df, function(x) return(!is.null(x)))]
  rownames(return_df) <- NULL

  if (length(return_df) == 0) {
    return("no CATE for such combination is available using this procedure")
  }
  else {
    return(return_df)
  }
}

find_match <- function(df, cov_name, cov_val) {

  CATE_df <- df[,colnames(df) %in% cov_name]
  # Sort cov_val based on the order of colnames in CATE_df
  cov_val_sort <- cov_val[match(colnames(CATE_df),cov_name)]
  cov_val_sort <- as.data.frame(t(cov_val_sort), stringsAsFactors=FALSE)
  colnames(cov_val_sort) <- colnames(CATE_df)

  L <- sapply(1:nrow(CATE_df), function(x)
    isTRUE(all.equal(cov_val_sort,CATE_df[x,],check.attributes = FALSE)))

  if (nrow(df[L,]) != 0) {
    rownames(df) <- NULL
    return(df[L,])
  }
}

#' Get the Size and CATE of Matched Group(s)
#'
#' \code{CATE} provides detailed information of the conditional average treatment effect (CATE)
#' and size of each matched group. First, given number of covariates used for matching,
#' \code{CATE(FLAME_object, num_covs = x)} returns the covariate values, size
#' and CATE of each matched group. Second, given a covariate combination,
#' \code{CATE(FLAME_object, num_covs = x, cov_name = c("x1", "x2", ...), cov_val = c(0,1,...))}
#' returns the CATE and size of the matched group. Third, if user would like to see all
#' matched groups given a specific covariate combination even when the FLAME
#' algorithm performs matching with more than the number of covariates
#' specified, \code{CATE(FLAME_object, cov_name = c("x1", "x2", ...), cov_val =
#' c(0,1,...))} returns all matched groups containing the covariate combination.
#'
#' @param FLAME_object object returned by applying the FLAME algorithm
#'   (\code{\link{FLAME_bit}} or \code{\link{FLAME_PostgreSQL}} or
#'   \code{\link{FLAME_SQLite}})
#' @param num_covs number of covariates used for matching
#' @param cov_name a vector of covariate names
#' @param cov_val  a vector of covariate values, where the value position should
#'   match cov_name position. In addition, it has to be in character R data
#'   type.
#' @examples
#' \donttest{
#' data(toy_data)
#' result <- FLAME::FLAME_bit(data = toy_data, holdout = toy_data)
#' FLAME::CATE(FLAME_object = result, num_covs  = 2)
#' FLAME::CATE(FLAME_object = result, num_covs  = 2, cov_name = c("X1", "X2"), cov_val = c("2", "2"))
#' FLAME::CATE(FLAME_object = result, cov_name = c("X1", "X2"), cov_val = c("2", "2"))
#' }
#' @return data frame with covariate values, CATE, and size of each matched
#'   group
#' @export

CATE <- function(FLAME_object, num_covs = NULL, cov_name = NULL, cov_val = NULL) {

  if (!is.null(cov_name) && typeof(cov_name) != "character") {
    stop("Please specify cov_name as a vector with character data type")
  }

  if (!is.null(cov_val) && typeof(cov_val) != "character") {
    stop("Please specify cov_val as a vector with character data type")
  }

  # If num_covs is provided but not cov_name or cov_val
  # ==> dataframe with all matched units given number of covariates used for matching
  # If num_covs, cov_name, and cov_val are provided
  # ==> dataframe with only units matching user provided cov_name and cov_val

  if (!is.null(num_covs)) {
    return(partial_CATE(FLAME_object, num_covs, cov_name, cov_val))
  }

  # If num_covs is not provided but cov_name and cov_val are provided
  # ==> all matched groups containing cov_name and cov_val

  if(is.null(num_covs)) {
    return(all_CATE(FLAME_object, cov_name, cov_val))
  }
}
