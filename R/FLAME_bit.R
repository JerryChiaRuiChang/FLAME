# Aggregate Function for getting number of times each value occurs
aggregate_table <- function(tab, list_val) {
  tab = unclass(tab)
  name = as.integer(names(tab))
  return(as.vector(tab[sapply(list_val, function(x) which(name ==  x))]))
}


# update_matched_bit takes a dataframe, a set of covariates to match on,
# the treatment indicator column and the matched indicator column.
# it returns the array indicating whether each unit is matched (the first return value),
# and a list of indices for the matched units (the second return value)

update_matched_bit <- function(data, cur_covs, covs_max_list) {

  data_wo_t = as.matrix(data[,cur_covs+1]) # the covariates values as a matrix

  # Compute b_u
  multiplier <- mapply(function(x,y) x^y, covs_max_list, seq(0,length(cur_covs)-1))
  b_u = sapply(1:nrow(data_wo_t), function(x) data_wo_t[x,] %*% multiplier)

  # Compute b_u+
  multiplier <- mapply(function(x,y) x^y, covs_max_list, seq(1,length(cur_covs)))
  b_u_plus = sapply(1:nrow(data_wo_t), function(x) data_wo_t[x,] %*% multiplier)
  b_u_plus = b_u_plus + data[,'treated']

  # Compute c_u
  c_u = aggregate_table(table(b_u), b_u)

  # Compute c_u+
  c_u_plus = aggregate_table(table(b_u_plus), b_u_plus)

  match_index = mapply(function(x,y) x != y, c_u, c_u_plus)
  index = b_u[match_index]
  return(list(match_index, index))
}

# Convert b_u to its original form

num_2_vector <- function(num, covs_max_list) {
  res = list()
  num = as.integer(num)
  for (i in rev(seq(0,length(covs_max_list) - 1))) {
    res = c(res, (num %/% as.integer(covs_max_list[i+1]^(i))))
    num = num %% as.integer(covs_max_list[i+1]^(i))
  }
  return(rev(unlist(res)))
}

#get_CATE function takes match_index and index (b_u values)
# and return dataframe that includes
#(1) list of covariates that are used to match at level l
#(1) conditional average treatment effect (effect)
#(2) size of each matched group (size)

get_CATE_bit <- function(data, match_index, index, cur_covs, covs_max_list, column, factor_level) {
  if (length(index) == 0) {
    CATE <- setNames(data.frame(matrix(ncol = length(cur_covs)+2, nrow = 0)),
                     c(column[(cur_covs + 1)],"effect","size"))
  }

  else {
    d = data[match_index,]
    d[,'b_u'] = index
    summary = data.frame(d %>% group_by(b_u,treated) %>% summarise(size = length(outcome), mean = mean(outcome)))
    summary = data.frame(summary %>% group_by(b_u) %>% summarize(size = sum(size), treated_lst = list(treated), mean_lst = list(mean)))
    CATE = data.frame(t(sapply(summary$b_u, num_2_vector, covs_max_list)))
    CATE$effect = mapply(function(x,y) x[which(y == 1)] - x[which(y == 0)], summary$mean_lst, summary$treated_lst)
    CATE$size = summary$size
    CATE <- CATE[order(CATE$effect),]
    colnames(CATE) = c(column[(cur_covs + 1)],"effect","size")
    CATE[,1:length(cur_covs)] <- mapply(function(x,y) factor_level[x,][CATE[,y]+1], cur_covs + 1, 1:length(cur_covs))
    rownames(CATE) = NULL
  }

  return(CATE)
}


#match_quality function takes holdout dataset, number of total covariates,
#list of current covariates, covariate c to temporily remove from, and trafeoff
#parameter as input. The function then computes Balancing Factor and Predictive Error,
#returning Match Quality.

match_quality_bit <- function(c, data, holdout, num_covs, cur_covs, covs_max_list, tradeoff,
                              PE_function, model, ridge_reg, lasso_reg, tree_depth) {

  # temporarly remove covariate c
  covs_to_match = cur_covs[cur_covs != c]
  covs_max_to_match = covs_max_list[-which(cur_covs == c)]

  # Calculate number of units unmatched (available)

  num_control = nrow(data[data[,'treated'] == 0,])
  num_treated = nrow(data[data[,'treated'] == 1,])

  # Number of matched units

  match_index = update_matched_bit(data, covs_to_match, covs_max_to_match)[[1]]

  num_control_matched = nrow(data[match_index & data[,'treated'] == 0,])
  num_treated_matched = nrow(data[match_index & data[,'treated'] == 1,])

  #Compute Predictive Error

  if (!is.null(PE_function)) {
    # Compute -PE based on user defined PE_function
    outcome_treated <- holdout[holdout[,'treated'] == 1,][,'outcome']
    outcome_control <- holdout[holdout[,'treated'] == 0,][,'outcome']
    covs_treated <- as.matrix(holdout[holdout[,'treated'] == 1,][,covs_to_match + 1])
    covs_control <- as.matrix(holdout[holdout[,'treated'] == 0,][,covs_to_match + 1])
    PE <- -PE_function(outcome_treated, outcome_control, covs_treated, covs_control)
  }

  else {
    if (!is.null(model)) {

      # Linear Regression
      if (model == "Linear") {
        source_python(system.file("Linear.py",package = "FLAME"))
        parameter = 0
      }

      # Ridge Regression
      if (model == "Ridge" && !is.null(ridge_reg)) {
        source_python(system.file("Ridge.py",package = "FLAME"))
        parameter = ridge_reg
      }

      # Lasso
      if (model == "Lasso" && !is.null(lasso_reg)) {
        source_python(system.file("Lasso.py",package = "FLAME"))
        parameter = lasso_reg
      }

      # Decision Tree
      if (model == "DecisionTree" && !is.null(tree_depth)) {
        source_python(system.file("DecisionTree.py",package = "FLAME"))
        parameter = tree_depth
      }
    } else {
      # Default Model is Ridge Regression with L2 Regularization = 0.1
      source_python(system.file("Ridge.py",package = "FLAME"))
      parameter = 0.1
    }

    #Compute PE based on Python scikit learn function
    if (length(covs_to_match) == 1) {
      PE <- predictive_error(r_to_py(holdout), as.integer(num_covs), list(covs_to_match), parameter)
    } else {
      PE <- predictive_error(r_to_py(holdout),as.integer(num_covs),covs_to_match, parameter)
    }
  }

  if (num_control == 0 | num_treated == 0) {
    return(PE)
  } else {
    BF <- num_control_matched/num_control + num_treated_matched/num_treated #Compute Balancing Factor
    return(tradeoff * BF + PE)
  }
}


#' Bit Vectors Implementation
#'
#' \code{FLAME_bit} applies FLAME matching algorithm based on bit vectors
#' implementation.
#'
#' @param data input data
#' @param holdout holdout training data
#' @param num_covs number of covariates
#' @param covs_max_list list indicates each covariate is binary/ternary/...
#' @param tradeoff tradeoff parameter to compute Matching Quality (default =
#'   0.1)
#' @param PE_function user defined function to compute predivtive error
#'   (optional)
#' @param model user defined model - Linear, Ridge, Lasso, or DecisionTree
#'   (optional)
#' @param ridge_reg L2 regularization parameter if model = Ridge (optional)
#' @param lasso_reg L1 regularization parameter if model = Lasso (optional)
#' @param tree_depth maximum depth of decision tree if model = DecisionTree
#'   (optional)
#' @return (1) list of covariates FLAME performs matching at each iteration, (2)
#' list of dataframe showing matched groups' sizes and conditional average
#' treatment effects (CATEs) at each iteration, (3) matching quality at each
#' iteration, and (4) the original data with additional column *matched*,
#' indicating the number of covariates each unit is matched. If a unit is never
#' matched, then *matched* will be 0.
#' @import dplyr
#' @import reticulate
#' @export

FLAME_bit <- function(data, holdout, num_covs, covs_max_list, tradeoff = 0.1, PE_function = NULL,
                             model = NULL, ridge_reg = NULL, lasso_reg = NULL, tree_depth = NULL) {

  if (Reduce("|", sapply(1:num_covs, function(x) !is.factor(data[,x] ))) |
      Reduce("|", sapply(1:num_covs, function(x) !is.factor(holdout[,x] )))) {
    stop("Covariates are not factor data type.")
  }

  if (!is.factor(data[,num_covs + 2]) | !is.factor(holdout[,num_covs + 2])) {
    stop("Treatment variable is not factor data type")
  }

  if (!is.numeric(data[,num_covs + 1]) | !is.numeric(holdout[,num_covs + 1])) {
    stop("Outcome variable is not numeric data type")
  }


  data$matched <- as.integer(0) #add column matched to input data
  column <- colnames(data)

  factor_level <- t(sapply(data[,1:num_covs], levels)) # Get levels of each factor

  # Convert each covariate and treated into type integer
  data[,c(1:num_covs)] <- sapply(data[,c(1:num_covs)], function(x) as.integer(x) - 1)
  data[,num_covs + 2] <- as.integer(levels(data[,num_covs+2])[data[,num_covs+2]])

  holdout[,c(1:num_covs)] <- sapply(holdout[,c(1:num_covs)],function(x) as.integer(x) - 1)
  holdout[,num_covs + 2] <- as.integer(levels(holdout[,num_covs+2])[holdout[,num_covs+2]])

  # Convert outcome variable to numeric
  data[,num_covs + 1] <- as.numeric(data[,num_covs + 1])
  holdout[,num_covs + 1] <- as.numeric(holdout[,num_covs + 1])

  #change input data and holdout training data column name
  colnames(data) <- c(paste("x",seq(0,num_covs-1), sep = ""),"outcome","treated","matched")
  colnames(holdout) <- c(paste("x",seq(0,num_covs-1), sep = ""),"outcome","treated")

  #Set up return objects

  covs_list = list() #list of covariates for matching at each level
  CATE = list() #list of dataframe that calculates conditional average treatment effect at each level
  SCORE = list()

  #Initialize the current covariates to be all covariates and set level to 1

  cur_covs = seq(0,num_covs - 1)
  level = 1

  # Get matched units without dropping anything

  return_match = update_matched_bit(data, cur_covs, covs_max_list)
  match_index = return_match[[1]]
  index = return_match[[2]]

  # Set matched = num_covs and get those matched units
  data[match_index,'matched'] = length(cur_covs)
  return_df = data[match_index,]

  covs_list[[level]] <- column[(cur_covs + 1)]
  CATE[[level]] <- get_CATE_bit(data, match_index, index, cur_covs, covs_max_list, column, factor_level)

  # Remove matched_units
  data = data[!match_index,]

  #while there are still covariates for matching

  while ((length(cur_covs) > 1) &&
         (sum(data[,'treated'] == 0) > 0) &&
         (sum(data[,'treated'] == 1) > 0)) {

    level = level + 1

    #Temporarily drop one covariate at a time to calculate Match Quality
    #Drop the covariate that returns highest Match Quality Score


    list_score <- unlist(lapply(cur_covs, match_quality_bit, data, holdout, num_covs, cur_covs, covs_max_list,
                                tradeoff, PE_function, model, ridge_reg, lasso_reg, tree_depth))
    quality <- max(list_score)

    cur_covs = cur_covs[-which(list_score == quality)] #Dropping one covariate
    covs_max_list = covs_max_list[-which(list_score == quality)]

    SCORE[[level-1]] <- quality
    covs_list[[level]] <- column[(cur_covs + 1)]

    # Update Match
    return_match = update_matched_bit(data, cur_covs, covs_max_list)
    match_index = return_match[[1]]
    index = return_match[[2]]

    # Set matched = num_covs and get those matched units
    data[match_index,'matched'] = length(cur_covs)
    return_df = rbind(return_df,data[match_index,])
    CATE[[level]] <- get_CATE_bit(data, match_index, index, cur_covs, covs_max_list, column, factor_level)

    # Remove matched_units
    data = data[!match_index,]
  }

  colnames(return_df) <- column
  rownames(return_df) <- NULL
  return_df[,1:num_covs] <- mapply(function(x,y) factor_level[x,][return_df[,y]+1], 1:num_covs, 1:num_covs)
  return(list(covs_list, CATE, unlist(SCORE), return_df))
}

#result_bit <- FLAME_bit(data, holdout, 15, covs_max_list, 0.1)

