#update_matched function takes list of covariates (cur_covs) to match
#and update column matched = 0 to matched = l (level) for matched units

update_matched_PostgreSQL <- function(db, cur_covs, compute_var) {

  #Convert column names to dynamic strings

  covariates <- toString(sprintf("x%s", cur_covs, cur_covs))
  equalcovariates <- paste(sprintf("S.x%s = data.x%s", cur_covs, cur_covs), collapse = " AND ")

  #Update Data

  if (compute_var) {
    dbExecute(db, gsub("[[:space:]]{2,}"," ",
                       sprintf("WITH tempgroups AS
                               (SELECT %s
                               FROM data
                               WHERE matched = 0
                               GROUP BY %s
                               HAVING SUM(treated) >= 2 AND COUNT(*) >= SUM(treated) + 2)
                               UPDATE data
                               SET matched = %s
                               WHERE EXISTS
                               (SELECT %s
                               FROM tempgroups S
                               WHERE %s)
                               AND matched = 0",covariates,covariates,length(cur_covs), covariates,equalcovariates)))
  }
  else {
    dbExecute(db, gsub("[[:space:]]{2,}"," ",
                       sprintf("WITH tempgroups AS
                               (SELECT %s
                               FROM data
                               WHERE matched = 0
                               GROUP BY %s
                               HAVING SUM(treated) > 0 AND SUM(treated) < COUNT(*))
                               UPDATE data
                               SET matched = %s
                               WHERE EXISTS
                               (SELECT %s
                               FROM tempgroups S
                               WHERE %s)
                               AND matched = 0",covariates,covariates,length(cur_covs), covariates,equalcovariates)))
  }

  num_matched <- as.integer(dbGetQuery(db, sprintf("SELECT count(*) FROM data WHERE matched = %s", length(cur_covs)))[1,1])
  message(paste("number of matched units = ", num_matched))
}

#get_CATE function takes list of covariates that are used to
#match at level l and return dataframe that includes
#(1) list of covariates that are used to match at level l
#(1) conditional average treatment effect (effect)
#(2) size of each matched group (size)

get_CATE_PostgreSQL <- function(db, cur_covs, column, factor_level, compute_var) {

  #Convert column names to dynamic strings

  covariates <- toString(sprintf("x%s", cur_covs, cur_covs))
  datacovariates <- toString(sprintf("control.x%s", cur_covs, cur_covs))
  equalcovariates <- paste(sprintf("control.x%s = treated.x%s", cur_covs, cur_covs), collapse = " AND ")

  #Get conditional average treatment effect

  if (compute_var) {

    CATE <- dbGetQuery(db, gsub("[[:space:]]{2,}"," ", sprintf(
      "WITH control AS
      (SELECT %s, AVG(outcome) AS conout, count(*) AS conc, VARIANCE(outcome) AS convar
      FROM data
      WHERE matched = %s AND treated = 0
      GROUP BY %s),
      treated AS
      (SELECT %s, AVG(outcome) AS treatout, count(*) AS treatc, VARIANCE(outcome) AS treatvar
      FROM data
      WHERE matched = %s AND treated = 1
      GROUP BY %s)
      SELECT %s, (treatout - conout) AS effect, (treatc + conc) AS size, (convar + treatvar) AS variance
      FROM
      (control INNER JOIN treated
      ON %s)",
      covariates,length(cur_covs),covariates,covariates,length(cur_covs),covariates,datacovariates,equalcovariates)))
  }

  else {
    CATE <- dbGetQuery(db, gsub("[[:space:]]{2,}"," ", sprintf(
      "WITH control AS
      (SELECT %s, AVG(outcome) AS conout, count(*) AS conc
      FROM data
      WHERE matched = %s AND treated = 0
      GROUP BY %s),
      treated AS
      (SELECT %s, AVG(outcome) AS treatout, count(*) AS treatc
      FROM data
      WHERE matched = %s AND treated = 1
      GROUP BY %s)
      SELECT %s, (treatout - conout) AS effect, (treatc + conc) AS size
      FROM
      (control INNER JOIN treated
      ON %s)",
      covariates,length(cur_covs),covariates,covariates,length(cur_covs),covariates,datacovariates,equalcovariates)))
  }

  if (compute_var) {
    #If the data frame to be returned is empty, convert its column names to covariates at current iteration
    #else, convert column names to back to its original column
    if (nrow(CATE) == 0) {
      CATE <- setNames(data.frame(matrix(ncol = length(cur_covs)+3, nrow = 0)),
                       c(column[(cur_covs + 1)],"effect","size", "variance"))
    } else {
      CATE <- data.frame(data.matrix(CATE)) # convert all columns into numeric
      CATE[,1:length(cur_covs)] <- mapply(function(x,y) factor_level[[x]][CATE[,y]], cur_covs + 1, 1:length(cur_covs))
      colnames(CATE) <- c(column[(cur_covs + 1)],"effect","size", "variance")
      CATE <- CATE[order(CATE$effect),]
      rownames(CATE) = NULL
    }
  }
  else {
    #If the data frame to be returned is empty, convert its column names to covariates at current iteration
    #else, convert column names to back to its original column
    if (nrow(CATE) == 0) {
      CATE <- setNames(data.frame(matrix(ncol = length(cur_covs)+2, nrow = 0)),
                       c(column[(cur_covs + 1)],"effect","size"))
    } else {
      CATE <- data.frame(data.matrix(CATE)) # convert all columns into numeric
      CATE[,1:length(cur_covs)] <- mapply(function(x,y) factor_level[[x]][CATE[,y]], cur_covs + 1, 1:length(cur_covs))
      colnames(CATE) <- c(column[(cur_covs + 1)],"effect","size")
      CATE <- CATE[order(CATE$effect),]
      rownames(CATE) = NULL
    }
  }
  return(CATE)
}

Regression_PE_PostgreSQL <- function(holdout_trt, holdout_ctl) {

  # MSE for treated

  model_lm <- lm(outcome ~ ., data = holdout_trt) # fit the data to lm model
  MSE_treated <- sum((holdout_trt$outcome - model_lm$fitted.values)^2)/length(holdout_trt$outcome) # compute mean squared error

  # MSE for control
  model_lm <- lm(outcome ~ ., data = holdout_ctl) # fit the data to lm model
  MSE_control <- sum((holdout_ctl$outcome - model_lm$fitted.values)^2)/length(holdout_ctl$outcome) # compute mean squared error

  return(MSE_treated + MSE_control)
}

#match_quality function takes holdout dataset, number of total covariates,
#list of current covariates, covariate c to temporily remove from, and trafeoff
#parameter as input. The function then computes Balancing Factor and Predictive Error,
#returning Match Quality.

match_quality_PostgreSQL <- function(c, db, holdout, num_covs, cur_covs, tradeoff,
                                     PE_function, model, ridge_reg, lasso_reg, tree_depth, compute_var, py_run) {

  #temporarly remove covariate c

  covs_to_match = cur_covs[cur_covs != c]

  #Convert column names to dynamic strings

  covariates <- toString(sprintf("x%s", covs_to_match, covs_to_match))
  equalcovariates <- paste(sprintf("S.x%s = data.x%s", covs_to_match, covs_to_match), collapse = " AND ")

  # Calculate number of units unmatched (available)

  num_control <- as.integer(dbGetQuery(db, "SELECT count(*) FROM data WHERE matched = 0 AND treated = 0")[1,1])
  num_treated <- as.integer(dbGetQuery(db, "SELECT count(*) FROM data WHERE matched = 0 AND treated = 1")[1,1])


  #get matched group for covariate list that exclude c

  if (compute_var) {

    match <- dbGetQuery(db, gsub("[[:space:]]{2,}"," ",
                                 sprintf("WITH tempgroups AS
                                         (SELECT %s
                                         FROM data
                                         WHERE matched = 0
                                         GROUP BY %s
                                         HAVING SUM(treated) >= 2 AND COUNT(*) >= SUM(treated) + 2)
                                         SELECT *
                                         FROM data
                                         WHERE EXISTS
                                         (SELECT *
                                         FROM tempgroups S
                                         WHERE %s)
                                         AND matched = 0",
                                         covariates,covariates,equalcovariates)))
  }

  else {

    match <- dbGetQuery(db, gsub("[[:space:]]{2,}"," ",
                                 sprintf("WITH tempgroups AS
                                         (SELECT %s
                                         FROM data
                                         WHERE matched = 0
                                         GROUP BY %s
                                         HAVING SUM(treated) > 0 AND SUM(treated) < COUNT(*))
                                         SELECT *
                                         FROM data
                                         WHERE EXISTS
                                         (SELECT *
                                         FROM tempgroups S
                                         WHERE %s)
                                         AND matched = 0",
                                         covariates,covariates,equalcovariates)))
  }

  match <- match[,-1] #Get rid of row.names
  dbWriteTable(db,"match",match, overwrite = TRUE) #write match dataframe into db

  if (nrow(match) == 0) {
    num_control_matched <- 0
    num_treated_matched <- 0
  } else {
    # Number of matched units
    num_control_matched <- as.integer(dbGetQuery(db, "SELECT count(*) FROM match WHERE treated = 0")[1,1])
    num_treated_matched <- as.integer(dbGetQuery(db, "SELECT count(*) FROM match WHERE treated = 1")[1,1])
  }

  #Compute Predictive Error

  if (!py_run) {
    holdout_trt <- holdout[holdout[,'treated'] == 1,]
    holdout_trt <- holdout_trt[,!(names(holdout_trt) %in% 'treated')]
    holdout_ctl <- holdout[holdout[,'treated'] == 0,]
    holdout_ctl <- holdout_trt[,!(names(holdout_ctl) %in% 'treated')]
    PE <- Regression_PE_PostgreSQL(holdout_trt, holdout_ctl)
  }

  if (!is.null(PE_function)) {
    # Compute -PE based on user defined PE_function
    outcome_treated <- holdout[holdout[,'treated'] == 1,][,'outcome']
    outcome_control <- holdout[holdout[,'treated'] == 0,][,'outcome']
    covs_treated <- as.matrix(holdout[holdout[,'treated'] == 1,][,covs_to_match + 1])
    covs_control <- as.matrix(holdout[holdout[,'treated'] == 0,][,covs_to_match + 1])
    PE <- -PE_function(outcome_treated, outcome_control, covs_treated, covs_control)
  }

  else {
    predictive_error  <- NULL
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

  #If the unmatched group does not have any control/treated units then return PE

  if (num_control == 0 | num_treated == 0) {
    return(PE)
  } else {
    BF <- num_control_matched/num_control + num_treated_matched/num_treated #Compute Balancing Factor
    return(tradeoff * BF + PE)
  }
}

#'PostgreSQL Database Implementation
#'
#'\code{FLAME_PostgreSQL} applies the FLAME algorithm based on PostgreSQL. If
#'your computer system does not have PostgreSQL installed, install from
#'\href{https://www.postgresql.org/download/}{here}. For setup of PostgreSQL
#'server, please refer to this
#'\href{http://www.postgresqltutorial.com/connect-to-postgresql-database/}{tutorial}.
#'User must connect to PostgreSQL server in R using the command
#'\code{dbConnect(dbDriver('PostgreSQL'), dbname="your_dbname",
#'host='your_localhost', port='your_port', user='your_username', password =
#''your_password')} and name the connection as \strong{db}
#'
#'@param db name of the database connection (\strong{must name the connection as
#'  db})
#'@param data input data
#'@param holdout holdout training data
#'@param compute_var indicator variable of computing variance (optional, default = FALSE)
#'@param tradeoff tradeoff parameter to compute Matching Quality (optional, default = 0.1)
#'@param PE_function user defined function to compute predivtive error
#'  (optional)
#'@param model user defined model - Linear, Ridge, Lasso, or DecisionTree
#'  (optional)
#'@param ridge_reg L2 regularization parameter if model = Ridge (optional)
#'@param lasso_reg L1 regularization parameter if model = Lasso (optional)
#'@param tree_depth maximum depth of decision tree if model = DecisionTree
#'  (optional)
#'@return (1) list of covariates FLAME performs matching at each iteration, (2)
#'  list of dataframe showing matched groups' sizes, conditional average treatment
#' effects (CATEs), and variance (if compute_var = TRUE) at each iteration, (3) matching quality at each
#'  iteration, and (4) the original data with additional column *matched*,
#'  indicating the number of covariates each unit is matched. If a unit is never
#'  matched, then *matched* will be 0.
#'@examples
#'\dontrun{
#'drv <- dbDriver('PostgreSQL')
#'
#'db <- dbConnect(drv, dbname="FLAME", host='localhost',
#'port=5432, user="postgres", password = 'new_password')
#'
#'FLAME_PostgreSQL(db = db, data = data, holdout = holdout)
#'
#'dbDisconnect(db)
#'}
#'@import RPostgreSQL
#'@import reticulate
#'@importFrom graphics boxplot
#'@importFrom stats rbinom rnorm runif setNames
#'@export

FLAME_PostgreSQL <- function(db, data, holdout, compute_var = FALSE, tradeoff = 0.1, PE_function = NULL,
                             model = NULL, ridge_reg = NULL, lasso_reg = NULL, tree_depth = NULL) {

  num_covs = ncol(data) - 2

  # If covariate(s) are not factor(s), then stop
  if (Reduce("|", sapply(1:num_covs, function(x) !is.factor(data[,x] ))) |
      Reduce("|", sapply(1:num_covs, function(x) !is.factor(holdout[,x] )))) {
    stop("Covariates are not factor data type")
  }

  # If treatment is not factor, then stop
  if (!is.factor(data[,num_covs + 2]) | !is.factor(holdout[,num_covs + 2])) {
    stop("Treatment variable is not factor data type")
  }

  # If outcome variable is not numeric, then stop
  if (!is.numeric(data[,num_covs + 1]) | !is.numeric(holdout[,num_covs + 1])) {
    stop("Outcome variable is not numeric data type")
  }

  if (!py_module_available("pandas")) {
    py_install("pandas")
    if (!py_module_available("pandas")) {
      warning("The package will use R’s default linear regression pandas module is not available. This will be VERY SLOW!
            For more information on how to attach Python module to R, please refer to https://rstudio.github.io/reticulate/reference/import.html.")
    }
  }

  if (!py_module_available("numpy")) {
    py_install("numpy")
    if (!py_module_available("numpy")) {
      warning("The package will use R’s default linear regression numpy module is not available. This will be VERY SLOW!
              For more information on how to attach Python module to R, please refer to https://rstudio.github.io/reticulate/reference/import.html.")
    }
  }

  if (!py_module_available("sklearn")) {
    py_install("sklearn")
    if (!py_module_available("sklearn")) {
      warning("The package will use R’s default linear regression numpy module is not available. This will be VERY SLOW!
              For more information on how to attach Python module to R, please refer to https://rstudio.github.io/reticulate/reference/import.html.")
    }
  }

  py_run = py_module_available("sklearn") && py_module_available("pandas") && py_module_available("numpy")

  #add column matched to input data
  data$matched <- as.integer(0)
  column <- colnames(data)

  factor_level <- lapply(data[,1:num_covs],levels) # Get levels of each factor

  # Convert each covariate and treated into type integer
  data[,c(1:num_covs)] <- sapply(data[,c(1:num_covs)], function(x) as.integer(x))
  data[,num_covs + 2] <- as.integer(levels(data[,num_covs+2])[data[,num_covs+2]])

  holdout[,c(1:num_covs)] <- sapply(holdout[,c(1:num_covs)],function(x) as.integer(x))
  holdout[,num_covs + 2] <- as.integer(levels(holdout[,num_covs+2])[holdout[,num_covs+2]])

  #change input data and holdout training data column name
  colnames(data) <- c(paste("x",seq(0,num_covs-1), sep = ""),"outcome","treated","matched")
  colnames(holdout) <- c(paste("x",seq(0,num_covs-1), sep = ""),"outcome","treated")

  #Write input data to database
  dbWriteTable(db,"data",data, overwrite = TRUE)

  #Set up return objects

  covs_list = list() #list of covariates for matching at each level
  CATE = list() #list of dataframe that calculates conditional average treatment effect at each level
  SCORE = list()

  #Initialize the current covariates to be all covariates and set level to 1

  cur_covs = seq(0,num_covs - 1)
  level = 1
  #Get matched units without dropping anything

  update_matched_PostgreSQL(db, cur_covs, compute_var)
  covs_list[[level]] <- column[(cur_covs + 1)]
  CATE[[level]] <- get_CATE_PostgreSQL(db, cur_covs, column, factor_level, compute_var)

  #while there are still covariates for matching

  while ((length(cur_covs) > 1) &&
         (dbGetQuery(db, "select count(*) from data where matched = 0 and treated = 0")[1,1] > 0) &&
         (dbGetQuery(db, "select count(*) from data where matched = 0 and treated = 1")[1,1] > 0)) {

    level = level + 1

    #Temporarily drop one covariate at a time to calculate Match Quality
    #Drop the covariate that returns highest Match Quality Score

    list_score <- unlist(lapply(cur_covs, match_quality_PostgreSQL, db, holdout, num_covs, cur_covs, tradeoff,
                                PE_function, model, ridge_reg, lasso_reg, tree_depth, compute_var, py_run))
    quality <- max(list_score)
    covs_to_drop <- cur_covs[which(list_score == quality)]

    cur_covs = cur_covs[! cur_covs %in% covs_to_drop]  #Dropping covariate(s)

    #Update Match
    SCORE[[level-1]] <- quality
    covs_list[[level]] <- column[(cur_covs + 1)]
    update_matched_PostgreSQL(db, cur_covs, compute_var)
    CATE[[level]] <- get_CATE_PostgreSQL(db, cur_covs, column, factor_level, compute_var)
  }


  return_df <- dbGetQuery(db, "SELECT * FROM data")[,-1]
  return_df[,1:num_covs] <- mapply(function(x,y) factor_level[[x]][return_df[,y]], 1:num_covs, 1:num_covs)
  colnames(return_df) <- column

  return_list = list(covs_list, CATE, unlist(SCORE), return_df)
  names(return_list) = c("covariate_lst", "CATE", "match_quality", "matched_data")

  return(return_list)
}


#data <- read.csv("/Users/Jerry/Desktop/flame_bit_breaks_on_this.csv")
#data[,c(1:22,24)] <- lapply(data[,c(1:22,24)], factor)
#holdout <- data

#drv <- dbDriver('PostgreSQL')
#db <- dbConnect(drv, dbname="FLAME", host='localhost',
#             port=5432, user="postgres", password = 'new_password')

#result_Postgres <- FLAME_PostgreSQL(db = db, data = data, holdout = holdout, compute_var = FALSE)
#dbDisconnect(db)

#compute_var = FALSE
#tradeoff = 0.1
#PE_function = NULL
#model = NULL
#ridge_reg = NULL
#lasso_reg = NULL
#tree_depth = NULL





