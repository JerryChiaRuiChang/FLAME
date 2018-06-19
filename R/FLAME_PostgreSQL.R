#update_matched function takes list of covariates (cur_covs) to match
#and update column matched = 0 to matched = l (level) for matched units

update_matched_PostgreSQL <- function(cur_covs, level) {

  #Convert column names to dynamic strings

  covariates <- toString(sprintf("x%s", cur_covs, cur_covs))
  equalcovariates <- paste(sprintf("S.x%s = data.x%s", cur_covs, cur_covs), collapse = " AND ")
  level <- toString(level)

  #Update Data

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
                             AND matched = 0",covariates,covariates,level, covariates,equalcovariates)))
}

#get_CATE function takes list of covariates that are used to
#match at level l and return dataframe that includes
#(1) list of covariates that are used to match at level l
#(1) conditional average treatment effect (effect)
#(2) size of each matched group (size)

get_CATE_PostgreSQL <- function(cur_covs, level, column) {

  #Convert column names to dynamic strings

  covariates <- toString(sprintf("x%s", cur_covs, cur_covs))
  datacovariates <- toString(sprintf("control.x%s", cur_covs, cur_covs))
  equalcovariates <- paste(sprintf("control.x%s = treated.x%s", cur_covs, cur_covs), collapse = " AND ")

  #Get conditional average treatment effect

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
    covariates,level,covariates,covariates,level,covariates,datacovariates,equalcovariates)))

  #If the data frame to be returned is empty, convert its column names to covariates at current iteration
  #else, convert column names to back to its original column
  if (nrow(CATE) == 0) {
    CATE <- setNames(data.frame(matrix(ncol = length(cur_covs)+2, nrow = 0)),
                     c(column[(cur_covs + 1)],"effect","size"))
  } else {
    colnames(CATE) <- c(column[(cur_covs + 1)],"effect","size")
  }

  return(CATE)
}

#match_quality function takes holdout dataset, number of total covariates,
#list of current covariates, covariate c to temporily remove from, and trafeoff
#parameter as input. The function then computes Balancing Factor and Predictive Error,
#returning Match Quality.

match_quality_PostgreSQL <- function(c, holdout, num_covs, cur_covs,tradeoff) {

  #temporarly remove covariate c

  covs_to_match = cur_covs[cur_covs != c]

  #Convert column names to dynamic strings

  covariates <- toString(sprintf("x%s", covs_to_match, covs_to_match))
  equalcovariates <- paste(sprintf("S.x%s = data.x%s", covs_to_match, covs_to_match), collapse = " AND ")

  #get matched group for covariate list that exclude c

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

  match <- match[,-1] #Get rid of row.names
  dbWriteTable(db,"match",match, overwrite = TRUE) #write match dataframe into db

  #get unmatched group for covariate list that exclude c
  if (nrow(match) == 0) {
    unmatch <- dbGetQuery(db, "SELECT * FROM data WHERE matched = 0")
  } else {
  unmatch <- dbGetQuery(db, gsub("[[:space:]]{2,}"," ",
                                 sprintf("SELECT *
                                         FROM data
                                         WHERE matched = 0
                                         AND NOT EXISTS
                                         (SELECT *
                                         FROM match S
                                         WHERE %s)",
                                         equalcovariates)))}

  unmatch <- unmatch[,-1] #Get rid of row.names
  dbWriteTable(db,"unmatch",unmatch, overwrite = TRUE) # write unmatch dataframe into db

  #Get number of units in each following group to calculate Balancing Factor
  #(1) units that have been matched and belong to control group (match_control)
  #(2) units that have been matched and belong to treated group (match_treated)
  #(3) units that haven't been matched and belong to control group (unmatch_control)
  #(4) units that haven't been matched and belong to treated group (unmatch_treated)

  if (nrow(match) == 0) {
    match_control <- 0
    match_treated <- 0
  } else {
  match_control <- as.integer(dbGetQuery(db, "SELECT count(*) FROM match WHERE treated = 0")[1,1])
  match_treated <- as.integer(dbGetQuery(db, "SELECT count(*) FROM match WHERE treated = 1")[1,1])
  }

  if (nrow(unmatch) == 0) {
    unmatch_control <- 0
    unmatch_treated <- 0
  } else {
  unmatch_control <- as.integer(dbGetQuery(db, "SELECT count(*) FROM unmatch WHERE treated = 0")[1,1])
  unmatch_treated <- as.integer(dbGetQuery(db, "SELECT count(*) FROM unmatch WHERE treated = 1")[1,1])
  }

  #Run python script PE.py to get Predictive Error
  source_python(system.file("PE.py",package = "FLAME"))

  #source_python("PE.py")
  #Compute Predictive Error
  if (length(covs_to_match) == 1) {
    PE <- predictive_error(r_to_py(holdout),seq(0,num_covs - 1),list(covs_to_match))
  } else {
    PE <- predictive_error(r_to_py(holdout),seq(0,num_covs - 1),covs_to_match)
  }

  #If the unmatched group does not have any control/treated units then return PE

  if (unmatch_control == 0 | unmatch_treated == 0) {
    return(PE)
  } else {
    BF <- match_control/unmatch_control + match_treated/unmatch_treated #Compute Balancing Factor
    return(tradeoff * BF + PE)
  }
}

#'PostgreSQL database implementation
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
#'@param num_covs number of covariates
#'@param tradeoff tradeoff parameter to compute Matching Quality
#'@return (1) list of covariates used for matching at each iteration (2) list of
#'  dataframe showing all matched units, size of each matched group, and its
#'  conditional average treatment effect (CATE).
#'@import reticulate
#'@import RPostgreSQL
#'@export

FLAME_PostgreSQL <- function(db,data,holdout,num_covs,tradeoff) {

  data <- data.frame(data) #Convert input data to data.frame if not already converted
  holdout <- data.frame(holdout) #Convert holdout data to data.frame if not already converted
  # Convert each covariate and treated into type integer
  data[,c(1:num_covs,num_covs+2)] <- sapply(data[,c(1:num_covs,num_covs+2)],as.integer)
  holdout[,c(1:num_covs,num_covs+2)] <- sapply(holdout[,c(1:num_covs,num_covs+2)],as.integer)

  column <- colnames(data)

  #change input data and holdout training data column name
  colnames(data) <- c(paste("x",seq(0,num_covs-1), sep = ""),"outcome","treated")
  colnames(holdout) <- c(paste("x",seq(0,num_covs-1), sep = ""),"outcome","treated")

  data$matched <- 0 #add column matched to input data

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

  update_matched_PostgreSQL(cur_covs,level)
  covs_list[[level]] <- column[(cur_covs + 1)]
  CATE[[level]] <- get_CATE_PostgreSQL(cur_covs,level,column)


  #while there are still covariates for matching

  while ((length(cur_covs) > 1) &&
         (dbGetQuery(db, "select count(*) from data where matched = 0 and treated = 0")[1,1] > 0) &&
         (dbGetQuery(db, "select count(*) from data where matched = 0 and treated = 1")[1,1] > 0)) {

    level = level + 1

    #Temporarily drop one covariate at a time to calculate Match Quality
    #Drop the covariate that returns highest Match Quality Score

    #quality = -Inf
    #covs_to_drop = NULL

    #for (c in cur_covs) {
    #  score = match_quality_PostgreSQL(c,holdout, num_covs, cur_covs, tradeoff)
    #  if (score > quality) {
    #    quality = score
    #    covs_to_drop = c
    #  }
    #}

    list_score <- unlist(lapply(cur_covs,match_quality_PostgreSQL,holdout, num_covs, cur_covs, tradeoff))
    quality <- max(list_score)
    covs_to_drop <- cur_covs[which(list_score == quality)]

    cur_covs = cur_covs[cur_covs != covs_to_drop] #Dropping one covariate

    #Update Match
    SCORE[[level-1]] <- quality
    covs_list[[level]] <- column[(cur_covs + 1)]
    update_matched_PostgreSQL(cur_covs,level)
    CATE[[level]] <- get_CATE_PostgreSQL(cur_covs,level,column)
  }

  return_list <- NULL
  return_list[[1]] <- covs_list
  return_list[[2]] <- CATE

  return(return_list)
}



#data <- data.frame(Data_Generation(100,100,10,0))
#holdout <- data
#num_covs <- 10
#tradeoff <- 0.1

# Connecting to RPostgreSQL

#drv <- dbDriver('PostgreSQL')
#db <- dbConnect(drv, dbname="FLAME", host='localhost',
#             port=5432, user="postgres", password = 'new_password')

#res_Postgres <- FLAME_PostgreSQL(db,data,data,10,0.1)

#dbDisconnect(db)

#db <- dbConnect(SQLite(),"tempdb")
#FLAME::FLAME_SQLite(db,data,holdout,num_covs,0.1)
#dbDisconnect(db)

#FLAME::FLAME_bit(data,holdout,10,100,100,rep(2,10),0.1)

