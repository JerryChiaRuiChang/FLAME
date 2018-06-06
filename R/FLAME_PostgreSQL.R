#update_matched function takes list of covariates (cur_covs) to match
#and update column matched = 0 to matched = l (level) for matched units

update_matched <- function(cur_covs, level) {

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

get_CATE <- function(cur_covs, level) {

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

  return(CATE)
}

#match_quality function takes holdout dataset, number of total covariates,
#list of current covariates, covariate c to temporily remove from, and trafeoff
#parameter as input. The function then computes Balancing Factor and Predictive Error,
#returning Match Quality.

match_quality <- function(holdout, num_covs, cur_covs, c, tradeoff) {

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

#' FLAME: Database Method
#'
#' @param db Name of Database
#' @param data Data Frame
#' @param holdout Holdout Training Data
#' @param num_covs Number of Covariates
#' @param tradeoff Tradeoff Parameter to compute Match Quality
#' @return List of covariates matched in each iteration, CATE for matched units
#' @import reticulate
#' @import RSQLite
#' @export

FLAME_PostgreSQL <- function(db,data,holdout,num_covs,tradeoff) {

  #Connect to database

  #Change dataframe column name and write it to db
  colnames(data) <- c(paste("x",seq(0,num_covs-1), sep = ""),"outcome","treated","matched")
  dbWriteTable(db,"data",data, overwrite = TRUE) #Write dataframe to database

  #Set up return objects

  covs_list = list() #list of covariates for matching at each level
  CATE = list() #list of dataframe that calculates conditional average treatment effect at each level
  SCORE = list()

  #Initialize the current covariates to be all covariates and set level to 1

  cur_covs = seq(0,num_covs - 1)
  level = 1

  #Get matched units without dropping anything

  update_matched(cur_covs,level)
  covs_list[[level]] <- cur_covs
  CATE[[level]] <- get_CATE(cur_covs,level)


  #while there are still covariates for matching

  while ((length(cur_covs) > 1) &&
         (dbGetQuery(db, "select count(*) from data where matched = 0 and treated = 0")[1,1] > 0) &&
         (dbGetQuery(db, "select count(*) from data where matched = 0 and treated = 1")[1,1] > 0)) {

    level = level + 1

    #Temporarily drop one covariate at a time to calculate Match Quality
    #Drop the covariate that returns highest Match Quality Score

    quality = -Inf
    covs_to_drop = NULL

    for (c in cur_covs) {
      score = match_quality(holdout, num_covs, cur_covs, c, tradeoff)
      if (score > quality) {
        quality = score
        covs_to_drop = c
      }
    }

    cur_covs = cur_covs[cur_covs != covs_to_drop] #Dropping one covariate

    #Update Match
    SCORE[[level-1]] <- quality
    covs_list[[level]] <- cur_covs
    update_matched(cur_covs,level)
    CATE[[level]] <- get_CATE(cur_covs,level)

  }

  return_list <- NULL
  return_list[[1]] <- covs_list
  return_list[[2]] <- CATE

  return(return_list)
}



#data <- data.frame(FLAME::Data_Generation(100,100,10,0))
#holdout <- data
#num_covs <- 10
#tradeoff <- 0.1


# Connecting to RPostgreSQL

#drv <- dbDriver('PostgreSQL')
#db<- dbConnect(drv, dbname="FLAME", host='localhost',
#               port=5432, user="postgres", password = 'new_password')

#dbDisconnect(db)
