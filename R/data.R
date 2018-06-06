#' 2010 Natality Dataset in the U.S.A.
#'
#' This dataset contains 2,190,410 units of pregnancy mother,
#' where each unit includes 86 attributes of a mother's conditions
#' during pregnancy and her health history. Treated units are mothers
#' who smoke cigarettes during pregnancy. 204,886 are treated units,
#' and 1,985,524 are control units. The outcome is whether the infant
#' is diagnosed with an abnormal condition after birth.
#'
#' @format A data frame with 2190410  rows and 86 covariates:
#' \describe{
#'   \item{X0}{Covariate}
#'   \item{X1}{Covariate}
#'   ...
#'   \item{outcome}{Outcome Variable}
#'   \item{treated}{Treated Unit = 1, Control Unit = 0}
#'   \item{matched}{Whether an unit is matched}
#' }
"natality"
