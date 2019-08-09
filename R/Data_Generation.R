#' Generate Synthetic Data
#'
#' \code{Data_Generation} generates synthetic data, where each covariate is a
#' binary variable.
#'
#' @param num_control number of samples in the control group
#' @param num_treated number of samples in the treated group
#' @param num_cov_dense number of important covariates
#' @param num_cov_unimportant number of unimportant covariates
#' @param non_linear number of non-linear covariates
#' @param U non-linear term coefficient
#' @return synthetic data
#' @examples
#' Data_Generation(10, 10, 10, 5, 5)
#' Data_Generation(10, 10, 5, 10, 5)
#'
#' @export

Data_Generation <- function(num_control, num_treated,
                            num_cov_dense, num_cov_unimportant, non_linear, U) {

  if (num_cov_dense != 0) {
    # Generate important x_{i} for both control and treated, where each x_{i} is bernoulli(0.5)
    xc = sapply(rep(0.5,num_cov_dense),function(p) rbinom(num_control,1,p))
    xt = sapply(rep(0.5,num_cov_dense),function(p) rbinom(num_treated,1,p))

    # Generate s ~ Uniform(-1,1) and alpha_{i} ~ N(10s,1)
    dense_sign = runif(num_cov_dense,-1,1)
    dense = sapply(dense_sign, function(s) rnorm(1,10*s,1))

    # Compute outcome for control units
    yc = xc %*% dense

    # Generate Beta_{i} for treated units, where Beta_{i} ~ N(1.5,0.15)
    treatment_eff_coef = rnorm(num_cov_dense,1.5, 0.15)
    treatment_effect = xt %*% treatment_eff_coef

    # Generate nonlinear term for treated units
    treatment_effect_second = rep(0,num_treated)
    if (U != 0) {
      xt_second = xt[,1:non_linear]
      for (i in 1:non_linear-1)
        for (j in (i+1):non_linear) {
          treatment_effect_second = treatment_effect_second + (xt_second[,i] * xt_second[,j])
        }
    }
    treatment_effect_second = matrix(treatment_effect_second)

    # Compute outcome for treated units
    yt = xt %*% dense + treatment_effect + U * treatment_effect_second

  }

  if (num_cov_unimportant != 0) {
    # Generate unimportant covariates for both control and treated
    xc2 = sapply(rep(0.5,num_cov_unimportant),function(p) rbinom(num_control,1,p))
    xt2 = sapply(rep(0.5,num_cov_unimportant),function(p) rbinom(num_treated,1,p))

    # Combine control and treated into dataframe
    df1 = cbind(xc,xc2,yc,rep(0,num_control))
    df2 = cbind(xt,xt2,yt,rep(1,num_treated))
  } else {
    # Combine control and treated into dataframe
    df1 = cbind(xc,yc,rep(0,num_control))
    df2 = cbind(xt,yt,rep(1,num_treated))
  }

  df = data.frame(rbind(df1,df2))
  colnames(df) <- c(paste("x",seq(1,num_cov_dense + num_cov_unimportant), sep = ""),"outcome","treated")

  num_covs = num_cov_dense + num_cov_unimportant

  # Convert each covariate and treated into type factor
  df[,c(1:num_covs,num_covs+2)] <- sapply(df[,c(1:num_covs,num_covs+2)], factor)
  df <- as.data.frame(unclass(df))

  return(df)
}


