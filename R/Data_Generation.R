#' Generate Synthetic Data
#'
#' \code{Data_Generation} generates synthetic data, where each covariate is a
#' binary variable.
#'
#' @param ncontrol number of samples in the control group
#' @param ntreated number of samples in the treated group
#' @param nimportant number of important covariates
#' @param ntrivial number of unimportant covariates
#' @param non_linear number of covariates with non-linear relationship
#' @param U coefficient of non-linear term
#' @return synthetic data
#' @examples
#' Data_Generation(10, 10, 10, 5, 5, 5)
#'
#' @export

Data_Generation <- function(ncontrol, ntreated, nimportant,
                            ntrivial, non_linear, U) {

  if (non_linear > nimportant) {
    stop("The number of covariates with non-linear relationship should be less than or equal to the number of important covariates.")
  }

  if (nimportant != 0) {
    # Generate important x_{i} for both control and treated, where each x_{i} is bernoulli(0.5)
    xc = sapply(rep(0.5,nimportant),function(p) rbinom(ncontrol,1,p))
    xt = sapply(rep(0.5,nimportant),function(p) rbinom(ntreated,1,p))

    # Generate s ~ Uniform(-1,1) and alpha_{i} ~ N(10s,1)
    sign = runif(nimportant,-1,1)
    alpha = sapply(sign, function(s) rnorm(1,10*s,1))

    # Compute outcome for control units
    yc = xc %*% alpha

    # Generate Beta_{i} for treated units, where Beta_{i} ~ N(1.5,0.15)
    treatment_eff_coef = rnorm(nimportant,1.5, 0.15)
    treatment_effect = xt %*% treatment_eff_coef

    # Generate nonlinear term for treated units
    treatment_effect_second = rep(0,ntreated)
    if (U != 0) {
      xt_second = xt[,1:non_linear]
      for (i in 1:(non_linear-1)) {
        for (j in (i+1):non_linear) {
          treatment_effect_second = treatment_effect_second + (xt_second[,i] * xt_second[,j])
        }
      }
    }
    treatment_effect_second = matrix(treatment_effect_second)

    # Compute outcome for treated units
    yt = xt %*% alpha + treatment_effect + U * treatment_effect_second

  }

  if (ntrivial != 0) {
    # Generate unimportant covariates for both control and treated
    xc2 = sapply(rep(0.5,ntrivial),function(p) rbinom(ncontrol,1,p))
    xt2 = sapply(rep(0.5,ntrivial),function(p) rbinom(ntreated,1,p))

    # Combine control and treated into dataframe
    df1 = cbind(xc,xc2,yc,rep(0,ncontrol))
    df2 = cbind(xt,xt2,yt,rep(1,ntreated))
  } else {
    # Combine control and treated into dataframe
    df1 = cbind(xc,yc,rep(0,ncontrol))
    df2 = cbind(xt,yt,rep(1,ntreated))
  }

  df = data.frame(rbind(df1,df2))
  colnames(df) <- c(paste("x",seq(1,nimportant + ntrivial), sep = ""),"outcome","treated")

  num_covs = ntrivial + nimportant

  # Convert each covariate and treated into type factor
  df[,c(1:num_covs,num_covs+2)] <- sapply(df[,c(1:num_covs,num_covs+2)], factor)
  df <- as.data.frame(unclass(df))

  return(df)
}




