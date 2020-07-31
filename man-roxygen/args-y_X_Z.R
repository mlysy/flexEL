#' @param y A numeric vector of resoponses of length \code{n_obs} where \code{n_obs} 
#'   is the number of observations.
#' @param X A numeric matrix of covariates in the location function of dimension 
#'   \code{n_obs} x \code{n_bet} where \code{n_obs} is the number of observations and 
#'   \code{b_bet} is the number of location function's coefficients (length of \code{beta}).
#' @param Z A numeric matrix of covariates in the scale function of dimension 
#'   \code{n-obs} x \code{n_gam} where \code{n_obs} is the number of observations and 
#'   \code{b_gam} is the number of scale function's coefficients (length of \code{gamma}).
