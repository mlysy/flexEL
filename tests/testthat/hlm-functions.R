#--- heteroscedastic linear model ------------------------------------------

#' Fit a heteroscedastic linear model.
#'
#' @param y Vector of response variables.
#' @param X Matrix of mean covariates.  Intercept must be included manually if desired.
#' @param W Matrix of log-variance covariates.  Intercept must be included manually if desired.
#' @param niter Maximum number of iterations.
#' @param tol Maximum relative tolerance.
#' @return A list containing elements:
#' \describe{
#'   \item{\code{beta}}{The MLE of beta.}
#'   \item{\code{gamma}}{The MLE of gamma.}
#'   \item{\code{iter}}{The number of iterations.}
#'   \item{\code{re}}{The maximum relative error between the last two iterations.}
#' }
#' @details The heteroscedastic linear model is given by
#' \deqn{
#' y_i \sim N( x_i' \beta, \exp(w_i' \gamma) ).
#' }
hlm.fit <- function(y, X, W, niter = 100, tol = 1e-5) {
  # preliminary calculations
  px <- ncol(X)
  pw <- ncol(W)
  n <- length(y)
  # relative error: 2 * |theta1-theta2|/|theta1+theta2|
  rel.err <- function(theta1, theta2) {
    abs(theta1-theta2)/abs((theta1+theta2)/2)
  }
  # initialize the algorithm
  beta0 <- coef(lm(y ~ X - 1))
  z <- (y - X %*% beta0)^2 # z ~ Gamma(1/2, 2 * exp(W %*% gamma))
  Mg <- glm(z ~ W - 1, family = Gamma(link = "log"))
  gamma0 <- coef(Mg)
  for(ii in 1:niter) {
    # update beta
    Wg <- predict(Mg, type = "link") # W %*% gamma
    # weighted regression: yi ~ N( xi' beta, s^2 / wi)
    Mb <- lm(y ~ X - 1, weights = exp(-Wg)) # y ~ N(X %*% beta, exp(Wg))
    beta1 <- coef(Mb)
    # update gamma
    z <- (y - predict(Mb))^2 # (y - X %*% beta)^2
    # gamma regression: zi ~ Gamma(1, exp(wi' gamma))
    Mg <- glm(z ~ W - 1, family = Gamma("log"))
    gamma1 <- coef(Mg)
    # old values: beta0, gamma0
    # new values: beta1, gamma1
    # calculate relative error
    re <- rel.err(c(beta0, gamma0), c(beta1, gamma1))
    # store new parameter values
    beta0 <- beta1
    gamma0 <- gamma1
    if(max(re) < tol) {
      # break out of for-loop
      break # faster than while-loop with ii <- ii+1 when niter is known.
    }
  }
  # check convergence
  if(max(re) > tol | ii == niter) {
    warning("Did not converge.  Try increasing number of iterations.")
  }
  # output
  list(beta = beta0, gamma = gamma0, iter = ii, re = max(re))
}

#' Loglikelihood of the heteroscedastic linear model.
#'
#' @param beta Vector of mean regression parameters.
#' @param gamma Vector of variance regression parameters.
#' @param y Vector of responses.
#' @param X Mean covariate matrix.
#' @param W Variance covariate matrix
#' @return The loglikelihood evaluated at \code{(beta, gamma)}.
hlm.loglik <- function(beta, gamma, y, X, W) {
  sum(dnorm(y, mean = X %*% beta,
            sd = exp(W %*% (.5*gamma)), log = TRUE))
}

#' Predictions with the heteroscedastic linear model.
#'
#' Calculates the conditional mean and prediction interval for responses \code{y} given the covariate matrices \code{X} and \code{W}.
#'
#' @param beta Vector of mean regression parameters.
#' @param gamma Vector of variance regression parameters.
#' @param X Mean covariate matrix.
#' @param W Variance covariate matrix.
#' @param level Coverage level of the prediction interval (e.g., 95%).
#' @return A 3-column matrix consisting of the conditional mean and lower/upper prediction bounds.
hlm.predict <- function(beta, gamma, X, W, level = .95) {
  yhat <- X %*% beta
  se <- - exp(W %*% (gamma/2)) * qnorm((1-level)/2)
  ans <- cbind(yhat, yhat - se, yhat + se)
  rownames(ans) <- 1:length(yhat)
  colnames(ans) <- c("fit", "lwr", "upr")
  ans
}
