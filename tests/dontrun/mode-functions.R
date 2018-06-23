## Note: the gradient in HMC paper only works for theta in the support of EL
## so if the initial point is far away this may not work.

# ---- shared functions ----

#'
#'@param lambda is length-\cote{p} vector
#'@param omegas is length-\code{nObs} vector
#'@param gradmat is \code{nObs x p} matrix returned by gradMat
#'@return value of the gradient of the logEL evaluated at theta
logELgrad_R <- function(omegas, lambda, gradlist) {
  grad <- sapply(gradlist, function(x) lambda %*% x)
  return(c(-grad %*% omegas))
}

# ---- mr functions ----

#'
#'@param X is \code{nObs x p} matrix
#'@param y is length-\code{nObs} vector
#'@param theta is length-\cote{p} vector
#'@details ...
#'@return a list of gradient values of length \code{nObs}
mr.deltaG_R <- function(y, X, theta) {
  lx <- split(X, row(X))
  dg <- lapply(lx, function(x) tcrossprod(x,x))
  return(dg)
}

#'
#'@param X is \code{nObs x p} matrix
#'@param y is length-\code{nObs} vector
#'@param beta is length-\cote{p} vector
#'@details ...
#'@return
mr.neglogEL_R <- function(y, X, beta) {
  G <- mr.evalG_R(y, X, beta)
  lout <- lambdaNR_R(G)
  if (lout$conv) {
    lambda <- lout$lambda
    omegas <- c(1/(1-t(lambda) %*% t(G)) / sum(1/(1-t(lambda) %*% t(G))))
    res <- -sum(log(omegas)) # negative logEL
    gradlist <- mr.deltaG_R(y, X, beta)
    grad <- -logELgrad_R(omegas, lambda, gradlist) # negative gradient
    attr(res, "gradient") <- grad
  }
  else {
    # TODO: if not converged, what should be the gradient...??
    res <- Inf
    attr(res, "gradient") <- Inf
  }
  return(res)
}

# ---- qr functions ----

# smoothed indicator function 1(x <= 0)
ind.smooth_R <- function(x, s=100) {
  return(1/(1+exp(s*x)))
}

# 1st derivative of ind.smooth
ind1.smooth_R <- function(x, s=100) {
  return(-s*exp(s*x)/(1+exp(s*x))^2)
}

# 2nd derivative of ind.smooth
ind2.smooth_R <- function(x, s=100) {
  return((-s*s*exp(s*x)*(1+exp(s*x))+2*s*s*exp(2*s*x))/(1+exp(s*x))^3)
}

# smoothed check function
rho.smooth_R <- function(x, tau, s=100) {
  return(x*(tau-ind.smooth_R(x, s)))
}

# 1st derivate of rho_smooth
rho1.smooth_R <- function(x, tau, s=100) {
  return(tau-ind.smooth_R(x,s)-x*ind1.smooth_R(x,s))
}

# 2nd derivative of rho_smooth
rho2.smooth_R <- function(x, tau, s=100) {
  retval <- -2 * ind1.smooth_R(x,s) - x*ind2.smooth_R(x,s)
  return(retval)
}

# # generlized Huber function
# huber_gen <- function(x,delta) {
#   
# }

#'
#'@param X is \code{nObs x p} matrix
#'@param y is length-\code{nObs} vector
#'@param tau is a scalar of quantile percentage in (0,1)
#'@param beta is length-\cote{p} vector
#'@details ...
#'@return a list of gradient values of length \code{nObs}
qr.deltaG_R <- function(y, X, tau, beta, s=100) {
  lx <- split(X, row(X))
  dg <- mapply(function(x,y) {
    eps <- y - x %*% beta
    tcrossprod(rho2.smooth_R(eps,tau,s)[1]*x,x)
  }, lx, y, SIMPLIFY = FALSE)
  return(dg)
}

# smoothed version of qr.evalG
qr.evalG.smooth_R <- function(y,X,tau,beta,s) {
  tX <- t(X)
  yXb <- y - c(beta %*% tX)
  pyXb <- rho1.smooth_R(yXb,tau,s)
  G <- sweep(tX, MARGIN = 2, pyXb, `*`)
  return(t(G))
}

#'
#'@param X is \code{nObs x p} matrix
#'@param y is length-\code{nObs} vector
#'@param beta is length-\cote{p} vector
#'@details ...
#'@return
qr.neglogEL_R <- function(y, X, tau, beta, s=100) {
  G <- qr.evalG.smooth_R(y, X, tau, beta, s)
  lout <- lambdaNR_R(G)
  if (lout$conv) {
    lambda <- lout$lambda
    omegas <- c(1/(1-t(lambda) %*% t(G)) / sum(1/(1-t(lambda) %*% t(G))))
    res <- -sum(log(omegas)) # negative logEL
    gradlist <- qr.deltaG_R(y, X, tau, beta, s)
    grad <- -logELgrad_R(omegas, lambda, gradlist) # negative gradient
    attr(res, "gradient") <- grad
  }
  else {
    # TODO: if not converged, what should be the gradient...??
    message("qr.neglogEL_R: lambdaNR not coverged.")
    res <- Inf
    attr(res, "gradient") <- Inf
  }
  return(res)
}