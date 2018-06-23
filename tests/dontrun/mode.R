# import the helper functions
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")

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

# --------------------------------------------------------------------------- 
source("../dontrun/gen_eps.R")

## mr ##

# 1-d problem 
n <- 100
mu0 <- 1
X <- matrix(rep(1,n), n, 1) # each row of X is one observation
eps <- rnorm(n) # N(0,1) error term
y <- c(X * mu0) + eps

mr.neglogEL_R(y,X,mu0)

numpoints <- 100
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints)
logel.seq <- rep(NA,numpoints)
grad.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  temp <- -mr.neglogEL_R(y,X,mu.seq[ii])
  logel.seq[ii] <- temp
  grad.seq[ii] <- attributes(temp)$gradient
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, mean(y), expression(mu))
plot(grad.seq,type='l')
abline(h=0,col='blue')

nlm(mr.neglogEL_R,1.2,y=y,X=X)

# 2-d problem
n <- 200
p <- 2
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)
eps <- rnorm(n) # N(0,1) error term
beta_intercept <- 1
beta_slope <- 1.5
y <- 1 + c(X1 %*% beta_slope) + eps 
beta0 <- c(beta_intercept, beta_slope)

mr.neglogEL_R(y,X,beta0)

nlm(mr.neglogEL_R,beta0*1.05,y=y,X=X)

# 3-d problem
n <- 200
p <- 2
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(2*n),n,2)
X <- cbind(X0,X1)
eps <- rnorm(n) # N(0,1) error term
beta_intercept <- 1
beta_slope <- c(1.5,-1.5)
y <- 1 + c(X1 %*% beta_slope) + eps 
beta0 <- c(beta_intercept, beta_slope)

mr.neglogEL_R(y,X,beta0)

nlm(mr.neglogEL_R,beta0*1.05,y=y,X=X)

## qr ##

# 1-d problem 
n <- 200
mu0 <- 1
X <- matrix(rep(1,n), n, 1) # each row of X is one observation
eps <- rnorm(n) # N(0,1) error term
y <- c(X * mu0) + eps
tau <- 0.75

qr.neglogEL_R(y,X,tau,mu0)

numpoints <- 100
mu.seq <- seq(-.5+mu0+qnorm(tau),.5+mu0+qnorm(tau),length.out = numpoints)
logel.seq <- rep(NA,numpoints)
grad.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  temp <- -qr.neglogEL_R(y,X,tau,mu.seq[ii],s=10)
  logel.seq[ii] <- temp
  grad.seq[ii] <- attributes(temp)$gradient
}
logelmode <- plotEL(mu.seq, logel.seq, mu0+qnorm(tau), quantile(y,alpha), expression(mu))
plot(grad.seq,type='l')
abline(h=0,col='blue')

nlm(qr.neglogEL_R,1,y=y,X=X,tau=tau,s=10)

# 2-d problem
n <- 200
p <- 2
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)
eps <- rnorm(n) # N(0,1) error term
beta_intercept <- 0.5
beta_slope <- 1
y <- 1 + c(X1 %*% beta_slope) + eps 
beta0 <- c(beta_intercept, beta_slope)
tau <- 0.75

qr.neglogEL_R(y,X,tau,beta0,s = 10)

nlm(qr.neglogEL_R,beta0*1.05,y=y,X=X,tau=tau)

# 3-d problem
n <- 200
p <- 2
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(2*n),n,2)
X <- cbind(X0,X1)
eps <- rnorm(n) # N(0,1) error term
beta_intercept <- 1
beta_slope <- c(1.5,-1.5)
y <- 1 + c(X1 %*% beta_slope) + eps 
beta0 <- c(beta_intercept, beta_slope)
tau <- 0.75

qr.neglogEL_R(y,X,tau,beta0)

nlm(qr.neglogEL_R,beta0*1.05,y=y,X=X,tau=tau)

# ---- try nlm ----
negnorm <- function(x) {
  return(-dnorm(x, mean=1))
}
nlm(negnorm,1)
