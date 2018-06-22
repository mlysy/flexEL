# import the helper functions
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")

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
#'@param lambda is length-\cote{p} vector
#'@param omegas is length-\code{nObs} vector
#'@param gradmat is \code{nObs x p} matrix returned by gradMat
#'@return value of the gradient of the logEL evaluated at theta
logELgrad_R <- function(omegas, lambda, gradlist) {
  grad <- sapply(gradlist, function(x) lambda %*% x)
  return(c(-grad %*% omegas))
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

# --------------------------------------------------------------------------- 
source("../dontrun/gen_eps.R")

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

## Note: the gradient in HMC paper only works for theta in the support of EL
## so if the initial point is far away this may not work.

# ---- try nlm ----
negnorm <- function(x) {
  return(-dnorm(x, mean=1))
}
nlm(negnorm,1)
