# import the helper functions
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
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
