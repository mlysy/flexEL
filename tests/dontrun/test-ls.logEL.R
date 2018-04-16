# check that logEL plot under location-scale models
require(bayesEL)
source("../testthat/el-utils.R")

#---- mean reg: X and Z both dim 1----
# dimensions
n <- 100 # number of observations
numpoints <- 100

p <- 1
q <- 1
# X <- matrix(rep(1,n),n,p)
X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q)
Z <- matrix(rnorm(n),n,q)

beta0 <- rnorm(p)
gamma0 <- rnorm(q)
# beta0 <- 2
# gamma0 <- -0.01
theta0 <- c(beta0,gamma0)
y <- c(X %*% beta0 + exp(Z %*% gamma0)*rnorm(n))

plot(y,cex=0.3)
beta.seq <- seq(-.5+beta0, .5+beta0, length.out = numpoints)
gamma.seq <- seq(-.5+gamma0, .5+gamma0, length.out = numpoints)
theta.seq <- rbind(beta.seq, gamma.seq)

logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,beta.seq[ii],gamma0)
  logel.seq[1,ii] <- logEL(G)
  G <- mrls.evalG(y,X,Z,beta0,gamma.seq[ii])
  logel.seq[2,ii] <- logEL(G)
}
plot(beta.seq,logel.seq[1,], type='l')
abline(v=beta0,col='red')
plot(gamma.seq, logel.seq[2,], type='l')
abline(v=gamma0,col='red')

# solve by hlm: 
theta.hat <- hlm.fit(y = y, X = X, W = Z)
# quick check
rbind(true = beta0, est = theta.hat$beta) # beta
rbind(true = gamma0, est = theta.hat$gamma) # gamma

# optim check
loglik <- function(theta) {
  hlm.loglik(beta = theta[1:p], gamma = theta[p+1:q],
             y = y, X = X, W = Z)
}
theta.names <- c(paste0("beta[", 1:p-1, "]"), paste0("gamma[", 1:q-1, "]"))
theta.names <- parse(text = theta.names) # convert to greek symbols
optim_proj(fun = loglik,
           xsol = c(theta.hat$beta, theta.hat$gamma),
           xnames = theta.names)


# ---- mean reg: X and Z both dim 2 ---- 
# dimensions
n <- 100 # number of observations
numpoints <- 100

p <- 2
q <- 2
# X <- matrix(rep(1,n),n,p)
# X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q)
# Z <- matrix(rnorm(n),n,q)

X <- cbind(1,rnorm(n))
# Z <- cbind(1,rnorm(n))
Z <- cbind(rnorm(n),rnorm(n))

# beta0 <- rnorm(p)
# gamma0 <- rnorm(q)
beta0 <- c(1,2)
gamma0 <- c(-0.8, -0.2)
y <- c(X %*% beta0 + exp(Z %*% gamma0)*rnorm(n))
plot(y,cex=0.3)
beta1.seq <- seq(-.5+beta0[1], .5+beta0[1], length.out = numpoints)
beta2.seq <- seq(-.5+beta0[2], .5+beta0[2], length.out = numpoints)
gamma1.seq <- seq(-.5+gamma0[1], .5+gamma0[1], length.out = numpoints)
gamma2.seq <- seq(-.5+gamma0[2], .5+gamma0[2], length.out = numpoints)

logel.seq <- matrix(rep(NA,p*q*numpoints),p*q,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,c(beta1.seq[ii],beta0[2]),gamma0)
  logel.seq[1,ii] <- logEL(G)
  G <- mrls.evalG(y,X,Z,c(beta0[1],beta2.seq[ii]),gamma0)
  logel.seq[2,ii] <- logEL(G)
  G <- mrls.evalG(y,X,Z,beta0,c(gamma1.seq[ii],gamma0[2]))
  logel.seq[3,ii] <- logEL(G)
  G <- mrls.evalG(y,X,Z,beta0,c(gamma0[1],gamma2.seq[ii]))
  logel.seq[4,ii] <- logEL(G)
}
plot(beta1.seq,logel.seq[1,], type='l')
abline(v=beta0[1],col='red')
plot(beta2.seq,logel.seq[2,], type='l')
abline(v=beta0[2],col='red')
plot(gamma1.seq, logel.seq[3,], type='l')
abline(v=gamma0[1],col='red')
plot(gamma2.seq, logel.seq[4,], type='l')
abline(v=gamma0[2],col='red')

# ---- test sampler ---- 
source("hlm-functions.R")
library(optimCheck)
theta.hat <- hlm.fit(y = y, X = X, W = Z)

# quick check
rbind(true = beta0, est = theta.hat$beta) # beta
rbind(true = gamma0, est = theta.hat$gamma) # gamma

loglik <- function(theta) {
  hlm.loglik(beta = theta[1:p], gamma = theta[p+1:q],
             y = y, X = X, W = Z)
}

# optim check
theta.names <- c(paste0("beta[", 1:p-1, "]"),
                 paste0("gamma[", 1:q-1, "]"))
theta.names <- parse(text = theta.names) # convert to greek symbols

optim_proj(fun = loglik,
           xsol = c(theta.hat$beta, theta.hat$gamma),
           xnames = theta.names)

nsamples <- 20000
nburn <- 3000
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
sigs <- rep(0.1,length(betaInit)+length(gammaInit))

system.time(
  postout <- mrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sigs)
)
theta_chain <- postout$theta_chain
theta_accept <- postout$paccept
theta_accept
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(gamma[0]), main='')
hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')
