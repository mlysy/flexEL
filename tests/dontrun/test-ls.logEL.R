# check that logEL plot under location-scale models
require(bayesEL)
require(optimCheck)
source("hlm-functions.R")
source("../testthat/el-utils.R")

#---- mean reg: X and Z both dim 1----
# dimensions
n <- 100 # number of observations
numpoints <- 100

p <- 1
q <- 1
X <- matrix(rep(1,n),n,p)
# X <- matrix(rnorm(n),n,p)
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

# calculate marginal posterior
Theta.seq <- as.matrix(expand.grid(beta.seq, gamma.seq))
logel.mat <- apply(Theta.seq, 1, function(bb) {
  G <- mrls.evalG(y,X,Z,bb[1],bb[2])
  logEL(G)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

# mcmc
nsamples <- 20000
nburn <- 3000
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
sigs <- c(0.25,0.23)

system.time(
  postout <- mrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sigs)
)
theta_chain <- postout$Theta_chain
theta_accept <- postout$paccept
theta_accept
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
lines(beta.seq, norm_pdf(logel.marg[,1], beta.seq),
      cex=0.1, col = 'red', type='l')
hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
lines(gamma.seq, norm_pdf(logel.marg[,2], beta.seq),
      cex=0.1, col = 'red', type='l')

# ---- mean reg: X and Z both dim 2 ---- 
# dimensions
n <- 100 # number of observations
numpoints <- 100

p <- 2
q <- 2
X <- cbind(1,rnorm(n))
# X <- matrix(rnorm(p*n),n,p)
Z <- matrix(rnorm(q*n),n,q)
# Z <- cbind(1,rnorm(n))

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

# solve by hlm
theta.hat <- hlm.fit(y = y, X = X, W = Z)
# quick check
rbind(true = beta0, est = theta.hat$beta) # beta
rbind(true = gamma0, est = theta.hat$gamma) # gamma

loglik <- function(theta) {
  hlm.loglik(beta = theta[1:p], gamma = theta[p+1:q],
             y = y, X = X, W = Z)
}

# optim check
theta.names <- c(paste0("beta[", 1:p-1, "]"), paste0("gamma[", 1:q-1, "]"))
theta.names <- parse(text = theta.names) # convert to greek symbols
optim_proj(fun = loglik,
           xsol = c(theta.hat$beta, theta.hat$gamma),
           xnames = theta.names)

# mcmc
nsamples <- 20000
nburn <- 3000
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
sigs <- c(0.15,0.1,0.23,0.23)

system.time(
  postout <- mrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sigs)
)
theta_chain <- postout$Theta_chain
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

# --------------------------------------------------------------------------- 
# test qr.evalG and qrls.evalG for multiple quantile levels
n <- 20
alphas <- c(0.25,0.5)
Beta <- matrix(c(0.5,1,1,1),2,2)
x <- rnorm(n)
y <- 1 + x + rnorm(n)
X <- cbind(1,x)
qr.evalG(y,X,alphas,Beta)
qr.evalG(y,X,alphas[1],Beta[,1])
qr.evalG(y,X,alphas[2],Beta[,2])

# Beta <- cbind(rnorm(p),rnorm(p))
# Gamma <- cbind(rnorm(q),rnorm(q))
# t(qrls.evalG(y,X,Z,alphas,Beta,Gamma))
# --------------------------------------------------------------------------- 

# ---- quant reg: X and Z both dim 1 with 1 quant level ----
# dimensions
n <- 100 # number of observations
numpoints <- 100

p <- 1
q <- 1
X <- matrix(rep(1,n),n,p)
# X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q)
Z <- matrix(rnorm(n),n,q)

beta0 <- rnorm(p)
gamma0 <- rnorm(q)
# beta0 <- 2
# gamma0 <- -0.01
theta0 <- c(beta0,gamma0)
y <- c(X %*% beta0 + exp(Z %*% gamma0)*rnorm(n))
alpha <- 0.5

plot(y,cex=0.3)
beta.seq <- seq(-.5+beta0, .5+beta0, length.out = numpoints)
gamma.seq <- seq(-.5+gamma0, .5+gamma0, length.out = numpoints)
theta.seq <- rbind(beta.seq, gamma.seq)

logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta.seq[ii],gamma0)
  logel.seq[1,ii] <- logEL(G)
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma.seq[ii])
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

# calculate marginal posterior
Theta.seq <- as.matrix(expand.grid(beta.seq, gamma.seq))
logel.mat <- apply(Theta.seq, 1, function(bb) {
  G <- qrls.evalG(y,X,Z,alpha,bb[1],bb[2])
  logEL(G)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

# mcmc
nsamples <- 20000
nburn <- 3000
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
sigs <- c(0.45,0.3)

system.time(
  postout <- qrls.post(y,X,Z,alpha,nsamples,nburn,betaInit,gammaInit,sigs)
)
theta_chain <- postout$Theta_chain
theta_accept <- postout$paccept
theta_accept
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
lines(beta.seq, norm_pdf(logel.marg[,1], beta.seq),
      cex=0.1, col = 'red', type='l')
abline(v=theta0[1], col='black')
abline(v=mean(theta_chain[1,]), col='blue')
legend('topright',legend=c(expression('true value'),
                           expression('grid plot'),
                           expression('sample mean')),
       lty = c(1,1,1), col = c('black','red','blue'), cex = 0.6)
hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma[0]), main='')
lines(gamma.seq, norm_pdf(logel.marg[,2], beta.seq),
      cex=0.1, col = 'red', type='l')
abline(v=theta0[2], col='black')
abline(v=mean(theta_chain[2,]), col='blue')
legend('topright',legend=c(expression('true value'),
                           expression('grid plot'),
                           expression('sample mean')),
       lty = c(1,1,1), col = c('black','red','blue'), cex = 0.6)

# TODO: not done yet
# ---- quant reg: X and Z both dim 1 with 2 quant level ----
library(quantreg)
# dimensions
n <- 100 # number of observations
numpoints <- 100

p <- 1
q <- 1
# X <- matrix(rep(1,n),n,p)
X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q)
Z <- matrix(rnorm(n),n,q)

# beta0 <- rnorm(p)
# gamma0 <- rnorm(q)
beta0 <- 2
gamma0 <- -0.01
theta0 <- c(beta0,gamma0)
y <- c(X %*% beta0 + exp(Z %*% gamma0)*rnorm(n))
alpha <- c(0.25,0.75)

theta.hat <- hlm.fit(y,X,Z)
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
BetaInit <- cbind(betaInit,betaInit*1.01)
GammaInit <- cbind(gammaInit,gammaInit*0.99)
Sigs <- rbind(c(0.2,0.2),c(0.2,0.2))
nsamples <- 20000
nburn <- 3000

system.time(
  postout <- qrls.post(y,X,Z,alpha,nsamples,nburn,BetaInit,GammaInit,Sigs)
)
theta_chain <- postout$Theta_chain
paccept <- postout$paccept
paccept
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma[0]), main='')
hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')

# G <- qrls.evalG(y,X,Z,alpha,BetaInit,GammaInit)
# G1 <- qrls.evalG(y,X,Z,alpha[1],BetaInit[,1],GammaInit[,1])
# G2 <- qrls.evalG(y,X,Z,alpha[2],BetaInit[,2],GammaInit[,2])

postout1 <- qrls.post(y,X,Z,alpha[1],nsamples,nburn,BetaInit[,1],GammaInit[,1],Sigs[,1])
postout2 <- qrls.post(y,X,Z,alpha[2],nsamples,nburn,BetaInit[,2],GammaInit[,2],Sigs[,2])
hist(postout1$Theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
hist(postout1$Theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma[0]), main='')
hist(postout2$Theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
hist(postout2$Theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')

