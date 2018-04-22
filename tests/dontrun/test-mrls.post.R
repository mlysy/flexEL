require(bayesEL)
require(optimCheck)
source("hlm-functions.R")
source("../testthat/el-utils.R")
#---- mean reg: X and Z both dim 1----
# dimensions
n <- 100 # number of observations
p <- 1
q <- 1
# X <- matrix(rep(1,n),n,p)
X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept
Z <- matrix(rnorm(n),n,q)

beta0 <- rnorm(p)
gamma0 <- rnorm(q)
theta0 <- c(beta0,gamma0)
y <- c(X %*% beta0 + exp(Z %*% gamma0)*rnorm(n))
# plot(X,y,cex=0.3)

numpoints <- 100
beta.seq <- seq(-.5+beta0, .5+beta0, length.out = numpoints)
gamma.seq <- seq(-.5+gamma0, .5+gamma0, length.out = numpoints)

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
theta.hat <- hlm.fit(y = y, X = X, W = Z) # solve by hlm:
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
sigs <- c(0.15,0.2)

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
abline(v=beta0,col='red')
abline(v=mean(theta_chain[1,]),col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
lines(gamma.seq, norm_pdf(logel.marg[,2], beta.seq),
      cex=0.1, col = 'red', type='l')
abline(v=gamma0,col='red')
abline(v=mean(theta_chain[2,]),col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

#---- mean reg: X and Z both dim 1----
# dimensions
n <- 100 # number of observations
p <- 2
q <- 2
X <- cbind(rep(1,n),rnorm(n))
# X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept
Z <- matrix(rnorm(q*n),n,q)

beta0 <- rnorm(p)
gamma0 <- rnorm(q)
y <- c(X %*% beta0 + exp(Z %*% gamma0)*rnorm(n))

numpoints <- 100
# beta1.seq <- seq(-.5+beta0[1], .5+beta0[1], length.out = numpoints)
# beta2.seq <- seq(-.5+beta0[2], .5+beta0[2], length.out = numpoints)
# gamma1.seq <- seq(-.5+gamma0[1], .5+gamma0[1], length.out = numpoints)
# gamma2.seq <- seq(-.5+gamma0[2], .5+gamma0[2], length.out = numpoints)

# mcmc
nsamples <- 20000
nburn <- 3000
# solve by hlm
theta.hat <- hlm.fit(y = y, X = X, W = Z)
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma
sigs <- c(0.15,0.2,0.23,0.23)

system.time(
  postout <- mrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sigs)
)
theta_chain <- postout$Theta_chain
theta_accept <- postout$paccept
theta_accept
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
abline(v=beta0[1],col='red')
abline(v=mean(theta_chain[1,]),col='blue')
hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
abline(v=beta0[2],col='red')
abline(v=mean(theta_chain[2,]),col='blue')
hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(gamma[0]), main='')
abline(v=gamma0[1],col='red')
abline(v=mean(theta_chain[3,]),col='blue')
hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')
abline(v=gamma0[2],col='red')
abline(v=mean(theta_chain[4,]),col='blue')
