require(bayesEL)
source("hlm-functions.R")
source("../testthat/el-utils.R")

#---- mean reg: X and Z both dim 1----
# dimensions
n <- 100 # number of observations
p <- 1
q <- 1

# covariates
X <- matrix(rep(1,n),n,p)
# X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept
Z <- matrix(rnorm(n),n,q)

# parameters
beta0 <- rnorm(p)
gamma0 <- rnorm(q)
sig20 <- abs(rnorm(1)) # NEW: scale param

# normal(0,1) error
eps <- rnorm(n)

# t error with mean 0 var 1
df <- 5
v <- df/(df-2)
eps <- rt(n, df=df)/sqrt(v)

# chi-sqr with mean 0 var 1
df <- 3
v <- 2*df
m <- df
eps <- (rchisq(n, df=df)-m)/sqrt(v)

# log-normal with mean 0 var 1
mn <- 0
sn <- 1
m <- exp(mn+sn^2/2)
v <- (exp(sn^2)-1)*exp(2*mn+sn^2)
eps <- (rlnorm(n,mn,sn)-m)/sqrt(v)

# response
y <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*eps)
plot(y,cex=0.3)
# plot(X,y,cex=0.3)

# # calculate marginal posterior
# numpoints <- 50
# beta.seq <- seq(-.5+beta0, .5+beta0, length.out = numpoints)
# gamma.seq <- seq(-.5+gamma0, .5+gamma0, length.out = numpoints)
# 
# Theta.seq <- as.matrix(expand.grid(beta.seq, gamma.seq))
# logel.mat <- apply(Theta.seq, 1, function(bb) {
#   G <- mrls.evalG(y,X,Z,bb[1],bb[2])
#   logEL(G)
# })
# logel.mat <- matrix(logel.mat, numpoints, numpoints)
# el.mat <- exp(logel.mat - max(logel.mat))
# logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

# calculate the marginal posterior distribution
# Note: this might take some time to calculate
numpoints <- 50
beta.seq <- seq(-1+beta0, 1+beta0, length.out = numpoints)
gamma.seq <- seq(-1+gamma0, 1+gamma0, length.out = numpoints)
sig2.seq <- seq(-1+sig20, 1+sig20, length.out = numpoints)

Theta.seq <- as.matrix(expand.grid(beta.seq, gamma.seq, sig2.seq))
logel.mat <- apply(Theta.seq, 1, function(bb) {
  G <- mrls.evalG(y,X,Z,bb[1],bb[2],bb[3])
  logEL(G)
})
logel.mat <- array(logel.mat, c(numpoints, numpoints,numpoints))
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg1 <- log(apply(el.mat, MARGIN=1, sum))
logel.marg2 <- log(apply(el.mat, MARGIN=2, sum))
logel.marg3 <- log(apply(el.mat, MARGIN=3, sum))

# mcmc
nsamples <- 20000
nburn <- 3000
theta.hat <- hlm.fit(y = y, X = X, W = cbind(1,Z)) # solve by hlm:
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma[2]
sig2Init <- exp(2*theta.hat$gamma[1])
mwgSd <- c(0.1,0.25,0.15)

system.time(
  # postout <- mrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,mwgSd)
  postout <- mrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sig2Init,mwgSd)
)
theta_chain <- postout$Theta_chain
theta_accept <- postout$paccept
theta_accept
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta), main='')
lines(beta.seq, norm_pdf(logel.marg1, beta.seq),
      cex=0.1, col = 'red', type='l')
# lines(beta.seq, norm_pdf(logel.marg[,1], beta.seq),
#       cex=0.1, col = 'red', type='l')
abline(v=beta0,col='red')
abline(v=mean(theta_chain[1,]),col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma), main='')
lines(gamma.seq, norm_pdf(logel.marg2, gamma.seq),
      cex=0.1, col = 'red', type='l')
# lines(gamma.seq, norm_pdf(logel.marg[,2], beta.seq),
#       cex=0.1, col = 'red', type='l')
abline(v=gamma0,col='red')
abline(v=mean(theta_chain[2,]),col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(sigma), main='')
lines(sig2.seq, norm_pdf(logel.marg3, sig2.seq),
      cex=0.1, col = 'red', type='l')
abline(v=sig20,col='red')
abline(v=mean(theta_chain[3,]),col='blue')
legend('topright',legend=c(expression('grid plot & true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

#---- mean reg: X and Z both dim 2 ----
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
sig20 <- abs(rnorm(1))
y <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*rnorm(n))

# solve by hlm
theta.hat <- hlm.fit(y = y, X = X, W = cbind(1,Z))
# mcmc
nsamples <- 20000
nburn <- 3000
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma[2:(q+1)]
sig2Init <- exp(2*theta.hat$gamma[1])
mwgSd <- c(0.2,0.2,0.23,0.23,0.35)
RvDoMcmc <- matrix(c(1,1,1,1,1),nrow=5,ncol=1)
system.time(
  postout <- mrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sig2Init,mwgSd)
)
theta_chain <- postout$Theta_chain
theta_accept <- postout$paccept
theta_accept
# mixing of the chains
plot(theta_chain[5,],type = 'l')
# histograms
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
hist(theta_chain[5,],breaks=50,freq=FALSE,
     xlab = expression(sigma^2), main='')
abline(v=sig20,col='red')
abline(v=mean(theta_chain[5,]),col='blue')

