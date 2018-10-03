require(bayesEL)
source("hlm-functions.R")
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("gen_eps.R")

#---- mean reg: X and Z both dim 1----
# dimensions
n <- 300 # number of observations
p <- 1
q <- 1

# covariates
X <- matrix(rep(1,n),n,p)
# X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept
Z <- matrix(rnorm(n),n,q)
# Z <- matrix(runif(n,min=1,max=5),n,q)
# Z <- matrix(rchisq(n,df=5),n,q)

# parameters
# beta0 <- rnorm(p)
# gamma0 <- rnorm(q)
# sig20 <- abs(rnorm(1)) # NEW: scale param
beta0 <- 1
gamma0 <- -0.5
sig20 <- 1

# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "t", df = 5)

# response
# Note: the 0.5 is consistent with hlm.fit but not the LS model in MCMC
y <- c(X %*% beta0 + sqrt(sig20)*exp(0.5 * Z %*% gamma0)*eps)
plot(y,cex=0.3)
# plot(X,y,cex=0.3)

# # calculate marginal posterior (old: no location param)
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
numpoints <- 30
beta.seq <- seq(-1+beta0, 1+beta0, length.out = numpoints)
gamma.seq <- seq(-1+gamma0, 1+gamma0, length.out = numpoints)
sig2.seq <- seq(-1+sig20, 1+sig20, length.out = numpoints)

Theta.seq <- as.matrix(expand.grid(beta.seq, gamma.seq, sig2.seq))
system.time(
  logel.mat <- apply(Theta.seq, 1, function(bb) {
    G <- mrls.evalG(y,X,Z,bb[1],bb[2],bb[3])
    omegas <- omega.hat(G)
    logEL(omegas)
  })
)
logel.mat <- array(logel.mat, c(numpoints, numpoints,numpoints))
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg1 <- log(apply(el.mat, MARGIN=1, sum))
logel.marg2 <- log(apply(el.mat, MARGIN=2, sum))
logel.marg3 <- log(apply(el.mat, MARGIN=3, sum))

# calculate the conditional posterior
numpoints <- 100
beta.seq <- seq(beta0-1,beta0+1,length.out = numpoints)
gamma.seq <- seq(gamma0/2-1,gamma0/2+1,length.out = numpoints)
sig2.seq <- seq(sig20-1,sig20+1,length.out = numpoints)
# Note: need to keep the sig2.seq range > 0 mostly
logel.seq <- matrix(rep(NA,3*numpoints),3,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,beta.seq[ii],gamma0/2,sig20)
  omegas <- omega.hat(G)
  logel.seq[1,ii] <- logEL(omegas)
  G <- mrls.evalG(y,X,Z,beta0,gamma.seq[ii],sig20)
  omegas <- omega.hat(G)
  logel.seq[2,ii] <- logEL(omegas)
  G <- mrls.evalG(y,X,Z,beta0,gamma0/2,sig2.seq[ii])
  omegas <- omega.hat(G)
  logel.seq[3,ii] <- logEL(omegas)
}
logelmode1 <- plotEL(beta.seq, logel.seq[1,], beta0, NA, expression(beta))
logelmode2 <- plotEL(gamma.seq, logel.seq[2,], gamma0/2, NA, expression(gamma))
logelmode3 <- plotEL(sig2.seq, logel.seq[3,], sig20, NA, expression(sigma^2))
cbind(logelmode1,logelmode2,logelmode3)

# mcmc
nsamples <- 20000
nburn <- 3000
theta.hat <- hlm.fit(y = y, X = X, W = cbind(1,Z)) # solve by hlm
# quick check 
rbind(true = beta0, est = theta.hat$beta) # beta
rbind(true = c(log(sig20),gamma0), est = theta.hat$gamma) # gamma

betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma[2]/2 # Note: divided by 2 since the models differ by this factor
sig2Init <- exp(theta.hat$gamma[1])
mwgSd <- c(0.1,0.15,0.2)
# RvDoMcmc <- c(0,0,1)
RvDoMcmc <- c(1,0,0)
system.time(
  # postout <- mrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sig2Init,mwgSd)
  # postout <- mrls.post(y,X,Z,nsamples,nburn,beta0,gamma0/2,sig2Init,mwgSd,RvDoMcmc)
  postout <- mrls.post(y,X,Z,nsamples,nburn,beta0,gamma0,sig20,mwgSd,RvDoMcmc)
)
theta_chain <- postout$theta_chain
theta_accept <- postout$paccept
theta_accept

plot(theta_chain[1,],type='l')

# overlay marginal / conditional gird plot to histogram
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta), main='')
# marginal line
lines(beta.seq, norm_pdf(logel.marg1, beta.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(beta.seq, norm_pdf(logel.seq[1,], beta.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0,col='red')
abline(v=mean(theta_chain[1,]),col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma), main='')
# marginal line
lines(gamma.seq, norm_pdf(logel.marg2, gamma.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(gamma.seq, norm_pdf(logel.seq[2,], gamma.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0/2,col='red') # Note: divided by 2 since the models differ by this factor
abline(v=mean(theta_chain[2,]),col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(sigma), main='')
# marginal line
lines(sig2.seq, norm_pdf(logel.marg3, sig2.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(sig2.seq, norm_pdf(logel.seq[3,], sig2.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=sig20,col='red')
abline(v=mean(theta_chain[3,]),col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

#---- mean reg: X and Z both dim 2 ----
# dimensions
n <- 500 # number of observations
p <- 2
q <- 2

# covariates
# X <- cbind(rep(1,n),rnorm(n))
X <- matrix(rnorm(p*n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept
Z <- matrix(rnorm(q*n),n,q)

# TODO: random params seem to work pretty badly
# beta0 <- rnorm(p)
# beta0
# gamma0 <- rnorm(q)
# gamma0
# sig20 <- abs(rnorm(1))
# sig20
beta0 <- c(0.5,1.5)
gamma0 <- c(-0.5,0.5)
sig20 <- 0.5

# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "t", df = 5)

# response
y <- c(X %*% beta0 + sqrt(sig20)*exp(0.5 * Z %*% gamma0)*eps)
plot(X[,2],y,cex=0.3)

# for plotting conditional curves
numpoints <- 100

beta.seq1 <- seq(-.5+beta0[1],.5+beta0[1],length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,c(beta.seq1[ii],beta0[2]),gamma0/2,sig20)
  omegas <- omega.hat(G)
  logel.seq1[ii] <- logEL(omegas)
}
logelmode1 <- plotEL(beta.seq1, logel.seq1, beta0[1], NA, expression(beta[0]))

beta.seq2 <- seq(-.5+beta0[2],.5+beta0[2],length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,c(beta0[1],beta.seq2[ii]),gamma0/2,sig20)
  omegas <- omega.hat(G)
  logel.seq2[ii] <- logEL(omegas)
}
logelmode2 <- plotEL(beta.seq2, logel.seq2, beta0[2], NA, expression(beta[1]))

gamma.seq1 <- seq(-.5+gamma0[1]/2,.5+gamma0[1]/2,length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,beta0,c(gamma.seq1[ii],gamma0[2]/2),sig20)
  omegas <- omega.hat(G)
  logel.seq3[ii] <- logEL(omegas)
}
logelmode3 <- plotEL(gamma.seq1, logel.seq3, gamma0[1]/2, NA, expression(gamma[1]))

gamma.seq2 <- seq(-.5+gamma0[2]/2,.5+gamma0[2]/2,length.out = numpoints)
logel.seq4 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,beta0,c(gamma0[1]/2,gamma.seq2[ii]),sig20)
  omegas <- omega.hat(G)
  logel.seq4[ii] <- logEL(omegas)
}
logelmode4 <- plotEL(gamma.seq2, logel.seq4, gamma0[2]/2, NA, expression(gamma[2]))

sig2.seq <- seq(-.5+sig20,.5+sig20,length.out = numpoints)
logel.seq5 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,beta0,gamma0/2,sig2.seq[ii])
  omegas <- omega.hat(G)
  logel.seq5[ii] <- logEL(omegas)
}
logelmode5 <- plotEL(sig2.seq, logel.seq5, sig20, NA, expression(sigma^2))

# solve by hlm
theta.hat <- hlm.fit(y = y, X = X, W = cbind(1,Z))
# quick check 
rbind(true = beta0, est = theta.hat$beta) # beta
rbind(true = c(log(sig20),gamma0), est = theta.hat$gamma) # gamma

# mcmc
nsamples <- 20000
nburn <- 3000
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma[2:(q+1)]/2 # Note: divided by 2 since the models differ by this factor
sig2Init <- exp(theta.hat$gamma[1])
mwgSd <- c(0.1,0.1,0.1,0.1,0.1)
RvDoMcmc <- matrix(c(0,0,0,0,1),nrow=5,ncol=1)

system.time(
  # postout <- mrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sig2Init,mwgSd,RvDoMcmc)
  postout <- mrls.post(y,X,Z,nsamples,nburn,beta0,gamma0/2,sig2Init,mwgSd,RvDoMcmc)
)
theta_chain <- postout$Theta_chain
theta_accept <- postout$paccept
theta_accept

# mixing of the chains
plot(theta_chain[5,],type = 'l')

# histograms
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]), main='')
# conditional line
lines(beta.seq1, norm_pdf(logel.seq1, beta.seq1),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0[1],col='red')
abline(v=mean(theta_chain[1,]),col='blue')

hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
# conditional line
lines(beta.seq2, norm_pdf(logel.seq2, beta.seq2),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0[2],col='red')
abline(v=mean(theta_chain[2,]),col='blue')

hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(gamma[0]), main='')
# conditional line
lines(gamma.seq1, norm_pdf(logel.seq3, gamma.seq1),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0[1]/2,col='red')
abline(v=mean(theta_chain[3,]),col='blue')

hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')
# conditional line
lines(gamma.seq2, norm_pdf(logel.seq4, gamma.seq2),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0[2]/2,col='red')
abline(v=mean(theta_chain[4,]),col='blue')

hist(theta_chain[5,],breaks=50,freq=FALSE,
     xlab = expression(sigma^2), main='')
# conditional line
lines(sig2.seq, norm_pdf(logel.seq5, sig2.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=sig20,col='red')
abline(v=mean(theta_chain[5,]),col='blue')

