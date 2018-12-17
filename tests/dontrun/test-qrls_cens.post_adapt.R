require(bayesEL)
source("hlm-functions.R")
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("gen_eps.R")
source("smoothEL.R")
source("mode-functions.R")
library(numDeriv)

# ---- quant reg: X and Z both dim 1 ----
# dimensions
n <- 200 # number of observations
p <- 1
q <- 1

# covariates
# X <- matrix(rep(1,n),n,p)
X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept
Z <- matrix(rnorm(n),n,q)

# parameters
beta0 <- 1
beta0
gamma0 <- -1
gamma0
sig20 <- 1.5
sig20

# quantile level
alpha <- 0.75

# dist is one of "norm","t","chisq","lnorm"
genout <- gen_eps(n, dist = "chisq", df = 5, tau = alpha)
eps <- genout$eps
nu0 <- genout$nu0

# resposes
# yy <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*eps)
# plot(yy,cex=0.3)
# # random censoring
# cc <- rnorm(n,mean=2,sd=1)
# deltas <- yy<=cc
# y <- yy
# sum(1-deltas)/n
# y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]
cc <- rnorm(n,mean=1.35,sd=1)
deltas <- eps<=cc
sum(1-deltas)/n
eps[!deltas] <- cc[!deltas]
y <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*eps)

# calculate the conditional posterior
numpoints <- 100
beta.seq <- seq(beta0-1,beta0+1,length.out = numpoints)
gamma.seq <- seq(gamma0-1,gamma0+1,length.out = numpoints)
sig2.seq <- seq(sig20-1,sig20+1,length.out = numpoints)
nu.seq <- seq(nu0-1,nu0+1,length.out = numpoints)
# Note: need to keep the sig2.seq range > 0 mostly
logel.seq <- matrix(rep(NA,4*numpoints),4,numpoints)
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- qrls.evalG(y,X,Z,alpha,beta.seq[ii],gamma0,sig20,nu0)
  epsilons <- evalEpsilonsLS_R(y,X,Z,beta.seq[ii],gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq[1,ii] <- logEL(omegas,epsilons,deltas)

  G <- qrls.evalG_R(y,X,Z,alpha,beta0,gamma.seq[ii],sig20,nu0)
  epsilons <- evalEpsilonsLS_R(y,X,Z,beta0,gamma.seq[ii],sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq[2,ii] <- logEL(omegas,epsilons,deltas)

  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0,sig2.seq[ii],nu0)
  epsilons <- evalEpsilonsLS_R(y,X,Z,beta0,gamma0,sig2.seq[ii])
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq[3,ii] <- logEL(omegas,epsilons,deltas)
  
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0,sig20,nu.seq[ii])
  epsilons <- evalEpsilonsLS_R(y,X,Z,beta0,gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq[4,ii] <- logEL(omegas,epsilons,deltas)
}
logelmode1 <- plotEL(beta.seq, logel.seq[1,], beta0, NA, expression(beta))
logelmode2 <- plotEL(gamma.seq, logel.seq[2,], gamma0, NA, expression(gamma))
logelmode3 <- plotEL(sig2.seq, logel.seq[3,], sig20, NA, expression(sigma^2))
logelmode4 <- plotEL(nu.seq, logel.seq[4,], nu0, NA, expression(nu))

# smoothed version
numpoints <- 100
beta.seq <- seq(beta0-1,beta0+1,length.out = numpoints)
gamma.seq <- seq(gamma0-1,gamma0+1,length.out = numpoints)
gamma.seq <- seq(-1.25,-1.05, length.out = numpoints)
sig2.seq <- seq(sig20-1,sig20+1,length.out = numpoints)
nu.seq <- seq(nu0-1,nu0+1,length.out = numpoints)
# Note: need to keep the sig2.seq range > 0 mostly
logel.seq <- matrix(rep(NA,4*numpoints),4,numpoints)
s <- 10
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  # logel.seq[1,ii] <- -qrls_cens.neglogEL.smooth_R(y,X,Z,deltas,alpha,theta = c(beta.seq[ii],gamma0,sig20,nu0),s)
  logel.seq[2,ii] <- -qrls_cens.neglogEL.smooth_R(y,X,Z,deltas,alpha,theta = c(beta0,gamma.seq[ii],sig20,nu0),s)
  # logel.seq[3,ii] <- -qrls_cens.neglogEL.smooth_R(y,X,Z,deltas,alpha,theta = c(beta0,gamma0,sig2.seq[ii],nu0),s)
  # logel.seq[4,ii] <- -qrls_cens.neglogEL.smooth_R(y,X,Z,deltas,alpha,theta = c(beta0,gamma0,sig20,nu.seq[ii]),s)

  # G <- qrls.evalG.smooth_R(y,X,Z,alpha,beta.seq[ii],gamma0,sig20,nu0,s)
  # epsilons <- evalEpsilonsLS(y,X,Z,beta.seq[ii],gamma0,sig20)
  # omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons,s)$omegas
  # logel.seq[1,ii] <- logEL.smooth_R(omegas,epsilons,deltas,s)
  # 
  # G <- qrls.evalG.smooth_R(y,X,Z,alpha,beta0,gamma.seq[ii],sig20,nu0,s)
  # epsilons <- evalEpsilonsLS(y,X,Z,beta0,gamma.seq[ii],sig20)
  # omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons,s)$omegas
  # logel.seq[2,ii] <- logEL.smooth_R(omegas,epsilons,deltas,s)
  # 
  # G <- qrls.evalG.smooth_R(y,X,Z,alpha,beta0,gamma0,sig2.seq[ii],nu0,s)
  # epsilons <- evalEpsilonsLS(y,X,Z,beta0,gamma0,sig2.seq[ii])
  # omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons,s)$omegas
  # logel.seq[3,ii] <- logEL.smooth_R(omegas,epsilons,deltas,s)
  # 
  # G <- qrls.evalG.smooth_R(y,X,Z,alpha,beta0,gamma0,sig20,nu.seq[ii],s)
  # epsilons <- evalEpsilonsLS(y,X,Z,beta0,gamma0,sig20)
  # omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons,s)$omegas
  # logel.seq[4,ii] <- logEL.smooth_R(omegas,epsilons,deltas,s)
}
logelmode1 <- plotEL(beta.seq, logel.seq[1,], beta0, NA, expression(beta))
logelmode2 <- plotEL(gamma.seq, logel.seq[2,], gamma0, NA, expression(gamma))
logelmode3 <- plotEL(sig2.seq, logel.seq[3,], sig20, NA, expression(sigma^2))
logelmode4 <- plotEL(nu.seq, logel.seq[4,], nu0, NA, expression(nu))

# accelerated smoothed version
numpoints <- 200
beta.seq <- seq(beta0-1,beta0+1,length.out = numpoints)
gamma.seq <- seq(gamma0-1,gamma0+1,length.out = numpoints)
sig2.seq <- seq(sig20-1,sig20+1,length.out = numpoints)
# sig2.seq <- seq(1.4,1.6,length.out = numpoints)
nu.seq <- seq(nu0-1,nu0+1,length.out = numpoints)
# nu.seq <- seq(0.45,0.6,length.out = numpoints)
# Note: need to keep the sig2.seq range > 0 mostly
logel.seq <- matrix(rep(NA,4*numpoints),4,numpoints)
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- qrls.evalG.smooth_R(y,X,Z,alpha,beta.seq[ii],gamma0,sig20,nu0,s=10)
  epsilons <- evalEpsilonsLS(y,X,Z,beta.seq[ii],gamma0,sig20)
  logel.seq[1,ii] <- logEL_EMAC.smooth_R(G,epsilons,deltas,s=10)

  G <- qrls.evalG.smooth_R(y,X,Z,alpha,beta0,gamma.seq[ii],sig20,nu0,s=10)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,gamma.seq[ii],sig20)
  logel.seq[2,ii] <- logEL_EMAC.smooth_R(G,epsilons,deltas,s=10)

  G <- qrls.evalG.smooth_R(y,X,Z,alpha,beta0,gamma0,sig2.seq[ii],nu0,s=10)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,gamma0,sig2.seq[ii])
  logel.seq[3,ii] <- logEL_EMAC.smooth_R(G,epsilons,deltas,s=10)
  
  G <- qrls.evalG.smooth_R(y,X,Z,alpha,beta0,gamma0,sig20,nu.seq[ii],s=10)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,gamma0,sig20)
  logel.seq[4,ii] <- logEL_EMAC.smooth_R(G,epsilons,deltas,s=10)
}
logelmode1 <- plotEL(beta.seq, logel.seq[1,], beta0, NA, expression(beta))
logelmode2 <- plotEL(gamma.seq, logel.seq[2,], gamma0, NA, expression(gamma))
logelmode3 <- plotEL(sig2.seq, logel.seq[3,], sig20, NA, expression(sigma^2))
logelmode4 <- plotEL(nu.seq, logel.seq[4,], nu0, NA, expression(nu))

# mcmc
nsamples <- 10000
nburn <- 2000
theta.hat <- hlm.fit(y = y, X = X, W = cbind(1,Z)) # solve by hlm
# quick check 
rbind(true = beta0, est = theta.hat$beta) # beta
rbind(true = c(log(sig20),gamma0), est = theta.hat$gamma) # gamma

betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma[2]/2 # Note: divided by 2 since the models differ by this factor
sig2Init <- exp(theta.hat$gamma[1])
eps_new <- c((y - X %*% theta.hat$beta)*exp(-0.5*Z %*% theta.hat$gamma[2])/sqrt(sig2Init))
nuInit <- quantile(eps_new,alpha)
mwgSd <- c(0.1,0.1,0.1,0.1)
RvDoMcmc <- c(1,1,1,1)
system.time(
  # mrout <- postCens_R(Gfun=qrls.evalG_R,nThe=4,nBet=1,nGam=1,
  #                       y=y,X=X,Z=Z,deltas=deltas,alpha=alpha,
  #                     thetaInit=c(betaInit,gammaInit,sig2Init,nuInit), 
  #                       nsamples=nsamples,nburn=nburn,
  #                       mwgSds=mwgSd,adjust=FALSE)
  # postout <- qrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sig2Init,mwgSd)
  postout <- qrls_cens.post_adapt(y,X,Z,deltas,alpha,
                                  nsamples,nburn,
                                  betaInit,gammaInit,sig2Init,nuInit,
                                  mwgSd,rvDoMcmc=RvDoMcmc)
)
theta_chain <- postout$theta_chain
theta_accept <- postout$paccept
theta_accept

# mixing of the chain
plot(theta_chain[4,],type = 'l')

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
lines(gamma.seq/2, norm_pdf(logel.marg2, gamma.seq/2),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(gamma.seq/2, norm_pdf(logel.seq[2,], gamma.seq/2),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0/2,col='red') # Note: divided by 2 since the models differ by this factor
abline(v=mean(theta_chain[2,]),col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(sigma^2), main='')
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

hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(nu), main='')
# marginal line
lines(nu.seq, norm_pdf(logel.marg4, nu.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(nu.seq, norm_pdf(logel.seq[4,], nu.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=nu0,col='red')
abline(v=mean(theta_chain[4,]),col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- X dim 2 Z dim 1 ----
n <- 200
tau <- 0.75
beta0 <- c(0.5,1)
gamma0 <- -0.5
sig20 <- 1.5

X <- cbind(rep(1,n),rnorm(n))
Z <- matrix(rnorm(1*n),n,1)
genout <- gen_eps(n, dist = "chisq", df = 5, tau = tau)
eps <- genout$eps
nu0 <- genout$nu0
theta0 <- c(beta0,gamma0,sig20,nu0)
plot(yy~X[,2], cex=0.3)
# random censoring
cc <- rnorm(n,mean=1.35,sd=1)
deltas <- eps<=cc
sum(1-deltas)/n
eps[!deltas] <- cc[!deltas]
y <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*eps)
hlmout <- hlm(y, deltas, X, cbind(1,Z))
hlmout$conv

beta.hlm <- hlmout$coef$beta
gamma.hlm <- hlmout$coef$gamma[2]*0.5
sig2.hlm <- exp(hlmout$coef$gamma[1])
nu.hlm <- quantile((y-X %*% beta.hlm)*exp(-Z %*% gamma.hlm)/sqrt(sig2.hlm),alpha)
c(beta.hlm,gamma.hlm,sig2.hlm,nu.hlm)

# for plotting conditional curves
numpoints <- 100

beta.seq1 <- seq(-1+beta0[1],1+beta0[1],length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,c(beta.seq1[ii],beta0[2]),gamma0,sig20,nu0)
  epsilons <- evalEpsilonsLS_R(y,X,Z,c(beta.seq1[ii],beta0[2]),gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq1[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode1 <- plotEL(beta.seq1, logel.seq1, beta0[1], NA, expression(beta[0]))

# smoothed version
beta.seq1 <- seq(-1+beta0[1],1+beta0[1],length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  logel.seq1[ii] <- -qrls_cens.neglogEL.smooth(y,X,Z,deltas,tau,
                                               c(beta.seq1[ii],beta0[2],gamma0,sig20,nu0))
}
logelmode1 <- plotEL(beta.seq1, logel.seq1, beta0[1], NA, expression(beta[0]))

beta.seq2 <- seq(-1+beta0[2],1+beta0[2],length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,c(beta0[1],beta.seq2[ii]),gamma0,sig20,nu0)
  epsilons <- evalEpsilonsLS_R(y,X,Z,c(beta0[1],beta.seq2[ii]),gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq2[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode2 <- plotEL(beta.seq2, logel.seq2, beta0[2], NA, expression(beta[1]))

# smoothed version
beta.seq2 <- seq(-1+beta0[2],1+beta0[2],length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  logel.seq2[ii] <- -qrls_cens.neglogEL.smooth(y,X,Z,deltas,tau,
                                               c(beta0[1],beta.seq2[ii],gamma0,sig20,nu0))
}
logelmode2 <- plotEL(beta.seq2, logel.seq2, beta0[2], NA, expression(beta[1]))

gamma.seq <- seq(-1+gamma0[1],1+gamma0[1],length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma.seq[ii],sig20,nu0)
  epsilons <- evalEpsilonsLS_R(y,X,Z,beta0,gamma.seq[ii],sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq3[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode3 <- plotEL(gamma.seq, logel.seq3, gamma0[1], NA, expression(gamma[0]))

# smoothed version
gamma.seq <- seq(-1+gamma0[1],1+gamma0[1],length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  logel.seq3[ii] <- -qrls_cens.neglogEL.smooth(y,X,Z,deltas,tau,
                                               c(beta0,gamma.seq[ii],sig20,nu0))
}
logelmode3 <- plotEL(gamma.seq, logel.seq3, gamma0[1], NA, expression(gamma[0]))

sig2.seq <- seq(-1+sig20,1+sig20,length.out = numpoints)
logel.seq4 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0,sig2.seq[ii],nu0)
  epsilons <- evalEpsilonsLS_R(y,X,Z,beta0,gamma0,sig2.seq[ii])
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq4[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode4 <- plotEL(sig2.seq, logel.seq4, sig20, NA, expression(sigma^2))

# smoothed version
sig2.seq <- seq(-1+sig20,1+sig20,length.out = numpoints)
logel.seq4 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  logel.seq4[ii] <- -qrls_cens.neglogEL.smooth(y,X,Z,deltas,tau,
                                                 c(beta0,gamma0,sig2.seq[ii],nu0))
}
logelmode4 <- plotEL(sig2.seq, logel.seq4, sig20, NA, expression(sigma^2))

nu.seq <- seq(-1+nu0,1+nu0,length.out = numpoints)
logel.seq5 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0,sig20,nu.seq[ii])
  epsilons <- evalEpsilonsLS_R(y,X,Z,beta0,gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq5[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode5 <- plotEL(nu.seq, logel.seq5, nu0, NA, expression(nu))

# smoothed version
nu.seq <- seq(-1+nu0,1+nu0,length.out = numpoints)
logel.seq5 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  logel.seq5[ii] <- -qrls_cens.neglogEL.smooth(y,X,Z,deltas,tau,
                                               c(beta0,gamma0,sig20,nu.seq[ii]))
}
logelmode5 <- plotEL(nu.seq, logel.seq5, nu0, NA, expression(nu))

betaInit <- beta.hlm
gammaInit <- gamma.hlm
sig2Init <- sig2.hlm
nuInit <- nu.hlm
# run adaptive MCMC here
nsamples <- 10000
nburn <- 5000
mwgSds <- rep(0.1,5)
postout <- qrls_cens.post_adapt(y,X,Z,deltas,alpha,nsamples,nburn,
                                beta0,gammaInit,sig20,nu0,
                                mwgSd = mwgSds, rvDoMcmc = c(0,0,1,0,0))
postout$paccept
theta_chain <- postout$theta_chain

# mixing of the chain
plot(theta_chain[3,],type='l')

# histograms and conditional loglikelihood overlays
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
     xlab = expression(gamma[1]), main='')
# conditional line
lines(gamma.seq, norm_pdf(logel.seq3, gamma.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0,col='red')
abline(v=mean(theta_chain[3,]),col='blue')

hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(nu), main='')
# conditional line
lines(sig2.seq, norm_pdf(logel.seq4, sig2.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=nu0,col='red')
abline(v=mean(theta_chain[4,]),col='blue')

hist(theta_chain[5,],breaks=50,freq=FALSE,
     xlab = expression(nu), main='')
# conditional line
lines(nu.seq, norm_pdf(logel.seq5, nu.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=nu0,col='red')
abline(v=mean(theta_chain[5,]),col='blue')

# ---- X and Z both dim 2 ----
# dimensions
n <- 100 # number of observations
p <- 2
q <- 2
X <- cbind(rep(1,n),rnorm(n))
# X <- matrix(rnorm(p*n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept
Z <- matrix(rnorm(q*n),n,q)

# parameters
# beta0 <- rnorm(p)
beta0 <- c(0.5,1)
# gamma0 <- rnorm(q)
gamma0 <- c(-1,-0.5)
# sig20 <- abs(rnorm(1,mean=1))
sig20 <- 1

# quantile level
alpha <- 0.75

# dist is one of "norm","t","chisq","lnorm"
genout <- gen_eps(n, dist = "norm", df = NULL, tau = alpha)
eps <- genout$eps
nu0 <- genout$nu0

# response
yy <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*eps)
plot(yy~X[,2], cex=0.3)
# random censoring
cc <- rnorm(n,mean=2,sd=1)
deltas <- yy<=cc
y <- yy
sum(1-deltas)/n
y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

# for plotting conditional curves
numpoints <- 100

beta.seq1 <- seq(-1+beta0[1],1+beta0[1],length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,c(beta.seq1[ii],beta0[2]),gamma0,sig20,nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,c(beta.seq1[ii],beta0[2]),gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq1[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode1 <- plotEL(beta.seq1, logel.seq1, beta0[1], NA, expression(beta[0]))

beta.seq2 <- seq(-1+beta0[2],1+beta0[2],length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,c(beta0[1],beta.seq2[ii]),gamma0,sig20,nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,c(beta0[1],beta.seq2[ii]),gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq2[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode2 <- plotEL(beta.seq2, logel.seq2, beta0[2], NA, expression(beta[1]))

gamma.seq1 <- seq(-1+gamma0[1],1+gamma0[1],length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,c(gamma.seq1[ii],gamma0[2]),sig20,nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,c(gamma.seq1[ii],gamma0[2]),sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq3[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode3 <- plotEL(gamma.seq1, logel.seq3, gamma0[1], NA, expression(gamma[0]))

gamma.seq2 <- seq(-1+gamma0[2],1+gamma0[2],length.out = numpoints)
logel.seq4 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,c(gamma0[1],gamma.seq2[ii]),sig20,nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,c(gamma.seq1[ii],gamma0[2]),sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq4[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode4 <- plotEL(gamma.seq2, logel.seq4, gamma0[2], NA, expression(gamma[1]))

sig2.seq <- seq(-1+sig20,1+sig20,length.out = numpoints)
logel.seq5 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0,sig2.seq[ii],nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,gamma0,sig2.seq[ii])
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq5[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode5 <- plotEL(sig2.seq, logel.seq5, sig20, NA, expression(sigma^2))

nu.seq <- seq(-1+nu0,1+nu0,length.out = numpoints)
logel.seq6 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0,sig20,nu.seq[ii])
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq6[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode6 <- plotEL(nu.seq, logel.seq6, nu0, NA, expression(nu))

# mcmc
nsamples <- 2000
nburn <- 1000
# solve by hlm
theta.hat <- hlm.fit(y = y, X = X, W = cbind(1,Z))
betaInit <- theta.hat$beta
gammaInit <- 0.5*theta.hat$gamma[2:(q+1)]
sig2Init <- exp(theta.hat$gamma[1])
nuInit <- quantile((y-X %*% betaInit)*exp(-Z %*% gammaInit)/sqrt(sig2Init),alpha)
mwgSd <- c(0.2,0.2,0.2,0.2,0.2,0.2)

# choose which parameters to update, and others are fixed
RvDoMcmc <- rep(1,6)
# RvDoMcmc[3] <- 1

system.time(
  # postout <- qrls.post(y,X,Z,alpha,nsamples,nburn,beta0,gamma0/2,nu0,mwgSd,RvDoMcmc)
  postout <- qrls_cens.post_adapt(y,X,Z,deltas,alpha,nsamples,nburn,
                                  betaInit,gammaInit,sig2Init,nuInit,
                                  mwgSd,RvDoMcmc)
)
theta_chain <- postout$theta_chain
theta_accept <- postout$paccept
theta_accept

# mixing of the chains
plot(theta_chain[6,],type = 'l')

# histograms and conditional loglikelihood overlays
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
# conditional line
lines(beta.seq1, norm_pdf(logel.seq1, beta.seq1),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0[1],col='red')
abline(v=mean(theta_chain[1,]),col='blue')

hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[2]), main='')
# conditional line
lines(beta.seq2, norm_pdf(logel.seq2, beta.seq2),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0[2],col='red')
abline(v=mean(theta_chain[2,]),col='blue')

hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')
# conditional line
lines(gamma.seq1, norm_pdf(logel.seq3, gamma.seq1),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0[1],col='red')
abline(v=mean(theta_chain[3,]),col='blue')

hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(gamma[2]), main='')
# conditional line
lines(gamma.seq2, norm_pdf(logel.seq4, gamma.seq2),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0[2],col='red')
abline(v=mean(theta_chain[4,]),col='blue')

hist(theta_chain[6,],breaks=50,freq=FALSE,
     xlab = expression(nu), main='')
# conditional line
lines(nu.seq, norm_pdf(logel.seq6, nu.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=nu0,col='red')
abline(v=mean(theta_chain[6,]),col='blue')

# ---- coverage prob ----
# dimensions
n <- 200 # number of observations
p <- 2
q <- 1
X <- cbind(rep(1,n),rnorm(n))
Z <- matrix(rnorm(q*n),n,q)
# X <- matrix(rnorm(p*n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept

# parameters
beta0 <- rnorm(p)
# beta0 <- c(1,0.5)
beta0 
gamma0 <- rnorm(q)
# gamma0 <- -0.5
gamma0
sig20 <- abs(rnorm(1,mean=1))
# sig20 <- 1
sig20

# quantile level
alpha <- 0.75

# dist is one of "norm","t","chisq","lnorm"
genout <- gen_eps(n, dist = "norm", df = NULL, tau = alpha)
eps <- genout$eps
nu0 <- genout$nu0

# response
yy <- c(X %*% beta0 + sqrt(sig20)*exp(0.5 * Z %*% gamma0)*eps)
plot(yy~X[,2], cex=0.3)
# random censoring
cc <- rnorm(n,mean=2,sd=1)
deltas <- yy<=cc
y <- yy
sum(1-deltas)/n
y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

# for plotting conditional curves
numpoints <- 100

beta.seq1 <- seq(-1+beta0[1],1+beta0[1],length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,c(beta.seq1[ii],beta0[2]),0.5*gamma0,sig20,nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,c(beta.seq1[ii],beta0[2]),0.5*gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq1[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode1 <- plotEL(beta.seq1, logel.seq1, beta0[1], NA, expression(beta[0]))

beta.seq2 <- seq(-1+beta0[2],1+beta0[2],length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,c(beta0[1],beta.seq2[ii]),0.5*gamma0,sig20,nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,c(beta0[1],beta.seq2[ii]),0.5*gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq2[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode2 <- plotEL(beta.seq2, logel.seq2, beta0[2], NA, expression(beta[1]))

gamma.seq <- seq(-1+0.5*gamma0[1],1+0.5*gamma0[1],length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma.seq[ii],sig20,nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,gamma.seq[ii],sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq3[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode3 <- plotEL(gamma.seq, logel.seq3, 0.5*gamma0[1], NA, expression(gamma[0]))

sig2.seq <- seq(-1+sig20,1+sig20,length.out = numpoints)
logel.seq4 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,0.5*gamma0,sig2.seq[ii],nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,0.5*gamma0,sig2.seq[ii])
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq4[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode4 <- plotEL(sig2.seq, logel.seq4, sig20, NA, expression(sigma^2))

nu.seq <- seq(-1+nu0,1+nu0,length.out = numpoints)
logel.seq5 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,0.5*gamma0,sig20,nu.seq[ii])
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,0.5*gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq5[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode5 <- plotEL(nu.seq, logel.seq5, nu0, NA, expression(nu))

# mcmc
nsamples <- 5000
nburn <- 2000
# solve by hlm
theta.hat <- hlm.fit(y = y, X = X, W = cbind(1,Z))
betaInit <- theta.hat$beta
gammaInit <- 0.5*theta.hat$gamma[2:(q+1)]
sig2Init <- exp(theta.hat$gamma[1])
nuInit <- quantile((y-X %*% betaInit)*exp(-Z %*% gammaInit)/sqrt(sig2Init),alpha)
mwgSd <- c(0.2,0.2,0.2,0.2,0.2)

# choose which parameters to update, and others are fixed
RvDoMcmc <- rep(1,5)

system.time(
  postout <- qrls_cens.post_adapt(y,X,Z,deltas,alpha,nsamples,nburn,
                                  betaInit,gammaInit,sig2Init,nuInit,
                                  mwgSd,RvDoMcmc)
)
theta_chain <- postout$theta_chain
theta_accept <- postout$paccept
theta_accept

# mixing of the chains
plot(theta_chain[5,],type = 'l')
# plot(theta_chain[3,],theta_chain[5,],cex=.3)
pairs(t(theta_chain),cex=.3)

# histograms and conditional loglikelihood overlays
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
# conditional line
lines(beta.seq1, norm_pdf(logel.seq1, beta.seq1),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0[1],col='red')
abline(v=mean(theta_chain[1,]),col='blue')

hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[2]), main='')
# conditional line
lines(beta.seq2, norm_pdf(logel.seq2, beta.seq2),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0[2],col='red')
abline(v=mean(theta_chain[2,]),col='blue')

hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(gamma), main='')
# conditional line
lines(gamma.seq, norm_pdf(logel.seq3, gamma.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0[1]/2,col='red')
abline(v=mean(theta_chain[3,]),col='blue')

hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(sigma^2), main='')
# conditional line
lines(sig2.seq, norm_pdf(logel.seq4, sig2.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=sig20,col='red')
abline(v=mean(theta_chain[4,]),col='blue')

hist(theta_chain[5,],breaks=50,freq=FALSE,
     xlab = expression(sigma^2), main='')
# conditional line
lines(nu.seq, norm_pdf(logel.seq5, nu.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=nu0,col='red')
abline(v=mean(theta_chain[5,]),col='blue')

# repeat the experiment and calculate coverage probilities
con <- file("covprob.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

ns <- c(100,200,300)
alpha <- 0.75
beta0 <- c(1,0.5)
gamma0 <- -0.5
sig20 <- 1
nu0 <- gen_eps(100,tau=alpha)$nu0
theta0 <- c(beta0,gamma0,sig20,nu0)

nsamples <- 10000
nburn <- 2000
mwgSds <- rep(0.1,5)
reptimes <- 500
coverProbs <- matrix(rep(0,5*length(ns)),5,length(ns))
lengthCI <- matrix(rep(0,5*length(ns)),5,length(ns))
# system.time(
for (ii in 1:length(ns)) {
  system.time({
    n <- ns[ii]
    message("---- n = ", n, " ----")
    covers <- rep(0,5)
    lengths <- rep(0,5)
    # covers <- c(101, 94, 90, 96, 104)
    # lengths <- c(23.96676, 24.00932, 39.95447, 35.39038, 23.68906)
    reptimes <- 500
    for (jj in 1:reptimes) {
      message("** jj = ", jj, " **")
      # cat ("Press [enter] to continue\n")
      # line <- readline()
      if (jj %% 50 == 0) message("jj = ", jj)
      X <- cbind(rep(1,n),rnorm(n))
      Z <- matrix(rnorm(1*n),n,1)
      eps <- gen_eps(n, dist = "norm", df = NULL, tau = alpha)$eps
      yy <- c(X %*% beta0 + sqrt(sig20)*exp(0.5 * Z %*% gamma0)*eps)
      # plot(yy~X[,2], cex=0.3)
      # random censoring
      cc <- rnorm(n,mean=2.5,sd=1)
      deltas <- yy<=cc
      y <- yy
      message("censored pct = ", sum(1-deltas)/n)
      if (sum(1-deltas)/n > 0.2 || sum(1-deltas)/n < 0.1) {
        reptimes <- reptimes+1
        next
      }
      y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]
      hlmout <- hlm(y, deltas, X, cbind(1,Z))
      message("hlm conv = ", hlmout$conv)
      if (!hlmout$conv) {
        reptimes <- reptimes+1
        next
      }
      betaInit <- hlmout$coef$beta
      gammaInit <- 0.5*hlmout$coef$gamma[2]
      sig2Init <- exp(hlmout$coef$gamma[1])
      nuInit <- quantile((y-X %*% betaInit)*exp(-Z %*% gammaInit)/sqrt(sig2Init),alpha)
      # run adaptive MCMC here
      postout <- qrls_cens.post_adapt(y,X,Z,deltas,alpha,nsamples,nburn,
                                      betaInit,gammaInit,sig2Init,nuInit,
                                      mwgSd = mwgSds)
      theta_chain <- postout$theta_chain
      if (all(theta_chain == 0.0)) {
        message("MCMC not converged.")
        reptimes <- reptimes+1
        next
      }
      cat("paccept = ", postout$paccept, "\n")
      if (any(postout$paccept < 0.4) || any(postout$paccept > 0.5) ) {
        reptimes <- reptimes+1
        next
      }
      qts <- rep(NA,5)
      for (kk in 1:5) {
        qts <- quantile(theta_chain[kk,],c(0.025,0.975))
        # gamma range need to time 2 since the estimate for gamma/2
        if (kk == 3) qts <- 2*qts
        if (theta0[kk] >= qts[1] && theta0[kk] <= qts[2]) covers[kk] <- covers[kk] + 1
        lengths[kk] <- lengths[kk] + (qts[2]-qts[1])
      }
      message("reptimes = ", reptimes)
      cat("covers = ", covers, "\n")
      cat("lengths = ", lengths, "\n")
    }
    for (ll in 1:5) {
      coverProbs[ll,ii] <- covers[ll]/500
      lengthCI[ll,ii] <- lengths[ll]/500
    }
  })
}

# Restore output to console
sink() 
sink(type="message")

coverProbs
lengthCI
