require(bayesEL)
source("hlm-functions.R")
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("gen_eps.R")

##### 1 quantile level #####

# ---- X and Z both dim 1 (works fine for all eps here) ----
# dimensions
n <- 10 # number of observations
p <- 1
q <- 1

# covariates
# X <- matrix(rep(1,n),n,p)
X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain intercept
Z <- matrix(rnorm(n),n,q)
# Z <- X

# parameters
# beta0 <- rnorm(p)
# beta0
# gamma0 <- rnorm(q)
# gamma0
# sig20 <- abs(rnorm(1))
beta0 <- 0.5
gamma0 <- -0.5
sig20 <- 1
alpha <- 0.75

# dist is one of "norm","t","chisq","lnorm"
genout <- gen_eps(n, dist = "t", df = 5, tau = alpha)
eps <- genout$eps
nu0 <- genout$nu0

# response
y <- c(X %*% beta0 + sqrt(sig20)*exp(0.5 * Z %*% gamma0)*eps)
plot(y~X, cex=0.3)

# calculate the marginal posterior distributions
numpoints <- 50
beta.seq <- seq(-1+beta0, 1+beta0, length.out = numpoints)
gamma.seq <- seq(-1+gamma0, 1+gamma0, length.out = numpoints)
sig2.seq <- seq(sig20-1,sig20+1,length.out = numpoints)
nu.seq <- seq(-1+nu0, 1+nu0, length.out = numpoints)

# Note: this may take a long while to calculate 
Theta.seq <- as.matrix(expand.grid(beta.seq, gamma.seq, sig2.seq, nu.seq))
logel.mat <- apply(Theta.seq, 1, function(bb) {
  G <- qrls.evalG(y,X,Z,alpha,bb[1],bb[2],bb[3],bb[4])
  omegas <- omega.hat(G)
  logEL(omegas)
})
logel.mat <- array(logel.mat, c(numpoints, numpoints,numpoints))
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg1 <- log(apply(el.mat, MARGIN=1, sum))
logel.marg2 <- log(apply(el.mat, MARGIN=2, sum))
logel.marg3 <- log(apply(el.mat, MARGIN=3, sum))

# calculate the conditional posterior distributions
numpoints <- 100
beta.seq <- seq(-1+beta0,1+beta0,length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta.seq[ii],gamma0/2,sig20,nu0)
  omegas <- omega.hat(G)
  logel.seq1[ii] <- logEL(omegas)
}
logelmode1 <- plotEL(beta.seq, logel.seq1, beta0, NA, expression(beta[0]))

gamma.seq <- seq(-1+gamma0/2,1+gamma0/2,length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma.seq[ii],sig20,nu0)
  omegas <- omega.hat(G)
  logel.seq2[ii] <- logEL(omegas)
}
logelmode2 <- plotEL(gamma.seq, logel.seq2, gamma0/2, NA, expression(gamma[0]))

sig2.seq <- seq(-.5+sig20,.5+sig20,length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0/2,sig2.seq[ii],nu0)
  omegas <- omega.hat(G)
  logel.seq3[ii] <- logEL(omegas)
}
logelmode3 <- plotEL(sig2.seq, logel.seq3, sig20, NA, expression(sigma^2))

nu.seq <- seq(-1+nu0,1+nu0,length.out = numpoints)
logel.seq4 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0/2,sig20,nu.seq[ii])
  omegas <- omega.hat(G)
  logel.seq4[ii] <- logEL(omegas)
}
logelmode4 <- plotEL(nu.seq, logel.seq4, nu0, NA, expression(nu))

# mcmc
nsamples <- 10000
nburn <- 3000
theta.hat <- hlm.fit(y = y, X = X, W = cbind(1,Z)) # solve by hlm
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma[2]/2
sig2Init <- exp(theta.hat$gamma[1])
eps_new <- c((y - X %*% betaInit)*exp(-Z %*% gammaInit)/sqrt(sig2Init))
nuInit <- quantile(eps_new,alpha)
sigs <- c(0.2,0.2,0.1,0.2)
RvDoMcmc <- c(0,0,1,0)
system.time(
  # postout <- qrls.post(y,X,Z,alpha,nsamples,nburn,beta0,gamma0/2,sig20,nu0,sigs,RvDoMcmc)
  postout <- qrls.post_adapt(y,X,Z,alpha,nsamples,nburn,beta0,gamma0/2,sig2Init,nu0, sigs,RvDoMcmc)
)
theta_chain <- postout$theta_chain
theta_accept <- postout$paccept
theta_accept

# mixing of the chain
plot(theta_chain[3,],type = 'l')

# histograms
hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta), main='')
# marginal line
lines(beta.seq, norm_pdf(logel.marg1, beta.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(beta.seq, norm_pdf(logel.seq1, beta.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0, col='red')
abline(v=mean(theta_chain[1,]), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(gamma), main='')
# marginal line
lines(gamma.seq, norm_pdf(logel.marg2, gamma.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(gamma.seq, norm_pdf(logel.seq2, gamma.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0/2, col='red')
abline(v=mean(theta_chain[2,]), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(nu), main='')
# marginal line
lines(sig2.seq, norm_pdf(logel.marg3, sig2.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(sig2.seq, norm_pdf(logel.seq3, sig2.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=sig20, col='red')
abline(v=mean(theta_chain[3,]), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(nu), main='')
# marginal line
lines(nu.seq, norm_pdf(logel.marg4, nu.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(nu.seq, norm_pdf(logel.seq4, nu.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=nu0, col='red')
abline(v=mean(theta_chain[4,]), col='blue')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- X and Z both dim 2 ----
# dimensions
n <- 300 # number of observations
p <- 2
q <- 2
X <- cbind(rep(1,n),rnorm(n))
# X <- matrix(rnorm(p*n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept
Z <- matrix(rnorm(q*n),n,q)

# parameters
beta0 <- rnorm(p)
beta0
gamma0 <- rnorm(q)
gamma0
sig20 <- abs(rnorm(1))
alpha <- 0.75

# dist is one of "norm","t","chisq","lnorm"
genout <- gen_eps(n, dist = "t", df = 5, tau = alpha)
eps <- genout$eps
nu0 <- genout$nu0

# response
y <- c(X %*% beta0 + sqrt(sig20)*exp(0.5 * Z %*% gamma0)*eps)
# plot(y~X[,2], cex=0.3)

# for plotting conditional curves
numpoints <- 100

beta.seq1 <- seq(-1+beta0[1],1+beta0[1],length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  # G <- qrls.evalG(y,X,Z,alpha,c(beta.seq1[ii],beta0[2]),gamma0/2,sig20,nu0)
  G <- qrls.evalG.smooth(y,X,Z,alpha,c(beta.seq1[ii],beta0[2]),gamma0/2,sig20,nu0,sp=10)
  omegas <- omega.hat(G)
  logel.seq1[ii] <- logEL(omegas)
}
logelmode1 <- plotEL(beta.seq1, logel.seq1, beta0[1], NA, expression(beta[0]))

beta.seq2 <- seq(-.1+beta0[2],.1+beta0[2],length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  # G <- qrls.evalG(y,X,Z,alpha,c(beta0[1],beta.seq2[ii]),gamma0/2,sig20,nu0)
  G <- qrls.evalG.smooth(y,X,Z,alpha,c(beta0[1],beta.seq2[ii]),gamma0/2,sig20,nu0,sp=10)
  omegas <- omega.hat(G)
  logel.seq2[ii] <- logEL(omegas)
}
logelmode2 <- plotEL(beta.seq2, logel.seq2, beta0[2], NA, expression(beta[1]))

gamma.seq1 <- seq(-1+gamma0[1]/2,1+gamma0[1]/2,length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  # G <- qrls.evalG(y,X,Z,alpha,beta0,c(gamma.seq1[ii],gamma0[2]/2),sig20,nu0)
  G <- qrls.evalG.smooth(y,X,Z,alpha,beta0,c(gamma.seq1[ii],gamma0[2]/2),sig20,nu0,sp=10)
  omegas <- omega.hat(G)
  logel.seq3[ii] <- logEL(omegas)
}
logelmode3 <- plotEL(gamma.seq1, logel.seq3, gamma0[1]/2, NA, expression(gamma[0]))

gamma.seq2 <- seq(-1+gamma0[2]/2,1+gamma0[2]/2,length.out = numpoints)
logel.seq4 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  # G <- qrls.evalG(y,X,Z,alpha,beta0,c(gamma0[1]/2,gamma.seq2[ii]),sig20,nu0)
  G <- qrls.evalG.smooth(y,X,Z,alpha,beta0,c(gamma0[1]/2,gamma.seq2[ii]),sig20,nu0,sp=10)
  omegas <- omega.hat(G)
  logel.seq4[ii] <- logEL(omegas)
}
logelmode4 <- plotEL(gamma.seq2, logel.seq4, gamma0[2]/2, NA, expression(gamma[1]))

sig2.seq <- seq(-.75+sig20,1+sig20,length.out = numpoints)
logel.seq5 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0/2,sig2.seq[ii],nu0)
  G <- qrls.evalG.smooth(y,X,Z,alpha,beta0,gamma0/2,sig2.seq[ii],nu0,sp=10)
  omegas <- omega.hat(G)
  logel.seq5[ii] <- logEL(omegas)
}
logelmode5 <- plotEL(sig2.seq, logel.seq5, sig20, NA, expression(nu))

nu.seq <- seq(-1+nu0,1+nu0,length.out = numpoints)
logel.seq6 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  # G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0/2,sig20,nu.seq[ii])
  G <- qrls.evalG.smooth(y,X,Z,alpha,beta0,gamma0/2,sig20,nu.seq[ii],sp=10)
  omegas <- omega.hat(G)
  logel.seq6[ii] <- logEL(omegas)
}
logelmode6 <- plotEL(nu.seq, logel.seq6, nu0, NA, expression(nu))

# mcmc
nsamples <- 20000
nburn <- 3000
# solve by hlm
theta.hat <- hlm.fit(y = y, X = X, W = cbind(1,Z))
betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma[2:(q+1)]/2
sig2Init <- exp(theta.hat$gamma[1])
nuInit <- quantile((y-X %*% betaInit)*exp(-0.5*Z %*% gammaInit)/sqrt(sig2Init),alpha)
mwgSd <- c(0.15,0.2,0.2,0.2,0.2,0.25)

# choose which parameters to update, and others are fixed
RvDoMcmc <- rep(0,6)
RvDoMcmc[5] <- 1

system.time(
  # postout <- qrls.post(y,X,Z,alpha,nsamples,nburn,beta0,gamma0/2,nu0,mwgSd,RvDoMcmc)
  postout <- qrls.post_adapt(y,X,Z,alpha,nsamples,nburn,
                             beta0,gamma0/2,sig2Init,nu0,mwgSd,RvDoMcmc)
)
theta_chain <- postout$theta_chain
theta_accept <- postout$paccept
theta_accept

# mixing of the chains
plot(theta_chain[5,],type = 'l')

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
abline(v=gamma0[1]/2,col='red')
abline(v=mean(theta_chain[3,]),col='blue')

hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(gamma[2]), main='')
# conditional line
lines(gamma.seq2, norm_pdf(logel.seq4, gamma.seq2),
      cex=0.1, col = 'blue', type='l')
abline(v=gamma0[2],col='red')
abline(v=mean(theta_chain[4,]),col='blue')

hist(theta_chain[5,],breaks=50,freq=FALSE,
     xlab = expression(sigma^2), main='')
# conditional line
lines(sig2.seq, norm_pdf(logel.seq5, sig2.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=sig20,col='red')
abline(v=mean(theta_chain[5,]),col='blue')

hist(theta_chain[6,],breaks=50,freq=FALSE,
     xlab = expression(nu), main='')
# conditional line
lines(nu.seq, norm_pdf(logel.seq6, nu.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=nu0,col='red')
abline(v=mean(theta_chain[5,]),col='blue')

