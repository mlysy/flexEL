require(bayesEL)
source("hlm-functions.R")
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("gen_eps.R")

# ---- quant reg: X and Z both dim 1 ----
# dimensions
n <- 300 # number of observations
p <- 1
q <- 1

# covariates
X <- matrix(rep(1,n),n,p)
# X <- matrix(rnorm(n),n,p)
# Z <- matrix(rep(1,n),n,q) # Z should not contain an intercept
Z <- matrix(rnorm(n),n,q)

# parameters
beta0 <- rnorm(p)
beta0
gamma0 <- rnorm(q)
gamma0
sig20 <- abs(rnorm(1,mean=1)) # NEW: scale param
sig20

# quantile level
alpha <- 0.75

# dist is one of "norm","t","chisq","lnorm"
genout <- gen_eps(n, dist = "norm", df = NULL, tau = alpha)
eps <- genout$eps
nu0 <- genout$nu0

# resposes
yy <- c(X %*% beta0 + sqrt(sig20)*exp(0.5 * Z %*% gamma0)*eps)
plot(yy,cex=0.3)
# random censoring
cc <- rnorm(n,mean=2.5,sd=1)
deltas <- yy<=cc
y <- yy
sum(1-deltas)/n
y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

# calculate the conditional posterior
numpoints <- 100
beta.seq <- seq(beta0-1,beta0+1,length.out = numpoints)
gamma.seq <- seq(gamma0-1,gamma0+1,length.out = numpoints)
sig2.seq <- seq(sig20-1,sig20+1,length.out = numpoints)
nu.seq <- seq(nu0-1,nu0+1,length.out = numpoints)
# Note: need to keep the sig2.seq range > 0 mostly
logel.seq <- matrix(rep(NA,4*numpoints),4,numpoints)
for (ii in 1:numpoints) {
  message("ii = ", ii)
  G <- qrls.evalG(y,X,Z,alpha,beta.seq[ii],gamma0/2,sig20,nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,beta.seq[ii],0.5*gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq[1,ii] <- logEL(omegas,epsilons,deltas)
  
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma.seq[ii]/2,sig20,nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,0.5*gamma.seq[ii],sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq[2,ii] <- logEL(omegas,epsilons,deltas)

  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0/2,sig2.seq[ii],nu0)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,0.5*gamma0,sig2.seq[ii])
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq[3,ii] <- logEL(omegas,epsilons,deltas)
  
  G <- qrls.evalG(y,X,Z,alpha,beta0,gamma0/2,sig20,nu.seq[ii])
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,0.5*gamma0,sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq[4,ii] <- logEL(omegas,epsilons,deltas)
}
logelmode1 <- plotEL(beta.seq, logel.seq[1,], beta0, NA, expression(beta))
logelmode2 <- plotEL(gamma.seq, logel.seq[2,], gamma0/2, NA, expression(gamma))
logelmode3 <- plotEL(sig2.seq, logel.seq[3,], sig20, NA, expression(sigma^2))
logelmode4 <- plotEL(nu.seq, logel.seq[4,], nu0, NA, expression(nu))

# mcmc
nsamples <- 5000
nburn <- 2000
theta.hat <- hlm.fit(y = y, X = X, W = cbind(1,Z)) # solve by hlm
# quick check 
rbind(true = beta0, est = theta.hat$beta) # beta
rbind(true = c(log(sig20),gamma0), est = theta.hat$gamma) # gamma

betaInit <- theta.hat$beta
gammaInit <- theta.hat$gamma[2]/2 # Note: divided by 2 since the models differ by this factor
sig2Init <- exp(theta.hat$gamma[1])
nuInit <- quantile(eps)
mwgSd <- c(0.1,0.1,0.1,0.1)
RvDoMcmc <- c(0,0,0,1)
system.time(
  # postout <- qrls.post(y,X,Z,nsamples,nburn,betaInit,gammaInit,sig2Init,mwgSd)
  postout <- qrls_cens.post_adapt(y,X,Z,deltas,alpha,
                                  nsamples,nburn,
                                  beta0,gamma0/2,sig20,nu0,
                                  mwgSd,RvDoMcmc)
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
lines(gamma.seq, norm_pdf(logel.marg2, gamma.seq),
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
