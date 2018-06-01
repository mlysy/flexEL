# source("../../R/hlm.R")
library(bayesEL)
library(numDeriv)
source("gen_eps.R")
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")

n <- 300
p <- 2
q <- 1
# X1 <- rnorm(n)
X1 <- sample(1:25,n,replace = TRUE)
# X2 <- rbinom(n,1,0.3)
# X <- cbind(X1,X2)
X <- cbind(1,X1)

# W <- matrix(rnorm(n*q),n,q)
beta0 <- c(1,1)
beta0
gamma0 <- .5
gamma0
sig20 <- 1.5
sig20
# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "norm", df = NULL)
yy <- c(X %*% beta0) + sqrt(sig20)*exp(0.5*X1*gamma0)*eps
cc <- rnorm(n, mean=25)
delta <- yy <= cc
sum(1-delta)/length(delta)
y <- yy
y[!delta] <- cc[!delta]

# --------------------------------------------------------------------------- 
# plot conditional likelihood
numpoints <- 100
Z <- matrix(X1)
deltas <- delta
adjust <- TRUE

beta.seq1 <- seq(-.5+beta0[1],.5+beta0[1],length.out = numpoints)
logel.seq1 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  message("ii = ", ii)
  G <- mrls.evalG(y,X,Z,c(beta.seq1[ii],beta0[2]),0.5*gamma0,sig20)
  if (adjust) G <- adjG_R(G)
  epsilons <- evalEpsilonsLS(y,X,Z,c(beta.seq1[ii],beta0[2]),0.5*gamma0,sig20)
  if (adjust) {
    omegas <- omega.hat_R(G,deltas,epsilons,adjust)
    logel.seq1[ii] <- logEL_R(omegas,epsilons,deltas,adjust)
  }
  else {
    omegas <- omega.hat(G,deltas,epsilons)
    logel.seq1[ii] <- logEL(omegas,epsilons,deltas)
  }
}
logelmode1 <- plotEL(beta.seq1, logel.seq1, beta0[1], NA, expression(beta[0]))

beta.seq2 <- seq(-.5+beta0[2],.5+beta0[2],length.out = numpoints)
logel.seq2 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  message("ii = ", ii)
  G <- mrls.evalG(y,X,Z,c(beta0[1],beta.seq2[ii]),0.5*gamma0,sig20)
  if (adjust) G <- adjG_R(G)
  epsilons <- evalEpsilonsLS(y,X,Z,c(beta0[1],beta.seq2[ii]),0.5*gamma0,sig20)
  if (adjust) {
    omegas <- omega.hat_R(G,deltas,epsilons,adjust)
    logel.seq2[ii] <- logEL_R(omegas,epsilons,deltas,adjust)
  }
  else {
    omegas <- omega.hat(G,deltas,epsilons)
    logel.seq2[ii] <- logEL(omegas,epsilons,deltas)
  }
}
logelmode2 <- plotEL(beta.seq2, logel.seq2, beta0[2], NA, expression(beta[1]))

gamma.seq <- seq(-.5+gamma0[1],.5+gamma0[1],length.out = numpoints)
logel.seq3 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,beta0,0.5*gamma.seq[ii],sig20)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,0.5*gamma.seq[ii],sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq3[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode3 <- plotEL(0.5*gamma.seq, logel.seq3, 0.5*gamma0, NA, expression(gamma[1]))

gamma.seq2 <- seq(-1+gamma0[2],1+gamma0[2],length.out = numpoints)
logel.seq4 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,beta0,0.5*c(gamma0[1],gamma.seq2[ii]),sig20)
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,0.5*c(gamma0[1],gamma.seq2[ii]),sig20)
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq4[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode4 <- plotEL(0.5*gamma.seq2, logel.seq4, 0.5*gamma0[2], NA, expression(gamma[2]))

sig2.seq <- seq(-1+sig20,1+sig20,length.out = numpoints)
logel.seq5 <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mrls.evalG(y,X,Z,beta0,0.5*gamma0,sig2.seq[ii])
  epsilons <- evalEpsilonsLS(y,X,Z,beta0,0.5*gamma0,sig2.seq[ii])
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq5[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode5 <- plotEL(sig2.seq, logel.seq5, sig20, NA, expression(sigma^2))

# ---- use hlm ----
hlmout <- hlm(y,delta,X,X,rel_tol = 1e-5)
hlmout$conv
# quick check 
rbind(true = beta0, est = hlmout$coef$beta) # beta
rbind(true = c(log(sig20),gamma0), est = hlmout$coef$gamma) # gamma
beta.hlm <- hlmout$coef$beta
gamma.hlm <- 0.5*hlmout$coef$gamma[2]
sig2.hlm <- exp(hlmout$coef$gamma[1])
# rbind(true = c(log(sig20),gamma0), est = hlmout$coef$gamma) # gamma

# quantile esimate of eps is the empirical quantile of the residuals
# nu.hlm <- quantile(evalEpsilonsLS(y,X,W,beta.hlm,gamma.hlm,sig2.hlm), prob=alpha)
# bootstrap the residuals for the CI of it

# ---- use bayesEL ----
nsamples <- 500
nburn <- 0
betaInit <- beta.hlm
gammaInit <- gamma.hlm
sig2Init <- sig2.hlm
mwgSd <- rep(0.1,p+q)
rvDoMcmc <- rep(1,p+q)
# belout <- mrls_cens.post_adapt(y,X,X1,delta,nsamples,nburn,
#                                betaInit,gammaInit,sig2Init,
#                                mwgSd,rvDoMcmc)
# theta_chain <- belout$theta_chain
# theta_accept <- belout$paccept
# theta_accept
mrlsout <- postCens_R(Gfun=mrls.evalG_R,nThe=4,nBet=2,nGam=1,
                      y=y,X=X,Z=matrix(X1),deltas=delta,
                      thetaInit=c(betaInit,gammaInit,sig2Init),
                      nsamples=nsamples,nburn=nburn,
                      mwgSds=mwgSd,adjust = TRUE)
theta_chain <- mrlsout$theta_chain
theta_accept <- mrlsout$paccept
theta_accept



# mixing of the chains
plot(theta_chain[4,],type = 'l')

hist(theta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]), main='')
abline(v=beta0[1],col='red')
abline(v=mean(theta_chain[1,]),col='blue')
abline(v=hlmout$coef$beta[1],col='green')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean'),
                           expression('hlm sol')),
       lty = c(1,1), col = c('red','blue','green'), cex = 0.6)

hist(theta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[2]), main='')
abline(v=beta0[2],col='red')
abline(v=mean(theta_chain[2,]),col='blue')
abline(v=hlmout$coef$beta[2],col='green')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean'),
                           expression('hlm sol')),
       lty = c(1,1), col = c('red','blue','green'), cex = 0.6)

hist(theta_chain[3,],breaks=50,freq=FALSE,
     xlab = expression(gamma[1]), main='')
abline(v=gamma0[1]/2,col='red')
abline(v=mean(theta_chain[3,]),col='blue')
abline(v=hlmout$coef$gamma[2]/2,col='green')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean'),
                           expression('hlm sol')),
       lty = c(1,1), col = c('red','blue','green'), cex = 0.6)

hist(theta_chain[4,],breaks=50,freq=FALSE,
     xlab = expression(gamma[2]), main='')
abline(v=gamma0[2]/2,col='red')
abline(v=mean(theta_chain[4,]),col='blue')
abline(v=hlmout$coef$gamma[3]/2,col='green')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean'),
                           expression('hlm sol')),
       lty = c(1,1), col = c('red','blue','green'), cex = 0.6)

hist(theta_chain[5,],breaks=50,freq=FALSE,
     xlab = expression(sigma^2), main='')
abline(v=sig20,col='red')
abline(v=mean(theta_chain[5,]),col='blue')
abline(v=exp(hlmout$coef$gamma[1]),col='green')
legend('topright',legend=c(expression('true param'),
                           expression('sample mean'),
                           expression('hlm sol')),
       lty = c(1,1), col = c('red','blue','green'), cex = 0.6)

# ---- plot the lines ----
plot(y=y,x=X1,cex=.3)
abline(a=beta0[1],b=beta0[2],col='red')
abline(a=mean(theta_chain[1,]),b=mean(theta_chain[2,]),col='blue')
abline(a=beta.hlm[1],b=beta.hlm[2],col='green')

