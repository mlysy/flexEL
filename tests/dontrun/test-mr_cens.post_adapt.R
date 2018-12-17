# tests for mean regression postSampler 
library(bayesEL)
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("gen_eps.R")
source("logEL_EMAC_R.R")
source("smoothEL.R")

# ---- 1-d problem ----
n <- 200
mu0 <- rnorm(1,0,1)
# X <- matrix(rnorm(n,0,1),n,1) # each row of X is one observation
X <- matrix(rep(1,n), n, 1)
# X <- matrix(sample(15:25,n,replace = TRUE),n,1)
# eps <- rnorm(n) # N(0,1) error term

# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "norm", df = NULL)

# yy <- c(X * mu0) + eps
# # random censoring
# cc <- rnorm(n,mean=1,sd=1)
# deltas <- yy<=cc
# y <- yy
# sum(1-deltas)/n
# y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]

# random censoring
cc <- rnorm(n,mean=1.35,sd=1)
deltas <- eps<=cc
sum(1-deltas)/n
eps[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]
y <- c(X * mu0) + eps

# fix censoring time
# pct <- 0.75
# idx <- order(yy)[round(n*pct)]
# cc <- yy[idx]
# deltas <- yy<=cc
# sum(1-deltas)/n # percentage censored
# y <- yy
# y[as.logical(1-deltas)] <- cc

# grid plot
numpoints <- 100
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints)
logel.seq <- rep(NA,numpoints)
adjust <- FALSE
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- mr.evalG(y,X,mu.seq[ii])
  if (adjust) G <- adjG_R(G)
  epsilons <- y - c(X * mu.seq[ii])
  if (adjust) {
    omegas <- omega.hat_R(G,deltas,epsilons,adjust = adjust)
    logel.seq[ii] <- logEL_R(omegas,epsilons,deltas,adjust = adjust)
  }
  else {
    omegas <- omega.hat(G,deltas,epsilons)
    logel.seq[ii] <- logEL(omegas,epsilons,deltas)
  }
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, mean(y), expression(mu))

# Acceleration with fitted lm
logel_acc.seq <- rep(NA,numpoints) 
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- mr.evalG(y,X,mu.seq[ii])
  epsilons <- y - c(X * mu.seq[ii])
  logel_acc.seq[ii] <- logEL_EMAC_R(G,epsilons,deltas)
}
logelmode_acc <- plotEL(mu.seq, logel_acc.seq, mu0, mean(y), expression(mu))

nsamples <- 10000
nburn <- 2000
betaInit <- c(lm(y ~ X-1)$coefficients) 
betaInit
sigs <- rep(0.1,1)
system.time(
  # mrout <- postCens_R(Gfun=mr.evalG_R,nThe=1,nBet=1,nGam=0,
  #                       y=y,X=X,deltas=deltas,thetaInit=betaInit,
  #                       nsamples=nsamples,nburn=nburn,
  #                       mwgSds=sigs,adjust = TRUE)
  mrout <- mr_cens.post_adapt(y, X, deltas, nsamples, nburn, betaInit, 
                              mwgSd = sigs, DoAdapt = 1)
)

mu_chain <- mrout$beta_chain
mu_paccept <- mrout$paccept 
mu_paccept
plot(mu_chain[1,], xlab = 'mu', ylab = 'EL', type='l')

# overlay grid plot to histogram
hist(mu_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(mu), main='')
lines(mu.seq, norm_pdf(logel.seq, mu.seq),
      cex=0.1, col = 'red', type='l')
abline(v=mu0, col='red')
abline(v=mean(mu_chain), col='blue')
legend('topright',legend=c(expression('true value'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- 1-d problem with 1 slope param ----
n <- 500
X <- matrix(rnorm(n),n,1)
beta_S <- 1
# eps <- rnorm(n) # N(0,1) error term

# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "norm", df = NULL)

# random censoring
cc <- rnorm(n,mean=1.35,sd=1)
delta <- eps<=cc
# delta <- rep(1,n)
sum(1-delta)/n
eps[as.logical(1-delta)] <- cc[as.logical(1-delta)]
y <- c(X %*% beta_S) + eps 
plot(X,y,cex=0.3)

# grid plot
numpoints <- 500
beta.seq <- seq(beta_S-.5,beta_S+.5,length.out = numpoints)
# beta.seq <- seq(0.98,1.07,length.out = numpoints)
logel.seq <- matrix(rep(NA,numpoints),1,numpoints)

for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- mr.evalG(y,X,beta.seq[ii])
  epsilons <- y - c(X %*% beta.seq[ii])
  omegas <- omega.hat(G,delta,epsilons)
  logel.seq[1,ii] <- logEL(omegas,epsilons,deltas)
}
logelmode <- plotEL(beta.seq, logel.seq[1,], beta_S, NA, expression(beta), legend.loc = 'topleft') # there is indeed jumps

# ---- 2-d problem (1 intercept, 1 slope) ----
n <- 200
p <- 2
X1 <- matrix(rnorm(n),n,1)
# X1 <- matrix(sample(15:25,n,replace = TRUE),n,1)
X <- cbind(1,X1)
# eps <- rnorm(n) # N(0,1) error term

# dist is one of "norm","t","chisq","lnorm"
eps <- gen_eps(n, dist = "norm", df = NULL)

beta_I <- 1
beta_S <- 1.5
# beta_I <- rnorm(1)
# beta_S <- rnorm(1)
beta0 <- c(beta_I, beta_S)

# random censoring
# yy <- beta_I + c(X1 %*% beta_S) + eps 
# cc <- rnorm(n,mean=3,sd=1)
# deltas <- yy<=cc
# y <- yy
# sum(1-deltas)/n
# y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]
# plot(X1,yy,cex=0.3)

# random censoring
cc <- rnorm(n,mean=1.35,sd=1)
deltas <- eps<=cc
sum(1-deltas)/n
eps[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]
y <- beta_I + c(X1 %*% beta_S) + eps 
plot(X1,y,cex=0.3)

# grid plot of conditionals: beta1|beta2 and beta2|beta1
numpoints <- 100
beta1.seq <- seq(beta0[1]-.5,beta0[1]+.5,length.out = numpoints)
# beta1.seq <- seq(1.1,1.15,length.out = numpoints)
beta2.seq <- seq(beta0[2]-.5,beta0[2]+.5,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
adjust <- FALSE
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- mr.evalG(y,X,c(beta1.seq[ii],beta0[2]))
  epsilons <- y - c(X %*% c(beta1.seq[ii],beta0[2]))
  omegas <- omega.hat(G,deltas,epsilons)
  # omegas <- omega.hat.EM_R(G,deltas,epsilons,rel_tol = 1e-3,verbose = TRUE)$omegas
  logel.seq[1,ii] <- logEL(omegas,epsilons,deltas)

  G <- mr.evalG(y,X,c(beta0[1],beta2.seq[ii]))
  epsilons <- y - c(X %*% c(beta0[1],beta2.seq[ii]))
  if (adjust) {
    G <- adjG_R(G)
    omegas <- omega.hat_R(G,deltas,epsilons,adjust = adjust)
    logel.seq[2,ii] <- logEL_R(omegas,epsilons,deltas,adjust = adjust)
  }
  else {
    omegas <- omega.hat(G,deltas,epsilons)
    # omegas <- omega.hat.EM_R(G,deltas,epsilons,rel_tol = 1e-3)$omegas
    logel.seq[2,ii] <- logEL(omegas,epsilons,deltas)
  }
}
par(mfrow=c(1,2))
logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]),
                     legend.loc = 'topright', cex.lab = 1.5, cex.axis = 1.2, cex = 1.2)
logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]),
                     legend.loc = 'topleft', cex.lab = 1.5, cex.axis = 1.2, cex = 1.2)
par(mfrow=c(1,1))

# smoothed censored logEL
numpoints <- 100
beta1.seq <- seq(beta0[1]-1,beta0[1]+1,length.out = numpoints)
# beta2.seq <- seq(1.45,1.54,length.out = numpoints)
beta2.seq <- seq(beta0[2]-1,beta0[2]+1,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
s <- 10
for (ii in 1:numpoints) {
  if (ii %% 1 == 0) message("ii = ", ii)
  
  logel.seq[1,ii] <- -mr_cens.neglogEL.smooth_R(y,X,deltas,c(beta1.seq[ii],beta0[2]),s)
  logel.seq[2,ii] <- -mr_cens.neglogEL.smooth_R(y,X,deltas,c(beta0[1],beta2.seq[ii]),s)
  
  # G <- mr.evalG(y,X,c(beta1.seq[ii],beta0[2]))
  # epsilons <- y - c(X %*% c(beta1.seq[ii],beta0[2]))
  # omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons)$omegas
  # logel.seq[1,ii] <- logEL.smooth_R(omegas,epsilons,deltas)
  # 
  # G <- mr.evalG(y,X,c(beta0[1],beta2.seq[ii]))
  # epsilons <- y - c(X %*% c(beta0[1],beta2.seq[ii]))
  # omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons,s=10)$omegas
  # logel.seq[2,ii] <- logEL.smooth_R(omegas,epsilons,deltas,s=10)
}
logelmode1.smooth <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2.smooth <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))

# check the fitted acc
logel_acc.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  message("ii = ", ii)
  G <- mr.evalG(y,X,c(beta1.seq[ii],beta0[2]))
  epsilons <- y - c(X %*% c(beta1.seq[ii],beta0[2]))
  logel_acc.seq[1,ii] <- logEL_EMAC_R(G,epsilons,deltas)
  
  G <- mr.evalG(y,X,c(beta0[1],beta2.seq[ii]))
  epsilons <- y - c(X %*% c(beta0[1],beta2.seq[ii]))
  logel_acc.seq[2,ii] <- logEL_EMAC_R(G,epsilons,deltas)
}
logelmode1_acc <- plotEL(beta1.seq, logel_acc.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2_acc <- plotEL(beta2.seq, logel_acc.seq[2,], beta0[2], NA, expression(beta[1]))

# calculate marginal posterior
# Note: this may take a while to calculate
Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
adjust <- FALSE
counter <- 0
logel.mat <- apply(Beta.seq, 1, function(bb) {
  counter <<- counter+1
  if (counter %% 10 == 0) message("counter = ", counter)
  G <- mr.evalG(y,X,c(bb[1],bb[2]))
  if (adjust) G <- adjG_R(G)
  epsilons <- y - c(X %*% c(bb[1],bb[2]))
  # omegas <- omega.hat(G,deltas,epsilons)
  # logEL(omegas,epsilons,deltas)
  if (adjust) {
    omegas <- omega.hat_R(G,deltas,epsilons,adjust)
    logEL_R(omegas,epsilons,deltas,adjust)
  }
  else {
    # omegas <- omega.hat(G,deltas,epsilons)
    # logEL(omegas,epsilons,deltas)
    omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons,s=10)$omegas
    logEL.smooth_R(omegas,epsilons,deltas,s=10)
  }
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))
plot(beta2.seq,logel.marg[,2],type='l')

nsamples <- 10000 # 20000 worked
nburn <- 5000
betaInit <- c(lm(y ~ X1)$coefficients)
betaInit
sigs <- c(0.1,0.02)
# Note: For matching the conditional grid plot and histogram, pay attention to 
#   their inital values, the fixed params need to be the same.
RvDoMcmc <- c(1,1)
k <- which(as.logical(RvDoMcmc))
betaInit[-k] <- beta0[-k] # fix other params at their true values
system.time(
  # mrout <- postCens_R(Gfun=mr.evalG_R,nThe=2,nBet=2,nGam=0,
  #                       y=y,X=X,deltas=deltas,thetaInit=betaInit,
  #                       nsamples=nsamples,nburn=nburn,
  #                       mwgSds=sigs,adjust = TRUE)
  mrout <- mr_cens.post_adapt(y, X, deltas, nsamples, nburn, betaInit, sigs, RvDoMcmc)
)
beta_chain <- mrout$beta_chain
beta_paccept <- mrout$paccept
beta_paccept
plot(beta_chain[1,], xlab = expression(beta[0]), ylab = 'EL', type='l')
plot(beta_chain[2,], xlab = expression(beta[1]), ylab = 'EL', type='l')

# overlay marginal grid plot to histogram
# intercept
hist(beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]),main='')
# marginal line
lines(beta1.seq, norm_pdf(logel.marg[,1], beta1.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(beta1.seq, norm_pdf(logel.seq[1,], beta1.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0[1],col='red')
abline(v=mean(beta_chain[1,]), col='blue')
legend('topright',legend=c(expression('true value'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
# slope
hist(beta_chain[2,],breaks=50,freq=FALSE,
     xlab = expression(beta[1]),main='')
# marginal line
lines(beta2.seq, norm_pdf(logel.marg[,2], beta2.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(beta2.seq, norm_pdf(logel.seq[2,], beta2.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0[2],col='red')
abline(v=mean(beta_chain[2,]), col='blue')
legend('topright',legend=c(expression('true value'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

