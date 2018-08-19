# tests for quantile regression postSampler 
library(bayesEL)
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("gen_eps.R")
source("smoothEL.R")
source("mode-functions.R")

# ---- 1-d problem ----
n <- 300
mu <- 1
alpha <- 0.5
X <- matrix(rep(1,n), n, 1) # each row of X is one observation

genout <- gen_eps(n, dist = "norm", df = NULL, tau = alpha)
eps <- genout$eps
nu0 <- genout$nu0

quantile(eps, alpha)
# mu0 <- mu + qnorm(alpha, lower.tail = TRUE) # true param assuming alpha-tail of eps is 0
mu0 <- mu+nu0
mu0

# random censoring
# cc <- rnorm(n,mean=2,sd=1)
# yy <- c(X * mu) + eps
# deltas <- yy<=cc
# y <- yy
# sum(1-deltas)/n
# y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]
cc <- rnorm(n,mean=1.5,sd=1)
deltas <- eps<=cc
sum(1-deltas)/n
deltas <-eps<=cc
eps[!deltas] <- cc[!deltas]
y <- c(X * mu) + eps

# plot(yy,cex=.3,ylim = range(c(yy,y)))
# points(y,cex=.3,col='red')

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
mu.seq <- seq(-1+mu0,1+mu0,length.out = numpoints)
# mu.seq <- seq(-1+quantile(yy,alpha),1+quantile(yy,alpha),length.out = numpoints)
logel.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- qr.evalG(y,X,alpha,mu.seq[ii])
  epsilons <- y - c(X * mu.seq[ii])
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq[ii] <- logEL(omegas,epsilons,deltas)
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, quantile(y,alpha), expression(mu))

# smoothed version
numpoints <- 100
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints)
# mu.seq <- seq(-1+quantile(yy,alpha),1+quantile(yy,alpha),length.out = numpoints)
logel.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- qr.evalG.smooth_R(y,X,alpha,mu.seq[ii],s=10)
  epsilons <- y - c(X * mu.seq[ii])
  omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons,s=10)$omegas
  logel.seq[ii] <- logEL.smooth_R(omegas,epsilons,deltas,s=10)
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, quantile(y,alpha), expression(mu))

nsamples <- 10000
nburn <- 2000
library(quantreg)
betaInit <- c(rq(y ~ 1, tau = alpha, method = 'fn')$coefficients)
betaInit
sigs <- rep(0.2,1)
system.time(
  qrout <- qr_cens.post_adapt(y, X, deltas, alpha, nsamples, nburn, betaInit, sigs)
)
mu_chain <- qrout$beta_chain
mu_paccept <- qrout$paccept 
mu_paccept

# 10000+2000 with sigs = 0.2
# user  system elapsed 
# 93.837   0.404  94.887 

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

# # ---- 1-d problem with 2 quantile levels ---- 
# n <- 100
# mu <- 1
# alphas <- c(0.5,0.75)
# X <- matrix(rep(1,n), n, 1) # each row of X is one observation
# eps <- rnorm(n) # N(0,1) error term
# quantile(eps, alphas[1])
# quantile(eps, alphas[2])
# y <- c(X * mu) + eps
# mu01 <- mu + qnorm(alphas[1], lower.tail = TRUE) # true param assuming alpha-tail of eps is 0
# mu01
# mu02 <- mu + qnorm(alphas[2], lower.tail = TRUE) # true param assuming alpha-tail of eps is 0
# mu02
# 
# # conditional posterior
# numpoints <- 100
# mu.seq1 <- seq(-1+mu01,1+mu01,length.out = numpoints)
# logel.seq1 <- rep(NA,numpoints)
# for (ii in 1:numpoints) {
#   G <- qr.evalG(y,X,alphas,cbind(mu.seq1[ii],mu02))
#   logel.seq1[ii] <- logEL(G = G)
# }
# # logelmode1 <- plotEL(mu.seq1, logel.seq1, mu01, quantile(y,alphas[1]), expression(mu))
# # logelmode1
# 
# mu.seq2 <- seq(-1+mu02,1+mu02,length.out = numpoints)
# logel.seq2 <- rep(NA,numpoints)
# for (ii in 1:numpoints) {
#   G <- qr.evalG(y,X,alphas,cbind(mu01,mu.seq2[ii]))
#   logel.seq2[ii] <- logEL(G = G)
# }
# # logelmode2 <- plotEL(mu.seq2, logel.seq2, mu02, quantile(y,alphas[2]), expression(mu))
# # logelmode2
# 
# # marginal posterior
# Mu.seq <- as.matrix(expand.grid(mu.seq1, mu.seq2))
# logel.mat <- apply(Mu.seq, 1, function(mm) {
#   G <- qr.evalG(y,X,alphas, cbind(mm[1],mm[2]))
#   logEL(G)
# })
# logel.mat <- matrix(logel.mat, numpoints, numpoints)
# el.mat <- exp(logel.mat - max(logel.mat))
# logel.marg <- log(cbind(mu1 = rowSums(el.mat), mu2 = colSums(el.mat)))
# 
# library(quantreg)
# betaInit1 <- c(rq(y ~ 1, tau = alphas[1], method = 'fn')$coefficients)
# betaInit1
# betaInit2 <- c(rq(y ~ 1, tau = alphas[2], method = 'fn')$coefficients)
# betaInit2
# BetaInit <- cbind(betaInit1,betaInit2)
# # BetaInit <- cbind(mu01,mu02)
# nsamples <- 20000
# nburn <- 3000
# sigs1 <- rep(0.22,1)
# sigs2 <- rep(0.18,1)
# Sigs <- cbind(sigs1,sigs2)
# # if one param is not updated, then the lines should be the conditional ones
# RvDoMcmc <- cbind(1,0)
# k <- which(as.logical(RvDoMcmc))
# beta0 <- cbind(mu01,mu02)
# BetaInit[1,-k] <- beta0[1,-k]
# system.time(
#   qrout <- qr.post(y, X, alphas, nsamples, nburn, BetaInit, Sigs, RvDoMcmc)
# )
# 
# mu_chain <- qrout$Beta_chain
# mu_paccept <- qrout$paccept 
# mu_paccept
# plot(mu_chain[1,], xlab = 'mu', ylab = 'EL', type='l')
# plot(mu_chain[2,], xlab = 'mu', ylab = 'EL', type='l')
# 
# # overlay grid plot to histogram
# hist(mu_chain[1,],breaks=50,freq=FALSE,
#      xlab = expression(mu), main='')
# # marginal line
# lines(mu.seq1, norm_pdf(logel.marg[,1], mu.seq1),
#       cex=0.1, col = 'red', type='l')
# # conditional line
# # lines(mu.seq1, norm_pdf(logel.seq1, mu.seq1),
# #       cex=0.1, col = 'blue', type='l')
# abline(v=mu01, col='red')
# abline(v=mean(mu_chain[1,]), col='blue')
# legend('topright',legend=c(expression('grid plot & true param'),
#                            expression('sample mean')),
#        lty = c(1,1), col = c('red','blue'), cex = 0.6)
# 
# hist(mu_chain[2,],breaks=50,freq=FALSE,
#      xlab = expression(mu), main='')
# # marginal line
# lines(mu.seq2, norm_pdf(logel.marg[,2], mu.seq2),
#       cex=0.1, col = 'red', type='l')
# # conditional line
# # lines(mu.seq2, norm_pdf(logel.seq2, mu.seq2),
# #       cex=0.1, col = 'blue', type='l')
# abline(v=mu02, col='red')
# abline(v=mean(mu_chain[2,]), col='blue')
# legend('topright',legend=c(expression('grid plot & true param'),
#                            expression('sample mean')),
#        lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- 2-d problem (1 intercept, 1 slope) ----
n <- 200
p <- 2
alpha <- 0.5
X0 <- matrix(rep(1,n),n,1)
# X1 <- matrix(seq(-2,2,length.out = n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)

# eps <- rnorm(n) # N(0,1) error term
genout <- gen_eps(n, dist = "norm", df = NULL, tau = alpha)
eps <- genout$eps
nu0 <- genout$nu0

# beta_I <- rnorm(1)
# beta_S <- rnorm(1)
beta_I  <- 0.5
beta_S <- 1
beta0 <- c(beta_I, beta_S)
# yy <- c(X1 %*% beta_S) + eps 
# plot(X1,yy,cex=0.3)
# beta0 <- c(beta_I, beta_S)
# beta0

# random censoring
# cc <- rnorm(n,mean=1.5,sd=1)
# deltas <- yy<=cc
# y <- yy
# sum(1-deltas)/n
# y[as.logical(1-deltas)] <- cc[as.logical(1-deltas)]
cc <- rnorm(n,mean=1.35,sd=1)
deltas <- eps<=cc
sum(1-deltas)/n
eps[!deltas] <- cc[!deltas]
y <- beta0 + c(X1 %*% beta1) + eps 

# grid plot (conditionals: beta1|beta2 and beta2|beta1)
numpoints <- 100
beta1.seq <- seq(beta0[1]-1,beta0[1]+1,length.out = numpoints)
beta2.seq <- seq(beta0[2]-1,beta0[2]+1,length.out = numpoints)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  G <- qr.evalG(y,X,alpha,c(beta1.seq[ii],beta0[2]))
  epsilons <- y - c(X %*% c(beta1.seq[ii],beta0[2]))
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq[1,ii] <- logEL(omegas,epsilons,deltas)
  
  G <- qr.evalG(y,X,alpha,c(beta0[1],beta2.seq[ii]))
  epsilons <- y - c(X %*% c(beta0[1],beta2.seq[ii]))
  omegas <- omega.hat(G,deltas,epsilons)
  logel.seq[2,ii] <- logEL(omegas,epsilons,deltas)
}
logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))

# smoothed version
numpoints <- 100
beta1.seq <- seq(beta0[1]-1,beta0[1]+1,length.out = numpoints)
beta2.seq <- seq(beta0[2]-1,beta0[2]+1,length.out = numpoints)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  if (ii %% 10 == 0) message("ii = ", ii)
  G <- qr.evalG.smooth_R(y,X,alpha,c(beta1.seq[ii],beta0[2]),s=10)
  epsilons <- y - c(X %*% c(beta1.seq[ii],beta0[2]))
  omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons,s=10)$omegas
  logel.seq[1,ii] <- logEL.smooth_R(omegas,epsilons,deltas,s=10)
  
  G <- qr.evalG.smooth_R(y,X,alpha,c(beta0[1],beta2.seq[ii]),s=10)
  epsilons <- y - c(X %*% c(beta0[1],beta2.seq[ii]))
  omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons,s=10)$omegas
  logel.seq[2,ii] <- logEL.smooth_R(omegas,epsilons,deltas,s=10)
}
logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))

# calculate marginal posterior
Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
logel.mat <- apply(Beta.seq, 1, function(bb) {
  G <- qr.evalG(y,X,alpha,matrix(c(bb[1],bb[2]),2,1))
  epsilons <- y - c(X %*% matrix(c(bb[1],bb[2]),2,1))
  omegas <- omega.hat(G,deltas,epsilons)
  logEL(omegas,epsilons,deltas)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

# smoothed version
counter <- 0
Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
logel.mat <- apply(Beta.seq, 1, function(bb) {
  counter <<- counter + 1
  if (counter %% 10 == 0) message("ii = ", counter)
  G <- qr.evalG.smooth_R(y,X,alpha,c(bb[1],bb[2]),s=10)
  epsilons <- y - c(X %*% matrix(c(bb[1],bb[2]),2,1))
  omegas <- omega.hat.EM.smooth_R(G,deltas,epsilons,s=10)$omegas
  logEL.smooth_R(omegas,epsilons,deltas,s=10)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

library(quantreg)
betaInit <- c(rq(y ~ X1, tau = alpha, method = 'fn')$coefficients)
betaInit
nsamples <- 10000
nburn <- 3000
sigs <- c(0.15,0.15)
RvDoMcmc <- rbind(1,1)
k <- which(as.logical(RvDoMcmc))
betaInit[-k] <- beta0[-k]
system.time(
  qrout <- qr_cens.post_adapt(y, X, deltas, alpha, nsamples, nburn, 
                              betaInit, sigs, RvDoMcmc)
  # c(betaInit[1],beta0[2])
)
beta_chain <- qrout$beta_chain
beta_paccept <- qrout$paccept
beta_paccept
plot(beta_chain[1,], xlab = expression(beta[0]), ylab = 'EL', type='l')
plot(beta_chain[2,], xlab = expression(beta[1]), ylab = 'EL', type='l')
# overlay grid plot to histogram
# intercept
hist(beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta[0]),main='')
# marginal line
lines(beta1.seq, norm_pdf(logel.marg[,1], beta1.seq),
      cex=0.1, col = 'red', type='l')
# conditional line
lines(beta1.seq, norm_pdf(logel.seq[1,], beta1.seq),
      cex=0.1, col = 'blue', type='l')
abline(v=beta0[1], col='red')
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
abline(v=beta0[2], col='red')
abline(v=mean(beta_chain[2,]), col='blue')
legend('topright',legend=c(expression('true value'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# # ---- 2-d problem with 2 quantile levels (1 intercept, 1 slope) ----
# n <- 200
# alphas <- c(0.75,0.9)
# X0 <- matrix(rep(1,n),n,1)
# # X1 <- matrix(seq(-2,2,length.out = n),n,1)
# X1 <- matrix(rnorm(n),n,1)
# X <- cbind(X0,X1)
# 
# # eps <- rnorm(n) # N(0,1) error term
# genout <- gen_eps(n, dist = "norm", df = NULL, tau = alpha)
# eps <- genout$eps
# nu0 <- genout$nu0
# 
# beta_S <- 1.5
# beta_I1 <- qnorm(alphas[1], lower.tail = TRUE)
# beta01 <- c(beta_I1, beta_S)
# beta01
# beta_I2 <- qnorm(alphas[2], lower.tail = TRUE)
# beta02 <- c(beta_I2, beta_S)
# beta02
# y <- c(X1 %*% beta_S) + eps 
# plot(X1,y,cex=0.3)
# 
# # grid plot of conditionals: beta1|beta2 and beta2|beta1
# numpoints <- 100
# beta1.seq1 <- seq(beta01[1]-.5,beta01[1]+.5,length.out = numpoints)
# beta2.seq1 <- seq(beta01[2]-.5,beta01[2]+.5,length.out = numpoints)
# beta1.seq2 <- seq(beta02[1]-.5,beta02[1]+.5,length.out = numpoints)
# beta2.seq2 <- seq(beta02[2]-.5,beta02[2]+.5,length.out = numpoints)
# logel.seq <- matrix(rep(NA,4*numpoints),4,numpoints)
# for (ii in 1:numpoints) {
#   G <- qr.evalG(y,X,alphas,cbind(c(beta1.seq1[ii],beta01[2]),beta02))
#   logel.seq[1,ii] <- logEL(G)
#   G <- qr.evalG(y,X,alphas,cbind(c(beta01[1],beta2.seq1[ii]),beta02))
#   logel.seq[2,ii] <- logEL(G)
#   G <- qr.evalG(y,X,alphas,cbind(beta01,c(beta1.seq2[ii],beta02[2])))
#   logel.seq[3,ii] <- logEL(G)
#   G <- qr.evalG(y,X,alphas,cbind(beta01,c(beta02[1],beta2.seq2[ii])))
#   logel.seq[4,ii] <- logEL(G)
# }
# logelmode1 <- plotEL(beta1.seq1, logel.seq[1,], beta01[1], NA, expression(beta[0]))
# logelmode2 <- plotEL(beta2.seq1, logel.seq[2,], beta01[2], NA, expression(beta[1]))
# logelmode3 <- plotEL(beta1.seq2, logel.seq[3,], beta02[1], NA, expression(beta[0]))
# logelmode4 <- plotEL(beta2.seq2, logel.seq[4,], beta02[2], NA, expression(beta[1]))
# 
# library(quantreg)
# betaInit1 <- c(rq(y ~ X1, tau = alphas[1], method = 'fn')$coefficients)
# betaInit1
# betaInit2 <- c(rq(y ~ X1, tau = alphas[2], method = 'fn')$coefficients)
# betaInit2
# BetaInit <- cbind(betaInit1,betaInit2)
# 
# # MCMC
# nsamples <- 20000
# nburn <- 5000
# Sigs <- cbind(c(0.2,0.12),c(0.3,0.4))
# RvDoMcmc <- cbind(c(0,0),c(0,1))
# k <- which(as.logical(RvDoMcmc))
# beta0 <- cbind(beta01,beta02)
# BetaInit[-k] <- beta0[-k]
# system.time(
#   qrout <- qr.post(y, X, alphas, nsamples, nburn, BetaInit, Sigs, RvDoMcmc)
# )
# beta_chain <- qrout$Beta_chain
# beta_paccept <- qrout$paccept
# beta_paccept
# plot(beta_chain[1,], xlab = expression(beta[0]), ylab = 'EL', type='l')
# plot(beta_chain[2,], xlab = expression(beta[1]), ylab = 'EL', type='l')
# plot(beta_chain[3,], xlab = expression(beta[0]), ylab = 'EL', type='l')
# plot(beta_chain[4,], xlab = expression(beta[1]), ylab = 'EL', type='l')
# 
# # overlay conditional grid plot to histogram
# # intercept
# hist(beta_chain[1,],breaks=50,freq=FALSE,
#      xlab = expression(beta[0]),main='')
# lines(beta1.seq1, norm_pdf(logel.seq[1,], beta1.seq1),
#       cex=0.1, col = 'blue', type='l')
# abline(v=beta01[1], col='red')
# abline(v=mean(beta_chain[1,]), col='blue')
# legend('topright',legend=c(expression('true param'),
#                            expression('sample mean')),
#        lty = c(1,1), col = c('red','blue'), cex = 0.6)
# # slope
# hist(beta_chain[2,],breaks=50,freq=FALSE,
#      xlab = expression(beta[1]),main='')
# lines(beta2.seq1, norm_pdf(logel.seq[2,], beta2.seq1),
#       cex=0.1, col = 'blue', type='l')
# abline(v=beta01[2], col='red')
# abline(v=mean(beta_chain[2,]), col='blue')
# legend('topright',legend=c(expression('true param'),
#                            expression('sample mean')),
#        lty = c(1,1), col = c('red','blue'), cex = 0.6)
# # intercept
# hist(beta_chain[3,],breaks=50,freq=FALSE,
#      xlab = expression(beta[0]),main='')
# lines(beta1.seq2, norm_pdf(logel.seq[3,], beta1.seq2),
#       cex=0.1, col = 'blue', type='l')
# abline(v=beta02[1], col='red')
# abline(v=mean(beta_chain[3,]), col='blue')
# legend('topright',legend=c(expression('true param'),
#                            expression('sample mean')),
#        lty = c(1,1), col = c('red','blue'), cex = 0.6)
# # slope
# hist(beta_chain[4,],breaks=50,freq=FALSE,
#      xlab = expression(beta[1]),main='')
# lines(beta2.seq2, norm_pdf(logel.seq[4,], beta2.seq2),
#       cex=0.1, col = 'blue', type='l')
# abline(v=beta02[2], col='red')
# abline(v=mean(beta_chain[4,]), col='blue')
# legend('topright',legend=c(expression('true param'),
#                            expression('sample mean')),
#        lty = c(1,1), col = c('red','blue'), cex = 0.6)
# 
