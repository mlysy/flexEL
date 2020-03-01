# tests for mean regression postSampler 
library(flexEL)
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("gen_eps.R")

# ---- 1-d problem ----
n <- 100
mu0 <- 1
X <- matrix(rep(1,n), n, 1) # each row of X is one observation
eps <- rnorm(n) # N(0,1) error term
y <- c(X * mu0) + eps

# gird plot
numpoints <- 100
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints)
logel.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mr.evalG(y,X,mu.seq[ii])
  omegas <- omega.hat(G)
  logel.seq[ii] <- logEL(omegas)
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, mean(y), expression(mu))

# adjusted G matrix (TODO: in R now)
numpoints <- 100
mu.seq <- seq(-.5+mu0,.5+mu0,length.out = numpoints)
logel.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mr.evalGadj_R(y,X,mu.seq[ii])
  omegas <- omega.hat(G)
  logel.seq[ii] <- logEL(omegas)
}
logelmode <- plotEL(mu.seq, logel.seq, mu0, mean(y), expression(mu))

nsamples <- 20000
nburn <- 5000
betaInit <- c(lm(y ~ 1)$coefficients) # TODO: does not work ???
# betaInit <- round(c(lm(y ~ 1)$coefficients), digits = 6)
betaInit
sigs <- rep(0.25,1)
qrout <- mr.post(y, X, nsamples, nburn, betaInit, sigs)
mu_chain <- qrout$Beta_chain
mu_paccept <- qrout$paccept 
mu_paccept
plot(mu_chain[1,], xlab = 'mu', ylab = 'EL', type='l')

# overlay gird plot to histogram
hist(mu_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(mu), main='')
lines(mu.seq, norm_pdf(logel.seq, mu.seq),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode, col='red')
abline(v=mean(mu_chain), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# check R version 
qrout.R <- post_R(Gfun = mr.evalG_R, nThe = 1, nBet = 1, nGam = 0, 
                  y, X, nsamples, nburn, matrix(betaInit), matrix(sigs))
mu_chain.R <- qrout.R$Theta_chain
mu_paccept.R <- qrout.R$paccept 
mu_paccept.R
plot(mu_chain.R[1,], xlab = 'mu', ylab = 'EL', type='l')
hist(mu_chain.R[1,],breaks=50,freq=FALSE,
     xlab = expression(mu), main='')
lines(mu.seq, norm_pdf(logel.seq, mu.seq),
      cex=0.1, col = 'red', type='l')
abline(v=logelmode, col='red')
abline(v=mean(mu_chain), col='blue')
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- 2-d problem (1 intercept, 1 slope) ----
n <- 500
p <- 2
X0 <- matrix(rep(1,n),n,1)
X1 <- matrix(rnorm(n),n,1)
X <- cbind(X0,X1)
eps <- rnorm(n) # N(0,1) error term
beta_intercept <- 1
beta_slope <- 1.5
y <- 1 + c(X1 %*% beta_slope) + eps 
beta0 <- c(beta_intercept, beta_slope)
plot(X1,y,cex=0.3)

# grid plot of conditionals: beta1|beta2 and beta2|beta1
numpoints <- 100
beta1.seq <- seq(beta0[1]-.5,beta0[1]+.5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-.5,beta0[2]+.5,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  G <- mr.evalG(y,X,c(beta1.seq[ii],beta0[2]))
  omegas <- omega.hat(G)
  logel.seq[1,ii] <- logEL(omegas)
  G <- mr.evalG(y,X,c(beta0[1],beta2.seq[ii]))
  omegas <- omega.hat(G)
  logel.seq[2,ii] <- logEL(omegas)
}
logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))

# use adjusted G (TODO: in R now)
numpoints <- 100
beta1.seq <- seq(beta0[1]-.5,beta0[1]+.5,length.out = numpoints)
beta2.seq <- seq(beta0[2]-.5,beta0[2]+.5,length.out = numpoints)
beta.seq <- cbind(beta1.seq,beta2.seq)
logel.seq <- matrix(rep(NA,2*numpoints),2,numpoints)
for (ii in 1:numpoints) {
  G <- mr.evalGadj_R(y,X,c(beta1.seq[ii],beta0[2]))
  omegas <- omega.hat(G)
  logel.seq[1,ii] <- logEL(omegas)
  G <- mr.evalGadj_R(y,X,c(beta0[1],beta2.seq[ii]))
  omegas <- omega.hat(G)
  logel.seq[2,ii] <- logEL(omegas)
}
logelmode1 <- plotEL(beta1.seq, logel.seq[1,], beta0[1], NA, expression(beta[0]))
logelmode2 <- plotEL(beta2.seq, logel.seq[2,], beta0[2], NA, expression(beta[1]))

# calculate marginal posterior
Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
logel.mat <- apply(Beta.seq, 1, function(bb) {
  G <- mr.evalG(y,X,c(bb[1],bb[2]))
  omegas <- omega.hat(G)
  logEL(omegas)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

# calculate marginal posterior with adjusted G matrix (TODO: in R now)
Beta.seq <- as.matrix(expand.grid(beta1.seq, beta2.seq))
logel.mat <- apply(Beta.seq, 1, function(bb) {
  G <- mr.evalGadj_R(y,X,c(bb[1],bb[2]))
  omegas <- omega.hat(G)
  logEL(omegas)
})
logel.mat <- matrix(logel.mat, numpoints, numpoints)
el.mat <- exp(logel.mat - max(logel.mat))
logel.marg <- log(cbind(beta1 = rowSums(el.mat), beta2 = colSums(el.mat)))

nsamples <- 10000
nburn <- 3000
betaInit <- c(lm(y ~ X1)$coefficients)
betaInit
sigs <- rep(0.12,2)
# Note: For matching the conditional grid plot and histogram, pay attention to 
#   their inital values, the fixed params need to be the same.
RvDoMcmc <- c(1,1)
# k <- which(as.logical(RvDoMcmc))
# betaInit[-k] <- beta0[-k] # fix other params at their true values
system.time(
  # qrout <- mr.post(y, X, nsamples, nburn, betaInit, sigs, RvDoMcmc)
  qrout <- post_R(Gfun = mr.evalGadj_R, nThe = 2, nBet = 2, nGam = 0,
                  y, X, nsamples, nburn, matrix(betaInit), matrix(sigs))
)
# beta_chain <- qrout$Beta_chain
beta_chain <- qrout$Theta_chain
beta_paccept <- qrout$paccept
beta_paccept
plot(beta_chain[1,], xlab = expression(beta[0]), ylab = 'EL', type='l')
plot(beta_chain[2,], xlab = expression(beta[1]), ylab = 'EL', type='l')

# overlay marginal gird plot to histogram
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
legend('topright',legend=c(expression('grid plot'),
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
legend('topright',legend=c(expression('grid plot & mode'),
                           expression('sample mean')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)

# ---- Check conditional posterior ---------------------------------------------
n <- 200
p <- sample(2:10,1) # choose a random number of params
beta <- rnorm(p)
X <- cbind(1,matrix(rnorm((p-1)*n),n,p-1))
y <- c(X %*% beta + rnorm(n))

numpoints <- 100
k <- sample(p,1) # choose a random param to update, fix all else
beta.seq <- seq(-.5+beta[k],.5+beta[k],length.out = numpoints)
logel.seq <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  beta.temp <- beta
  beta.temp[k] <- beta.seq[ii]
  G <- mr.evalG(y,X,beta.temp)
  logel.seq[ii] <- logEL(G = G)
}
# logelmode <- plotEL(beta.seq, logel.seq, beta[k], NA, expression(beta))

betaInit <- c(lm(y ~ X-1)$coefficients)
betaInit
nsamples <- 20000
nburn <- 3000
sigs <- rep(0.2,p)
RvDoMcmc <- rep(0,p)
RvDoMcmc[k] <- 1
betaInit[-k] <- beta[-k] # fix other params at their true values
system.time(
  qrout <- mr.post(y, X, nsamples, nburn, betaInit, sigs, RvDoMcmc)
)
beta_chain <- qrout$Beta_chain
beta_paccept <- qrout$paccept
beta_paccept
hist(beta_chain[k,],breaks=50,freq=FALSE,
     xlab = expression(beta),main='')
lines(beta.seq, norm_pdf(logel.seq, beta.seq),
      cex=0.1, col = 'red', type='l')

# another way to get conditional posterior 
# (NOTE: not the right way to go since G would be different)
y.tilda <- y - c(X[,-k] %*% beta[-k])
numpoints <- 100
beta.seq.tilda <- seq(-.5+beta[k],.5+beta[k],length.out = numpoints)
logel.seq.tilda <- rep(NA,numpoints)
for (ii in 1:numpoints) {
  G <- mr.evalG(y.tilda,matrix(X[,k]),beta.seq.tilda[ii])
  logel.seq.tilda[ii] <- logEL(G = G)
}
# logelmode <- plotEL(beta.seq.tilda, logel.seq.tilda, beta[k], NA)

betaInit <- c(lm(y ~ X[,k]-1)$coefficients)
betaInit
nsamples <- 20000
nburn <- 3000
sigs <- 0.2
# Note: here G has only one column
system.time(
  qrout <- mr.post(y.tilda, matrix(X[,k]), nsamples, nburn, betaInit, sigs)
)
beta_chain <- qrout$Beta_chain
beta_paccept <- qrout$paccept
beta_paccept
hist(beta_chain[1,],breaks=50,freq=FALSE,
     xlab = expression(beta),main='')
lines(beta.seq.tilda, norm_pdf(logel.seq.tilda, beta.seq.tilda),
      cex=0.1, col = 'blue', type='l')

# the two "different" conditional curves are not the same
plot(beta.seq, norm_pdf(logel.seq, beta.seq), type='l', col='red')
lines(beta.seq.tilda, norm_pdf(logel.seq.tilda, beta.seq.tilda), type='l',col='blue')
legend('topright',legend=c(expression('use y'),
                           expression('use y.tilda')),
       lty = c(1,1), col = c('red','blue'), cex = 0.6)
