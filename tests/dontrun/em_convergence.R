library(bayesEL)
library(numDeriv)
source("hlm-functions.R")
source("gen_eps.R")
source("../testthat/el-utils.R")
source("../testthat/el-rfuns.R")
source("../testthat/el-model.R")
source("mode-functions.R")
source("logEL_EMAC_R.R")

# with random data
n <- 200
p <- 2
max_iter <- 500
rel_tol <- 1e-7
G <- matrix(rnorm(n*p), n, p)
deltas <- rep(1,n)
numcens <- sample(round(n/2),1)
censinds <- sample(n,numcens)
deltas[censinds] <- 0
1-sum(deltas)/n
epsilons <- rnorm(n)
omegahat <- omega.hat.EM_R(G, deltas, epsilons, adjust = FALSE, 
                           max_iter = max_iter, rel_tol = rel_tol, 
                           verbose = FALSE, dbg= TRUE)

plot(omegahat$logels, cex=.3)
len <- length(omegahat$logels)
len

logEL(omegahat$omegas,epsilons,deltas)

# ok, since we do not know what's the value l_i converges to, this does not work well
# idx <- len
# els <- omegahat$logels[1:(idx-1)]
# plot(els, cex=.3)
# optval <- omegahat$logels[idx]
# dif1 <- diff(log(abs(els-optval)))
# plot(exp(dif1), cex=.3)

# use another way, this shows the correct plot
# dif2 <- diff(omegahat$logels)[2:(len-1)]/diff(omegahat$logels)[1:(len-2)]
# plot(dif2, cex=.3) # coverge to a const => linear convergence
dif2 <- diff(log(diff(omegahat$logels)))
cseq <- exp(dif2)
plot(cseq, cex=.5) # coverge to a const => linear convergence 
# Obs: the higer pct of censoring, the larger the covergence rate c
# idx <- 1/(1:length(cseq))
# mod <- lm(cseq~idx+I(idx^2), weights = 1:length(cseq))
# points(predict(mod), cex=.5, col="red")

# pred <- predict(lm(cseq[1:49]~idx))
# plot(exp(pred), cex=.5)
# points(cseq[1:49],col='red', cex=.5)

omegahat.acc <- logEL_EMAC_R(G,epsilons,deltas,rel_tol = 1e-3, dbg = TRUE)
idx <- 1/(1:length(omegahat.acc$logels))
length(idx)
# fit <- lm(omegahat.acc$logels~idx+I(idx^2))
# pred <- predict(fit)
# plot(omegahat.acc$logels,ylim = range(c(pred,omegahat.acc$logels)))
# points(pred,col="blue")

nfit <- ceiling(1/2*length(idx))
if (nfit < 5) nfit <- 5
idx_fit <- idx[1:nfit]
idx_pre <- idx[(nfit+1):length(idx)]
fit <- lm(omegahat.acc$logels[1:nfit]~idx_fit+I(idx_fit^2), weights = 1:nfit)
# ws <- 1:nfit
# ws <- ws/sum(ws)
sres <- sum(resid(fit)^2)
sres
sres < 1e-2
pred <- predict(fit, newdata = data.frame(idx_fit=idx))
plot(pred, ylim=range(c(omegahat.acc$logels,pred)))
points(omegahat.acc$logels, col='red')
legend("bottomright", fill=c("black","red"), legend = c("pred","true"))

# compare the values
cbind(coef(fit)[1],omegahat.acc$logel,logEL(omegahat$omegas,epsilons,deltas))
coef(fit)[1]-logEL(omegahat$omegas,epsilons,deltas)
coef(fit)[1]+coef(fit)[2]*1/nfit-logEL(omegahat$omegas,epsilons,deltas)
omegahat.acc$logel-logEL(omegahat$omegas,epsilons,deltas)

# Let's generate a perfectly linear convergence sequecence, see how it looks like
c <- 0.5
n <- 100
lseq <- rep(NA,n)
lseq[1] <- 1
lseq[2] <- 2
for (ii in 3:n) {
  lseq[ii] <- (c+1)*lseq[ii-1]-c*lseq[ii-2]
}
plot(lseq, cex=.5)

# check the fitted acceleration 
accout <- logEL_EMAC2_R(G, epsilons, deltas,
                        max_iter = 50, rel_tol = 1e-3, fit_iter = 5, dbg=TRUE)
plot(accout$logels, cex=.5)
xx <- 1:length(log(diff(accout$logels)))
cest <- coef(lm(log(diff(accout$logels)) ~ xx-1))
exp(cest)
accout$logels[1] + 1/(1-exp(cest))*(accout$logels[2]-accout$logels[1])

dif2 <- diff(log(diff(accout$logels)))
cseq <- exp(dif2)
plot(cseq, cex=.5) 

plot(accout$logels.acc, cex=.5)
plot(accout$ccs, cex=.5)
