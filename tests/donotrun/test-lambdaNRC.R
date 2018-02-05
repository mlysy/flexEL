#--- check that lambdaNRC.R is working properly ---------------------------------
require(bayesEL)
source("~/bayesEL/tests/testthat/el-utils.R")
source("~/bayesEL/tests/donotrun/mle-check.R")

#---- ## checking correctness of lambdaNRCens ## ----

# 1-d problem 
p <- 1
N <- 200
mu <- 1
z <- rnorm(N, mean = mu) # true lifetime variable
c <- rnorm(N, mean = 2*mu) # censoring variable 
delta <- z <= c
sum(!delta)/N # censored percentage
y <- z 
y[!delta] <- c[!delta] # observed lifetime
X <- matrix(rep(p,N),p,N) # p x n matrix
ws0 <- rep(1/N,N)
# order data here
ord <- order(y)
y <- y[ord]
X <- matrix(X[,ord],p,N)
G <- mr.evalG_R(y, X, mu)
ws0 <- ws0[ord]
qs <- get_qs(ws0)

# optimization in C++
lambdahat <- lambdaNRC(G,qs)

# optimization in R
lambdahatout <- lambdaNRC_R(G,qs)
lambdahat_R <- lambdahatout$lambda

# difference of the two implementation 
lambdahat - lambdahat_R 

par(mfrow=c(1,1))
curve(sapply(x, QfunCens, G = G), from = lambdahat-1, to = lambdahat+1)
abline(v = lambdahat, col='red')
abline(v = lambdahat_R, col='blue')

# 2-dim problem : 
N <- 200
beta <- c(1,2)
X <- rbind(rep(1,N),rnorm(N))
z <- beta %*% X
c <- rnorm(N, mean = 2*mean(z)) # censoring variable 
delta <- z <= c
sum(!delta)/N # censored percentage
y <- z 
y[!delta] <- c[!delta] # observed lifetime
ws0 <- rep(1/N,N)
# order data here
ord <- order(y)
y <- y[ord]
X <- matrix(X[,ord],p,N)
G <- mr.evalG_R(y, X, beta)
ws0 <- ws0[ord]
qs <- get_qs(ws0)

# optimization in C++
lambdahat <- lambdaNRC(G,qs)

# optimization in R
lambdahatout <- lambdaNRC_R(G,qs)
lambdahat_R <- lambdahatout$lambda

# difference of the two implementation 
lambdahat - lambdahat_R 

# derivative test 
require(numDeriv)

# grad Q(lambdahat) = 0 (sometimes unreliable)
Qfcens <- function(lambda) QfunCens(lambda, G)
grad(Qfcens, lambdahat)
grad(Qfcens, lambdahat_R)

# visual test
mle.check(loglik = Qfcens, theta.mle = lambdahat)
mle.check(loglik = Qfcens, theta.mle = lambdahat_R)
