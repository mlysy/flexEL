#--- check that lambdaNRC.R is working properly ---------------------------------
require(bayesEL)
source("~/bayesEL/tests/testthat/el-utils.R")
source("~/bayesEL/tests/donotrun/mle-check.R")

# TODO: modify the following for lambdaNRC

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

# TODO: modify the following multi-dim problem : 
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

# wraper
Qfcens <- function(lambda) QfunCens(lambda, G)

mle.check(loglik = Qfcens, theta.mle = lambdahat)

mle.check(loglik = Qfcens, theta.mle = lambdahat_R)

#---- Remove later: the following are copied from test-lambdaNR.R ----

# 1-d problem
N <- 10 # number of observations
m <- 1 # number of dimensions
y <- rnorm(N*m)
X <- matrix(rep(1,N))
G <- t(matrix(y, m, N))

# optimization in C++
lambda0 <- rnorm(m)
lambdahat <- lambdaNR(G = G)

# by implementation in R
lambdahat_Rout <- lambdaNR_R(G = G)
lambdahat_R <- lambdahat_Rout$lambda

# difference of the two implementation 
lambdahat - lambdahat_R 

par(mfrow=c(1,1))
curve(exp(sapply(x, Qfun, G = G)), from = lambdahat-1, to = lambdahat+1)
abline(v = lambdahat, col='red')
abline(v = lambdahat_R, col='blue')

# multi-dimensional problem
N <- 10 # number of observations
m <- 3 # number of dimensions
y <- rnorm(N*m)
X <- matrix(rep(1,N))
G <- t(matrix(rnorm(N*m), m, N))

# optimization in C++
lambda0 <- rnorm(m)
lambdahat <- lambdaNR(G = G)

# by implementation in R
lambdahat_Rout <- lambdaNR_R(G = G)
lambdahat_R <- lambdahat_Rout$lambda # output is a row vector

# difference of the two implementation 
lambdahat - lambdahat_R 

# derivative test 
require(numDeriv)

# grad Q(lambdahat) = 0 (sometimes unreliable)
Qf <- function(lambda) Qfun(lambda, G)
grad(Qf, lambdahat)
grad(Qf, lambdahat_R)

# visual mode check (more reliable)
mle.check(loglik = Qf, theta.mle = lambdahat)

# visual mode check (more reliable)
mle.check(loglik = Qf, theta.mle = lambdahat_R)

