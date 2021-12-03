## simulate some data ##
n <- 20
p <- 2
q <- 2
X <- replicate(p-1, rnorm(n))
X <- cbind(1, X) # X may have an intercept term, if so, it should be explicit
Z <- replicate(q, rnorm(n)) # Z shall not have an intercept term
beta0 <- c(1, 2)
gamma0 <- c(0.5, 0.25)
sig20 <- 0.5
y <- c(X %*% beta0 + sqrt(sig20)*exp(Z %*% gamma0)*rnorm(n)) # with N(0,1) error term

## calculate G matrix given data and certain parameter values ##
# a single quantile level (with continuity correction)
alpha <- 0.5
beta <- c(1, 2)
gamma <- c(0.5, 0.25)
sig2 <- 0.5
nu <- 0
qrls_evalG(y, X, Z, alpha, beta, gamma, sig2, nu, s = 1)
# multiple quantile levels (with continuity correction)
alpha <- c(0.25, 0.75)
beta <- c(1, 2)
gamma <- c(0.5, 0.25)
sig2 <- 0.5
nu <- c(-0.5, 0.5)
qrls_evalG(y, X, Z, alpha, beta, gamma, sig2, nu, s = 1)
