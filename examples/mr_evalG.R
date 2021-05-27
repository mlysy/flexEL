## simulate some data ##
n <- 20
p <- 2
X <- replicate(p-1, rnorm(n))
X <- cbind(1, X)
beta0 <- c(1, 2)
y <- c(X %*% beta0) + rnorm(n) # with N(0,1) error term

## calculate G matrix given data and certain parameter values ##
beta <- c(1, 2)
mr_evalG(y, X, beta)
