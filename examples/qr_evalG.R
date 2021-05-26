## simulate some data ##
n <- 20
p <- 2
X <- replicate(p-1, rnorm(n))
X <- cbind(1, X)
beta0 <- c(1, 2)
y <- c(X %*% beta0) + rnorm(n)

## calculate G matrix given data and certain parameter values ##
# a single quantile level (with continuity correction)
alpha <- 0.5
beta <- c(1, 2)
qr_evalG(y, X, alpha, beta, s = 1)

# multiple quantile levels (with continuity correction)
alpha <- c(0.25, 0.75)
Beta <- cbind(c(0.5, 2), c(1.5, 2)) # each column corresponds to one quantile level
qr_evalG(y, X, alpha, Beta, s = 1)
