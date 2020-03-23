n <- 20
p <- 4
X <- replicate(p, rnorm(n))
beta0 <- rnorm(p)
y <- c(X %*% beta0) + rnorm(n) # with N(0,1) error term
mr_evalG(y, X, beta0)
