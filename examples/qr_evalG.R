n <- 20
p <- 4
X <- replicate(p, rnorm(n))
X[1,] <- rep(1,p)
beta <- rnorm(p)
alpha <- runif(1)
y <- c(X %*% beta) + rnorm(n) # with N(0,1) error term
qr_evalG(y, X, alpha, beta)
