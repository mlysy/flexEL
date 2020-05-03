n <- 20
p <- 5
q <- 3
X <- replicate(p, rnorm(n))
Z <- replicate(q, rnorm(n))
beta <- rnorm(p)
gamma <- rnorm(q)
sig2 <- 1.5
y <- c(X %*% beta + exp(Z %*% gamma)) + rnorm(n) # with N(0,1) error term
max_iter <- 100
rel_tol <- 1e-5
mrls_evalG(y,X,Z,beta,gamma,sig2)
