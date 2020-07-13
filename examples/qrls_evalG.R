n <- 20
p <- 5
q <- 3
X <- replicate(p, rnorm(n))
Z <- replicate(q, rnorm(n))
alpha <- runif(1)
beta <- rnorm(p)
gamma <- rnorm(q)
nu <- rnorm(1)
sig2 <- 1.5
y <- c(X %*% beta + sqrt(sig2)*exp(Z %*% gamma)*rnorm(n)) # with multiplicative N(0,1) error

# no continuity correction
qrls_evalG(y, X, Z, alpha, beta, gamma, sig2, nu)

# with continuity correction
sp <- 10
qrls_evalG(y, X, Z, alpha, beta, gamma, sig2, nu, sp)
