n <- 20
p <- 5
q <- 3
X <- replicate(p, rnorm(n))
Z <- replicate(q, rnorm(n))
alpha <- 0.5
beta <- rnorm(p)
gamma <- rnorm(q)
nu <- 0.5
sig2 <- 1.5
y <- c(X %*% beta + sqrt(sig2)*exp(Z %*% gamma)*rnorm(n)) # with multiplicative N(0,1) error
sp <- 10
qrls_evalG_smooth(y,X,Z,alpha,beta,gamma,sig2,nu,sp)
