n <- 1000
p <- 2
# X <- replicate(p-1, rnorm(n))
X <- sample(seq(from = 1, to = 10, length.out = 100), n, replace = TRUE)
X <- cbind(1, X) # n_obs x n_eqs
beta <- c(1, 2)
alpha <- c(0.25,0.75)
y <- c(X %*% beta) + exp(0.1*X[,2])*rnorm(n) # with N(0,1) error term
G <- qr_evalG(y, X, alpha, cbind(beta, beta), s = 10)
G

lfun <- function(beta, y, X, alpha, gel) {
  G <- qr_evalG(y, X, alpha, matrix(beta, ncol = length(alpha)), s = 1)
  -gel$logel(G)
}

gel <- GenEL$new(n, 2*p)
gel$supp_adj <- TRUE

tmp <- nlm(lfun, c(0.5, 1.25, 0.75, 1.55), y, X, alpha, gel)
tmp
Beta_hat <- matrix(tmp$estimate, ncol = length(alpha))

plot(x = X[,2], y = y, cex = .3)
lines(x = X[,2], y = X %*% Beta_hat[,1], col = 'red')
lines(x = X[,2], y = X %*% Beta_hat[,2], col = 'red')

# ----- right-censored logel with qr ----

## simulate some data ##
n <- 200
p <- 2
X <- replicate(p-1, rnorm(n))
X <- cbind(1, X)
beta0 <- c(1, 2)
y <- c(X %*% beta0) + rnorm(n)
c <- rnorm(2*mean(y), n = n)
delta <- y <= c
sum(!delta)
y_obs <- y
y_obs[!delta] <- c[!delta]

## calculate G matrix given data and certain parameter values ##
# a single quantile level (with continuity correction)
alpha <- 0.5
beta <- c(1, 2)
qr_evalG(y, X, alpha, beta, s = 1)

qr_neglogcel <- function(beta, y_obs, X, delta, alpha, cel) {
  G <- flexEL::qr_evalG(y, X, alpha, beta, s = 1)
  eps <- y - X %*% beta
  -cel$logel(G, delta, eps)
}

cel <- flexEL::CensEL$new(n, p)
cel$supp_adj <- TRUE
cel$smooth <- TRUE
cel$smooth_s <- 1

qr_neglogcel(c(1, 2), y_obs, X, delta, 0.5, cel)
nlm(qr_neglogcel, c(1,2), y_obs, X, delta, 0.5, cel)

# ---- right-censored logel with mr -----

mr_neglogcel <- function(beta, y_obs, X, cel) {
  G <- flexEL::mr_evalG(y, X, beta)
  -gel$logel(G)
}

cel <- flexEL::CensEL$new(n, p)
cel$supp_adj <- TRUE
cel$smooth <- TRUE
cel$smooth_s <- 1

nlm(mr_neglogcel, c(1,2), y_obs, X, cel)
