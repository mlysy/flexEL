# ---- smoothed check functions ----

#' Equation (3) in Zheng 2011
#' @param x Numeric vector.
#' @param tau Quantile (probabilities).
#' @param alpha Tuning parameter for smoothing.
sfun <- function(x, tau, alpha) {
  tau*x + alpha*log(1+exp(-x/alpha))
}

#' First derivative of \code{sfun}
sfun1 <- function(x, tau, alpha) {
  tau - 1/(exp(x/alpha) + 1)
}

#' Shifted \code{sfun} so that the argmin is at 0
zfun <- function(x, tau, alpha) {
  temp <- alpha*log(1/tau-1)
  tau*(x+temp) + alpha*log(1+exp(-(x+temp)/alpha))
}

#' First derivative of \code{zfun}
zfun1 <- function(x, tau, alpha) {
  temp <- x + alpha*log(1/tau-1)
  tau - 1/(exp(temp/alpha) + 1)
}

# ---- visualize and test the check functions ----

tau <- 0.75
alpha <- 1
curve(rho_alpha(x, tau), from = -10, to = 10, ylab ='')
curve(sfun(x, tau, alpha), col = 'red', add = TRUE)
curve(rho_smooth_R(x, tau, s = 1/alpha), add = TRUE, col = 'green')
curve(zfun(x, tau, alpha), add = TRUE, col = 'blue')
legend('topleft', lty = c(1,1,1,1), col = c('black', 'red', 'green', 'blue'),
       legend = c("Koenker & Bassett", "Zheng", "flexEL", "Shifted Zheng"),
       title = paste0("alpha = ", alpha, " tau = ", tau))

# check that the minimum of the two functions are at x0 and 0 respectively
x0 <- alpha*log(1/tau-1)
x0

nlm(sfun, p = 0, tau = tau, alpha = alpha)
nlm(zfun, p = x0, tau = tau, alpha = alpha)

# check that zfun is shifted sfun by x0
x <- seq(from = -1, to = 1, length.out=10)
testthat::expect_equal(sfun(x+x0, tau, alpha), zfun(x, tau, alpha))

# check that sfun1 is the correct derivative of sfun
x <- seq(from = -1, to = 1, length.out=10)
testthat::expect_equal(numDeriv::grad(sfun, x, tau = tau, alpha = alpha), 
                       sfun1(x, tau, alpha))

# check that zfun1 is the correct derivative of zfun
x <- seq(from = -1, to = 1, length.out=10)
testthat::expect_equal(numDeriv::grad(zfun, x, tau = tau, alpha = alpha), 
                       zfun1(x, tau, alpha))

# ---- calculate G matrix ----

qr_evalG_sfun <- function(y, X, beta, tau, alpha) {
  tX <- t(X) # tX is nEqs x nObs
  yXb <- y - c(beta %*% tX)
  pyXb <- sfun1(yXb, tau, alpha)
  G <- sweep(tX, MARGIN = 2, pyXb, `*`)
  return(t(G))
}

qr_evalG_zfun <- function(y, X, beta, tau, alpha) {
  tX <- t(X) # tX is nEqs x nObs
  yXb <- y - c(beta %*% tX)
  pyXb <- zfun1(yXb, tau, alpha)
  G <- sweep(tX, MARGIN = 2, pyXb, `*`)
  return(t(G))
}

n <- 500
p <- 2
X <- replicate(p, rnorm(n))
X[,1] <- rep(1,n)
beta <- c(1, 2)
eps <- rnorm(n)
y <- c(X %*% beta) + eps # with N(0,1) error term
plot(x = X[,2], y = y, cex = .3)
lines(x = X[,2], y = c(X %*% beta))

x <- rnorm(n, mean = 0, sd = 1)
eps <- rnorm(n)
y <- 2*x + eps

plot(x, y, cex = .3)
lines(x, y = 2*x)

b1 <- seq(from = 1, to = 3, length.out = 100)
# ll <- sapply(bb, function(b) {neg_logel_qr_sfun(c(1, b), y, X, 0.75, 10)})
ll <- sapply(b1, function(b) {neg_logel_qr_sfun(c(0.6745, b), y, cbind(1,x), 0.75, 1)})
plot(x = b1, y = ll, type = 'l')
abline(v = b1[which.min(ll)], col = 'red')
b1[which.min(ll)]

b0 <- seq(from = -1, to = 2, length.out = 100)
# ll <- sapply(b0, function(b) {neg_logel_qr_sfun(c(b, 2), y, X, 0.75, 10)})
ll <- sapply(b0, function(b) {neg_logel_qr_sfun(c(b, 2), y, cbind(1,x), 0.75, 10)})
plot(x = b0, y = ll, type = 'l')
abline(v = b0[which.min(ll)], col = 'red')
b0[which.min(ll)]

res <- nlm(neg_logel_qr_sfun, p = c(0.8,1.2), X = X, y = y, tau = 0.25, alpha = 10)
res

b1 <- seq(from = 1, to = 3, length.out = 100)
ll <- sapply(b1, function(b) {neg_logel_qr_zfun(c(0.6745, b), y, cbind(1,x), 0.75, 10)})
plot(x = b1, y = ll, type = 'l')
abline(v = b1[which.min(ll)], col = 'red')
b1[which.min(ll)]

b0 <- seq(from = 0, to = 3, length.out = 100)
ll <- sapply(bb, function(b) {neg_logel_qr_zfun(c(b, 2), y, cbind(1,x), 0.75, 100)})
plot(x = bb, y = ll, type = 'l')
abline(v = bb[which.min(ll)], col = 'red')
bb[which.min(ll)]

res <- nlm(neg_logel_qr_zfun, p = beta, X = X, y = y, tau = 0.25, alpha = 10)
res

b1 <- seq(from = 1, to = 3, length.out = 100)
# ll <- sapply(bb, function(b) {neg_logel_qr_sfun(c(1, b), y, X, 0.75, 10)})
ll <- sapply(b1, function(b) {neg_logel_qr_smooth(c(0.6745, b), y, cbind(1,x), 0.75, 10)})
plot(x = b1, y = ll, type = 'l')
abline(v = b1[which.min(ll)], col = 'red')
b1[which.min(ll)]

b0 <- seq(from = -1, to = 2, length.out = 100)
# ll <- sapply(b0, function(b) {neg_logel_qr_sfun(c(b, 2), y, X, 0.75, 10)})
ll <- sapply(b0, function(b) {neg_logel_qr_smooth(c(b, 2), y, cbind(1,x), 0.75, 100)})
plot(x = b0, y = ll, type = 'l')
abline(v = b0[which.min(ll)], col = 'red')
b0[which.min(ll)]

res <- nlm(neg_logel_qr_smooth_orig, p = beta, X = X, y = y, tau = 0.25, alpha = 10)
res

lines(x = x, y = X %*% beta)
lines(x = X[,2], y = X %*% res$estimate, col = 'blue')


