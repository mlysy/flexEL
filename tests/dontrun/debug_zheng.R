n <- 1e4
x <- rep(1, n)
y <- rnorm(n)
tau <- 0.1
alpha <- 1
quantile(y, tau)

tt <- seq(from = -2, to = 2, length.out = 500)
ll <- sapply(tt, function(t) {neg_logel_qr_sfun(t, y, cbind(x), tau, alpha)})
plot(x = tt, y = ll, cex = .3)
tt[which.min(ll)]
abs(tt[which.min(ll)]-quantile(y, tau))

tt <- seq(from = -2, to = 2, length.out = 500)
ll <- sapply(tt, function(t) {neg_logel_qr_zfun(t, y, cbind(x), tau, alpha)})
plot(x = tt, y = ll, cex = .3)
tt[which.min(ll)]
abs(tt[which.min(ll)]-quantile(y, tau))

tt <- seq(from = -2, to = 2, length.out = 500)
ll <- sapply(tt, function(t) {neg_logel_qr_smooth(t, y, cbind(x), tau, 1/alpha)})
plot(x = tt, y = ll, cex = .3)
tt[which.min(ll)]
abs(tt[which.min(ll)]-quantile(y, tau))

tt <- seq(from = -2, to = 2, length.out = 500)
ll <- sapply(tt, function(t) {neg_logel_qr_orig(t, y, cbind(x), tau)})
plot(x = tt, y = ll, cex = .3)
tt[which.min(ll)]
abs(tt[which.min(ll)]-quantile(y, tau))

loss <- rep(NA, length(y))
for (ii in 1:length(y)) {
  loss[ii] <- sum(rho_alpha(y-y[ii], tau))
}
plot(x = y, y = loss, cex = .3)
y[which.min(loss)]
quantile(y, tau)
