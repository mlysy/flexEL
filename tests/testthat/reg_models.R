# ---- mean regression ----

# location model
# mr_evalG is the same for non-censoring and censoring for mean regression
# X: nObs x nEqs matrix
# return: a nObs x nEqs matrix G
mr_evalG_R <- function(y, X, beta) {
  tX <- t(X)
  # yXb <- y - c(beta %*% tX)
  yXb <- y - c(X %*% beta)
  G <- sweep(tX, MARGIN = 2, yXb, `*`)
  return(t(G))
}

# location-scale model
mrls_evalG_R <- function(y, X, Z, beta, gamma, sig2) {
  nObs <- nrow(X)
  nBeta <- length(beta)
  nGamma <- length(gamma)
  G <- matrix(NaN, nObs, nBeta + nGamma + 1)
  eZg <- c(exp(-Z %*% gamma)) # e^{-z'gamma}
  yXbeZg <- c((y - X %*% beta)*eZg) # (y-x'beta)e^{-z'gamma}
  yXbeZg2 <- yXbeZg * yXbeZg # (y-x'beta)^2*e^{-2z'gamma}
  G[,1:nBeta] <- yXbeZg * eZg * X
  G[,nBeta+1:nGamma] <- (1-yXbeZg2/sig2) * Z
  G[,nBeta+nGamma+1] <- 1/sig2 * yXbeZg2 - 1;
  return(G)
}

# ---- quantile regression ----

# check function
rho_alpha <- function(u, alpha) {
  u * (alpha - (u <= 0))
}

# 1st derivative of check function
phi_alpha <- function(u, alpha) {
  # (u <= 0) - alpha
  alpha - (u <= 0)
}

# location model
# qr_evalG is the same for non-censoring and censoring for quantile regression
# X: nObs x nEqs matrix
# return: a nObs x nEqs matrix G
qr_evalG_R <- function(y, X, alpha, beta) {
  tX <- t(X) # tX is nEqs x nObs
  yXb <- y - c(beta %*% tX)
  pyXb <- phi_alpha(yXb,alpha)
  G <- sweep(tX, MARGIN = 2, pyXb, `*`)
  return(t(G))
}

# location-scale model
# Note: only works for 1 quantile now
qrls_evalG_R <- function(y, X, Z, alpha, beta, gamma, sig2, nu) {
  if (length(nu) > 1) stop("Only supports 1 quantile level at the moment.")
  nObs <- nrow(X)
  nBeta <- length(beta)
  nGamma <- length(gamma)
  G <- matrix(NaN, nObs, nBeta + nGamma + 2)
  eZg <- c(exp(-Z %*% gamma)) # e^{-z'gamma}
  yXbeZg <- c((y - X %*% beta)*eZg) # (y-x'beta)e^{-z'gamma}
  yXbeZg2 <- yXbeZg * yXbeZg # (y-x'beta)^2*e^{-2z'gamma}
  G[,1:nBeta] <- yXbeZg * eZg * X
  G[,nBeta+1:nGamma] <- (1-yXbeZg2/sig2) * Z
  G[,nBeta+nGamma+1] <- 1/sig2 * yXbeZg2 - 1;
  G[,nBeta+nGamma+2] <- phi_alpha(yXbeZg/sqrt(sig2)-nu, alpha)
  return(G)
}

# ---- quantile regression (smooth) ----

# smoothed indicator function 1(x <= 0)
ind_smooth_R <- function(x, s=10) {
  return(1/(1+exp(s*x)))
}

# 1st derivative of ind_smooth
ind1_smooth_R <- function(x, s=10) {
  return(-s*exp(s*x)/(1+exp(s*x))^2)
}

# 2nd derivative of ind_smooth
ind2_smooth_R <- function(x, s=10) {
  return((-s*s*exp(s*x)*(1+exp(s*x))+2*s*s*exp(2*s*x))/((1+exp(s*x))^3))
}

# smoothed check function
rho_smooth_R <- function(x, tau, s=10) {
  return(x*(tau-ind_smooth_R(x, s)))
}

# 1st derivate of rho_smooth
rho1_smooth_R <- function(x, tau, s=10) {
  return(tau-ind_smooth_R(x,s)-x*ind1_smooth_R(x,s))
}

# 2nd derivative of rho_smooth
rho2_smooth_R <- function(x, tau, s=10) {
  retval <- -2 * ind1_smooth_R(x,s) - x*ind2_smooth_R(x,s)
  return(retval)
}

# smoothed version of qr_evalG
qr_evalG_smooth_R <- function(y,X,tau,beta,s=10) {
  tX <- t(X)
  yXb <- y - c(beta %*% tX)
  pyXb <- rho1_smooth_R(yXb,tau,s)
  G <- sweep(tX, MARGIN = 2, pyXb, `*`)
  return(t(G))
}


qrls_evalG_smooth_R <- function(y, X, Z, tau, beta, gamma, sig2, nu, s=10) {
  nObs <- nrow(X)
  nBeta <- length(beta)
  nGamma <- length(gamma)
  G <- matrix(NaN, nObs, nBeta + nGamma + 2)
  eZg <- c(exp(-Z %*% gamma)) # e^{-z'gamma}
  yXbeZg <- c((y - X %*% beta)*eZg) # (y-x'beta)e^{-z'gamma}
  yXbeZg2 <- yXbeZg * yXbeZg # (y-x'beta)^2*e^{-2z'gamma}
  G[,1:nBeta] <- yXbeZg * eZg * X
  # G[,nBeta+1:nGamma] <- yXbeZg2 * Z
  # G[,nBeta+1:nGamma] <- (1-yXbeZg2) * Z
  G[,nBeta+1:nGamma] <- (1-yXbeZg2/sig2) * Z
  G[,nBeta+nGamma+1] <- 1/sig2 * yXbeZg2 - 1;
  G[,nBeta+nGamma+2] <- rho1_smooth_R(yXbeZg/sqrt(sig2)-nu, tau, s)
  return(G)
}

