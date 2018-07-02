#' Fit the normal heteroscedastic linear model with censoring.
#'
#' @name hlm
#' @aliases hlm_loglik
#' @param y Response vector.
#' @param delta Censoring indicator vector (1: observed, 0: censored).
#' @param X Covariate matrix for the mean term.  Must include the intercept if required.
#' @param W Covariate matrix for the variance term.  Must include the intercept if required.
#' @param max_iter Maximum number of iterations of ECM algorithm (see Details).
#' @param rel_tol Maximum relative parameter difference between successive iterations (see Details).
#' @param multi_cycle Whether or not to compute the E-step after both CM steps, or only after the second one.  The former converges in fewer steps, but each step is longer.
#' @return For \code{hlm_loglik}, the loglikelihood function.  For \code{hlm}, a list with elements
#' \describe{
#'   \item{\code{coef}}{The vector \code{(beta_hat, gamma_h)} giving the MLE estimates for the model.}
#'   \item{\code{vcov}}{The estimated variance for the MLE.}
#'   \item{\code{loglik}}{The value of the loglikelihood function evaluated at the MLE.}
#'   \item{\code{niter}}{The number of ECM iterations.}
#' }
#' @details The HLM model is given by
#' \preformatted{
#' z_i = x_i'beta + exp(.5 w_i'gamma) * eps_i,    eps_i ~iid N(0,1),
#'
#' y_i = min(z_i, c_i),   delta_i = 1[z_i <= c_i].
#' }
#' Maximum likelihood estimates of \code{(beta, gamma)} are obtained by an ECM method extending Smyth (1989).  That is, without censoring the parameters can be sequentially updated using \code{lm} for \code{beta}, and \code{glm(family = "Gamma", link = "log")} for \code{gamma}.
#'
#' The function \code{hlm} does the fitting, and \code{hlm_loglik} computes the loglikelihood function.
#' @export
hlm <- function(y, delta, X, W,
                max_iter = 1e3, rel_tol = 1e-5, multi_cycle = FALSE) {
  # Helper functions
  # E[Z|Z>a]
  f <- function(a) {
    return(dnorm(a)/pnorm(-a))
  }
  # E[Z^2|Z>a]
  g <- function(a) {
    return(1+a*dnorm(a)/pnorm(-a))
  }
  # relative error
  rel_err <- function(theta1, theta2) {
    abs(theta1-theta2)/abs((theta1+theta2)/2)
  }
  # loglikelihood function
  loglik <- function(theta) {
    hlm_loglik(beta = theta[1:px], gamma = theta[px+1:pw], y, delta, X, W)
  }
  ## loglik_orig <- function(theta) {
  ##   t <- exp(y)
  ##   beta <- theta[1:px]
  ##   gamma <- theta[px+1:pw]
  ##   mu <- c(X%*%beta)
  ##   sigma <- exp(c(W%*%gamma)/2)
  ##   ll_orig <- ifelse(delta, dlnorm(t, meanlog = mu, sdlog = sigma, log = TRUE),
  ##                     plnorm(t, meanlog = mu, sdlog = sigma, lower.tail = FALSE, log.p = TRUE))
  ##   return(sum(ll_orig))
  ## }
  # some constants
  X <- as.matrix(X)
  W <- as.matrix(W)
  # covariate dimensions
  px <- ncol(X)
  pw <- ncol(W)
  # covariate names
  if(is.null(colnames(X))) colnames(X) <- 1:px
  if(is.null(colnames(W))) colnames(W) <- 1:pw
  XW_names <- c(paste0("X", colnames(X)), paste0("W", colnames(W)))
  # covariates for censored observations
  Xc <- as.matrix(X[!delta,])
  Wc <- as.matrix(W[!delta,])
  # (added: in case that R drops dimention if there's only one cencored obs.)
  if (sum(!delta)==1) {
    Xc <- t(Xc)
    Wc <- t(Wc)
  }
  yc <- y[!delta]
  # T_tilde and U_tilde.  some of these values will never be updated
  T <- y
  U <- y^2
  ## if(debug) browser()
  # initialize the model parameters
  # (removed)
  # betat <- coef(lm.fit(x = X, y = y))
  # (added)
  glm_conv <- TRUE
  betat <- tryCatch(expr = coef(lm.fit(x = X, y = y)),
                    warning = function(w) {glm_conv <<- FALSE; coef(lm.fit(x = X, y = y))},
                    error = function(e) {rep(NA,px)})
  # (removed)
  # gammat <- coef(glm.fit(x = W, y = (y-X%*%betat)^2,
  #                        family = Gamma(link = "log")))
  # (added)
  gammat <- tryCatch(expr = coef(glm.fit(x = W, 
                                         y = (y-X%*%betat)^2, 
                                         family = Gamma(link = "log"))),
                     warning = function(w) {glm_conv <<- FALSE; coef(glm.fit(x = W, 
                                                                            y = (y-X%*%betat)^2, 
                                                                            family = Gamma(link = "log")))},
                     error = function(e) {rep(NA,pw)})
  # (added)
  if (anyNA(gammat)) {
    out <- list(conv = FALSE,
                coef = list(beta = betat, gamma = gammat),
                niter = 0)
    return(out)
  }
  # (added)
  conv <- TRUE
  # main loop
  for(ii in 1:max_iter) {
    # E-step
    mut <- c(Xc %*% betat)
    Wg <- c(W %*% gammat)
    sigmat <- exp(Wg[!delta]/2)
    zt <- (yc - mut)/sigmat
    sft <- sigmat * f(zt)
    T[!delta] <- sft + mut
    U[!delta] <- sigmat^2 * g(zt) + 2*mut * sft + mut^2
    # M-step: beta
    beta <- coef(lm.wfit(x = X, y = T, w = exp(-Wg)))
    if(multi_cycle) {
      # recompute E-step
      mut <- c(Xc %*% beta)
      zt <- (yc - mut)/sigmat
      sft <- sigmat * f(zt)
      T[!delta] <- sft + mut
      U[!delta] <- sigmat^2 * g(zt) + 2*mut * sft + mut^2
    }
    # M-step: gamma
    Xb <- c(X %*% beta)
    R <- U - 2*T*Xb + Xb^2
    # (removed)
    # gamma <- coef(glm.fit(x = W, y = R, family = Gamma(link = "log")))
    # (added: tryCatch)
    glm_conv <- TRUE
    gamma <- tryCatch(expr = coef(glm.fit(x = W, y = R, family = Gamma(link = "log"))),
                      warning = function(w) {glm_conv <<- FALSE; coef(glm.fit(x = W, y = R, family = Gamma(link = "log")))},
                      error = function(e) {rep(NA,pw)})
    # (added: stop when glm returns an error)
    if (anyNA(gamma)) break
    # relative error
    theta_rel <- rel_err(c(beta, gamma), c(betat, gammat))
    # update
    betat <- beta
    gammat <- gamma
    if(max(theta_rel) < rel_tol) break
  }
  niter <- ii # (added)
  # (added: to handle glm error or warning)
  if((!glm_conv) || anyNA(gamma) || (ii == max_iter && max(theta_rel) > rel_tol)) {
  # (removed)
  # if(ii == max_iter && max(theta_rel) > rel_tol) {
    # (removed: for the ease of optimCheck)
    # stop("ECM algorithm did not converge.")
    # (added)
    conv <- FALSE
    out <- list(conv = conv,
                coef = list(beta = betat, gamma = gammat),
                niter = niter)
    return(out)
  }
  # return some relevant values
  theta_hat <- c(betat, gammat)
  names(theta_hat) <- XW_names
  ll_max <- loglik(theta_hat)
  ## ll_max_orig <- loglik_orig(theta_hat)
  var_hat <- solve(-hessian(func=loglik,x=theta_hat))
  # (added: tryCatch)
  # var_hat <- tryCatch(expr=solve(-hessian(func=loglik,x=theta_hat)),
  #                     error=function(e){matrix(NA,ncol=length(theta_hat),
  #                                              nrow=length(theta_hat))})
  colnames(var_hat) <- XW_names
  rownames(var_hat) <- XW_names
  # (added)
  out <- list(conv = conv, 
              coef = list(beta = betat, gamma = gammat),
              vcov = var_hat, loglik = ll_max, niter = niter)
  # (removed: grad_hat removed since it is not assigned)
  # out <- list(coef = list(beta = betat, gamma = gammat),
  #             vcov = var_hat, score = grad_hat,
  #             loglik = ll_max, niter)
  return(out)
}


#' @rdname hlm
#' @export
hlm_loglik <- function(beta, gamma, y, delta, X, W) {
  ## t <- exp(y)
  ## beta <- theta[1:px]
  ## gamma <- theta[px+1:pw]
  X <- as.matrix(X)
  W <- as.matrix(W)
  mu <- c(X%*%beta)
  sigma <- exp(c(W%*%gamma)/2)
  ll <- rep(NA, length(y))
  ll[delta] <- dnorm(y[delta], mean = mu[delta], sd = sigma[delta],
                     log = TRUE)
  ll[!delta] <- pnorm(y[!delta], mean = mu[!delta], sd = sigma[!delta],
                      lower.tail = FALSE, log.p = TRUE)
  ## ll    <- ifelse(delta, dnorm(y, mean = mu, sd = sigma, log = TRUE),
  ##                 pnorm(y, mean = mu, sd = sigma,
  ##                       lower.tail = FALSE, log.p = TRUE))
  return(sum(ll))
}
