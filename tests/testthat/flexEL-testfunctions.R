#--- R version of C++ functions ------------------------------------------------

# setup for GenEL tests
GenEL_setup <- function(set_opts, supp_adj, weighted, check_conv,
                        n_obs, n_eqs, max_iter) {
  # setup
  if(missing(n_obs)) n_obs <- sample(10:20, 1)
  if(missing(n_eqs)) n_eqs <- sample(1:3, 1)
  gel <- GenEL$new(n_obs, n_eqs)
  G <- matrix(rnorm(n_obs * n_eqs), n_obs , n_eqs)
  weights <- runif(n_obs, 0, 2)
  if(missing(max_iter)) max_iter <- sample(c(1, 10, 50, 100), 1)
  rel_tol <- runif(1, 1e-6, 1e-2)
  supp_adj_a <- runif(1, 1, 5)
  weight_adj <- runif(1, 0, 2)
  lambda0 <- rnorm(n_eqs)
  # list of arguments to GenEL methods
  el_args <- list(G = G, check_conv = check_conv)
  # el_opts will contain all options to give lambda_nr.
  # if set_opts == TRUE, it will be passed to gel$set_opts
  el_opts <- list(max_iter = max_iter, rel_tol = rel_tol, lambda0 = lambda0,
                  supp_adj = supp_adj, supp_adj_a = supp_adj_a,
                  weight_adj = weight_adj)
  # options will be passed directly either way
  gel$max_iter <- max_iter
  gel$rel_tol <- rel_tol
  gel$lambda0 <- lambda0
  gel$supp_adj <- supp_adj
  gel$supp_adj_a <- supp_adj_a
  gel$weight_adj <- weight_adj
  if(weighted) {
    el_args <- c(el_args, list(weights = weights))
  }
  if(set_opts) {
    do.call(gel$set_opts, el_opts)
  }
  list(gel = gel, el_opts = el_opts, el_args = el_args)
}

# setup for CensEL tests
CensEL_setup <- function(set_opts, supp_adj, check_conv,
                         n_obs, n_eqs, max_iter_nr, max_iter_em) {
  gel_setup <- GenEL_setup(set_opts = TRUE, supp_adj = supp_adj,
                           weighted = FALSE, check_conv = check_conv,
                           n_obs = n_obs, n_eqs = n_eqs,
                           max_iter = max_iter_nr)
  n_obs <- gel_setup$gel$n_obs
  n_eqs <- gel_setup$gel$n_eqs
  supp_adj <- gel_setup$gel$supp_adj
  delta <- sample(0:1, size = n_obs, replace = TRUE)
  epsilon <- rnorm(n_obs)
  omega <- runif(n_obs+supp_adj, 0, 1)
  omega <- omega/sum(omega)
  smooth_s <- runif(1, 0, 1)
  if(missing(max_iter_em)) max_iter_em <- sample(c(100, 200, 500), 1)
  abs_tol <- runif(1, 1e-4, 1e-3)
  # list of arguments to CensEL methods
  cel_args <- list(G = gel_setup$el_args$G,
                   delta = delta, epsilon = epsilon,
                   omega = omega,
                   check_conv = check_conv)
  cel <- CensEL$new(n_obs, n_eqs)
  gel_opts <- gel_setup$el_opts
  cel_opts <- list(max_iter = max_iter_em, abs_tol = abs_tol,
                   smooth_s = smooth_s)
  # set options via active bindings
  cel$max_iter <- max_iter_em
  cel$abs_tol <- abs_tol
  cel$smooth_s <- smooth_s
  if(set_opts) {
    # set options via set_opts
    do.call(cel$set_opts, cel_opts)
  }
  # only way to set gel options
  do.call(cel$gel$set_opts, gel_opts)
  list(cel = cel,
       cel_opts = cel_opts, gel_opts = gel_opts,
       cel_args = cel_args)
}

# maximum relative error
max_rel_err <- function(lambda_new, lambda_old) {
  # Note: for stability (using the same tolerance as in C++ implementation)
  max(abs(lambda_new - lambda_old)/(.1 + abs(lambda_new + lambda_old)))
}

# support adjustment
adj_G <- function(G, supp_adj_a) {
  n_obs <- nrow(G)
  gbar <- colMeans(G)
  if (missing(supp_adj_a)) supp_adj_a <- max(1,0.5*log(n_obs))
  gadd <- -supp_adj_a*gbar
  return(unname(rbind(G,gadd)))
}

# weighted log_star and its derivatives
log_sharp <- function(x, q) {
  if (length(x) != length(q)) stop("x and q must be of the same length.")
  cond <- x >= q
  ans <- rep(NaN,length(x))
  ans[cond] <- log(x[cond])
  ans[!cond] <- -1/(2*q[!cond]^2)*x[!cond]^2 + 2/q[!cond]*x[!cond] - 3/2 + log(q[!cond])
  return(ans)
}

# 1st derivative of log_sharp
log_sharp1 <- function(x, q) {
  cond <- x >= q
  ans <- rep(NaN,length(x))
  ans[cond] <- 1/(x[cond])
  ans[!cond] <- -1/(q[!cond]^2)*x[!cond] + 2/q[!cond]
  return(ans)
}

# 2nd derivative of log_sharp
log_sharp2 <- function(x, q) {
  cond <- x >= q
  ans <- rep(NaN,length(x))
  ans[cond] <- -1/(x[cond]^2)
  ans[!cond] <- -1/q[!cond]^2
  return(ans)
}

# perform support adjustment on G, weights, and optionally delta and epsilon
# if supp_adj = FALSE just return inputs, except weights which gets converted
# from missing to ones.
adjust_support <- function(G, weights, delta, epsilon,
                           supp_adj, supp_adj_a, weight_adj) {
  n_obs <- nrow(G)
  n_eqs <- ncol(G)
  n_obs2 <- n_obs + supp_adj
  if(missing(weights)) {
    weights <- rep(1, n_obs)
  }
  if(supp_adj) {
    G <- adj_G(G, supp_adj_a)
    weights <- c(weights, weight_adj)
  }
  out <- list(G = G, weights = weights)
  if(!missing(delta)) {
    if(supp_adj) delta <- c(delta, FALSE)
    out <- c(out, list(delta = delta))
  }
  if(!missing(epsilon)) {
    if(supp_adj) epsilon <- c(epsilon, -Inf)
    out <- c(out, list(epsilon = epsilon))
  }
  out
}

# weighted newton-raphson
# returns: lambda and has_converged flag
lambda_nr <- function(G, weights, max_iter, rel_tol, lambda0,
                      supp_adj, supp_adj_a, weight_adj,
                      check_conv) {
  ## if(missing(weights)) {
  ##   weights <- rep(1, n_obs)
  ## }
  ## if(supp_adj) {
  ##   G <- adj_G(G, supp_adj_a)
  ##   weights <- c(weights, weight_adj)
  ## }
  ## G <- t(G)
  el_args <- adjust_support(G = G, weights = weights,
                            supp_adj = supp_adj, supp_adj_a = supp_adj_a,
                            weight_adj = weight_adj)
  weights <- el_args$weights
  G <- t(el_args$G)
  n_obs2 <- ncol(G)
  n_eqs <- nrow(G)
  lambda_old <- lambda0
  norm_weights <- weights/sum(weights)
  # newton-raphson loop
  for (ii in 1:max_iter) {
    # Q1 and Q2
    Glambda <- t(lambda_old) %*% G
    Glambda <- sum(norm_weights) - Glambda
    rho <- rep(NaN,n_obs2)
    Q2 <- matrix(rep(0,n_eqs*n_eqs), n_eqs, n_eqs)
    for (jj in 1:n_obs2) {
      rho[jj] <- log_sharp1(Glambda[jj], norm_weights[jj])
      Q2 <- Q2 - norm_weights[jj]*log_sharp2(Glambda[jj], norm_weights[jj])*(G[,jj] %*% t(G[,jj]))
    }
    Q1 <- G %*% (rho * norm_weights)
    lambda <- lambda_old - solve(Q2,Q1)
    err <- max_rel_err(lambda, lambda_old) # maximum relative error
    if (err < rel_tol) {
      break
    }
    lambda_old <- lambda # complete cycle
  }
  has_converged <- !(ii == max_iter && err > rel_tol)
  if(check_conv && !has_converged) {
    lambda <- rep(NaN, n_eqs)
  }
  list(lambda = c(lambda), has_converged = has_converged)
}

# weighted omega_hat for newton-raphson algorithm
omega_hat_nr <- function(G, weights, max_iter, rel_tol, lambda0,
                         supp_adj, supp_adj_a, weight_adj,
                         check_conv) {
  lambda_out <- lambda_nr(G = G, weights = weights,
                          max_iter = max_iter, rel_tol = rel_tol,
                          lambda0 = lambda0,
                          supp_adj = supp_adj, supp_adj_a = supp_adj_a,
                          weight_adj = weight_adj, check_conv = check_conv)
  el_args <- adjust_support(G = G, weights = weights,
                            supp_adj = supp_adj, supp_adj_a = supp_adj_a,
                            weight_adj = weight_adj)
  G <- el_args$G
  weights <- el_args$weights
  if (check_conv & !lambda_out$has_converged) {
    omega <- rep(NaN, nrow(G))
  }
  else {
    lambda <- lambda_out$lambda
    norm_weights <- weights / sum(weights)
    omega <- c(norm_weights/(1-lambda %*% t(G)))
  }
  # needs to be normalized for supp_adj
  omega / sum(omega)
}

# weighted logel and its gradient
logel_grad <- function(G, weights, max_iter, rel_tol, lambda0,
                      supp_adj, supp_adj_a, weight_adj,
                      check_conv) {
  lambda_out <- lambda_nr(G = G, weights = weights,
                          max_iter = max_iter, rel_tol = rel_tol,
                          lambda0 = lambda0,
                          supp_adj = supp_adj, supp_adj_a = supp_adj_a,
                          weight_adj = weight_adj, check_conv = check_conv)
  if(check_conv && !lambda_out$has_converged) {
    out <- list(logel = -Inf, grad = matrix(NaN, nrow(G), ncol(G)))
  } else {
    omega <- omega_hat_nr(G = G, weights = weights,
                          max_iter = max_iter, rel_tol = rel_tol,
                          lambda0 = lambda0,
                          supp_adj = supp_adj, supp_adj_a = supp_adj_a,
                          weight_adj = weight_adj, check_conv = check_conv)
    el_args <- adjust_support(G = G, weights = weights,
                              supp_adj = supp_adj, supp_adj_a = supp_adj_a,
                              weight_adj = weight_adj)
    weights <- el_args$weights
    lambda <- lambda_out$lambda
    logel <- if(all(omega > 0)) sum(weights * log(omega)) else NaN
    ## logel <- sum(weights * log(omega))
    n_obs <- nrow(G)
    grad <- t(lambda %*% t(omega[1:n_obs]))
    if(supp_adj) {
      grad <- grad - rep(omega[n_obs+1] * supp_adj_a/n_obs, n_obs) %*% t(lambda)
    }
    grad <- grad * sum(weights)
    out <- list(logel = logel, grad = grad)
  }
  out
}

# smoothed indicator function
smooth_indicator <- function(eps1, eps2, s) {
  1/(1 + exp(s * (eps1 - eps2)))
}

# expected weights from E-step of EM algorithm
expected_weights <- function(delta, epsilon, omega, smooth_s) {
  n_obs <- length(delta)
  n_obs2 <- length(omega)
  supp_adj <- n_obs2 == n_obs + 1
  if(supp_adj) {
    epsilon <- c(epsilon, -Inf)
    delta <- c(as.logical(delta), FALSE)
  }
  weights <- rep(0, n_obs2)
  for(ii in 1:n_obs2) {
    if(delta[ii]) {
      weights[ii] <- weights[ii] + 1
    } else {
      omega_tilde <- smooth_indicator(epsilon[ii], epsilon, smooth_s)
      omega_tilde[ii] <- .5
      omega_tilde <- omega_tilde * omega
      weights <-  weights + omega_tilde / sum(omega_tilde)
    }
  }
  ## weights <- rep(NA, n_obs2)
  ## for(ii in 1:n_obs2) {
  ##   omega_tilde <- smooth_indicator(epsilon[ii], epsilon, smooth_s)
  ##   omega_tilde[ii] <- .5
  ##   omega_tilde <- omega_tilde * omega
  ##   omega_tilde <- omega_tilde / sum(omega_tilde)
  ##   weights[ii] <- delta[ii] + sum((1-delta) * omega_tilde)
  ## }
  weights
}

# omega_hat from EM algorithm
# em_opts: list with elements smooth_s, max_iter, abs_tol
# nr_opts: list with elements max_iter, rel_tol, lambda0, supp_adj, supp_adj_a, weight_adj
# returns: omega_hat and has_converged
omega_hat_em <- function(G, delta, epsilon,
                         em_opts, nr_opts,
                         check_conv) {
  # setup
  cel <- CensEL$new(nrow(G), ncol(G))
  do.call(cel$set_opts, em_opts)
  do.call(cel$gel$set_opts, nr_opts)
  # initialize quantities
  n_obs <- cel$n_obs
  supp_adj <- cel$gel$supp_adj
  n_obs2 <- n_obs + supp_adj
  has_converged_em <- TRUE
  omega <- rep(1/n_obs2, n_obs2)
  weights <- cel$expected_weights(omega = omega,
                                  delta = delta, epsilon = epsilon)
  logel_old <- cel$gel$logel(G, check_conv = FALSE)
  for(ii in 1:cel$max_iter) {
    # M-step
    if(supp_adj) cel$gel$set_opts(weight_adj = weights[n_obs2])
    # lots of redundant calculations here
    lambda <- cel$gel$lambda_nr(G, weights = weights[1:n_obs],
                                check_conv = FALSE)
    omega <- cel$gel$omega_hat(G, weights = weights[1:n_obs],
                               check_conv = FALSE)
    logel <- cel$gel$logel(G,
                           weights = weights[1:n_obs],
                           check_conv = FALSE)
    # now establish convergence
    has_converged_nr <- cel$gel$lambda_nr(G, weights = weights[1:n_obs],
                                          check_conv = TRUE)
    has_converged_nr <- all(is.nan(has_converged_nr))
    em_err <- abs(logel - logel_old) # might be NaN if lambda_nr DNC
    if(has_converged_em && has_converged_nr) has_converged_em <- FALSE
    if(is.nan(em_err) || em_err < cel$abs_tol) break
    ## if(isTRUE(em_err < cel$abs_tol)) break
    logel_old <- logel
    cel$gel$set_opts(lambda0 = lambda) # update starting point
    # E-step
    weights <- cel$expected_weights(omega = omega,
                                    delta = delta, epsilon = epsilon)
  }
  has_converged_em <- has_converged_em &&
    !(ii == cel$max_iter && em_err > cel$abs_tol)
  if(check_conv && !has_converged_em) {
    omega <- rep(NaN, n_obs2)
  }
  list(omega = omega, has_converged = has_converged_em)
}

