# ---- testing R and C++ implementations of omega_hat are equal and optimality ----

source("test_utils.R")
source("el_rfuns.R")
source("reg_models.R")

# library(testthat) # not loaded automatically
context("omega_hat")

ntest <- 50

# Non-censored case:
# checking R and C++ implementations are equal
test_that("no censoring: omegahat_R == omegahat_cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p) # random G here
    omegahat_cpp <- omega_hat(G = G, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    omegahat_R <- omega_hat_R(G = G, max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    expect_equal(omegahat_cpp, omegahat_R)
  }
})

# checking optimality of the solution from C++
test_that("no censoring: omegahat_cpp is optimal", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(5, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p) # random G here
    omegahat_cpp <- omega_hat(G = G, max_iter = max_iter, rel_tol = rel_tol, 
                              support = FALSE, verbose = FALSE)
    # check optimality by optim_proj if converged
    if (!any(is.nan(omegahat_cpp))) {
      # optimal of omega_check should be at xsol = 1s if omegahat is optimal
      ocheck <- optim_proj(xsol = rep(1,n-p),
                           xrng = 0.05,
                           fun = function(x) {omega_check(x, omegahat_cpp, G)},
                           plot = FALSE)
      expect_lt(max_xdiff(ocheck),0.01)
    }
  }
})

# checking optimality of the solution from C++ (with support correction) 
test_that("no censoring: omegahat_R == omegahat_cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p),n,p) # random G here
    omegahat_cpp <- omega_hat(G = G, max_iter = max_iter, rel_tol = rel_tol, 
                              support = TRUE, verbose = FALSE)
    omegahat_R <- omega_hat_R(G = adjG_R(G), max_iter = max_iter, rel_tol = rel_tol, verbose = FALSE)
    expect_equal(omegahat_cpp, omegahat_R)
  }
})

# Censored case:
test_that("under censoring: omegahatC_R == omegahatC_cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    abs_tol <- runif(1, 1e-5, 1e-3)
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    omegahat_cpp <- omega_hat(G, deltas, epsilons, max_iter = max_iter, 
                              rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)
    omegahat_R <- omega_hat_R(G, deltas, epsilons, max_iter = max_iter, 
                              rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)$omegas
    if (!any(is.nan(omegahat_cpp)) && any(is.nan(omegahat_R))) {
      message("R version did not converge but C++ does.")
    }
    else {
      expect_equal(omegahat_cpp, omegahat_R)
    }
  }
})

test_that("under censoring: omegahatC_R == omegahatC_cpp (with support correction) ", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    abs_tol <- runif(1, 1e-5, 1e-3)
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    omegahat_cpp <- omega_hat(G, deltas, epsilons, max_iter = max_iter, 
                              rel_tol = rel_tol, abs_tol = abs_tol, support = TRUE, verbose = FALSE)
    omegahat_R <- omega_hat_EM_R(adjG_R(G), deltas, epsilons, adjust = TRUE, max_iter = max_iter, 
                              rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)$omegas
    if (!any(is.nan(omegahat_cpp)) && any(is.nan(omegahat_R))) {
      message("R version did not converge but C++ does.")
    }
    else {
      expect_equal(omegahat_cpp, omegahat_R, tolerance = 1e-4) # TODO: is this tolerance ok??
    }
  }
})

# checking optimality of the solution from C++
test_that("under censoring: omegahat_cpp is optimal", {
  for(ii in 1:ntest) {
    n <- sample(10:30,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    # max_iter <- 200
    rel_tol <- runif(1, 1e-8, 1e-6)
    abs_tol <- runif(1, 1e-5, 1e-3)
    # rel_tol <- 1e-7
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    omegahat_cpp <- omega_hat(G, deltas, epsilons, max_iter = max_iter, 
                              rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)
    if (!any(is.nan(omegahat_cpp))) {
      idx0 <- (abs(omegahat_cpp) < 1e-5 & !deltas)
      if (n-p-sum(idx0) > 0) {
        ocheck <- optim_proj(xsol = rep(1,n-p-sum(idx0)),
                             xrng = 0.01,
                             npts = 201, 
                             fun = function(x) {omega_pcheck(x, omegahat_cpp, G, deltas, epsilons, idx0, rel_tol)},
                             plot = FALSE)
        expect_lt(max_xdiff(ocheck), 0.01)
      }
    }
  }
})

# Censored case + smooth:
test_that("under censoring: omegahatCS_R == omegahatCS_cpp", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    s <- sample(1:100,1)
    support <- FALSE
    omegahat_cpp <- omega_hat_EM_smooth(G, deltas, epsilons, s, 
                                        max_iter = max_iter, rel_tol = rel_tol, abs_tol = rel_tol, 
                                        support = support, verbose = FALSE)
    omegahat_R <- omega_hat_EM_smooth_R(G, deltas, epsilons, s, max_iter = max_iter, rel_tol = rel_tol, abs_tol = rel_tol, verbose = FALSE)$omegas
    if (!any(is.nan(omegahat_cpp)) && any(is.nan(omegahat_R))) {
      message("R version did not converge but C++ does.")
    }
    else {
      expect_equal(omegahat_cpp, omegahat_R)
    }
  }
})

test_that("under censoring: omegahatCS_R == omegahatCS_cpp (with support correction)", {
  for(ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    epsilons <- rnorm(n)
    s <- sample(1:100,1)
    support <- TRUE
    omegahat_cpp <- omega_hat_EM_smooth(G, deltas, epsilons, s, 
                                        max_iter = max_iter, rel_tol = rel_tol, abs_tol = rel_tol, 
                                        support = support, verbose = FALSE)
    omegahat_R <- omega_hat_EM_smooth_R(adjG_R(G), deltas, epsilons, s, 
                                        max_iter = max_iter, rel_tol = rel_tol, abs_tol = rel_tol,
                                        adjust = support, verbose = FALSE)$omegas
    
    if (!any(is.nan(omegahat_cpp)) && any(is.nan(omegahat_R))) {
      message("R version did not converge but C++ does.")
    }
    else {
      expect_equal(omegahat_cpp, omegahat_R, tolerance = 1e-4)
    }
  }
})

# checking optimality of the solution from C++
test_that("under censoring: omegahat_cpp is optimal", {
  for(ii in 1:ntest) {
    n <- sample(10:30,1)
    p <- sample(1:(n-2), 1)
    # max_iter <- sample(c(2, 10, 100), 1)
    max_iter <- sample(c(10, 100, 500), 1)
    # max_iter <- 200
    rel_tol <- runif(1, 1e-8, 1e-6)
    # rel_tol <- 1e-7
    G <- matrix(rnorm(n*p), n, p)
    deltas <- rep(1,n)
    numcens <- sample(round(n/4),1)
    censinds <- sample(n,numcens)
    deltas[censinds] <- 0
    1-sum(deltas)/n
    epsilons <- rnorm(n)
    s <- sample(1:100,1)
    support <- FALSE
    omegahat_cpp <- omega_hat_EM_smooth(G, deltas, epsilons, s, max_iter = max_iter, 
                                          rel_tol = 1e-5, abs_tol = 1e-5, support = support, verbose = FALSE)
    # omegahat_cpp <- omega_hat_EM_smooth_R(G, deltas, epsilons, s, max_iter = max_iter,
    #                                       rel_tol = 1e-3, verbose = FALSE)$omegas
    if (!any(is.nan(omegahat_cpp))) {
      idx0 <- (abs(omegahat_cpp) < 1e-5 & !deltas)
      if (n-p-sum(idx0) > 0) {
        ocheck <- optim_proj(xsol = rep(1,n-p-sum(idx0)),
                             xrng = 0.01,
                             npts = 201,
                             fun = function(x) {omega_smooth_pcheck(x, omegahat_cpp, G, deltas, epsilons, idx0, s)},
                             plot = FALSE)
        # ocheck <- optim_proj(xsol = rep(1,n-p),
        #                      xrng = 0.001,
        #                      npts = 201, 
        #                      fun = function(x) {omega_smooth_check(x, omegahat_cpp, G, deltas, epsilons, s)},
        #                      plot = FALSE)
        # print(ocheck)
        expect_lt(max_xdiff(ocheck), 0.01)
      }
    }
  }
})

