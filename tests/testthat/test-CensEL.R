library(testthat)
library(flexEL)
context("CensEL")

source("flexEL-testfunctions.R")

#--- smooth_indicator ----------------------------------------------------------

test_that("R and C++ versions of `flexEL:::smooth_indicator()` are the same.", {
  job_descr <- expand.grid(has_inf = c("eps1", "eps2", "none"))
  n_jobs <- nrow(job_descr)
  for(ii in 1:n_jobs) {
    has_inf <- job_descr$has_inf[ii]
    n <- sample(10:20, 1)
    eps1 <- runif(1, -1, 1)
    eps2 <- runif(n, -1, 1)
    s <- runif(1, 0, 20)
    if(has_inf == "eps1") {
      eps1 <- -Inf
    } else if(has_inf == "eps2") {
      eps2[sample(n, 5)] <- -Inf
    }
    smooth_r <- smooth_indicator(eps1, eps2, s)
    smooth_cpp <- flexEL:::smooth_indicator(eps1, eps2, s)
    expect_equal(smooth_r, smooth_cpp)
  }
})

#--- expected_weights ----------------------------------------------------------

test_that("R and C++ versions of `CensEL$expected_weights()` are the same.", {
  n_test <- 50
  for (ii in 1:n_test) {
    n_obs <- sample(10:20,1)
    n_eqs <- sample(1:(n_obs-2), 1)
    supp_adj <- as.logical(rbinom(1, 1, .5))
    delta <- sample(0:1, size = n_obs, replace = TRUE)
    epsilon <- rnorm(n_obs)
    omega <- runif(n_obs+supp_adj, 0, 1)
    omega <- omega/sum(omega)
    smooth_s <- runif(1, 0, 1)
    gel <- GenEL$new(n_obs, n_eqs)
    cel <- CensEL$new(gel, smooth_s = smooth_s)
    weights_cpp <- cel$expected_weights(delta, epsilon, omega)
    weights_r <- expected_weights(delta, epsilon, omega, smooth_s)
    expect_equal(weights_cpp, weights_r)
  }
})


# ---- eval_weights ----

test_that("eval_weights with default settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    cel <- CensEL$new(n, p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    omega <- runif(n, 0, 1)
    omega <- omega/sum(omega)
    weights_cpp <- cel$eval_weights(delta, epsilon, omega)
    # weights_cpp
    weights_R <- evalWeights_R(delta, omega, epsilon)
    # weights_R
    # flexEL::.EvalWeights(omega, delta, epsilon, FALSE)
    expect_equal(weights_cpp, weights_R)
  }
})

test_that("eval_weights with support correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    cel <- CensEL$new(n, p)
    cel$supp_adj <- TRUE
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    omega <- runif(n + 1, 0, 1)
    omega <- omega/sum(omega)
    weights_cpp <- cel$eval_weights(delta, epsilon, omega)
    weights_R <- evalWeights_R(c(delta, 0), omega, c(epsilon, -Inf))
    # flexEL::.EvalWeights(omega, delta, epsilon, FALSE)
    expect_equal(weights_cpp, weights_R)
  }
})

test_that("eval_weights with continuity correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    cel <- CensEL$new(n, p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    omega <- runif(n, 0, 1)
    omega <- omega/sum(omega)
    s <- runif(1, 10, 100)
    cel$smooth <- TRUE
    cel$smooth_s <- s
    weights_cpp <- cel$eval_weights(delta, epsilon, omega)
    # weights_cpp
    weights_R <- evalWeights_smooth_R(delta, omega, epsilon, s)
    # weights_R
    # flexEL:::.EvalWeightsSmooth(omega, delta, epsilon, s, FALSE)
    expect_equal(weights_cpp, weights_R)
  }
})

test_that("eval_weights with support and continuity correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    cel <- CensEL$new(n, p)
    cel$supp_adj <- TRUE
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    omega <- runif(n + 1, 0, 1)
    omega <- omega/sum(omega)
    s <- runif(1, 10, 100)
    cel$smooth <- TRUE
    cel$smooth_s <- s
    weights_cpp <- cel$eval_weights(delta, epsilon, omega)
    # weights_cpp
    weights_R <- evalWeights_smooth_R(c(delta, 0), omega, c(epsilon, -Inf), s, TRUE)
    # weights_R
    # flexEL:::.EvalWeightsSmooth(omega, delta, epsilon, s, TRUE)
    expect_equal(weights_cpp, weights_R)
  }
})

# ---- omega_hat -----

# nconv <- 0
test_that("omega_hat with default settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    cel <- CensEL$new(n, p)
    omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    omega_cpp
    omega_R_lst <- omega_hat_EM_R(G, deltas = delta, epsilons = epsilon)
    omega_R_lst$omegas
    # range(omega_cpp-omega_R_lst$omegas)
    # omega_R_lst$conv
    if (!all(is.na(omega_cpp)) & omega_R_lst$conv) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R_lst$omegas, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("omega_hat with given convergence settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-7, 1e-5)
    abs_tol <- runif(1, 1e-4, 1e-3)
    cel <- CensEL$new(n, p)
    cel$max_iter_nr <- max_iter
    cel$rel_tol <- rel_tol
    cel$max_iter_em <- max_iter
    cel$abs_tol <- abs_tol
    omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    # omega_cpp
    omega_R_lst <- omega_hat_EM_R(G, deltas = delta, epsilons = epsilon,
                                  max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol)
    # omega_R_lst$omegas
    # range(omega_cpp-omega_R_lst$omegas)
    # omega_R_lst$conv
    if (!all(is.na(omega_cpp)) & omega_R_lst$conv) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R_lst$omegas, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("omega_hat with support correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    adj_a <- runif(1, 1, 5)
    cel <- CensEL$new(n, p)
    cel$supp_adj <- TRUE
    cel$supp_adj_a <- adj_a
    omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    # omega_cpp
    omega_R_lst <- omega_hat_EM_R(adjG_R(G, adj_a), deltas = delta, epsilons = epsilon, adjust = TRUE)
    # omega_R_lst$omegas
    # range(omega_cpp-omega_R_lst$omegas)
    # omega_R_lst$conv
    if (!all(is.na(omega_cpp)) & omega_R_lst$conv) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R_lst$omegas, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("omega_hat with continuity correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    s <- runif(1, 1, min(n, 100))
    cel <- CensEL$new(n, p)
    cel$smooth <- TRUE
    cel$smooth_s <- s
    omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    omega_cpp
    omega_R_lst <- omega_hat_EM_smooth_R(G, deltas = delta, epsilons = epsilon, s = s)
    omega_R_lst$omegas
    # range(omega_cpp-omega_R_lst$omegas)
    # omega_R_lst$conv
    if (!all(is.na(omega_cpp)) & omega_R_lst$conv) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R_lst$omegas, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("omega_hat with support and continuity correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    adj_a <- runif(1, 1, 5)
    s <- runif(1, 1, min(n, 100))
    cel <- CensEL$new(n, p)
    cel$supp_adj <- TRUE
    cel$supp_adj_a <- adj_a
    cel$smooth <- TRUE
    cel$smooth_s <- s
    omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    omega_cpp
    omega_R_lst <- omega_hat_EM_smooth_R(adjG_R(G, adj_a), deltas = delta,
                                         epsilons = epsilon, s = s, adjust = TRUE)
    omega_R_lst$omegas
    # range(omega_cpp-omega_R_lst$omegas)
    # omega_R_lst$conv
    if (!all(is.na(omega_cpp)) & omega_R_lst$conv) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R_lst$omegas, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("omega_hat with given convergence settings and support correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-7, 1e-5)
    abs_tol <- runif(1, 1e-4, 1e-3)
    adj_a <- runif(1, 1, 5)
    cel <- CensEL$new(n, p)
    cel$max_iter_nr <- max_iter
    cel$rel_tol <- rel_tol
    cel$max_iter_em <- max_iter
    cel$abs_tol <- abs_tol
    cel$supp_adj <- TRUE
    cel$supp_adj_a <- adj_a
    omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    # omega_cpp
    omega_R_lst <- omega_hat_EM_R(adjG_R(G, adj_a), deltas = delta, epsilons = epsilon, adjust = TRUE,
                                  max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol)
    # omega_R_lst$omegas
    range(omega_cpp-omega_R_lst$omegas)
    # omega_R_lst$conv
    if (!all(is.na(omega_cpp)) & omega_R_lst$conv) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R_lst$omegas, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("omega_hat with given convergence settings and continuity correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-7, 1e-5)
    abs_tol <- runif(1, 1e-4, 1e-3)
    s <- runif(1, 1, min(n, 100))
    cel <- CensEL$new(n, p)
    cel$max_iter_nr <- max_iter
    cel$rel_tol <- rel_tol
    cel$max_iter_em <- max_iter
    cel$abs_tol <- abs_tol
    cel$smooth <- TRUE
    cel$smooth_s <- s
    omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    # omega_cpp
    omega_R_lst <- omega_hat_EM_smooth_R(G, deltas = delta, epsilons = epsilon, s = s, adjust = FALSE,
                                         max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol)
    # omega_R_lst$omegas
    # range(omega_cpp-omega_R_lst$omegas)
    # omega_R_lst$conv
    if (!all(is.na(omega_cpp)) & omega_R_lst$conv) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R_lst$omegas, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("omega_hat with given convergence settings, support and continuity correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-7, 1e-5)
    abs_tol <- runif(1, 1e-4, 1e-3)
    s <- runif(1, 1, min(n, 100))
    adj_a <- runif(1, 1, 5)
    cel <- CensEL$new(n, p)
    cel$max_iter_nr <- max_iter
    cel$rel_tol <- rel_tol
    cel$max_iter_em <- max_iter
    cel$abs_tol <- abs_tol
    cel$smooth <- TRUE
    cel$smooth_s <- s
    cel$supp_adj <- TRUE
    cel$supp_adj_a <- adj_a
    omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    # omega_cpp
    omega_R_lst <- omega_hat_EM_smooth_R(adjG_R(G, adj_a), deltas = delta, epsilons = epsilon, s = s, adjust = TRUE,
                                         max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol)
    # omega_R_lst$omegas
    # range(omega_cpp-omega_R_lst$omegas)
    # omega_R_lst$conv
    if (!all(is.na(omega_cpp)) & omega_R_lst$conv) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R_lst$omegas, tolerance = 1e-3)
    }
  }
})
# nconv

# ---- logel ----

# converged result
check_res <- function(x) {
  all(is.finite(x) & !is.na(x))
}

# nconv <- 0
test_that("logel with default settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    cel <- CensEL$new(n, p)
    # cel$omega_hat(G, delta = delta, epsilon = epsilon)
    logel_cpp <- cel$logel(G, delta = delta, epsilon = epsilon)
    # logel_cpp
    omega_R_lst <- omega_hat_EM_R(G, deltas = delta, epsilons = epsilon)
    # omega_R_lst$omegas
    logel_R <- logEL_R(omega_R_lst$omegas, epsilon, delta)
    # logel_R
    # range(logel_cpp-logel_R)
    # omega_R_lst$conv
    if (check_res(logel_cpp) & check_res(logel_R)) {
      # nconv <<- nconv + 1
      expect_equal(logel_cpp, logel_R, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("logel with given convergence settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-7, 1e-5)
    abs_tol <- runif(1, 1e-4, 1e-3)
    cel <- CensEL$new(n, p)
    cel$set_opts(max_iter_nr = max_iter, rel_tol = rel_tol,
                 max_iter_em = max_iter, abs_tol = abs_tol)
    # cel$omega_hat(G, delta = delta, epsilon = epsilon)
    logel_cpp <- cel$logel(G, delta = delta, epsilon = epsilon)
    # logel_cpp
    omega_R_lst <- omega_hat_EM_R(G, deltas = delta, epsilons = epsilon,
                                  max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol)
    # omega_R_lst$omegas
    logel_R <- logEL_R(omega_R_lst$omegas, epsilon, delta)
    # logel_R
    # logel_cpp-logel_R
    # omega_R_lst$conv
    if (check_res(logel_cpp) & check_res(logel_R)) {
      # nconv <<- nconv + 1
      expect_equal(logel_cpp, logel_R, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("logel with support correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    adj_a <- runif(1, 1, 5)
    cel <- CensEL$new(n, p)
    cel$set_supp_adj(supp_adj = TRUE, supp_adj_a = adj_a)
    # omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    # omega_cpp
    logel_cpp <- cel$logel(G, delta, epsilon)
    omega_R_lst <- omega_hat_EM_R(adjG_R(G, adj_a), deltas = delta, epsilons = epsilon, adjust = TRUE)
    # omega_R_lst$omegas
    logel_R <- logEL_R(omega_R_lst$omegas, epsilon, delta, adjust = TRUE)
    # logel_cpp - logel_R
    if (check_res(logel_cpp) & check_res(logel_R)) {
      # nconv <<- nconv + 1
      expect_equal(logel_cpp, logel_R, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("logel with continuity correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    s <- runif(1, 1, min(n, 100))
    cel <- CensEL$new(n, p)
    cel$set_smooth(smooth = TRUE, smooth_s = s)
    # omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    # omega_cpp
    omega_R_lst <- omega_hat_EM_smooth_R(G, deltas = delta, epsilons = epsilon, s = s)
    # omega_R_lst$omegas
    # range(omega_cpp-omega_R_lst$omegas)
    logel_cpp <- cel$logel(G, delta, epsilon)
    logel_R <- logEL_smooth_R(omega_R_lst$omegas, epsilon, delta, s)
    logel_cpp-logel_R
    # omega_R_lst$conv
    if (check_res(logel_cpp) & check_res(logel_R)) {
      # nconv <<- nconv + 1
      expect_equal(logel_cpp, logel_R, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("logel with support and continuity correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    adj_a <- runif(1, 1, 5)
    s <- runif(1, 1, min(n, 100))
    cel <- CensEL$new(n, p)
    cel$set_supp_adj(supp_adj = TRUE, supp_adj_a = adj_a)
    cel$set_smooth(smooth = TRUE, smooth_s = s)
    # omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    # omega_cpp
    logel_cpp <- cel$logel(G, delta, epsilon)
    # logel_cpp
    omega_R_lst <- omega_hat_EM_smooth_R(adjG_R(G, adj_a), deltas = delta, epsilons = epsilon, s = s, adjust = TRUE)
    # omega_R_lst$omegas
    logel_R <- logEL_smooth_R(omega_R_lst$omegas, epsilon, delta, s = s, adjust = TRUE)
    # logel_R
    # logel_cpp - logel_R
    if (check_res(logel_cpp) & check_res(logel_R)) {
      # nconv <<- nconv + 1
      expect_equal(logel_cpp, logel_R, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("logel with given convergence settings and support correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-7, 1e-5)
    abs_tol <- runif(1, 1e-4, 1e-3)
    adj_a <- runif(1, 1, 5)
    cel <- CensEL$new(n, p)
    cel$set_opts(max_iter_nr = max_iter, rel_tol = rel_tol,
                 max_iter_em = max_iter, abs_tol = abs_tol)
    cel$set_supp_adj(supp_adj = TRUE, supp_adj_a = adj_a)
    # omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    logel_cpp <- cel$logel(G, delta = delta, epsilon = epsilon)
    # logel_cpp
    omega_R_lst <- omega_hat_EM_R(adjG_R(G, adj_a), deltas = delta, epsilons = epsilon,
                                  max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, adjust = TRUE)
    # omega_R_lst$omegas
    # range(omega_cpp - omega_R_lst$omegas)
    logel_R <- logEL_R(omega_R_lst$omegas, epsilon, delta, adjust = TRUE)
    # logel_R
    # logel_cpp-logel_R
    # omega_R_lst$conv
    if (check_res(logel_cpp) & check_res(logel_R)) {
      # nconv <<- nconv + 1
      expect_equal(logel_cpp, logel_R, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("logel with given convergence settings and continuity correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-7, 1e-5)
    abs_tol <- runif(1, 1e-4, 1e-3)
    s <- runif(1, 1, min(n, 100))
    cel <- CensEL$new(n, p)
    cel$set_opts(max_iter_nr = max_iter, rel_tol = rel_tol,
                 max_iter_em = max_iter, abs_tol = abs_tol)
    cel$set_smooth(smooth = TRUE, smooth_s = s)
    # omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    # omega_cpp
    omega_R_lst <- omega_hat_EM_smooth_R(G, deltas = delta, epsilons = epsilon, s = s,
                                         max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol)
    # omega_R_lst$omegas
    # range(omega_cpp-omega_R_lst$omegas)
    logel_cpp <- cel$logel(G, delta, epsilon)
    logel_R <- logEL_smooth_R(omega_R_lst$omegas, epsilon, delta, s)
    # logel_cpp-logel_R
    # omega_R_lst$conv
    if (check_res(logel_cpp) & check_res(logel_R)) {
      # nconv <<- nconv + 1
      expect_equal(logel_cpp, logel_R, tolerance = 1e-3)
    }
  }
})
# nconv

# nconv <- 0
test_that("logel with given convergence settings, support and continuity correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    G <- matrix(rnorm(n*p),n,p)
    delta <- rep(1,n)
    numcens <- sample(round(n/2),1)
    censinds <- sample(n,numcens)
    delta[censinds] <- 0
    epsilon <- rnorm(n)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-7, 1e-5)
    abs_tol <- runif(1, 1e-4, 1e-3)
    adj_a <- runif(1, 1, 5)
    s <- runif(1, 1, min(n, 100))
    cel <- CensEL$new(n, p)
    cel$set_opts(max_iter_nr = max_iter, rel_tol = rel_tol,
                 max_iter_em = max_iter, abs_tol = abs_tol)
    cel$set_supp_adj(supp_adj = TRUE, supp_adj_a = adj_a)
    cel$set_smooth(smooth = TRUE, smooth_s = s)
    # omega_cpp <- cel$omega_hat(G, delta = delta, epsilon = epsilon)
    # omega_cpp
    logel_cpp <- cel$logel(G, delta, epsilon)
    # logel_cpp
    omega_R_lst <- omega_hat_EM_smooth_R(adjG_R(G, adj_a),
                                         deltas = delta, epsilons = epsilon, s = s,
                                         adjust = TRUE, max_iter = max_iter,
                                         rel_tol = rel_tol, abs_tol = abs_tol)
    # omega_R_lst$omegas
    logel_R <- logEL_smooth_R(omega_R_lst$omegas, epsilon, delta, s = s, adjust = TRUE)
    # logel_R
    # logel_cpp - logel_R
    if (check_res(logel_cpp) & check_res(logel_R)) {
      # nconv <<- nconv + 1
      expect_equal(logel_cpp, logel_R, tolerance = 1e-3)
    }
  }
})
# nconv
