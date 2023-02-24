## library(testthat)
## library(flexEL)

source("flexEL-testfunctions.R")

# TODO: use CensEL_setup() in CensEL$expected_weights() test.

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

#--- CensEL methods ------------------------------------------------------------

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
    cel <- CensEL$new(n_obs, n_eqs, smooth_s = smooth_s,
                      gel_opts = list(supp_adj = supp_adj))
    weights_cpp <- cel$expected_weights(delta, epsilon, omega)
    weights_r <- expected_weights(delta, epsilon, omega, smooth_s)
    expect_equal(weights_cpp, weights_r)
  }
})

test_that("R and C++ versions of `CensEL$omega_hat()` are the same.", {
  test_cases <- expand.grid(
    set_opts = c(FALSE, TRUE),
    supp_adj = c(FALSE, TRUE),
    check_conv = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  n_test <- nrow(test_cases)
  for (ii in 1:n_test) {
    # setup
    ## ii <- 1
    ## check_conv <- FALSE
    ## check_conv <- sample(c(TRUE, FALSE), 1)
    setup <- CensEL_setup(set_opts = test_cases$set_opts[ii],
                          supp_adj = test_cases$supp_adj[ii],
                          check_conv = test_cases$check_conv[ii])
    arg_names <- c("G", "delta", "epsilon", "check_conv")
    out_cpp <- do.call(setup$cel$omega_hat,
                       args = setup$cel_args[arg_names])
    out_r <- do.call(omega_hat_em,
                     c(setup$cel_args[arg_names],
                       list(em_opts = setup$cel_opts,
                            nr_opts = setup$gel_opts)))
    expect_equal(out_cpp, out_r$omega)
  }
})

