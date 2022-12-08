## library(testthat)
## library(flexEL)
context("GenEL")

source("flexEL-testfunctions.R")

test_that("R and C++ implementation of `GenEL$lambda_nr()` are the same.", {
  test_cases <- expand.grid(
    set_opts = c(FALSE, TRUE),
    supp_adj = c(FALSE, TRUE),
    weighted = c(FALSE, TRUE)
  )
  n_test <- nrow(test_cases)
  for(ii in 1:n_test) {
    # setup
    check_conv <- sample(c(TRUE, FALSE), 1)
    setup <- GenEL_setup(set_opts = test_cases$set_opts[ii],
                         supp_adj = test_cases$supp_adj[ii],
                         weighted = test_cases$weighted[ii],
                         check_conv = check_conv)
    out_cpp <- do.call(setup$gel$lambda_nr, setup$el_args)
    out_r <- do.call(lambda_nr, c(setup$el_args, setup$el_opts))$lambda
    expect_equal(out_cpp, out_r)
  }
})

test_that("R and C++ implementation of `GenEL$omega_hat()` are the same.", {
  test_cases <- expand.grid(
    set_opts = c(FALSE, TRUE),
    supp_adj = c(FALSE, TRUE),
    weighted = c(FALSE, TRUE)
  )
  n_test <- nrow(test_cases)
  for(ii in 1:n_test) {
    # setup
    check_conv <- sample(c(TRUE, FALSE), 1)
    setup <- GenEL_setup(set_opts = test_cases$set_opts[ii],
                         supp_adj = test_cases$supp_adj[ii],
                         weighted = test_cases$weighted[ii],
                         check_conv = check_conv)
    out_cpp <- do.call(setup$gel$omega_hat, setup$el_args)
    out_r <- do.call(omega_hat, c(setup$el_args, setup$el_opts))
    expect_equal(out_cpp, out_r)
  }
})

test_that("R and C++ implementation of `GenEL$logel()` are the same.", {
  test_cases <- expand.grid(
    set_opts = c(FALSE, TRUE),
    supp_adj = c(FALSE, TRUE),
    weighted = c(FALSE, TRUE)
  )
  n_test <- nrow(test_cases)
  for(ii in 1:n_test) {
    # setup
    check_conv <- sample(c(TRUE, FALSE), 1)
    setup <- GenEL_setup(set_opts = test_cases$set_opts[ii],
                         supp_adj = test_cases$supp_adj[ii],
                         weighted = test_cases$weighted[ii],
                         check_conv = check_conv)
    out_cpp <- do.call(setup$gel$logel, setup$el_args)
    out_r <- do.call(logel_grad, c(setup$el_args, setup$el_opts))$logel
    expect_equal(out_cpp, out_r)
  }
})

test_that("R and C++ implementation of `GenEL$logel_grad()` are the same.", {
  # TODO: check against numerical gradient of `GenEL$logel()`
  test_cases <- expand.grid(
    set_opts = c(FALSE, TRUE),
    supp_adj = c(FALSE, TRUE),
    weighted = c(FALSE, TRUE)
  )
  n_test <- nrow(test_cases)
  for(ii in 1:n_test) {
    # setup
    check_conv <- sample(c(TRUE, FALSE), 1)
    setup <- GenEL_setup(set_opts = test_cases$set_opts[ii],
                         supp_adj = test_cases$supp_adj[ii],
                         weighted = test_cases$weighted[ii],
                         check_conv = check_conv)
    out_cpp <- do.call(setup$gel$logel_grad, setup$el_args)
    out_r <- do.call(logel_grad, c(setup$el_args, setup$el_opts))
    expect_equal(out_cpp, out_r)
  }
})
