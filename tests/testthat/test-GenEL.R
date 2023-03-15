## library(testthat)
## library(flexEL)

source("flexEL-testfunctions.R")

# TODO: add test for support adjustment

test_that("R and C++ implementation of `GenEL$lambda_nr()` are the same.", {
  test_cases <- expand.grid(
    set_opts = c(FALSE, TRUE),
    supp_adj = c(FALSE, TRUE),
    weighted = c(FALSE, TRUE),
    check_conv = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  n_test <- nrow(test_cases)
  for(ii in 1:n_test) {
    # setup
    setup <- GenEL_setup(set_opts = test_cases$set_opts[ii],
                         supp_adj = test_cases$supp_adj[ii],
                         weighted = test_cases$weighted[ii],
                         check_conv = test_cases$check_conv[ii])
    out_cpp <- do.call(setup$gel$lambda_nr, setup$el_args)
    out_r <- do.call(lambda_nr, c(setup$el_args, setup$el_opts))$lambda
    expect_equal(out_cpp, out_r)
  }
})

test_that("R and C++ implementation of `GenEL$omega_hat()` are the same.", {
  test_cases <- expand.grid(
    set_opts = c(FALSE, TRUE),
    supp_adj = c(FALSE, TRUE),
    weighted = c(FALSE, TRUE),
    check_conv = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  n_test <- nrow(test_cases)
  for(ii in 1:n_test) {
    # setup
    setup <- GenEL_setup(set_opts = test_cases$set_opts[ii],
                         supp_adj = test_cases$supp_adj[ii],
                         weighted = test_cases$weighted[ii],
                         check_conv = test_cases$check_conv[ii])
    out_cpp <- do.call(setup$gel$omega_hat, setup$el_args)
    out_r <- do.call(omega_hat_nr, c(setup$el_args, setup$el_opts))
    expect_equal(out_cpp, out_r)
  }
})

test_that("R and C++ implementation of `GenEL$logel()` are the same.", {
  test_cases <- expand.grid(
    set_opts = c(FALSE, TRUE),
    supp_adj = c(FALSE, TRUE),
    weighted = c(FALSE, TRUE),
    check_conv = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  n_test <- nrow(test_cases)
  for(ii in 1:n_test) {
    # setup
    setup <- GenEL_setup(set_opts = test_cases$set_opts[ii],
                         supp_adj = test_cases$supp_adj[ii],
                         weighted = test_cases$weighted[ii],
                         check_conv = test_cases$check_conv[ii])
    out_cpp <- do.call(setup$gel$logel, setup$el_args)
    out_r <- do.call(logel_grad, c(setup$el_args, setup$el_opts))$logel
    expect_equal(out_cpp, out_r)
  }
})

test_that("R and C++ implementation of `GenEL$logel(full_out = TRUE)` are the same.", {
  test_cases <- expand.grid(
    set_opts = c(FALSE, TRUE),
    supp_adj = c(FALSE, TRUE),
    weighted = c(FALSE, TRUE),
    check_conv = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  n_test <- nrow(test_cases)
  for(ii in 1:n_test) {
    # setup
    setup <- GenEL_setup(set_opts = test_cases$set_opts[ii],
                         supp_adj = test_cases$supp_adj[ii],
                         weighted = test_cases$weighted[ii],
                         check_conv = test_cases$check_conv[ii],
                         full_out = TRUE)
    out_cpp <- do.call(setup$gel$logel, setup$el_args)
    out_r <- do.call(logel_grad, c(setup$el_args, setup$el_opts))
    out_r <- out_r[c("logel", "omega", "lambda")]
    expect_equal(out_cpp, out_r)
  }
})

test_that("R and C++ implementation of `GenEL$logel_grad()` are the same.", {
  # TODO: check against numerical gradient of `GenEL$logel()`
  test_cases <- expand.grid(
    set_opts = c(FALSE, TRUE),
    supp_adj = c(FALSE, TRUE),
    weighted = c(FALSE, TRUE),
    check_conv = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  n_test <- nrow(test_cases)
  for(ii in 1:n_test) {
    # setup
    setup <- GenEL_setup(set_opts = test_cases$set_opts[ii],
                         supp_adj = test_cases$supp_adj[ii],
                         weighted = test_cases$weighted[ii],
                         check_conv = test_cases$check_conv[ii])
    out_cpp <- do.call(setup$gel$logel_grad, setup$el_args)
    out_r <- do.call(logel_grad, c(setup$el_args, setup$el_opts))
    out_r <- out_r[c("logel", "grad")]
    ## if(any(sapply(out_r, anyNA))) browser()
    expect_equal(out_cpp, out_r)
  }
})
