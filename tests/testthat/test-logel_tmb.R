## library(testthat)
## library(flexEL)

source("flexEL-testfunctions.R")

test_that("R and TMB implementation of `logel()` are the same.", {
  skip_if_not_installed(pkg = "TMB")
  # compile TMB file
  model <- "tmb_el"
  flexel_path <- system.file("include", package = "flexEL")
  TMB::compile(paste0(model, ".cpp"),
               PKG_CXXFLAGS = paste0("-std=gnu++11 -I",
                                     flexel_path))
  dyn.load(TMB::dynlib(model))
  # test itself
  test_cases <- expand.grid(
    set_opts = c(FALSE, TRUE),
    supp_adj = c(FALSE, TRUE),
    weighted = c(FALSE),
    check_conv = c(TRUE),
    stringsAsFactors = FALSE
  )
  n_test <- nrow(test_cases)
  for(ii in 1:n_test) {
    # setup
    setup <- GenEL_setup(set_opts = test_cases$set_opts[ii],
                         supp_adj = test_cases$supp_adj[ii],
                         weighted = test_cases$weighted[ii],
                         check_conv = test_cases$check_conv[ii])
    setup$el_args <- setup$el_args["G"]
    setup$el_opts <- c(setup$el_opts[c("max_iter", "rel_tol", "supp_adj")])
    out_r <- do.call(flexEL::logEL,
                     c(setup$el_args, setup$el_opts, list(grad = TRUE)))
    out_tmb <- do.call(logel_tmb,
                       c(list(model = model), setup$el_args, setup$el_opts))
    expect_equal(out_tmb, out_r)
  }
})
