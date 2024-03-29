# library(testthat)
# library(flexEL)
context("GenEL")

source("el_rfuns.R")

test_cases <- expand.grid(
  set_opts = c(FALSE, TRUE),
  supp_adj = c(FALSE, TRUE),
  weights = c(FALSE, TRUE),
  method = c("lambda_nr", "omega_hat", "logel", "logel_grad")
)
ntest <- nrow(test_cases)

test_that("R and C++ implementations of methods are the same.", {
  for(ii in 1:ntest) {
    # setup
    n <- sample(20:50,1)
    p <- sample(1:3, 1)
    gel <- GenEL$new(n, p)
    G <- matrix(rnorm(n*p),n,p)
    weights <- runif(n, 0, 2)
    # construct argument lists from test cases
    el_opts <- NULL # list of arguments to set_opts
    el_args <- list(G = G) # list of arguments to other methods
    if(test_cases$set_opts[ii]) {
      # set nr parameters
      max_iter <- sample(c(1, 10, 50, 100), 1)
      rel_tol <- runif(1, 1e-6, 1e-2)
      el_opts <- c(el_opts, list(max_iter = max_iter, rel_tol = rel_tol))
      # check active members
      gel$max_iter <- max_iter
      gel$rel_tol <- rel_tol
    }
    if(test_cases$supp_adj[ii]) {
      # set support correction
      adj_a <- runif(1, 1, 5)
      el_opts <- c(el_opts, list(supp_adj = TRUE, supp_adj_a = adj_a))
    }
    if(test_cases$weights[ii]) {
      # weighted version
      el_args <- c(el_args, list(weights = weights))
    }
    # R and C++ calculations
    if(length(el_opts) > 0) {
      do.call(gel$set_opts, args = el_opts)
    }
    if(test_cases$method[ii] == "lambda_nr") {
      lambda_cpp <- do.call(gel$lambda_nr, el_args)
      lambda_R_lst <- do.call(lambda_nr_R, args = c(el_opts, el_args))
      if (!lambda_R_lst$convergence) {
        expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
      } else {
        expect_equal(lambda_cpp, lambda_R_lst$lambda)
      }
    }
    if(test_cases$method[ii] == "omega_hat") {
      omega_cpp <- do.call(gel$omega_hat, el_args)
      omega_R <- do.call(omega_hat_R, args = c(el_opts, el_args))
      if (!any(is.na(omega_cpp)) & !any(is.na(omega_R))) {
        expect_equal(omega_cpp, omega_R, tolerance = 1e-6)
      } else {
        expect_equal(any(is.na(omega_cpp)), any(is.na(omega_R)))
      }
    }
    if(test_cases$method[ii] == "logel") {
      lel_cpp <- do.call(gel$logel, el_args)
      lel_R <- do.call(logel_R, args = c(el_opts, el_args))
      expect_equal(lel_cpp, lel_R)
    }
    if(test_cases$method[ii] == "logel_grad") {
      dldG_cpp <- do.call(gel$logel_grad, el_args)$grad
      dldG_R <- do.call(logel_grad_R, args = c(el_opts, el_args))
      expect_equal(dldG_cpp, dldG_R, tolerance = 1e-4)
    }
  }
})

#--- old tests -----------------------------------------------------------------

## ntest <- 3

## # nconv <- 0
## test_that("lambda_nr with default options", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     gel <- GenEL$new(n, p)
##     G <- matrix(rnorm(n*p),n,p)
##     lambda_cpp <- gel$lambda_nr(G)
##     # lambda_cpp
##     lambda_R_lst <- lambdaNR_R(G)
##     # lambda_R_lst$lambda
##     if (!lambda_R_lst$convergence) {
##       expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
##     }
##     else {
##       # nconv <<- nconv + 1
##       expect_equal(lambda_cpp, lambda_R_lst$lambda)
##     }
##   }
## })

## # nconv <- 0
## test_that("lambda_nr with given convergence settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(1, 10, 50, 100), 1)
##     rel_tol <- runif(1, 1e-6, 1e-2)
##     gel <- GenEL$new(n, p)
##     gel$max_iter <- max_iter
##     gel$rel_tol <- rel_tol
##     G <- matrix(rnorm(n*p),n,p)
##     lambda_cpp <- gel$lambda_nr(G)
##     lambda_R_lst <- lambdaNR_R(G, max_iter = max_iter, rel_tol = rel_tol)
##     if (!lambda_R_lst$convergence) {
##       expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
##     }
##     else {
##       # nconv <<- nconv + 1
##       expect_equal(lambda_cpp, lambda_R_lst$lambda)
##     }
##   }
## })

## # nconv <- 0
## test_that("lambda_nr with given convergence settings and support correction", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(1, 10, 50, 100), 1)
##     rel_tol <- runif(1, 1e-6, 1e-2)
##     adj_a <- runif(1, 1, 5)
##     gel <- GenEL$new(n, p)
##     gel$set_opts(max_iter = max_iter, rel_tol = rel_tol,
##                  supp_adj = TRUE, supp_adj_a = adj_a)
##     G <- matrix(rnorm(n*p),n,p)
##     lambda_cpp <- gel$lambda_nr(G)
##     lambda_cpp
##     lambda_R_lst <- lambdaNR_R(adjG_R(G, adj_a), max_iter = max_iter, rel_tol = rel_tol)
##     lambda_R_lst$lambda
##     if (!lambda_R_lst$convergence) {
##       expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
##     }
##     else {
##       # nconv <<- nconv + 1
##       expect_equal(lambda_cpp, lambda_R_lst$lambda)
##     }
##   }
## })

## # ---- lambda_nr with weights ----

## # nconv <- 0
## test_that("weighted lambda_nr with default options", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     gel <- GenEL$new(n, p)
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     lambda_cpp <- gel$lambda_nr(G, weights)
##     # lambda_cpp
##     lambda_R_lst <- lambdaNRC_R(G, weights/sum(weights))
##     # lambda_R_lst$lambda
##     if (!lambda_R_lst$convergence) {
##       expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
##     }
##     else {
##       # nconv <<- nconv + 1
##       expect_equal(lambda_cpp, lambda_R_lst$lambda)
##     }
##   }
## })

## # nconv <- 0
## test_that("weighted lambda_nr with given convergence settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(1, 10, 50, 100), 1)
##     rel_tol <- runif(1, 1e-6, 1e-4)
##     gel <- GenEL$new(n, p)
##     gel$max_iter <- max_iter
##     gel$rel_tol <- rel_tol
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     lambda_cpp <- gel$lambda_nr(G, weights)
##     # lambda_cpp
##     lambda_R_lst <- lambdaNRC_R(G, weights/sum(weights),
##                                 max_iter = max_iter, rel_tol = rel_tol)
##     # lambda_R_lst$lambda
##     if (!lambda_R_lst$convergence) {
##       expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
##     }
##     else {
##       # nconv <<- nconv + 1
##       expect_equal(lambda_cpp, lambda_R_lst$lambda)
##     }
##   }
## })

## # nconv <- 0
## test_that("weighted lambda_nr with given convergence settings and support correction", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(1, 10, 50, 100), 1)
##     rel_tol <- runif(1, 1e-6, 1e-2)
##     adj_a <- runif(1, 1, 5)
##     gel <- GenEL$new(n, p)
##     gel$set_opts(max_iter = max_iter, rel_tol = rel_tol,
##                  supp_adj = TRUE, supp_adj_a = adj_a)
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     lambda_cpp <- gel$lambda_nr(G, weights)
##     # lambda_cpp
##     weights_R <- c(weights, 1)
##     lambda_R_lst <- lambdaNRC_R(adjG_R(G, adj_a), weights_R/sum(weights_R),
##                                 max_iter = max_iter, rel_tol = rel_tol)
##     # lambda_R_lst$lambda
##     if (!lambda_R_lst$convergence) {
##       expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
##     }
##     else {
##       # nconv <<- nconv + 1
##       expect_equal(lambda_cpp, lambda_R_lst$lambda)
##     }
##   }
## })

## # ---- omega_hat ----

## # nconv <- 0
## test_that("omega_hat with default settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     gel <- GenEL$new(n, p)
##     G <- matrix(rnorm(n*p),n,p)
##     omega_cpp <- gel$omega_hat(G)
##     omega_R <- omega_hat_R(G)
##     if (!any(is.na(omega_cpp)) & !any(is.na(omega_R))) {
##       # nconv <<- nconv + 1
##       expect_equal(omega_cpp, omega_R)
##     }
##     else {
##       expect_equal(any(is.na(omega_cpp)), any(is.na(omega_R)))
##     }
##   }
## })

## # nconv <- 0
## test_that("omega_hat with given convergence settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(1, 10, 50, 100), 1)
##     rel_tol <- runif(1, 1e-6, 1e-2)
##     gel <- GenEL$new(n, p)
##     gel$max_iter <- max_iter
##     gel$rel_tol <- rel_tol
##     G <- matrix(rnorm(n*p),n,p)
##     omega_cpp <- gel$omega_hat(G)
##     omega_R <- omega_hat_R(G, max_iter = max_iter, rel_tol = rel_tol)
##     if (!any(is.na(omega_cpp)) & !any(is.na(omega_R))) {
##       # nconv <<- nconv + 1
##       expect_equal(omega_cpp, omega_R)
##     }
##     else {
##       expect_equal(any(is.na(omega_cpp)), any(is.na(omega_R)))
##     }
##   }
## })

## # nconv <- 0
## test_that("omega_hat with default settings and support correction", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(1, 10, 50, 100), 1)
##     rel_tol <- runif(1, 1e-6, 1e-2)
##     adj_a <- runif(1, 1, 5)
##     gel <- GenEL$new(n, p)
##     gel$set_opts(max_iter = max_iter, rel_tol = rel_tol,
##                  supp_adj = TRUE, supp_adj_a = adj_a)
##     G <- matrix(rnorm(n*p),n,p)
##     omega_cpp <- gel$omega_hat(G)
##     omega_R <- omega_hat_R(adjG_R(G, adj_a), adjust = TRUE,
##                            max_iter = max_iter, rel_tol = rel_tol)
##     if (!any(is.na(omega_cpp)) & !any(is.na(omega_R))) {
##       # nconv <<- nconv + 1
##       expect_equal(omega_cpp, omega_R)
##     }
##     else {
##       expect_equal(any(is.na(omega_cpp)), any(is.na(omega_R)))
##     }
##   }
## })

## # ---- omega_hat with weights ----

## # nconv <- 0
## test_that("weighted omega_hat with default settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     gel <- GenEL$new(n, p)
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     weights <- weights/sum(weights)
##     omega_cpp <- gel$omega_hat(G, weights)
##     # omega_cpp
##     omega_R <- weighted_omega_hat_R(G, weights)
##     # omega_R
##     if (!any(is.na(omega_cpp)) & !any(is.na(omega_R))) {
##       # nconv <<- nconv + 1
##       expect_equal(omega_cpp, omega_R)
##     }
##     else {
##       expect_equal(any(is.na(omega_cpp)), any(is.na(omega_R)))
##     }
##   }
## })

## # nconv <- 0
## test_that("weighted omega_hat with given convergence settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(1, 10, 50, 100), 1)
##     rel_tol <- runif(1, 1e-6, 1e-2)
##     gel <- GenEL$new(n, p)
##     gel$max_iter <- max_iter
##     gel$rel_tol <- rel_tol
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     omega_cpp <- gel$omega_hat(G, weights)
##     # omega_cpp
##     omega_R <- weighted_omega_hat_R(G, weights/sum(weights),
##                                     max_iter = max_iter, rel_tol = rel_tol)
##     # omega_R
##     if (!any(is.na(omega_cpp)) & !any(is.na(omega_R))) {
##       # nconv <<- nconv + 1
##       expect_equal(omega_cpp, omega_R, tolerance = 1e-4)
##     }
##     else {
##       expect_equal(any(is.na(omega_cpp)), any(is.na(omega_R)))
##     }
##   }
## })

## # nconv <- 0
## test_that("weighted omega_hat with default settings and support correction", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(1, 10, 50, 100), 1)
##     rel_tol <- runif(1, 1e-6, 1e-2)
##     adj_a <- runif(1, 1, 5)
##     gel <- GenEL$new(n, p)
##     gel$set_opts(max_iter = max_iter, rel_tol = rel_tol,
##                  supp_adj = TRUE, supp_adj_a = adj_a)
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     omega_cpp <- gel$omega_hat(G, weights)
##     # omega_cpp
##     weights_R <- c(weights, 1)
##     omega_R <- weighted_omega_hat_R(adjG_R(G, adj_a), weights_R/sum(weights_R),
##                                     adjust = TRUE,
##                                     max_iter = max_iter, rel_tol = rel_tol)
##     # omega_R
##     if (!any(is.na(omega_cpp)) & !any(is.na(omega_R))) {
##       # nconv <<- nconv + 1
##       expect_equal(omega_cpp, omega_R, tolerance = 1e-4)
##     }
##     else {
##       expect_equal(any(is.na(omega_cpp)), any(is.na(omega_R)))
##     }
##   }
## })

## # ---- logel ----

## # converged result
## check_res <- function(x) {
##   all(is.finite(x) & !is.na(x))
## }

## # nconv <- 0
## test_that("logel with default settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     gel <- GenEL$new(n, p)
##     G <- matrix(rnorm(n*p),n,p)
##     logel_cpp <- gel$logel(G)
##     omegahat_R <- omega_hat_R(G)
##     logel_R <- logEL_R(omegahat_R)
##     ## if (check_res(logel_R) & check_res(logel_cpp)) {
##       # nconv <<- nconv + 1
##       expect_equal(logel_cpp, logel_R)
##     ## }
##   }
## })

## # nconv <- 0
## test_that("logel with given convergence settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(1, 10, 50, 100), 1)
##     rel_tol <- runif(1, 1e-6, 1e-2)
##     gel <- GenEL$new(n, p)
##     gel$max_iter <- max_iter
##     gel$rel_tol <- rel_tol
##     G <- matrix(rnorm(n*p),n,p)
##     logel_cpp <- gel$logel(G)
##     # logel_cpp
##     omegahat_R <- omega_hat_R(G, max_iter = max_iter, rel_tol = rel_tol)
##     logel_R <- logEL_R(omegahat_R)
##     # logel_R
##     ## if (check_res(logel_R) & check_res(logel_cpp)) {
##       # nconv <<- nconv + 1
##       expect_equal(logel_cpp, logel_R)
##     ## }
##   }
## })

## # nconv <- 0
## test_that("logel with given convergence settings and support correction", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(100, 200, 500), 1)
##     rel_tol <- runif(1, 1e-6, 1e-2)
##     adj_a <- runif(1, 1, 5)
##     gel <- GenEL$new(n, p)
##     gel$set_opts(max_iter = max_iter, rel_tol = rel_tol,
##                  supp_adj = TRUE, supp_adj_a = adj_a)
##     G <- matrix(rnorm(n*p),n,p)
##     logel_cpp <- gel$logel(G)
##     omegahat_R <- omega_hat_R(adjG_R(G, adj_a), max_iter = max_iter, rel_tol = rel_tol)
##     logel_R <- logEL_R(omegahat_R, adjust = TRUE)
##     ## if (check_res(logel_R) & check_res(logel_cpp)) {
##       # nconv <<- nconv + 1
##       expect_equal(logel_cpp, logel_R)
##     ## }
##   }
## })

## # ---- logel with weights ----

## # nconv <- 0
## test_that("weighted logel with default settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     gel <- GenEL$new(n, p)
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     # weights <- weights/sum(weights)
##     logel_cpp <- gel$logel(G, weights)
##     # logel_cpp
##     logel_R <- weighted_logEL_G_R(G, weights)
##     # as.numeric(logel_R)
##     if (check_res(logel_R) & check_res(logel_cpp)) {
##       # nconv <<- nconv + 1
##       expect_equal(logel_cpp, as.numeric(logel_R))
##     }
##   }
## })

## # nconv <- 0
## test_that("weighted logel with given convergence settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(100, 200, 500), 1)
##     rel_tol <- runif(1, 1e-6, 1e-4)
##     gel <- GenEL$new(n, p)
##     gel$max_iter <- max_iter
##     gel$rel_tol <- rel_tol
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     # weights <- weights/sum(weights)
##     logel_cpp <- gel$logel(G, weights)
##     # logel_cpp
##     logel_R <- weighted_logEL_G_R(G, weights, max_iter, rel_tol)
##     # as.numeric(logel_R)
##     if (check_res(logel_R) & check_res(logel_cpp)) {
##       # nconv <<- nconv + 1
##       expect_equal(logel_cpp, as.numeric(logel_R), tolerance = 1e-4)
##     }
##   }
## })

## # nconv <- 0
## test_that("weighted logel with given convergence settings and support correction", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:(n-2), 1)
##     max_iter <- sample(c(100, 200, 500), 1)
##     rel_tol <- runif(1, 1e-6, 1e-4)
##     adj_a <- runif(1, 1, 5)
##     gel <- GenEL$new(n, p)
##     gel$set_opts(max_iter = max_iter, rel_tol = rel_tol,
##                  supp_adj = TRUE, supp_adj_a = adj_a)
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     # weights <- weights/sum(weights)
##     logel_cpp <- gel$logel(G, weights)
##     logel_cpp
##     logel_R <- weighted_logEL_G_R(adjG_R(G, adj_a), c(weights, 1),
##                                   max_iter = max_iter, rel_tol = rel_tol)
##     as.numeric(logel_R)
##     if (check_res(logel_R) & check_res(logel_cpp)) {
##       # nconv <<- nconv + 1
##       expect_equal(logel_cpp, as.numeric(logel_R))
##     }
##   }
## })

## # ---- logel_grad ----

## # nconv <- 0
## test_that("dldG with default settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:2, 1)
##     gel <- GenEL$new(n, p)
##     G <- matrix(rnorm(n*p),n,p)
##     dldG_cpp <- gel$logel_grad(G)$dldG
##     dldG_nd <- matrix(tryCatch(numDeriv::grad(logEL_G_R, G, method.args = list(eps = gel$rel_tol)),
##                                error = function(e) {
##                                  rep(NA, nrow(G) * ncol(G))
##                                }), nrow = nrow(G), ncol = ncol(G))
##     if (check_res(dldG_cpp) & check_res(dldG_nd)) {
##       # nconv <<- nconv + 1
##       expect_equal(dldG_cpp, dldG_nd, tolerance = 1e-3)
##     }
##   }
## })

## # nconv <- 0
## test_that("dldG with given convergence settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:2, 1)
##     max_iter <- sample(c(100, 200, 500), 1)
##     rel_tol <- runif(1, 1e-6, 1e-2)
##     gel <- GenEL$new(n, p)
##     gel$max_iter <- max_iter
##     gel$rel_tol <- rel_tol
##     G <- matrix(rnorm(n*p),n,p)
##     dldG_cpp <- gel$logel_grad(G)$dldG
##     logEL_G_R_use <- function(G) {
##       logEL_G_R(G, max_iter = max_iter, rel_tol = rel_tol)
##     }
##     dldG_nd <- matrix(tryCatch(numDeriv::grad(logEL_G_R_use, G, method.args = list(eps = gel$rel_tol)),
##                                error = function(e) {
##                                  rep(NA, nrow(G) * ncol(G))
##                                }), nrow = nrow(G), ncol = ncol(G))
##     if (check_res(dldG_cpp) & check_res(dldG_nd)) {
##       # nconv <<- nconv + 1
##       expect_equal(dldG_cpp, dldG_nd, tolerance = 1e-3)
##     }
##   }
## })

## # nconv <- 0
## test_that("dldG with given convergence settings and support correction", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:2, 1)
##     max_iter <- sample(c(100, 200, 500), 1)
##     rel_tol <- runif(1, 1e-5, 1e-3)
##     adj_a <- runif(1, 1, 5)
##     gel <- GenEL$new(n, p)
##     gel$set_opts(max_iter = max_iter, rel_tol = rel_tol,
##                  supp_adj = TRUE, supp_adj_a = adj_a)
##     G <- matrix(rnorm(n*p),n,p)
##     dldG_cpp <- gel$logel_grad(G)$dldG
##     # dldG_cpp
##     logEL_adjG_R_use <- function(G) {
##       logEL_adjG_R(G, max_iter = max_iter, rel_tol = rel_tol, an = adj_a)
##     }
##     dldG_nd <- matrix(tryCatch(numDeriv::grad(logEL_adjG_R_use, G,
##                                               method.args = list(eps = gel$rel_tol)),
##                                error = function(e) {
##                                  rep(NA, (nrow(G) + 1) * ncol(G))
##                                }), nrow = nrow(G), ncol = ncol(G))
##     # dldG_nd
##     range(dldG_cpp-dldG_nd)
##     if (check_res(dldG_cpp) & check_res(dldG_nd)) {
##       # nconv <<- nconv + 1
##       expect_equal(dldG_cpp, dldG_nd, tolerance = 1e-4)
##     }
##   }
## })

## # ---- logel_grad with weights ----

## # nconv <- 0
## test_that("weights dldG with default settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:2, 1)
##     gel <- GenEL$new(n, p)
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     # weights <- rep(1,n)
##     dldG_cpp <- gel$logel_grad(G, weights)$dldG
##     weighted_logEL_G_R_use <- function(G) {
##       weighted_logEL_G_R(G, weights)
##     }
##     dldG_nd <- matrix(tryCatch(numDeriv::grad(weighted_logEL_G_R_use, G,
##                                               method.args = list(eps = gel$rel_tol)),
##                                error = function(e) {
##                                  rep(NA, nrow(G) * ncol(G))
##                                }), nrow = nrow(G), ncol = ncol(G))
##     if (check_res(dldG_cpp) & check_res(dldG_nd)) {
##       # nconv <<- nconv + 1
##       expect_equal(dldG_cpp, dldG_nd, tolerance = 1e-3)
##     }
##   }
## })
## # nconv

## # nconv <- 0
## test_that("weighted dldG with given convergence settings", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:2, 1)
##     max_iter <- sample(c(100, 200, 500), 1)
##     rel_tol <- runif(1, 1e-6, 1e-2)
##     gel <- GenEL$new(n, p)
##     gel$max_iter <- max_iter
##     gel$rel_tol <- rel_tol
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     dldG_cpp <- gel$logel_grad(G, weights)$dldG
##     weighted_logEL_G_R_use <- function(G) {
##       weighted_logEL_G_R(G, weights, max_iter, rel_tol)
##     }
##     dldG_nd <- matrix(tryCatch(numDeriv::grad(weighted_logEL_G_R_use, G,
##                                               method.args = list(eps = gel$rel_tol)),
##                                error = function(e) {
##                                  rep(NA, nrow(G) * ncol(G))
##                                }), nrow = nrow(G), ncol = ncol(G))
##     if (check_res(dldG_cpp) & check_res(dldG_nd)) {
##       # nconv <<- nconv + 1
##       expect_equal(dldG_cpp, dldG_nd, tolerance = 1e-3)
##     }
##   }
## })

## # nconv <- 0
## test_that("weighted dldG with given convergence settings and support correction", {
##   for (ii in 1:ntest) {
##     n <- sample(10:20,1)
##     p <- sample(1:2, 1)
##     max_iter <- sample(c(100, 200, 500), 1)
##     rel_tol <- runif(1, 1e-5, 1e-3)
##     adj_a <- runif(1, 1, 5)
##     gel <- GenEL$new(n, p)
##     gel$set_opts(max_iter = max_iter, rel_tol = rel_tol,
##                  supp_adj = TRUE, supp_adj_a = adj_a)
##     G <- matrix(rnorm(n*p),n,p)
##     weights <- runif(n, 0, 2)
##     dldG_cpp <- gel$logel_grad(G, weights)$dldG
##     weighted_logEL_G_R_use <- function(G) {
##       weighted_logEL_G_R(adjG_R(G, adj_a), c(weights, 1), max_iter, rel_tol)
##     }
##     dldG_nd <- matrix(tryCatch(numDeriv::grad(weighted_logEL_G_R_use, G,
##                                               method.args = list(eps = gel$rel_tol)),
##                                error = function(e) {
##                                  rep(NA, nrow(G) * ncol(G))
##                                }), nrow = nrow(G), ncol = ncol(G))
##     range(dldG_cpp-dldG_nd)
##     if (check_res(dldG_cpp) & check_res(dldG_nd)) {
##       # nconv <<- nconv + 1
##       expect_equal(dldG_cpp, dldG_nd, tolerance = 1e-4)
##     }
##   }
## })
