
context("GenEL")

source("el_rfuns.R")

ntest <- 5

# ---- lambda_nr -----

# nconv <- 0
test_that("lambda_nr with default options", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    gel <- GenEL$new(n, p)
    G <- matrix(rnorm(n*p),n,p)
    lambda_cpp <- gel$lambda_nr(G)
    lambda_R_lst <- lambdaNR_R(G)
    if (!lambda_R_lst$convergence) {
      expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
    }
    else {
      # nconv <<- nconv + 1
      expect_equal(lambda_cpp, lambda_R_lst$lambda)
    }
  }
})

# nconv <- 0
test_that("lambda_nr with given convergence settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-2)
    gel <- GenEL$new(n, p)
    gel$max_iter <- max_iter
    gel$rel_tol <- rel_tol
    G <- matrix(rnorm(n*p),n,p)
    lambda_cpp <- gel$lambda_nr(G, verbose = FALSE)
    lambda_R_lst <- lambdaNR_R(G, max_iter = max_iter, rel_tol = rel_tol)
    if (!lambda_R_lst$convergence) {
      expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
    }
    else {
      # nconv <<- nconv + 1
      expect_equal(lambda_cpp, lambda_R_lst$lambda)
    }
  }
})

# nconv <- 0
test_that("lambda_nr with given convergence settings and support correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-2)
    adj_a <- runif(1, 1, 5)
    gel <- GenEL$new(n, p)
    gel$max_iter <- max_iter
    gel$rel_tol <- rel_tol
    gel$supp_adj <- TRUE
    gel$supp_adj_a <- adj_a
    G <- matrix(rnorm(n*p),n,p)
    lambda_cpp <- gel$lambda_nr(G, verbose = FALSE)
    lambda_R_lst <- lambdaNR_R(adjG_R(G, adj_a), max_iter = max_iter, rel_tol = rel_tol)
    if (!lambda_R_lst$convergence) {
      expect_equal(all(is.na(lambda_cpp)), all(is.na(lambda_R_lst$lambda)))
    }
    else {
      # nconv <<- nconv + 1
      expect_equal(lambda_cpp, lambda_R_lst$lambda)
    }
  }
})

# ---- omega_hat ----

# nconv <- 0
test_that("omega_hat with default settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    gel <- GenEL$new(n, p)
    G <- matrix(rnorm(n*p),n,p)
    omega_cpp <- gel$omega_hat(G)
    omega_R <- omega_hat_R(G)
    if (!any(is.na(omega_cpp)) & !any(is.na(omega_R))) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R)
    }
    else {
      expect_equal(any(is.na(omega_cpp)), any(is.na(omega_R)))
    }
  }
})

# nconv <- 0
test_that("omega_hat with given convergence settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-2)
    gel <- GenEL$new(n, p)
    gel$max_iter <- max_iter
    gel$rel_tol <- rel_tol
    G <- matrix(rnorm(n*p),n,p)
    omega_cpp <- gel$omega_hat(G)
    omega_R <- omega_hat_R(G, max_iter = max_iter, rel_tol = rel_tol)
    if (!any(is.na(omega_cpp)) & !any(is.na(omega_R))) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R)
    }
    else {
      expect_equal(any(is.na(omega_cpp)), any(is.na(omega_R)))
    }
  }
})

# nconv <- 0
test_that("omega_hat with default settings and support correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-2)
    adj_a <- runif(1, 1, 5)
    gel <- GenEL$new(n, p)
    gel$max_iter <- max_iter
    gel$rel_tol <- rel_tol
    gel$supp_adj <- TRUE
    gel$supp_adj_a <- adj_a
    G <- matrix(rnorm(n*p),n,p)
    omega_cpp <- gel$omega_hat(G)
    omega_R <- omega_hat_R(adjG_R(G, adj_a), adjust = TRUE,
                           max_iter = max_iter, rel_tol = rel_tol)
    if (!any(is.na(omega_cpp)) & !any(is.na(omega_R))) {
      # nconv <<- nconv + 1
      expect_equal(omega_cpp, omega_R)
    }
    else {
      expect_equal(any(is.na(omega_cpp)), any(is.na(omega_R)))
    }
  }
})

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
    gel <- GenEL$new(n, p)
    G <- matrix(rnorm(n*p),n,p)
    logel_cpp <- gel$logel(G)
    omegahat_R <- omega_hat_R(G)
    logel_R <- logEL_R(omegahat_R)
    if (check_res(logel_R) & check_res(logel_cpp)) {
      # nconv <<- nconv + 1
      expect_equal(logel_cpp, logel_R)
    }
  }
})

# nconv <- 0
test_that("logel with given convergence settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-2)
    gel <- GenEL$new(n, p)
    gel$max_iter <- max_iter
    gel$rel_tol <- rel_tol
    G <- matrix(rnorm(n*p),n,p)
    logel_cpp <- gel$logel(G)
    omegahat_R <- omega_hat_R(G, max_iter = max_iter, rel_tol = rel_tol)
    logel_R <- logEL_R(omegahat_R)
    if (check_res(logel_R) & check_res(logel_cpp)) {
      # nconv <<- nconv + 1
      expect_equal(logel_cpp, logel_R)
    }
  }
})

# nconv <- 0
test_that("logel with given convergence settings and support correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-2)
    adj_a <- runif(1, 1, 5)
    gel <- GenEL$new(n, p)
    gel$max_iter <- max_iter
    gel$rel_tol <- rel_tol
    gel$supp_adj <- TRUE
    gel$supp_adj_a <- adj_a
    G <- matrix(rnorm(n*p),n,p)
    logel_cpp <- gel$logel(G)
    omegahat_R <- omega_hat_R(adjG_R(G, adj_a), max_iter = max_iter, rel_tol = rel_tol)
    logel_R <- logEL_R(omegahat_R, adjust = TRUE)
    if (check_res(logel_R) & check_res(logel_cpp)) {
      # nconv <<- nconv + 1
      expect_equal(logel_cpp, logel_R)
    }
  }
})

# ---- logel_grad ----

# nconv <- 0
test_that("dldG with default settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    gel <- GenEL$new(n, p)
    G <- matrix(rnorm(n*p),n,p)
    dldG_cpp <- gel$logel_grad(G)$dldG
    dldG_nd <- matrix(tryCatch(numDeriv::grad(logEL_G_R, G, method.args = list(eps = gel$rel_tol)),
                               error = function(e) {
                                 rep(NA, nrow(G) * ncol(G))
                               }), nrow = nrow(G), ncol = ncol(G))
    if (check_res(dldG_cpp) & check_res(dldG_nd)) {
      # nconv <<- nconv + 1
      expect_equal(dldG_cpp, dldG_nd, tolerance = 1e-4)
    }
  }
})

# nconv <- 0
test_that("dldG with given convergence settings", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-6, 1e-2)
    gel <- GenEL$new(n, p)
    gel$max_iter <- max_iter
    gel$rel_tol <- rel_tol
    G <- matrix(rnorm(n*p),n,p)
    dldG_cpp <- gel$logel_grad(G)$dldG
    logEL_G_R_use <- function(G) {
      logEL_G_R(G, max_iter = max_iter, rel_tol = rel_tol)
    }
    dldG_nd <- matrix(tryCatch(numDeriv::grad(logEL_G_R_use, G, method.args = list(eps = gel$rel_tol)),
                               error = function(e) {
                                 rep(NA, nrow(G) * ncol(G))
                               }), nrow = nrow(G), ncol = ncol(G))
    if (check_res(dldG_cpp) & check_res(dldG_nd)) {
      # nconv <<- nconv + 1
      expect_equal(dldG_cpp, dldG_nd, tolerance = 1e-4)
    }
  }
})

# nconv <- 0
test_that("dldG with given convergence settings and support correction", {
  for (ii in 1:ntest) {
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    max_iter <- sample(c(100, 200, 500), 1)
    rel_tol <- runif(1, 1e-5, 1e-3)
    adj_a <- runif(1, 1, 5)
    gel <- GenEL$new(n, p)
    gel$max_iter <- max_iter
    gel$rel_tol <- rel_tol
    gel$supp_adj <- TRUE
    gel$supp_adj_a <- adj_a
    G <- matrix(rnorm(n*p),n,p)
    dldG_cpp <- gel$logel_grad(G)$dldG
    # dldG_cpp
    logEL_adjG_R_use <- function(G) {
      logEL_adjG_R(G, max_iter = max_iter, rel_tol = rel_tol, an = adj_a)
    }
    dldG_nd <- matrix(tryCatch(numDeriv::grad(logEL_adjG_R_use, G,
                                              method.args = list(eps = gel$rel_tol)),
                               error = function(e) {
                                 rep(NA, (nrow(G) + 1) * ncol(G))
                               }), nrow = nrow(G), ncol = ncol(G))
    # dldG_nd
    if (check_res(dldG_cpp) & check_res(dldG_nd)) {
      # nconv <<- nconv + 1
      expect_equal(dldG_cpp, dldG_nd, tolerance = 1e-3)
    }
  }
})


