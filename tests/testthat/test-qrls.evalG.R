# ---- testing qrls.evalG implementations in R and C++ are equal ----

source("el-utils.R")

n <- 5
p <- 2
q <- 1
X <- matrix(rchisq(p*n,df=2),n,p)
Z <- matrix(rchisq(q*n,df=1),n,q)
beta <- rnorm(p)
gamma <- rnorm(q)
eps <- rnorm(n)
y <- X %*% beta + exp(Z %*% gamma) + eps 
alpha <- 0.5
G_R <- qrls.evalG_R(y,X,Z,beta,gamma,alpha)
G_R
G_cpp <- qrls.evalG(y,X,Z,beta,gamma,alpha)

# ---- testing qrls.evalG implementations in R and C++ are equal ----

source("el-utils.R")
# library(testthat)

context("mrls.evalG")

ntest <- 50

test_that("mr.evalG.R == mr.evalG.cpp", {
  for(ii in 1:ntest) {
    # Location-scale model + mean regression
    n <- sample(10:20,1)
    p <- sample(1:(n-2), 1)
    q <- sample(1:(n-2), 1)
    X <- replicate(p, rnorm(n))
    Z <- replicate(q, rnorm(n))
    beta <- rnorm(p)
    gamma <- rnorm(q)
    y <- c(X %*% beta + exp(Z %*% gamma)) + rnorm(n) # with N(0,1) error term
    max_iter <- sample(c(2, 10, 100), 1)
    rel_tol <- runif(1, 1e-6, 1e-5)
    # checking G matrix from cpp and R
    G.cpp <- mrls.evalG(y,X,Z,beta,gamma)
    G.R <- mrls.evalG_R(y,X,Z,beta,gamma)
    expect_equal(G.R, G.cpp)
  }
})


