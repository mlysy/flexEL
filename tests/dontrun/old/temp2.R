source("../testthat/el-utils.R")
load("G1.RData")
lambdaNR_R(G1, verbose = T)
rep(1,nrow(G1)) %*% G1
