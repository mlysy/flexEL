library(TMB)
library(numDeriv)
library(flexEL)

model <- "el_tmb" # name of cpp file
flexel_path <- system.file("include", package = "flexEL")
TMB::compile(paste0(model, ".cpp"),
             PKG_CXXFLAGS = paste0("-std=gnu++11 -I",
                                   flexel_path))

dyn.load(dynlib(model))

n_obs <- 10
n_eqs <- 3

logel_adf <- MakeADFun(
  data = list(),
  parameters = list(G = matrix(0, n_eqs, n_obs)),
  silent = TRUE,
  DLL = model
)

logel_r <- function(G) logEL(G)
logel_grad_r <- function(G) logEL(G, grad = TRUE)$grad


G <- matrix(rnorm(n_obs * n_eqs), n_obs, n_eqs)

logel_r(G) - logel_adf$f(t(G))

logel_grad_r(G) - t(matrix(logel_adf$g(t(G)), n_eqs, n_obs))
