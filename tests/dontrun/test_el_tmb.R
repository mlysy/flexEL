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
max_iter <- 1
rel_tol <- 1e-10
supp_adj <- FALSE

logel_adf <- MakeADFun(
  data = list(
    max_iter = max_iter,
    rel_tol = rel_tol,
    supp_adj = supp_adj
  ),
  parameters = list(G = matrix(0, n_eqs, n_obs)),
  silent = TRUE,
  DLL = model
)

## logel_r <- function(G) {
##   logEL(G, max_iter = max_iter,
##         rel_tol = rel_tol, supp_adj = supp_adj)
## }
## logel_grad_r <- function(G) logEL(G, grad = TRUE)$grad

G <- matrix(rnorm(n_obs * n_eqs), n_obs, n_eqs)

ll_r <- logEL(G, max_iter = max_iter, rel_tol = rel_tol,
              supp_adj = supp_adj)
ll_tmb <- logel_adf$f(t(G))
ll_r - ll_tmb

dldg_r <- logEL(G, max_iter = max_iter, rel_tol = rel_tol,
                supp_adj = supp_adj, grad = TRUE)$grad
dldg_tmb <- t(matrix(logel_adf$g(t(G)), n_eqs, n_obs))
dldg_r - dldg_tmb

## numDeriv::grad(func = logel_r, x = G)
