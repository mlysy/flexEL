# no censoring
n <- 20
p <- 5
max_iter <- 200
rel_tol <- 1e-4
G <- matrix(rnorm(n*p),n,p) # random G here
flexEL::logEL(G = G, support = FALSE,
              max_iter = max_iter, rel_tol = rel_tol, abs_tol = 1e-3,
              return_omega = FALSE, verbose = FALSE)

# no censoring, get the derivative of log EL w.r.t. G
flexEL::logEL(G = G, support = FALSE,
              max_iter = max_iter, rel_tol = rel_tol, abs_tol = 1e-3,
              return_omega = FALSE, return_dldG = TRUE, verbose = FALSE)$dldG

# right censoring
n <- 20
p <- 5
max_iter <- 200
rel_tol <- 1e-4
abs_tol <- 1e-3
G <- matrix(rnorm(n*p), n, p)
deltas <- rep(1,n)
numcens <- sample(round(n/2),1)
censinds <- sample(n,numcens)
deltas[censinds] <- 0
epsilons <- rnorm(n)
flexEL::logEL(G = G, delta = deltas, eps = epsilons, support = FALSE, 
              max_iter = max_iter, rel_tol = rel_tol, abs_tol = abs_tol, verbose = FALSE)
