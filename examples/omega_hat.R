# no censoring
n <- 20
p <- 7
G <- matrix(rnorm(n*p),n,p) # random G here
max_iter <- 100
rel_tol <- 1e-5
omega_hat(G = G, max_iter = max_iter, rel_tol = rel_tol, 
          support = FALSE, verbose = FALSE) # without support correction
omega_hat(G = G, max_iter = max_iter, rel_tol = rel_tol, 
          support = TRUE, verbose = FALSE) # with support correction

# right-censoring
n <- 20
p <- 7
G <- matrix(rnorm(n*p),n,p) # random G here
max_iter <- 100
rel_tol <- 1e-5
abs_tol <- 1e-3
deltas <- rep(1,n)
numcens <- 5
censinds <- sample(n,numcens)
deltas[censinds] <- 0
epsilons <- rnorm(n)
omega_hat(G, deltas, epsilons, max_iter = max_iter, 
          rel_tol = rel_tol, abs_tol = abs_tol, support = FALSE, verbose = FALSE) # without support correction
omega_hat(G, deltas, epsilons, max_iter = max_iter, 
          rel_tol = rel_tol, abs_tol = abs_tol, support = TRUE, verbose = FALSE) # with support correction
