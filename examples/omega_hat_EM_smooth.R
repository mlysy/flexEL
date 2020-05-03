n <- 20
p <- 7
G <- matrix(rnorm(n*p), n, p)
deltas <- rep(1,n)
numcens <- sample(round(n/2),1)
censinds <- sample(n,numcens)
deltas[censinds] <- 0
epsilons <- rnorm(n)
s <- 10 # parameter for continuity correction
max_iter <- 500
rel_tol <- 1e-5
omega_hat_EM_smooth(G, deltas, epsilons, s, 
                    max_iter = max_iter, rel_tol = rel_tol, abs_tol = rel_tol, 
                    support = FALSE, verbose = FALSE) # without support correction
omega_hat_EM_smooth(G, deltas, epsilons, s, 
                    max_iter = max_iter, rel_tol = rel_tol, abs_tol = rel_tol, 
                    support = TRUE, verbose = FALSE) # with support correction
