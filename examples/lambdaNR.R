n <- 20
p <- 8
G <- matrix(rnorm(n*p),n,p) # randomly generated G
max_iter <- 100
rel_tol <- 1e-5

# without support correction
lambdaNR(G = G,
         max_iter = max_iter, rel_tol = rel_tol, support = FALSE, 
         verbose = FALSE)

# with support correction
lambdaNR(G = G,
         max_iter = max_iter, rel_tol = rel_tol, support = TRUE, 
         verbose = FALSE)
