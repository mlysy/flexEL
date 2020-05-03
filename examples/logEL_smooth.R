n <- 20 # total number of observations
s <- 10 # parameter for continuity correction
omegas <- abs(rnorm(n))
omegas <- omegas/sum(omegas) # a probability vector
epsilons <- rnorm(n)
deltas <- rep(1,n)
numcens <- 5
deltas[sample(1:n,numcens)] <- 0 # censoring indicators
logEL_smooth(omegas,epsilons,deltas,s)
