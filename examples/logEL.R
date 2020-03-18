
# no censoring
n <- 20
omegas <- abs(rnorm(n))
omegas <- omegas/sum(omegas) # create a random probability vector
logEL(omegas)

# right censoring
n <- 20 # total number of observations
k <- 4 # number of censored observations
idx <- sample(n,k) # indices of censored observations
deltas <- rep(1,n)
deltas[idx] <- 0 # censoring indicators
omegas <- abs(rnorm(n)) 
omegas <- omegas/sum(omegas) # create a random probability vector
epsilons <- rnorm(n) # time to event for each observation
logEL(omegas,epsilons,deltas)
