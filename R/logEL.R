#' log empirical likelihood
#'
#' @param omegas \code{nObs} probability vector. 
#' @param G \code{nObs x nEqs} matrix of constraints.
#' @param deltas \code{nObs} vector of censoring indicators. 
#' @param epsilons \code{nObs} vector of residuals. 
#' @param max_iter Maximum number of Newton-Raphson steps.
#' @param rel_tol Relative tolerance of Newton-Raphson convergence.
#' @return Log empirical likelihood of the input beta
#' @details ...
#' @export
logEL <- function(G, deltas, epsilons, 
                  max_iter = 100, rel_tol = 1e-7, verbose = FALSE) {
    # non-censoing case
    if (missing(deltas) && missing(epsilons)) {
        if (!is.vector(omegas)) {
            stop("omegas must be a vector")
        }
        if(nrow(G) != length(omegas)) {
            stop("G and deltas have inconsistent dimensions.")
        }
        .logEL(G = t(G), maxIter = max_iter, relTol = rel_tol, verbose = verbose)
    }
    # censoring case
    else {
        # input checks
        if(nrow(G) != length(deltas)) {
            stop("G and deltas have inconsistent dimensions.")
        }
        if(nrow(G) != length(epsilons)) {
            stop("G and epsilons have inconsistent dimensions.")
        }
        if(length(deltas) != length(epsilons)) {
            stop("deltas and epsilons have inconsistent lengths.")
        }
        .logELC(G = t(G), deltas = deltas, epsilons = epsilons, 
               maxIter = max_iter, relTol = rel_tol, verbose = verbose)
    }
}