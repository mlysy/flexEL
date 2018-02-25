#' Solve the inner loop optimization of an EL function.
#'
#' @param y Lenght-\code{nObs} vector of responses.
#' @param X \code{nObs x n Eqs} matrix of covariate matrix. 
#' @param deltas Lenght-\code{nObs} vector of censoring indicators.
#' @param omegas Lenght-\code{nObs} vector, empirical distribution (probability vector).
#' @return Length-\code{nObs} vector of weights for the weighted NR algorithm.
#' @details ... 
#' @export
evalWeights <- function(y, X, deltas, omegas, beta) {
    # TODO: do some checking here
    .evalWeights(y, t(X), deltas, omegas, beta)
}