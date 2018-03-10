#' Solve the inner loop optimization of an EL function.
#'
#' @param deltas Lenght-\code{nObs} vector of censoring indicators.
#' @param omegas Lenght-\code{nObs} vector, empirical distribution (probability vector).
#' @return Length-\code{nObs} vector of weights for the weighted NR algorithm.
#' @details ... 
#' @export
getWeights <- function(deltas, omegas, epsilons) {
    .getWeights(deltas, omegas, epsilons)
}
# getWeights <- function(y, X, deltas, omegas, beta) {
#     .getWeights(y, t(X), deltas, omegas, beta)
# }