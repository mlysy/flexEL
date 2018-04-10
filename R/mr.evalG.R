#' Evaluate the G matrix for mean regression (location model)
#'
#' @param y Length-\code{nObs} vector of responses.
#' @param X \code{nObs x nEqs} matrix of covariates.
#' @param beta Length-\code{nEqs} vector of coefficients in location model.
#' @return G matrix for location mean regression model. 
#' @details ...
#' @export mr.evalG
mr.evalG <- function(y, X, beta) { 
    if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
    if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
    G <- .MeanReg_evalG(y,t(X),beta)
    return(t(G))
}