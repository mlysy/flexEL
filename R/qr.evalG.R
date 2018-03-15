#' Evaluate the G matrix for quantile regression (location model)
#'
#' @param y Length-\code{nObs} vector of responses.
#' @param X \code{nObs x nEqns} matrix of covariates.
#' @param beta Length-\code{nEqs} vector of coefficients in location model.
#' @details Returns the G matrix for location mean regression model. 
#' @export
qr.evalG <- function(y, X, alpha, beta) { 
    if (!is.vector(y)) stop("y should be a vector.") # TODO: allow y to be 1d matrix too
    if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
    G <- .QuantReg_evalG(y,t(X),alpha,beta)
    return(t(G))
}