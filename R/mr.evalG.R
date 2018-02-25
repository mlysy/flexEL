#' Solve the inner loop optimization of an EL function.
#'
#' @param y Length-\code{nObs} vector of responses.
#' @param X \code{nObs x nEqns} matrix of covariates.
#' @param beta Length-\code{nEqs} vector of coefficients in location model.
#' @return G \code{nObs x nEqns} matrix of constraints.
#' @details Returns the G matrix for location mean regression model. 
#' @export
mr.evalG <- function(y, X, beta) { 
    if (!is.vector(y)) stop("y should be a vector.")
    if (nrow(X) != length(y)) stop("y and X have inconsistent dimensions.")
    G <- .MeanReg_evalG(y,t(X),beta)
    return(t(G))
}