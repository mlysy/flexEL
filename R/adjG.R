#' Adjusted G matrix.
#'
#' @param G A matrix of dimension \code{nObs} x \code{nEqs}.
#' @param a A tuning parameter. 
#' @return A matrix of dimension \code{nObs+1} x \code{nEqs}.
#' @details ...
#' @export adjG
adjG <- function(G, a) {
  nObs <- nrow(G)
  if (missing(a)) a <- max(1,0.5*log(nObs))
  return(t(.adjG(t(G),a)))
}