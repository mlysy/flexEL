#' Calculates adjusted G matrix.
#'
#' @param G A matrix of dimension `n_obs x n_eq`.
#' @param a A number as tuning parameter.
#' @return A matrix of dimension `(n_obs+1) x n_eq`.
#' @noRd
adjG <- function(G, a) {
  nObs <- nrow(G)
  if (missing(a)) a <- max(1,0.5*log(nObs))
  return(t(.adjG(t(G),a)))
}
