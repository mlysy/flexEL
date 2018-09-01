#' Solve the inner loop optimization of an EL function.
#'
#' @param x Scalar or vector of real values.
#' @param s Scalar or vector of the same length as \code{x}, parameter of the smoothed indicator function.
#' @return Scalar or vector of the same length as \code{x} where the values are between 0 and 1 evaluated by the smoothed indicator function.
#' @details ...
#' @export ind.smooth
ind.smooth <- function(x, s=10) {
  if (any(s <= 0)) stop("Choose a positive s.")
  if (length(x) != 1 && length(s) == 1) s <- rep(s,length(x))
  return(.ind.smooth(x,s))
}
