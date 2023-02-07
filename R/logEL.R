#' Calculate the log empirical likelihood for a given `G` matrix.
#'
#' A simplified wrapper to `GenEL$logEL()` and `GenEL$logEL_grad()`.  See [`GenEL`] documentation for more options.
#'
#' @param G A numeric matrix of size `n_obs x n_eqs`.
#' @param max_iter A positive integer controlling the maximum number of iterations.
#' @param rel_tol A small positive number controlling accuracy at convergence.
#' @param supp_adj A boolean indicating whether to conduct support correction or not.
#' @param grad Whether or not to return the gradient of the log empirical likelihood with respect to `G`.
#' @return The log empirical likelihood evaluated at `G`, or if `grad == TRUE`, a list with elements `logel` and `grad`, with the latter being a matrix of size `n_obs x n_eqs`.
#' @export
logEL <- function(G, max_iter = 100, rel_tol = 1e-7, supp_adj = FALSE,
                  grad = FALSE) {
  el <- flexEL::GenEL$new(nrow(G), ncol(G))
  el$set_opts(max_iter = max_iter, rel_tol = rel_tol, supp_adj = supp_adj)
  if(grad) {
    out <- el$logel_grad(G)
  } else {
    out <- el$logel(G)
  }
  out
}
