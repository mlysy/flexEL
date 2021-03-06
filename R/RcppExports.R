# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

.adjG <- function(G, a) {
    .Call(`_flexEL_adjG`, G, a)
}

.LambdaNR <- function(G, max_iter, rel_tol, support, verbose) {
    .Call(`_flexEL_LambdaNR`, G, max_iter, rel_tol, support, verbose)
}

.OmegaHat <- function(G, max_iter, rel_tol, support, verbose) {
    .Call(`_flexEL_OmegaHat`, G, max_iter, rel_tol, support, verbose)
}

.LogEL <- function(omegas, support) {
    .Call(`_flexEL_LogEL`, omegas, support)
}

.LogELGrad <- function(G, max_iter, rel_tol, support = FALSE, verbose = FALSE) {
    .Call(`_flexEL_LogELGrad`, G, max_iter, rel_tol, support, verbose)
}

.EvalWeights <- function(omegas, deltas, epsilons, support) {
    .Call(`_flexEL_EvalWeights`, omegas, deltas, epsilons, support)
}

.LambdaNRCens <- function(G, weights, max_iter, rel_tol, support, verbose) {
    .Call(`_flexEL_LambdaNRCens`, G, weights, max_iter, rel_tol, support, verbose)
}

.OmegaHatEM <- function(G, omegas_init, deltas, epsilons, max_iter, rel_tol, abs_tol, support, verbose) {
    .Call(`_flexEL_OmegaHatEM`, G, omegas_init, deltas, epsilons, max_iter, rel_tol, abs_tol, support, verbose)
}

.LogELCens <- function(omegas, deltas, epsilons, support) {
    .Call(`_flexEL_LogELCens`, omegas, deltas, epsilons, support)
}

.LogELSmooth <- function(omegas, deltas, epsilons, sp, support) {
    .Call(`_flexEL_LogELSmooth`, omegas, deltas, epsilons, sp, support)
}

.EvalWeightsSmooth <- function(omegas, deltas, epsilons, sp, support) {
    .Call(`_flexEL_EvalWeightsSmooth`, omegas, deltas, epsilons, sp, support)
}

.OmegaHatEMSmooth <- function(G, omegas_init, deltas, epsilons, sp, max_iter, rel_tol, abs_tol, support, verbose) {
    .Call(`_flexEL_OmegaHatEMSmooth`, G, omegas_init, deltas, epsilons, sp, max_iter, rel_tol, abs_tol, support, verbose)
}

.MeanRegEvalG <- function(y, X, beta) {
    .Call(`_flexEL_MeanReg_evalG`, y, X, beta)
}

.MeanRegLSEvalG <- function(y, X, Z, beta, gamma, sig2) {
    .Call(`_flexEL_MeanRegLS_EvalG`, y, X, Z, beta, gamma, sig2)
}

.QuantRegEvalG <- function(y, X, tauArr, beta) {
    .Call(`_flexEL_QuantRegEvalG`, y, X, tauArr, beta)
}

.QuantRegLSEvalG <- function(y, X, Z, tauArr, beta, gamma, sig2, nu) {
    .Call(`_flexEL_QuantRegLSEvalG`, y, X, Z, tauArr, beta, gamma, sig2, nu)
}

.QuantRegLSEvalGSmooth <- function(y, X, Z, tauArr, beta, gamma, sig2, nu, s) {
    .Call(`_flexEL_QuantRegLSEvalGSmooth`, y, X, Z, tauArr, beta, gamma, sig2, nu, s)
}

