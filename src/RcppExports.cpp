// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// evalWeights
Eigen::VectorXd evalWeights(Eigen::VectorXd deltas, Eigen::VectorXd omegas, Eigen::VectorXd epsilons);
RcppExport SEXP _bayesEL_evalWeights(SEXP deltasSEXP, SEXP omegasSEXP, SEXP epsilonsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type deltas(deltasSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type omegas(omegasSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type epsilons(epsilonsSEXP);
    rcpp_result_gen = Rcpp::wrap(evalWeights(deltas, omegas, epsilons));
    return rcpp_result_gen;
END_RCPP
}
// lambdaNRC
Eigen::VectorXd lambdaNRC(Eigen::MatrixXd G, Eigen::VectorXd weights, int maxIter, double relTol, bool verbose);
RcppExport SEXP _bayesEL_lambdaNRC(SEXP GSEXP, SEXP weightsSEXP, SEXP maxIterSEXP, SEXP relTolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(lambdaNRC(G, weights, maxIter, relTol, verbose));
    return rcpp_result_gen;
END_RCPP
}
// omegaHatEM
Eigen::VectorXd omegaHatEM(Eigen::VectorXd omegasInit, Eigen::MatrixXd G, Eigen::VectorXd deltas, Eigen::VectorXd epsilons, int maxIter, double relTol, bool verbose);
RcppExport SEXP _bayesEL_omegaHatEM(SEXP omegasInitSEXP, SEXP GSEXP, SEXP deltasSEXP, SEXP epsilonsSEXP, SEXP maxIterSEXP, SEXP relTolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type omegasInit(omegasInitSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type deltas(deltasSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type epsilons(epsilonsSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(omegaHatEM(omegasInit, G, deltas, epsilons, maxIter, relTol, verbose));
    return rcpp_result_gen;
END_RCPP
}
// logELC
double logELC(Eigen::VectorXd omegas, Eigen::MatrixXd G, Eigen::VectorXd deltas, Eigen::VectorXd epsilons);
RcppExport SEXP _bayesEL_logELC(SEXP omegasSEXP, SEXP GSEXP, SEXP deltasSEXP, SEXP epsilonsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type omegas(omegasSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type deltas(deltasSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type epsilons(epsilonsSEXP);
    rcpp_result_gen = Rcpp::wrap(logELC(omegas, G, deltas, epsilons));
    return rcpp_result_gen;
END_RCPP
}
// lambdaNR
Eigen::VectorXd lambdaNR(Eigen::MatrixXd G, int maxIter, double relTol, bool verbose);
RcppExport SEXP _bayesEL_lambdaNR(SEXP GSEXP, SEXP maxIterSEXP, SEXP relTolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(lambdaNR(G, maxIter, relTol, verbose));
    return rcpp_result_gen;
END_RCPP
}
// omegaHat
Eigen::VectorXd omegaHat(Eigen::MatrixXd G, Eigen::VectorXd lambda);
RcppExport SEXP _bayesEL_omegaHat(SEXP GSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(omegaHat(G, lambda));
    return rcpp_result_gen;
END_RCPP
}
// logEL
double logEL(Eigen::VectorXd omegas);
RcppExport SEXP _bayesEL_logEL(SEXP omegasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type omegas(omegasSEXP);
    rcpp_result_gen = Rcpp::wrap(logEL(omegas));
    return rcpp_result_gen;
END_RCPP
}
// MeanReg_evalG
Eigen::MatrixXd MeanReg_evalG(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd beta);
RcppExport SEXP _bayesEL_MeanReg_evalG(SEXP ySEXP, SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(MeanReg_evalG(y, X, beta));
    return rcpp_result_gen;
END_RCPP
}
// MeanRegLS_evalG
Eigen::MatrixXd MeanRegLS_evalG(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, Eigen::VectorXd beta, Eigen::VectorXd gamma, Eigen::VectorXd sig2);
RcppExport SEXP _bayesEL_MeanRegLS_evalG(SEXP ySEXP, SEXP XSEXP, SEXP ZSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP sig2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type sig2(sig2SEXP);
    rcpp_result_gen = Rcpp::wrap(MeanRegLS_evalG(y, X, Z, beta, gamma, sig2));
    return rcpp_result_gen;
END_RCPP
}
// MeanReg_post
Rcpp::List MeanReg_post(Eigen::VectorXd y, Eigen::MatrixXd X, int nsamples, int nburn, Eigen::MatrixXd BetaInit, Eigen::MatrixXd mwgSd, Eigen::MatrixXd RvDoMcmc, int maxIter, double relTol);
RcppExport SEXP _bayesEL_MeanReg_post(SEXP ySEXP, SEXP XSEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP BetaInitSEXP, SEXP mwgSdSEXP, SEXP RvDoMcmcSEXP, SEXP maxIterSEXP, SEXP relTolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type BetaInit(BetaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type mwgSd(mwgSdSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type RvDoMcmc(RvDoMcmcSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    rcpp_result_gen = Rcpp::wrap(MeanReg_post(y, X, nsamples, nburn, BetaInit, mwgSd, RvDoMcmc, maxIter, relTol));
    return rcpp_result_gen;
END_RCPP
}
// MeanReg_post_adapt
Rcpp::List MeanReg_post_adapt(Eigen::VectorXd y, Eigen::MatrixXd X, int nsamples, int nburn, Eigen::VectorXd betaInit, Eigen::VectorXd mwgSd, Eigen::VectorXd rvDoMcmc, int maxIter, double relTol);
RcppExport SEXP _bayesEL_MeanReg_post_adapt(SEXP ySEXP, SEXP XSEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP betaInitSEXP, SEXP mwgSdSEXP, SEXP rvDoMcmcSEXP, SEXP maxIterSEXP, SEXP relTolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mwgSd(mwgSdSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type rvDoMcmc(rvDoMcmcSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    rcpp_result_gen = Rcpp::wrap(MeanReg_post_adapt(y, X, nsamples, nburn, betaInit, mwgSd, rvDoMcmc, maxIter, relTol));
    return rcpp_result_gen;
END_RCPP
}
// MeanRegLS_post
Rcpp::List MeanRegLS_post(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, int nsamples, int nburn, Eigen::MatrixXd BetaInit, Eigen::MatrixXd GammaInit, Eigen::VectorXd Sig2Init, Eigen::MatrixXd mwgSd, Eigen::MatrixXd RvDoMcmc, int maxIter, double relTol);
RcppExport SEXP _bayesEL_MeanRegLS_post(SEXP ySEXP, SEXP XSEXP, SEXP ZSEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP BetaInitSEXP, SEXP GammaInitSEXP, SEXP Sig2InitSEXP, SEXP mwgSdSEXP, SEXP RvDoMcmcSEXP, SEXP maxIterSEXP, SEXP relTolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type BetaInit(BetaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type GammaInit(GammaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type Sig2Init(Sig2InitSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type mwgSd(mwgSdSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type RvDoMcmc(RvDoMcmcSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    rcpp_result_gen = Rcpp::wrap(MeanRegLS_post(y, X, Z, nsamples, nburn, BetaInit, GammaInit, Sig2Init, mwgSd, RvDoMcmc, maxIter, relTol));
    return rcpp_result_gen;
END_RCPP
}
// MeanRegLS_post_adapt
Rcpp::List MeanRegLS_post_adapt(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, int nsamples, int nburn, Eigen::VectorXd betaInit, Eigen::VectorXd gammaInit, Eigen::VectorXd sig2Init, Eigen::VectorXd mwgSd, Eigen::VectorXd rvDoMcmc, int maxIter, double relTol);
RcppExport SEXP _bayesEL_MeanRegLS_post_adapt(SEXP ySEXP, SEXP XSEXP, SEXP ZSEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP betaInitSEXP, SEXP gammaInitSEXP, SEXP sig2InitSEXP, SEXP mwgSdSEXP, SEXP rvDoMcmcSEXP, SEXP maxIterSEXP, SEXP relTolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type gammaInit(gammaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type sig2Init(sig2InitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mwgSd(mwgSdSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type rvDoMcmc(rvDoMcmcSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    rcpp_result_gen = Rcpp::wrap(MeanRegLS_post_adapt(y, X, Z, nsamples, nburn, betaInit, gammaInit, sig2Init, mwgSd, rvDoMcmc, maxIter, relTol));
    return rcpp_result_gen;
END_RCPP
}
// MeanRegCens_post
Rcpp::List MeanRegCens_post(Eigen::VectorXd omegasInit, Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd deltas, int nsamples, int nburn, Eigen::VectorXd betaInit, Eigen::VectorXd mwgSd, Eigen::VectorXd rvDoMcmc, int maxIter, double relTol);
RcppExport SEXP _bayesEL_MeanRegCens_post(SEXP omegasInitSEXP, SEXP ySEXP, SEXP XSEXP, SEXP deltasSEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP betaInitSEXP, SEXP mwgSdSEXP, SEXP rvDoMcmcSEXP, SEXP maxIterSEXP, SEXP relTolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type omegasInit(omegasInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type deltas(deltasSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mwgSd(mwgSdSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type rvDoMcmc(rvDoMcmcSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    rcpp_result_gen = Rcpp::wrap(MeanRegCens_post(omegasInit, y, X, deltas, nsamples, nburn, betaInit, mwgSd, rvDoMcmc, maxIter, relTol));
    return rcpp_result_gen;
END_RCPP
}
// QuantReg_evalG
Eigen::MatrixXd QuantReg_evalG(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd alphaArr, Eigen::MatrixXd Beta);
RcppExport SEXP _bayesEL_QuantReg_evalG(SEXP ySEXP, SEXP XSEXP, SEXP alphaArrSEXP, SEXP BetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type alphaArr(alphaArrSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Beta(BetaSEXP);
    rcpp_result_gen = Rcpp::wrap(QuantReg_evalG(y, X, alphaArr, Beta));
    return rcpp_result_gen;
END_RCPP
}
// QuantRegLS_evalG
Eigen::MatrixXd QuantRegLS_evalG(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, Eigen::VectorXd alphaArr, Eigen::MatrixXd Beta, Eigen::MatrixXd Gamma, Eigen::VectorXd Nu);
RcppExport SEXP _bayesEL_QuantRegLS_evalG(SEXP ySEXP, SEXP XSEXP, SEXP ZSEXP, SEXP alphaArrSEXP, SEXP BetaSEXP, SEXP GammaSEXP, SEXP NuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type alphaArr(alphaArrSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Beta(BetaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type Nu(NuSEXP);
    rcpp_result_gen = Rcpp::wrap(QuantRegLS_evalG(y, X, Z, alphaArr, Beta, Gamma, Nu));
    return rcpp_result_gen;
END_RCPP
}
// QuantReg_post
Rcpp::List QuantReg_post(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd alphaArr, int nsamples, int nburn, Eigen::MatrixXd BetaInit, Eigen::MatrixXd Sigs, Eigen::MatrixXd RvDoMcmc, int maxIter, double relTol);
RcppExport SEXP _bayesEL_QuantReg_post(SEXP ySEXP, SEXP XSEXP, SEXP alphaArrSEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP BetaInitSEXP, SEXP SigsSEXP, SEXP RvDoMcmcSEXP, SEXP maxIterSEXP, SEXP relTolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type alphaArr(alphaArrSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type BetaInit(BetaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Sigs(SigsSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type RvDoMcmc(RvDoMcmcSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    rcpp_result_gen = Rcpp::wrap(QuantReg_post(y, X, alphaArr, nsamples, nburn, BetaInit, Sigs, RvDoMcmc, maxIter, relTol));
    return rcpp_result_gen;
END_RCPP
}
// QuantReg_post_adapt
Rcpp::List QuantReg_post_adapt(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd alphaArr, int nsamples, int nburn, Eigen::VectorXd betaInit, Eigen::VectorXd mwgSd, Eigen::VectorXd rvDoMcmc, int maxIter, double relTol);
RcppExport SEXP _bayesEL_QuantReg_post_adapt(SEXP ySEXP, SEXP XSEXP, SEXP alphaArrSEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP betaInitSEXP, SEXP mwgSdSEXP, SEXP rvDoMcmcSEXP, SEXP maxIterSEXP, SEXP relTolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type alphaArr(alphaArrSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mwgSd(mwgSdSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type rvDoMcmc(rvDoMcmcSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    rcpp_result_gen = Rcpp::wrap(QuantReg_post_adapt(y, X, alphaArr, nsamples, nburn, betaInit, mwgSd, rvDoMcmc, maxIter, relTol));
    return rcpp_result_gen;
END_RCPP
}
// QuantRegLS_post
Rcpp::List QuantRegLS_post(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, Eigen::VectorXd alphaArr, int nsamples, int nburn, Eigen::MatrixXd BetaInit, Eigen::MatrixXd GammaInit, Eigen::VectorXd NuInit, Eigen::MatrixXd Sigs, Eigen::MatrixXd RvDoMcmc, int maxIter, double relTol);
RcppExport SEXP _bayesEL_QuantRegLS_post(SEXP ySEXP, SEXP XSEXP, SEXP ZSEXP, SEXP alphaArrSEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP BetaInitSEXP, SEXP GammaInitSEXP, SEXP NuInitSEXP, SEXP SigsSEXP, SEXP RvDoMcmcSEXP, SEXP maxIterSEXP, SEXP relTolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type alphaArr(alphaArrSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type BetaInit(BetaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type GammaInit(GammaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type NuInit(NuInitSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Sigs(SigsSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type RvDoMcmc(RvDoMcmcSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    rcpp_result_gen = Rcpp::wrap(QuantRegLS_post(y, X, Z, alphaArr, nsamples, nburn, BetaInit, GammaInit, NuInit, Sigs, RvDoMcmc, maxIter, relTol));
    return rcpp_result_gen;
END_RCPP
}
// QuantRegLS_post_adapt
Rcpp::List QuantRegLS_post_adapt(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, Eigen::VectorXd alphaArr, int nsamples, int nburn, Eigen::VectorXd betaInit, Eigen::VectorXd gammaInit, Eigen::VectorXd nuInit, Eigen::VectorXd mwgSd, Eigen::VectorXd rvDoMcmc, int maxIter, double relTol);
RcppExport SEXP _bayesEL_QuantRegLS_post_adapt(SEXP ySEXP, SEXP XSEXP, SEXP ZSEXP, SEXP alphaArrSEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP betaInitSEXP, SEXP gammaInitSEXP, SEXP nuInitSEXP, SEXP mwgSdSEXP, SEXP rvDoMcmcSEXP, SEXP maxIterSEXP, SEXP relTolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type alphaArr(alphaArrSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type gammaInit(gammaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type nuInit(nuInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mwgSd(mwgSdSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type rvDoMcmc(rvDoMcmcSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    rcpp_result_gen = Rcpp::wrap(QuantRegLS_post_adapt(y, X, Z, alphaArr, nsamples, nburn, betaInit, gammaInit, nuInit, mwgSd, rvDoMcmc, maxIter, relTol));
    return rcpp_result_gen;
END_RCPP
}
// QuantRegCens_post
Rcpp::List QuantRegCens_post(Eigen::VectorXd omegasInit, Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd deltas, Eigen::VectorXd alphaArr, int nsamples, int nburn, Eigen::VectorXd betaInit, Eigen::VectorXd mwgSd, Eigen::VectorXd rvDoMcmc, int maxIter, double relTol);
RcppExport SEXP _bayesEL_QuantRegCens_post(SEXP omegasInitSEXP, SEXP ySEXP, SEXP XSEXP, SEXP deltasSEXP, SEXP alphaArrSEXP, SEXP nsamplesSEXP, SEXP nburnSEXP, SEXP betaInitSEXP, SEXP mwgSdSEXP, SEXP rvDoMcmcSEXP, SEXP maxIterSEXP, SEXP relTolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type omegasInit(omegasInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type deltas(deltasSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type alphaArr(alphaArrSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mwgSd(mwgSdSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type rvDoMcmc(rvDoMcmcSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type relTol(relTolSEXP);
    rcpp_result_gen = Rcpp::wrap(QuantRegCens_post(omegasInit, y, X, deltas, alphaArr, nsamples, nburn, betaInit, mwgSd, rvDoMcmc, maxIter, relTol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayesEL_evalWeights", (DL_FUNC) &_bayesEL_evalWeights, 3},
    {"_bayesEL_lambdaNRC", (DL_FUNC) &_bayesEL_lambdaNRC, 5},
    {"_bayesEL_omegaHatEM", (DL_FUNC) &_bayesEL_omegaHatEM, 7},
    {"_bayesEL_logELC", (DL_FUNC) &_bayesEL_logELC, 4},
    {"_bayesEL_lambdaNR", (DL_FUNC) &_bayesEL_lambdaNR, 4},
    {"_bayesEL_omegaHat", (DL_FUNC) &_bayesEL_omegaHat, 2},
    {"_bayesEL_logEL", (DL_FUNC) &_bayesEL_logEL, 1},
    {"_bayesEL_MeanReg_evalG", (DL_FUNC) &_bayesEL_MeanReg_evalG, 3},
    {"_bayesEL_MeanRegLS_evalG", (DL_FUNC) &_bayesEL_MeanRegLS_evalG, 6},
    {"_bayesEL_MeanReg_post", (DL_FUNC) &_bayesEL_MeanReg_post, 9},
    {"_bayesEL_MeanReg_post_adapt", (DL_FUNC) &_bayesEL_MeanReg_post_adapt, 9},
    {"_bayesEL_MeanRegLS_post", (DL_FUNC) &_bayesEL_MeanRegLS_post, 12},
    {"_bayesEL_MeanRegLS_post_adapt", (DL_FUNC) &_bayesEL_MeanRegLS_post_adapt, 12},
    {"_bayesEL_MeanRegCens_post", (DL_FUNC) &_bayesEL_MeanRegCens_post, 11},
    {"_bayesEL_QuantReg_evalG", (DL_FUNC) &_bayesEL_QuantReg_evalG, 4},
    {"_bayesEL_QuantRegLS_evalG", (DL_FUNC) &_bayesEL_QuantRegLS_evalG, 7},
    {"_bayesEL_QuantReg_post", (DL_FUNC) &_bayesEL_QuantReg_post, 10},
    {"_bayesEL_QuantReg_post_adapt", (DL_FUNC) &_bayesEL_QuantReg_post_adapt, 10},
    {"_bayesEL_QuantRegLS_post", (DL_FUNC) &_bayesEL_QuantRegLS_post, 13},
    {"_bayesEL_QuantRegLS_post_adapt", (DL_FUNC) &_bayesEL_QuantRegLS_post_adapt, 13},
    {"_bayesEL_QuantRegCens_post", (DL_FUNC) &_bayesEL_QuantRegCens_post, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayesEL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
