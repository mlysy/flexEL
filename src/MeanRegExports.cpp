// port some of MeanReg functions to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
#include "MeanRegModel.h"
#include "dVecTobArr.h"

// [[Rcpp::export(".MeanReg_evalG")]]
Eigen::MatrixXd MeanReg_evalG(Eigen::VectorXd y, Eigen::MatrixXd X, 
                              Eigen::VectorXd beta) {
    // InnerEL<MeanRegModel> MR(y, X, NULL); // instantiate
    InnerEL<MeanRegModel> MR;
    MR.setData(y,X,NULL); 
    MR.evalG(beta);
    Eigen::MatrixXd G = MR.getG(); // G is nEqs x nObs
    return(G); 
}

// [[Rcpp::export(".MeanRegLS_evalG")]]
Eigen::MatrixXd MeanRegLS_evalG(Eigen::VectorXd y, 
                                Eigen::MatrixXd X, Eigen::MatrixXd Z,
                                Eigen::VectorXd beta, Eigen::VectorXd gamma) {
  // InnerEL<MeanRegModel> MR(y, X, NULL); // instantiate
  InnerEL<MeanRegModel> MR;
  MR.setData(y,X,Z,NULL); 
  MR.evalG(beta,gamma);
  Eigen::MatrixXd G = MR.getG(); // G is nEqs x nObs
  return(G); 
}

// Old code:
// // [[Rcpp::export(".MeanReg_logEL")]]
// double MeanReg_logEL(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd beta,
//                      int maxIter = 100, double relTol = 1e-7) {
//     // InnerEL<MeanRegModel> MR(y, X, NULL); // instantiate
//     InnerEL<MeanRegModel> MR;
//     MR.setData(y,X,NULL); 
//     MR.evalG(beta);
//     // int nIter;
//     // double maxErr; 
//     double logELmean = MR.logEL(); 
//     // TODO: check convergence here
//     return(logELmean);
// }

// [[Rcpp::export(".MeanReg_post")]]
Rcpp::List MeanReg_post(Eigen::VectorXd y, Eigen::MatrixXd X, 
                        int nsamples, int nburn, 
                        Eigen::MatrixXd BetaInit, Eigen::MatrixXd Sigs, 
                        int maxIter = 100, double relTol = 1e-7) {
  InnerEL<MeanRegModel> MR;
  MR.setData(y,X,NULL); 
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  MR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept; 
  Eigen::MatrixXd Beta_chain = MR.postSample(nsamples, nburn, 
                                             BetaInit, Sigs, paccept, BetaInit.rows());
  // return(Beta_chain);
  Rcpp::List retlst; 
  retlst["Beta_chain"] = Beta_chain;
  retlst["paccept"] = paccept; 
  return(retlst); 
}

// TODO: rvDoMcmc uses 0/1 double converted to bool for now
// [[Rcpp::export(".MeanReg_post_adapt")]]
Rcpp::List MeanReg_post_adapt(Eigen::VectorXd y, Eigen::MatrixXd X,
                               int nsamples, int nburn,
                               Eigen::VectorXd betaInit, Eigen::VectorXd mwgSd,
                               Eigen::VectorXd rvDoMcmc,
                               int maxIter = 100, double relTol = 1e-7) {
  InnerEL<MeanRegModel> MR;
  MR.setData(y,X,NULL);
  MR.setTol(maxIter, relTol);
  Eigen::VectorXd paccept;
  bool *rvdomcmc = new bool[betaInit.size()];
  dVec_to_bArr(rvDoMcmc, rvdomcmc);
  Eigen::MatrixXd beta_chain = MR.postSampleAdapt(nsamples, nburn,
                                                  betaInit, mwgSd.data(),
                                                  rvdomcmc, paccept);
  delete[] rvdomcmc;
  Rcpp::List retlst;
  retlst["beta_chain"] = beta_chain;
  retlst["paccept"] = paccept;
  return(retlst);
}

// Note: BetaInit and GammaInit must have the same number of columns
// [[Rcpp::export(".MeanRegLS_post")]]
Rcpp::List MeanRegLS_post(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, 
                        int nsamples, int nburn, 
                        Eigen::MatrixXd BetaInit, Eigen::MatrixXd GammaInit, 
                        Eigen::MatrixXd Sigs, 
                        int maxIter = 100, double relTol = 1e-7) {
  InnerEL<MeanRegModel> MR;
  MR.setData(y,X,Z,NULL);
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  MR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept;
  // concatenate (stack) them together
  Eigen::MatrixXd ThetaInit(BetaInit.rows()+GammaInit.rows(), BetaInit.cols());
  ThetaInit << BetaInit, GammaInit;
  // std::cout << "MeanRegExports: thetaInit = " << thetaInit.transpose() << std::endl;
  Eigen::MatrixXd Theta_chain = MR.postSample(nsamples, nburn,
                                              ThetaInit, Sigs, paccept, BetaInit.rows());
  Rcpp::List retlst; 
  retlst["Theta_chain"] = Theta_chain;
  retlst["paccept"] = paccept;
  return(retlst); 
}

/*
// [[Rcpp::export(".MeanReg_PostSample")]]
Eigen::MatrixXd MeanReg_PostSample(Eigen::VectorXd y, Eigen::MatrixXd X, 
                             int nObs, int nEqs, Eigen::VectorXd lambda0, 
                             int nsamples, int nburn, Eigen::VectorXd betaInit,
                             Eigen::VectorXd sigs, int maxIter = 100, double eps = 1e-7)  {
    InnerEL<MeanRegModel> MR(y, X, nObs, nEqs, lambda0); // instantiate
    Eigen::MatrixXd beta_chain = MR.PostSample(nsamples, nburn, betaInit, sigs, maxIter, eps);
    return(beta_chain);
}
*/
