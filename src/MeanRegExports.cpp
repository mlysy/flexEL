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

// Note: sig2 should actually be scaler, but have to put in vector to be 
// consistent with the implementation for quantreg.
// [[Rcpp::export(".MeanRegLS_evalG")]]
Eigen::MatrixXd MeanRegLS_evalG(Eigen::VectorXd y, 
                                Eigen::MatrixXd X, Eigen::MatrixXd Z,
                                Eigen::VectorXd beta, Eigen::VectorXd gamma, 
                                Eigen::VectorXd sig2) {
  // InnerEL<MeanRegModel> MR(y, X, NULL); // instantiate
  InnerEL<MeanRegModel> MR;
  MR.setData(y,X,Z,NULL); 
  MR.evalG(beta,gamma,sig2);
  // MR.evalG(beta,gamma,Eigen::VectorXd::Zero(0));
  Eigen::MatrixXd G = MR.getG(); // G is nEqs x nObs
  return(G); 
}

// [[Rcpp::export(".MeanReg_post")]]
Rcpp::List MeanReg_post(Eigen::VectorXd y, Eigen::MatrixXd X, 
                        int nsamples, int nburn, 
                        Eigen::MatrixXd BetaInit, Eigen::MatrixXd mwgSd, 
                        Eigen::MatrixXd RvDoMcmc, 
                        int maxIter = 100, double relTol = 1e-7) {
  InnerEL<MeanRegModel> MR;
  MR.setData(y,X,NULL); 
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  MR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept; 
  Eigen::MatrixXd Beta_chain = MR.postSample(nsamples, nburn, 
                                             BetaInit, mwgSd, 
                                             RvDoMcmc, paccept);
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
  // Eigen::MatrixXd paccept;
  Eigen::VectorXd paccept;
  // int nTheta = BetaInit.rows();
  // int numTheta = BetaInit.cols();
  // bool *rvdomcmc = new bool[nTheta*numTheta];
  int nTheta = betaInit.size();
  bool *rvdomcmc = new bool[nTheta];
  dVec_to_bArr(rvDoMcmc, rvdomcmc);
  // for (int ii=0; ii<numTheta; ii++) {
  //   dVec_to_bArr(rvDoMcmc.col(ii), &rvdomcmc[ii*nTheta]);
  // }
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
                        Eigen::VectorXd Sig2Init, Eigen::MatrixXd mwgSd, 
                        Eigen::MatrixXd RvDoMcmc, 
                        int maxIter = 100, double relTol = 1e-7) {
  InnerEL<MeanRegModel> MR;
  MR.setData(y,X,Z,NULL);
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  MR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept;
  
  // concatenate (stack) them together
  // Eigen::MatrixXd ThetaInit(BetaInit.rows()+GammaInit.rows(), BetaInit.cols());
  // ThetaInit << BetaInit, GammaInit;
  
  // concatenate (stack) them together
  Eigen::MatrixXd ThetaInit(BetaInit.rows()+GammaInit.rows()+1, BetaInit.cols());
  ThetaInit.topRows(BetaInit.rows()) = BetaInit;
  ThetaInit.block(BetaInit.rows(),0,GammaInit.rows(),GammaInit.cols()) = GammaInit;
  ThetaInit.bottomRows(1) = Sig2Init.transpose();
  
  // std::cout << "MeanRegExports: thetaInit = " << thetaInit.transpose() << std::endl;
  Eigen::MatrixXd Theta_chain = MR.postSample(nsamples, nburn,
                                              ThetaInit, mwgSd,
                                              RvDoMcmc, paccept);
                                              // BetaInit.rows(), GammaInit.rows());
  Rcpp::List retlst; 
  retlst["Theta_chain"] = Theta_chain;
  retlst["paccept"] = paccept;
  return(retlst); 
}
