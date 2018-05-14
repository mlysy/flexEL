// port some of QuantReg to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
#include "InnerELC.h"
#include "QuantRegModel.h"
#include "dVecTobArr.h"

// double alpha
// [[Rcpp::export(".QuantReg_evalG")]]
Eigen::MatrixXd QuantReg_evalG(Eigen::VectorXd y, Eigen::MatrixXd X,
                               Eigen::VectorXd alphaArr, Eigen::MatrixXd Beta) {
  // std::cout << "print col of vector? " << Beta.col(0).transpose() << std::endl;
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate
  InnerEL<QuantRegModel> QR;
  // QR.setData(y,X,&alpha); 
  QR.setData(y,X,alphaArr.data());
  QR.evalG(Beta);
  Eigen::MatrixXd G = QR.getG(); // G is nEqs x nObs
  return(G); 
}

// [[Rcpp::export(".QuantRegLS_evalG")]]
Eigen::MatrixXd QuantRegLS_evalG(Eigen::VectorXd y, 
                                 Eigen::MatrixXd X, Eigen::MatrixXd Z, 
                                 Eigen::VectorXd alphaArr, 
                                 Eigen::MatrixXd Beta, 
                                 Eigen::MatrixXd Gamma,
                                 Eigen::VectorXd Nu) {
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate
  InnerEL<QuantRegModel> QR;
  // QR.setData(y,X,Z,&alpha); 
  QR.setData(y,X,Z,alphaArr.data());
  QR.evalG(Beta,Gamma,Nu);
  Eigen::MatrixXd G = QR.getG(); // G is nEqs x nObs
  return(G); 
}

// [[Rcpp::export(".QuantReg_post")]]
Rcpp::List QuantReg_post(Eigen::VectorXd y, Eigen::MatrixXd X, 
                         Eigen::VectorXd alphaArr, int nsamples, int nburn, 
                         Eigen::MatrixXd BetaInit, Eigen::MatrixXd Sigs, 
                         Eigen::MatrixXd RvDoMcmc, 
                         int maxIter = 100, double relTol = 1e-7) {
  InnerEL<QuantRegModel> QR;
  // QR.setData(y,X,&alpha); 
  QR.setData(y,X,alphaArr.data());
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  QR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept; 
  Eigen::MatrixXd Beta_chain = QR.postSample(nsamples, nburn, 
                                             BetaInit, Sigs, 
                                             RvDoMcmc, paccept); 
  // BetaInit.rows(), 0);
  // return(Beta_chain);
  Rcpp::List retlst; 
  retlst["Beta_chain"] = Beta_chain;
  retlst["paccept"] = paccept; 
  return(retlst); 
}

// for location models
// TODO: rvDoMcmc uses 0/1 double converted to bool for now
// [[Rcpp::export(".QuantReg_post_adapt")]]
Rcpp::List QuantReg_post_adapt(Eigen::VectorXd y, Eigen::MatrixXd X, 
                               Eigen::VectorXd alphaArr, int nsamples, int nburn, 
                               Eigen::VectorXd betaInit, Eigen::VectorXd mwgSd, 
                               Eigen::VectorXd rvDoMcmc,  
                               int maxIter = 100, double relTol = 1e-7) {
  InnerEL<QuantRegModel> QR;
  // QR.setData(y,X,&alpha); 
  QR.setData(y,X,alphaArr.data());
  QR.setTol(maxIter, relTol);
  // Eigen::MatrixXd paccept;
  Eigen::VectorXd paccept;
  int nTheta = betaInit.size();
  bool *rvdomcmc = new bool[nTheta];
  dVec_to_bArr(rvDoMcmc, rvdomcmc);
  // int nTheta = BetaInit.rows();
  // int numTheta = BetaInit.cols();
  // bool *rvdomcmc = new bool[nTheta*numTheta];
  // for (int ii=0; ii<numTheta; ii++) {
  //   dVec_to_bArr(rvDoMcmc.col(ii), &rvdomcmc[ii*nTheta]);
  //   // std::cout << "rvdomcmc[ii*nTheta] = " << rvdomcmc[ii*nTheta] << std::endl;
  // }
  Eigen::MatrixXd beta_chain = QR.postSampleAdapt(nsamples, nburn,
                                                  betaInit, mwgSd.data(),
                                                  rvdomcmc, paccept);
  delete[] rvdomcmc; 
  Rcpp::List retlst;
  retlst["beta_chain"] = beta_chain;
  retlst["paccept"] = paccept;
  return(retlst);
}

// Note: BetaInit and GammaInit must have the same number of columns
// [[Rcpp::export(".QuantRegLS_post")]]
Rcpp::List QuantRegLS_post(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, 
                           Eigen::VectorXd alphaArr, int nsamples, int nburn, 
                           Eigen::MatrixXd BetaInit, Eigen::MatrixXd GammaInit, 
                           Eigen::VectorXd NuInit, Eigen::MatrixXd Sigs, 
                           Eigen::MatrixXd RvDoMcmc, 
                           int maxIter = 100, double relTol = 1e-7) {
  InnerEL<QuantRegModel> QR;
  QR.setData(y,X,Z,alphaArr.data());
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  QR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept;
  // concatenate (stack) them together
  Eigen::MatrixXd ThetaInit(BetaInit.rows()+GammaInit.rows()+1, BetaInit.cols());
  ThetaInit.topRows(BetaInit.rows()) = BetaInit;
  ThetaInit.block(BetaInit.rows(),0,GammaInit.rows(),GammaInit.cols()) = GammaInit;
  ThetaInit.bottomRows(1) = NuInit.transpose();
  // std::cout << "QuantRegExports: ThetaInit = \n" << ThetaInit << std::endl;
  Eigen::MatrixXd Theta_chain = QR.postSample(nsamples, nburn,
                                              ThetaInit, Sigs, 
                                              RvDoMcmc, paccept);
                                              // BetaInit.rows(), GammaInit.rows());
  Rcpp::List retlst; 
  retlst["Theta_chain"] = Theta_chain;
  retlst["paccept"] = paccept;
  return(retlst); 
}

// Note: BetaInit and GammaInit must have the same number of columns
// [[Rcpp::export(".QuantRegLS_post_adapt")]]
Rcpp::List QuantRegLS_post_adapt(Eigen::VectorXd y, Eigen::MatrixXd X,
                                 Eigen::MatrixXd Z, Eigen::VectorXd alphaArr, 
                                 int nsamples, int nburn,
                                 Eigen::VectorXd betaInit, Eigen::VectorXd gammaInit, 
                                 Eigen::VectorXd nuInit,
                                 Eigen::VectorXd mwgSd, Eigen::VectorXd rvDoMcmc,
                                 int maxIter = 100, double relTol = 1e-7) {
  InnerEL<QuantRegModel> QR;
  QR.setData(y,X,Z,alphaArr.data());
  QR.setTol(maxIter, relTol);
  // Eigen::MatrixXd paccept;
  Eigen::VectorXd paccept;
  // int nTheta = BetaInit.rows();
  // int numTheta = BetaInit.cols();
  // bool *rvdomcmc = new bool[nTheta*numTheta];
  
  // concatenate (stack) them together
  Eigen::VectorXd thetaInit(betaInit.size()+gammaInit.size()+1, betaInit.size());
  thetaInit.head(betaInit.size()) = betaInit;
  thetaInit.segment(betaInit.size(),gammaInit.size()) = gammaInit;
  thetaInit.tail(1) = nuInit;
  
  int nTheta = thetaInit.size();
  bool *rvdomcmc = new bool[nTheta];
  dVec_to_bArr(rvDoMcmc, rvdomcmc);
  // for (int ii=0; ii<numTheta; ii++) {
  //   dVec_to_bArr(rvDoMcmc.col(ii), &rvdomcmc[ii*nTheta]);
  // }
  Eigen::MatrixXd theta_chain = QR.postSampleAdapt(nsamples, nburn,
                                                   thetaInit, mwgSd.data(),
                                                   rvdomcmc, paccept);
  delete[] rvdomcmc;
  Rcpp::List retlst;
  retlst["theta_chain"] = theta_chain;
  retlst["paccept"] = paccept;
  return(retlst);
}

// [[Rcpp::export(".QuantRegCens_post")]]
Rcpp::List QuantRegCens_post(Eigen::VectorXd omegasInit,
                             Eigen::VectorXd y, Eigen::MatrixXd X,
                             Eigen::VectorXd deltas,
                             Eigen::VectorXd alphaArr,
                             int nsamples, int nburn,
                             Eigen::VectorXd betaInit, Eigen::VectorXd mwgSd,
                             Eigen::VectorXd rvDoMcmc,
                             int maxIter = 100, double relTol = 1e-7) {
  InnerELC<QuantRegModel> QRC;
  QRC.setData(y,X,deltas,alphaArr.data());
  QRC.setTol(maxIter, relTol);
  Eigen::VectorXd paccept;
  QRC.setOmegas(omegasInit); // set initial value for first EM
  Eigen::MatrixXd beta_chain = QRC.postSample(nsamples, nburn,
                                              betaInit, mwgSd,
                                              rvDoMcmc, paccept);
  Rcpp::List retlst;
  retlst["beta_chain"] = beta_chain;
  retlst["paccept"] = paccept;
  return(retlst);
}

