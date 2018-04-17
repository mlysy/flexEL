// port some of QuantReg to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
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
                                 Eigen::MatrixXd Beta, Eigen::MatrixXd Gamma) {
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate
  InnerEL<QuantRegModel> QR;
  // QR.setData(y,X,Z,&alpha); 
  QR.setData(y,X,Z,alphaArr.data());
  QR.evalG(Beta,Gamma);
  Eigen::MatrixXd G = QR.getG(); // G is nEqs x nObs
  return(G); 
}

// Old code:
// // [[Rcpp::export(".QuantReg_logEL")]]
// double QuantReg_logEL(Eigen::VectorXd y, Eigen::MatrixXd X, 
//                  double alpha, Eigen::VectorXd beta, 
//                  int maxIter = 100, double relTol = 1e-7) {
//     // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate
//     InnerEL<QuantRegModel> QR;
//     QR.setData(y,X,&alpha); 
//     QR.evalG(beta); 
//     // int nIter;
//     // double maxErr;
//     double logELquant = QR.logEL();
//     // TODO: check convergence 
//     return(logELquant);
// }

// now can only sample for sinle quantile
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
  Eigen::VectorXd paccept;
  bool *rvdomcmc = new bool[betaInit.size()];
  dVec_to_bArr(rvDoMcmc, rvdomcmc);
  Eigen::MatrixXd beta_chain = QR.postSampleAdapt(nsamples, nburn,
                                                  betaInit, mwgSd.data(),
                                                  rvdomcmc, paccept);
  delete[] rvdomcmc; 
  Rcpp::List retlst;
  retlst["beta_chain"] = beta_chain;
  retlst["paccept"] = paccept;
  return(retlst);
}

// now can only sample for single quantile
// [[Rcpp::export(".QuantReg_post")]]
Rcpp::List QuantReg_post(Eigen::VectorXd y, Eigen::MatrixXd X, 
                         Eigen::VectorXd alphaArr, int nsamples, int nburn, 
                         Eigen::MatrixXd BetaInit, Eigen::MatrixXd Sigs, 
                         int maxIter = 100, double relTol = 1e-7) {
  InnerEL<QuantRegModel> QR;
  // QR.setData(y,X,&alpha); 
  QR.setData(y,X,alphaArr.data());
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  QR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept; 
  Eigen::MatrixXd Beta_chain = QR.postSample(nsamples, nburn, 
                                             BetaInit, Sigs, paccept, BetaInit.rows());
  // return(Beta_chain);
  Rcpp::List retlst; 
  retlst["Beta_chain"] = Beta_chain;
  retlst["paccept"] = paccept; 
  return(retlst); 
}

// Note: BetaInit and GammaInit must have the same number of columns
// [[Rcpp::export(".QuantRegLS_post")]]
Rcpp::List QuantRegLS_post(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, 
                           Eigen::VectorXd alphaArr, int nsamples, int nburn, 
                           Eigen::MatrixXd BetaInit, Eigen::MatrixXd GammaInit, 
                           Eigen::MatrixXd Sigs, 
                           int maxIter = 100, double relTol = 1e-7) {
  InnerEL<QuantRegModel> QR;
  QR.setData(y,X,Z,alphaArr.data());
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  QR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept;
  // concatenate (stack) them together
  Eigen::MatrixXd ThetaInit(BetaInit.rows()+GammaInit.rows(), BetaInit.cols());
  ThetaInit << BetaInit, GammaInit;
  // std::cout << "MeanRegExports: thetaInit = " << thetaInit.transpose() << std::endl;
  Eigen::MatrixXd Theta_chain = QR.postSample(nsamples, nburn,
                                              ThetaInit, Sigs, paccept, BetaInit.rows());
  Rcpp::List retlst; 
  retlst["Theta_chain"] = Theta_chain;
  retlst["paccept"] = paccept;
  return(retlst); 
}
