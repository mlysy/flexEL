// port some of QuantReg to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
#include "QuantRegModel.h"

// [[Rcpp::export(".QuantReg_evalG")]]
Eigen::MatrixXd QuantReg_evalG(Eigen::VectorXd y, Eigen::MatrixXd X,
                              double alpha, Eigen::VectorXd beta) {
    // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate
    InnerEL<QuantRegModel> QR;
    QR.setData(y,X,&alpha); 
    QR.evalG(beta);
    Eigen::MatrixXd G = QR.getG(); // G is nEqs x nObs
    return(G); 
}

// [[Rcpp::export(".QuantRegLS_evalG")]]
Eigen::MatrixXd QuantRegLS_evalG(Eigen::VectorXd y, 
                                 Eigen::MatrixXd X, Eigen::MatrixXd Z, 
                                 double alpha, 
                                 Eigen::VectorXd beta, Eigen::VectorXd gamma) {
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate
  InnerEL<QuantRegModel> QR;
  QR.setData(y,X,Z,&alpha); 
  QR.evalG(beta,gamma);
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

// convert Eigen VectorXd dVec to C++ bool array bArr
// bArr must have the same length as dVec
void dVec_to_bArr(Eigen::VectorXd dVec, bool *bArr) {
  // std::cout << "dVec.size() = " << dVec.size() << std::endl;
  for (int ii=0; ii<dVec.size(); ii++) {
    if (dVec(ii) == 0) bArr[ii] = false;
    else bArr[ii] = true;
  }
}

// TODO: rvDoMcmc uses 0/1 int for now
// [[Rcpp::export(".QuantReg_post_adapt")]]
Rcpp::List QuantReg_post_adapt(Eigen::VectorXd y, Eigen::MatrixXd X, 
                         double alpha, int nsamples, int nburn, 
                         Eigen::VectorXd betaInit, Eigen::VectorXd mwgSd, 
                         Eigen::VectorXd rvDoMcmc,  
                         int maxIter = 100, double relTol = 1e-7) {
  InnerEL<QuantRegModel> QR;
  QR.setData(y,X,&alpha); 
  QR.setTol(maxIter, relTol);
  Eigen::VectorXd paccept;
  bool *rvdomcmc = new bool[betaInit.size()];
  dVec_to_bArr(rvDoMcmc, rvdomcmc);
  // // BEGIN DEBUG
  // std::cout << "betaInit = " << betaInit.transpose() << std::endl;
  // for (int ii=0; ii<betaInit.size(); ii++) {
  //   std::cout << "rvdomcmc[" << ii << "] = " << rvdomcmc[ii] << std::endl;
  // }
  // for (int ii=0; ii<betaInit.size(); ii++) {
  //   std::cout << "mwgSd.data()[" << ii << "] = " << mwgSd.data()[ii] << std::endl;
  // }
  // // END DEBUG
  // Eigen::MatrixXd beta_chain = QR.postSample(nsamples, nburn,
  //                                            betaInit, mwgSd, paccept);
  Eigen::MatrixXd beta_chain = QR.postSampleAdapt(nsamples, nburn,
                                                  betaInit, mwgSd.data(),
                                                  rvdomcmc, paccept);
  delete[] rvdomcmc; 
  Rcpp::List retlst;
  retlst["beta_chain"] = beta_chain;
  retlst["paccept"] = paccept;
  return(retlst);
}

// [[Rcpp::export(".QuantReg_post")]]
Rcpp::List QuantReg_post(Eigen::VectorXd y, Eigen::MatrixXd X, 
                             double alpha, int nsamples, int nburn, 
                             Eigen::VectorXd betaInit, Eigen::VectorXd sigs, 
                             int maxIter = 100, double relTol = 1e-7) {
  InnerEL<QuantRegModel> QR;
  QR.setData(y,X,&alpha); 
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  QR.setTol(maxIter, relTol);
  Eigen::VectorXd paccept; 
  Eigen::MatrixXd beta_chain = QR.postSample(nsamples, nburn, 
                                             betaInit, sigs, paccept);
  // return(beta_chain);
  Rcpp::List retlst; 
  retlst["beta_chain"] = beta_chain;
  retlst["paccept"] = paccept; 
  return(retlst); 
}
