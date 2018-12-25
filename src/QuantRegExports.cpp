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

// double tau
// [[Rcpp::export(".QuantReg_evalG")]]
Eigen::MatrixXd QuantReg_evalG(Eigen::VectorXd y, Eigen::MatrixXd X,
                               Eigen::VectorXd tauArr, Eigen::VectorXd beta) {
  // InnerEL<QuantRegModel> QR(y, X, &tau); // instantiate
  el::InnerEL<QuantRegModel> QR;
  // QR.setData(y,X,&tau); 
  QR.setData(y,X,tauArr.data());
  QR.evalG(beta);
  Eigen::MatrixXd G = QR.getG(); // G is nEqs x nObs
  return(G); 
}

// [[Rcpp::export(".QuantRegLS_evalG")]]
Eigen::MatrixXd QuantRegLS_evalG(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, 
                                 Eigen::VectorXd tauArr, Eigen::VectorXd beta, 
                                 Eigen::VectorXd gamma, double sig2, Eigen::VectorXd Nu) {
  // InnerEL<QuantRegModel> QR(y, X, &tau); // instantiate
  el::InnerEL<QuantRegModel> QR;
  // QR.setData(y,X,Z,&tau); 
  QR.setData(y,X,Z,tauArr.data());
  QR.evalG(beta,gamma,sig2,Nu);
  Eigen::MatrixXd G = QR.getG(); // G is nEqs x nObs
  return(G); 
}

// [[Rcpp::export(".QuantRegLS_evalGSmooth")]]
Eigen::MatrixXd QuantRegLS_evalGSmooth(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, 
                                       Eigen::VectorXd tauArr, Eigen::VectorXd beta, 
                                       Eigen::VectorXd gamma, double sig2, Eigen::VectorXd Nu, double s) {
  // InnerEL<QuantRegModel> QR(y, X, &tau); // instantiate
  el::InnerEL<QuantRegModel> QR;
  // QR.setData(y,X,Z,&tau); 
  QR.setData(y,X,Z,tauArr.data());
  QR.evalGSmooth(beta,gamma,sig2,Nu,s);
  Eigen::MatrixXd G = QR.getG(); // G is nEqs x nObs
  return(G); 
}


/*
// [[Rcpp::export(".QuantReg_post")]]
Rcpp::List QuantReg_post(Eigen::VectorXd y, Eigen::MatrixXd X, 
                         Eigen::VectorXd alphaArr, 
                         int nsamples, int nburn, 
                         Eigen::MatrixXd BetaInit, Eigen::MatrixXd MwgSds, 
                         Eigen::MatrixXd RvDoMcmc, 
                         int maxIter = 100, double relTol = 1e-7) {
  InnerEL<QuantRegModel> QR;
  QR.setData(y,X,alphaArr.data());
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  QR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept; 
  Eigen::MatrixXd beta_chain = QR.postSample(nsamples, nburn, 
                                             BetaInit, MwgSds, 
                                             RvDoMcmc, paccept); 
  Rcpp::List retlst; 
  retlst["beta_chain"] = beta_chain;
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
  int nThe = betaInit.size();
  bool *rvdomcmc = new bool[nThe];
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

// Note: BetaInit and GammaInit must have the same number of columns
// [[Rcpp::export(".QuantRegLS_post")]]
Rcpp::List QuantRegLS_post(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, 
                           Eigen::VectorXd alphaArr, int nsamples, int nburn, 
                           Eigen::MatrixXd betaInit, Eigen::MatrixXd gammaInit, 
                           Eigen::VectorXd sig2Init, Eigen::MatrixXd nuInit, 
                           Eigen::MatrixXd mwgSds, Eigen::MatrixXd rvDoMcmc, 
                           int maxIter = 100, double relTol = 1e-7) {
  InnerEL<QuantRegModel> QR;
  QR.setData(y,X,Z,alphaArr.data());
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  QR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept;
  // concatenate (stack) them together
  Eigen::MatrixXd thetaInit(betaInit.rows()+gammaInit.rows()+1, betaInit.cols());
  thetaInit.topRows(betaInit.rows()) = betaInit;
  thetaInit.block(betaInit.rows(),0,gammaInit.rows(),gammaInit.cols()) = gammaInit;
  thetaInit.block(betaInit.rows()+gammaInit.rows(),0,1,gammaInit.cols()) = sig2Init;
  thetaInit.bottomRows(1) = nuInit.transpose();
  // std::cout << "QuantRegExports: thetaInit = \n" << thetaInit << std::endl;
  Eigen::MatrixXd theta_chain = QR.postSample(nsamples, nburn,
                                              thetaInit, mwgSds, 
                                              rvDoMcmc, paccept);
                                              // betaInit.rows(), gammaInit.rows());
  Rcpp::List retlst; 
  retlst["theta_chain"] = theta_chain;
  retlst["paccept"] = paccept;
  return(retlst); 
}

// Note: BetaInit and GammaInit must have the same number of columns
// [[Rcpp::export(".QuantRegLS_post_adapt")]]
Rcpp::List QuantRegLS_post_adapt(Eigen::VectorXd y, Eigen::MatrixXd X,
                                 Eigen::MatrixXd Z, Eigen::VectorXd alphaArr, 
                                 int nsamples, int nburn,
                                 Eigen::VectorXd betaInit, Eigen::VectorXd gammaInit, 
                                 Eigen::VectorXd sig2Init, Eigen::VectorXd nuInit,
                                 Eigen::VectorXd mwgSd, Eigen::VectorXd rvDoMcmc,
                                 int maxIter = 100, double relTol = 1e-7) {
  InnerEL<QuantRegModel> QR;
  QR.setData(y,X,Z,alphaArr.data());
  QR.setTol(maxIter, relTol);
  Eigen::VectorXd paccept;
  
  // concatenate (stack) them together
  Eigen::VectorXd thetaInit(betaInit.size()+gammaInit.size()+2, betaInit.size());
  thetaInit.head(betaInit.size()) = betaInit;
  thetaInit.segment(betaInit.size(),gammaInit.size()) = gammaInit;
  thetaInit.segment(betaInit.size()+gammaInit.size(),1) = sig2Init;
  thetaInit.tail(1) = nuInit;
  // std::cout << "thetaInit = " << thetaInit.transpose() << std::endl;
  
  int nThe = thetaInit.size();
  bool *rvdomcmc = new bool[nThe];
  dVec_to_bArr(rvDoMcmc, rvdomcmc);
  Eigen::MatrixXd theta_chain = QR.postSampleAdapt(nsamples, nburn,
                                                   thetaInit, mwgSd.data(),
                                                   rvdomcmc, paccept);
  // Eigen::MatrixXd theta_chain = Eigen::MatrixXd::Zero(1,1);
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
                             int maxIter = 100, double relTol = 1e-7, double absTol = 1e-3) {
  InnerELC<QuantRegModel> QRC;
  QRC.setData(y,X,deltas,alphaArr.data());
  QRC.setTol(maxIter, relTol, absTol);
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

// [[Rcpp::export(".QuantRegCens_post_adapt")]]
Rcpp::List QuantRegCens_post_adapt(Eigen::VectorXd omegasInit, 
                                   Eigen::VectorXd y, Eigen::MatrixXd X,
                                   Eigen::VectorXd deltas, Eigen::VectorXd alphaArr,
                                   int nsamples, int nburn,
                                   Eigen::VectorXd betaInit, Eigen::VectorXd mwgSd,
                                   Eigen::VectorXd rvDoMcmc, Eigen::VectorXd doAdapt,
                                   int maxIter = 100, double relTol = 1e-7, double absTol = 1e-3) {
  InnerELC<QuantRegModel> QRC;
  QRC.setData(y,X,deltas,alphaArr.data());
  QRC.setTol(maxIter, relTol, absTol);
  int nTheta = betaInit.size();
  Eigen::VectorXd paccept(nTheta);
  bool *doadapt = new bool[nTheta];
  dVec_to_bArr(doAdapt, doadapt);
  QRC.setOmegas(omegasInit); // set initial value for first EM
  Eigen::MatrixXd beta_chain = QRC.postSampleAdapt(nsamples, nburn,
                                                   betaInit, mwgSd.data(),
                                                   rvDoMcmc, doadapt, paccept);
  // Eigen::MatrixXd beta_chain = Eigen::MatrixXd::Zero(1,1);
  delete[] doadapt;
  Rcpp::List retlst;
  retlst["beta_chain"] = beta_chain;
  retlst["paccept"] = paccept;
  return(retlst);
}

// [[Rcpp::export(".QuantRegCensLS_post_adapt")]]
Rcpp::List QuantRegCensLS_post_adapt(Eigen::VectorXd omegasInit, 
                                     Eigen::VectorXd y, Eigen::MatrixXd X,
                                     Eigen::MatrixXd Z, Eigen::VectorXd deltas, 
                                     Eigen::VectorXd alphaArr,
                                     int nsamples, int nburn,
                                     Eigen::VectorXd betaInit, 
                                     Eigen::VectorXd gammaInit, 
                                     double sig2Init,
                                     Eigen::VectorXd nuInit, Eigen::VectorXd mwgSd, 
                                     Eigen::VectorXd rvDoMcmc, Eigen::VectorXd doAdapt, 
                                     int maxIter = 100, double relTol = 1e-7, double absTol = 1e-3) {
  InnerELC<QuantRegModel> QRC;
  QRC.setData(y,X,Z,deltas,alphaArr.data());
  QRC.setTol(maxIter, relTol, absTol);
  QRC.setOmegas(omegasInit); // set initial value for first EM
  Eigen::VectorXd paccept;
  
  // concatenate (stack) them together
  Eigen::VectorXd thetaInit(betaInit.size()+gammaInit.size()+2, betaInit.size());
  thetaInit.head(betaInit.size()) = betaInit;
  thetaInit.segment(betaInit.size(),gammaInit.size()) = gammaInit;
  thetaInit.segment(betaInit.size()+gammaInit.size(),1)(0) = sig2Init;
  // thetaInit.tail(1) = nuInit;
  // TODO: length of nuInit should be equal to alphaArr(0)
  thetaInit.tail(alphaArr(0)) = nuInit; 
  // std::cout << thetaInit << std::endl;
  
  int nTheta = thetaInit.size();
  bool *doadapt = new bool[nTheta];
  dVec_to_bArr(doAdapt, doadapt);
  Eigen::MatrixXd theta_chain = QRC.postSampleAdapt(nsamples, nburn,
                                                    thetaInit, mwgSd.data(),
                                                    rvDoMcmc, doadapt, paccept);
  // Eigen::MatrixXd theta_chain = Eigen::MatrixXd::Zero(1,1);
  delete[] doadapt;
  Rcpp::List retlst;
  retlst["theta_chain"] = theta_chain;
  retlst["paccept"] = paccept;
  return(retlst);
}
*/

// [[Rcpp::export(".rho1.smooth")]]
double rho1Smooth(double u, double tau, double s) {
  el::InnerELC<QuantRegModel> QRC;
  return(QRC.phi_tau_smooth(u,tau,s));
}
