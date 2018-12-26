/**
 * @file MeanRegExports.cpp
 * 
 * @brief Export MeanRegModel functions to R.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "InnerEL.h"
#include "InnerELC.h"
#include "MeanRegModel.h"
#include "dVecTobArr.h"

//[[Rcpp::depends("RcppEigen")]]

using namespace Rcpp;
using namespace Eigen;

/* ------ G matrix construction ------ */

// [[Rcpp::export(".MeanReg_evalG")]]
Eigen::MatrixXd MeanReg_evalG(Eigen::VectorXd y, Eigen::MatrixXd X, 
                              Eigen::VectorXd beta) {
  el::InnerEL<MeanRegModel> MR;
  MR.setData(y,X);
  MR.evalG(beta);
  Eigen::MatrixXd G = MR.getG(); // G is nEqs x nObs
  return(G);
}

// [[Rcpp::export(".MeanRegLS_evalG")]]
Eigen::MatrixXd MeanRegLS_evalG(Eigen::VectorXd y, 
                                Eigen::MatrixXd X, Eigen::MatrixXd Z,
                                Eigen::VectorXd beta, Eigen::VectorXd gamma, 
                                double sig2) {
  el::InnerEL<MeanRegModel> MR;
  MR.setData(y,X,Z); 
  MR.evalG(beta,gamma,sig2);
  Eigen::MatrixXd G = MR.getG(); // G is nEqs x nObs
  return(G); 
}

/* ------ posterior samplers ------ */
/*
// [[Rcpp::export(".MeanReg_post")]]
Rcpp::List MeanReg_post(Eigen::VectorXd y, Eigen::MatrixXd X, 
                        int nsamples, int nburn, 
                        Eigen::MatrixXd BetaInit, 
                        Eigen::MatrixXd MwgSd, Eigen::MatrixXd RvDoMcmc, 
                        int maxIter = 100, double relTol = 1e-7) {
  // InnerEL<QuantRegModel> QR(y, X, &alpha); // instantiate QR(nObs, nEqs, alpha, lambda0);
  InnerEL<MeanRegModel> MR;
  MR.setData(y,X,NULL); 
  MR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept; 
  Eigen::MatrixXd beta_chain = MR.postSample(nsamples, nburn, 
                                             BetaInit, MwgSd, 
                                             RvDoMcmc, paccept);
  Rcpp::List retlst; 
  retlst["beta_chain"] = beta_chain;
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
  MR.setTol(maxIter,relTol);
  Eigen::VectorXd paccept;
  int nThe = betaInit.size();
  bool *rvdomcmc = new bool[nThe];
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
                          Eigen::MatrixXd betaInit, Eigen::MatrixXd gammaInit, 
                          Eigen::VectorXd sig2Init, 
                          Eigen::MatrixXd mwgSd, Eigen::MatrixXd rvDoMcmc, 
                          int maxIter = 100, double relTol = 1e-7) {
  InnerEL<MeanRegModel> MR;
  MR.setData(y,X,Z,NULL);
  MR.setTol(maxIter, relTol);
  Eigen::MatrixXd paccept;
  
  // concatenate (stack) them together
  Eigen::MatrixXd thetaInit(betaInit.rows()+gammaInit.rows()+1, betaInit.cols());
  thetaInit.topRows(betaInit.rows()) = betaInit;
  thetaInit.block(betaInit.rows(),0,gammaInit.rows(),gammaInit.cols()) = gammaInit;
  thetaInit.bottomRows(1) = sig2Init.transpose();
  
  // std::cout << "MeanRegExports: thetaInit = " << thetaInit.transpose() << std::endl;
  Eigen::MatrixXd theta_chain = MR.postSample(nsamples, nburn,
                                              thetaInit, mwgSd,
                                              rvDoMcmc, paccept);
                                              // betaInit.rows(), gammaInit.rows());
  Rcpp::List retlst; 
  retlst["theta_chain"] = theta_chain;
  retlst["paccept"] = paccept;
  return(retlst); 
}


// Note: BetaInit and GammaInit must have the same number of columns
// [[Rcpp::export(".MeanRegLS_post_adapt")]]
Rcpp::List MeanRegLS_post_adapt(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z,
                                int nsamples, int nburn,
                                Eigen::VectorXd betaInit, Eigen::VectorXd gammaInit, 
                                double sig2Init,
                                Eigen::VectorXd mwgSd, Eigen::VectorXd rvDoMcmc,
                                int maxIter = 100, double relTol = 1e-7) {
  InnerEL<MeanRegModel> MR;
  MR.setData(y,X,Z,NULL);
  MR.setTol(maxIter, relTol);
  Eigen::VectorXd paccept;
  
  // concatenate (stack) them together
  Eigen::VectorXd thetaInit(betaInit.size()+gammaInit.size()+1, betaInit.size());
  thetaInit.head(betaInit.size()) = betaInit;
  thetaInit.segment(betaInit.size(),gammaInit.size()) = gammaInit;
  thetaInit.tail(1)(0) = sig2Init;
  
  int nThe = thetaInit.size();
  bool *rvdomcmc = new bool[nThe];
  dVec_to_bArr(rvDoMcmc, rvdomcmc);
  Eigen::MatrixXd theta_chain = MR.postSampleAdapt(nsamples, nburn,
                                                  thetaInit, mwgSd.data(),
                                                  rvdomcmc, paccept);
  delete[] rvdomcmc;
  Rcpp::List retlst;
  retlst["theta_chain"] = theta_chain;
  retlst["paccept"] = paccept;
  return(retlst);
}

// [[Rcpp::export(".MeanRegCens_post")]]
Rcpp::List MeanRegCens_post(Eigen::VectorXd omegasInit, 
                            Eigen::VectorXd y, Eigen::MatrixXd X,
                            Eigen::VectorXd deltas,
                            int nsamples, int nburn,
                            Eigen::VectorXd betaInit, Eigen::VectorXd mwgSd,
                            Eigen::VectorXd rvDoMcmc,  Eigen::VectorXd doAdapt,
                            int maxIter = 100, double relTol = 1e-7, double absTol = 1e-3) {
  InnerELC<MeanRegModel> MRC;
  MRC.setData(y,X,deltas,NULL);
  MRC.setTol(maxIter, relTol, absTol);
  Eigen::VectorXd paccept(betaInit.size());
  MRC.setOmegas(omegasInit); // set initial value for first EM
  Eigen::MatrixXd beta_chain = MRC.postSample(nsamples, nburn,
                                              betaInit, mwgSd,
                                              rvDoMcmc, paccept);
  Rcpp::List retlst;
  retlst["beta_chain"] = beta_chain;
  retlst["paccept"] = paccept;
  return(retlst);
}

// [[Rcpp::export(".MeanRegCens_post_adapt")]]
Rcpp::List MeanRegCens_post_adapt(Eigen::VectorXd omegasInit, 
                                  Eigen::VectorXd y, Eigen::MatrixXd X,
                                  Eigen::VectorXd deltas,
                                  int nsamples, int nburn,
                                  Eigen::VectorXd betaInit, Eigen::VectorXd mwgSd,
                                  Eigen::VectorXd rvDoMcmc, Eigen::VectorXd doAdapt,
                                  int maxIter = 100, double relTol = 1e-7, double absTol = 1e-3) {
  InnerELC<MeanRegModel> MRC;
  MRC.setData(y,X,deltas,NULL);
  MRC.setTol(maxIter, relTol, absTol);
  int nTheta = betaInit.size();
  Eigen::VectorXd paccept(nTheta);
  bool *doadapt = new bool[nTheta];
  dVec_to_bArr(doAdapt, doadapt);
  MRC.setOmegas(omegasInit); // set initial value for first EM
  Eigen::MatrixXd beta_chain = MRC.postSampleAdapt(nsamples, nburn,
                                                   betaInit, mwgSd.data(),
                                                   rvDoMcmc, doadapt, paccept);
  // MatrixXd beta_chain = MatrixXd::Zero(1,1);
  // std::cout << "betaInit = " << betaInit.transpose() << std::endl;
  // std::cout << "mwgSd = " << mwgSd.transpose() << std::endl;
  // std::cout << "rvDoMcmc = " << rvDoMcmc.transpose() << std::endl;
  delete[] doadapt;
  Rcpp::List retlst;
  retlst["beta_chain"] = beta_chain;
  retlst["paccept"] = paccept;
  return(retlst);
}

// Note: BetaInit and GammaInit must have the same number of columns
// [[Rcpp::export(".MeanRegCensLS_post_adapt")]]
Rcpp::List MeanRegCensLS_post_adapt(Eigen::VectorXd omegasInit, 
                                    Eigen::VectorXd y, Eigen::MatrixXd X,
                                    Eigen::MatrixXd Z, Eigen::VectorXd deltas, 
                                    int nsamples, int nburn,
                                    Eigen::VectorXd betaInit, Eigen::VectorXd gammaInit, 
                                    double sig2Init, Eigen::VectorXd mwgSd, 
                                    Eigen::VectorXd rvDoMcmc, Eigen::VectorXd doAdapt,
                                    int maxIter = 100, double relTol = 1e-7, double absTol = 1e-3) {
  InnerELC<MeanRegModel> MRC;
  MRC.setData(y,X,Z,deltas,NULL);
  MRC.setTol(maxIter, relTol, absTol);
  MRC.setOmegas(omegasInit); // set initial value for first EM
  Eigen::VectorXd paccept;
  
  // concatenate (stack) them together
  Eigen::VectorXd thetaInit(betaInit.size()+gammaInit.size()+1, betaInit.size());
  thetaInit.head(betaInit.size()) = betaInit;
  thetaInit.segment(betaInit.size(),gammaInit.size()) = gammaInit;
  thetaInit.tail(1)(0) = sig2Init;
  
  int nTheta = thetaInit.size();
  bool *doadapt = new bool[nTheta];
  dVec_to_bArr(doAdapt, doadapt);
  Eigen::MatrixXd theta_chain = MRC.postSampleAdapt(nsamples, nburn,
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