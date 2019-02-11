/**
 * @file InnerElcExports.cpp
 * 
 * @brief Export InnerELC functions to R.
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "InnerELC.h"
#include "MeanRegModel.h"

//[[Rcpp::depends("RcppEigen")]]

using namespace Rcpp;
using namespace Eigen;

/*
// Note: for testing purpose only
// [[Rcpp::export(".evalEpsilonsLS")]]
Eigen::VectorXd evalEpsilonsLS(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z,
                               Eigen::VectorXd beta, Eigen::VectorXd gamma, double sig2) {
  // InnerELC<MeanRegModel> ILC;
  int nObs = G.cols();
  int nEqs = G.rows();
  InnerELC<MeanRegModel> ILC(nObs,nEqs);
  Eigen::VectorXd deltas = VectorXd::Zero(y.size()).array()+1.0;
  // ILC.setData(y,X,Z,deltas,NULL);
  ILC.setData(y,X,Z,deltas);
  ILC.evalEpsilons(beta,gamma,sig2);
  Eigen::VectorXd epsilons = ILC.getEpsilons();
  return(epsilons);
}
*/

// returns weights for the weighted maximum log EL
// Note: for testing purpose only
// [[Rcpp::export(".evalWeights")]]
Eigen::VectorXd evalWeights(Eigen::VectorXd deltas, Eigen::VectorXd omegas, 
                           Eigen::VectorXd epsilons, bool support) {
  int nObs = epsilons.size(); // TODO: this is problematic -- which nObs to use
  el::InnerELC<MeanRegModel> ILC(nObs,1);
  ILC.setOpts(support);
  ILC.setDeltas(deltas);
  ILC.setOmegas(omegas);
  ILC.setEpsilons(epsilons); 
  ILC.evalWeights(); 
  Eigen::VectorXd weights = ILC.getWeights(); 
  return(weights);
}

// [[Rcpp::export(".lambdaNRC")]]
Eigen::VectorXd lambdaNRC(Eigen::MatrixXd G, Eigen::VectorXd weights, 
                          int maxIter, double relTol, bool support, bool verbose) {
  int nObs = G.cols();
  int nEqs = G.rows();
  el::InnerELC<MeanRegModel> ILC(nObs,nEqs);
  ILC.setOpts(maxIter, relTol, support);
  ILC.setG(G); // assign a given G
  ILC.setWeights(weights);
  // ILC.setTol(maxIter, relTol);
  
  // initialize variables for output here 
  int nIter;
  double maxErr;
  ILC.lambdaNR(nIter, maxErr);
  VectorXd lambda = ILC.getLambda(); // output

  // check convergence
  bool not_conv = (nIter == maxIter) && (maxErr > relTol);
  if(verbose) {
    Rprintf("nIter = %i, maxErr = %f\n", nIter, maxErr);
  }

  // fill in NaN if not converged
  if (not_conv) {
    for (int ii=0; ii<lambda.size(); ii++) {
      lambda(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  return lambda;
  
  // return the status and value
  // return Rcpp::List::create(_["lambda"] = lambda,
  //                           _["convergence"] = !not_conv);
}

// G: m x N matrix
// lambda0: m-vector of starting values
// [[Rcpp::export(".omega.hat.EM")]]
Eigen::VectorXd omegaHatEM(Eigen::VectorXd omegasInit, 
                           Eigen::MatrixXd G, Eigen::VectorXd deltas,
                           Eigen::VectorXd epsilons, 
                           int maxIter, double relTol, double absTol, bool support, bool verbose) {
  int nObs = G.cols();
  int nEqs = G.rows();
  el::InnerELC<MeanRegModel> ILC(nObs,nEqs); 
  ILC.setOpts(maxIter, relTol, absTol, support);
  ILC.setDeltas(deltas);
  ILC.setG(G); // assign a given G
  ILC.setEpsilons(epsilons);
  // ILC.setTol(maxIter, relTol, absTol);
  ILC.setOmegas(omegasInit); // set initial omegas from uncensored omega.hat
  ILC.evalOmegas();
  VectorXd omegasnew = ILC.getOmegas(); // output
  return omegasnew;
  // return Eigen::VectorXd::Zero(nObs);
}

// [[Rcpp::export(".logELC")]]
double logELC(Eigen::VectorXd omegas, Eigen::VectorXd epsilons, 
              Eigen::VectorXd deltas, bool support) {
  int nObs = omegas.size();
  el::InnerELC<MeanRegModel> ILC(nObs,1);
  ILC.setOpts(support);
  ILC.setDeltas(deltas);
  ILC.setEpsilons(epsilons); 
  ILC.setOmegas(omegas);
  double logel = ILC.logEL();
  return logel; 
}

// // [[Rcpp::export(".evalPsos.smooth")]]
// double evalPsosSmooth(int ii, Eigen::VectorXd omegas, 
//                       Eigen::VectorXd epsilons, double s, bool support) {
//   int nObs = epsilons.size();
//   el::InnerELC<MeanRegModel> ILC(nObs,1);
//   ILC.setOpts(support);
//   ILC.setEpsilons(epsilons);
//   ILC.setOmegas(omegas);
//   return ILC.evalPsosSmooth(ii-1,s); // ii is 1 larger in R than in C++
// }

// [[Rcpp::export(".logEL.smooth")]]
double logELSmooth(Eigen::VectorXd omegas, 
                   Eigen::VectorXd epsilons, 
                   Eigen::VectorXd deltas, double s, bool support) {
  int nObs = omegas.size();
  el::InnerELC<MeanRegModel> ILC(nObs,1);
  ILC.setOpts(support);
  ILC.setOmegas(omegas);
  ILC.setEpsilons(epsilons);
  ILC.setDeltas(deltas);
  return ILC.logELSmooth(s);
}

// [[Rcpp::export(".evalWeights.smooth")]]
Eigen::VectorXd evalWeightsSmooth(Eigen::VectorXd deltas, 
                                  Eigen::VectorXd omegas, 
                                  Eigen::VectorXd epsilons, double s, bool support) {
  int nObs = epsilons.size();
  el::InnerELC<MeanRegModel> ILC(nObs,1);
  ILC.setOpts(support);
  ILC.setOmegas(omegas);
  ILC.setEpsilons(epsilons);
  ILC.setDeltas(deltas);
  ILC.evalWeightsSmooth(s);
  Eigen::VectorXd weights = ILC.getWeights(); 
  return(weights);
}

// [[Rcpp::export(".omega.hat.EM.smooth")]]
Eigen::VectorXd omegaHatEMSmooth(Eigen::VectorXd omegasInit,
                                 Eigen::MatrixXd G, Eigen::VectorXd deltas,
                                 Eigen::VectorXd epsilons, double s, 
                                 int maxIter, double relTol, double absTol, 
                                 bool support, bool verbose) {
  int nObs = G.cols();
  int nEqs = G.rows();
  el::InnerELC<MeanRegModel> ILC(nObs,nEqs); 
  ILC.setOpts(maxIter, relTol, absTol, support);
  ILC.setDeltas(deltas);
  ILC.setG(G); // assign a given G
  ILC.setEpsilons(epsilons); 
  // ILC.setTol(maxIter, relTol, absTol);
  ILC.setOmegas(omegasInit); // set initial omegas from uncensored omega.hat
  ILC.evalOmegasSmooth(s);
  VectorXd omegasnew = ILC.getOmegas(); // output
  return omegasnew; 
}
