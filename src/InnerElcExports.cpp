// port some of InnerELC to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerELC.h"
#include "MeanRegModel.h"

// [[Rcpp::export(".evalEpsilonsLS")]]
Eigen::VectorXd evalEpsilonsLS(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, 
                               Eigen::VectorXd beta, Eigen::VectorXd gamma, double sig2) {
  InnerELC<MeanRegModel> ILC; 
  Eigen::VectorXd deltas = VectorXd::Zero(y.size()).array()+1.0;
  ILC.setData(y,X,Z,deltas,NULL); 
  ILC.evalEpsilons(beta,gamma,sig2);
  Eigen::VectorXd epsilons = ILC.getEpsilons();
  return(epsilons);
}

// TODO: getWeights does not need y X, only depends on deltas, epsilons and omegas, 
//       remove y X and use pesudo input to init
// returns weights for the weighted NR algorithm 
// [[Rcpp::export(".evalWeights")]]
Eigen::VectorXd evalWeights(Eigen::VectorXd deltas, Eigen::VectorXd omegas, 
                           Eigen::VectorXd epsilons) {
    // TODO: pseudo-input, actually can have setG to allocate the space but do this for now 
    int nObs = omegas.size();
    VectorXd y = VectorXd::Zero(nObs);
    MatrixXd X = MatrixXd::Zero(1,nObs);
    // InnerELC<MeanRegModel> ILC(y, X, deltas, NULL); // instantiate
    InnerELC<MeanRegModel> ILC; 
    ILC.setData(y,X,deltas,NULL); 
    ILC.setOmegas(omegas);
    ILC.setEpsilons(epsilons); 
    ILC.evalWeights(); 
    Eigen::VectorXd weights = ILC.getWeights(); 
    return(weights);
}
// Eigen::VectorXd evalWeights(Eigen::VectorXd y, Eigen::MatrixXd X, 
//                            Eigen::VectorXd deltas, Eigen::VectorXd omegas, 
//                            Eigen::VectorXd epsilons) {
//     InnerELC<MeanRegModel> ILC(y, X, deltas, NULL); // instantiate
//     ILC.omegas = omegas;
//     ILC.setEpsilons(epsilons); 
//     ILC.evalWeights(); 
//     Eigen::VectorXd weights = ILC.weights; 
//     return(weights);
// }

// Note: y, X are not actually needed here but instantiating an InnerEL object needs them 
// G: m x N matrix
// lambda0: m-vector of starting values
// [[Rcpp::export(".lambdaNRC")]]
Eigen::VectorXd lambdaNRC(Eigen::MatrixXd G, Eigen::VectorXd weights, 
                          int maxIter, double relTol, bool verbose) {
  // TODO: pseudo-input, actually can have setG to allocate the space but do this for now 
  int nObs = G.cols();
  int nEqs = G.rows();
  VectorXd y = VectorXd::Zero(nObs);
  MatrixXd X = MatrixXd::Zero(nEqs, nObs);
  VectorXd deltas = VectorXd::Zero(nObs);
  // InnerELC<MeanRegModel> ILC(y, X, delta, NULL); // instantiate
  InnerELC<MeanRegModel> ILC; 
  ILC.setData(y,X,deltas,NULL); 
  ILC.setG(G); // assign a given G
  ILC.setWeights(weights); 
  // initialize variables for output here 
  int nIter;
  double maxErr;
  VectorXd lambda;
  bool not_conv;
  ILC.setTol(maxIter, relTol, 0.0);
  ILC.lambdaNR(nIter, maxErr);
  lambda = ILC.getLambda(); // output
  // check convergence
  not_conv = (nIter == maxIter) && (maxErr > relTol);
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
}

/* return a list -- value & convergence status
Rcpp::List lambdaNRC(Eigen::MatrixXd G, Eigen::VectorXd weights, 
                     int maxIter, double relTol, bool verbose) {
    // TODO: pseudo-input, actually can have setG to allocate the space but do this for now 
    int nObs = G.cols();
    int nEqs = G.rows();
    VectorXd y = VectorXd::Zero(nObs);
    MatrixXd X = MatrixXd::Zero(nEqs, nObs);
    VectorXd deltas = VectorXd::Zero(nObs);
    // InnerELC<MeanRegModel> ILC(y, X, delta, NULL); // instantiate
    InnerELC<MeanRegModel> ILC; 
    ILC.setData(y,X,deltas,NULL); 
    ILC.setG(G); // assign a given G
    ILC.setWeights(weights); 
    // initialize variables for output here 
    int nIter;
    double maxErr;
    VectorXd lambda;
    bool not_conv;
    ILC.lambdaNR(nIter, maxErr, maxIter, relTol);
    lambda = ILC.getLambda(); // output
    // check convergence
    not_conv = (nIter == maxIter) && (maxErr > relTol);
    if(verbose) {
        Rprintf("nIter = %i, maxErr = %f\n", nIter, maxErr);
    }
    return Rcpp::List::create(_["lambda"] = lambda,
                              _["convergence"] = !not_conv);
}
*/

// G: m x N matrix
// lambda0: m-vector of starting values
// [[Rcpp::export(".omega.hat.EM")]]
Eigen::VectorXd omegaHatEM(Eigen::VectorXd omegasInit, 
                           Eigen::MatrixXd G, Eigen::VectorXd deltas,
                           Eigen::VectorXd epsilons, 
                           int maxIter, double relTol, double absTol, bool verbose) {
    // TODO: pseudo-input, actually can have setG to allocate the space but do this for now 
    int nObs = G.cols();
    int nEqs = G.rows();
    VectorXd y = VectorXd::Zero(nObs);
    MatrixXd X = MatrixXd::Zero(nEqs, nObs);
    // InnerELC<MeanRegModel> ILC(y, X, deltas, NULL); // instantiate
    InnerELC<MeanRegModel> ILC; 
    ILC.setData(y,X,deltas,NULL); 
    ILC.setG(G); // assign a given G
    ILC.setEpsilons(epsilons); 
    ILC.setTol(maxIter, relTol, absTol); // TODO: change this
    // Note: set initial omegas from uncensored omega.hat
    ILC.setOmegas(omegasInit);
    // std::cout << "omegasInit = " << omegasInit.transpose() << std::endl;
    ILC.evalOmegas();
    VectorXd omegasnew = ILC.getOmegas(); // output
    return omegasnew; 
}

// [[Rcpp::export(".logELC")]]
double logELC(Eigen::VectorXd omegas, Eigen::VectorXd epsilons, 
              Eigen::VectorXd deltas) {
  int nObs = omegas.size();
  VectorXd y = VectorXd::Zero(nObs);
  MatrixXd X = MatrixXd::Zero(1, nObs);
  InnerELC<MeanRegModel> ILC;
  ILC.setData(y,X,deltas,NULL);
  ILC.setEpsilons(epsilons); 
  // set omegas which is obtained from EM algorithm 
  ILC.setOmegas(omegas);
  double logel = ILC.logEL();
  return logel; 
}

// [[Rcpp::export(".ind.smooth")]]
Eigen::VectorXd indSmooth(Eigen::VectorXd x, Eigen::VectorXd s) {
  InnerELC<MeanRegModel> ILC;
  return ILC.indSmooth(x,s);
}

// [[Rcpp::export(".evalPsos.Smooth")]]
double evalPsosSmooth(int ii, Eigen::VectorXd omegas, 
                      Eigen::VectorXd epsilons, double s) {
  InnerELC<MeanRegModel> ILC;
  ILC.setEpsilons(epsilons);
  ILC.setOmegas(omegas);
  return ILC.evalPsosSmooth(ii,s);
}

// [[Rcpp::export(".logEL.smooth")]]
double logELSmooth(Eigen::VectorXd omegas, 
                   Eigen::VectorXd epsilons, 
                   Eigen::VectorXd deltas, double s) {
  InnerELC<MeanRegModel> ILC;
  ILC.setOmegas(omegas);
  ILC.setEpsilons(epsilons);
  ILC.setDeltas(deltas);
  return ILC.logELSmooth(s);
}

// [[Rcpp::export(".evalWeights.smooth")]]
Eigen::VectorXd evalWeightsSmooth(Eigen::VectorXd deltas, 
                                  Eigen::VectorXd omegas, 
                                  Eigen::VectorXd epsilons, double s) {
  InnerELC<MeanRegModel> ILC;
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
                                 bool verbose) {
  // TODO: pseudo-input, actually can have setG to allocate the space but do this for now 
  int nObs = G.cols();
  int nEqs = G.rows();
  VectorXd y = VectorXd::Zero(nObs);
  MatrixXd X = MatrixXd::Zero(nEqs, nObs);
  // InnerELC<MeanRegModel> ILC(y, X, deltas, NULL); // instantiate
  InnerELC<MeanRegModel> ILC; 
  ILC.setData(y,X,deltas,NULL); 
  ILC.setG(G); // assign a given G
  ILC.setEpsilons(epsilons); 
  ILC.setTol(maxIter, relTol, absTol);
  // Note: set initial omegas from uncensored omega.hat
  ILC.setOmegas(omegasInit);
  // std::cout << "omegasInit = " << omegasInit.transpose() << std::endl;
  ILC.evalOmegasSmooth(s);
  VectorXd omegasnew = ILC.getOmegas(); // output
  return omegasnew; 
}

// // Returns the maximized log empirical likelihood given G, 
// //      or -Inf if no omegas satisfies G.
// // Note: omegas should be the returned value from uncensored omega.hat
// // [[Rcpp::export(".logELC")]]
// double logELC(Eigen::VectorXd omegas, Eigen::MatrixXd G, 
//               Eigen::VectorXd deltas, Eigen::VectorXd epsilons) {
//     int nObs = G.cols();
//     int nEqs = G.rows();
//     VectorXd y = VectorXd::Zero(nObs);
//     MatrixXd X = MatrixXd::Zero(nEqs, nObs);
//     InnerELC<MeanRegModel> ILC;
//     ILC.setData(y,X,deltas,NULL);
//     ILC.setG(G); // set the given G
//     ILC.setEpsilons(epsilons); 
//     // set omegas which is obtained from EM algorithm 
//     ILC.setOmegas(omegas);
//     double logel = ILC.logEL();
//     return logel; 
// }
