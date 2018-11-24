// port some of InnerEL to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
#include "MeanRegModel.h"

// [[Rcpp::export(".lambdaNR")]]
Eigen::VectorXd lambdaNR(Eigen::MatrixXd G, int maxIter, double relTol, bool verbose) {
  int nObs = G.cols();
  int nEqs = G.rows();
  InnerEL<MeanRegModel> IL(nObs,nEqs);
  IL.setG(G); // assign the given G
  
  // initialize variables for output here 
  int nIter;
  double maxErr;
  bool not_conv;
  IL.setTol(maxIter, relTol);
  IL.lambdaNR(nIter, maxErr);
  VectorXd lambda = IL.getLambda(); // output (could be not converged)
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
  
  // return the status and value
  // return Rcpp::List::create(_["lambda"] = lambda,
  //                           _["convergence"] = !not_conv);
}

// Eigen::VectorXd y, Eigen::MatrixXd X not needed here since there is no ordering 
// as in the censored case 
// [[Rcpp::export(".omega.hat")]]
Eigen::VectorXd omegaHat(Eigen::MatrixXd G, Eigen::VectorXd lambda) {
  // TODO: pseudo-input, actually can have setG to allocate the space but do this for now 
  int nObs = G.cols();
  int nEqs = G.rows();
  InnerEL<MeanRegModel> IL(nObs,nEqs);
  IL.setG(G); // assign the given G
  IL.setLambda(lambda); 
  IL.evalOmegas(); // calculate omegas
  VectorXd omegasnew = IL.getOmegas(); // get omegas
  return omegasnew;
}

// Calculate logEL given omegas
// [[Rcpp::export(".logEL")]]
double logEL(Eigen::VectorXd omegas) {
  InnerEL<MeanRegModel> IL;
  IL.setOmegas(omegas); 
  double logel = IL.logEL();
  return logel; 
}

