// port some of InnerEL to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
#include "MeanRegModel.h"

/**
 * @brief      simulate one draw x ~ p(x | theta)
 *
 * @param[in]  G        A matrix of dimension \code{nEqs} x \code{nObs}
 * @param[in]  maxIter  Maximum number of iterations allowed
 * @param[in]  relTol   Relative tolerance for convergence
 * @param[in]  verbose  Whether to print number of iteration and maximum error of the return value
 */
// [[Rcpp::export(".lambdaNR")]]
Eigen::VectorXd lambdaNR(Eigen::MatrixXd G, int maxIter, double relTol, bool verbose) {
  int nObs = G.cols();
  int nEqs = G.rows();
  // // Note: y, X are not actually needed here but instantiating an InnerEL object needs them
  // VectorXd y = VectorXd::Zero(nObs);
  // MatrixXd X = MatrixXd::Zero(nEqs, nObs);
  // InnerEL<MeanRegModel> IL(y, X, NULL); // instantiate
  // InnerEL<MeanRegModel> IL; // Here init with an MeanRegModel, but it is the same using QuantRegModel
  // IL.setData(y, X, NULL);
  InnerEL<MeanRegModel> IL(nObs,nEqs);
  IL.setG(G); // assign the given G
  // initialize variables for output here 
  int nIter;
  double maxErr;
  // VectorXd lambda;
  bool not_conv;
  IL.setTol(maxIter, relTol); // New
  IL.lambdaNR(nIter, maxErr);
  // IL.lambdaNR(nIter, maxErr, maxIter, relTol);
  // lambda = IL.getLambda(); // output (could be not converged)
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
}

/* Old: return a list -- value & convergence status
Rcpp::List lambdaNR(Eigen::MatrixXd G, 
                    int maxIter, double relTol, bool verbose) {
    // TODO: pseudo-input, actually can have setG to allocate the space but do this for now 
    int nObs = G.cols();
    int nEqs = G.rows();
    VectorXd y = VectorXd::Zero(nObs);
    MatrixXd X = MatrixXd::Zero(nEqs, nObs);
    // Here init with an MeanRegModel, but it is the same using QuantRegModel
    // InnerEL<MeanRegModel> IL(y, X, NULL); // instantiate
    InnerEL<MeanRegModel> IL;
    IL.setData(y, X, NULL);
    IL.setG(G); // assign the given G
    // initialize variables for output here 
    int nIter;
    double maxErr;
    VectorXd lambda;
    bool not_conv;
    IL.lambdaNR(nIter, maxErr, maxIter, relTol);
    lambda = IL.getLambda(); // output
    // check convergence
    not_conv = (nIter == maxIter) && (maxErr > relTol);
    if(verbose) {
        Rprintf("nIter = %i, maxErr = %f\n", nIter, maxErr);
    }
    return Rcpp::List::create(_["lambda"] = lambda,
                              _["convergence"] = !not_conv);
}
*/ 

// Eigen::VectorXd y, Eigen::MatrixXd X not needed here since there is no ordering 
// as in the censored case 
// [[Rcpp::export(".omega.hat")]]
Eigen::VectorXd omegaHat(Eigen::MatrixXd G, Eigen::VectorXd lambda) {
  // TODO: pseudo-input, actually can have setG to allocate the space but do this for now 
  int nObs = G.cols();
  int nEqs = G.rows();
  // VectorXd y = VectorXd::Zero(nObs);
  // MatrixXd X = MatrixXd::Zero(nEqs, nObs);
  // Here init with an MeanRegModel, but it is the same using QuantRegModel
  // InnerEL<MeanRegModel> IL(y, X, NULL); // instantiate
  // InnerEL<MeanRegModel> IL;
  // IL.setData(y,X,NULL);
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

