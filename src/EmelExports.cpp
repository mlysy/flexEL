// port some of InnerELC to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerELC.h"
#include "MeanRegModel.h"

// returns weights for the weighted NR algorithm 
// [[Rcpp::export(".evalWeights")]]
Eigen::VectorXd evalWeights(Eigen::VectorXd y, Eigen::MatrixXd X, 
                           Eigen::VectorXd deltas, Eigen::VectorXd omegas, 
                           Eigen::VectorXd beta) {
    InnerELC<MeanRegModel> IL(y, X, deltas, NULL); // instantiate
    IL.omegas = omegas;
    IL.evalWeights(beta); 
    Eigen::VectorXd weights = IL.weights; 
    return(weights);
}


// Note: y, X are not actually needed here but instantiating an InnerEL object needs them 
// G: m x N matrix
// lambda0: m-vector of starting values
// [[Rcpp::export(".EMEL")]]
Rcpp::List EMEL(Eigen::VectorXd beta, Eigen::MatrixXd G, 
                Eigen::VectorXd deltas,Eigen::VectorXd omegas0, 
                int maxIter, double eps, bool verbose) {
    int nObs = G.cols();
    int nEqs = G.rows();
    VectorXd y = VectorXd::Zero(nObs);
    MatrixXd X = MatrixXd::Zero(nEqs, nObs);
    InnerELC<MeanRegModel> IL(y, X, deltas, NULL); // instantiate
    IL.G = G; // assign a given G
    IL.omegas = omegas0; 
    // initialize variables for output here 
    int nIter;
    double maxErr;
    bool not_conv;
    IL.EMEL(beta, nIter, maxErr, maxIter, eps);
    VectorXd omegasnew = IL.omegas; // output
    // check convergence
    not_conv = (nIter == maxIter) && (maxErr > eps);
    if(verbose) {
        Rprintf("nIter = %i, maxErr = %f\n", nIter, maxErr);
    }
    return Rcpp::List::create(_["omegas"] = omegasnew,
                              _["convergence"] = !not_conv);
}
