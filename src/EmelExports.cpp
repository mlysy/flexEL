// port some of InnerELC to R

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerELC.h"
#include "MeanRegModel.h"

// Note: y, X are not actually needed here but instantiating an InnerEL object needs them 
// G: m x N matrix
// lambda0: m-vector of starting values
// [[Rcpp::export(".EMEL")]]
Rcpp::List EMEL(Eigen::MatrixXd G, Eigen::VectorXd delta, Eigen::VectorXd ws0, 
                     int maxIter, double eps, bool verbose) {
    int nObs = G.cols();
    int nEqs = G.rows();
    VectorXd y = VectorXd::Zero(nObs);
    MatrixXd X = MatrixXd::Zero(nEqs, nObs);
    InnerELC<MeanRegModel> IL(y, X, delta, NULL); // instantiate
    IL.G = G; // assign a given G
    IL.ws = ws0; 
    // initialize variables for output here 
    int nIter;
    double maxErr;
    bool not_conv;
    IL.EMEL(nIter, maxErr, maxIter, eps);
    VectorXd wsnew = IL.ws; // output
    // check convergence
    not_conv = (nIter == maxIter) && (maxErr > eps);
    if(verbose) {
        Rprintf("nIter = %i, maxErr = %f\n", nIter, maxErr);
    }
    return Rcpp::List::create(_["ws"] = wsnew,
                              _["convergence"] = !not_conv);
}
