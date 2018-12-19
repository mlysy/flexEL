#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "AdjG.h"

// [[Rcpp::export(".adjG")]]
Eigen::MatrixXd adjG(Eigen::MatrixXd G, double a) {
  return adj_G(G,a);
}
