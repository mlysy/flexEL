#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "IndSmooth.h"

// // [[Rcpp::export(".ind.smooth")]]
// Eigen::VectorXd indSmooth(Eigen::VectorXd x, Eigen::VectorXd s) {
//   return ind_smooth(x,s);
// }
// 
// // [[Rcpp::export(".ind1.smooth")]]
// Eigen::VectorXd ind1Smooth(Eigen::VectorXd x, Eigen::VectorXd s) {
//   return ind1_smooth(x,s);
// }

// [[Rcpp::export(".ind.smooth")]]
double indSmooth(double x, double s) {
  return ind_smooth(x,s);
}

// [[Rcpp::export(".ind1.smooth")]]
double ind1Smooth(double x, double s) {
  return ind1_smooth(x,s);
}