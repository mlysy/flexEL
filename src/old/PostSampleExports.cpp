// posterior sampler for mean regression with location model

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
using namespace Eigen;
#include "InnerEL.h"
#include "InnerELC.h"
#include "MeanRegModel.h"
#include "MeanRegLSModel.h"
#include "PostSample.h"


// [[Rcpp::export(".MeanRegPost")]]
Eigen::MatrixXd MeanRegPost(int nsample, int nburn, 
                            Eigen::MatrixXd &theta_chain, 
                            const Ref<const Eigen::VectorXd>& theta_init, 
                            const Ref<const Eigen::VectorXd>& sigs_init,
                            const Ref<const Rcpp::LogicalVector>& do_rv_mcmc) {
    

}

// // [[Rcpp::export(".MeanRegPostLS")]]
// Eigen::MatrixXd MeanRegPostLS() {
//     
// }