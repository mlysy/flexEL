#ifndef POSTSAMPLE_h
#define POSTSAMPLE_h

// class for posterior sampling 

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
#include "mwgAdapt.h"
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

template <typename elOpt>
class PostSample { 
private:
    template <typename elModel>
    elOpt<elModel> EL; 
public:
    template <typename elModel>
    PostSample(elOpt<elModel> EL);
    MatrixXd Sample(int nsample, int nburn, const Ref<const VectorXd>& thetaInit, 
                    const Ref<const VectorXd>& sigsInit); 
}

template<typename elOpt>
template<typename elModel>
inline PostSample<elOpt>::do_mh(elOpt<elModel> EL) {
    EL = EL;
}

void do_mh(VectorXd &theta_cur, double &sig, bool &acc) {
    return;
}

void store_theta(VectorXd &theta_cur) {
    return;
}

void Sample(int nsample, int nburn, VectorXd thetaInit, 
            const Ref<const VectorXd>& sigsInit) {
    int nrv = thetaInit.rows(); 
    Eigen::VectorXd mwgSd(nrv);
    for (int ii=0; ii<nrv; ii++) {
        mwgSd(ii) = R::unif_rand(); // starting values of random walk jump sizes
    }
    Eigen::VectorXd theta_cur = thetaInit;
    // default constructor: all parameters are updated
    mwgAdapt tuneMCMC(nrv);
    // debug mode constructor: some of the parameters are fixed to investigate
    // MCMC convergence issues (using fewer parameters)
    // do_rv_mcmc is a bool (or int) vector of true/false (or 0/1)
    // mwgAdapt tuneMCMC(number_of_rvs, do_rv_mcmc);
    // do_rv_mcmc was passed to C++ from R as type
    // Rcpp::LogicalVector (or Rcpp::IntegerVector.  
    // R always stores bool arrays internally as ints,
    // because in R the array could contain NA, 
    //  but in C++ bool[] can only contain T/F)
    // mwgAdapt tuneMCMC(number_of_rvs, LOGICAL(do_rv_mcmc));
    
    // MCMC loop
    for(int ii=-nburn; ii<nsample; ii++) {
        for(int jj=0; jj<nrv; jj++) {
            if(do_rv_mcmc(jj)) {
                // single component mwg update,
                // takes current state of MCMC and changes only one component
                // isAccepted is boolean for whether or not proposal was accepted
                theta_curr(jj) = do_mh(theta_curr, mwgSd(jj), isAccepted(jj));
            }
        }
        // R equivalent of Theta_out[ii,] <- theta_curr
        store_theta(theta_chain, theta_curr, ii);
        tuneMCMC.adapt(mwgSd, isAccepted);
    }
    return;
}