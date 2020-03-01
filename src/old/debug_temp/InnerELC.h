/**
 * @file InnerELC.h
 * 
 * @brief Inner optimization routine for empirial likelihood problems under right-censoring.
 */

#ifndef INNERELC_h
#define INNERELC_h

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <Rmath.h> // for random number 
#include <cmath>  // for abs on scalars
#include "SortOrder.h"
#include "IndSmooth.h" // for smoothed indicator function
#include "BlockOuter.h"
#include "AdjG.h" // for support correction
// #include "MwgAdapt.h" // for adaptive mcmc

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

/* --------------------------------------------------------------------------- */

/**
 * @brief EL namespace
 * 
 * Wrap the exported library components into a namespace called \b el to avoid potential naming conflicts with other libraries or user-defined headers.
 */
namespace el {
  
  /**
   * @file       InnerELC.h
   *
   * @class      InnerELC
   *
   * @brief      A template class for empirical likelihood inner optimization calculation with right-censored responses.
   */
  template <typename ELModel>
  class InnerELC : public ELModel { 
  private:
    
    // required members in ELModel
    using ELModel::nObs_;
    using ELModel::nEqs_;
    
    // placeholders for EM algorithm
    VectorXd deltas_; // censoring indicators
    VectorXd weights_; // weights in weighted maximum log EL
    VectorXd omegas_; // empirical distribution
    VectorXd omegasInit_; // initial value for omegas, in case of reset
    VectorXd omegasps_; // partial sum of omegas
    VectorXd epsilons_; // residuals used for ordering in EM
    VectorXi epsOrd_; // order of epsilons
    VectorXd psots_; // partial sum of omegatildas
    VectorXd psoss_; // partial sum of omegas (another s just to distinguish with the name of another double variable)
    ArrayXd logels_; // logel values
  
    // placeholders for lambdaNR
    MatrixXd G_; // G matrix for the estimating equations
    MatrixXd GUsed_;
    VectorXd lambda0_;
    VectorXd lambdaOld_;
    VectorXd lambdaNew_;
    MatrixXd GGt_;
    MatrixXd GGtUsed_;
    VectorXd Glambda_;
    ArrayXd Gl11_;
    VectorXd lGq_;
    VectorXd Q1_;
    MatrixXd Q2_;
    LDLT<MatrixXd> Q2ldlt_;
    VectorXd rho_;
    VectorXd relErr_; // relative difference between two vectors (element-wise absolute difference divided by sum)
    double absErr_; // absolute difference between two scalars

    // for support modification
    int nObs1_; // nObs1_ = nObs_+1: for initial space allocation
    int nObs2_; // nObs2_ = nObs_+support: used nObs_ in calculation
    
    // tolerance values for lambdaNR and EM
    double relTol_;
    double absTol_;
    int maxIter_;
    bool support_;
    double supa_; // tuning parameter for support currection
    
    // maximum relative error in lambda: same for cens / non-cens
    double maxRelErr();
    
    // logsharp and its derivatives for the EL dual problem
    double logsharp(double x, double q); // A support-refined log function
    double logsharp1(double x, double q); // First derivative of logsharp
    double logsharp2(double x, double q); // Second derivative of logsharp
    
    // Helper function for evalWeights: calculate partial sum of omegas.
    double evalPsos(const int ii); // Partial sum of omegas_jj s.t. eps_jj >= eps_ii
    double evalPsosSmooth(const int ii, const double s); // smoothed version
    
  public:
    
    // constructors
    InnerELC();
    InnerELC(int nObs);
    InnerELC(int nObs, int nEqs);
    
    // set options functions
    void setOpts(const int& maxIter, const double& relTol, 
                 const double& absTol, const bool& support, const double& supa, 
                 const Ref<const VectorXd>& lambda0);
    void setOpts(const int& maxIter, const double& relTol, 
                 const double& absTol, const bool& support, 
                 const Ref<const VectorXd>& lambda0);
    void setOpts(const int& maxIter, const double& relTol, 
                 const double& absTol, const bool& support, const double& supa);
    void setOpts(const int& maxIter, const double& relTol, 
                 const double& absTol, const bool& support);
    void setOpts(const int& maxIter, const double& relTol, 
                 const bool& support, const double& supa);
    void setOpts(const int& maxIter, const double& relTol, 
                 const bool& support);
    void setOpts(const bool& support);
    // void setTol(const int& maxIter, const double& relTol, const double& absTol);
    // void setTol(const int& maxIter, const double& relTol);
    
    // core calculations
    void lambdaNR(int& nIter, double& maxErr); // dual solution by Newton-Raphson algorithm
    void evalWeights(); // calculate weights according to epsilons
    void evalOmegas(); // empirical distribution
    double logEL(); // log empirical likelihood
    void evalWeightsSmooth(const double s);
    void evalOmegasSmooth(const double s);
    double logELSmooth(const double s);
    
    // set and get functions 
    void setDeltas(const Ref<const VectorXd>& deltas);
    void setLambda(const Ref<const VectorXd>& lambda); // assigned to lambdaNew
    void setWeights(const Ref<const VectorXd>& weights); 
    void setOmegas(const Ref<const VectorXd>& omegas);
    void setEpsilons(const Ref<const VectorXd>& epsilons);
    void setG(const Ref<const MatrixXd>& G);
    VectorXd getLambda(); 
    VectorXd getWeights(); 
    VectorXd getOmegas(); 
    VectorXd getEpsilons();
    MatrixXd getG();

    // // posterior sampler
    // MatrixXd postSample(int nsamples,int nburn, VectorXd betaInit, 
    //                     const Ref<const VectorXd>& sigs, 
    //                     VectorXd &RvDoMcmc,VectorXd &paccept);
    // void mwgStep(VectorXd &thetaCur, const int &idx, const double &mwgsd,
    //              bool &accept, double &logELCur);
    // MatrixXd postSampleAdapt(int nsamples, int nburn, VectorXd thetaInit,
    //                          double *mwgSd, VectorXd &rvDoMcmc, bool *doAdapt, 
    //                          VectorXd &paccept);
  };
}

/* --------------------------------------------------------------------------- */

// private functions 

// maximum relative error in lambda
template<typename ELModel>
inline double el::InnerELC<ELModel>::maxRelErr() {
  // TODO: added for numerical stability, what is a good tolerance to use ?
  if ((lambdaNew_ - lambdaOld_).array().abs().maxCoeff() < 1e-10) return(0);
  relErr_ = ((lambdaNew_ - lambdaOld_).array() / (lambdaNew_ + lambdaOld_).array()).abs();
  return(relErr_.maxCoeff());
}

// logsharp
template<typename ELModel>
inline double el::InnerELC<ELModel>::logsharp(double x, double q) {
  double xq;
  if(x >= q) {
    return(log(x));
  } else {
    xq = x/q;
    return(-0.5*xq*xq + 2.0*xq - 1.5 + log(q));
  }
}

// d logsharp(x)/dx
template<typename ELModel>
inline double el::InnerELC<ELModel>::logsharp1(double x, double q) {
  if(x >= q) {
    return(1.0/x);
  } else {
    return(-1.0/(q*q)*x + 2.0/q);
  }
}

// d^2 logsharp(x)/dx^2
template<typename ELModel>
inline double el::InnerELC<ELModel>::logsharp2(double x, double q) {
  if(x >= q) {
    return(-1.0/(x*x));
  } else {
    return(-1.0/(q*q));
  }
}

template<typename ELModel>
inline double el::InnerELC<ELModel>::evalPsos(const int ii) {
  double psos = 0;
  int kk;
  for (int jj=nObs2_-1; jj>=0; jj--) {
    // Note: epsOrd corresponds to epsilons acesndingly 
    kk = epsOrd_(jj); // kk is the index of the jj-th largest epsilon
    psos += omegas_(kk);
    if (kk == ii) break; // until (backwardsly) reaching ii-th largest epsilon 
  }
  return psos;
}

template<typename ELModel>
inline double el::InnerELC<ELModel>::evalPsosSmooth(const int ii, const double s) {
  double psos_smooth = 0;
  if (support_ && ii == (nObs2_-1)) {
    psos_smooth = omegas_.head(nObs2_-1).sum() + 0.5*omegas_(nObs2_-1);
  }
  else {
    for (int jj=0; jj<nObs2_; jj++) {
      psos_smooth += ind_smooth(epsilons_(ii)-epsilons_(jj),s)*omegas_(jj);
    }
  }
  return(psos_smooth);
  // VectorXd ss = ArrayXd::Zero(omegas.size()) + s;
  // return (ind_smooth(epsilons[ii]-epsilons.array(),ss).array()*omegas.array()).sum();
}

/* --------------------------------------------------------------------------- */

// public functions

// constructors
/**
 * @brief Default constructor for InnerELC.
 */
template<typename ELModel>
inline el::InnerELC<ELModel>::InnerELC() {}

// ctor with dimensions as input
/**
 * @brief Constructor for InnerELC with number of observations only as input for memory allocation.
 * 
 * @param nObs    Number of observations.
 */
template<typename ELModel>
inline el::InnerELC<ELModel>::InnerELC(int nObs) {
  nObs_ = nObs;
  // support correction
  support_ = false;
  nObs1_ = nObs_+1;
  nObs2_ = nObs_+support_;
  omegas_ = VectorXd::Zero(nObs1_).array() + 1.0/(double)nObs1_; // Initialize omegas to 1/nObs?
}

// ctor with dimensions as input
/**
 * @brief Constructor for InnerELC with dimensions of G matrix as inputs for memory allocation.
 * 
 * @param nObs    Number of observations.
 * @param nEqs    Number of estimating equations.
 */
template<typename ELModel>
inline el::InnerELC<ELModel>::InnerELC(int nObs, int nEqs): ELModel(nObs, nEqs){
  // support correction
  support_ = false;
  nObs1_ = nObs_+1;
  nObs2_ = nObs_+support_;
  // space allocation
  omegasInit_ = VectorXd::Zero(nObs1_).array() + 1.0/(double)nObs1_; 
  omegas_ = VectorXd::Zero(nObs1_).array() + 1.0/(double)nObs1_; 
  G_ = MatrixXd::Zero(nEqs_,nObs1_); // NEW: JAN 1
  GGt_ = MatrixXd::Zero(nEqs_,nObs1_*nEqs_);
  lambda0_ = VectorXd::Zero(nEqs_); 
  lambdaOld_ = VectorXd::Zero(nEqs_); // Initialize to all 0's
  lambdaNew_ = VectorXd::Zero(nEqs_);
  Q1_ = VectorXd::Zero(nEqs_);
  Q2_ = MatrixXd::Zero(nEqs_,nEqs_);
  Glambda_ = VectorXd::Zero(nObs1_);
  Gl11_ = ArrayXd::Zero(nObs1_);
  rho_ = VectorXd::Zero(nObs1_);
  relErr_ = VectorXd::Zero(nEqs_);
  Q2ldlt_.compute(MatrixXd::Identity(nEqs_,nEqs_));
  deltas_ = VectorXd::Zero(nObs1_);
  weights_ = VectorXd::Zero(nObs1_);
  epsilons_ = VectorXd::Zero(nObs1_);
  epsOrd_ = VectorXi::Zero(nObs1_);
  psoss_ = VectorXd::Zero(nObs1_);
  logels_ = ArrayXd::Zero(nObs1_);
  lGq_ = VectorXd::Zero(nObs1_);
}

// set options
/**
 * @brief Set tolerance values for NR, support correction option and tuning parameter, and initial value of lambda.
 * 
 * @param maxIter    Maximum number of iterations.
 * @param relTol     Relative tolerance (for Newton-Raphson algorithm).
 * @param absTol     Absolute tolerance (for EM algorithm).
 * @param support    Whether to have support correction.
 * @param supa       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 * @param lambda0    Initial value for lambda.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setOpts(const int& maxIter, 
                                           const double& relTol, const double& absTol, 
                                           const bool& support, const double& supa,
                                           const Ref<const VectorXd>& lambda0) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  absTol_ = absTol;
  support_ = support;
  supa_ = supa;
  nObs2_ = nObs_+support_;
  lambda0_ = lambda0;
  if (support_) {
    epsilons_.tail(1)(0) = -INFINITY;
    deltas_.tail(1)(0) = 0;
  }
}

/**
 * @brief Set tolerance values for NR, support correction option and tuning parameter, and initial value of lambda.
 * 
 * @param maxIter    Maximum number of iterations.
 * @param relTol     Relative tolerance (for Newton-Raphson algorithm).
 * @param absTol     Absolute tolerance (for EM algorithm).
 * @param support    Whether to have support correction.
 * @param lambda0    Initial value for lambda.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setOpts(const int& maxIter, 
                                           const double& relTol, const double& absTol, 
                                           const bool& support, 
                                           const Ref<const VectorXd>& lambda0) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  absTol_ = absTol;
  support_ = support;
  supa_ = std::max(1.0,0.5*log(nObs_));
  nObs2_ = nObs_+support_;
  lambda0_ = lambda0;
  if (support_) {
    epsilons_.tail(1)(0) = -INFINITY;
    deltas_.tail(1)(0) = 0;
  }
}

/**
 * @brief Set tolerance values for NR, support correction option and tuning parameter, and initial value of lambda.
 * 
 * @param maxIter    Maximum number of iterations.
 * @param relTol     Relative tolerance (for Newton-Raphson algorithm).
 * @param absTol     Absolute tolerance (for EM algorithm).
 * @param support    Whether to have support correction.
 * @param supa       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setOpts(const int& maxIter, 
                                           const double& relTol, const double& absTol, 
                                           const bool& support, const double& supa) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  absTol_ = absTol;
  support_ = support;
  supa_ = supa;
  nObs2_ = nObs_+support_;
  if (support_) {
    epsilons_.tail(1)(0) = -INFINITY;
    deltas_.tail(1)(0) = 0;
  }
}

/**
 * @brief Set tolerance values for NR, support correction option and tuning parameter, and initial value of lambda.
 * 
 * @param maxIter    Maximum number of iterations.
 * @param relTol     Relative tolerance (for Newton-Raphson algorithm).
 * @param absTol     Absolute tolerance (for EM algorithm).
 * @param support    Whether to have support correction.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setOpts(const int& maxIter, 
                                           const double& relTol, const double& absTol, 
                                           const bool& support) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  absTol_ = absTol;
  support_ = support;
  supa_ = std::min((double)nObs_,0.5*log(nObs_));
  nObs2_ = nObs_+support_;
  if (support_) {
    epsilons_.tail(1)(0) = -INFINITY;
    deltas_.tail(1)(0) = 0;
  }
}

/**
 * @brief Set tolerance values for NR, support correction option and tuning parameter, and initial value of lambda.
 * 
 * @param maxIter    Maximum number of iterations.
 * @param relTol     Relative tolerance (for Newton-Raphson algorithm).
 * @param support    Whether to have support correction.
 * @param supa       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setOpts(const int& maxIter, const double& relTol,
                                           const bool& support, const double& supa) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  support_ = support;
  supa_ = supa;
  nObs2_ = nObs_+support_;
  if (support_) {
    epsilons_.tail(1)(0) = -INFINITY;
    deltas_.tail(1)(0) = 0;
  }
}

/**
 * @brief Set tolerance values for NR, support correction option and tuning parameter, and initial value of lambda.
 * 
 * @param maxIter    Maximum number of iterations.
 * @param relTol     Relative tolerance (for Newton-Raphson algorithm).
 * @param support    Whether to have support correction.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setOpts(const int& maxIter, const double& relTol,
                                           const bool& support) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  support_ = support;
  supa_ = std::min((double)nObs_,0.5*log(nObs_));
  nObs2_ = nObs_+support_;
  if (support_) {
    epsilons_.tail(1)(0) = -INFINITY;
    deltas_.tail(1)(0) = 0;
  }
}

/**
 * @brief Set tolerance values for NR, support correction option and tuning parameter, and initial value of lambda.
 * 
 * @param support    Whether to have support correction.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setOpts(const bool& support) {
  support_ = support;
  supa_ = std::min((double)nObs_,0.5*log(nObs_));
  nObs2_ = nObs_+support_;
  if (support_) {
    epsilons_.tail(1)(0) = -INFINITY;
    deltas_.tail(1)(0) = 0;
  }
}

/**
 * @brief Find the optimal lambda by a Newton-Raphson algorithm (with right-censored EL).
 * 
 * @param[out] nIter    Number of iterations to achieve convergence.
 * @param[out] maxErr   Maximum relative error among entires in lambda at the last step.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::lambdaNR(int& nIter, double& maxErr) {

  lambdaOld_ = lambda0_;
  lambdaNew_.fill(0.0);
  
  // Note: these two cannot be preallocate untill `support` is set
  // but has to do it this way since .block does not return a MatrixXd (TODO: alternative way?)
  GGtUsed_ = GGt_.block(0,0,nEqs_,nObs2_*nEqs_);
  GUsed_ = G_.block(0,0,nEqs_,nObs2_);
  
  block_outer(GGtUsed_,GUsed_);
  
  // To avoid numerical problem (TODO: test this again and see if it is ok without it?)
  VectorXd weights_save = weights_;
  for (int ll=0;ll<nObs2_; ll++) {
    if (weights_save(ll) < relTol_*nObs_) weights_save(ll) = 0;
  }
  
  // newton-raphson loop
  int ii;
  for(ii=0; ii<maxIter_; ii++) {
    // Q1 and Q2
    Glambda_.noalias() = lambdaOld_.transpose() * GUsed_;
    Glambda_ = weights_.head(nObs2_).sum() + Glambda_.array();
    Q2_.fill(0.0); 
    for(int jj=0; jj<nObs2_; jj++) {
      rho_(jj) = logsharp1(Glambda_(jj), weights_(jj));
      
      // To avoid numerical problem
      if (weights_(jj) > relTol_*nObs_) {
        Q2_ += weights_(jj) * logsharp2(Glambda_(jj), weights_(jj)) * GGtUsed_.block(0,jj*nEqs_,nEqs_,nEqs_);
      }
    }
    // To avoid numerical problem
    for (int kk=0; kk<nObs2_; kk++) {
      if (isinf(rho_(kk)) && weights_save(kk) == 0) rho_(kk) = 0;
    }
    // update lambda
    Q1_ = GUsed_ * (rho_.array()*weights_save.head(nObs2_).array()).matrix();
    Q2ldlt_.compute(Q2_);
    lambdaNew_.noalias() = lambdaOld_ - Q2ldlt_.solve(Q1_);
    maxErr = maxRelErr();
    if (maxErr < relTol_) {
      break;
    }
    lambdaOld_ = lambdaNew_; // complete cycle
  }
  nIter = ii; // output lambda and also nIter and maxErr
  return;
}


/**
 * @brief Calculate weights for weighted log EL in EM according to epsilons.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::evalWeights() {
  // find the indices for increasing order of epsilons 
  psots_.fill(0.0);
  int kk;
  double psos;
  for (int ii=0; ii<nObs2_; ii++) {
    for (int jj=0; jj<nObs2_; jj++) {
      kk = epsOrd_(jj);
      if (deltas_(kk) == 0) {
        psos = evalPsos(kk);
        // to prevent dividing by 0
        if (abs(psos) >= 1e-10) psots_(ii) += omegas_(ii)/psos;
        else if (omegas_(ii) >= 1e-10 && evalPsos(kk) < 1e-10) {
          // TODO: this means a problem
          std::cout << "evalWeights: dividing by 0 problem." << std::endl;
        }
      }
      if (kk == ii) break;
    }
  }
  weights_.array() = deltas_.array() + psots_.array();
}


/**
* @brief Evaluate omegas using an EM algorithm.
*/
template<typename ELModel>
inline void el::InnerELC<ELModel>::evalOmegas() {
  // Need to have a valid starting value if the last one is not valid in MCMC
  if (omegas_ != omegas_) {
    // std::cout << "evalOmegas: resetting omegas_." << std::endl;
    omegas_ = omegasInit_;
  }
  int nIter;
  double maxErr;
  // lGq_(nObs2_);
  double logelOld = logEL();
  double logel = logelOld;
  int ii;
  for(ii=0; ii<maxIter_; ii++) {
    // E-step:
    evalWeights(); // assigns weights according to epsilons
    // M-step:
    lambdaNR(nIter, maxErr);
    lGq_.head(nObs2_) = ((lambdaNew_.transpose() * G_.block(0,0,nEqs_,nObs2_)).array() + weights_.head(nObs2_).sum()).transpose();
    omegas_.head(nObs2_).array() = weights_.head(nObs2_).array() / lGq_.head(nObs2_).array();
    omegas_.head(nObs2_).array() = omegas_.head(nObs2_).array() / (omegas_.head(nObs2_).array().sum()); // normalize
    logel = logEL();
    maxErr = abs(logel-logelOld); // absolute error in log EL
    if (maxErr < absTol_) break;
    logelOld = logel;
  }
  nIter = ii;
  if (nIter == maxIter_ && maxErr > absTol_) {
    // TODO: maybe should assign nan elsewhere
    for (int ii=0; ii<nObs2_; ii++) {
      omegas_(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  return;
}

/**
* @brief Calculate logEL using omegas and deltas.
*/
template<typename ELModel>
inline double el::InnerELC<ELModel>::logEL() {
  // evalOmegas(maxIter,relTol); 
  if (omegas_ != omegas_) return -INFINITY; // (NaN is not equal to themselves)
  else if ((omegas_.head(nObs_).array() < -1e-10/nObs2_).any()) return -INFINITY;
  else {
    omegas_ = omegas_.array().abs();
    for (int ii=0; ii<nObs2_; ii++) {
      psoss_(ii) = evalPsos(ii);
    }
    for (int jj=0; jj<nObs2_; jj++) {
      if (deltas_(jj) == 1) logels_(jj) = log(omegas_(jj));
      else  logels_(jj) = log(psoss_(jj));
    }
    return(logels_.head(nObs2_).sum());
  }
}

/**
 * @brief Calculate weights for weighted log EL in EM according to epsilons.
 * 
 * @param s        Tuning parameter for smoothing.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::evalWeightsSmooth(const double s) {
  psots_.fill(0.0); // TODO initialize psots somewhere else???
  for (int ii=0; ii<nObs2_; ii++) {
    psoss_(ii) = evalPsosSmooth(ii,s);
  }
  if (support_) {
    for (int jj=0; jj<nObs2_; jj++) {
      for (int kk=0; kk<nObs2_; kk++) {
        if (jj == nObs2_-1 && kk == nObs2_-1) {
          psots_(jj) += (1-deltas_(kk))*ind_smooth(0.0,s)*omegas_(jj)/psoss_(kk);
        }
        else psots_(jj) += (1-deltas_(kk))*ind_smooth(epsilons_(kk)-epsilons_(jj),s)*omegas_(jj)/psoss_(kk);
      }
    }
  }
  else {
    for (int jj=0; jj<nObs2_; jj++) {
      for (int kk=0; kk<nObs2_; kk++) {
        psots_(jj) += (1-deltas_(kk))*ind_smooth(epsilons_(kk)-epsilons_(jj),s)*omegas_(jj)/psoss_(kk);
      }
    }
  }
  weights_ = deltas_.array()+psots_.array();
}

/**
 * @brief Evaluate omegas using an EM algorithm with continuity correction.
 * 
 * @param s        Tuning parameter for smoothing.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::evalOmegasSmooth(const double s) {
  if (omegas_ != omegas_) {
    // std::cout << "evalOmegas: resetting omegas." << std::endl;
    omegas_ = omegasInit_;
  }
  int nIter;
  double maxErr;
  double logelOld = logELSmooth(s);
  double logel = logelOld;
  int ii;
  for(ii=0; ii<maxIter_; ii++) {
    // E-step:
    evalWeightsSmooth(s); // assigns weights according to epsilons
    // M-step:
    lambdaNR(nIter, maxErr);
    lGq_.head(nObs2_) = ((lambdaNew_.transpose() * G_.block(0,0,nEqs_,nObs2_)).array() + weights_.head(nObs2_).sum()).transpose();
    omegas_.head(nObs2_).array() = weights_.head(nObs2_).array() / lGq_.head(nObs2_).array();
    omegas_.head(nObs2_).array() = omegas_.head(nObs2_).array() / (omegas_.head(nObs2_).array().sum()); // normalize
    logel = logELSmooth(s);
    absErr_ = abs(logel-logelOld);
    if (absErr_ < absTol_) break;
    logelOld = logel;
  }
  nIter = ii; 
  if (nIter == maxIter_ && absErr_ > absTol_) {
    // TODO: maybe should assign nan elsewhere 
    for (int ii=0; ii<nObs2_; ii++) {
      omegas_(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  return;
}

/**
 * @brief Calculate logEL using omegas and deltas.
 * 
 * @param s        Tuning parameter for smoothing.
 */
template<typename ELModel>
inline double el::InnerELC<ELModel>::logELSmooth(const double s) {
  if (omegas_ != omegas_) return -INFINITY; // (NaN is not equal to themselves)
  else if ((omegas_.array() < -1e-10/omegas_.size()).any()) return -INFINITY;
  else {
    omegas_ = omegas_.array().abs();
    // psoss_ = VectorXd::Zero(nObs2_);
    for (int ii=0; ii<nObs2_; ii++) {
      psoss_(ii) = evalPsosSmooth(ii,s);
    }
    // logels_ = ArrayXd::Zero(nObs2_);
    for (int jj=0; jj<nObs2_; jj++) {
      if (deltas_(jj) == 1) logels_(jj) = log(omegas_(jj));
      else  logels_(jj) = log(psoss_(jj));
    }
    return(logels_.head(nObs2_).sum());
  }
}

// setters

/**
 * @brief Set the value of epsilons.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setEpsilons(const Ref<const VectorXd>& epsilons) {
  epsilons_.head(nObs_) = epsilons;
  epsOrd_.head(nObs2_) = sort_inds(epsilons_.head(nObs2_)); 
}

/**
 * @brief Set the value of lambda.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setLambda(const Ref<const VectorXd>& lambda) {
  lambdaNew_ = lambda;
}

/**
 * @brief Set the value of weights.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setWeights(const Ref<const VectorXd>& weights) {
  weights_.head(nObs_) = weights; 
}

/**
 * @brief Set the value of omegas.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setOmegas(const Ref<const VectorXd>& omegas) {
  // nObs = _omegas.size(); // TODO: where to set nObs
  omegasInit_.head(nObs2_) = omegas; // new
  omegas_.head(nObs2_) = omegas; 
}

/**
 * @brief Set the value of deltas.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setDeltas(const Ref<const VectorXd>& deltas) {
  psots_ = VectorXd::Zero(nObs2_); // TODO: initialization maybe should be done at better places
  deltas_.head(nObs_) = deltas;
}

/**
 * @brief Set the value of G.
 */
template<typename ELModel>
inline void el::InnerELC<ELModel>::setG(const Ref<const MatrixXd>& G) {
  G_.block(0,0,nEqs_,nObs_) = G;
  if (support_) adj_G(G_,supa_);
}

// getters

/**
 * @brief Get the value of lambda.
 */
template<typename ELModel>
inline VectorXd el::InnerELC<ELModel>::getLambda() {
  return(lambdaNew_);
}

/**
 * @brief Get the value of weights.
 */
template<typename ELModel>
inline VectorXd el::InnerELC<ELModel>::getWeights() {
  return(weights_.head(nObs2_));
}

/**
 * @brief Get the value of omegas.
 */
template<typename ELModel>
inline VectorXd el::InnerELC<ELModel>::getOmegas() {
  return(omegas_.head(nObs2_));
}

/**
 * @brief Get the value of epsilons.
 */
template<typename ELModel>
inline VectorXd el::InnerELC<ELModel>::getEpsilons() {
  return(epsilons_.head(nObs2_));
}

/**
 * @brief Get the value of G.
 */
template<typename ELModel>
inline MatrixXd el::InnerELC<ELModel>::getG() {
  return(G_.block(0,0,nEqs_,nObs2_));
}

/* TAKE OUT ALL MCMC SAMPLERS FOR NOW:
// posterior sampler: location model, single quantile case
// Note: omegasInit needed as the starting value for EM?
template<typename ELModel>
inline MatrixXd InnerELC<ELModel>::postSample(int nsamples, int nburn,
                                              VectorXd betaInit, 
                                              const Ref<const VectorXd>& sigs,
                                              VectorXd &RvDoMcmc, VectorXd &paccept) {
  std::cout << "--------------------------- ELC postSample ----------------------------" << std::endl;
  VectorXd betaOld = betaInit;
  VectorXd betaNew = betaInit;
  VectorXd betaProp = betaInit;
  int betalen = betaInit.size();
  MatrixXd beta_chain(betalen,nsamples); 
  // beta_chain.fill(0.0); // debug
  ELModel::evalG(betaInit);
  evalEpsilons(betaInit);
  evalOmegas(); // omegasInit should have been assigned
  double logELOld = logEL();
  // std::cout << "after calling logEL()." << std::endl;
  double logELProp;
  bool satisfy;
  double u;
  double a;
  double ratio;
  int nIter;
  double maxErr;
  // paccept.fill(0.0);
  paccept = VectorXd::Zero(betalen);
  // std::cout << "sigs = " << sigs.transpose() << std::endl;
  // std::cout << "paccept = " << paccept.transpose() << std::endl;
  // std::cout << "RvDoMcmc = " << RvDoMcmc.transpose() << std::endl;

  for (int ii=-nburn; ii<nsamples; ii++) {
    if (ii % 200 == 0) {
      std::cout << "ii = " << ii << std::endl;
    }
    for (int jj=0; jj<betalen; jj++) {
      if (RvDoMcmc(jj)) {
        betaProp = betaOld;
        betaProp(jj) += sigs(jj)*R::norm_rand();
        // std::cout << "betaProp = " << betaProp.transpose() << std::endl;
        // check if proposed beta satisfies the constraint
        satisfy = false;
        ELModel::evalG(betaProp);
        evalEpsilons(betaProp);
        evalWeights();
        // std::cout << "weights = " << weights.transpose() << std::endl;
        // check whether the constrains are satisfied or not
        // std::cout << "old lambdaNew = " << lambdaNew.transpose() << std::endl;
        lambdaNR(nIter, maxErr);
        // std::cout << "new lambdaNew = " << lambdaNew << std::endl;
        if (nIter < maxIter) satisfy = true;
        // if does not satisfy, keep the old beta
        if (satisfy == false) break;
        // if does satisfy, decide whether to accept
        u = R::unif_rand();
        // use the lambda calculate just now to get the logEL for Prop
        // to avoid an extra call of lambdaNR
        // VectorXd logomegahat = log(weights.array() /
        //   (weights.transpose() + lambdaNew.transpose() * G).array().sum());
        // logELProp = (weights.array()*logomegahat.array()).sum();
        evalOmegas(); // lambdaNR is in here as well
        // std::cout << "omegas = " << omegas.transpose() << std::endl;
        logELProp = logEL();
        // std::cout << "logELOld = " << logELOld << std::endl;
        // std::cout << "logELProp = " << logELProp << std::endl;
        ratio = exp(logELProp-logELOld);
        // std::cout << "ratio = " << ratio << std::endl;
        a = std::min(1.0,ratio);
        if (u < a) { // accepted
          paccept(jj) += 1;
          betaNew = betaProp;
          betaOld = betaNew;
          // betaOld = betaProp;
          logELOld = logELProp; // store the new one
        }
      }
    }
    if (ii >= 0) {
      beta_chain.col(ii) = betaNew;
      // beta_chain.col(ii) = betaProp;
    }
  }
  paccept /= (nsamples+nburn);
  return(beta_chain);
}

template<typename ELModel>
inline void InnerELC<ELModel>::mwgStep(VectorXd &thetaCur,
                                       const int &idx,
                                       const double &mwgsd,
                                       bool &accept, 
                                       double &logELCur) {
  // std::cout << "in mwgStep: logELCur = " << logELCur << std::endl;
  int nTheta = thetaCur.size();
  accept = false;
  VectorXd thetaProp = thetaCur;
  thetaProp(idx) += mwgsd*R::norm_rand();
  // sig2 has to be positive
  if (idx == nBet+nGam && thetaProp(idx) < 0) {
    // std::cout << "negative sig2 = " << thetaProp(idx) << std::endl;
    return;
  }
  
  if (nTheta == nBet) {
    ELModel::evalG(thetaProp);
    evalEpsilons(thetaProp);
  }
  else {
    if (nTheta == nBet + nGam + 1){
      ELModel::evalG(thetaProp.head(nBet), 
                     thetaProp.segment(nBet,nGam), 
                     thetaProp.tail(1)(0),
                     VectorXd::Zero(0));
    }
    else {
      ELModel::evalG(thetaProp.head(nBet), 
                     thetaProp.segment(nBet,nGam), 
                     thetaProp.segment(nBet+nGam,1)(0),
                     thetaProp.tail(nQts));
    }
    evalEpsilons(thetaProp.head(nBet),
                 thetaProp.segment(nBet,nGam),
                 thetaProp.segment(nBet+nGam,1)(0));
  }
  int nIter;
  double maxErr;
  evalWeights();
  lambdaNR(nIter, maxErr);
  bool satisfy = false;
  if (nIter < maxIter || maxErr <= relTol) satisfy = true;
  // if does not satisfy, keep the old theta
  if (satisfy == false) {
    // std::cout << "lambdaNR not converged so keep the old one."<< std::endl;
    // std::cout << "nIter = " << nIter << std::endl;
    // std::cout << "maxErr = " << maxErr << std::endl;
    return;
  }
  // if does satisfy, flip a coin
  double u = R::unif_rand();
  // std::cout << "before evalOmegas" << std::endl;
  evalOmegas();
  // std::cout << "omegas = " << omegas.transpose() << std::endl;
  double logELProp = logEL();
  // std::cout << "logELProp = " << logELProp << std::endl;
  double ratio = exp(logELProp-logELCur);
  double a = std::min(1.0,ratio);
  if (u < a) { // accepted
    accept = true;
    thetaCur = thetaProp;
    logELCur = logELProp;
  }
}

// template<typename ELModel>
// inline MatrixXd InnerELC<ELModel>::postSampleAdapt(int nsamples, int nburn, 
//                                                    VectorXd thetaInit,
//                                                    double *mwgSd, bool *rvDoMcmc, 
//                                                    VectorXd &paccept) {
template<typename ELModel>
inline MatrixXd InnerELC<ELModel>::postSampleAdapt(int nsamples, int nburn, 
                                                   VectorXd thetaInit,
                                                   double *mwgSd,
                                                   VectorXd &rvDoMcmc,
                                                   bool *doAdapt, 
                                                   VectorXd &paccept) {
  int nTheta = thetaInit.size();
  MatrixXd theta_chain(nTheta,nsamples);
  paccept = VectorXd::Zero(nTheta);
  // theta_chain.fill(0.0); // debug
  // paccept.fill(0.0);
  MwgAdapt tuneMCMC(nTheta, doAdapt);
  bool *isAccepted = new bool[nTheta];
  for (int ii=0; ii<nTheta; ii++) {
    isAccepted[ii] = false;
  }
  VectorXd thetaCur = thetaInit;
  if (nTheta == nBet) {
    std::cout << "location model." << std::endl;
    ELModel::evalG(thetaCur);
    evalEpsilons(thetaCur);
  }
  else {
    if (nTheta == nBet + nGam + 1) {
      std::cout << "location-scale model for mean regression." << std::endl;
      ELModel::evalG(thetaCur.head(nBet),
                     thetaCur.segment(nBet,nGam),
                     thetaCur.tail(1)(0),
                     VectorXd::Zero(0));
    }
    else {
      std::cout << "location-scale model for quantile regression." << std::endl;
      ELModel::evalG(thetaCur.head(nBet),
                     thetaCur.segment(nBet,nGam),
                     thetaCur.segment(nBet+nGam,1)(0),
                     thetaCur.tail(nQts));
    }
    evalEpsilons(thetaCur.head(nBet),
                 thetaCur.segment(nBet,nGam),
                 thetaCur.segment(nBet+nGam,1)(0)); // Note: have to take the element here
  }
  evalOmegas();
  // std::cout << "initial omegas = " << omegas.transpose() << std::endl;
  double logELCur = logEL();
  // MCMC loop
  for(int ii=-nburn; ii<nsamples; ii++) {
    // if (ii % 200 == 0) {
    //   std::cout << "ii = " << ii << std::endl;
    //   
    //   // DEBUG begins
    //   std::cout << "mwgSd = ";
    //   for (int ll=0; ll<nTheta; ll++) {
    //     std::cout << mwgSd[ll] << ' ';
    //   }
    //   std::cout << std::endl;
    // 
    //   std::cout << "isAccepted = ";
    //   for (int ll=0; ll<nTheta; ll++) {
    //     std::cout << isAccepted[ll] << ' ';
    //   }
    //   std::cout << std::endl;
    //   // DEBUG ends
    // }
    for(int jj=0; jj<nTheta; jj++) {
      if(rvDoMcmc[jj]) {
        // modifies thetaCur's jj-th entry
        mwgStep(thetaCur,jj,mwgSd[jj],isAccepted[jj],logELCur);
        // std::cout << "in sampler loop: thetaCur = " << thetaCur.transpose() << std::endl;
        // std::cout << "mwgSd[jj] = " << mwgSd[jj] << std::endl;
        if (isAccepted[jj]) paccept(jj) += 1; // add 1 to paccept if accepted
      }
    }
    if (ii >= 0) {
      theta_chain.col(ii) = thetaCur;
    }
    tuneMCMC.adapt(mwgSd, isAccepted);
    // TODO: temp to check if the chain is converging..
    if (ii > 0 && ii % 1000 == 0) {
      int noconv = 0;
      VectorXd pacctemp = paccept.array() / (nburn+ii);
      std::cout << "pactemp = "; 
      for (int kk=0; kk<nTheta; kk++) {
        std::cout << pacctemp(kk) << ' ';
        if (pacctemp(kk) < 0.4) noconv += 1;
      }
      std::cout << std::endl;
      if (noconv == nTheta) {
        theta_chain.fill(0.0);
        break;
      }
    }
  }
  paccept /= (nsamples+nburn);
  delete[] isAccepted; // deallocate memory
  return(theta_chain);
}
*/


#endif
