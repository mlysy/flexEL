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
// #include "MwgAdapt.h" // for adaptive mcmc

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

/**
 * @brief el namespace
 * 
 * Wrap the exported library components into a namespace called \b el to avoid potential naming conflicts with other libraries or user-defined headers.
 */
namespace el {
  
  template <typename ELModel>
  class InnerELC : public ELModel { 
  private:
    
    // required members in ELModel
    // using ELModel::G_; // REMOVED: JAN 1
    using ELModel::nObs_;
    using ELModel::nEqs_;
    
    // placeholders for EM algorithm
    VectorXd deltas_; /**< censoring indicators */
    VectorXd weights_; /**< weights in weighted maximum log EL */
    VectorXd omegas_; /**< empirical distribution */
    VectorXd omegasInit_; /**< initial value for omegas, in case of reset */
    VectorXd omegasps_; /**< partial sum of omegas */
    VectorXd epsilons_; /**< residuals used for ordering in EM */
    VectorXi epsOrd_; /**< order of epsilons */
    VectorXd psots_; /**< partial sum of omegatildas */
  
    // placeholders for lambdaNR
    MatrixXd G_; // NEW: JAN 1
    VectorXd lambda0_; // TODO: not used yet!
    VectorXd lambdaOld_;
    VectorXd lambdaNew_;
    MatrixXd GGt_;
    VectorXd Glambda_;
    ArrayXd Gl11_;
    VectorXd Q1_;
    MatrixXd Q2_;
    LDLT<MatrixXd> Q2ldlt_;
    VectorXd rho_;
    VectorXd relErr_; /**< relative difference between two vectors (element-wise absolute difference divided by sum) */
    double absErr_; /**< absolute difference between two scalars */

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
    
    /**
     * @brief      Helper function for evalWeights: calculate partial sum of omegas.
     * 
     * @return     Partial sum of omegas_jj s.t. eps_jj >= eps_ii.
     */
    double evalPsos(const int ii);
    
  public:
    
    /**
     * @brief Default constructor for InnerELC.
     */
    InnerELC(); // default ctor
    
    /**
     * @brief Constructor for InnerELC with number of observations only as input for memory allocation.
     * 
     * @param nObs    Number of observations.
     */
    InnerELC(int nObs);
    
    /**
     * @brief Constructor for InnerELC with dimensions of G matrix as inputs for memory allocation.
     * 
     * @param nObs    Number of observations.
     * @param nEqs    Number of estimating equations.
     */
    InnerELC(int nObs, int nEqs);
    
    // logsharp and its derivatives for the EL dual problem
    /**
     * @brief A support-refined log function.
     */
    double logsharp(double x, double q);
    
    /**
     * @brief First derivative of logsharp.
     */
    double logsharp1(double x, double q);
    
    /**
     * @brief Second derivative of logsharp.
     */
    double logsharp2(double x, double q);
    
    // Set options
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
    
    // Newton-Raphson algorithm
    /**
     * @brief Set tolerance values for NR and EM (EM uses absTol to determine the convergence of omegas).
     */
    void setTol(const int& maxIter, const double& relTol, const double& absTol);
    
    /**
     * @brief Set tolerance values for NR.
     */
    void setTol(const int& maxIter, const double& relTol);
    
    /**
     * @brief Find the optimal lambda by a Newton-Raphson algorithm (with right-censored EL).
     * 
     * @param[out] nIter    Number of iterations to achieve convergence.
     * @param[out] maxErr   Maximum relative error among entires in lambda at the last step.
     */
    void lambdaNR(int& nIter, double& maxErr);
    
    // eval functions 
    /**
     * @brief Calculate weights for weighted log EL in EM according to epsilons.
     */
    void evalWeights(); // calculate weights according to epsilons 
    
    /**
     * @brief Evaluate omegas using an EM algorithm.
     */
    void evalOmegas();
  
    /**
     * @brief Calculate logEL using omegas and deltas.
     */
    double logEL();
    
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
    
    // smoothed version functions
    // VectorXd indSmooth(VectorXd x, VectorXd s);
    double evalPsosSmooth(const int ii, const double s); // any way to test while set as private?
    double logELSmooth(const double s);
    void evalWeightsSmooth(const double s);
    void evalOmegasSmooth(const double s);
    
    // void evalEpsilons(const Ref<const VectorXd>& beta); // for location model 
    // void evalEpsilons(const Ref<const VectorXd>& beta,
    //                   const Ref<const VectorXd>& gamma,
    //                   const double& sig2); // for location-scale model
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

// default ctor
template<typename ELModel>
inline el::InnerELC<ELModel>::InnerELC() {}

// ctor with dimensions as input
template<typename ELModel>
inline el::InnerELC<ELModel>::InnerELC(int nObs) {
  nObs_ = nObs;
  omegas_ = VectorXd::Zero(nObs_).array() + 1.0/(double)nObs_; // Initialize omegas to 1/nObs
  // support correction
  support_ = false;
  nObs1_ = nObs_+1;
  nObs2_ = nObs_+support_;
}

// ctor with dimensions as input
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
}

// set options
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
}

template<typename ELModel>
inline void el::InnerELC<ELModel>::setOpts(const int& maxIter, 
                                           const double& relTol, const double& absTol, 
                                           const bool& support, 
                                           const Ref<const VectorXd>& lambda0) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  absTol_ = absTol;
  support_ = support;
  supa_ = std::min((double)nObs_,0.5*log(nObs_));
  nObs2_ = nObs_+support_;
  lambda0_ = lambda0;
}

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
}

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
}

template<typename ELModel>
inline void el::InnerELC<ELModel>::setOpts(const int& maxIter, const double& relTol,
                                           const bool& support, const double& supa) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  support_ = support;
  supa_ = supa;
  nObs2_ = nObs_+support_;
}

template<typename ELModel>
inline void el::InnerELC<ELModel>::setOpts(const int& maxIter, const double& relTol,
                                           const bool& support) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  support_ = support;
  supa_ = std::min((double)nObs_,0.5*log(nObs_));
  nObs2_ = nObs_+support_;
}

template<typename ELModel>
inline void el::InnerELC<ELModel>::setOpts(const bool& support) {
  support_ = support;
  supa_ = std::min((double)nObs_,0.5*log(nObs_));
  nObs2_ = nObs_+support_;
}

// set tolerance for NR
template<typename ELModel>
inline void el::InnerELC<ELModel>::setTol(const int& maxIter, const double& relTol) {
  maxIter_ = maxIter;
  relTol_ = relTol;
}

// set tolerance for NR and EM
template<typename ELModel>
inline void el::InnerELC<ELModel>::setTol(const int& maxIter,
                                          const double& relTol, const double& absTol) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  absTol_ = absTol;
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

// maximum relative error in lambda
template<typename ELModel>
inline double el::InnerELC<ELModel>::maxRelErr() {
  // TODO: added for numerical stability, what is a good tolerance to use ?
  if ((lambdaNew_ - lambdaOld_).array().abs().maxCoeff() < 1e-10) return(0);
  relErr_ = ((lambdaNew_ - lambdaOld_).array() / (lambdaNew_ + lambdaOld_).array()).abs();
  return(relErr_.maxCoeff());
}

// Newton-Raphson algorithm
// template<typename ELModel>
// inline void el::InnerELC<ELModel>::lambdaNR(int& nIter, double& maxErr) {
//   // prevent bad starting values
//   lambdaNew_.fill(0.0);
//   lambdaOld_.fill(0.0);
//   block_outer(GGt_,G_);
//   int ii;
//   
//   // To avoid numerical problem
//   VectorXd weights_save = weights_;
//   for (int ll=0;ll<nObs_; ll++) {
//     if (weights_save(ll) < relTol_*nObs_) weights_save(ll) = 0;
//   }
//   
//   // newton-raphson loop
//   for(ii=0; ii<maxIter_; ii++) {
//     // Q1 and Q2
//     Glambda_.noalias() = lambdaOld_.transpose() * G_;
//     Glambda_ = weights_.sum() + Glambda_.array();
//     Q2_.fill(0.0); 
//     for(int jj=0; jj<nObs_; jj++) {
//       rho_(jj) = logsharp1(Glambda_(jj), weights_(jj));
//       
//       // To avoid numerical problem
//       if (weights_(jj) > relTol_*nObs_) {
//         Q2_ += weights_(jj) * logsharp2(Glambda_(jj), weights_(jj)) * GGt_.block(0,jj*nEqs_,nEqs_,nEqs_);
//       }
//     }
//     
//     // To avoid numerical problem
//     for (int kk=0; kk<nObs_; kk++) {
//       if (isinf(rho_(kk)) && weights_save(kk) == 0) rho_(kk) = 0;
//     }
//     
//     Q1_ = G_ * (rho_.array()*weights_save.array()).matrix();
//     // update lambda
//     Q2ldlt_.compute(Q2_);
//     lambdaNew_.noalias() = lambdaOld_ - Q2ldlt_.solve(Q1_);
//     // maxErr = maxRelErr(lambdaNew_, lambdaOld_); // maximum relative error
//     maxErr = maxRelErr();
//     if (maxErr < relTol_) {
//         break;
//     }
//     lambdaOld_ = lambdaNew_; // complete cycle
//   }
//   nIter = ii; // output lambda and also nIter and maxErr
//   return;
// }
template<typename ELModel>
inline void el::InnerELC<ELModel>::lambdaNR(int& nIter, double& maxErr) {
  // prevent bad starting values
  lambdaOld_ = lambda0_;
  lambdaNew_.fill(0.0);
  
  MatrixXd GGtUsed_ = GGt_.block(0,0,nEqs_,nObs2_*nEqs_);
  MatrixXd GUsed_ = G_.block(0,0,nEqs_,nObs2_);
  block_outer(GGtUsed_,GUsed_);
  
  // To avoid numerical problem
  VectorXd weights_save = weights_;
  for (int ll=0;ll<nObs2_; ll++) {
    if (weights_save(ll) < relTol_*nObs_) weights_save(ll) = 0;
  }
  
  // newton-raphson loop
  int ii;
  for(ii=0; ii<maxIter_; ii++) {
    // Q1 and Q2
    Glambda_.head(nObs2_).noalias() = lambdaOld_.transpose() * GUsed_;
    Glambda_.head(nObs2_) = weights_.head(nObs2_).sum() + Glambda_.head(nObs2_).array();
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
    
    Q1_ = GUsed_ * (rho_.array()*weights_save.head(nObs2_).array()).matrix();
    // update lambda
    Q2ldlt_.compute(Q2_);
    lambdaNew_.noalias() = lambdaOld_ - Q2ldlt_.solve(Q1_);
    // maxErr = maxRelErr(lambdaNew_, lambdaOld_); // maximum relative error
    maxErr = maxRelErr();
    if (maxErr < relTol_) {
      break;
    }
    lambdaOld_ = lambdaNew_; // complete cycle
  }
  nIter = ii; // output lambda and also nIter and maxErr
  return;
}

template<typename ELModel>
inline VectorXd el::InnerELC<ELModel>::getLambda() {
  return(lambdaNew_);
}

// returns partial sum of omegas_ according to the epsilons(jj) that are larger  
//   than epsilons(ii)
// Note: epsOrd must have been assigned, i.e. have called sort_inds(epsilons); 
// template<typename ELModel>
// inline double el::InnerELC<ELModel>::evalPsos(const int ii) {
//   double psos = 0;
//   int kk;
//   for (int jj=nObs_-1; jj>=0; jj--) {
//     // Note: epsOrd corresponds to epsilons acesndingly 
//     kk = epsOrd_(jj); // kk is the index of the jj-th largest epsilon
//     psos += omegas_(kk);
//     if (kk == ii) break; // until (backwardsly) reaching ii-th largest epsilon 
//   }
//   return psos;
// }
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
inline void el::InnerELC<ELModel>::setEpsilons(const Ref<const VectorXd>& epsilons) {
  epsilons_ = epsilons;
  epsOrd_ = sort_inds(epsilons_); 
}

// // evaluate epsilons location model
// template<typename ELModel>
// inline void InnerELC<ELModel>::evalEpsilons(const Ref<const VectorXd>& beta) {
//   epsilons = (y.transpose() - beta.transpose() * X).transpose(); 
//   epsOrd = sort_inds(epsilons);
// }

// // evaluate epsilons location-scale model
// template<typename ELModel>
// inline void InnerELC<ELModel>::evalEpsilons(const Ref<const VectorXd>& beta,
//                                             const Ref<const VectorXd>& gamma,
//                                             const double& sig2) {
//   if (sig2 < 0.0) {
//     std::cout << "negative sig2 in evalEpsilons." << std::endl;
//   }
//   epsilons = (y.transpose() - beta.transpose() * X).transpose(); 
//   epsilons.transpose().array() *= (-gamma.transpose()*Z).array().exp()/sqrt(sig2);
//   epsOrd = sort_inds(epsilons);
// }

// Note: epsilons must have been assigned
// template<typename ELModel>
// inline void el::InnerELC<ELModel>::evalWeights() {
//   // find the indices for increasing order of epsilons 
//   psots_.fill(0.0);
//   int kk;
//   double psos;
//   for (int ii=0; ii<nObs_; ii++) {
//     for (int jj=0; jj<nObs_; jj++) {
//       kk = epsOrd_(jj);
//       if (deltas_(kk) == 0) {
//         psos = evalPsos(kk);
//         // to prevent dividing by 0
//         if (abs(psos) >= 1e-10) psots_(ii) += omegas_(ii)/psos;
//         else if (omegas_(ii) >= 1e-10 && evalPsos(kk) < 1e-10) {
//           // TODO: this means a problem
//           std::cout << "evalWeights: dividing by 0 problem." << std::endl;
//         }
//       }
//       if (kk == ii) break;
//     }
//   }
//   // assigned weights are still in the order of original data
//   weights_.array() = deltas_.array() + psots_.array();
// }
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
  // assigned weights are still in the order of original data
  weights_.array() = deltas_.array() + psots_.array();
}

// exported as `omega.hat.EM`
// Note: convergence of EM depends on change in the value of logEL
// Note: epsilons must have been assigned
// Note: maxIter and relTol are the same as in lambdaNR, and should have been assigned
// template<typename ELModel>
// inline void el::InnerELC<ELModel>::evalOmegas() {
//   // if (omegas != omegas) return; // if initial value nan, stop
//   // Need to have a valid starting value if the last one is not valid in MCMC
//   if (omegas_ != omegas_) {
//     // std::cout << "evalOmegas: resetting omegas_." << std::endl;
//     omegas_ = omegasInit_;
//   }
//   int ii;
//   int nIter;
//   double maxErr;
//   VectorXd lGq;
//   double logelOld = logEL();
//   double logel = logelOld;
//   for(ii=0; ii<maxIter_; ii++) {
//     // E-step:
//     evalWeights(); // assigns weights according to epsilons
//     // M-step:
//     lambdaNR(nIter, maxErr);
//     lGq = ((lambdaNew_.transpose() * G_).array() + weights_.sum()).transpose();
//     omegas_.array() = weights_.array() / lGq.array();
//     omegas_.array() = omegas_.array() / (omegas_.array().sum()); // normalize
//     logel = logEL();
//     maxErr = abs(logel-logelOld); // absolute error in log EL
//     if (maxErr < absTol_) break;
//     logelOld = logel;
//   }
//   nIter = ii; 
//   if (nIter == maxIter_ && maxErr > absTol_) {
//     // TODO: maybe should assign nan elsewhere 
//     for (int ii=0; ii<nObs_; ii++) {
//       omegas_(ii) = std::numeric_limits<double>::quiet_NaN();
//     }
//   }
//   return;
// }
template<typename ELModel>
inline void el::InnerELC<ELModel>::evalOmegas() {
  // if (omegas != omegas) return; // if initial value nan, stop
  // Need to have a valid starting value if the last one is not valid in MCMC
  if (omegas_ != omegas_) {
    // std::cout << "evalOmegas: resetting omegas_." << std::endl;
    omegas_ = omegasInit_;
  }
  int nIter;
  double maxErr;
  VectorXd lGq(nObs2_);
  double logelOld = logEL();
  double logel = logelOld;
  int ii;
  for(ii=0; ii<maxIter_; ii++) {
    // E-step:
    evalWeights(); // assigns weights according to epsilons
    // M-step:
    lambdaNR(nIter, maxErr);
    lGq = ((lambdaNew_.transpose() * G_.block(0,0,nEqs_,nObs2_)).array() + weights_.head(nObs2_).sum()).transpose();
    omegas_.head(nObs2_).array() = weights_.head(nObs2_).array() / lGq.array();
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

template<typename ELModel>
inline void el::InnerELC<ELModel>::setLambda(const Ref<const VectorXd>& lambda) {
  lambdaNew_ = lambda;
}

// template<typename ELModel>
// inline void el::InnerELC<ELModel>::setWeights(const Ref<const VectorXd>& weights) {
//     weights_ = weights; 
// }
template<typename ELModel>
inline void el::InnerELC<ELModel>::setWeights(const Ref<const VectorXd>& weights) {
  weights_.head(nObs_) = weights; 
}

// template<typename ELModel>
// inline void el::InnerELC<ELModel>::setOmegas(const Ref<const VectorXd>& omegas) {
//   // nObs = _omegas.size(); // TODO: where to set nObs
//   omegasInit_. = omegas; // new
//   omegas_ = omegas; 
// }
template<typename ELModel>
inline void el::InnerELC<ELModel>::setOmegas(const Ref<const VectorXd>& omegas) {
  // nObs = _omegas.size(); // TODO: where to set nObs
  omegasInit_.head(nObs_) = omegas; // new
  omegas_.head(nObs_) = omegas; 
}

// template<typename ELModel>
// inline void el::InnerELC<ELModel>::setDeltas(const Ref<const VectorXd>& deltas) {
//   psots_ = VectorXd::Zero(nObs_); // TODO: initialization maybe should be done at better places
//   deltas_ = deltas; 
// }
template<typename ELModel>
inline void el::InnerELC<ELModel>::setDeltas(const Ref<const VectorXd>& deltas) {
  psots_ = VectorXd::Zero(nObs2_); // TODO: initialization maybe should be done at better places
  deltas_.head(nObs_) = deltas; 
}

// template<typename ELModel>
// inline VectorXd el::InnerELC<ELModel>::getWeights() {
//   return(weights_);
// }
template<typename ELModel>
inline VectorXd el::InnerELC<ELModel>::getWeights() {
    return(weights_.head(nObs2_));
}

// template<typename ELModel>
// inline VectorXd el::InnerELC<ELModel>::getOmegas() {
//   return(omegas_);
// }
template<typename ELModel>
inline VectorXd el::InnerELC<ELModel>::getOmegas() {
  return(omegas_.head(nObs2_));
}

// template<typename ELModel>
// inline VectorXd el::InnerELC<ELModel>::getEpsilons() {
//   return(epsilons_);
// }
template<typename ELModel>
inline VectorXd el::InnerELC<ELModel>::getEpsilons() {
  return(epsilons_.head(nObs2_));
}

// Note: omegas must have been assigned before calling 
// i.e., in InnerElcExports.cpp, assign omegas first 
// template<typename ELModel>
// inline double el::InnerELC<ELModel>::logEL() {
//   // evalOmegas(maxIter,relTol); 
//   if (omegas_ != omegas_) return -INFINITY; // (NaN is not equal to themselves)
//   else if ((omegas_.array() < -1e-10/omegas_.size()).any()) return -INFINITY;
//   else {
//     omegas_ = omegas_.array().abs();
//     VectorXd psos(nObs_); 
//     for (int ii=0; ii<nObs_; ii++) {
//       psos(ii) = evalPsos(ii);
//     }
//     ArrayXd logel(nObs_);
//     for (int jj=0; jj<nObs_; jj++) {
//       if (deltas_(jj) == 1) logel(jj) = log(omegas_(jj));
//       else  logel(jj) = log(psos(jj));
//     }
//     return(logel.sum());
//     // return((deltas.array()*omegas.array().log()
//     //           + (1-deltas.array())*psos.array().log()).sum());
//   }
// }
template<typename ELModel>
inline double el::InnerELC<ELModel>::logEL() {
  // evalOmegas(maxIter,relTol); 
  if (omegas_ != omegas_) return -INFINITY; // (NaN is not equal to themselves)
  else if ((omegas_.array() < -1e-10/omegas_.size()).any()) return -INFINITY;
  else {
    if (support_) {
      omegas_.tail(1)(0) = -INFINITY;
      deltas_.tail(1)(0) = 0;
    }
    omegas_ = omegas_.array().abs();
    VectorXd psos(nObs2_); 
    for (int ii=0; ii<nObs2_; ii++) {
      psos(ii) = evalPsos(ii);
    }
    ArrayXd logel(nObs2_);
    for (int jj=0; jj<nObs2_; jj++) {
      if (deltas_(jj) == 1) logel(jj) = log(omegas_(jj));
      else  logel(jj) = log(psos(jj));
    }
    return(logel.sum());
    // return((deltas.array()*omegas.array().log()
    //           + (1-deltas.array())*psos.array().log()).sum());
  }
}

// set function for G matrix
template<typename ELModel>
inline void el::InnerELC<ELModel>::setG(const Ref<const MatrixXd>& G) {
  G_.block(0,0,nEqs_,nObs2_) = G;
  // G_ = G; 
}

// get function for G matrix
template<typename ELModel>
inline MatrixXd el::InnerELC<ELModel>::getG() {
  // return(G_);
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

// returns partial sum of omegas according to the epsilons(jj) that are larger  
//   than epsilons(ii)
// Note: epsOrd must have been assigned, i.e. have called sort_inds(epsilons); 
// Note: s must be nonnegative
// template<typename ELModel>
// inline double el::InnerELC<ELModel>::evalPsosSmooth(const int ii, const double s) {
//   double psos_smooth = 0;
//   for (int jj=0; jj<nObs_; jj++) {
//     psos_smooth += ind_smooth(epsilons_(ii)-epsilons_(jj),s)*omegas_(jj);
//   }
//   return(psos_smooth);
//   // VectorXd ss = ArrayXd::Zero(omegas.size()) + s;
//   // return (ind_smooth(epsilons[ii]-epsilons.array(),ss).array()*omegas.array()).sum();
// }
template<typename ELModel>
inline double el::InnerELC<ELModel>::evalPsosSmooth(const int ii, const double s) {
  double psos_smooth = 0;
  for (int jj=0; jj<nObs2_; jj++) {
    psos_smooth += ind_smooth(epsilons_(ii)-epsilons_(jj),s)*omegas_(jj);
  }
  return(psos_smooth);
  // VectorXd ss = ArrayXd::Zero(omegas.size()) + s;
  // return (ind_smooth(epsilons[ii]-epsilons.array(),ss).array()*omegas.array()).sum();
}

// template<typename ELModel>
// inline double el::InnerELC<ELModel>::logELSmooth(const double s) {
//   if (omegas_ != omegas_) return -INFINITY; // (NaN is not equal to themselves)
//   else if ((omegas_.array() < -1e-10/omegas_.size()).any()) return -INFINITY;
//   else {
//     omegas_ = omegas_.array().abs();
//     VectorXd psos(nObs_); 
//     for (int ii=0; ii<nObs_; ii++) {
//       psos(ii) = evalPsosSmooth(ii,s);
//     }
//     ArrayXd logel(nObs_);
//     for (int jj=0; jj<nObs_; jj++) {
//       if (deltas_(jj) == 1) logel(jj) = log(omegas_(jj));
//       else  logel(jj) = log(psos(jj));
//     }
//     return(logel.sum());
//   }
// }
template<typename ELModel>
inline double el::InnerELC<ELModel>::logELSmooth(const double s) {
  if (omegas_ != omegas_) return -INFINITY; // (NaN is not equal to themselves)
  else if ((omegas_.array() < -1e-10/omegas_.size()).any()) return -INFINITY;
  else {
    if (support_) {
      omegas_.tail(1)(0) = -INFINITY;
      deltas_.tail(1)(0) = 0;
    }
    omegas_ = omegas_.array().abs();
    VectorXd psos(nObs2_); 
    for (int ii=0; ii<nObs2_; ii++) {
      psos(ii) = evalPsosSmooth(ii,s);
    }
    ArrayXd logel(nObs2_);
    for (int jj=0; jj<nObs2_; jj++) {
      if (deltas_(jj) == 1) logel(jj) = log(omegas_(jj));
      else  logel(jj) = log(psos(jj));
    }
    return(logel.sum());
  }
}

// template<typename ELModel>
// inline void el::InnerELC<ELModel>::evalWeightsSmooth(const double s) {
//   psots_.fill(0.0); // TODO initialize psots somewhere else???
//   VectorXd psoss = VectorXd::Zero(nObs_);
//   for (int ii=0; ii<nObs_; ii++) {
//     psoss(ii) = evalPsosSmooth(ii,s);
//   }
//   for (int jj=0; jj<nObs_; jj++) {
//     for (int kk=0; kk<nObs_; kk++) {
//       psots_(jj) += (1-deltas_(kk))*ind_smooth(epsilons_(kk)-epsilons_(jj),s)*omegas_(jj)/psoss(kk);
//     }
//   }
//   weights_ = deltas_.array()+psots_.array();
// }
template<typename ELModel>
inline void el::InnerELC<ELModel>::evalWeightsSmooth(const double s) {
  psots_.fill(0.0); // TODO initialize psots somewhere else???
  VectorXd psoss = VectorXd::Zero(nObs2_);
  for (int ii=0; ii<nObs2_; ii++) {
    psoss(ii) = evalPsosSmooth(ii,s);
  }
  for (int jj=0; jj<nObs2_; jj++) {
    for (int kk=0; kk<nObs2_; kk++) {
      psots_(jj) += (1-deltas_(kk))*ind_smooth(epsilons_(kk)-epsilons_(jj),s)*omegas_(jj)/psoss(kk);
    }
  }
  weights_ = deltas_.array()+psots_.array();
}

// template<typename ELModel>
// inline void el::InnerELC<ELModel>::evalOmegasSmooth(const double s) {
//   if (omegas_ != omegas_) {
//     // std::cout << "evalOmegas: resetting omegas." << std::endl;
//     omegas_ = omegasInit_;
//   }
//   int ii;
//   int nIter;
//   double maxErr;
//   VectorXd lGq;
//   double logelOld = logELSmooth(s);
//   double logel = logelOld;
//   for(ii=0; ii<maxIter_; ii++) {
//     // E-step:
//     evalWeightsSmooth(s); // assigns weights according to epsilons
//     // M-step:
//     lambdaNR(nIter, maxErr);
//     lGq = ((lambdaNew_.transpose() * G_).array() + weights_.sum()).transpose();
//     omegas_.array() = weights_.array() / lGq.array();
//     omegas_.array() = omegas_.array() / (omegas_.array().sum()); // normalize
//     logel = logELSmooth(s);
//     absErr_ = abs(logel-logelOld);
//     if (absErr_ < absTol_) break;
//     logelOld = logel;
//   }
//   nIter = ii; 
//   if (nIter == maxIter_ && absErr_ > absTol_) {
//     // TODO: maybe should assign nan elsewhere 
//     for (int ii=0; ii<nObs_; ii++) {
//       omegas_(ii) = std::numeric_limits<double>::quiet_NaN();
//     }
//   }
//   return;
// }
template<typename ELModel>
inline void el::InnerELC<ELModel>::evalOmegasSmooth(const double s) {
  if (omegas_ != omegas_) {
    // std::cout << "evalOmegas: resetting omegas." << std::endl;
    omegas_ = omegasInit_;
  }
  int nIter;
  double maxErr;
  VectorXd lGq(nObs2_);
  double logelOld = logELSmooth(s);
  double logel = logelOld;
  int ii;
  for(ii=0; ii<maxIter_; ii++) {
    // E-step:
    evalWeightsSmooth(s); // assigns weights according to epsilons
    // M-step:
    lambdaNR(nIter, maxErr);
    lGq = ((lambdaNew_.transpose() * G_.block(0,0,nEqs_,nObs2_)).array() + weights_.head(nObs2_).sum()).transpose();
    omegas_.head(nObs2_).array() = weights_.head(nObs2_).array() / lGq.array();
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

#endif
