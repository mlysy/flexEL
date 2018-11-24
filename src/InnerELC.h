/**
 * @file InnerELC.h
 * @brief Inner optimization routine for empirial likelihood problems under right-censoring.
 */

#ifndef INNERELC_h
#define INNERELC_h

#include <math.h>
#include <Rmath.h> // for random number 
#include <cmath>  // for abs on scalars
#include "SortOrder.h"
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;
#include "IndSmooth.h" // for smoothed indicator function
#include "BlockOuter.h"
// #include "MwgAdapt.h" // for adaptive mcmc

// [[Rcpp::depends(RcppEigen)]]

template <typename ELModel>
class InnerELC : public ELModel { 
private:
  
  // required members in ELModel
  using ELModel::G;
  using ELModel::nObs;
  using ELModel::nEqs;
  
  // placeholders for EM algorithm
  VectorXd deltas; /**< censoring indicators */
  VectorXd weights; /**< weights in weighted maximum log EL */
  VectorXd omegas; /**< empirical distribution */
  VectorXd omegasInit; /**< initial value for omegas, in case of reset */
  VectorXd omegasps; /**< partial sum of omegas */
  VectorXd epsilons; /**< residuals used for ordering in EM */
  VectorXi epsOrd; /**< order of epsilons */
  VectorXd psots; /**< partial sum of omegatildas */

  // placeholders for lambdaNR
  VectorXd lambdaOld;
  VectorXd lambdaNew;
  MatrixXd GGt;
  VectorXd Glambda;
  ArrayXd Gl11;
  VectorXd Q1;
  MatrixXd Q2;
  LDLT<MatrixXd> Q2ldlt;
  VectorXd rho;
  VectorXd relErr; /**< relative difference between two vectors (element-wise absolute difference divided by sum) */
  double absErr; /**< absolute difference between two scalars */
  
  // tolerance values for lambdaNR and EM
  double relTol;
  double absTol;
  int maxIter;
  
  // maximum relative error in lambda: same for cens / non-cens
  double maxRelErr(const Ref<const VectorXd>& lambdaNew,
                   const Ref<const VectorXd>& lambdaOld);
  
  /**
   * @brief      Helper function for evalWeights: calculate partial sum of omegas
   * @return     Partial sum of omegas_jj s.t. eps_jj >= eps_ii
   */
  double evalPsos(const int ii);
  
public:
  
  /**
   * @brief Default constructor for InnerELC.
   */
  InnerELC(); // default ctor
  
  /**
   * @brief Constructor for InnerELC with number of observations only as input for memory allocation.
   * @param _nObs    Number of observations.
   */
  InnerELC(int _nObs);
  
  /**
   * @brief Constructor for InnerELC with dimensions of G matrix as inputs for memory allocation.
   * @param _nObs    Number of observations.
   * @param _nEqs    Number of estimating equations.
   */
  InnerELC(int _nObs, int _nEqs);
  
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
  
  // Newton-Raphson algorithm
  /**
   * @brief Set tolerance values for NR and EM (EM uses absTol to determine the convergence of omegas).
   */
  void setTol(const int& _maxIter, const double& _relTol, const double& _absTol);
  
  /**
   * @brief Set tolerance values for NR.
   */
  void setTol(const int& _maxIter, const double& _relTol);
  
  /**
   * @brief Find the optimal lambda by a Newton-Raphson algorithm (with right-censored EL).
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
  void setDeltas(const Ref<const VectorXd>& _deltas);
  void setLambda(const Ref<const VectorXd>& _lambda); // assigned to lambdaNew
  void setWeights(const Ref<const VectorXd>& _weights); 
  void setOmegas(const Ref<const VectorXd>& _omegas);
  void setEpsilons(const Ref<const VectorXd>& _epsilons);
  VectorXd getLambda(); 
  VectorXd getWeights(); 
  VectorXd getOmegas(); 
  VectorXd getEpsilons();
  
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

// default ctor
template<typename ELModel>
inline InnerELC<ELModel>::InnerELC() {}

// ctor with dimensions as input
template<typename ELModel>
inline InnerELC<ELModel>::InnerELC(int _nObs) {
  nObs = _nObs;
  omegas = VectorXd::Zero(nObs).array() + 1.0/(double)nObs; // Initialize omegas to 1/nObs
}

// ctor with dimensions as input
template<typename ELModel>
inline InnerELC<ELModel>::InnerELC(int _nObs, int _nEqs): ELModel(_nObs, _nEqs){
  // Initialize omegas to 1/nObs
  omegas = VectorXd::Zero(nObs).array() + 1.0/(double)nObs; 
  // Newton-Raphson initialization
  GGt = MatrixXd::Zero(nEqs,nObs*nEqs);
  lambdaOld = VectorXd::Zero(nEqs); // Initialize to all 0's
  lambdaNew = VectorXd::Zero(nEqs);
  Q1 = VectorXd::Zero(nEqs);
  Q2 = MatrixXd::Zero(nEqs,nEqs);
  Glambda = VectorXd::Zero(nObs);
  Gl11 = ArrayXd::Zero(nObs);
  rho = VectorXd::Zero(nObs);
  relErr = VectorXd::Zero(nEqs);
  Q2ldlt.compute(MatrixXd::Identity(nEqs,nEqs));
}

// set tolerance for NR
template<typename ELModel>
inline void InnerELC<ELModel>::setTol(const int& _maxIter, const double& _relTol) {
  maxIter = _maxIter;
  relTol = _relTol;
}

// set tolerance for NR and EM
template<typename ELModel>
inline void InnerELC<ELModel>::setTol(const int& _maxIter,
                                      const double& _relTol, const double& _absTol) {
  maxIter = _maxIter;
  relTol = _relTol;
  absTol = _absTol;
}

// logsharp
template<typename ELModel>
inline double InnerELC<ELModel>::logsharp(double x, double q) {
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
inline double InnerELC<ELModel>::logsharp1(double x, double q) {
    if(x >= q) {
        return(1.0/x);
    } else {
        return(-1.0/(q*q)*x + 2.0/q);
    }
}

// d^2 logsharp(x)/dx^2
template<typename ELModel>
inline double InnerELC<ELModel>::logsharp2(double x, double q) {
    if(x >= q) {
        return(-1.0/(x*x));
    } else {
        return(-1.0/(q*q));
    }
}

// maximum relative error in lambda
template<typename ELModel>
inline double InnerELC<ELModel>::maxRelErr(const Ref<const VectorXd>& lambdaNew,
                                           const Ref<const VectorXd>& lambdaOld) {
  // TODO: added for numerical stability, what is a good tolerance to use ?
  if ((lambdaNew - lambdaOld).array().abs().maxCoeff() < 1e-10) return(0);
  relErr = ((lambdaNew - lambdaOld).array() / (lambdaNew + lambdaOld).array()).abs();
  return(relErr.maxCoeff());
}

// Newton-Raphson algorithm
template<typename ELModel>
inline void InnerELC<ELModel>::lambdaNR(int& nIter, double& maxErr) {
  // prevent bad starting values
  lambdaNew.fill(0.0);
  lambdaOld.fill(0.0);
  block_outer(GGt,G);
  int ii;
  
  // To avoid numerical problem
  VectorXd weights_save = weights;
  for (int ll=0;ll<nObs; ll++) {
    if (weights_save(ll) < relTol*nObs) weights_save(ll) = 0;
  }
  
  // newton-raphson loop
  for(ii=0; ii<maxIter; ii++) {
    // Q1 and Q2
    Glambda.noalias() = lambdaOld.transpose() * G;
    Glambda = weights.sum() + Glambda.array();
    Q2.fill(0.0); 
    for(int jj=0; jj<nObs; jj++) {
      rho(jj) = logsharp1(Glambda(jj), weights(jj));
      
      // To avoid numerical problem
      if (weights(jj) > relTol*nObs) {
        Q2 += weights(jj) * logsharp2(Glambda(jj), weights(jj)) * GGt.block(0,jj*nEqs,nEqs,nEqs);
      }
    }
    
    // To avoid numerical problem
    for (int kk=0; kk<nObs; kk++) {
      if (isinf(rho(kk)) && weights_save(kk) == 0) rho(kk) = 0;
    }
    
    Q1 = G * (rho.array()*weights_save.array()).matrix();
    // update lambda
    Q2ldlt.compute(Q2);
    lambdaNew.noalias() = lambdaOld - Q2ldlt.solve(Q1);
    maxErr = maxRelErr(lambdaNew, lambdaOld); // maximum relative error
    if (maxErr < relTol) {
        break;
    }
    lambdaOld = lambdaNew; // complete cycle
  }
  nIter = ii; // output lambda and also nIter and maxErr
  return;
}

template<typename ELModel>
inline VectorXd InnerELC<ELModel>::getLambda() {
  return(lambdaNew);
}

// returns partial sum of omegas according to the epsilons(jj) that are larger  
//   than epsilons(ii)
// Note: epsOrd must have been assigned, i.e. have called sort_inds(epsilons); 
template<typename ELModel>
inline double InnerELC<ELModel>::evalPsos(const int ii) {
  double psos = 0;
  int kk;
  for (int jj=nObs-1; jj>=0; jj--) {
    // Note: epsOrd corresponds to epsilons acesndingly 
    kk = epsOrd(jj); // kk is the index of the jj-th largest epsilon
    psos += omegas(kk);
    if (kk == ii) break; // until (backwardsly) reaching ii-th largest epsilon 
  }
  return psos;
}

template<typename ELModel>
inline void InnerELC<ELModel>::setEpsilons(const Ref<const VectorXd>& _epsilons) {
  epsilons = _epsilons;
  epsOrd = sort_inds(epsilons); 
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
template<typename ELModel>
inline void InnerELC<ELModel>::evalWeights() {
  // find the indices for increasing order of epsilons 
  psots.fill(0.0);
  int kk;
  double psos;
  for (int ii=0; ii<nObs; ii++) {
    for (int jj=0; jj<nObs; jj++) {
      kk = epsOrd(jj);
      if (deltas(kk) == 0) {
        psos = evalPsos(kk);
        // to prevent dividing by 0
        if (abs(psos) >= 1e-10) psots(ii) += omegas(ii)/psos;
        else if (omegas(ii) >= 1e-10 && evalPsos(kk) < 1e-10) {
          // TODO: this means a problem
          std::cout << "evalWeights: dividing by 0 problem." << std::endl;
        }
      }
      if (kk == ii) break;
    }
  }
  // assigned weights are still in the order of original data
  weights.array() = deltas.array() + psots.array();
}

// exported as `omega.hat.EM`
// Note: convergence of EM depends on the value of logEL
// Note: epsilons must have been assigned
// Note: maxIter and relTol are the same as in lambdaNR, and should have been assigned
template<typename ELModel>
inline void InnerELC<ELModel>::evalOmegas() {
  // if (omegas != omegas) return; // if initial value nan, stop
  // Need to have a valid starting value if the last one is not valid in MCMC
  if (omegas != omegas) {
    // std::cout << "evalOmegas: resetting omegas." << std::endl;
    omegas = omegasInit;
  }
  int ii;
  int nIter;
  double maxErr;
  VectorXd lGq;
  double logelOld = logEL();
  double logel = logelOld;
  for(ii=0; ii<maxIter; ii++) {
    // E-step:
    evalWeights(); // assigns weights according to epsilons
    // M-step:
    lambdaNR(nIter, maxErr);
    lGq = ((lambdaNew.transpose() * G).array() + weights.sum()).transpose();
    omegas.array() = weights.array() / lGq.array();
    omegas.array() = omegas.array() / (omegas.array().sum()); // normalize
    logel = logEL();
    maxErr = abs(logel-logelOld); // absolute error in log EL
    if (maxErr < absTol) break;
    logelOld = logel;
  }
  nIter = ii; 
  if (nIter == maxIter && maxErr > absTol) {
    // TODO: maybe should assign nan elsewhere 
    for (int ii=0; ii<nObs; ii++) {
      omegas(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  return;
}

template<typename ELModel>
inline void InnerELC<ELModel>::setLambda(const Ref<const VectorXd>& _lambda) {
  lambdaNew = _lambda;
}

template<typename ELModel>
inline void InnerELC<ELModel>::setWeights(const Ref<const VectorXd>& _weights) {
    weights = _weights; 
}

template<typename ELModel>
inline void InnerELC<ELModel>::setOmegas(const Ref<const VectorXd>& _omegas) {
  // nObs = _omegas.size(); // TODO: where to set nObs
  omegasInit = _omegas; // new
  omegas = _omegas; 
}

template<typename ELModel>
inline void InnerELC<ELModel>::setDeltas(const Ref<const VectorXd>& _deltas) {
  psots = VectorXd::Zero(nObs); // TODO: initialization maybe should be done at better places
  deltas = _deltas; 
}

template<typename ELModel>
inline VectorXd InnerELC<ELModel>::getWeights() {
    return(weights);
}

template<typename ELModel>
inline VectorXd InnerELC<ELModel>::getOmegas() {
  return(omegas);
}

template<typename ELModel>
inline VectorXd InnerELC<ELModel>::getEpsilons() {
  return(epsilons);
}

// Note: omegas must have been assigned before calling 
// i.e., in InnerElcExports.cpp, assign omegas first 
template<typename ELModel>
inline double InnerELC<ELModel>::logEL() {
  // evalOmegas(maxIter,relTol); 
  if (omegas != omegas) return -INFINITY; // (NaN is not equal to themselves)
  else if ((omegas.array() < -1e-10/omegas.size()).any()) return -INFINITY;
  else {
    omegas = omegas.array().abs();
    VectorXd psos(nObs); 
    for (int ii=0; ii<nObs; ii++) {
      psos(ii) = evalPsos(ii);
    }
    ArrayXd logel(nObs);
    for (int jj=0; jj<nObs; jj++) {
      if (deltas(jj) == 1) logel(jj) = log(omegas(jj));
      else  logel(jj) = log(psos(jj));
    }
    return(logel.sum());
    // return((deltas.array()*omegas.array().log()
    //           + (1-deltas.array())*psos.array().log()).sum());
  }
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
template<typename ELModel>
inline double InnerELC<ELModel>::evalPsosSmooth(const int ii, const double s) {
  double psos_smooth = 0;
  for (int jj=0; jj<nObs; jj++) {
    psos_smooth += ind_smooth(epsilons(ii)-epsilons(jj),s)*omegas(jj);
  }
  return(psos_smooth);
  // VectorXd ss = ArrayXd::Zero(omegas.size()) + s;
  // return (ind_smooth(epsilons[ii]-epsilons.array(),ss).array()*omegas.array()).sum();
}

template<typename ELModel>
inline double InnerELC<ELModel>::logELSmooth(const double s) {
  if (omegas != omegas) return -INFINITY; // (NaN is not equal to themselves)
  else if ((omegas.array() < -1e-10/omegas.size()).any()) return -INFINITY;
  else {
    omegas = omegas.array().abs();
    VectorXd psos(nObs); 
    for (int ii=0; ii<nObs; ii++) {
      psos(ii) = evalPsosSmooth(ii,s);
    }
    ArrayXd logel(nObs);
    for (int jj=0; jj<nObs; jj++) {
      if (deltas(jj) == 1) logel(jj) = log(omegas(jj));
      else  logel(jj) = log(psos(jj));
    }
    return(logel.sum());
  }
}

template<typename ELModel>
inline void InnerELC<ELModel>::evalWeightsSmooth(const double s) {
  psots.fill(0.0); // TODO initialize psots somewhere else???
  VectorXd psoss = VectorXd::Zero(nObs);
  for (int ii=0; ii<nObs; ii++) {
    psoss(ii) = evalPsosSmooth(ii,s);
  }
  for (int jj=0; jj<nObs; jj++) {
    for (int kk=0; kk<nObs; kk++) {
      psots(jj) += (1-deltas(kk))*ind_smooth(epsilons(kk)-epsilons(jj),s)*omegas(jj)/psoss(kk);
    }
  }
  weights = deltas.array()+psots.array();
}

template<typename ELModel>
inline void InnerELC<ELModel>::evalOmegasSmooth(const double s) {
  if (omegas != omegas) {
    // std::cout << "evalOmegas: resetting omegas." << std::endl;
    omegas = omegasInit;
  }
  int ii;
  int nIter;
  double maxErr;
  VectorXd lGq;
  double logelOld = logELSmooth(s);
  double logel = logelOld;
  for(ii=0; ii<maxIter; ii++) {
    // E-step:
    evalWeightsSmooth(s); // assigns weights according to epsilons
    // M-step:
    lambdaNR(nIter, maxErr);
    lGq = ((lambdaNew.transpose() * G).array() + weights.sum()).transpose();
    omegas.array() = weights.array() / lGq.array();
    omegas.array() = omegas.array() / (omegas.array().sum()); // normalize
    logel = logELSmooth(s);
    absErr = abs(logel-logelOld);
    if (absErr < absTol) break;
    logelOld = logel;
  }
  nIter = ii; 
  if (nIter == maxIter && absErr > absTol) {
    // TODO: maybe should assign nan elsewhere 
    for (int ii=0; ii<nObs; ii++) {
      omegas(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  return;
}

#endif
