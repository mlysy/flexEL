#ifndef INNERELC_h
#define INNERELC_h

// class for computing inner optimization of EL likelihood with censoring

// using namespace std;
#include <math.h>
#include <Rmath.h> // for random number 
#include <cmath>  // for abs on scalars
// #include <math.h> // for isnan
#include "SortOrder.h"
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;
#include "MwgAdapt.h" // for adaptive mcmc
// [[Rcpp::depends(RcppEigen)]]

// main class begins here
template <typename ELModel>
class InnerELC : public ELModel { 
private:
  using ELModel::nObs;
  using ELModel::nEqs;
  using ELModel::nBet;
  using ELModel::nGam;
  using ELModel::nQts;
  // need access to y, X, Z in evalEpsilons
  using ELModel::y;
  using ELModel::X;
  using ELModel::Z;
  using ELModel::G;
  VectorXd deltas; 
  VectorXd weights; 
  VectorXd omegasInit; // new
  VectorXd omegas; 
  VectorXd omegasps; // partial sum of omegas 
  // constants for logsharp calculations
  // double trunc, aa, bb, cc; 
  // temporary storage for Newton-Raphson
  VectorXd lambdaOld;
  VectorXd lambdaNew;
  MatrixXd GGt;
  VectorXd Glambda;
  ArrayXd Gl11;
  VectorXd Q1;
  MatrixXd Q2;
  LDLT<MatrixXd> Q2ldlt;
  VectorXd rho;
  // VectorXd relErr;
  double relErr; // new
  double absErr; // new
  // tolerance for Newton-Raphson lambdaNR and evalOmegas (New)
  int maxIter;
  double relTol;
  double absTol;
  VectorXd epsilons; // error term used to order omegas in evalWeights
  // vector<size_t> epsOrd; // C vector of indicies of ordered epsilons
  VectorXi epsOrd; // vector of indicies of ordered epsilons
  VectorXd psots; // partial sum of omegatildas
  // columnwise outer product (see below)
  void blockOuter(void);
  // maximum relative error in lambda: same for cens / non-cens
  double maxRelErr(const Ref<const VectorXd>& lambdaNew,
                   const Ref<const VectorXd>& lambdaOld);
  double maxRelErr(const double& valnew,
                   const double& valold);
  // helper function for evalWeights: calculate partial sum of omegas
  // partial sum of omegas_jj s.t. eps_jj >= eps_ii
  double evalPsos(const int ii);
public:
  // constructor for regression-like problems
  // InnerELC(const Ref<const VectorXd>& _y, const Ref<const MatrixXd>& _X, 
  //          const Ref<const VectorXd>& _deltas,
  //          void* params);
  InnerELC(); // default ctor
  void setData(const Ref<const VectorXd>& _y, 
               const Ref<const MatrixXd>& _X, 
               const Ref<const VectorXd>& _deltas, 
               void* params); 
  void setData(const Ref<const VectorXd>& _y, 
               const Ref<const MatrixXd>& _X, 
               const Ref<const MatrixXd>& _Z,
               const Ref<const VectorXd>& _deltas, 
               void* params); 
  // logsharp and its derivatives
  double logsharp(double x, double q);
  double logsharp1(double x, double q);
  double logsharp2(double x, double q);
  // Newton-Raphson algorithm
  void setTol(const int& _maxIter, const double& _relTol, const double& _absTol);
  void lambdaNR(int& nIter, double& maxErr);
  // void lambdaNR(int& nIter, double& maxErr, 
  //               int maxIter, double relTol);
  // eval functions 
  void evalEpsilons(const Ref<const VectorXd>& beta); // for location model 
  void evalEpsilons(const Ref<const VectorXd>& beta,
                    const Ref<const VectorXd>& gamma,
                    const double& sig2); // for location-scale model
  void evalWeights(); // calculate weights according to epsilons 
  // void evalOmegas(int& nIter, double& maxErr, int maxIter, double relTol);
  // void evalOmegas(int maxIter, double relTol);
  void evalOmegas();
  // log empirical likelihood calculation with given omegas and G
  double logEL();
  // double logEL(int maxIter, double relTol); 
  // set and get functions 
  void setLambda(const Ref<const VectorXd>& _lambda); // assigned to lambdaNew
  void setWeights(const Ref<const VectorXd>& _weights); 
  void setOmegas(const Ref<const VectorXd>& _omegas);
  void setEpsilons(const Ref<const VectorXd>& _epsilons);
  void setDeltas(const Ref<const VectorXd>& _deltas);
  VectorXd getLambda(); 
  VectorXd getWeights(); 
  VectorXd getOmegas(); 
  VectorXd getEpsilons();
  // posterior sampler
  MatrixXd postSample(int nsamples,int nburn, VectorXd betaInit, 
                      const Ref<const VectorXd>& sigs, 
                      VectorXd &RvDoMcmc,VectorXd &paccept);
  void mwgStep(VectorXd &thetaCur, const int &idx, const double &mwgsd,
               bool &accept, double &logELCur);
  MatrixXd postSampleAdapt(int nsamples, int nburn, VectorXd thetaInit,
                           double *mwgSd, VectorXd &rvDoMcmc, bool *doAdapt, 
                           VectorXd &paccept);
  // smooth functions
  VectorXd indSmooth(VectorXd x, VectorXd s);
  double evalPsosSmooth(const int ii, const double s); // any way to test while set as private?
  double logELSmooth(const double s);
  void evalWeightsSmooth(const double s);
  void evalOmegasSmooth(const double s);
};

/*
// constructor for mean regression (without alpha)
template<typename ELModel>
inline InnerELC<ELModel>::InnerELC(const Ref<const VectorXd>& _y,
                                   const Ref<const MatrixXd>& _X,
                                   const Ref<const VectorXd>& _deltas,
                                   void* params) : ELModel(_y, _X, params) {
    // std::cout << nObs << std::endl;
    // std::cout << nEqs << std::endl;
    epsilons = VectorXd::Zero(nObs);
    psots = VectorXd::Zero(nObs);
    // epsOrd = vector<size_t>(nObs); // Note: epsOrd is a C vector not Eigen VectorXd
    epsOrd = VectorXi::Zero(nObs);
    // Newton-Raphson initialization
    deltas = _deltas;
    omegas = VectorXd::Zero(nObs).array() + 1.0/(double)nObs; // Initialize to 1/nObs
    omegasps = VectorXd::Zero(nObs);
    W = MatrixXd::Zero(nObs,nObs);
    weights = VectorXd::Zero(nObs); // Initialize with the current omegas?
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
 */

// default ctor
template<typename ELModel>
inline InnerELC<ELModel>::InnerELC() {}

// set data with default ctor
template<typename ELModel>
inline void InnerELC<ELModel>::setData(const Ref<const VectorXd>& _y,
                                       const Ref<const MatrixXd>& _X,
                                       const Ref<const VectorXd>& _deltas,
                                       void* params) {
    ELModel::setData(_y,_X,params); 
    epsilons = VectorXd::Zero(nObs);
    psots = VectorXd::Zero(nObs);
    epsOrd = VectorXi::Zero(nObs);
    // Newton-Raphson initialization
    deltas = _deltas;
    omegas = VectorXd::Zero(nObs).array() + 1.0/(double)nObs; // Initialize to 1/nObs
    omegasps = VectorXd::Zero(nObs);
    weights = VectorXd::Zero(nObs); // Initialize with the current omegas?
    GGt = MatrixXd::Zero(nEqs,nObs*nEqs);
    lambdaOld = VectorXd::Zero(nEqs); // Initialize to all 0's
    lambdaNew = VectorXd::Zero(nEqs);
    Q1 = VectorXd::Zero(nEqs);
    Q2 = MatrixXd::Zero(nEqs,nEqs);
    Glambda = VectorXd::Zero(nObs);
    Gl11 = ArrayXd::Zero(nObs);
    rho = VectorXd::Zero(nObs);
    // relErr = VectorXd::Zero(nEqs);
    Q2ldlt.compute(MatrixXd::Identity(nEqs,nEqs));
}

template<typename ELModel>
inline void InnerELC<ELModel>::setData(const Ref<const VectorXd>& _y,
                                       const Ref<const MatrixXd>& _X,
                                       const Ref<const MatrixXd>& _Z,
                                       const Ref<const VectorXd>& _deltas,
                                       void* params) {
  ELModel::setData(_y,_X,_Z,params); 
  epsilons = VectorXd::Zero(nObs);
  psots = VectorXd::Zero(nObs);
  epsOrd = VectorXi::Zero(nObs);
  // Newton-Raphson initialization
  deltas = _deltas;
  omegas = VectorXd::Zero(nObs).array() + 1.0/(double)nObs; // Initialize to 1/nObs
  omegasps = VectorXd::Zero(nObs);
  weights = VectorXd::Zero(nObs); // Initialize with the current omegas?
  GGt = MatrixXd::Zero(nEqs,nObs*nEqs);
  lambdaOld = VectorXd::Zero(nEqs); // Initialize to all 0's
  lambdaNew = VectorXd::Zero(nEqs);
  Q1 = VectorXd::Zero(nEqs);
  Q2 = MatrixXd::Zero(nEqs,nEqs);
  Glambda = VectorXd::Zero(nObs);
  Gl11 = ArrayXd::Zero(nObs);
  rho = VectorXd::Zero(nObs);
  // relErr = VectorXd::Zero(nEqs);
  Q2ldlt.compute(MatrixXd::Identity(nEqs,nEqs));
}

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

// for an (m x N) matrix G = [g1 ... gN], returns the (m x mN) matrix
// GGt = [g1 g1' ... gN gN']
template<typename ELModel>
inline void InnerELC<ELModel>::blockOuter(void) {
  // for each row of G, compute outer product and store as block
    for(int ii=0; ii<nObs; ii++) {
        GGt.block(0,ii*nEqs,nEqs,nEqs).noalias() = G.col(ii) * G.col(ii).transpose();
    }
    return;
}

// maximum relative error in lambda
template<typename ELModel>
inline double InnerELC<ELModel>::maxRelErr(const Ref<const VectorXd>& lambdaNew,
                                           const Ref<const VectorXd>& lambdaOld) {
  // TODO: have to find another way to check convergence of omegas
  // maybe write another function abs err or rel error 
  // since the omegas could be really small and thus hard to determine convergence
  
  // TODO: added for numerical stability, what is a good tolerance to use ?
  if ((lambdaNew - lambdaOld).array().abs().maxCoeff() < 1e-10) return(0);
  relErr = ((lambdaNew - lambdaOld).array() / (lambdaNew + lambdaOld).array()).abs().maxCoeff();
  // return(relErr.maxCoeff());
  return(relErr);

  // absErr = (lambdaNew - lambdaOld).array().abs().maxCoeff();
  // relErr = ((lambdaNew - lambdaOld).array() / (lambdaNew + lambdaOld).array()).abs().maxCoeff();
  // return(std::min(absErr,relErr));
}

template<typename ELModel>
inline double InnerELC<ELModel>::maxRelErr(const double& valnew,
                                           const double& valold) {
  absErr = std::abs(valnew - valold);
  relErr = std::abs((valnew - valold) / (valnew + valold));
  return(std::min(absErr,relErr));
}

// Newton-Raphson algorithm
template<typename ELModel>
inline void InnerELC<ELModel>::lambdaNR(int& nIter, double& maxErr) {
  // prevent bad starting values
  lambdaNew.fill(0.0);
  lambdaOld.fill(0.0);
  blockOuter(); // initialize GGt according to epsilons order (epsOrd)
  int ii;
  
  // TODO: to avoid numerical problem
  VectorXd weights_save = weights;
  for (int ll=0;ll<nObs; ll++) {
    if (weights_save(ll) < relTol*nObs) weights_save(ll) = 0;
  }
  
  // newton-raphson loop
  for(ii=0; ii<maxIter; ii++) {
    // std::cout << "lambdaOld = " << lambdaOld.transpose() << std::endl;
    // Q1 and Q2
    Glambda.noalias() = lambdaOld.transpose() * G;
    Glambda = weights.sum() + Glambda.array();
    Q2.fill(0.0); 
    for(int jj=0; jj<nObs; jj++) {
      rho(jj) = logsharp1(Glambda(jj), weights(jj));
      
      // TODO: to avoid numerical problem
      if (weights(jj) > relTol*nObs) {
        Q2 += weights(jj) * logsharp2(Glambda(jj), weights(jj)) * GGt.block(0,jj*nEqs,nEqs,nEqs);
      }
    }
    
    // TODO: to avoid numerical problem
    for (int kk=0; kk<nObs; kk++) {
      if (isinf(rho(kk)) && weights_save(kk) == 0) rho(kk) = 0;
    }
    
    Q1 = G * (rho.array()*weights_save.array()).matrix();
    // Q1 = G * (rho.array()*weights.array()).matrix();
    // std::cout << "Q2 = \n" << Q2 << std::endl;
    // std::cout << "Q1 = \n" << Q1.transpose() << std::endl;
    // update lambda
    Q2ldlt.compute(Q2);
    lambdaNew.noalias() = lambdaOld - Q2ldlt.solve(Q1);
    maxErr = maxRelErr(lambdaNew, lambdaOld); // maximum relative error
    // std::cout << "In lambdaNR: maxErr = " << maxErr << std::endl;
    if (maxErr < relTol) {
        break;
    }
    lambdaOld = lambdaNew; // complete cycle
  }
  nIter = ii; // output lambda and also nIter and maxErr
  // if (nIter == maxIter) {
  //   std::cout << "Q2ldlt.solve(Q1) = " << Q2ldlt.solve(Q1).transpose() << std::endl;
  //   std::cout << "lambdaNew = " << lambdaNew.transpose() << std::endl;
  // }
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
      // epsOrd corresponds to epsilons acesndingly 
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

// evaluate epsilons location model
template<typename ELModel>
inline void InnerELC<ELModel>::evalEpsilons(const Ref<const VectorXd>& beta) {
  epsilons = (y.transpose() - beta.transpose() * X).transpose(); 
  epsOrd = sort_inds(epsilons);
}

// evaluate epsilons location-scale model
template<typename ELModel>
inline void InnerELC<ELModel>::evalEpsilons(const Ref<const VectorXd>& beta,
                                            const Ref<const VectorXd>& gamma,
                                            const double& sig2) {
  if (sig2 < 0.0) {
    std::cout << "negative sig2 in evalEpsilons." << std::endl;
  }
  epsilons = (y.transpose() - beta.transpose() * X).transpose(); 
  epsilons.transpose().array() *= (-gamma.transpose()*Z).array().exp()/sqrt(sig2);
  epsOrd = sort_inds(epsilons);
}

// Note: epsilons must have been assigned
template<typename ELModel>
inline void InnerELC<ELModel>::evalWeights() {
  // std::cout << "---- In evalWeights ----" << std::endl;
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
        // psots(ii) += omegas(ii)/evalPsos(kk); // old code
      }
      if (kk == ii) break;
    }
  }
  // assigned weights are still in the order of original data
  weights.array() = deltas.array() + psots.array();
  // std::cout << "evalWeighes: weights = " << weights.transpose() << std::endl;
  // std::cout << "evalWeights: deltas = " << deltas.transpose() << std::endl;
  // std::cout << "evalWeights: epsilons = " << epsilons.transpose() << std::endl;
  // std::cout << "evalWeights: psots = " << psots.transpose() << std::endl;
}

// exported as `omega.hat.EM`
// Note: epsilons must have been assigned
// Note: maxIter and relTol are the same as in lambdaNR, and should have been assigned
template<typename ELModel>
inline void InnerELC<ELModel>::evalOmegas() {
  // inline void InnerELC<ELModel>::evalOmegas(int maxIter, double relTol) {
  // std::cout << "**** In evalOmegas ****" << std::endl;
  // std::cout << "prelim: omegas = " << omegas.transpose() << std::endl;
  // if (omegas != omegas) return; // if initial value nan, stop
  // TODO: need to have a valid starting value if the last one is not valid
  if (omegas != omegas) {
    // std::cout << "evalOmegas: resetting omegas." << std::endl;
    omegas = omegasInit;
  }
  int ii;
  int nIter;
  double maxErr;
  VectorXd lGq;
  // omegas.fill(1.0/nObs);
  // VectorXd omegasOld = omegas;
  double logelOld = logEL();
  double logel = logelOld;
  // std::cout << "Initial lambdaNew = " << lambdaNew.transpose() << std::endl;
  // std::cout << "Initial lambdaOld = " << lambdaOld.transpose() << std::endl;
  // std::cout << "Initial omegas = " << omegas.transpose() << std::endl;
  for(ii=0; ii<maxIter; ii++) {
    // E-step:
    evalWeights(); // assigns weights according to epsilons
    // std::cout << "weights = " << weights.transpose() << std::endl;
    // M-step:
    // std::cout << "lambdaOld = " << lambdaOld.transpose() << std::endl;
    // lambdaNR(nIter, maxErr, maxIter, relTol); 
    lambdaNR(nIter, maxErr);
    // std::cout << "lambdaNew = " << lambdaNew.transpose() << std::endl;
    // TODO: Check convergence of NR here ???
    // if (nIter == maxIter && maxErr > relTol) {
    //   std::cout << "evalOmegas: lambdaNRC did not converge in EM." << std::endl;
    //   std::cout << "nIter = " << nIter << std::endl;
    //   std::cout << "maxErr = " << maxErr << std::endl;
    // }
    lGq = ((lambdaNew.transpose() * G).array() + weights.sum()).transpose();
    // std::cout << "lGq = " << lGq.transpose() << std::endl;
    // std::cout << "lambdaOld = " << lambdaOld.transpose() << std::endl;
    omegas.array() = weights.array() / lGq.array();
    // std::cout << "In evalOmegas before normalize: omegas = \n" << omegas.transpose() << std::endl;
    omegas.array() = omegas.array() / (omegas.array().sum()); // normalize
    // std::cout << "In evalOmegas: omegas = " << omegas.transpose() << std::endl;
    // maxErr = maxRelErr(omegas, omegasOld); 
    logel = logEL();
    maxErr = maxRelErr(logel, logelOld);
    // std::cout << "In evalOmegas: maxErr = " << maxErr << std::endl;
    if (maxErr < relTol) break;
    // omegasOld = omegas;
    logelOld = logel;
  }
  nIter = ii; 
  // std::cout << "nIter = " << nIter << std::endl;
  if (nIter == maxIter && maxErr > relTol) {
    // std::cout << "evalOmegas not converged. maxErr = " << maxErr << std::endl;
    // std::cout << "the omegas now = " << omegas.transpose() << std::endl;
    // TODO: maybe should assign nan elsewhere 
    for (int ii=0; ii<nObs; ii++) {
      omegas(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  // std::cout << "---- End of evalOmegas ----" << std::endl;
  return;
}

/* Old code
template<typename ELModel>
inline void InnerELC<ELModel>::evalOmegas() {
// inline void InnerELC<ELModel>::evalOmegas(int maxIter, double relTol) {
  std::cout << "**** In evalOmegas ****" << std::endl;
  // Note: have to save this o.w. lambdaOld got modified
  VectorXd lambdaOldEM = lambdaOld; 
  int ii;
  int nIter;
  double maxErr;
  VectorXd lGq;
  for(ii=0; ii<maxIter; ii++) {
      // E-step:
      evalWeights(); // assigns weights according to epsilons
      // std::cout << "weights = " << weights.transpose() << std::endl;
      // M-step:
      // std::cout << "lambdaOldEM = " << lambdaOldEM.transpose() << std::endl;
      // lambdaNR(nIter, maxErr, maxIter, relTol); 
      lambdaNR(nIter, maxErr);
      // Check convergence of NR here
      if (nIter == maxIter & maxErr > relTol) {
        std::cout << "lambdaNRC did not converge in EM" << std::endl;
      }
      lGq = ((lambdaNew.transpose() * G).array() + weights.sum()).transpose();
      // std::cout << "lambdaNew = " << lambdaNew.transpose() << std::endl;
      // std::cout << "lambdaOld = " << lambdaOld.transpose() << std::endl;
      // std::cout << "lambdaOldEM = " << lambdaOldEM.transpose() << std::endl;
      omegas.array() = weights.array() / lGq.array();
      // std::cout << "In evalOmegas before normalize: omegas = \n" << omegas.transpose() << std::endl;
      omegas.array() = omegas.array() / (omegas.array().sum()); // normalize
      // std::cout << "In evalOmegas: omegas = \n" << omegas.transpose() << std::endl;
      maxErr = maxRelErr(lambdaNew, lambdaOldEM); 
      // std::cout << "In evalOmegas: maxErr = " << maxErr << std::endl;
      if (maxErr < relTol) break;
      lambdaOldEM = lambdaNew;
  }
  nIter = ii; 
  if (nIter == maxIter & maxErr > relTol) {
    // omegas = VectorXd::Zero(nObs); 
    std::cout << "evalOmegas not converged." << std::endl;
    for (int ii=0; ii<nObs; ii++) {
      omegas(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  std::cout << "last lambdaNew = " << lambdaNew.transpose() << std::endl;
  std::cout << "G = \n" << G << std::endl;
  std::cout << "lGq = \n" << lGq.transpose() << std::endl;
  std::cout << "weights = \n" << weights.transpose() << std::endl;
  std::cout << "In evalOmegas: omegas = \n" << omegas.transpose() << std::endl;
  return;
}
*/

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
  nObs = _omegas.size(); // TODO: where to set nObs
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

// smoothed indicator function 
// Note: length of x and s must be the same, and s must be nonnegative to be consistent with
// the original indicator function x <= 0.
template<typename ELModel>
inline VectorXd InnerELC<ELModel>::indSmooth(VectorXd x, VectorXd s) {
  return(1.0/(1.0+(x.array()*s.array()).exp()));
}

// returns partial sum of omegas according to the epsilons(jj) that are larger  
//   than epsilons(ii)
// Note: epsOrd must have been assigned, i.e. have called sort_inds(epsilons); 
// Note: s must be nonnegative
template<typename ELModel>
inline double InnerELC<ELModel>::evalPsosSmooth(const int ii, const double s) {
  // psos <- sum(ind.smooth_R(epsilons[ii]-epsilons,s)*omegas)
  VectorXd ss = ArrayXd::Zero(omegas.size()) + s;
  return (indSmooth(epsilons[ii]-epsilons.array(),ss).array()*omegas.array()).sum();
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
    // std::cout << "logel = " << logel.transpose() << std::endl;
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
    VectorXd ss = ArrayXd::Zero(omegas.size()) + s;
    psots(jj) = ((1-deltas.array())*indSmooth(epsilons.array()-epsilons(jj),ss).array()*omegas(jj)/psoss.array()).sum();
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
    // std::cout << "weights = " << weights.transpose() << std::endl;
    // M-step:
    // std::cout << "lambdaOld = " << lambdaOld.transpose() << std::endl;
    // lambdaNR(nIter, maxErr, maxIter, relTol); 
    lambdaNR(nIter, maxErr);
    // std::cout << "lambdaNew = " << lambdaNew.transpose() << std::endl;
    // TODO: Check convergence of NR here ???
    // if (nIter == maxIter && maxErr > relTol) {
    //   std::cout << "evalOmegas: lambdaNRC did not converge in EM." << std::endl;
    //   std::cout << "nIter = " << nIter << std::endl;
    //   std::cout << "maxErr = " << maxErr << std::endl;
    // }
    lGq = ((lambdaNew.transpose() * G).array() + weights.sum()).transpose();
    // std::cout << "lGq = " << lGq.transpose() << std::endl;
    // std::cout << "lambdaOld = " << lambdaOld.transpose() << std::endl;
    omegas.array() = weights.array() / lGq.array();
    // std::cout << "In evalOmegas before normalize: omegas = \n" << omegas.transpose() << std::endl;
    omegas.array() = omegas.array() / (omegas.array().sum()); // normalize
    // std::cout << "In evalOmegas: omegas = " << omegas.transpose() << std::endl;
    // maxErr = maxRelErr(omegas, omegasOld); 
    logel = logELSmooth(s);
    // maxErr = maxRelErr(logel, logelOld);
    absErr = abs(logel-logelOld);
    // std::cout << "In evalOmegas: maxErr = " << maxErr << std::endl;
    if (absErr < absTol) break;
    // omegasOld = omegas;
    logelOld = logel;
  }
  nIter = ii; 
  // std::cout << "nIter = " << nIter << std::endl;
  if (nIter == maxIter && absErr > absTol) {
    // std::cout << "evalOmegas not converged. maxErr = " << maxErr << std::endl;
    // std::cout << "the omegas now = " << omegas.transpose() << std::endl;
    // TODO: maybe should assign nan elsewhere 
    for (int ii=0; ii<nObs; ii++) {
      omegas(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  // std::cout << "---- End of evalOmegas ----" << std::endl;
  return;
}

#endif
