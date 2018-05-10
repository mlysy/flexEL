#ifndef INNERELC_h
#define INNERELC_h

// class for computing inner optimization of EL likelihood with censoring

// using namespace std;
#include <math.h>
#include <Rmath.h> // for random number 
#include "SortOrder.h"
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]


// main class begins here
template <typename ELModel>
class InnerELC : public ELModel { 
private:
    using ELModel::nObs;
    using ELModel::nEqs;
    using ELModel::X; // need access to X, y in evalWeights
    using ELModel::y;
    using ELModel::G;
    VectorXd deltas; 
    VectorXd weights; 
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
    VectorXd relErr;
    // tolerance for Newton-Raphson lambdaNR and evalOmegas (New)
    int maxIter;
    double relTol;
    VectorXd epsilons; // error term used to order omegas in evalWeights
    // vector<size_t> epsOrd; // C vector of indicies of ordered epsilons
    VectorXi epsOrd; // vector of indicies of ordered epsilons
    VectorXd psots; // partial sum of omegatildas
    // columnwise outer product (see below)
    void blockOuter(void);
    // maximum relative error in lambda: same for cens / non-cens
    double maxRelErr(const Ref<const VectorXd>& lambdaNew,
                     const Ref<const VectorXd>& lambdaOld);
    // helper function for evalWeights: calculate partial sum of omegas
    // partial sum of omegas_jj s.t. eps_jj >= eps_ii
    double evalPsos(const int ii);
public:
    // constructor for regression-like problems
    // InnerELC(const Ref<const VectorXd>& _y, const Ref<const MatrixXd>& _X, 
    //          const Ref<const VectorXd>& _deltas,
    //          void* params);
    InnerELC(); // default ctor
    void setData(const Ref<const VectorXd>& _y, const Ref<const MatrixXd>& _X, 
                 const Ref<const VectorXd>& _deltas, void* params); 
    // logsharp and its derivatives
    double logsharp(double x, double q);
    double logsharp1(double x, double q);
    double logsharp2(double x, double q);
    // Newton-Raphson algorithm
    void setTol(const int& _maxIter, const double& _relTol);
    void lambdaNR(int& nIter, double& maxErr);
    // void lambdaNR(int& nIter, double& maxErr, 
    //               int maxIter, double relTol);
    // eval functions 
    void evalEpsilons(const Ref<const VectorXd>& beta); // for location model 
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
    VectorXd getLambda(); 
    VectorXd getWeights(); 
    VectorXd getOmegas(); 
    // log empirical likelihood calculation 
    // double logEL(const Ref<const VectorXd>& theta, int maxIter, double relTol); 
    // posterior sampler
    MatrixXd postSample(int nsamples,int nburn,VectorXd betaInit, 
                        const Ref<const VectorXd>& sigs, 
                        VectorXd &RvDoMcmc,VectorXd &paccept);
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
    // epsOrd = vector<size_t>(nObs); // Note: epsOrd is a C vector not Eigen VectorXd
    epsOrd = VectorXi::Zero(nObs);
    // Newton-Raphson initialization
    deltas = _deltas;
    omegas = VectorXd::Zero(nObs).array() + 1.0/(double)nObs; // Initialize to 1/nObs
    omegasps = VectorXd::Zero(nObs);
    // W = MatrixXd::Zero(nObs,nObs);
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

template<typename ELModel>
inline void InnerELC<ELModel>::setTol(const int& _maxIter,
                                     const double& _relTol) {
  maxIter = _maxIter;
  relTol = _relTol;
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
  // TODO: added for numerical stability, what is a good tolerance to use ?
  if ((lambdaNew - lambdaOld).array().abs().maxCoeff() < 1e-10) return(0);
  
  // std::cout << "In maxRelErr: lambdaNew = " << lambdaNew << ", lambdaOld = " << lambdaOld << std::endl;
  relErr = ((lambdaNew - lambdaOld).array() / (lambdaNew + lambdaOld).array()).abs();
  return(relErr.maxCoeff());
}

// Newton-Raphson algorithm
template<typename ELModel>
inline void InnerELC<ELModel>::lambdaNR(int& nIter, double& maxErr) {
    int ii, jj;
    blockOuter(); // initialize GGt according to epsilons order (epsOrd)
    // newton-raphson loop
    for(ii=0; ii<maxIter; ii++) {
        // Q1 and Q2
        Glambda.noalias() = lambdaOld.transpose() * G;
        Glambda = weights.sum() + Glambda.array();
        Q2.fill(0.0); 
        for(jj=0; jj<nObs; jj++) {
            rho(jj) = logsharp1(Glambda(jj), weights(jj));
            Q2 += weights(jj) * logsharp2(Glambda(jj), weights(jj)) * GGt.block(0,jj*nEqs,nEqs,nEqs);
        }
        Q1 = G * (rho.array()*weights.array()).matrix();
        // update lambda
        Q2ldlt.compute(Q2);
        // std::cout << "Q2invQ1 = " << Q2ldlt.solve(Q1).transpose() << std::endl;
        lambdaNew.noalias() = lambdaOld - Q2ldlt.solve(Q1);
        maxErr = maxRelErr(lambdaNew, lambdaOld); // maximum relative error
        // std::cout << "maxErr = " << maxErr << std::endl;
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

// evaluate epsilons
template<typename ELModel>
inline void InnerELC<ELModel>::evalEpsilons(const Ref<const VectorXd>& beta) {
  epsilons = (y.transpose() - beta.transpose() * X).transpose(); 
}

// Note: epsilons must have been assigned
template<typename ELModel>
inline void InnerELC<ELModel>::evalWeights() {
  // std::cout << "---- In evalWeights ----" << std::endl;
  // find the indices for increasing order of epsilons 
  epsOrd = sort_inds(epsilons);  // TODO: might remove it from here
  // std::cout << "omegas = " << omegas.transpose() << std::endl;
  psots.fill(0.0);
  int kk;
  for (int ii=0; ii<nObs; ii++) {
    for (int jj=0; jj<nObs; jj++) {
      kk = epsOrd(jj);
      if (deltas(kk) == 0) {
        // std::cout << "evalPsos(kk) = " << evalPsos(kk) << std::endl;
        psots(ii) += omegas(ii)/evalPsos(kk);
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
  int ii;
  int nIter;
  double maxErr;
  VectorXd lGq;
  // omegas.fill(1.0/nObs);
  VectorXd omegasOld = omegas;
  // std::cout << "omegas = " << omegas.transpose() << std::endl;
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
    maxErr = maxRelErr(omegas, omegasOld); 
    // std::cout << "In evalOmegas: maxErr = " << maxErr << std::endl;
    if (maxErr < relTol) break;
    omegasOld = omegas;
  }
  nIter = ii; 
  if (nIter == maxIter & maxErr > relTol) {
    // omegas = VectorXd::Zero(nObs); 
    // std::cout << "evalOmegas not converged." << std::endl;
    for (int ii=0; ii<nObs; ii++) {
      omegas(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  // std::cout << "last lambdaNew = " << lambdaNew.transpose() << std::endl;
  // std::cout << "G = \n" << G << std::endl;
  // std::cout << "lGq = \n" << lGq.transpose() << std::endl;
  // std::cout << "weights = \n" << weights.transpose() << std::endl;
  // std::cout << "In evalOmegas: omegas = \n" << omegas.transpose() << std::endl;
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
    omegas = _omegas; 
}

template<typename ELModel>
inline VectorXd InnerELC<ELModel>::getWeights() {
    return(weights);
}

template<typename ELModel>
inline VectorXd InnerELC<ELModel>::getOmegas() {
    return(omegas);
}

// Note: omegas must have been assigned before calling 
// i.e., in InnerElcExports.cpp, assign omegas first 
template<typename ELModel>
inline double InnerELC<ELModel>::logEL() {
  // evalOmegas(maxIter,relTol); 
  if (omegas != omegas) return -INFINITY; // (NaN is not equal to themselves)
  else {
    VectorXd psos(nObs); 
    for (int ii=0; ii<nObs; ii++) {
      psos(ii) = evalPsos(ii);
    }
    // std::cout << "epsOrd = " << epsOrd.transpose() << std::endl;
    // std::cout << "psos = " << psos.transpose() << std::endl;
    return((deltas.array()*omegas.array().log()
              + (1-deltas.array())*psos.array().log()).sum());
  }
}

// template<typename ELModel>
// inline double InnerELC<ELModel>::logEL() {
//     // TODO: (sort it here?) 
//     // sort epsilons decendingly and record its order 
//     epsOrd = sort_inds(epsilons); 
//     // std::cout << "epsOrd = " << epsOrd.transpose() << std::endl; 
//     // check if omega is not a feasible solution with G return -Inf
//     VectorXd rhs = G * omegas; // if feasible, should be approx a vector of 0s
//     // std::cout << "rhs = " << rhs.transpose() << std::endl; 
//     // TODO: what tolarance?  
//     if (rhs.array().sum() > 0.01) return -INFINITY;
//     else {
//         VectorXd psos(nObs); 
//         for (int ii=0; ii<nObs; ii++) {
//             psos(ii) = evalPsos(ii); 
//         }
//         // std::cout << "psos = " << psos.transpose() << std::endl; 
//         return((deltas.array()*omegas.array().log()
//                + (1-deltas.array())*psos.array().log()).sum());
//     }
// }

// posterior sampler: location model, single quantile case
// Note: omegasInit needed as the starting value for EM?
template<typename ELModel>
inline MatrixXd InnerELC<ELModel>::postSample(int nsamples, int nburn,
                                              VectorXd betaInit, 
                                              const Ref<const VectorXd>& sigs,
                                              VectorXd &RvDoMcmc, VectorXd &paccept) {
  // std::cout << "--------------------------- In postSample ----------------------------" << std::endl;
  VectorXd betaOld = betaInit;
  VectorXd betaNew = betaOld;
  VectorXd betaProp = betaOld;
  int betalen = betaInit.size();
  MatrixXd beta_chain(betalen,nsamples);
  ELModel::evalG(betaOld);
  evalEpsilons(betaInit);
  // evalWeights();
  int nIter;
  double maxErr;
  // lambdaNR(nIter, maxErr);
  // if (nIter == maxIter && maxErr > relTol) {
  //   // TODO: what to do ??
  //   std::cout << "betaInit not valid." << std::endl;
  //   // return NULL;
  // }
  std::cout << "Initial omegas = " << omegas.transpose() << std::endl;
  evalOmegas(); // omegasInit should have been assigned
  std::cout << "First omegas = " << omegas.transpose() << std::endl;
  double logELOld = logEL();
  double logELProp;
  bool satisfy;
  double u;
  double a;
  double ratio;
  paccept = VectorXd::Zero(betaInit.size());

  for (int ii=-nburn; ii<nsamples; ii++) {
    std::cout << "#### ii = " << ii << " ####" << std::endl;
    for (int jj=0; jj<betalen; jj++) {
      if (RvDoMcmc(jj)) {
        betaProp = betaOld;
        betaProp(jj) += sigs(jj)*R::norm_rand();
        // std::cout << "betaProp = " << betaProp << std::endl;
        // check if proposed beta satisfies the constraint
        bool satisfy = false;
        ELModel::evalG(betaProp);
        evalEpsilons(betaProp);
        evalWeights();
        // std::cout << "weights = " << weights.transpose() << std::endl;
        // check whether the constrains are satisfied or not
        lambdaNR(nIter, maxErr);
        if (nIter < maxIter) satisfy = true;
        // if does not satisfy, keep the old beta
        if (satisfy == false) break;
        // if does satisfy, flip a coin
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
        // TODO: suspect  the calcluation of logELProp is not correct (it is not consistent with how logELOld is obtained) check this
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

#endif
