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
    void lambdaNR(int& nIter, double& maxErr, 
                  int maxIter, double relTol);
    // eval functions 
    void evalWeights(); // calculate weights according to epsilons 
    // void evalOmegas(int& nIter, double& maxErr, int maxIter, double relTol);
    void evalOmegas(int maxIter, double relTol);
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
    // // posterior sampler
    // MatrixXd PostSample(int nsamples, int nburn, VectorXd betaInit,
    //                     const Ref<const VectorXd>& sigs,
    // 		      int maxIter, double relTol);
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
    // std::cout << "In maxRelErr: lambdaNew = " << lambdaNew << ", lambdaOld = " << lambdaOld << std::endl;
    relErr = ((lambdaNew - lambdaOld).array() / (lambdaNew + lambdaOld).array()).abs();
    return(relErr.maxCoeff());
}

// Newton-Raphson algorithm
template<typename ELModel>
inline void InnerELC<ELModel>::lambdaNR(int& nIter, double& maxErr,
                                        int maxIter, double relTol) {
    int ii, jj;
    blockOuter(); // initialize GGt according to epsilons order (epsOrd)
    // newton-raphson loop
    for(ii=0; ii<maxIter; ii++) {
        // Q1 and Q2
        Glambda.noalias() = lambdaOld.transpose() * G;
        Glambda = weights.sum() + Glambda.array();
        Q2.fill(0.0); // TODO: needed? already init to 0's 
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

// Note: epsilons must have been assigned
template<typename ELModel>
inline void InnerELC<ELModel>::evalWeights() {
    // find the indices for increasing order of epsilons 
    epsOrd = sort_inds(epsilons);  // TODO: might remove it from here
    psots.noalias() = VectorXd::Zero(nObs); // clear previous values
    int kk;
    for (int ii=0; ii<nObs; ii++) {
        for (int jj=0; jj<nObs; jj++) {
            kk = epsOrd(jj);
            if (deltas(kk) == 0) {
                psots(ii) += omegas(ii)/evalPsos(kk);
            }
            if (kk == ii) break;
        }
    }
    // assigned weights are still in the order of original data
    weights.array() = deltas.array() + psots.array();
    // std::cout << "evalWeights: deltas = " << deltas.transpose() << std::endl;
    // std::cout << "evalWeights: epsilons = " << epsilons.transpose() << std::endl;
    // std::cout << "evalWeights: psots = " << psots.transpose() << std::endl;
}

// exported as `omega.hat.EM`
// Note: epsilons must have been assigned
template<typename ELModel>
inline void InnerELC<ELModel>::evalOmegas(int maxIter, double relTol) {
  // Problem: have to save this o.w. lambdaOld got modified
  VectorXd lambdaOldEM = lambdaOld; 
    int ii; 
    int nIter; 
    double maxErr;
    for(ii=0; ii<maxIter; ii++) {
        // E-step:
        evalWeights(); // assigns weights according to epsilons
        // std::cout << "weights = " << weights.transpose() << std::endl;
        // M-step:
        // std::cout << "lambdaOldEM = " << lambdaOldEM.transpose() << std::endl;
        lambdaNR(nIter, maxErr, maxIter, relTol); 
        // Check convergence of NR here
        if (nIter == maxIter & maxErr > relTol) {
          std::cout << "lambdaNRC did not converge in EM" << std::endl;
        }
        VectorXd lGq = ((lambdaNew.transpose() * G).array() + weights.sum()).transpose();
        // std::cout << "lGq = \n" << lGq.transpose() << std::endl;
        // std::cout << "lambdaNew = " << lambdaNew.transpose() << std::endl;
        // std::cout << "lambdaOld = " << lambdaOld.transpose() << std::endl;
        // std::cout << "lambdaOldEM = " << lambdaOldEM.transpose() << std::endl;
        // if (ii == 1) break;
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

// // posterior sampler
// inline MatrixXd InnerELC::PostSample(int nsamples, int nburn,
//                                     VectorXd y, MatrixXd X, VectorXd betaInit,
//                                     VectorXd sigs, int maxIter, double relTol) {
//   VectorXd betaOld = betaInit;
//   VectorXd betaNew = betaOld;
//   VectorXd betaProp = betaOld;
//   int betalen = betaInit.size();
//   MatrixXd beta_chain(betaInit.size(),nsamples);
//   double logELOld = logEL(y, X, betaOld, maxIter, relTol); // NEW: chache old
//   
//   for (int ii=-nburn; ii<nsamples; ii++) {
//     for (int jj=0; jj<betalen; jj++) {
//       betaProp = betaOld;
//       betaProp(jj) = betaOld(jj) + sigs(jj)*R::norm_rand();
//       // check if proposed beta satisfies the constraint
//       bool satisfy = false;
//       int nIter = 0;
//       double maxErr;
//       // had a BUG here?! Didn't change G!!!
//       Gfun(y,X,betaProp); // NEW: change G with betaProp
//       InnerELC::lambdaNR(nIter, maxErr, maxIter, relTol);
//       if (nIter < maxIter) satisfy = true;
//       // if does not satisfy, keep the old beta
//       if (satisfy == false) break;
//       // if does satisfy, flip a coin
//       double u = R::unif_rand();
//       // use the lambda calculate just now to get the logEL for Prop
//       // to avoid an extra call of lambdaNR
//       VectorXd logomegahat = log(1/(1-(lambdaNew.transpose()*G).array())) - 
//         log((1/(1-(lambdaNew.transpose()*G).array())).sum());
//       double logELProp = logomegahat.sum();
//       double ratio = exp(logELProp-logELOld);
//       // double ratio = exp(logEL(y, X, betaProp, maxIter, relTol) -
//       //                    logEL(y, X, betaOld, maxIter, relTol));
//       double a = std::min(1.0,ratio);
//       if (u < a) { // accepted
//         betaNew = betaProp;
//         betaOld = betaNew;
//         logELOld = logELProp; // NEW: store the new one
//       }
//     }
//     if (ii >= 0) {
//       beta_chain.col(ii) = betaNew;
//     }
//   }
//   return(beta_chain);
// }

#endif
