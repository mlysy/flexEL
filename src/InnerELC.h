#ifndef INNERELC_h
#define INNERELC_h

// class for computing inner optimization of EL likelihood with censoring

#include <vector> // for the sorting to work 
using namespace std;
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

// helper function for sorting indices according to values in a C vector
template <typename T>
class compare_acc_vec {
    const T& vec;
public:
    compare_acc_vec(const T& vec): vec(vec) { }
    bool operator () (size_t lind, size_t rind) const {
        return vec[lind] > vec[rind]; // < to sort ascendingly; > to sort descendingly
    }
};

// return indices after sorting ascendingly
// Note: have to change return type to vector<int> if want to convert to
//       Eigen Vector types, i.e., VectorXi
template <typename T>
vector<size_t> sort_inds(const vector<T> &vec) {
    vector<size_t> inds(vec.size());
    iota(inds.begin(), inds.end(), 0);
    sort(inds.begin(), inds.end(), compare_acc_vec <decltype(vec)> (vec));
    return inds;
}


// main class begins here
template <typename elModel>
class InnerELC : public elModel { 
private:
    using elModel::nObs;
    using elModel::nEqs;
    using elModel::X; // need access to X, y in evalWeights
    using elModel::y;
    // constants for logsharp calculations
    // double trunc, aa, bb, cc; 
    // temporary storage for Newton-Raphson
    MatrixXd GGt;
    VectorXd Glambda;
    ArrayXd Gl11;
    VectorXd Q1;
    MatrixXd Q2;
    LDLT<MatrixXd> Q2ldlt;
    VectorXd rho;
    VectorXd relErr;
    VectorXd epsilons; // error term used to order omegas in evalWeights
    vector<size_t> epsOrd; // order of epsilons
    VectorXd psots; // partial sum of omegatildas
    // columnwise outer product (see below)
    void BlockOuter(void);
    // maximum relative error in lambda: same for cens / non-cens
    double MaxRelErr(const Ref<const VectorXd>& lambdaNew,
                     const Ref<const VectorXd>& lambdaOld);
    // helper function for evalWeights: calculate partial sum of omegas
    // partial sum of omegas_jj s.t. eps_jj >= eps_ii
    double evalPsos(const int ii); 
public:
    // TODO: public for testing for now
    VectorXd deltas; 
    VectorXd omegas; 
    VectorXd omegasps; // partial sum of omegas 
    // TODO: W for calculation in evalWeights, 
    // it is only upper triangular, space!!!! what's better??
    MatrixXd W; 
    VectorXd weights;  
    // constructor for regression-like problems
    InnerELC(const Ref<const VectorXd>& y, const Ref<const MatrixXd>& X, 
             const Ref<const VectorXd>& deltas,
             void* params);
    // logsharp and its derivatives
    double logsharp(double x, double q);
    double logsharp1(double x, double q);
    double logsharp2(double x, double q);
    // TODO: should create 'set' functions for them 
    using elModel::G;
    VectorXd lambdaOld;
    VectorXd lambdaNew;
    // Newton-Raphson algorithm
    void LambdaNR(int& nIter, double& maxErr,
                  int maxIter, double relTol);
    void evalWeights(const Ref<const VectorXd>& beta); // calculate weights
    void EMEL(VectorXd beta, int& nIter, double& maxErr,int maxIter, double relTol);
    // log empirical likelihood calculation 
    // double logEL(const Ref<const VectorXd>& theta, int maxIter, double relTol); 
    // // posterior sampler
    // MatrixXd PostSample(int nsamples, int nburn, VectorXd betaInit,
    //                     const Ref<const VectorXd>& sigs,
    // 		      int maxIter, double relTol);
};

// constructor
// constructor for mean regression (without alpha)
template<typename elModel>
inline InnerELC<elModel>::InnerELC(const Ref<const VectorXd>& _y,
                                   const Ref<const MatrixXd>& _X,
                                   const Ref<const VectorXd>& _deltas,
                                   void* params) : elModel(_y, _X, params) {
    // std::cout << nObs << std::endl;
    // std::cout << nEqs << std::endl;
    epsilons = VectorXd::Zero(nObs);
    psots = VectorXd::Zero(nObs); 
    epsOrd = vector<size_t>(nObs); // Note: epsOrd is a C vector not Eigen VectorXd
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

// logsharp
template<typename elModel>
inline double InnerELC<elModel>::logsharp(double x, double q) {
    double xq;
    if(x >= q) {
        return(log(x));
    } else {
        xq = x/q;
        return(-0.5*xq*xq + 2.0*xq - 1.5 + log(q));
    }
}

// d logsharp(x)/dx
template<typename elModel>
inline double InnerELC<elModel>::logsharp1(double x, double q) {
    if(x >= q) {
        return(1.0/x);
    } else {
        return(-1.0/(q*q) + 2.0/q);
    }
}

// d^2 logsharp(x)/dx^2
template<typename elModel>
inline double InnerELC<elModel>::logsharp2(double x, double q) {
    if(x >= q) {
        return(-1.0/(x*x));
    } else {
        return(-1.0/(q*q));
    }
}

// for an (m x N) matrix G = [g1 ... gN], returns the (m x mN) matrix
// GGt = [g1 g1' ... gN gN']
template<typename elModel>
inline void InnerELC<elModel>::BlockOuter(void) {
  // for each row of G, compute outer product and store as block
    for(int ii=0; ii<nObs; ii++) {
        GGt.block(0,ii*nEqs,nEqs,nEqs).noalias() = G.col(ii) * G.col(ii).transpose();
    }
    return;
}

// maximum relative error in lambda
template<typename elModel>
inline double InnerELC<elModel>::MaxRelErr(const Ref<const VectorXd>& lambdaNew,
                                           const Ref<const VectorXd>& lambdaOld) {
    // std::cout << "In MaxRelErr: lambdaNew = " << lambdaNew << ", lambdaOld = " << lambdaOld << std::endl;
    relErr = ((lambdaNew - lambdaOld).array() / (lambdaNew + lambdaOld).array()).abs();
    return(relErr.maxCoeff());
}

// Newton-Raphson algorithm
template<typename elModel>
inline void InnerELC<elModel>::LambdaNR(int& nIter, double& maxErr,
                                        int maxIter, double relTol) {
    int ii, jj;
    BlockOuter(); // initialize GGt
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
        lambdaNew.noalias() = lambdaOld - Q2ldlt.solve(Q1);
        maxErr = MaxRelErr(lambdaNew, lambdaOld); // maximum relative error
        if (maxErr < relTol) {
            break;
        }
        lambdaOld = lambdaNew; // complete cycle
    }
    nIter = ii; // output lambda and also nIter and maxErr
    return;
}

// wrong way not ordering by epsilons
// template<typename elModel>
// inline void InnerELC<elModel>::evalWeights() {
//     for (int ii=0; ii<nObs; ii++) {
//         // omegas is a col vector of length nObs
//         omegasps(ii) = omegas.block(ii,0,nObs-ii,1).sum();
//     }
//     // std::cout << "omegas = \n" << omegas << std::endl;
//     // std::cout << "omegasps = \n" << omegasps << std::endl;
//     for (int jj=0; jj<nObs; jj++) {
//         // std::cout << omegas.block(jj,0,nObs-jj,1).array() / omegasps(jj) << std::endl;
//         // std::cout << "W.block(jj,0,1,nObs-jj) = \n" << W.block(jj,0,1,nObs-jj) << std::endl;
//         W.block(jj,jj,1,nObs-jj) = (omegas.block(jj,0,nObs-jj,1).array() / omegasps(jj)).transpose(); 
//     }
//     // std::cout << "W = \n" << W << std::endl;
//     // std::cout << "(1-deltas.array()).matrix().transpose() * W = " << (1-deltas.array()).matrix().transpose() * W << std::endl;
//     weights = deltas + (1-deltas.array()).matrix().transpose() * W; 
// }

template<typename elModel>
inline double InnerELC<elModel>::evalPsos(const int ii) {
    double psos = 0;
    int kk;
    for (int jj=0; jj <nObs; jj++) {
        kk = epsOrd.at(jj); // index of jj-th epsilon (ordered decreasingly)
        psos += omegas(kk);
        if (kk == ii) break;
    }
    return psos;
}

template<typename elModel>
inline void InnerELC<elModel>::evalWeights(const Ref<const VectorXd>& beta) {
    // find the indices for decreasing order of epsilons 
    epsilons.noalias() = y.transpose() - beta.transpose() * X;
    // need to work with vector version 
    vector<double> epsVec(epsilons.data(), epsilons.data()+epsilons.size());
    epsOrd = sort_inds(epsVec);
    int kk;
    for (int ii=0; ii<nObs; ii++) {
        for (int jj=nObs-1; jj>=0; jj--) { // here is backwards (smaller ones)
            kk = epsOrd.at(jj);
            if (deltas(kk) == 0) {
                psots(ii) += omegas(ii)/evalPsos(kk);
            }
            if (kk == ii) break;
        }
    }
    weights.array() = deltas.array() + psots.array();
}

// exported as omega.hat.EM
template<typename elModel>
inline void InnerELC<elModel>::EMEL(VectorXd beta, int& nIter, double& maxErr,
                                    int maxIter, double relTol) {
    // Problem: have to save this o.w. lambdaOld got modified
    VectorXd lambdaOldEM = lambdaOld; 
    int ii; 
    for(ii=0; ii<maxIter; ii++) {
        // E-step:
        evalWeights(beta);
        std::cout << "weights = " << weights << std::endl;
        // M-step:
        // std::cout << "lambdaOldEM = " << lambdaOldEM << std::endl;
        LambdaNR(nIter, maxErr, maxIter, relTol);
        std::cout << "lambdaNew = " << lambdaNew << std::endl;
        // std::cout << "G = \n" << G << std::endl;
        VectorXd lGq = ((lambdaNew.transpose() * G).array() + weights.sum()).transpose();
        std::cout << "lGq = \n" << lGq << std::endl;
        // TODO: no element-wise dividion???
        // for (int jj=0; jj<nObs; jj++){
        //     omegas(jj) = weights(jj)/lGq(jj);
        // }
        omegas.array() = weights.array() / lGq.array();
        std::cout << "In EMEL before normalize: omegas = \n" << omegas << std::endl;
        omegas = omegas.array() / (omegas.array().sum()); // normalize
        // std::cout << "In EMEL: omegas = \n" << omegas << std::endl; 
        maxErr = MaxRelErr(lambdaNew, lambdaOldEM); 
        // std::cout << "In EMEL: maxErr = " << maxErr << std::endl;
        if (maxErr < relTol) break;
        lambdaOldEM = lambdaNew;
    }
    nIter = ii; 
    return;
}

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
//       InnerELC::LambdaNR(nIter, maxErr, maxIter, relTol);
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
