#ifndef INNERELC_h
#define INNERELC_h

// class for computing inner optimization of EL likelihood with censoring

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

template <typename elModel>
class InnerELC : public elModel { 
private:
    using elModel::nObs;
    using elModel::nEqs;
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
    // columnwise outer product (see below)
    void BlockOuter(void);
    // maximum relative error in lambda: same for cens / non-cens
    double MaxRelErr(const Ref<const VectorXd>& lambdaNew,
                     const Ref<const VectorXd>& lambdaOld);
public:
    // constructor for regression-like problems
    InnerELC(const Ref<const VectorXd>& y, const Ref<const MatrixXd>& X, 
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
    void LambdaNR(VectorXd qs, int& nIter, double& maxErr,
                  int maxIter, double tolEps);
    // log empirical likelihood calculation 
    // double logEL(const Ref<const VectorXd>& theta, int maxIter, double tolEps); 
    // // posterior sampler
    // MatrixXd PostSample(int nsamples, int nburn, VectorXd betaInit,
    //                     const Ref<const VectorXd>& sigs,
    // 		      int maxIter, double tolEps);
};

// constructor
// constructor for mean regression (without alpha)
template<typename elModel>
inline InnerELC<elModel>::InnerELC(const Ref<const VectorXd>& y,
                                   const Ref<const MatrixXd>& X,
                                   void* params) : elModel(y, X, params) {
    // std::cout << nObs << std::endl;
    // std::cout << nEqs << std::endl;
    // logstar constants
    // aa = -.5 * nObs*nObs;
    // bb = 2.0 * nObs;
    // cc = -1.5 - log(nObs);
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

// logsharp
template<typename elModel>
inline double InnerELC<elModel>::logsharp(double x, double q) {
    if(x >= q) {
        return(log(x));
    } else {
        return(-1.0/(2.0*q*q)*x*x + 2/q*x - 3.0/2.0 + log(q));
    }
}

// d logsharp(x)/dx
template<typename elModel>
inline double InnerELC<elModel>::logsharp1(double x, double q) {
    if(x >= q) {
        return(1.0/x);
    } else {
        return(-1.0/(q*q) + 2/q);
    }
}

// d^2 logsharp(x)/dx^2
template<typename elModel>
inline double InnerELC<elModel>::logsharp2(double x, double q) {
    if(x >= q) {
        return(-1.0/(x*x));
    } else {
        return(-1/(q*q));
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
  relErr = ((lambdaNew - lambdaOld).array() / (lambdaNew + lambdaOld).array()).abs();
  return(relErr.maxCoeff());
}

// Newton-Raphson algorithm
template<typename elModel>
inline void InnerELC<elModel>::LambdaNR(VectorXd qs, 
                               int& nIter, double& maxErr,
                               int maxIter, double tolEps) {
    int ii, jj;
    BlockOuter(); // initialize GGt
    // newton-raphson loop
    for(ii=0; ii<maxIter; ii++) {
        // Q1 and Q2
        Glambda.noalias() = lambdaOld.transpose() * G;
        Glambda = qs.sum() + Glambda.array();
        Q2.fill(0.0); // TODO: needed? already init to 0's 
        for(jj=0; jj<nObs; jj++) {
            rho(jj) = logsharp1(Glambda(jj), qs(jj));
            Q2 += qs(jj) * logsharp2(Glambda(jj), qs(jj)) * GGt.block(0,jj*nEqs,nEqs,nEqs);
        }
        Q1 = G * (rho.array()*qs.array()).matrix();
        // update lambda
        Q2ldlt.compute(Q2);
        lambdaNew.noalias() = lambdaOld - Q2ldlt.solve(Q1);
        maxErr = MaxRelErr(lambdaNew, lambdaOld); // maximum relative error
        if (maxErr < tolEps) {
            break;
        }
        lambdaOld = lambdaNew; // complete cycle
    }
    nIter = ii; // output lambda and also nIter and maxErr
    return;
}

// // posterior sampler
// inline MatrixXd InnerELC::PostSample(int nsamples, int nburn,
//                                     VectorXd y, MatrixXd X, VectorXd betaInit,
//                                     VectorXd sigs, int maxIter, double tolEps) {
//   VectorXd betaOld = betaInit;
//   VectorXd betaNew = betaOld;
//   VectorXd betaProp = betaOld;
//   int betalen = betaInit.size();
//   MatrixXd beta_chain(betaInit.size(),nsamples);
//   double logELOld = logEL(y, X, betaOld, maxIter, tolEps); // NEW: chache old
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
//       InnerELC::LambdaNR(nIter, maxErr, maxIter, tolEps);
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
//       // double ratio = exp(logEL(y, X, betaProp, maxIter, tolEps) -
//       //                    logEL(y, X, betaOld, maxIter, tolEps));
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