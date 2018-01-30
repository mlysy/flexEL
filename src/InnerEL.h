#ifndef INNEREL_h
#define INNEREL_h

// class for computing inner optimization of EL likelihood

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

template <typename elModel>
class InnerEL : public elModel { 
private:
  // constants for logstar calculations
  double trunc, aa, bb, cc; 
  // temporary storage for Newton-Raphson
  // MatrixXd GGt;
  VectorXd Q1;
  MatrixXd Q2;
  LDLT<MatrixXd> Q2ldlt;
  VectorXd Glambda;
  VectorXd rho;
  VectorXd relErr;
  // columnwise outer product (see below)
  void BlockOuter(void);
  // maximum relative error in lambda: same for cens / non-cens
  double MaxRelErr(const Ref<const VectorXd>& lambdaNew,
                   const Ref<const VectorXd>& lambdaOld);
public:
  InnerEL(VectorXd, MatrixXd, int, int, VectorXd); // constructor for mean regression 
  InnerEL(VectorXd, MatrixXd, double, int, int, VectorXd); // constructor for quantile regression 
  // logstar and its derivatives
  double logstar(double x);
  double logstar1(double x);
  double logstar2(double x);
  // TODO: should create 'set' functions for them 
  VectorXd lambdaOld;
  VectorXd lambdaNew;
  // MatrixXd G;
  // // TODO: for accessibility, remove if not used 
  // int getnObs() const;
  // int getnEqs() const;
  // Newton-Raphson algorithm
  void LambdaNR(int& nIter, double& maxErr,
                int maxIter, double tolEps);
  // log empirical likelihood calculation 
  double logEL(VectorXd beta, int maxIter, double eps); 
  // posterior sampler
  MatrixXd PostSample(int nsamples, int nburn, VectorXd betaInit,
                      VectorXd sigs, int maxIter, double eps);
};

// constructor for mean regression (without alpha)
template<typename elModel>
inline InnerEL<elModel>::InnerEL(VectorXd y, MatrixXd X, int nObs, int nEqs, VectorXd lambda0) : elModel(y, X, nObs,nEqs) {
  // logstar constants
  trunc = 1.0/nObs;
  aa = -.5 * nObs*nObs;
  bb = 2.0 * nObs;
  cc = -1.5 - log(nObs);
  // Newton-Raphson initialization
  // G = MatrixXd::Zero(nEqs,nObs);
  // GGt = MatrixXd::Zero(nEqs,nObs*nEqs);
  // lambdaOld = VectorXd::Zero(nEqs); 
  lambdaOld.noalias() = lambda0; // TODO: is this a copy or reference? 
  lambdaNew = VectorXd::Zero(nEqs);
  Q1 = VectorXd::Zero(nEqs);
  Q2 = MatrixXd::Zero(nEqs,nEqs);
  Glambda = VectorXd::Zero(nObs);
  rho = VectorXd::Zero(nObs);
  relErr = VectorXd::Zero(nEqs);
  Q2ldlt.compute(MatrixXd::Identity(nEqs,nEqs));
}

// constructor for quantile regression (with alpha)
// TODO: too much copying???
template<typename elModel>
inline InnerEL<elModel>::InnerEL(VectorXd y, MatrixXd X, double alpha, 
                                 int nObs, int nEqs, VectorXd lambda0) : elModel(y, X, alpha, nObs,nEqs) {
  // logstar constants
  trunc = 1.0/nObs;
  aa = -.5 * nObs*nObs;
  bb = 2.0 * nObs;
  cc = -1.5 - log(nObs);
  // Newton-Raphson initialization
  // G = MatrixXd::Zero(nEqs,nObs);
  // GGt = MatrixXd::Zero(nEqs,nObs*nEqs);
  // lambdaOld = VectorXd::Zero(nEqs); 
  lambdaOld.noalias() = lambda0; // TODO: is this a copy or reference? 
  lambdaNew = VectorXd::Zero(nEqs);
  Q1 = VectorXd::Zero(nEqs);
  Q2 = MatrixXd::Zero(nEqs,nEqs);
  Glambda = VectorXd::Zero(nObs);
  rho = VectorXd::Zero(nObs);
  relErr = VectorXd::Zero(nEqs);
  Q2ldlt.compute(MatrixXd::Identity(nEqs,nEqs));
}

// // for accessiblity
// template<typename elModel>
// inline int InnerEL<elModel>::getnObs() const {return nObs;}
// template<typename elModel>
// inline int InnerEL<elModel>::getnEqs() const {return nEqs;} 

// logstar
template<typename elModel>
inline double InnerEL<elModel>::logstar(double x) {
  if(x >= trunc) {
    return(log(x));
  } else {
    return((aa*x + bb)*x + cc);
  }
}

// d logstar(x)/dx
template<typename elModel>
inline double InnerEL<elModel>::logstar1(double x) {
  if(x >= trunc) {
    return(1.0/x);
  } else {
    // return(aa*x + bb);
    return(2.0*aa*x + bb); // TODO: should be (2*aa*x + bb) ?
  }
}

// d^2 logstar(x)/dx^2
template<typename elModel>
inline double InnerEL<elModel>::logstar2(double x) {
  if(x >= trunc) {
    return(-1.0/(x*x));
  } else {
    // return(aa);
    return(2.0*aa); // TODO: should be 2aa ?
  }
}

// for an (m x N) matrix G = [g1 ... gN], returns the (m x mN) matrix
// GGt = [g1 g1' ... gN gN']
template<typename elModel>
inline void InnerEL<elModel>::BlockOuter(void) {
  // for each row of G, compute outer product and store as block
  for(int ii=0; ii<this->nObs; ii++) {
    elModel::GGt.block(0,ii*this->nEqs,this->nEqs,this->nEqs).noalias() = elModel::G.col(ii) * elModel::G.col(ii).transpose();
  }
  return;
}

// maximum relative error in lambda
template<typename elModel>
inline double InnerEL<elModel>::MaxRelErr(const Ref<const VectorXd>& lambdaNew,
                                 const Ref<const VectorXd>& lambdaOld) {
  relErr = ((lambdaNew - lambdaOld).array() / (lambdaNew + lambdaOld).array()).abs();
  return(relErr.maxCoeff());
}

// Newton-Raphson algorithm
template<typename elModel>
inline void InnerEL<elModel>::LambdaNR(int& nIter, double& maxErr,
                              int maxIter, double tolEps) {
  int ii, jj;
  BlockOuter(); // initialize GGt
  //lambdaOld = lambdaIn; // initialize lambda
  // newton-raphson loop
  for(ii=0; ii<maxIter; ii++) {
    // Q1 and Q2
    Glambda.noalias() = lambdaOld.transpose() * elModel::G;
    Glambda = 1.0 - Glambda.array();
    Q2.fill(0.0);
    for(jj=0; jj<this->nObs; jj++) {
      rho(jj) = logstar1(Glambda(jj));
      Q2 += logstar2(Glambda(jj)) * elModel::GGt.block(0,jj*this->nEqs,this->nEqs,this->nEqs);
    }
    Q1 = -elModel::G * rho;
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

// log empirical likelihood for linear regression Y = X'beta + eps
// beta is ncol(X) x 1 = p x 1 vector
template<typename elModel>
inline double InnerEL<elModel>::logEL(VectorXd beta, int maxIter, double eps) {
  int nIter;
  double maxErr;
  elModel::evalG(beta);
  LambdaNR(nIter, maxErr, maxIter, eps);
  VectorXd logomegahat = log(1/(1-(lambdaNew.transpose()*elModel::G).array())) - 
    log((1/(1-(lambdaNew.transpose()*elModel::G).array())).sum());
  return(logomegahat.sum());
}

// posterior sampler
template<typename elModel>
inline MatrixXd InnerEL<elModel>::PostSample(int nsamples, int nburn, VectorXd betaInit,
                                    VectorXd sigs, int maxIter, double eps) {
  VectorXd betaOld = betaInit;
  VectorXd betaNew = betaOld;
  VectorXd betaProp = betaOld;
  int betalen = betaInit.size();
  MatrixXd beta_chain(betaInit.size(),nsamples);
  double logELOld = logEL(betaOld, maxIter, eps); // NEW: chache old
  
  for (int ii=-nburn; ii<nsamples; ii++) {
    for (int jj=0; jj<betalen; jj++) {
      betaProp = betaOld;
      betaProp(jj) = betaOld(jj) + sigs(jj)*R::norm_rand();
      // check if proposed beta satisfies the constraint
      bool satisfy = false;
      int nIter = 0;
      double maxErr;
      // had a BUG here?! Didn't change G!!!
      elModel::evalG(betaProp); // NEW: change G with betaProp
      LambdaNR(nIter, maxErr, maxIter, eps);
      if (nIter < maxIter) satisfy = true;
      // if does not satisfy, keep the old beta
      if (satisfy == false) break;
      // if does satisfy, flip a coin
      double u = R::unif_rand();
      // use the lambda calculate just now to get the logEL for Prop
      // to avoid an extra call of lambdaNR
      VectorXd logomegahat = log(1/(1-(lambdaNew.transpose()*elModel::G).array())) - 
        log((1/(1-(lambdaNew.transpose()*elModel::G).array())).sum());
      double logELProp = logomegahat.sum();
      double ratio = exp(logELProp-logELOld);
      // double ratio = exp(logEL(y, X, betaProp, maxIter, eps) -
      //                    logEL(y, X, betaOld, maxIter, eps));
      double a = std::min(1.0,ratio);
      if (u < a) { // accepted
        betaNew = betaProp;
        betaOld = betaNew;
        logELOld = logELProp; // NEW: store the new one
      }
    }
    if (ii >= 0) {
      beta_chain.col(ii) = betaNew;
    }
  }
  return(beta_chain);
}

#endif