#ifndef InnerEL_h
#define InnerEL_h

// class for computing inner optimization of EL likelihood
// this should be inherited by specific EL problems, which then
// specify how G is computed.  Same for lambda.

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

class InnerEL {
private:
      int nObs; // number of observations (= N)
      int nDims; // length of lambda (= m)
      // constants relating to logstar calculations
      double trunc, aa, bb, cc;
      // temporary storage for Newton-Raphson
      MatrixXd GGt;
      VectorXd Q1;
      MatrixXd Q2;
      LDLT<MatrixXd> Q2ldlt;
      VectorXd Glambda;
      VectorXd rho;
      VectorXd relErr;
      // maximum relative error in lambda
      double MaxRelErr(const Ref<const VectorXd>& lambdaNew,
                       const Ref<const VectorXd>& lambdaOld);
      // columnwise outer product (see below)
      void BlockOuter(void);
public:
      // TODO: change ctor of InnerEL to assign lambdaOld ??
      InnerEL(int, int, VectorXd); // constructor
      // for accessibility
      int getnObs() const;
      int getnDims() const;
      // logstar and its derivatives
      double logstar(double x);
      double logstar1(double x);
      double logstar2(double x);
      // these should get filled by subclass
      // TODO: should create 'set' functions for them -- not directly access
      VectorXd lambdaOld;
      MatrixXd G;
      VectorXd lambdaNew;
      // Newton-Raphson algorithm
      void LambdaNR(int& nIter, double& maxErr,
                    int maxIter, double tolEps);
      // NEW: form G
      virtual void Gfun(VectorXd y, MatrixXd X, VectorXd beta) = 0; 
      // NEW: logEL declaration
      virtual double logEL(VectorXd y, MatrixXd X, VectorXd beta, int maxIter, double eps) = 0;
      // NEW: posterior sampler
      MatrixXd PostSample(int nsamples, int nburn,
                    VectorXd y, MatrixXd X, VectorXd betaInit,
                    VectorXd sigs, int maxIter, double eps);
};

// constructor
inline InnerEL::InnerEL(int N, int m, VectorXd lambda0) {
      nObs = N;
      nDims = m;
      // logstar constants
      trunc = 1.0/nObs;
      aa = -.5 * nObs*nObs;
      bb = 2.0 * nObs;
      cc = -1.5 - log(nObs);
      // Newton-Raphson initialization
      G = MatrixXd::Zero(nDims,nObs);
      GGt = MatrixXd::Zero(nDims,nObs*nDims);
      // lambdaOld = VectorXd::Zero(nDims); 
      lambdaOld = lambda0; // NEW: initialize lambdaOld here
      lambdaNew = VectorXd::Zero(nDims);
      Q1 = VectorXd::Zero(nDims);
      Q2 = MatrixXd::Zero(nDims,nDims);
      Glambda = VectorXd::Zero(nObs);
      rho = VectorXd::Zero(nObs);
      relErr = VectorXd::Zero(nDims);
      Q2ldlt.compute(MatrixXd::Identity(nDims,nDims));
}

// for accessiblity
inline int InnerEL::getnObs() const {return nObs;}
inline int InnerEL::getnDims() const {return nDims;} 

// logstar
inline double InnerEL::logstar(double x) {
      if(x >= trunc) {
            return(log(x));
      } else {
            return((aa*x + bb)*x + cc);
      }
}

// d logstar(x)/dx
inline double InnerEL::logstar1(double x) {
      if(x >= trunc) {
            return(1.0/x);
      } else {
            // return(aa*x + bb);
            return(2.0*aa*x + bb); // TODO: should be (2*aa*x + bb) ?
      }
}

// d^2 logstar(x)/dx^2
inline double InnerEL::logstar2(double x) {
      if(x >= trunc) {
            return(-1.0/(x*x));
      } else {
            // return(aa);
            return(2.0*aa); // TODO: should be 2aa ?
      }
}

// for an (m x N) matrix G = [g1 ... gN], returns the (m x mN) matrix
// GGt = [g1 g1' ... gN gN']
inline void InnerEL::BlockOuter(void) {
      // for each row of G, compute outer product and store as block
      for(int ii=0; ii<nObs; ii++) {
            GGt.block(0,ii*nDims,nDims,nDims).noalias() = G.col(ii) * G.col(ii).transpose();
      }
      return;
}

// maximum relative error in lambda
inline double InnerEL::MaxRelErr(const Ref<const VectorXd>& lambdaNew,
                                 const Ref<const VectorXd>& lambdaOld) {
      relErr = ((lambdaNew - lambdaOld).array() / (lambdaNew + lambdaOld).array()).abs();
      return(relErr.maxCoeff());
}

// Newton-Raphson algorithm
inline void InnerEL::LambdaNR(int& nIter, double& maxErr,
                              int maxIter, double tolEps) {
      int ii, jj;
      BlockOuter(); // initialize GGt
      //lambdaOld = lambdaIn; // initialize lambda
      // newton-raphson loop
      for(ii=0; ii<maxIter; ii++) {
            // Q1 and Q2
            Glambda.noalias() = lambdaOld.transpose() * G;
            Glambda = 1.0 - Glambda.array();
            Q2.fill(0.0);
            for(jj=0; jj<nObs; jj++) {
                  rho(jj) = logstar1(Glambda(jj));
                  Q2 += logstar2(Glambda(jj)) * GGt.block(0,jj*nDims,nDims,nDims);
            }
            Q1 = -G * rho;
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

// posterior sampler
inline MatrixXd InnerEL::PostSample(int nsamples, int nburn,
                               VectorXd y, MatrixXd X, VectorXd betaInit,
                               VectorXd sigs, int maxIter, double eps) {
      VectorXd betaOld = betaInit;
      VectorXd betaNew = betaOld;
      VectorXd betaProp = betaOld;
      int betalen = betaInit.size();
      MatrixXd beta_chain(betaInit.size(),nsamples);
      double logELOld = logEL(y, X, betaOld, maxIter, eps); // NEW: chache old
      
      for (int ii=-nburn; ii<nsamples; ii++) {
        for (int jj=0; jj<betalen; jj++) {
          betaProp = betaOld;
          betaProp(jj) = betaOld(jj) + sigs(jj)*R::norm_rand();
          // check if proposed beta satisfies the constraint
          bool satisfy = false;
          int nIter = 0;
          double maxErr;
          // had a BUG here?! Didn't change G!!!
          Gfun(y,X,betaProp); // NEW: change G with betaProp
          InnerEL::LambdaNR(nIter, maxErr, maxIter, eps);
          if (nIter < maxIter) satisfy = true;
          // if does not satisfy, keep the old beta
          if (satisfy == false) break;
          // if does satisfy, flip a coin
          double u = R::unif_rand();
          // use the lambda calculate just now to get the logEL for Prop
          // to avoid an extra call of lambdaNR
          VectorXd logomegahat = log(1/(1-(lambdaNew.transpose()*G).array())) - 
            log((1/(1-(lambdaNew.transpose()*G).array())).sum());
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
