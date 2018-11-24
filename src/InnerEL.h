/**
* @file InnerEL.h
* @brief Inner optimization for empirial likelihood problems.
*/

#ifndef INNEREL_h
#define INNEREL_h

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;
#include "BlockOuter.h" // columnwise outer product
// #include "MwgAdapt.h" // for adaptive mcmc

// [[Rcpp::depends(RcppEigen)]]

template <typename ELModel>
class InnerEL : public ELModel {
private:
  
  // required members in ELModel
  using ELModel::G; /**< matrix of estimating equations */
  using ELModel::nObs; /**< number of observations (number of columns of G) */
  using ELModel::nEqs; /**< number of estimating equations (number of rows of G) */
  
  // constants for logstar
  double trunc, aa, bb, cc; /**< constants in logstar functions */
  
  // placeholders for lambdaNR
  VectorXd lambdaOld; /**< old lambda in Newton-Raphson iterations */
  VectorXd lambdaNew; /**< new lambda in Newton-Raphson iterations */
  VectorXd omegas; /**< empirical distribution */
  MatrixXd GGt; /**< G times G transpose */
  VectorXd Glambda; /**< G times lambda */
  ArrayXd Gl11; /**< values of omegas before standardization */
  VectorXd Q1; /**< placeholder in lambdaNR */
  MatrixXd Q2; /**< placeholder in lambdaNR */
  LDLT<MatrixXd> Q2ldlt; /**< placeholder in lambdaNR */
  VectorXd rho; /**< placeholder in lambdaNR */
  VectorXd relErr; /**< relative error */

  // tolerance values for lambdaNR
  int maxIter; /**< maximum number of iterations */
  double relTol; /**< relative tolerance */

  // maximum relative error in lambda
  /**
   * @brief      Calculating the maximum relative error between \p lambdaNew and \p lambdaOld
   * @return     maximum relative error
   */
  double maxRelErr(const Ref<const VectorXd>& lambdaNew,
                   const Ref<const VectorXd>& lambdaOld);
public:
  
  // constructors
  /**
    * @brief Default constructor for InnerEL.
    */
  InnerEL();
  
  /**
   * @brief Constructor for InnerEL with dimensions of G matrix as inputs for memory allocation.
   * @param _nObs    Number of observations.
   * @param _nEqs    Number of estimating equations.
   */
  InnerEL(int _nObs, int _nEqs);
    
  // logstar and its derivatives for the EL dual problem
  /**
   * @brief A support-refined log function.
   */
  double logstar(double x);
  
  /**
   * @brief First derivative of logstar.
   */
  double logstar1(double x);
  
  /**
   * @brief Second derivative of logstar.
   */
  double logstar2(double x);

  // Newton-Raphson algorithm
  /**
   * @brief Set tolerance values for NR.
   */
  void setTol(const int& _maxIter, const double& _relTol);
  
  /**
   * @brief Find the optimal lambda by a Newton-Raphson algorithm.
   * @param[out] nIter    Number of iterations to achieve convergence.
   * @param[out] maxErr   Maximum relative error among entires in lambda at the last step.
   */
  void lambdaNR(int& nIter, double& maxErr); // Note: relTol and maxIter must be set before calling
  
  // EL evaluation
  /**
   * @brief Evaluate omegas based on G and lambdaNew.
   */
  void evalOmegas();
  
  /**
   * @brief Calculate logEL using omegas.
   */
  double logEL();
  
  // set and get functions 
  void setLambda(const Ref<const VectorXd>& _lambda); // assigned to lambdaNew
  void setOmegas(const Ref<const VectorXd>& _omegas); 
  void setG(const Ref<const MatrixXd>& _G);
  VectorXd getLambda(); 
  VectorXd getOmegas();
  MatrixXd getG(); 
    
    // nBet, nGam and nQts are FOR THE MCMC SAMPELERS
    // using ELModel::nBet;
    // using ELModel::nGam;
    // using ELModel::nQts;
    // // posterior samplers:
    // void mwgStep(VectorXd &thetaCur, const int &idx, const double &mwgsd,
    //              bool &accept, double &logELCur);
    // MatrixXd postSample(int nsamples, int nburn, MatrixXd ThetaInit,
    //                     const Ref<const MatrixXd>& MwgSds, 
    //                     MatrixXd &RvDoMcmc, MatrixXd &Paccept);
    // MatrixXd postSampleAdapt(int nsamples, int nburn, VectorXd thetaInit,
    //                          double *mwgSd, bool *rvDoMcmc, VectorXd &paccept);
};

// default ctor 
template<typename ELModel>
inline InnerEL<ELModel>::InnerEL(){}

// ctor with dimensions as input
template<typename ELModel>
inline InnerEL<ELModel>::InnerEL(int _nObs, int _nEqs): ELModel(_nObs, _nEqs){
  omegas = VectorXd::Zero(nObs).array() + 1.0/(double)nObs; // Initialize to 1/nObs
  trunc = 1.0 / nObs;
  aa = -.5 * nObs*nObs;
  bb = 2.0 * nObs;
  cc = -1.5 - log(nObs);
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

// logstar
template<typename ELModel>
inline double InnerEL<ELModel>::logstar(double x) {
  if(x >= trunc) {
    return(log(x));
  } else {
    return((aa*x + bb)*x + cc);
  }
}

// d logstar(x)/dx
template<typename ELModel>
inline double InnerEL<ELModel>::logstar1(double x) {
  if(x >= trunc) {
    return(1.0/x);
  } 
  else {
    // return(aa*x + bb);
    return(2.0*aa*x + bb); // TODO: should be (2*aa*x + bb) ?
  }
}

// d^2 logstar(x)/dx^2
template<typename ELModel>
inline double InnerEL<ELModel>::logstar2(double x) {
  if(x >= trunc) {
    return(-1.0/(x*x));
  } else {
    // return(aa);
    return(2.0*aa); // TODO: should be 2aa ?
  }
}

// maximum relative error in lambda
template<typename ELModel>
inline double InnerEL<ELModel>::maxRelErr(const Ref<const VectorXd>& lambdaNew,
                                 const Ref<const VectorXd>& lambdaOld) {
  // TODO: added for numerical stability, what is a good tolerance to use ?
  if ((lambdaNew - lambdaOld).array().abs().maxCoeff() < 1e-10) return(0);
  
  relErr = ((lambdaNew - lambdaOld).array() / (lambdaNew + lambdaOld).array()).abs();
  return(relErr.maxCoeff());
}

// set tolerance values for lambdaNR
template<typename ELModel>
inline void InnerEL<ELModel>::setTol(const int& _maxIter, const double& _relTol) {
  maxIter = _maxIter;
  relTol = _relTol;
}

// Note: maxIter and relTol must have been assigned 
// Newton-Raphson algorithm
template<typename ELModel>
inline void InnerEL<ELModel>::lambdaNR(int& nIter, double& maxErr) {
  lambdaOld.fill(0.0);
  lambdaNew.fill(0.0);
  int ii, jj;
  // blockOuter(); // initialize GGt
  block_outer(GGt,G);
  // newton-raphson loop
  for(ii=0; ii<maxIter; ii++) {
    // Q1 and Q2
    Glambda.noalias() = lambdaOld.transpose() * G;
    Glambda = 1.0 - Glambda.array();
    Q2.fill(0.0);
    for(jj=0; jj<nObs; jj++) {
        rho(jj) = logstar1(Glambda(jj));
        Q2 += logstar2(Glambda(jj)) * GGt.block(0,jj*nEqs,nEqs,nEqs);
    }
    Q1 = -G * rho;
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

// Note: G and lambdaNew must have been assigned before calling
// i.e., in InnerElExports.cpp, assign lambdaNew first 
template<typename ELModel>
inline void InnerEL<ELModel>::evalOmegas() {
  // G and lambdaNew must have been assigned
  if (lambdaNew != lambdaNew) { // if lambdaNew is NaN 
    for (int ii=0; ii<nObs; ii++) {
      omegas(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  else {
    Glambda.noalias() = lambdaNew.transpose() * G;
    Gl11 = 1.0/(1.0-Glambda.array());
    omegas.array() = Gl11.array() / Gl11.sum(); // in fact, Gl11.sum() should be equal to nObs
  }
}

// Note: omegas must have been assigned before calling 
// i.e., in InnerElExports.cpp, evaluate & assign omegas first 
template<typename ELModel>
inline double InnerEL<ELModel>::logEL() {
  // if omegas are NaN, return -Inf
  if (omegas != omegas) return -INFINITY;
  else return(omegas.array().log().sum());
}

// setter for lamba
template<typename ELModel>
inline void InnerEL<ELModel>::setLambda(const Ref<const VectorXd>& _lambda) {
  lambdaNew = _lambda;
}

// getter for lamba
template<typename ELModel>
inline VectorXd InnerEL<ELModel>::getLambda() {
  return(lambdaNew);
}

// setter for omegas
template<typename ELModel>
inline void InnerEL<ELModel>::setOmegas(const Ref<const VectorXd>& _omegas) {
    omegas = _omegas; 
}

// getter for omegas
template<typename ELModel>
inline VectorXd InnerEL<ELModel>::getOmegas() {
    return(omegas);
}

// set function for G matrix
template<typename ELModel>
inline void InnerEL<ELModel>::setG(const Ref<const MatrixXd>& _G) {
  G = _G; 
}

// get function for G matrix
template<typename ELModel>
inline MatrixXd InnerEL<ELModel>::getG() {
  return(G);
}

/* TAKE OUT ALL MCMC SAMPLERS FOR NOW:

// This works for location models but also for multiple quantile levels
template<typename ELModel>
inline MatrixXd InnerEL<ELModel>::postSample(int nsamples, int nburn,
                MatrixXd ThetaInit, const Ref<const MatrixXd>& MwgSds,
                MatrixXd &RvDoMcmc, MatrixXd &Paccept) {
  int nThe = ThetaInit.rows(); // dimension of Theta
  int numThe = ThetaInit.cols(); // numer of Thetas
  MatrixXd Theta_chain(nThe*numThe,nsamples);
  MatrixXd ThetaOld = ThetaInit;
  MatrixXd ThetaNew = ThetaOld;
  MatrixXd ThetaProp = ThetaOld;
  ELModel::evalG(ThetaOld);
  int nIter;
  double maxErr;
  lambdaNR(nIter, maxErr);
  // TODO: what if not converged ?
  if (nIter == maxIter && maxErr > relTol) {
    std::cout << "ThetaInit not valid." << std::endl;
    // return NULL;
  }
  evalOmegas();
  double logELOld = logEL(); 
  double logELProp;
  bool satisfy;
  double u;
  double a;
  double ratio;
  Paccept = MatrixXd::Zero(ThetaInit.rows(),ThetaInit.cols()); // TODO: need to initialize to 0 ?
  
  bool go_next;
  for (int ii=-nburn; ii<nsamples; ii++) {
    go_next = false;
    for (int kk=0; kk<numThe; kk++) {  
      if (go_next == true) break;
      for (int jj=0; jj<nThe; jj++) {
        if (RvDoMcmc(jj,kk)) {
          ThetaProp = ThetaOld;
          ThetaProp(jj,kk) += MwgSds(jj,kk)*R::norm_rand();
          // check if proposed Theta satisfies the constraint
          satisfy = false;
          ELModel::evalG(ThetaProp); // NEW: change G with ThetaProp
          lambdaNR(nIter, maxErr);
          if (nIter < maxIter) satisfy = true;
          // if does not satisfy, keep the old Theta
          if (satisfy == false) {
            go_next = true; // break out two loops
            break;
          }
          // if does satisfy
          u = R::unif_rand();
          // use the lambda calculate just now to get the logEL for Prop
          // to avoid an extra call of lambdaNR
          VectorXd logomegahat = 1/(1-(lambdaNew.transpose()*G).array());
          logomegahat = log(logomegahat.array()) - log(logomegahat.sum());
          // VectorXd logomegahat = log(1/(1-(lambdaNew.transpose()*ELModel::G).array())) -
          //   log((1/(1-(lambdaNew.transpose()*ELModel::G).array())).sum());
          logELProp = logomegahat.sum();
          ratio = exp(logELProp-logELOld);
          a = std::min(1.0,ratio);
          if (u < a) { // accepted
            Paccept(jj,kk) += 1; 
            ThetaNew = ThetaProp;
            ThetaOld = ThetaNew;
            logELOld = logELProp; // NEW: store the new one
          }
        }
      }
    }
    if (ii >= 0) {
      // Theta_chain.col(ii) = ThetaNew;
      Theta_chain.col(ii) = Map<VectorXd>(ThetaNew.data(), ThetaNew.size());
    }
  }
  Paccept /= (nsamples+nburn); 
  return(Theta_chain);
}

// mwgStep updates the idx entry of thetaCur
template<typename ELModel>
inline void InnerEL<ELModel>::mwgStep(VectorXd &thetaCur,
                                      const int &idx,
                                      const double &mwgsd,
                                      bool &accept, 
                                      double &logELCur) {
  int nThe = thetaCur.size();
  accept = false;
  VectorXd thetaProp = thetaCur;
  thetaProp(idx) += mwgsd*R::norm_rand();
  // sig2 has to be positive
  if (idx == nBet+nGam && thetaProp(idx) < 0) return;
  
  if (nThe == nBet) {
    // location mean regression and quantile regression (single quantile)
    ELModel::evalG(thetaProp);
  }
  else if (nThe == nBet + nGam + 1){
    // location-scale mean regression
    ELModel::evalG(thetaProp.head(nBet), 
                   thetaProp.segment(nBet,nGam), 
                   thetaProp.tail(1)(0),
                   VectorXd::Zero(0));
  }
  else {
    // location-scale quantile regression (single or multiple quantile)
    // using ELModel::nQts;
    ELModel::evalG(thetaProp.head(nBet), 
                   thetaProp.segment(nBet,nGam), 
                   thetaProp.segment(nBet+nGam,1)(0),
                   thetaProp.tail(nQts));
  }
  int nIter;
  double maxErr;
  lambdaNR(nIter, maxErr);
  bool satisfy = false;
  if (nIter < maxIter || maxErr <= relTol) satisfy = true;
  // if does not satisfy, keep the old theta
  if (satisfy == false) return;
  // if does satisfy, might accept new theta
  double u = R::unif_rand();
  VectorXd logomegahat = log(1/(1-(lambdaNew.transpose()*ELModel::G).array())) -
    log((1/(1-(lambdaNew.transpose()*ELModel::G).array())).sum());
  double logELProp = logomegahat.sum();
  double ratio = exp(logELProp-logELCur);
  double a = std::min(1.0,ratio);
  if (u < a) { // accepted
    accept = true;
    thetaCur = thetaProp;
    logELCur = logELProp;
  }
}

template<typename ELModel>
inline MatrixXd InnerEL<ELModel>::postSampleAdapt(int nsamples, int nburn,
                                                  VectorXd thetaInit,
                                                  double *mwgSd, bool *rvDoMcmc,
                                                  VectorXd &paccept) {
  
  int nThe = thetaInit.size();
  MwgAdapt tuneMCMC(nThe, rvDoMcmc);
  bool *isAccepted = new bool[nThe];
  for (int ii=0; ii<nThe; ii++) {
    isAccepted[ii] = false;
  }
  MatrixXd theta_chain(nThe,nsamples);
  paccept = VectorXd::Zero(nThe);
  VectorXd thetaCur = thetaInit;
  if (nThe == nBet) {
    ELModel::evalG(thetaCur);
  }
  else if (nThe == nBet + nGam + 1){
    ELModel::evalG(thetaCur.head(nBet), 
                   thetaCur.segment(nBet,nGam), 
                   thetaCur.tail(1)(0),
                   VectorXd::Zero(0));
  }
  else {
    // using ELModel::nQts;
    ELModel::evalG(thetaCur.head(nBet), 
                   thetaCur.segment(nBet,nGam), 
                   thetaCur.segment(nBet+nGam,1)(0),
                   thetaCur.tail(nQts));
    // std::cout << "G = \n" << G << std::endl;
  }
  int nIter;
  double maxErr;
  lambdaNR(nIter, maxErr);
  // TODO: throw an error ??
  if (nIter == maxIter && maxErr > relTol) {
    std::cout << "thetaInit not valid." << std::endl;
  }
  evalOmegas();
  double logELCur = logEL();
  // MCMC loop
  for(int ii=-nburn; ii<nsamples; ii++) {
    for(int jj=0; jj<nThe; jj++) {
      if(rvDoMcmc[jj]) {
        // std::cout << "isAccepted[" << jj << "] = " << isAccepted[jj] << std::endl;
        // modifies thetaCur's jj-th entry
        mwgStep(thetaCur,jj,mwgSd[jj],isAccepted[jj],logELCur);
        if (isAccepted[jj]) paccept(jj) += 1; // add 1 to paccept if accepted
      }
    }
    if (ii >= 0) {
      theta_chain.col(ii) = thetaCur;
    }
    tuneMCMC.adapt(mwgSd, isAccepted);
  }
  paccept /= (nsamples+nburn);
  delete[] isAccepted; // deallocate memory
  return(theta_chain);
}

*/

#endif
