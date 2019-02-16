/**
 * @file InnerEL.h
 * 
 * @brief Inner optimization for empirial likelihood problems.
 */

#ifndef INNEREL_h
#define INNEREL_h

#include <Rcpp.h>
#include <RcppEigen.h>
#include "BlockOuter.h" // columnwise outer product
#include "AdjG.h" // for support correction
// #include "MwgAdapt.h" // for adaptive mcmc

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

/* --------------------------------------------------------------------------- */

/**
 * @brief EL namespace
 * 
 * Wrap the exported library components into a namespace called \b el to avoid potential naming conflicts with other libraries or user-defined headers.
 */
namespace el {
  
  /**
   * @file       InnerEL.h
   *
   * @class      InnerEL
   *
   * @brief      A template class for empirical likelihood inner optimization calculation with fully observed responses.
   */
  template <typename ELModel>
  class InnerEL : public ELModel {
    
  private:
    
    // required members in ELModel
    using ELModel::nObs_; 
    using ELModel::nEqs_;
    
    // constants for logstar
    double trunc_, aa_, bb_, cc_; //constants in logstar functions
  
    // placeholders for lambdaNR
    MatrixXd G_; // G matrix for the estimating equations
    MatrixXd GUsed_;
    VectorXd lambda0_; // initial lambda in Newton-Raphson iterations
    VectorXd lambdaOld_; // old lambda in Newton-Raphson iterations
    VectorXd lambdaNew_; // new lambda in Newton-Raphson iterations
    VectorXd omegas_; //  empirical distribution
    MatrixXd GGt_; //  G times G transpose
    MatrixXd GGtUsed_;
    VectorXd Glambda_; // G times lambda
    ArrayXd Gl11_; // values of omegas before standardization
    VectorXd Q1_; // 1st derivative of the objective function
    MatrixXd Q2_; // 2nd derivative of the objective function
    LDLT<MatrixXd> Q2ldlt_; // robust Cholesky decomposition of Q2
    VectorXd rho_; // placeholder in lambdaNR
    VectorXd relErr_; //relative error
    
    // tolerance values for lambdaNR and support correction
    int maxIter_; // maximum number of iterations
    double relTol_; // relative tolerance
    bool support_; // whether do support correction or not
    double supa_; // tuning parameter for support currection

    // for support modification
    int nObs1_; // nObs1_ = nObs_+1: for initial space allocation
    int nObs2_; // nObs2_ = nObs_+support: used nObs_ in calculation
    
    // maximum relative error in lambda
    double maxRelErr();
    // double maxRelErr(const Ref<const VectorXd>& lambdaNew,
    //                  const Ref<const VectorXd>& lambdaOld);
    
    // logstar and its derivatives for the EL dual problem
    double logstar(double x); // A support-refined log function
    double logstar1(double x); // First derivative of logstar
    double logstar2(double x); // Second derivative of logstar
    
  public:
    
    // constructors
    InnerEL();
    InnerEL(int nObs, int nEqs);
    
    // set options functions
    void setOpts(const int& maxIter, const double& relTol,
                 const bool& support, const double& supa, 
                 const Ref<const VectorXd>& lambda0);
    void setOpts(const int& maxIter, const double& relTol,
                 const bool& support, 
                 const Ref<const VectorXd>& lambda0);
    void setOpts(const int& maxIter, const double& relTol,
                 const bool& support);
    void setOpts(const int& maxIter, const double& relTol,
                 const bool& support, const double& supa);
    void setOpts(const bool& support, const double& supa);
    void setOpts(const bool& support);
    // void setTol(const int& maxIter, const double& relTol);
    
    // core caculations
    void lambdaNR(int& nIter, double& maxErr); // Note: relTol and maxIter must be set before calling
    void evalOmegas(); // empirical distribution
    double logEL(); // log empirical likelihood
    
    // set and get functions 
    void setLambda(const Ref<const VectorXd>& lambda); // assigned to lambdaNew
    void setOmegas(const Ref<const VectorXd>& omegas); 
    void setG(const Ref<const MatrixXd>& G);
    VectorXd getLambda(); 
    VectorXd getOmegas();
    MatrixXd getG(); // TODO: is it better to return a reference?..
    
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
}

/* --------------------------------------------------------------------------- */

// private functions

// maximum relative error in lambda
template<typename ELModel>
inline double el::InnerEL<ELModel>::maxRelErr() {
  // TODO: added for numerical stability, what is a good tolerance to use ?
  if ((lambdaNew_ - lambdaOld_).array().abs().maxCoeff() < 1e-10) return(0);
  
  relErr_ = ((lambdaNew_ - lambdaOld_).array() / (lambdaNew_ + lambdaOld_).array()).abs();
  return(relErr_.maxCoeff());
}

// logstar
template<typename ELModel>
inline double el::InnerEL<ELModel>::logstar(double x) {
  if(x >= trunc_) {
    return(log(x));
  } else {
    return((aa_*x + bb_)*x + cc_);
  }
}

// d logstar(x)/dx
template<typename ELModel>
inline double el::InnerEL<ELModel>::logstar1(double x) {
  if(x >= trunc_) {
    return(1.0/x);
  } 
  else {
    // return(aa*x + bb);
    return(2.0*aa_*x + bb_); // TODO: should be (2*aa*x + bb) ?
  }
}

// d^2 logstar(x)/dx^2
template<typename ELModel>
inline double el::InnerEL<ELModel>::logstar2(double x) {
  if(x >= trunc_) {
    return(-1.0/(x*x));
  } else {
    // return(aa);
    return(2.0*aa_); // TODO: should be 2aa ?
  }
}

/* --------------------------------------------------------------------------- */

// public functions

// constructors
/**
 * @brief Default constructor for InnerEL.
 */
template<typename ELModel>
inline el::InnerEL<ELModel>::InnerEL(){}

/**
 * @brief Constructor for InnerEL with dimensions of G matrix as inputs for memory allocation.
 * 
 * @param nObs    Number of observations.
 * @param nEqs    Number of estimating equations.
 */
template<typename ELModel>
inline el::InnerEL<ELModel>::InnerEL(int nObs, int nEqs): ELModel(nObs, nEqs) {
  // support correction
  support_ = false;
  nObs1_ = nObs_+1;
  nObs2_ = nObs_+support_;
  // initialization of constants
  trunc_ = 1.0 / nObs_;
  aa_ = -.5 * nObs_*nObs_;
  bb_ = 2.0 * nObs_;
  cc_ = -1.5 - log(nObs_);
  // space allocation
  omegas_ = VectorXd::Zero(nObs1_).array() + 1.0/(double)nObs1_; // Initialize to 1/nObs_
  G_ = MatrixXd::Zero(nEqs_,nObs1_); // NEW: JAN 1
  GGt_ = MatrixXd::Zero(nEqs_,nObs1_*nEqs_);
  lambda0_ = VectorXd::Zero(nEqs_); // default initial value of lambda
  lambdaOld_ = VectorXd::Zero(nEqs_); // Initialize to all 0's
  lambdaNew_ = VectorXd::Zero(nEqs_);
  Q1_ = VectorXd::Zero(nEqs_);
  Q2_ = MatrixXd::Zero(nEqs_,nEqs_);
  Glambda_ = VectorXd::Zero(nObs1_);
  Gl11_ = ArrayXd::Zero(nObs1_);
  rho_ = VectorXd::Zero(nObs1_);
  relErr_ = VectorXd::Zero(nEqs_);
  Q2ldlt_.compute(MatrixXd::Identity(nEqs_,nEqs_));
}

// set options functions

/**
 * @brief Set tolerance values for NR, support correction option and tuning parameter, and initial value of lambda.
 * 
 * @param maxIter    Maximum number of iterations.
 * @param relTol     Relative tolerance.
 * @param support    Whether to have support correction.
 * @param supa       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 * @param lambda0    Initial value for lambda.
 */
template<typename ELModel>
inline void el::InnerEL<ELModel>::setOpts(const int& maxIter, const double& relTol, 
                                          const bool& support, const double& supa, 
                                          const Ref<const VectorXd>& lambda0) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  support_ = support;
  supa_ = supa;
  nObs2_ = nObs_+support_;
  lambda0_ = lambda0;
}

/**
 * @brief Set tolerance values for NR, support correction option, and initial value of lambda.
 * 
 * @param maxIter    Maximum number of iterations.
 * @param relTol     Relative tolerance.
 * @param lambda0    Initial value for lambda.
 */
template<typename ELModel>
inline void el::InnerEL<ELModel>::setOpts(const int& maxIter, const double& relTol, 
                                          const bool& support, 
                                          const Ref<const VectorXd>& lambda0) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  support_ = support;
  supa_ = std::max(1.0,0.5*log(nObs_));
  nObs2_ = nObs_+support_;
  lambda0_ = lambda0;
}

/**
 * @brief Set tolerance values for NR, support correction option. Initial value of lambda is omitted and is default to be a vector of zeros.
 * 
 * @param maxIter    Maximum number of iterations.
 * @param relTol     Relative tolerance.
 * @param support    Whether to have support correction.
 * @param supa       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 */
template<typename ELModel>
inline void el::InnerEL<ELModel>::setOpts(const int& maxIter, const double& relTol, 
                                          const bool& support, const double& supa) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  support_ = support;
  supa_ = supa;
  nObs2_ = nObs_+support_;
}

/**
 * @brief Set tolerance values for NR, support correction option. Initial value of lambda is omitted and is default to be a vector of zeros.
 * 
 * @param maxIter    Maximum number of iterations.
 * @param relTol     Relative tolerance.
 * @param support    Whether to have support correction.
 */
template<typename ELModel>
inline void el::InnerEL<ELModel>::setOpts(const int& maxIter, const double& relTol, 
                                          const bool& support) {
  maxIter_ = maxIter;
  relTol_ = relTol;
  support_ = support;
  supa_ = std::max(1.0,0.5*log(nObs_));
  nObs2_ = nObs_+support_;
}

/**
 * @brief Set support correction option and tuning parameter only. 
 * 
 * @param support    Whether to have support correction.
 * @param supa       Tuning parameter for support correction (referred as "a" in chen-et-al08).
 */
template<typename ELModel>
inline void el::InnerEL<ELModel>::setOpts(const bool& support, const double& supa) {
  support_ = support;
  supa_ = supa;
  nObs2_ = nObs_+support_;
}

/**
 * @brief Set support correction option only. 
 * 
 * @param support    Whether to have support correction.
 */
template<typename ELModel>
inline void el::InnerEL<ELModel>::setOpts(const bool& support) {
  support_ = support;
  supa_ = std::max(1.0,0.5*log(nObs_));
  nObs2_ = nObs_+support_;
}

/**
 * @brief Find the optimal lambda by a Newton-Raphson algorithm.
 * 
 * @param[out] nIter    Number of iterations to achieve convergence.
 * @param[out] maxErr   Maximum relative error among entires in lambda at the last step.
 */
template<typename ELModel>
inline void el::InnerEL<ELModel>::lambdaNR(int& nIter, double& maxErr) {
  
  lambdaOld_ = lambda0_; // set to initial value
  lambdaNew_.fill(0.0); // may not be needed here..
  
  // Note: these two cannot be preallocate untill `support` is set
  GGtUsed_ = GGt_.block(0,0,nEqs_,nObs2_*nEqs_);
  GUsed_ = G_.block(0,0,nEqs_,nObs2_);
  
  block_outer(GGtUsed_,GUsed_);
  
  // newton-raphson loop
  int ii, jj;
  for(ii=0; ii<maxIter_; ii++) {
    // Q1 and Q2
    Glambda_.noalias() = lambdaOld_.transpose() * GUsed_;
    Glambda_= 1.0 - Glambda_.array();
    Q2_.fill(0.0);
    for(jj=0; jj<nObs2_; jj++) {
      rho_(jj) = logstar1(Glambda_(jj));
      Q2_ += logstar2(Glambda_(jj)) * GGtUsed_.block(0,jj*nEqs_,nEqs_,nEqs_);
    }
    // update lambda
    Q1_ = -GUsed_ * rho_;
    Q2ldlt_.compute(Q2_);
    lambdaNew_.noalias() = lambdaOld_ - Q2ldlt_.solve(Q1_);
    maxErr = maxRelErr();
    if (maxErr < relTol_) {
      break;
    }
    lambdaOld_ = lambdaNew_; // complete cycle
  }
  nIter = ii; // output lambda and also nIter and maxErr
  return;
}

/**
 * @brief Evaluate omegas based on G and lambdaNew.
 */
template<typename ELModel>
inline void el::InnerEL<ELModel>::evalOmegas() {
  // G and lambdaNew must have been assigned
  if (lambdaNew_ != lambdaNew_) { // if lambdaNew is NaN 
    for (int ii=0; ii<nObs2_; ii++) {
      omegas_(ii) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  else {
    // if (support_) adj_G(G_,supa_);
    MatrixXd GUsed_ = G_.block(0,0,nEqs_,nObs2_); // TODO: new allocation right now..
    Glambda_.head(nObs2_).noalias() = lambdaNew_.transpose() * GUsed_;
    Gl11_.head(nObs2_) = 1.0/(1.0-Glambda_.head(nObs2_).array());
    omegas_.head(nObs2_).array() = Gl11_.head(nObs2_).array() / Gl11_.head(nObs2_).sum(); // in fact, Gl11.sum() should be equal to nObs
  }
}

/**
 * @brief Calculate logEL using omegas.
 * 
 * @return log empirical likelihood.
 */
template<typename ELModel>
inline double el::InnerEL<ELModel>::logEL() {
  // if omegas are NaN, return -Inf
  if (omegas_.head(nObs2_) != omegas_.head(nObs2_)) return -INFINITY;
  else return(omegas_.head(nObs2_).array().log().sum());
}

// setters

/**
 * @brief Set the value of lambda (e.g. to be used directly to calculate omegas).
 */
template<typename ELModel>
inline void el::InnerEL<ELModel>::setLambda(const Ref<const VectorXd>& lambda) {
  // lambdaOld_ = lambda;
  lambdaNew_ = lambda;
}

/**
 * @brief Set the value of omegas (e.g. to be used directly to calculate log EL).
 */
template<typename ELModel>
inline void el::InnerEL<ELModel>::setOmegas(const Ref<const VectorXd>& omegas) {
  omegas_.head(nObs2_) = omegas; 
}

/**
 * @brief Set the value of G (e.g. to be used directly to calculate lambda or log EL).
 */
template<typename ELModel>
inline void el::InnerEL<ELModel>::setG(const Ref<const MatrixXd>& G) {
  G_.block(0,0,nEqs_,nObs_) = G;
  if (support_) adj_G(G_,supa_);
}

// getters

/**
 * @brief Get the value of lambda.
 */
template<typename ELModel>
inline VectorXd el::InnerEL<ELModel>::getLambda() {
  return(lambdaNew_);
}

/**
 * @brief Get the value of omegas.
 */
template<typename ELModel>
inline VectorXd el::InnerEL<ELModel>::getOmegas() {
    return(omegas_.head(nObs2_));
}

/**
 * @brief Get the value of G.
 */
template<typename ELModel>
inline MatrixXd el::InnerEL<ELModel>::getG() {
  return(G_.block(0,0,nEqs_,nObs2_));
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
