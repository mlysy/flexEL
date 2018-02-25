/*

Metropolis-within-Gibbs Adaptive MCMC

Martin Lysy, Feb 2018

// basic usage:

double mwgSd = sd_init(); // starting values of random walk jump sizes
// default constructor: all parameters are updated
mwgAdapt tuneMCMC(number_of_rvs);
// debug mode constructor: some of the parameters are fixed to investigate
// MCMC convergence issues (using fewer parameters)
// do_rv_mcmc is a bool (or int) vector of true/false (or 0/1)
mwgAdapt tuneMCMC(number_of_rvs, do_rv_mcmc);
// do_rv_mcmc was passed to C++ from R as type
// Rcpp::LogicalVector (or Rcpp::IntegerVector.  
// R always stores bool arrays internally as ints,
// because in R the array could contain NA, 
//  but in C++ bool[] can only contain T/F)
mwgAdapt tuneMCMC(number_of_rvs, LOGICAL(do_rv_mcmc));

// MCMC loop
for(int ii=0; ii<number_of_iterations; ii++) {
  for(int jj=0; jj<number_of_rvs; jj++) {
    if(do_rv_mcmc[jj]) {
      // single component mwg update,
      // takes current state of MCMC and changes only one component
      // isAccepted is boolean for whether or not proposal was accepted
      theta_curr[jj] = do_mh(theta_curr, mwgSd[jj], isAccepted[jj]);
    }
  }
  // R equivalent of Theta_out[ii,] <- theta_curr
  store_theta(Theta_out, theta_curr, ii);
  tuneMCMC.adapt(mwgSd, isAccepted);
}

*/

// -----------------------------------------------------------------------------

#ifndef mwgAdapt_h
#define mwgAdapt_h

class mwgAdapt {
 private:
  int nRV; // number of components
  bool *doAdapt; // whether or not to update each component
  double *adaptMax, *adaptRate; // adaptation max and rate
  int nIter; // number of iterations
  int *nAccept; // number of accepted proposals per component
  static const double targAcc = 0.44; // ideal acceptance rate
  double acc, lsig, delta; // temporaries for adaptation function
  void initialize(int nrv); // allocate space the same way for all constructors
 public:
  mwgAdapt(int nrv); // default adaptation parameters
  // specify max/rates for each variable
  mwgAdapt(int nrv, double *amax, double *arate);
  // same as before, but can set some of the parameters not to adapt.
  // this is mainly for debugging, e.g.,
  // when some of the parameters are fixed at known values.
  mwgAdapt(int nrv, double *amax, double *arate, bool *adapt);
  // same thing, but when adapt is integer-valued
  // this is because boolean inputs from R are actually ints, so can't be
  // converted to boolean arrays without memory allocation.
  mwgAdapt(int nrv, double *amax, double *arate, int *adapt);
  // partial parameter updates, default adaptation parameters
  mwgAdapt(int nrv, bool *adapt);
  mwgAdapt(int nrv, int *adapt);
  ~mwgAdapt;
  void adapt(double *mwgSd, bool *accept); // adapts the standard deviations
  void reset(); // resets the iteration count to zero
  // TODO: implement copy/move constuctors.
  // for now, mwgAdapt is neither copyable nor movable.
  mwgAdapt(const mwgAdapt&) = delete;
  mwgAdapt& operator=(const mwgAdapt&) = delete;
};

inline void mwgAdapt::adapt(double *mwgSd, bool *accept) {
  nIter++; // add one iteration
  for(int ii=0; ii<nRV; ii++) {
    if(doAdapt[ii]) {
      nAccept[ii] += (int) accept[ii];
      acc = (double) nAccept[ii] / (double) nIter;
      delta = pow((double) nIter, -adaptRate[ii]);
      if(delta > adaptMax[ii]) delta = adaptMax[ii];
      lsig = log(mwgSd[ii]);
      lsig += acc < targAcc ? -delta : delta;
      mwgSd[ii] = exp(lsig);
    }
  }
  return;
}

inline void mwgAdapt::reset() {
  nIter = 0;
  for(int ii=0; ii<nRV; ii++) {
    nAccept = 0;
  }
  return;
}

// allocate memory and default initialization
inline void mwgAdapt::initialize(int nrv) {
  nRV = nrv; // number of components
  reset(); // initialize nIter and nAccept to 0
  // memory allocation
  adaptMax = new double[nRV];
  adaptRate = new double[nRV];
  doAdapt = new bool[nRV];
  for(int ii=0; ii<nRV; ii++) {
    adaptMax[ii] = 0.01;
    adaptRate[ii] = 0.5;
    doAdapt[ii] = true;
  }
  return;
}

// deallocate memory
inline mwgAdapt::~mwgAdapt {
  delete [] adaptMax;
  delete [] adaptRate;
  delete [] doAdapt;
}

// --- various constructors ----------------------------------------------------

inline void mwgAdapt::mwgAdapt(int nrv) {
  initialize(nrv);
}

inline void mwgAdapt::mwgAdapt(int nrv, double *amax, double *arate) {
  initialize(nrv);
  for(int ii=0; ii<nRV; ii++) {
    adaptMax[ii] = amax[ii];
    adaptRate[ii] = arate[ii];
  }
}

inline void mwgAdapt::mwgAdapt(int nrv, double *amax, double *arate,
			       bool *adapt) {
  initialize(nrv);
  for(int ii=0; ii<nRV; ii++) {
    adaptMax[ii] = amax[ii];
    adaptRate[ii] = arate[ii];
    doAdapt[ii] = adapt[ii];
  }
}

inline void mwgAdapt::mwgAdapt(int nrv, double *amax, double *arate,
			       int *adapt) {
  initialize(nrv);
  for(int ii=0; ii<nRV; ii++) {
    adaptMax[ii] = amax[ii];
    adaptRate[ii] = arate[ii];
    doAdapt[ii] = (adapt[ii] != 0);
  }
}

inline void mwgAdapt::mwgAdapt(int nrv, bool *adapt) {
  initialize(nrv);
  for(int ii=0; ii<nRV; ii++) {
    doAdapt[ii] = adapt[ii];
  }
}

inline void mwgAdapt::mwgAdapt(int nrv, int *adapt) {
  initialize(nrv);
  for(int ii=0; ii<nRV; ii++) {
    doAdapt[ii] = (adapt[ii] != 0);
  }
}

#endif
