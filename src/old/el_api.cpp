/// @file el_api.cpp

// Usage 1: no support adjustment
MeanRegModel mrG(X, y); // X = MatrixXd, y = VectorXd
InnerEL el(mrG.get_nObs(), mrG.get_nEq()); // memory allocation

bool support = true; // if does support adjustment
el.setOpts(support);
// el contains G_, which has nObs + 1 rows (columns?) for support adjustment

for(int ii=0; ii<nTheta; ii++) {
  mrG.set_theta(Theta.col(ii));
  // set el.G_ to mrG.evalG()
  // method 1: this has an extra copy
  mrG.evalG(tmpG); // R side: tmpG <- mrG_eval(X, y, Theta[,ii])
  el.set_G(tmpG);
  loglik[ii] = el.logEL();
  // method 2: this exposes internal storage directly (a bit unsafe)
  mrG.evalG(el.get_Gptr());
  loglik[ii] = el.logEL();
  // method 3: doesn't work well for Gadj...
  mrG.evalG(tmpG);
  loglik[ii] = el.logEL(tmpG);
}

// can we do something similar with InnerELC?
