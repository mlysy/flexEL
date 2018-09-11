# Summary of Code

## Code Structure

### Base Classes

The base classes correspond to the specific regression models, currently include mean regression (`MeanRegModel`) and quantile regression (`QuantRegModel`). The main functionalities of these classes are to evaluate the `G` matrix in the empirical likelihood framework. 

### Template Classes

The template classes inherit from a base class, and the main functionalities include a Newton-Raphson solver, evaluating log empirical likeihood, and posterior samplers. Currently the template classes include one for uncensored data (`InnerEL`) and one for censored data (`InnerELC`).

## Naming Conventions

### C++ Class Names

Class names are all in upper Camel Case (a.k.a PascalCase). E.g., `InnerEL`, `InnerELC`, `MeanRegModel`, and `QuantRegModel`.

### C++ Class Methods and Member Variable Names

1. Functions in C++ as methods in a class are in lower Camel Case. E.g., `evalG`, `logEL`.

2. Helper functions in C++ which are in separate header files are named with underscore joining words in lower case only. E.g., `ind_smooth`, `sort_inds`. [TODO: `logstar`, `logsharp` and related methods are not following this convention right now.]

3. Member variables are in lower case or lower Camel Case if composed of multiple words. E.g., `weights`, `psots`, `nObs`, `nEqs`.

### R Function Argument / Local Variable Names

1. Names of data inputs are in English letters, if it is a matrix, its name is a single upper case letter, otherwise, it is a single lower case letter. E.g., `X`, `Z`, `y`, `G`. Exception: the vector of censoring indicators is currently called `deltas`, the residual vector is called `epsilons`. [TODO: call them `delta` and `eps`?]

2. Names of model parameters are in greek letters and lower Camel Case, if it is a vector, currently it has an "s" at the end. E.g., `beta` in `mr.evalG`, `betaInit` in `mr.post`, `alpha` in `qr.evalG` and `omegas` in `logEL`. [TODO: rename `alpha` to `tau`?] [TODO: call it `omega`?] [TODO: qr.post takes multiple quantile values are in upper Camel Case, E.g., `BetaInit`]

3. Names of dimension inputs are in lower case only. E.g., `nsamples`, `nburn`. [TODO: call it `nsamp`?]

4. Names of precision control are in lower case with underscore. E.g., `max_iter`, `rel_tol`, and `abs_tol`. [TODO: in C++ code they are called `maxIter`, `relTol`, and `absTol`, make them the same?]

5. Other local variable names are normally in lower case only, but if it is composed of multiple words, the name is in lower Camel Case. E.g., `thetaInit`, `relTol`, `mwgSd`, `rvDoMCMC`, `doAdapt`. [TODO: qr.post has argument `Sigs` which should be renamed to `mwgSd`.]

### List of Exported R Functions and Purposes

* `lambdaNR.R`: a wraper function for two Newton-Raphson solvers (uncensored, censored).
* `[mr,qr,mrls,qrls].evalG.R`: functions to evaluate G matrix corresponding to a specific regression model.
* `[qr,qrls].evalG.smooth.R`: functions to evaluate G matrix for smoothed quantile regression models. [TODO: qr.evalG.smooth.R not implemented yet]
* `logEL.R`: a wraper function for log empirical likelihood evaluations (uncensored, censored) given the probability weights `omegas`, and if censoring `epsilons` and `deltas` as well. There is no optimization in this function.
* `logEL.smooth.R`: to evaluate smoothed log empirical likelihood for censored data given `omegas`, `epsilons` and `deltas`.
* `evalEpsilons.R`: a wraper function for evaluating the residuals (location model, location-scale model). [TODO: currently it is only for location-scale model.]
* `omega.hat.R`: a wraper function to calculate the probability vector omega in a empirical likelihood setting (uncensored, censored). With censoring, an EM algorithm is used. 
* `omega.hat.EM.smooth.R`: to calculate the probability vector omega using the smoothed version of censored empirical likelihood.
* `[mr,qr,mrls,qrls,mr_cens,qr_cens].post.R`: regular MWG posterior samplers. [TODO: combine these into the post_adapt samplers with `doAdapt` argument set to false.] [TODO: these qr models are able to take multiple quantile levels at the moment, remove this feature?]
* `[mr,qr,mrls,qrls,mr_cens,qr_cens,mrls_cens,qrls_cens].post_adapt.R`: adaptive MWG posterior samplers.
* `[mr_cens,qrls_cens].neglogEL.smooth.R`: smoothed censored empirical likeihood which can be directly used for optimization routine in R, such as `nlm` and `optim`. [TODO: other regression models.]
* `hlm.R`: the HLM estimator by wang-et-al15 which works for censored data. 
* ~~`evalPsos.smooth.R`~~
* ~~`evalWeights.R`~~
* ~~`evalWeights.smooth.R`~~
* ~~`ind.smooth.R`~~
* ~~`ind1.smooth.R`~~
