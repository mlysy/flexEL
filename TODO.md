# TODO List

- [x] Use Markdown in roxygen documentation. 

- [ ] Please be consistent with roxygen punctuation!  The rule is: *everything* is sentence case, i.e., only(ish) first word is capitalized and ends with period.  This is the simplest rule to remember so it's easy to get it right :)

- [ ] Don't export internal functions.  For example, it looks like `adjG` isn't used for anything except to test C++ code.  If this is the case, the test can be conducted by replacing `adjG` by `flexEL:::.adjG`.   If you do want to export `adjG` (and it might well be useful to do so), then please clean up its documentation.  Currently the `details` section is blank :)

- [ ] `qrls_*`: What is the `nu` parameter?  I can't understand from the provided description.  In general, the models should be described in the `details` section.  Please use `\preformatted{}` (or with Markdown, triple backticks) to describe models rather than LaTeX with `\eqn{}` and `\deqn{}`.  The reason is that the latter formats differently depending on HTML or PDF output and takes a lot more work to get right.  If the model is too complicated to explain with `\preformatted{}`, then it should be documented somewhere in a vignette and the `details` section can provide a link to this.

- [ ] `qrls_evalG_smooth`: Is a separate function really necessary?  How about only `qrls_evalG` but with an `sp = 10` default (probably better to smooth :), and explaining that `sp = 0` means no smoothing.  Inside the R function you can call the corresponding C++ function `.QuantRegLSEvalGSmooth` or `.QuantRegLSEvalG`.

- [ ] `lambdaNR_cens`: Seems to me like this only relates to censoring indirectly through `weights`.  If so, I think it's better to have just `lambdaNR` with an optional `weights` argument.  Might also want to rename `.LambdaNRC` to something more informative.

- [ ] `logEL` and `logEL_smooth`: To me the distinction is not between smooth and unsmoothm, but between general EL and the specific case of EM with right censoring.  So I would drop the `delta` argument from `logEL`, and replace `logEL_smooth` by `logEL_cens`.  Where once again `sp = 0` means no smoothing.  Once again, please consider renaming C++ functions if necessary :)

- [ ] `qrls_evalG` and `qrls_evalG_smooth`: Please combine into one function.

- [ ] `omega_hat` and `omega_hat_EM_smooth`.  I would call the second `omega_hat_cens`, separating the general case from the censoring case which also includes smoothing.

	In fact, what if instead the functions `logEL` and `logEL_cens` were parametrized wrt `G`, `delta`, `eps`, etc, and a new argument `omega_out = TRUE` returns a list with elements `logEL` and `omega_hat` (as opposed to just the scalar `logEL`)?  My reasoning is that users typically won't be using `omega_hat`, so having `logEL` compute it internally makes their parameter optimization routines 1-2 lines shorter/easier to write.  On the other hand, I'm not seeing a situation where they would want only `omega_hat` and the small overhead of computing `logEL` with it makes a real difference.  Or am I missing something?  
	
	If you do go with this option, I would still make the change `omega_hat_EM_smooth -> omega_hat_cens`, but `omega_hat` and `omega_hat_cens` would be internal functions.  That means you can remove the `@export` tag, and add the `@noRd` tag to prevent `man` files from being generated.
	
- [x] Doxygen documentation.  Please use my package [**rdoxygen**](https://github.com/mlysy/rdoxygen) for this.  In particular, the Doxygen documentation would be located in its own folder instead of `src` (easier to see what is what).  

	On the `README` page for **rdoxygen** you'll see that it's possible to add the Doxygen documentation to the package as a vignette.  Upon careful reflection, I suggest we do not do this.  The reason is that Doxygen documentation considerably increases the size of the package, which CRAN doesn't like.  So right now my preference is to have the Doxydoc only in the GitHub repo, but not in the R package (use `.Rbuildignore` to exclude it).  More generally, many people (including myself) find the Doxygen HTML output fairly ugly.  I've been working on an R vignette format for this which (1) takes up way less space and (2) looks much nicer :) Unfortunately though I have zero time to spend on this project...

- [x] Move portable C++ code to `inst/include/flexEL`.  This way other R packages can easily include our C++ code (otherwise it's essentially impossible without copy-pasting).

- [ ] Make sure build passes `devtools::check()`.  This catches all sorts of issues that need to be resolved before submitting to CRAN.  It should return zero warnings and errors.  Ideal zero notes as well, but some of these are benign (e.g., "first package submission").  After this passes, please test on CRAN's Windows machines as well by running `devtools::check_win_devel()`.  You will have to set yourself as the package maintainer in the DESCRIPTION (`role = "cre"`) to send the test output to your email, not mine.

- [ ] C++ naming conventions.  I suggest you use only one naming convention for methods (snake case :)  Also, please check include guards / doxygen `@filename`s and make them consistent with the current names of the files.  Please finish the Doxygen documentation with `@param` and `@return` arguments if needed (the more thoroughly documented, the better :) Also, you can Doxygen document Rcpp functions in `src` as well!

- [ ] Vignette.  I think the organization could be improved for readers who are not familiar with EL.  How about this:

	- The EL Framework.  Rather than provide the estimator, please provide the model first.  The general model is covariate-free, i.e., you have iid observations $\bm{X}_1, \ldots, \bm{X}_n$, $\bm{X}_i \in \mathbb{R}^d$, such that $E[g_k(\bm{X}, \bm{\theta})] = 0$ for $k = 1,\ldots, m$.  The EL estimator of $\tth$ is then given by ...
	
	- Example 1: Mean Regression.  Now you can walk through the R code for this.
	
	- Other Models.  Perhaps here you can state the L/LS and MR/QR pairs.  This is also a good place to talk about smoothing.

	- Censored Regression.  The general regression model is of the form $\eps_i = g(y_i, \bm{x}_i, \tth) \iid F(\epsilon)$.  However, in survival analysis it often occurs that $y_i$ is right-censored.  Please write down the corresponding logEL optimization problem, then explain that **flexEL** solves this using an EM algorithm (cite the original Zhou paper and your thesis).  Then provide an example with e.g., `qrls`.

	- C++ API.  This is great.  What you can do now though is actually compile the C++ code calling our library!  That is, if you put
	
	```cpp
	//[Rcpp::depends("flexEL")
	```
	
	Anywhere in the C++ code, **Rcpp** will add `flexEL/include` to the source file includes.  Thus you can put
	
	```cpp
	#include <flexEL/inner_el.h>
	```
	
	And **Rcpp** will automatically know which file to use! (Prefixed by `flexEL` here because I put the source code in `inst/include/flexEL`.  The `inst` gets dropped when the package is installed.)

	- Future Work.  Excellent.
