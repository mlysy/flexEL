# flexEL: Flexible Empirical Likelihood Methods for Regression

*Shimeng Huang, Martin Lysy*

---

### Description

Various tools for implementing and calibrating empirical likelihood models.  In particular, provides the loglikelihood and gradient functions for arbitrary moment constraint matrices. The inner optimization problem is efficiently computed in C++ using the [**Eigen**](http://eigen.tuxfamily.org/index.php?title=Main_Page) linear algebra library.   Also provides functions for implementing right-censored regression models, where the inner optimization is conducted via expectation-maximation.  Users may interface with the library through R or directly through C++, as the underlying C++ code is exposes as a standalone header-only library.

### Installation

Install the R package [**devtools**](https://CRAN.R-project.org/package=devtools) and run
```r
devtools::install_github("mlysy/flexEL")
```

### Usage

Please see package [vignette](https://htmlpreview.github.io/?https://github.com/mlysy/flexEL/blob/master/doc/flexEL.html).
