#--- create package ------------------------------------------------------------

require(Rcpp) # since package has C++ code

# nice book on writing packages: http://r-pkgs.had.co.nz/

# only do this once
# Rcpp.package.skeleton(name = "bayesEL", example_code = FALSE)
# delete man/bayesEL-package.Rd
# add RcppEigen to DESCRIPTION: LinkingTo

# step 1: put C++ code into the package
# in this case, C++ code is NewtonRaphson.cpp
# any time you add C++ code, run this:
Rcpp::compileAttributes()

# to "compile" package
require(devtools)
devtools::install()

# after installing, quit + restart R

# QUESTION: seems need to create the bayesELnew-package.R file by hand before document()

# step 2: wrap C++ code in an R function and document it.
# to automate documentation with roxygen2, need to delete NAMESPACE,
# making sure package still loads C++ code properly, etc.
devtools::document() # create documentation

# step 3: test code immediately

# create a tests folder and do things informally.
# we'll do it more formally later. 


#--- next steps: ---------------------------------------------------------------

# 1.  write C++ code for EL with mean regression, quantile regression.
#     do in separate header files, with superclass InnerEL.
# 2.  EL with censoring.  Model is y = mu(x,theta) + sigma(x,theta) * eps
#     and arbitrary E[g(X, theta)] = 0

# use template class
