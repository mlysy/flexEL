# Summary of Code

## Code Structure

### Base Classes

The base classes correspond to the specific regression models, currently they are mean regression (`MeanRegModel`) and quantile regression (`QuantRegModel`). The main functionalities of these classes are to evaluate the `G` matrix in the empirical likelihood framework. 


### Template Classes

The template classes inherit from a base class, and the main functionalities include a Newton-Raphson solver, evaluating log empirical likeihood, and posterior samplers. Currently the template classes include one for uncensored data (`InnerEL`) and one for censored data (`InnerELC`).

## Naming Conventions

### Class Names

Class names are all in upper camel case (a.k.a PascalCase). E.g., `InnerEL`.

### Function Names

1. Functions as methods in a class are all in lower camel case. E.g., `evalG`.

2. Helper functions in separate header files are named with underscore joining words in lower case only. E.g., `ind_smooth`.

### Function Argument Names

1. Names of data inputs are in English letters, if it is a matrix, its name is a single upper case letter, otherwise, it is a single lower case letter. E.g., `X`, `Z`, `y`. Exception: the vector of censoring indicators is currently called `deltas`.

2. Names of model parameters are in greek letters, if it is a vector, currently it has an "s" at the end. E.g., `omegas`. [TODO: plan to change this to without "s"]

3. Other variable names are normally in lower case only, but if it is composed of multiple words, the name is in lower camel case. E.g., `thetaInit`, `relTol`.

### List of Exported R Functions and Purpose


