/**
* @file mainpage.h
* @brief Main Page documentation.
*/

/**
* @mainpage About
*
* @authors Shimeng Huang, Martin Lysy
*
* @section intro Introduction
* flexEL is a computationally efficient library in C++, which has a flexible framework for users to solve any type of regression problems with minimum programming effort. 
* 
* More details of the framework can be found <a href="https://htmlpreview.github.io/?https://github.com/mlysy/flexEL/doc/flexEL.html">here</a>.
* 
* @section struc Code Structure
* flexEL is currently capable to work with regression problems given fully observed covariates and fully observed or right-censored time-to-event resposes. The corresponding computations are handled separately by the template classes `flexEL::InnerEL` and `flexEL::InnerELC` respectively. These two classes are designed as base classes that will inherit from a regression model class which could be user defined.
* Two examples of regression problems, mean regression and quantile regression, are provided in this package, defined by the classes `MeanRegModel` and `QuantRegModel` respectively.
*/
