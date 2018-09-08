# Summary of Code

## Code Structure

### Base Classes

The base classes correspond to the specific regression models, currently they are mean regression (`MeanRegModel`) and quantile regression (`QuantRegModel`). The main functionalities of these classes are to evaluate the `G` matrix in the empirical likelihood framework. 


### Template Classes

The template classes inherit from a base class, and the main functionalities include a Newton-Raphson solver, evaluating log empirical likeihood, and posterior samplers. Currently the template classes include one for uncensored data (`InnerEL`) and one for censored data (`InnerELC`).


## Naming Conventions

### Class Names

### Function Names

### Variables Names
