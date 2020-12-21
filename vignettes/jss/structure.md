
# Introduction

- Motivation: What is EL?  

	- Partially specified model via moment conditions (estimand).  
	- An estimator that's "in some sense" as efficient as full parametric likelihood (refs).

- Gap: Much theory, but relatively little software (should list the ones we know of).  

	- Bulk of computations done via a convex optimization problem dual to eq (??) (describe here or relegate to methods section).  
	- Existing libraries are either written in a high-level programming language for which inner optimization is not efficient, or provide efficient implementation but for e.g., pecific regression problems (like \pkg{emplik} in \proglang{R}).
	
- Contribution: A framework for EL researchers to develop fast/efficient implementations of their own EL models and related methods.

	- This is achieved with a low-level C++ implementation of the NR method using Eigen for linear algebra.
	
	- Other EL techniques provided include: support correction, censoring, gradients.  (is this too technical for intro?)
	
- How to address theory/existing literature?

	- Should definitely cite as much existing literature as possible, theoretical or otherwise.
	
	- As for explanations, I think in many cases just the "equation" is good enough (for example, for support correction can just show how to add the observation, and explain that this is known to make NR have a unique solution for every `\bm{\theta}`.
	
