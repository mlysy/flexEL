---
title:
  formatted: "\\pkg{flexEL}: A Fast and Flexible Framework for Empirical Likelihood Modeling"
  plain: "flexEL: A Fast and Flexible Framework for Empirical Likelihood Modeling"
author: 
  - name: Shimeng Huang
    affiliation: University of Waterloo
  - name: Martin Lysy
    affiliation: University of Waterloo
    address:
      - 200 University Avenue West
      - Ontario, Canada
    email: \email{mlysy@uwaterloo.ca}
abstract: >
  Here is the abstract.
keywords:
  # at least one keyword must be supplied
  formatted: [keywords, not capitalized, "\\proglang{Java}"]
  plain:     [keywords, not capitalized, Java]
date: "2020-09-16"
documentclass: jss
classoption: article
header-includes:
  - \usepackage{caption}
  - \captionsetup[table]{skip=.1em}
output: 
  bookdown::pdf_book:
    toc: false
    template: jss-template.tex
    keep_tex: true
    keep_md: true
    highlight: tango
    # highlight_bw: false
    md_extensions: +tex_math_dollars
    # latex_engine: xelatex
  html_document:
    keep_md: true
---



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
	
# Methodology

## Basics

- Dual problem.

- Uncertainty.  Either chi-square inversion for single parameters, or turns out that mode-quadrature (like MLE) also gives asymptotic normality, cf Qin-Lawless94.

- Might wish to provide gradient here, for optimization over theta.

## Support Correction

## Censoring

- What this is and how we estimate it with \pkg{flexEL} (EM algorithm, give the steps).

# Illustrations

## Mean Regression

- Show users how to create their own `Gfun` and use this from \proglang{R}.  

- Can include gradient-based optimization here as well.  I suggest to use `optim(method = "BFGS")` for optimization, not `nlm()` (it's just more complicated/annoying).  I would say use `nlminb()` which is way faster than `optim()`, but for some reason R Core suggests to use the other optimizers instead...

## Quantile Regression

- What it is, smoothness correction, present built-in function for this.

## Built-in Models

- We have a bunch of these: quantile vs mean regression, and location vs location-scale.

- Censoring as well.

- Should we put all this in the Methods section?  Have separate illustrations for some or all?

## Bayesian EL

- Good to illustrate with \proglang{Stan} NUTS algorithm, because we've implemented it already and \proglang{Stan} is a really good autodiff engine + best NUTS implementation.

- Can illustrate with the example from the JRSSB EL-HMC paper.

# Conclusion

- Everything below here are just tests for using JSS Markdown.

# Numbering

## Equations

Here's a reference to equation \@ref(eq:xyz).
\begin{equation}
x + y = z.
(\#eq:xyz)
\end{equation}
Here's an unnumbered equation:
\[
a/b = c.
\]
Here's some inline math: $y = \exp(x^2 - 1)$.

## Figures

Reference to Figure \@ref(fig:plot).  Also an example of embedding \proglang{R} code.

```r
R> x <- {function(y) {
+    y + 1:10
+  }}(3)
R> x
```

```
 [1]  4  5  6  7  8  9 10 11 12 13
```

```r
R> plot(x, pch = x, col = x)
```

![(\#fig:plot)A simple plot.](jss-draft_files/figure-latex/plot-1.pdf) 

## Tables

Reference to Table \@ref(tab:mtcars).
\begin{table}

\caption{(\#tab:mtcars)A table of the first 5 rows of the mtcars data.}
\centering
\begin{tabular}[t]{lrrrrrrrr}
\toprule
  & mpg & cyl & disp & hp & drat & wt & qsec & vs\\
\midrule
Mazda RX4 & 21.0 & 6 & 160 & 110 & 3.90 & 2.620 & 16.46 & 0\\
Mazda RX4 Wag & 21.0 & 6 & 160 & 110 & 3.90 & 2.875 & 17.02 & 0\\
Datsun 710 & 22.8 & 4 & 108 & 93 & 3.85 & 2.320 & 18.61 & 1\\
Hornet 4 Drive & 21.4 & 6 & 258 & 110 & 3.08 & 3.215 & 19.44 & 1\\
Hornet Sportabout & 18.7 & 8 & 360 & 175 & 3.15 & 3.440 & 17.02 & 0\\
\bottomrule
\end{tabular}
\end{table}

# Embedded Files

<!-- ## C++ File -->

\subsection[C++ File alpha + beta]{C++ File $\alpha + \beta$}

- LaTeX in heading doesn't seem to work.  [Here's](https://github.com/jgm/pandoc/issues/3555) how far I got with this.
- In fact, seems to be a pure LaTeX problem, in that `latexmk` compile of tex output takes several tries to get it right.
- Working solution: Use `\subsection[plaintext]{latex}` instead of Markdown section.
- New problem: `\tightlist` not defined...

```cpp
// test document embedding

// basic setter/getter for scalar double
class foo {
 private:
  double x;
 public:
  foo(double _x);
  void set_x(double _x);
  double get_x();
};

// constructor
inline foo::foo(double _x) {
  set_x(_x);
}

// getter
double get_x() {
  return x;
}

// setter
void set_x(double _x) {
  x = _x;
  return;
}
```

## R Markdown File

```bash
---
title: "Embedding R Markdown Files"
author: "Martin Lysy"
date: "`r Sys.Date()`"
output: 
  html_document:
    keep_md: true
---

## Introduction

Testing.
```
