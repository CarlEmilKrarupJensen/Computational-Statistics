Computational-Statistics
This repository contains solutions to the four exam assignments for the course on Computational Statistics at the University of Copenhagen.

Each project is organized in its own folder:

```text
assignment_1/
assignment_2/
assignment_3/
assignment_4/
```

The projects combine statistical methodology with computational implementation, simulation studies, validation, benchmarking, profiling, and comparison with existing R functionality.

---

## Repository structure

Each assignment folder contains the project code, reports, supporting functions, tests, and presentation material for that assignment. The Quarto reports contain the numerical investigations, computational comparisons, and figures.

The general workflow across the projects is:

- implementation of the statistical method,
- correctness checks and testing,
- simulation or empirical evaluation,
- comparison with alternative methods or R implementations,
- profiling and benchmarking,
- discussion of accuracy, speed, and memory trade-offs.

---

## Assignment 1 — Kernel Density Estimation and Bandwidth Selection

This project implements a kernel density estimator using the Epanechnikov kernel together with several bandwidth-selection approaches.

The project includes:

- AMISE plug-in bandwidth selection,
- likelihood cross-validation,
- least-squares cross-validation,
- simulation studies for normal, skewed, and multimodal densities,
- comparison with oracle bandwidths,
- comparison with selected built-in R bandwidth selectors,
- correctness testing with `testthat`,
- profiling of computational bottlenecks,
- benchmarking of runtime and memory usage,
- alternative KDE implementations using matrix, chunked, and binned approaches,
- comparison with R's `density()` implementation,
- analysis of speed-memory-accuracy trade-offs.


---

## Assignment 2 — Project B: Importance Sampling for Finite-Horizon Ruin

This project estimates the probability that a random walk crosses a ruin threshold within a fixed time horizon, using exponential tilting to make rare events easier to simulate.

The project includes:

- crude Monte Carlo, full-path importance sampling, and stopped-path importance sampling,
- likelihood-ratio weights and unbiasedness arguments for the importance sampling estimators,
- pilot selection of the tilting parameter and independent validation,
- comparison of estimator variance, Monte Carlo uncertainty, and relative efficiency,
- analysis of variance-minimizing versus computing-cost-adjusted tilts,
- early stopping and survivor-index bookkeeping to reduce simulation work,
- correctness testing, profiling, and runtime benchmarking,
- comparison of the computing time required to achieve a target precision.

### Key questions

- How much can exponential tilting improve precision for a rare ruin event?
- How does stopping at ruin affect runtime and estimator variance?
- Does the variance-minimizing tilt also minimize the time needed for a given precision?

---

## Assignment 3 — Project A: Student's t Estimation with EM

This project estimates the location and scale of a Student's t distribution with fixed degrees of freedom, using its latent-variable representation to derive an expectation-maximization (EM) algorithm.

The project includes:

- derivation of the complete-data likelihood and EM updates,
- maximum-likelihood estimation using EM, BFGS, and Newton methods,
- validation of likelihoods, derivatives, information calculations, and optimization routines,
- contamination studies examining robustness and optimizer reliability,
- Fisher information, standard errors, and confidence intervals,
- Monte Carlo evaluation of estimation uncertainty and interval coverage,
- comparison of convergence, likelihood agreement, and runtime across datasets,
- profiling and targeted computational improvements.

### Key questions

- How do the EM weights explain robustness to outlying observations?
- How do EM, BFGS, and Newton compare in speed and convergence reliability?
- How well does Fisher information predict uncertainty and confidence-interval coverage?

---

## Assignment 4 — Project B: Gradient Methods for Dose-Response Curves

This project fits a four-parameter logistic dose-response curve by nonlinear least squares and compares deterministic and stochastic gradient methods.

The project includes:

- analytical gradients and Hessians, with derivative checks and local geometry diagnostics,
- full-gradient descent with Armijo backtracking,
- stochastic gradient descent, Polyak and Nesterov momentum, and AdaGrad,
- sensitivity studies for starting values, learning rates, sample size, and local conditioning,
- comparison of fitted curves and objective values against a multi-start BFGS reference,
- accuracy-versus-time comparisons, profiling, and memory measurements,
- exact compression of repeated-dose data,
- comparison of R and compiled C++ implementations of stochastic gradient descent.

### Key questions

- When do momentum and adaptive scaling improve convergence?
- How do starting values, learning rates, and local geometry affect numerical accuracy?
- How do repeated-dose compression and compiled implementations change runtime and memory costs?

---

## Tools

The projects are primarily implemented in R and use tools such as:

- Quarto for reports and presentation material,
- `testthat` for correctness testing,
- `bench` for benchmarking,
- `Rprof()` and related profiling tools for identifying bottlenecks,
- `ggplot2` for visualization,
- Rcpp and C++ for selected performance-critical routines.

The emphasis throughout the repository is not only on obtaining statistically correct results, but also on understanding the computational cost and trade-offs of different implementations.








