# Computational-Statistics
This repository contains solutions to the four exam assignments for the course on Computational Statistics at the University of Copenhagen.

Each project is organized in its own folder:

```text
assignment1/
assignment2/
assignment3/
assignment4/
```

The projects combine statistical methodology with computational implementation, simulation studies, validation, benchmarking, profiling, and comparison with existing R functionality.

---

## Repository structure

Each assignment folder contains the code, report, figures, and supporting functions and presentation material.

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

## Assignment 2 — [Project title]

**Topic:** [One-sentence description of the project]

Main components:

- [method or model implemented],
- [simulation or empirical study],
- [computational aspect],
- [comparison with alternatives],
- [main statistical or computational question].

### Key questions

- [Question 1]
- [Question 2]
- [Question 3]

---

## Assignment 3 — [Project title]

**Topic:** [One-sentence description of the project]

Main components:

- [method or model implemented],
- [simulation or empirical study],
- [computational aspect],
- [comparison with alternatives],
- [main statistical or computational question].

### Key questions

- [Question 1]
- [Question 2]
- [Question 3]

---

## Assignment 4 — [Project title]

**Topic:** [One-sentence description of the project]

Main components:

- [method or model implemented],
- [simulation or empirical study],
- [computational aspect],
- [comparison with alternatives],
- [main statistical or computational question].

### Key questions

- [Question 1]
- [Question 2]
- [Question 3]

---

## Tools

The projects are primarily implemented in R and use tools such as:

- Quarto for reports,
- `testthat` for correctness testing,
- `bench` for benchmarking,
- `Rprof()` and related profiling tools for identifying bottlenecks,
- `ggplot2` for visualization.

The emphasis throughout the repository is not only on obtaining statistically correct results, but also on understanding the computational cost and trade-offs of different implementations.






