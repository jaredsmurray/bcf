# Bayesian Causal Forests

Welcome to the `bcf` site! This page provides hands-on examples of how to conduct Bayesian causal forest (BCF) analyses. You can find methodological details on the underlying modeling approach in the original BCF paper: [Hahn, Murray, and Carvalho 2020](https://projecteuclid.org/journals/bayesian-analysis/volume-15/issue-3/Bayesian-Regression-Tree-Models-for-Causal-Inference--Regularization-Confounding/10.1214/19-BA1195.full).

## Why BCF?

BCF is a cutting-edge model for causal inference that builds on Bayesian Additive Regression Trees (BART, [Chipman, George, and McCulloch 2010](https://projecteuclid.org/euclid.aoas/1273584455)). BART and BCF both combine Bayesian regularization with regression trees to provide a highly flexible response surface that, thanks to regularization from prior distributions, does not overfit to the training data. BCF extends BART's flexibility by specifying different models for relationships between (1) covariates and the outcome and (2) covariates and the treatment effect, and regularizing the treatment effect directly.

BCF performs remarkably well in simulation and has led the pack at recent rigorous causal inference competitions, such as those held at the Atlantic Causal Inference Conference (see, for example, [Dorie et al. 2019](https://projecteuclid.org/euclid.ss/1555056030)).

## Getting Started

If you are just getting started with `bcf`, we recommend beginning with the tutorial vignettes.

## Installation

This package requires compilation, so make sure you have Rtools properly installed if you are on Windows -- see [this site](https://cran.r-project.org/bin/windows/Rtools/) for details.

Install the latest release from CRAN:

```{r}
install.packages("bcf")
```

Install the latest development version from GitHub:

```{r}
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("jaredsmurray/bcf")
```