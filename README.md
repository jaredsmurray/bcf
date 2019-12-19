# Bayesian Causal Forests

Welcome to the `BCF` site! This page provides hands-on examples of how to conduct Bayesian causal forest (BCF) analyses. You can find methodological details on the underlying modeling approach in the original BCF paper: [Hahn, Murray, and Carvalho 2017](https://arxiv.org/pdf/1706.09523.pdf).

## Why BCF?

BCF is a cutting-edge causal inference methodology that builds on Bayesian Additive Regression Trees (BART, [Chipman, George, and McCulloch 2010](https://projecteuclid.org/euclid.aoas/1273584455)). More immediately, it builds on the application of BART to causal inference ([Hill 2011](https://www.tandfonline.com/doi/abs/10.1198/jcgs.2010.08162)).  BART and BCF both combine Bayesian regularization with regression trees to provide a highly flexible response surface that, thanks to the Bayesian regularizing priors, is not overfit to the training data.  BCF further extends BART's flexibility by specifying different models for relationships between (a) covariates and the outcome and (b) covariates and the treatment effect.

BCF performs remarkably well in simulation and has led the pack at recent rigorous causal inference competitions, such as those held at the Atlantic Causal Inference Conference. 

[MMF: the following list of bullets doesn't seem like it belongs here, on the front page for the whole package, does it? If you agree, OK to move this to the changelog?]This implementation further extends existing BCF functionality by:

- allowing for heteroskedastic error
- automating multi-chain, multi-core implementations
- providing a suite of convergence diagnostic functions via the `coda` package
- accelerating some underlying computations, resulting in shorter runtimes

## Getting Started

If you are just getting started with BCF, we recommend starting with the tutorial vignettes on the "Articles" page of this website.

## Installation

This package requires compilation, so make sure you have Rtools properly installed -- details [here](https://cran.r-project.org/bin/windows/Rtools/).

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

## Example

[MMF: This example seems redundant with, but less clear than, the vignette. Cut? If we don't cut it, I'll go back through and review more closely.]

```{r}
set.seed(1)

p <- 3 # two control variables and one effect moderator
n <- 1000
n_burn <- 1000
n_sim <- 1000

x <- matrix(rnorm(n*p), nrow=n)


# create targeted selection, whereby a unit's likelihood of joining the intervention 
# (pi) is related to its expected outcome (mu)
q <- -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2])) -0.1

# generate treatment variable
pi <- pnorm(q)
z <- rbinom(n,1,pi)

# tau is the true treatment effect, which varies across units as a function of
# X3, the effect moderator
tau <-  1/(1 + exp(-x[,3]))

mu <- q

# generate the response using q, tau and z
y_noiseless <- mu + tau*z

# set the noise level relative to the expected mean function of Y
sigma <- diff(range(mu + tau*pi))/8

# draw the response variable with additive error
y <- y_noiseless + sigma*rnorm(n)

weights <- 1000.0*rep(1, n) # [MMF: doesn't affect sigma --> not needed?]

bcf_out <- bcf2::bcf(y            = y,
                 z                = z,
                 x_control        = x,
                 x_moderate       = x,
                 pihat            = pi,
                 nburn            = n_burn,
                 nsim             = n_sim,
                 w                = weights,
                 n_chains         = 4,
                 n_chain_clusters = 2,
                 random_seed      = 1,
                 update_interval  = 1)

```

```{r}
bcf2::summarise_bcf(bcf_out)

# [MMF: Even if we keep this example on the home page, I don't think we need the following checks here, do we? Seems like TMI, since we haven't yet clarified what mu/tau/etc are...]
mean_square_error <- function (x,y){
  mean((x-y)^2)
}

assess_closeness <- function(x,y, title){
  cat("Assessing Cloesness of ", title, "\n")
  
  mse = mean_square_error(x,y)
  
  print("MSE")
  print(mse)
  
  print("Error")
  print(sqrt(mse)/abs(mean(x)))
}

mu_correct = bcf_out$yhat - t(t(bcf_out$tau)*z)

print("Y Mean")
print(mean(y))

print("Tau Mean")
print(mean(tau))

print("mu Mean")
print(mean(mu))

assess_closeness(bcf_out$mu,mu_correct,'mu_compare')
```
