# Bayesian Causal Forests

Welcome to the BCF site! This page provides more details on this implementation of Bayesian causal forests (BCF).

## Getting Started

If you are just getting started with BCF, we recommend starting with the tutorial vignette and the examples throughout the package documentation.

Something else that may be useful is the original BCF paper here: [https://arxiv.org/pdf/1706.09523.pdf](https://arxiv.org/pdf/1706.09523.pdf).

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

## Examples

```{r}
set.seed(1)

p <- 3 # two control variables and one effect moderator
n <- 1000
n_burn <- 1000
n_sim <- 1000

x <- matrix(rnorm(n*p), nrow=n)


# create targeted selection, whereby a practice's likelihood of joining the intervention (pi) is related to their expected outcome (mu)
q <- -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2])) -0.1

# generate treatment variable
pi <- pnorm(q)
z <- rbinom(n,1,pi)

# tau is the true treatment effect. It varies across practices as a function of
# X3, the effect moderator
tau <-  1/(1 + exp(-x[,3]))

mu <- q

# generate the response using q, tau and z
y_noiseless <- mu + tau*z

# set the noise level relative to the expected mean function of Y
sigma <- diff(range(mu + tau*pi))/8

# draw the response variable with additive error
y <- y_noiseless + sigma*rnorm(n)

weights <- 1000.0*rep(1, n)

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

mean_square_error <- function (x,y){
  mean((x-y)^2)
}
```

```{r}
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

## Latest package updates: weights

The original version of {bcf} does not allow for weights. 
In the new version, the model is specified as follows:

`y_i∼N(μ(x_i)+τ(x_i ) z_i,σ^2/w_i )`

We changed several parts of the code to incorporate these weights:

* Code to calculate sufficient statistics
* Code to update leaf node means
* Code to update variance across leaf node means
