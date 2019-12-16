# bcf 2.0.0 (bcf2)

This page provides more details on the latest updates to bcf.

## Major changes

### Weights

The original version of {bcf} does not allow for weights, which we often use in practical applications to account for heteroskedasticity. Where the original BCF model was specified as:

$y_i \sim N(\mu(x_i) + \tau(x_i) z_i, \sigma^2)$,

which assumes that all outcomes $y_i$ have the same variance $\sigma^2$, in the extended version we can relax this assumption to allow the variance to reflect the uncertainty in the $y_i$:

$y_i∼N(μ(x_i)+τ(x_i ) z_i,σ^2/w_i )$

We changed several parts of the code to incorporate these weights:

* Code to calculate sufficient statistics
* Code to update leaf node means
* Code to update variance across leaf node means

### Automating multichain processing

It is useful in Bayesian analysis to produce different runs of the same model, with different starting values, as a way of assessing convergence; if the different runs produce drastically different posterior distributions, it is a sign that the model has not converged fully.  In this version of {bcf} we have automated multichain processing and incorporated key MCMC diagnostics from the {coda} package, including effective sample sizes and the Gelman-Rubin statistic ("R hat").

*Question for Peter*: is the multichain processing also multicore?  If so, we should talk it up. 

### Speed-ups

Finally, our implementation parallelizes some steps of the sampling procedure to maximize efficiency in a multicore environment.  Our testing shows that these enhancements have reduced runtimes by [*Peter please fill in here*].
