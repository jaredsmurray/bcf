# bcf 2.0.0 (bcf2)

This page provides more details on the latest updates to bcf.

## Major changes

This implementation further extends existing BCF functionality by:

- allowing for heteroskedastic error
- automating multi-chain, multi-core implementations
- providing a suite of convergence diagnostic functions via the `coda` package
- accelerating some underlying computations, resulting in shorter runtimes

### Weights

The original version of {bcf} does not allow for weights, which we often use in practical applications to account for heteroskedasticity. Where the original BCF model was specified as:

y<sub>i</sub> &sim; N(&mu;(x<sub>i</sub>) + &tau;(x<sub>i</sub>) z<sub>i</sub>, &sigma;<sup>2</sup>),

which assumes that all outcomes y<sub>i</sub> have the same variance &sigma.<sup>2</sup>, in the extended version we can relax this assumption to allow the variance to reflect the uncertainty in the y<sub>i</sub>:

y<sub>i</sub> &sim; N(&mu;(x<sub>i</sub>) + &tau;(x<sub>i</sub>) z<sub>i</sub>, &sigma;<sup>2</sup>/w<sub>i</sub>)

We changed several parts of the code to incorporate these weights:

* Code to calculate sufficient statistics
* Code to update leaf node means
* Code to update variance across leaf node means

### Implementing a prediction method

In this version of the package, we have incorporated code from Jared Murray (one of the original package authors) to predict values of the outcome based on a new set of covariates. Once users have produced a satisfactory BCF run (using training data), they are able to use this run to predict on a new set of (test) data. This is possible even with runs that have multiple chains.

### Automating multichain processing in parallel

It is useful in Bayesian analysis to produce different runs of the same model, with different starting values, as a way of assessing convergence; if the different runs produce drastically different posterior distributions, it is a sign that the model has not converged fully.  In this version of {bcf} we have automated multichain processing and incorporated key MCMC diagnostics from the {coda} package, including effective sample sizes and the Gelman-Rubin statistic ("R hat"). In addition, all runs happen in parallel, on different cores, so as to provide all these extra benefits without much cost to the timing of the runs.

### Within chain parallel processing

Finally, our implementation parallelizes some steps of the sampling procedure to maximize efficiency in a multicore environment.  Our testing shows that these enhancements have reduced runtimes by around 2x, across a variety of experimental conditions.
