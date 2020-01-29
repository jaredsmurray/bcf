# bcf 2.0.0

This page provides more details on the latest updates to `bcf`.

## Major changes

This implementation extends existing `bcf` functionality by:

- allowing for heteroskedastic errors
- automating multi-chain, multi-core implementations
- providing a suite of convergence diagnostic functions via the `coda` package
- accelerating some underlying computations, resulting in shorter runtimes
- providing a function to predict treatment effects based on an existing model using new data

### Weights

The original version of `bcf` does not allow for weights, which are often used in practical applications to account for heteroskedasticity. Where the original BCF model was specified as:

y<sub>i</sub> &sim; N(&mu;(x<sub>i</sub>) + &tau;(x<sub>i</sub>) z<sub>i</sub>, &sigma;<sup>2</sup>),

which assumes that all outcomes y<sub>i</sub> have the same variance &sigma;<sup>2</sup>, in the extended version we can relax this assumption to allow for heteroskedasticity in y<sub>i</sub>:

y<sub>i</sub> &sim; N(&mu;(x<sub>i</sub>) + &tau;(x<sub>i</sub>) z<sub>i</sub>, &sigma;<sup>2</sup>/w<sub>i</sub>)

Incorporating weights impacts several parts of the code, including the computation of:

* sufficient statistics
* leaf node means
* leaf node means variance
* error variance (sigma)

### Automating multichain processing in parallel

In Bayesian analysis, it is useful to produce different runs of the same model -- with different starting values -- as a way of assessing convergence. If the different runs produce drastically different posterior distributions, it is a sign that the model has not converged fully.  In this version of `bcf` we have automated multichain processing and incorporated key MCMC diagnostics from the `coda` package, including effective sample sizes and the Gelman-Rubin statistic ("R hat"). In addition, all runs happen in parallel, on different cores, to provide all these extra benefits without increasing runtime much.

### Within-chain parallel processing

Finally, our implementation conducts some steps of the sampling procedure in parallel to maximize computational efficiency. Our testing shows that these enhancements have reduced runtimes by around 50%, across various experimental conditions.

### Implementing a prediction method

In this version of the package, we have incorporated code from Jared Murray (one of the original package authors) to predict the treatment effect for a new set of units. Once users have produced a satisfactory `bcf` run (using training data), they can use this fitted `bcf` object to predict on a new set of test data. This is possible even with runs that have multiple chains.
