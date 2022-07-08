## R CMD check results
Duration: 10m 43.7s

> checking files in ‘vignettes’ ... NOTE
  The following directory looks like a leftover from 'knitr':
    ‘figure’
  Please remove from your package.

> checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

0 errors v | 0 warnings v | 2 notes x

R CMD check succeeded

### > checking for GNU extensions in Makefiles ... NOTE GNU make is a SystemRequirements
    Because our package uses RcppParallel, we inherit a dependence on GNU make. We can’t eliminate this note. 
    
### > checking files in ‘vignettes’ ... NOTE The following directory looks like a leftover from 'knitr': ‘figure’ Please remove from your package.
    This directory and these files are cached figures for rebuilding vignettes with long run times. It is not leftover.
    
