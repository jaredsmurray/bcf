set.seed(1)
library(rbenchmark)

library(profvis)

p <- 10 # two control variables and one effect moderator
n <- 100
n_burn <- 100
n_sim <- 100

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

bcf_run_name = function(n_threads, save_trees){
  sprintf('t%d_%s_tree',n_threads, save_trees)
}

bcf_wrapper = function(n_threads, save_trees){
  save_tree_directory <- if(save_trees=='yes') "." else NULL
  
  bcf_out <- bcf::bcf(y                   = y,
                      z                   = z,
                      w                   = weights,
                      x_control           = x,
                      x_moderate          = x,
                      pihat               = pi,
                      nburn               = n_burn,
                      nsim                = n_sim,
                      save_tree_directory = save_tree_directory,
                      n_chains            = 1,
                      n_cores             = 1,
                      n_threads           = n_threads)
  
  
}

# bcf_wrapper(1,'no') 

# p = profvis({ bcf_wrapper(4,'no') })
  

tests = list()
tests[[bcf_run_name(8,'yes')]] = expression(bcf_wrapper(8,'yes'))
tests[[bcf_run_name(8,'no')]]  = expression(bcf_wrapper(8,'no'))
tests[[bcf_run_name(4,'yes')]] = expression(bcf_wrapper(4,'yes'))
tests[[bcf_run_name(4,'no')]]  = expression(bcf_wrapper(4,'no'))
tests[[bcf_run_name(2,'yes')]] = expression(bcf_wrapper(2,'yes'))
tests[[bcf_run_name(2,'no')]]  = expression(bcf_wrapper(2,'no'))
tests[[bcf_run_name(1,'yes')]] = expression(bcf_wrapper(1,'yes'))
tests[[bcf_run_name(1,'no')]]  = expression(bcf_wrapper(1,'no'))
 
df = do.call(benchmark, c(tests, list(replications=1)))

show(df)

p = profvis({ bcf_wrapper(4,'no') })


