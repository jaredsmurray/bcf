set.seed(1)

p <- 3 # two control variables and one effect moderator
n <- 10
n_burn <- 10
n_sim <- 15


x <- matrix(rnorm(n*p), nrow=n)

weights <- 1.0*rep(1, n)


# create targeted selection, whereby a practice's likelihood of joining the intervention (pi) is related to their expected outcome (mu)
q <- -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2])) -0.1

# generate treatment variable
pi <- pnorm(q)
z <- rbinom(n,1,pi)

# tau is the true treatment effect. It varies across practices as a function of
# X3, the effect moderator
tau <- 1/(1 + exp(-x[,3]))

# generate the response using q, tau and z
mu <- (q + tau*z)

# set the noise level relative to the expected mean function of Y
sigma <- diff(range(q + tau*pi))/8

# draw the response variable with additive error
y <- mu + sigma*rnorm(n)



out2 <- bcf2::bcf(y          = y,
                  z          = z,
                  x_control  = x,
                  x_moderate = x,
                  pihat      = pi,
                  nburn      = n_burn,
                  nsim       = n_sim,
                  w          = weights, 
                  update_interval = 100)

                  
cat("BCF run complete\n")

