set.seed(1)

p <- 3 # two control variables and one effect moderator
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

bcf_out <- bcf::bcf(y            = y,
                 z                = z,
                 x_control        = x,
                 x_moderate       = x,
                 pihat            = pi,
                 nburn            = n_burn,
                 nsim             = n_sim,
                 w                = weights,
                 n_chains         = 4)

# saveRDS(bcf_out, file = "examples/my_data.rds")

