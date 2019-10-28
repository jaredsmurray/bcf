set.seed(1)

p <- 3 # two control variables and one effect moderator
n <- 1000
n_burn <- 100
n_sim <- 150


x <- matrix(rnorm(n*p), nrow=n)



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

weights <- 1000.0*rep(1, n)

set.seed(1)
out <- bcf2::bcf(y               = y,
                  z               = z,
                  x_control       = x,
                  x_moderate      = x,
                  pihat           = pi,
                  nburn           = n_burn,
                  nsim            = n_sim,
                  w               = weights,
                  random_seed     = 1,
                  update_interval = 100)

originatOut = readRDS("examples/output_original.rds")


mean_square_error <- function (x,y){
  mean((x-y)^2)
}

assess_closeness <- function(x,y, title){
  cat("Assessing Cloesness of ", title, "\n")
  print("Correlation")
  print(cor(x,y))
  print("MSE")
  print(mean_square_error(x,y))
  plot(x, y, col = z + 1, main=title)
  abline(a=0, b=1)
}

assess_closeness(colMeans(originatOut$yhat), colMeans(out$yhat),'yhat')

assess_closeness(colMeans(originatOut$tau), colMeans(out$tau),'tau')

assess_closeness(colMeans(originatOut$mu), colMeans(out$mu),'mu')



