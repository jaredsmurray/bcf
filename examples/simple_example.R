set.seed(1)

p <- 3 # two control variables and one effect moderator
n <- 10
n_burn <- 15
n_sim <- 15


x <- matrix(rnorm(n*p), nrow=n)



# create targeted selection, whereby a practice's likelihood of joining the intervention (pi) is related to their expected outcome (mu)
q <- -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2])) -0.1

# generate treatment variable
pi <- pnorm(q)
z <- rbinom(n,1,pi)

# tau is the true treatment effect. It varies across practices as a function of
# X3, the effect moderator
tau <- 1 + 1/(1 + exp(-x[,3]))

# generate the response using q, tau and z
mu <- (10 + 10*q + tau*z)

# set the noise level relative to the expected mean function of Y
sigma <- diff(range(q + tau*pi))/8

# draw the response variable with additive error
y <- mu + sigma*rnorm(n)

weights <- 1000.0*rep(1, n)

bcf_out <- bcf2::bcf(y          = y,
                 z          = z,
                 x_control  = x,
                 x_moderate = x,
                 pihat      = pi,
                 nburn      = n_burn,
                 nsim       = n_sim,
                 w          = weights,
                 n_chains = 1,
                 random_seed = 1,
                 update_interval = 100)

cat("BCF Run Complete \n")

# bcf2::summarise_bcf(bcf_out)

# coda::traceplot(bcf_out$chains)


mean_square_error <- function (x,y){
  mean((x-y)^2)
}

assess_closeness <- function(x,y, title){
  cat("Assessing Cloesness of ", title, "\n")
  # print("Correlation")
  # print(cor(x,y))
  
  mse = mean_square_error(x,y)
  
  print("MSE")
  print(mse)
  
  print("Error")
  print(sqrt(mse)/abs(mean(x)))
  
  # plot(x, y, col = z + 1, main=title)
  # abline(a=0, b=1)
}


Tm = t( t(bcf_out$tau) * (1/ (bcf_out$sdy*(bcf_out$b1 - bcf_out$b0))) )

Tc = t( t(bcf_out$mu - bcf_out$muy) * (1/(bcf_out$sdy*bcf_out$mu_scale)) )


mu_theory = bcf_out$muy + bcf_out$sdy*( t(t(Tc)*bcf_out$mu_scale) + t(t(Tm)*bcf_out$b0) )

mu_correct = bcf_out$yhat - t(t(bcf_out$tau)*z)


assess_closeness(mu_theory,mu_correct,'mu_compare')


df <-        data.frame("z"           = z,
                        "y"           = y,
                        "y_hat"       = colMeans(bcf_out$yhat),
                        "tau"         = tau,
                        "tau_bar"     = colMeans(bcf_out$tau),
                        "mu_correct"  = colMeans(mu_correct),
                        "mu_theory"   = colMeans(mu_theory),
                        "mu_diff"     = colMeans(mu_correct) - colMeans(mu_theory),
                        "Tm"          = colMeans(Tm),
                        "Tc"          = colMeans(Tc))

# print(round(t(df)))

print(t(df))

diff = mu_theory - mu_correct
