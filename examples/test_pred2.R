library(bcf)
library(tidyverse)
set.seed(1)

p <- 3 # two control variables and one effect moderator
n <- 100
n_burn <- 100
n_sim <- 150


x <- matrix(rnorm(n*p), nrow=n)

weights <- 1.0*rep(1, n)


# create targeted selection, whereby a practice's likelihood of joining the intervention (pi) is related to their expected outcome (mu)
q <- -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2])) -0.1

# generate treatment variable
pi <- pnorm(q)
z <- rbinom(n,1,pi)

# tau is the true treatment effect. It varies across practices as a function of
# X3, the effect moderator
tau <- (1/(1 + exp(-x[,3])))

mu <- q

# generate the response using q, tau and z
y_noiseless <- mu + tau*z

# set the noise level relative to the expected mean function of Y
sigma <- diff(range(mu + tau*pi))/8

# draw the response variable with additive error
y <- y_noiseless + sigma*rnorm(n)

out2 <- bcf::bcf(y               = y,
                  z               = z,
                  x_control       = x,
                  x_moderate      = x,
                  pihat           = pi,
                  nburn           = n_burn,
                  nsim            = n_sim,
                  w               = weights,
                  n_chains        = 2,
                  nthin           = 3,
                  update_interval = 1)

cat("BCF run complete\n")



cat("Starting Prediction \n")

pred_out = predict(object=out2,
                   x_predict_control=x,
                   x_predict_moderate=x,
                   pi_pred=pi,
                   z_pred=z,
                   save_tree_directory = '..')



cat("Predictions Compelete\n")


mean_square_error <- function (x,y){
  mean((x-y)^2)
}

assess_closeness <- function(x,y, title){
  cat("Assessing Cloesness of ", title, "\n")
  print("Correlation")
  print(cor(x,y))
  
  mse = mean_square_error(x,y)
  
  print("MSE")
  print(mse)
  
  print("Error")
  print(sqrt(mse)/abs(mean(x)))
  plot(x, y, col = z + 1, main=title)
  abline(a=0, b=1)
}

print("Y Mean")
print(mean(y))

print("Tau Mean")
print(mean(tau))

print("mu Mean")
print(mean(mu))

assess_closeness(colMeans(pred_out$yhat), colMeans(out2$yhat),'yhat')
assess_closeness(colMeans(pred_out$tau), colMeans(out2$tau),'tau')
assess_closeness(colMeans(pred_out$mu), colMeans(out2$mu),'mu')

# summarise_bcf(out2)
# summarise_bcf(pred_out)


