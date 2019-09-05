library(bcf2)
library(tidyverse)
set.seed(1)

p <- 3 # two control variables and one effect moderator
n <- 1000
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
tau <- 1/(1 + exp(-x[,3]))

# generate the response using q, tau and z
mu <- (q + tau*z)

# set the noise level relative to the expected mean function of Y
sigma <- diff(range(q + tau*pi))/8

# draw the response variable with additive error
y <- mu + sigma*rnorm(n)

out2 <- bcf2::bcf(y               = y,
                  z               = z,
                  x_control       = x,
                  x_moderate      = x,
                  pihat           = pi,
                  z_pred          = z,
                  x_pred_moderate = x,
                  x_pred_control  = x,
                  pi_pred         = pi,
                  nburn           = n_burn,
                  nsim            = n_sim,
                  w               = weights, 
                  update_interval = 1,
                  ntree_moderate  = 3,
                  ntree_control   = 3,
                  verbose         = TRUE,
                  use_muscale     = TRUE,
                  use_tauscale    = TRUE)

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

assess_closeness(colMeans(out2$y_preds), colMeans(out2$yhat),'yhat')

assess_closeness(colMeans(out2$tau_preds), colMeans(out2$tau),'tau')

assess_closeness(colMeans(out2$mu_preds), colMeans(out2$mu),'mu')

normal_mu = colMeans(out2$mu)

df <- data.frame(mu_preds = colMeans(out2$mu_preds))

mod = lm(normal_mu ~ mu_preds, data = df)
print(mod)

assess_closeness(colMeans(out2$mu_preds)+0.1348, colMeans(out2$mu),'mu_mod')

yhat_preds_2 = colMeans(out2$mu_preds) + colMeans(out2$tau_preds)*z

assess_closeness(yhat_preds_2, colMeans(out2$yhat),'yhat_2')

z_matrix = matrix(1,length(out2$tau_scale),1)%*%z

yhat_preds_3 = out2$mu_preds + out2$tau_preds*z_matrix

assess_closeness(colMeans(yhat_preds_3), colMeans(out2$yhat),'yhat_3')

A = matrix( 
    c(2, 4, 3, 1, 5, 7), # the data elements 
     nrow=2,              # number of rows 
     ncol=3,              # number of columns 
     byrow = TRUE) 

x = c(1,2,0)

b = c(1,2,3,4)