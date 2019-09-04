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

# Load posterior samples of the trees & generate predictions
cor(colMeans(out2$y_preds), colMeans(out2$yhat))
plot(colMeans(out2$y_preds), colMeans(out2$yhat), col = z + 1)
abline(a=0, b=1)

cor(colMeans(out2$tau_preds), colMeans(out2$tau))
plot(colMeans(out2$tau_preds), colMeans(out2$tau), col = z + 1)
abline(a=0, b=1)

cor(colMeans(out2$mu_preds), colMeans(out2$mu))
plot(colMeans(out2$mu_preds), colMeans(out2$mu), col = z + 1)
abline(a=0, b=1)

yhat <- colMeans(out2$mu) + colMeans(out2$tau)*z
cor(yhat, colMeans(out2$yhat))
plot(yhat, colMeans(out2$yhat))

### Let's try adding mu_scale and tau_scale



# tau_preds_tmp <- out2$tau_preds %>% 
#   as.data.frame() %>% 
#   mutate_all(~./out2$mu_scale) %>% 
#   as.matrix()
#   
# ggplot(NULL, aes(x = colMeans(out2$y_preds), y = colMeans(out2$yhat), color = as.factor(z))) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_abline(slope = -1.03, intercept = -2.26) +
#   geom_abline(slope = -0.957, intercept = -1.57) +
#   geom_abline(slope = 0.618, intercept = 0.143)

# tbl <- data.frame(preds = colMeans(out2$y_preds), y = colMeans(out2$yhat), t = as.factor(z)) %>% 
#   mutate(group2 = case_when(preds < -0.75 ~ 1,
#                             preds > -0.75 & preds < 1 ~ 2,
#                             preds > 1 ~ 3))
# tbl %>% 
#   filter(t == 0) %>% 
#   group_by(group2) %>% 
#   do({
#     mod = lm(y ~ preds, data = .)
#     data.frame(Intercept = coef(mod)[1],
#                Slope = coef(mod)[2])
#   })

# tbl %>% 
#   filter(t == 1) %>% 
#   do({
#     mod = lm(y ~ preds, data = .)
#     data.frame(Intercept = coef(mod)[1],
#                Slope = coef(mod)[2])
#   })
