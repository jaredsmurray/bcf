par(mfrow=c(2,3))

# sample size
n <- 200

# prognostic function mu
mu <- function(x){
result = 0.5 - 0.25*x + 0.25*(cos(5*x) - sin(5*x+1))*dnorm(x,-0.25,0.2)
result = 20*result
return(result)}

xgrid = seq(-1,1,length.out=200)
plot(xgrid,mu(xgrid),lwd=2,type='l')


# confounder variables
x1 <- runif(n,-1/2,1/2)  
x2 <- runif(n, -1/2,1/2)


# treatment effect function
tau <- function(x){
  return(3 + 3*log(6*x+3.01)/10)
}

xgrid = seq(-0.5,0.5,length.out=200)
plot(xgrid,tau(xgrid),lwd=2,type='l')


# treatment variable
lambda <- 1
muhat = mu(x1 + x2) - mean(mu(x1+x2)) + 0.1*rnorm(n)
intercept = 1
pi = 1/(1 + exp(intercept-lambda*muhat + 0*x2))
pi = 0.1 + 0.8*pi

hist(pi,col='gray')

z <-  rbinom(n,1,pi)  

EY = mu(x1 + x2) + tau(x1)*z
sig = 1*sd(tau(x1))

# response variable
error <- sig*rnorm(n)
y <- EY + error 

hist(error,col='gray')

# fit bcf, not using pi
 fit_default = bcf(y, z, x_control = cbind(x1,x2), 
                  nburn = 500, nsim = 500, pihat = matrix(pi),update_interval = 100000, include_pi = "none")



# extract posterior mean of CATEs
tauhat <- colMeans(fit_default$tau)


# fit "just another covariate" with BART, not using pi
 dbfit <- bart(cbind(x1,x2,z),y,cbind(c(x1,x1),c(x2,x2),c(rep(1,n),rep(0,n))),verbose=FALSE)

# extract posterior mean of CATEs
dbtauhat = dbfit$yhat.test.mean[1:n] - dbfit$yhat.test.mean[1:n + n]

# not using pi
  bf1 <- bart(cbind(x1[z==1],x2[z==1]),y[z==1],cbind(x1,x2),verbose=FALSE)
  bf0 <- bart(cbind(x1[z==0],x2[z==0]),y[z==0],cbind(x1,x2),verbose=FALSE)

# extract posterior mean of CATEs
btauhat = bf1$yhat.test.mean[1:n] - bf0$yhat.test.mean[1:n]

# plot estimates
ests = c(btauhat,tauhat,dbtauhat)
plot(x1,tau(x1),ylim=c(2.5,4))
points(x1,btauhat,col='gray',pch=20)
points(x1,tauhat,pch=20,col='green')
points(x1,dbtauhat,pch=20,col='orange')


rmse.spline <- sqrt(mean((tau(x1)-btauhat)^2))
rmse.bart <- sqrt(mean((tau(x1)-dbtauhat)^2))
rmse.bcf <- sqrt(mean((tau(x1)-tauhat)^2))

print("Without the propensity transformation")
print(cbind(rmse.spline,rmse.bart,rmse.bcf))

# fit bcf, using true propensity function
fit_default = bcf(y, z, x_control = cbind(x1,x2), 
                  nburn = 500, nsim = 500, pihat = matrix(pi),update_interval = 100000)


# extract posterior mean of CATEs
tauhat <- colMeans(fit_default$tau)


# fit "just another covariate" with BART, using pi
dbfit <- bart(cbind(x1,x2,pi,z),y,cbind(c(x1,x1),c(pi,pi),c(x2,x2),c(rep(1,n),rep(0,n))), verbose = FALSE)

# extract posterior mean of CATEs
dbtauhat = dbfit$yhat.test.mean[1:n] - dbfit$yhat.test.mean[1:n + n]

# fit "separate regressions" with BART, using pi
bf1 <- bart(cbind(x1[z==1],x2[z==1],pi[z==1]),y[z==1],cbind(x1,x2,pi),verbose=FALSE)
bf0 <- bart(cbind(x1[z==0],x2[z==0],pi[z==0]),y[z==0],cbind(x1,x2,pi),verbose=FALSE)

# extract posterior mean of CATEs
btauhat = bf1$yhat.test.mean[1:n] - bf0$yhat.test.mean[1:n]

# plot estimates
ests = c(btauhat,tauhat,dbtauhat)
plot(x1,tau(x1),ylim=c(2.5,4))
points(x1,btauhat,col='gray',pch=20)
points(x1,tauhat,pch=20,col='green')
points(x1,dbtauhat,pch=20,col='orange')



rmse.spline <- sqrt(mean((tau(x1)-btauhat)^2))
rmse.bart <- sqrt(mean((tau(x1)-dbtauhat)^2))
rmse.bcf <- sqrt(mean((tau(x1)-tauhat)^2))
print("With the propensity transformation")
print(cbind(rmse.spline,rmse.bart,rmse.bcf))


