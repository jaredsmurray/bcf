.ident <- function(...){
# courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
  args <- c(...)
  if( length( args ) > 2L ){
    #  recursively call ident()
    out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
  }else{
    out <- identical( args[1] , args[2] )
  }
  return( all( out ) )
}

.cp_quantile = function(x, num=10000, cat_levels=8){
  nobs = length(x)
  nuniq = length(unique(x))

  if(nuniq==1) {
    ret = x[1]
    warning("A supplied covariate contains a single distinct value.")
  } else if(nuniq < cat_levels) {
    xx = sort(unique(x))
    ret = xx[-length(xx)] + diff(xx)/2
  } else {
    q = approxfun(sort(x),quantile(x,p = 0:(nobs-1)/nobs))
    ind = seq(min(x),max(x),length.out=num)
    ret = q(ind)
  }

  return(ret)
}

#' Fit Bayesian Causal Forests
#'
#' @references Hahn, Murray, and Carvalho(2017). Bayesian regression tree models for causal inference: regularization, confounding, and heterogeneous effects.
#'  https://arxiv.org/abs/1706.09523. (Call citation("bcf") from the
#' command line for citation information in Bibtex format.)
#'
#' @details Fits the Bayesian Causal Forest model (Hahn et. al. 2018): For a response
#' variable y, binary treatment z, and covariates x,
#' \deqn{y_i = \mu(x_i, \pi_i) + \tau(x_i, \pi_i)z_i + \epsilon_i}
#' where \eqn{\pi_i} is an (optional) estimate of the propensity score \eqn{\Pr(Z_i=1 | X_i=x_i)} and
#' \eqn{\epsilon_i \sim N(0,\sigma^2)}
#'
#' Some notes:
#' \itemize{
#'    \item x_control and x_moderate must be numeric matrices. See e.g. the makeModelMatrix function in the
#'    dbarts package for appropriately constructing a design matrix from a data.frame
#'    \item sd_control and sd_moderate are the prior SD(mu(x)) and SD(tau(x)) at a given value of x (respectively). If
#'    use_muscale = FALSE, then this is the parameter \eqn{\sigma_\mu} from the original BART paper, where the leaf parameters
#'    have prior distribution \eqn{N(0, \sigma_\mu/m)}, where m is the number of trees.
#'    If use_muscale=TRUE then sd_control is the prior median of a half Cauchy prior for SD(mu(x)). If use_tauscale = TRUE,
#'    then sd_moderate is the prior median of a half Normal prior for SD(tau(x)).
#'    \item By default the prior on \eqn{\sigma^2} is calibrated as in Chipman, George and McCulloch (2008).
#'
#'
#' }
#' @param y Response variable
#' @param z Treatment variable
#' @param x_control Design matrix for the "prognostic" function mu(x)
#' @param x_moderate Design matrix for the covariate-dependent treatment effects tau(x)
#' @param pihat Length n estimates of
#' @param nburn Number of burn-in MCMC iterations
#' @param nsim Number of MCMC iterations to save after burn-in
#' @param nthin Save every nthin'th MCMC iterate. The total number of MCMC iterations will be nsim*nthin + nburn.
#' @param update_interval Print status every update_interval MCMC iterations
#' @param ntree_control Number of trees in mu(x)
#' @param sd_control SD(mu(x)) marginally at any covariate value (or its prior median if use_muscale=TRUE)
#' @param base_control Base for tree prior on mu(x) trees (see details)
#' @param power_control Power for the tree prior on mu(x) trees
#' @param ntree_moderate Number of trees in tau(x)
#' @param sd_moderate SD(tau(x)) marginally at any covariate value (or its prior median if use_tauscale=TRUE)
#' @param base_moderate Base for tree prior on tau(x) trees (see details)
#' @param power_moderate Power for the tree prior on tau(x) trees (see details)
#' @param nu Degrees of freedom in the chisq prior on \eqn{sigma^2}
#' @param lambda Scale parameter in the chisq prior on \eqn{sigma^2}
#' @param sigq Calibration quantile for the chisq prior on \eqn{sigma^2}
#' @param sighat Calibration estimate for the chisq prior on \eqn{sigma^2}
#' @param include_pi Takes values "control", "moderate", "both" or "none". Whether to
#' include pihat in mu(x) ("control"), tau(x) ("moderate"), both or none. Values of "control"
#' or "both" are HIGHLY recommended with observational data.
#' @param use_muscale Use a half-Cauchy hyperprior on the scale of mu.
#' @param use_tauscale Use a half-Normal prior on the scale of tau.
#' @return A list with elements
#' \item{tau}{\code{nsim} by \code{n} matrix of posterior samples of individual treatment effects}
#' \item{mu}{\code{nsim} by \code{n} matrix of posterior samples of individual treatment effects}
#' \item{sigma}{Length \code{nsim} vector of posterior samples of sigma}
#' @examples
#'\donttest{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#' #
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' # If you didn't know pi, you would estimate it here
#' pihat = pnorm(q)
#'
#' bcf_fit = bcf(y, z, x, x, pihat, nburn=2000, nsim=2000)
#'
#' # Get posterior of treatment effects
#' tau_post = bcf_fit$tau
#' tauhat = colMeans(tau_post)
#' plot(tau, tauhat); abline(0,1)
#'
#'}
#'\dontshow{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#' #
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' pihat = pnorm(q)
#'
#' bcf_fit = bcf(y, z, x, x, pihat, nburn=10, nsim=10)
#'
#' # Get posterior of treatment effects
#' tau_post = bcf_fit$tau
#' tauhat = colMeans(tau_post)
#' plot(tau, tauhat); abline(0,1)
#'
#'}
bcf <- function(y, z, x_control, x_moderate=x_control, pihat,
                nburn, nsim, nthin = 1, update_interval = 100,
                ntree_control = 200,
                sd_control = 2*sd(y),
                base_control = 0.95,
                power_control = 2,
                ntree_moderate = 50,
                sd_moderate = sd(y),
                base_moderate = 0.25,
                power_moderate = 3,
                nu = 3, lambda = NULL, sigq = .9, sighat = NULL,
                include_pi = "control", use_muscale=TRUE, use_tauscale=TRUE
) {

  pihat = as.matrix(pihat)
  if( !.ident(length(y),
              length(z),
              nrow(x_control),
              nrow(x_moderate),
              nrow(pihat)
  )
  ) {

    stop("Data size mismatch. The following should all be equal:
         length(y): ", length(y), "\n",
         "length(z): ", length(z), "\n",
         "nrow(x_control): ", nrow(x_control), "\n",
         "nrow(x_moderate): ", nrow(x_moderate), "\n",
         "nrow(pihat): ", nrow(pihat),"\n"
    )
  }

  if(any(is.na(y))) stop("Missing values in y")
  if(any(is.na(z))) stop("Missing values in z")
  if(any(is.na(x_control))) stop("Missing values in x_control")
  if(any(is.na(x_moderate))) stop("Missing values in x_moderate")
  if(any(is.na(pihat))) stop("Missing values in pihat")

  if(any(!is.finite(y))) stop("Non-numeric values in y")
  if(any(!is.finite(z))) stop("Non-numeric values in z")
  if(any(!is.finite(x_control))) stop("Non-numeric values in x_control")
  if(any(!is.finite(x_moderate))) stop("Non-numeric values in x_moderate")
  if(any(!is.finite(pihat))) stop("Non-numeric values in pihat")

  if(!all(sort(unique(z)) == c(0,1))) stop("z must be a vector of 0's and 1's, with at least one of each")

  if(length(unique(y))<5) warning("y appears to be discrete")

  if(nburn<0) stop("nburn must be positive")
  if(nsim<0) stop("nsim must be positive")
  if(nthin<0) stop("nthin must be positive")
  if(nthin>nsim+1) stop("nthin must be < nsim")
  if(nburn<100) warning("A low (<100) value for nburn was supplied")

  ### TODO range check on parameters

  ###
  x_c = matrix(x_control, ncol=ncol(x_control))
  x_m = matrix(x_moderate, ncol=ncol(x_moderate))
  if(include_pi=="both" | include_pi=="control") {
    x_c = cbind(x_control, pihat)
  }
  if(include_pi=="both" | include_pi=="moderate") {
    x_m = cbind(x_moderate, pihat)
  }
  cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_quantile(x_c[,i]))
  cutpoint_list_m = lapply(1:ncol(x_m), function(i) .cp_quantile(x_m[,i]))

  yscale = scale(y)
  sdy = sd(y)
  muy = mean(y)

  if(is.null(lambda)) {
    if(is.null(sighat)) {
      lmf = lm(yscale~z+as.matrix(x_c))
      sighat = summary(lmf)$sigma #sd(y) #summary(lmf)$sigma
    }
    qchi = qchisq(1.0-sigq,nu)
    lambda = (sighat*sighat*qchi)/nu
  }

  dir = tempdir()

  perm = order(z, decreasing=TRUE)

  fitbcf = bcfoverparRcppClean(yscale[perm], z[perm], t(x_c[perm,]), t(x_m[perm,,drop=FALSE]), t(x_m[1,,drop=FALSE]),
                        cutpoint_list_c, cutpoint_list_m,
                        random_des = matrix(1),
                        random_var = matrix(1),
                        random_var_ix = matrix(1),
                        random_var_df = 3,
                        nburn, nsim, nthin,
                        ntree_moderate, ntree_control, lambda, nu,
                        con_sd = ifelse(abs(2*sdy - sd_control)<1e-6, 2, sd_control/sdy),
                        mod_sd = ifelse(abs(sdy - sd_moderate)<1e-6, 1, sd_moderate/sdy)/ifelse(use_tauscale,0.674,1), # if HN make sd_moderate the prior median
                        base_moderate, power_moderate, base_control, power_control,
                        "tmp", status_interval = update_interval,
                        use_mscale = use_muscale, use_bscale = use_tauscale, b_half_normal = TRUE)

  #B = drop(fit$post_B)
  #B0 = fit$b0
  #EYs = fit$post_yhat

  #return(List::create(_["m_post"] = m_post, _["b_post"] = b_post, _["b_est_post"] = b_est_post,
  #                     _["sigma"] = sigma_post, _["msd"] = msd_post, _["bsd"] = bsd_post,
  #                     _["gamma"] = gamma_post, _["random_var_post"] = random_var_post

  m_post = muy + sdy*fitbcf$m_post[,order(perm)]
  tau_post = sdy*fitbcf$b_post[,order(perm)]
  #yhat_post = muy + sdy*fitbcf$m_post
  #yhat_post[,z==1] = yhat_post[,z==1] + fitbcf$b_post
  #yhat_post = (muy + sdy*fitbcf$m_post)
  #yhat_post[,z[perm]==1] = yhat_post[,z[perm]==1] + sdy*fitbcf$b_post
  #yhat_post = yhat_post[,order(perm)]

  list(sigma = sdy*fitbcf$sigma,
       yhat = muy + sdy*fitbcf$yhat_post[,order(perm)],
#       mu  = m_post,
       tau = tau_post,
       mu_scale = fitbcf$msd*sdy,
       tau_scale = fitbcf$bsd*sdy,
       perm = perm
  )

}
