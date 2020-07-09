#' @importFrom stats approxfun lm qchisq quantile sd
#' @importFrom RcppParallel RcppParallelLibs
Rcpp::loadModule(module = "TreeSamples", TRUE)

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

.get_chain_tree_files = function(tree_path, chain_id){
  out <- list("con_trees" = paste0(tree_path,'/',"con_trees.", chain_id, ".txt"), 
              "mod_trees" = paste0(tree_path,'/',"mod_trees.", chain_id, ".txt"))
  
  return(out) 
}

.get_do_type = function(n_cores){
  if(n_cores>1){
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    `%doType%`  <- foreach::`%dopar%`
  } else {
    cl <- NULL
    `%doType%`  <- foreach::`%do%`
  }
  
  do_type_config <- list('doType'  = `%doType%`,
                         'n_cores' = n_cores,
                         'cluster' = cl)
  
  return(do_type_config)
}

.cleanup_after_par = function(do_type_config){
  if(do_type_config$n_cores>1){
    parallel::stopCluster(do_type_config$cluster)
  }
}

#' Fit Bayesian Causal Forests
#'
#' @references Hahn, Murray, and Carvalho(2017). Bayesian regression tree models for causal inference: regularization, confounding, and heterogeneous effects.
#'  https://arxiv.org/abs/1706.09523. (Call citation("bcf") from the
#' command line for citation information in Bibtex format.)
#'
#' @details Fits the Bayesian Causal Forest model (Hahn et. al. 2018): For a response
#' variable y, binary treatment z, and covariates x, we return estimates of mu, tau, and sigma in
#' the model
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
#' @param x_control Design matrix for the prognostic function mu(x)
#' @param x_moderate Design matrix for the covariate-dependent treatment effects tau(x)
#' @param pihat Length n estimates of propensity score
#' @param w An optional vector of weights. When present, BCF fits a model \eqn{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function.
#' @param random_seed A random seed passed to R's set.seed
#' @param n_chains  An optional integer of the number of MCMC chains to run
#' @param n_cores An optional integer of the number of cores to run your MCMC chains on
#' @param n_threads An optional integer of the number of threads to parallelize within chain bcf operations on
#' @param nburn Number of burn-in MCMC iterations
#' @param nsim Number of MCMC iterations to save after burn-in. The chain will run for nsim*nthin iterations after burn-in
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
#' @param save_tree_directory Specify where trees should be saved. Keep track of this for predict(). Defaults to working directory.
#' @param nu Degrees of freedom in the chisq prior on \eqn{sigma^2}
#' @param lambda Scale parameter in the chisq prior on \eqn{sigma^2}
#' @param sigq Calibration quantile for the chisq prior on \eqn{sigma^2}
#' @param sighat Calibration estimate for the chisq prior on \eqn{sigma^2}
#' @param include_pi Takes values "control", "moderate", "both" or "none". Whether to
#' include pihat in mu(x) ("control"), tau(x) ("moderate"), both or none. Values of "control"
#' or "both" are HIGHLY recommended with observational data.
#' @param use_muscale Use a half-Cauchy hyperprior on the scale of mu.
#' @param use_tauscale Use a half-Normal prior on the scale of tau.
#' @param verbose logical, whether to print log of MCMC iterations, defaults to FALSE.
#' @return A fitted bcf object that is a list with elements
#' \item{tau}{\code{nsim} by \code{n} matrix of posterior samples of individual-level treatment effect estimates}
#' \item{mu}{\code{nsim} by \code{n} matrix of posterior samples of prognostic function E(Y|Z=0, x=x) estimates}
#' \item{sigma}{Length \code{nsim} vector of posterior samples of sigma}
#' @examples
#'\donttest{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#' 
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
#' pihat = pnorm(q)
#'
#' # nburn and nsim should be much larger, at least a few thousand each
#' # The low values below are for CRAN.
#' bcf_fit = bcf(y, z, x, x, pihat, nburn=100, nsim=10)
#'
#' # Get posterior of treatment effects
#' tau_post = bcf_fit$tau
#' tauhat = colMeans(tau_post)
#' plot(tau, tauhat); abline(0,1)
#'}
#'
#' @useDynLib bcf
#' @export
bcf <- function(y, z, x_control, x_moderate=x_control, pihat, w = NULL, 
                random_seed = sample.int(.Machine$integer.max, 1),
                n_chains = 4,
                n_cores  = n_chains,
                n_threads = max((RcppParallel::defaultNumThreads()-2)/n_cores,1), #max number of threads, minus a arbitrary holdback, over the number of cores
                nburn, nsim, nthin = 1, update_interval = 100,
                ntree_control = 250,
                sd_control = 2*sd(y),
                base_control = 0.95,
                power_control = 2,
                ntree_moderate = 50,
                sd_moderate = sd(y),
                base_moderate = 0.25,
                power_moderate = 3,
                save_tree_directory = '..',
                nu = 3, lambda = NULL, sigq = .9, sighat = NULL,
                include_pi = "control", use_muscale=TRUE, use_tauscale=TRUE, verbose=FALSE
) {

  
  if(is.null(w)){
    w <- matrix(1, ncol = 1, nrow = length(y))
    }

  pihat = as.matrix(pihat)
  if(!.ident(length(y),
             length(z),
             length(w),
             nrow(x_control),
             nrow(x_moderate),
             nrow(pihat))
    ) {
    stop("Data size mismatch. The following should all be equal:
         length(y): ", length(y), "\n",
         "length(z): ", length(z), "\n",
         "length(w): ", length(w), "\n",
         "nrow(x_control): ", nrow(x_control), "\n",
         "nrow(x_moderate): ", nrow(x_moderate), "\n",
         "nrow(pihat): ", nrow(pihat),"\n"
    )
  }

  if(any(is.na(y))) stop("Missing values in y")
  if(any(is.na(z))) stop("Missing values in z")
  if(any(is.na(w))) stop("Missing values in w")
  if(any(is.na(x_control))) stop("Missing values in x_control")
  if(any(is.na(x_moderate))) stop("Missing values in x_moderate")
  if(any(is.na(pihat))) stop("Missing values in pihat")
  if(any(!is.finite(y))) stop("Non-numeric values in y")
  if(any(!is.finite(z))) stop("Non-numeric values in z")
  if(any(!is.finite(w))) stop("Non-numeric values in w")
  if(any(!is.finite(x_control))) stop("Non-numeric values in x_control")
  if(any(!is.finite(x_moderate))) stop("Non-numeric values in x_moderate")
  if(any(!is.finite(pihat))) stop("Non-numeric values in pihat")
  if(!all(sort(unique(z)) == c(0,1))) stop("z must be a vector of 0's and 1's, with at least one of each")

  if(length(unique(y))<5) warning("y appears to be discrete")

  if(nburn<0) stop("nburn must be positive")
  if(nsim<0) stop("nsim must be positive")
  if(nthin<0) stop("nthin must be positive")
  if(nthin>nsim+1) stop("nthin must be < nsim")
  if(nburn<1000) warning("A low (<1000) value for nburn was supplied")
  if(nsim*nburn<1000) warning("A low (<1000) value for total iterations after burn-in was supplied")

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

  sdy = sqrt(Hmisc::wtd.var(y, w))
  muy = stats::weighted.mean(y, w)
  yscale = (y-muy)/sdy


  if(is.null(lambda)) {
    if(is.null(sighat)) {
      lmf = lm(yscale~z+as.matrix(x_c), weights = w)
      sighat = summary(lmf)$sigma #sd(y) #summary(lmf)$sigma
    }
    qchi = qchisq(1.0-sigq,nu)
    lambda = (sighat*sighat*qchi)/nu
  }

  dir = tempdir()

  perm = order(z, decreasing=TRUE)

  con_sd = ifelse(abs(2*sdy - sd_control)<1e-6, 2, sd_control/sdy)
  mod_sd = ifelse(abs(sdy - sd_moderate)<1e-6, 1, sd_moderate/sdy)/ifelse(use_tauscale,0.674,1) # if HN make sd_moderate the prior median

  RcppParallel::setThreadOptions(numThreads=n_threads)
  
  do_type_config <- .get_do_type(n_cores)
  `%doType%` <- do_type_config$doType
  
  chain_out <- foreach::foreach(iChain=1:n_chains) %doType% {
    
    this_seed = random_seed + iChain - 1
    
    cat("Calling bcfoverparRcppClean From R\n")
    set.seed(this_seed)
    
    tree_files = .get_chain_tree_files(save_tree_directory, iChain)
    
    print(tree_files)

    fitbcf = bcfoverparRcppClean(y_ = yscale[perm], z_ = z[perm], w_ = w[perm],
                                 x_con_ = t(x_c[perm,,drop=FALSE]), x_mod_ = t(x_m[perm,,drop=FALSE]), 
                                 x_con_info_list = cutpoint_list_c, 
                                 x_mod_info_list = cutpoint_list_m,
                                 random_des = matrix(1),
                                 random_var = matrix(1),
                                 random_var_ix = matrix(1),
                                 random_var_df = 3,
                                 burn = nburn, nd = nsim, thin = nthin,
                                 ntree_mod = ntree_moderate, ntree_con = ntree_control, 
                                 lambda = lambda, nu = nu,
                                 con_sd = con_sd,
                                 mod_sd = mod_sd, # if HN make sd_moderate the prior median
                                 mod_alpha = base_moderate, 
                                 mod_beta = power_moderate, 
                                 con_alpha = base_control, 
                                 con_beta = power_control,
                                 treef_con_name_ = tree_files$con_trees, 
                                 treef_mod_name_ = tree_files$mod_trees, 
                                 status_interval = update_interval,
                                 use_mscale = use_muscale, use_bscale = use_tauscale, 
                                 b_half_normal = TRUE, verbose_sigma=verbose)
    
    cat("bcfoverparRcppClean returned to R\n")

    ac = fitbcf$m_post[,order(perm)]

    Tm = fitbcf$b_post[,order(perm)] * (1.0/ (fitbcf$b1 - fitbcf$b0))

    Tc = ac * (1.0/fitbcf$msd) 

    tau_post = sdy*fitbcf$b_post[,order(perm)]

    mu_post  = muy + sdy*(Tc*fitbcf$msd + Tm*fitbcf$b0)
    
    list(sigma = sdy*fitbcf$sigma,
         yhat = muy + sdy*fitbcf$yhat_post[,order(perm)],
         sdy = sdy,
         con_sd = con_sd,
         mod_sd = mod_sd,
         muy = muy,
         mu  = mu_post,
         tau = tau_post,
         mu_scale = fitbcf$msd,
         tau_scale = fitbcf$bsd,
         b0 = fitbcf$b0,
         b1 = fitbcf$b1,
         perm = perm,
         include_pi = include_pi,
         random_seed=this_seed
    )
    
  }


  all_sigma = c()
  all_mu_scale = c()
  all_tau_scale = c()

  all_b0 = c()
  all_b1 = c()
  
  all_yhat = c()
  all_mu   = c()
  all_tau  = c()
  
  chain_list=list()

  n_iter = length(chain_out[[1]]$sigma)
  
  
  for (iChain in 1:n_chains){
    sigma        <- chain_out[[iChain]]$sigma
    mu_scale     <- chain_out[[iChain]]$mu_scale
    tau_scale    <- chain_out[[iChain]]$tau_scale
    
    b0          <- chain_out[[iChain]]$b0
    b1          <- chain_out[[iChain]]$b1

    yhat         <- chain_out[[iChain]]$yhat
    tau          <- chain_out[[iChain]]$tau
    mu           <- chain_out[[iChain]]$mu

    # -----------------------------    
    # Support Old Output
    # -----------------------------
    all_sigma       = c(all_sigma,     sigma)
    all_mu_scale    = c(all_mu_scale,  mu_scale)
    all_tau_scale   = c(all_tau_scale, tau_scale)
    all_b0 = c(all_b0, b0)
    all_b1 = c(all_b1, b1)

    all_yhat = rbind(all_yhat, yhat)
    all_mu   = rbind(all_mu,   mu)
    all_tau  = rbind(all_tau,  tau)

    # -----------------------------    
    # Make the MCMC Object
    # -----------------------------

    scalar_df <- data.frame("sigma"     = sigma,
                            "tau_bar"   = matrixStats::rowWeightedMeans(tau, w),
                            "mu_bar"    = matrixStats::rowWeightedMeans(mu, w),
                            "yhat_bar"  = matrixStats::rowWeightedMeans(yhat, w),
                            "mu_scale"  = mu_scale, 
                            # "tau_scale" = tau_scale,
                            "b0"  = b0, 
                            "b1"  = b1)
    
    # y_df <- as.data.frame(chain$yhat)
    # colnames(y_df) <- paste0('y',1:ncol(y_df))
    # 
    # mu_df <- as.data.frame(chain$mu)
    # colnames(mu_df) <- paste0('mu',1:ncol(mu_df))
    # 
    # tau_df <- as.data.frame(chain$tau)
    # colnames(tau_df) <- paste0('tau',1:ncol(tau_df))
    
    chain_list[[iChain]] <- coda::as.mcmc(scalar_df)
    # -----------------------------    
    # Sanity Check Constants Accross Chains
    # -----------------------------
    if(chain_out[[iChain]]$sdy        != chain_out[[1]]$sdy)        stop("sdy not consistent between chains for no reason")
    if(chain_out[[iChain]]$con_sd     != chain_out[[1]]$con_sd)     stop("con_sd not consistent between chains for no reason")
    if(chain_out[[iChain]]$mod_sd     != chain_out[[1]]$mod_sd)     stop("mod_sd not consistent between chains for no reason")
    if(chain_out[[iChain]]$muy        != chain_out[[1]]$muy)        stop("muy not consistent between chains for no reason")
    if(chain_out[[iChain]]$include_pi != chain_out[[1]]$include_pi) stop("include_pi not consistent between chains for no reason")
    if(any(chain_out[[iChain]]$perm   != chain_out[[1]]$perm))      stop("perm not consistent between chains for no reason")
  }

  fitObj <- list(sigma = all_sigma,
                 yhat = all_yhat,
                 sdy = chain_out[[1]]$sdy,
                 muy = chain_out[[1]]$muy,
                 mu  = all_mu,
                 tau = all_tau,
                 mu_scale = all_mu_scale,
                 tau_scale = all_tau_scale,
                 b0 = all_b0,
                 b1 = all_b1,
                 perm = perm,
                 include_pi = chain_out[[1]]$include_pi,
                 random_seed = chain_out[[1]]$random_seed,
                 coda_chains = coda::as.mcmc.list(chain_list),
                 raw_chains = chain_out)
  
  attr(fitObj, "class") <- "bcf"
  
  .cleanup_after_par(do_type_config)
  
  return(fitObj)
}

#' Takes a fitted bcf object produced by bcf() and produces summary stats and MCMC diagnostics.
#' This function is built using the coda package and meant to mimic output from rstan::print.stanfit().
#' It includes, for key parameters, posterior summary stats, effective sample sizes, 
#' and Gelman and Rubin's convergence diagnostics. 
#' By default, those parameters are: sigma (the error standard deviation when the weights
#' are all equal), tau_bar (the estimated sample average treatment effect), mu_bar
#' (the average outcome under control/z=0 across all observations in the sample), and
#' yhat_bat (the average outcome under the realized treatment assignment across all
#' observations in the sample).
#' 
#' We strongly suggest updating the coda package to our 
#' Github version, which uses the Stan effective size computation. 
#' We found the native coda effective size computation to be overly optimistic in some situations
#' and are in discussions with the coda package authors to change it on CRAN.
#' @param object output from a BCF predict run.
#' @param ... additional arguments affecting the summary produced.
#' @param params_2_summarise parameters to summarise.
#' @examples
#'\donttest{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#' 
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
#' # Get model fit diagnostics
#' summary(bcf_fit)
#'
#'}
#' @export
summary.bcf <- function(object,
                        ..., 
                        params_2_summarise = c('sigma','tau_bar','mu_bar','yhat_bar')){

  chains_2_summarise <- object$coda_chains[,params_2_summarise]

  message("Summary statistics for each Markov Chain Monte Carlo run")
  print(summary(chains_2_summarise))

  cat("\n----\n\n")


  message("Effective sample size for summary parameters")
  
  ef = function(e) {
    if(e$message == "unused argument (crosschain = TRUE)") {
      cat("Reverting to coda's default ESS calculation. See ?summary.bcf for details.\n\n")
      print(coda::effectiveSize(chains_2_summarise))
    } else {
      stop(e)
    }
  }
  tryCatch(print(coda::effectiveSize(chains_2_summarise, crosschain = TRUE)),
           error = ef) 
  cat("\n----\n\n")
  
  
  if (length(chains_2_summarise) > 1){
    message("Gelman and Rubin's convergence diagnostic for summary parameters")
    print(coda::gelman.diag(chains_2_summarise, autoburnin = FALSE))
    cat("\n----\n\n")
    
  }
  
}


#' Takes a fitted bcf object produced by bcf() and produces predictions for a new set of covariate values
#' 
#' This function takes in an existing BCF model fit and uses it to predict estimates for new data.
#' It is important to note that this function requires that you indicate where the trees from the model fit are saved.
#' You can do so using the save_tree_directory argument in bcf(). Otherwise, they will be saved in the working directory.
#' The bcf() function automatically saves those in the same directory as the 
#' @param object output from a BCF predict run
#' @param ... additional arguments affecting the predictions produced.
#' @param x_predict_control matrix of covariates for the "prognostic" function mu(x) for predictions (optional)
#' @param x_predict_moderate matrix of covariates for the covariate-dependent treatment effects tau(x) for predictions (optional)
#' @param z_pred Treatment variable for predictions (optional except if x_pre is not empty)
#' @param pi_pred propensity score for prediction
#' @param save_tree_directory directory where the trees have been saved
#' @param n_cores An optional integer of the number of cores to run your MCMC chains on
#' @examples
#'\donttest{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
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
#' bcf_fit = bcf(y               = y,
#'               z               = z,
#'               x_control       = x,
#'               x_moderate      = x,
#'               pihat           = pihat,
#'               nburn           = n_burn,
#'               nsim            = n_sim,
#'               n_chains        = 2,
#'               update_interval = 1,
#'               save_tree_directory = './trees')
#'
#' # Predict using new data
#' 
#' x_pred = matrix(rnorm(n*p), nrow=n)
#' 
#' pred_out = predict(bcf_out=bcf_fit,
#'                    x_predict_control=x_pred,
#'                    x_predict_moderate=x_pred,
#'                    pi_pred=pihat,
#'                    z_pred=z,
#'                    save_tree_directory = './trees')
#'
#'}
#' @export
predict.bcf <- function(object, 
                        x_predict_control,
                        x_predict_moderate,
                        pi_pred,
                        z_pred, 
                        save_tree_directory,
                        n_cores=2,
                        ...) {
                        
    if(any(is.na(x_predict_moderate))) stop("Missing values in x_predict_moderate")
    if(any(is.na(x_predict_control))) stop("Missing values in x_predict_control")
    if(any(is.na(z_pred))) stop("Missing values in z_pred")
    if(any(!is.finite(x_predict_moderate))) stop("Non-numeric values in x_pred_moderate")
    if(any(!is.finite(x_predict_control))) stop("Non-numeric values in x_pred_control")
    if(any(!is.finite(pi_pred))) stop("Non-numeric values in pi_pred")
    if(!all(sort(unique(z_pred)) == c(0,1))) stop("z_pred must be a vector of 0's and 1's, with at least one of each")

    if((is.null(x_predict_moderate) & !is.null(x_predict_control)) | (!is.null(x_predict_moderate) & is.null(x_predict_control))) {
        stop("If you want to predict, you need to add values to both x_pred_control and x_pred_moderate")
    }

    pi_pred = as.matrix(pi_pred)
    if(!.ident(length(z_pred),
                nrow(x_predict_moderate),
                nrow(x_predict_control),
                nrow(pi_pred))
        ) {
        stop("Data size mismatch. The following should all be equal:
            length(z_pred): ", length(z_pred), "\n",
            "nrow(x_pred_moderate): ", nrow(x_predict_moderate), "\n",
            "nrow(x_pred_control): ", nrow(x_predict_control), "\n",
            "nrow(pi_pred): ", nrow(pi_pred), "\n"
        )
    }


    cat("Initializing BCF Prediction\n")
    x_pm = matrix(x_predict_moderate, ncol=ncol(x_predict_moderate))
    x_pc = matrix(x_predict_control, ncol=ncol(x_predict_control))

    if(object$include_pi=="both" | object$include_pi=="control") {
        x_pc = cbind(x_predict_control, pi_pred)
    }
    if(object$include_pi=="both" | object$include_pi=="moderate") {
        x_pm = cbind(x_predict_moderate, pi_pred)
    }


    cat("Starting Prediction \n")

    n_chains = length(object$coda_chains)
    
    do_type_config <- .get_do_type(n_cores)
    `%doType%` <- do_type_config$doType
    
    chain_out <- foreach::foreach(iChain=1:n_chains) %doType% {
      
      tree_files = .get_chain_tree_files(save_tree_directory, iChain)

      cat("Starting to Predict Chain ", iChain, "\n")
      
      mods = TreeSamples$new()
      mods$load(tree_files$mod_trees)
      Tm = mods$predict(t(x_pm))
      
      cons = TreeSamples$new()
      cons$load(tree_files$con_trees)
      Tc = cons$predict(t(x_pc))
      
      
      list(Tm = Tm,
           Tc = Tc)
    }
    
    all_yhat = c()
    all_mu   = c()
    all_tau  = c()
    
    chain_list=list()

    muy = object$muy
      
    sdy = object$sdy
    
    for (iChain in 1:n_chains){
      
      
        # Extract Chain Specific Information
    
        Tm = chain_out[[iChain]]$Tm
        Tc = chain_out[[iChain]]$Tc
        
        this_chain_bcf_out = object$raw_chains[[iChain]]
        
        b1 = this_chain_bcf_out$b1
        b0 = this_chain_bcf_out$b0
        mu_scale = this_chain_bcf_out$mu_scale
        


        # Calculate, tau, y, and mu

        
        mu  = muy + sdy*(Tc*mu_scale + Tm*b0)
        tau = sdy*(b1 - b0)*Tm
        yhat = mu + t(t(tau)*z_pred)
        
        
        # Package Output up
        all_yhat = rbind(all_yhat, yhat)
        all_mu   = rbind(all_mu,   mu)
        all_tau  = rbind(all_tau,  tau)
        
        
        
        scalar_df <- data.frame("tau_bar"   = matrixStats::rowWeightedMeans(tau, w=NULL),
                                "mu_bar"    = matrixStats::rowWeightedMeans(mu, w=NULL),
                                "yhat_bar"  = matrixStats::rowWeightedMeans(yhat, w=NULL))

        chain_list[[iChain]] <- coda::as.mcmc(scalar_df)
    }

   .cleanup_after_par(do_type_config)


    list(tau = all_tau,
         mu = all_mu,
         yhat = all_yhat,
         coda_chains = coda::as.mcmc.list(chain_list))
}
