############################################################################################
# Cox Proportional Hazards Models - Comparison Models for BMA Approach
############################################################################################


### Function to calculate the log posterior density of regression effects using the
###   first Cox proportional hazards model (global treatment effect only)
# theta: scalar global treatment effect
# mu0: prior mean scalar for theta
# Sig0: prior variance for theta
# eta0: gamma shape hyperparameter for baseline hazard
# phi0: gamma rate hyperparameter for baseline hazard
# y.mat: S x 2 matrix of event times for all region-by-treatment combinations
# nu.mat: S x 2 matrix of statuses (1 if failed, 0 if censored) for all region-by-trt combinations
# z.mat: S x 2 matrix of trt indicators (1 if trt, 0 otherwise) for all region-by-trt combinations
# eta.star: piece of posterior density (scalar)
log.ph.dens1 <- function( theta, mu0, Sig0, eta0, phi0, y.mat, nu.mat, z.mat, eta.star ){
  
  # Calculate piece of posterior density (not on log scale)
  phi.star <- phi0 + sum( y.mat * exp(z.mat * theta) )
  
  # Log of full conditional of theta_l
  log.val <- -eta.star * log(phi.star) + sum( z.mat * nu.mat * theta ) -
    .5 * t(theta - mu0) %*% solve(Sig0) %*% (theta - mu0)
  
  return(log.val)
  
}



### Function to calculate the log posterior density of regression effects using the
###   second Cox proportional hazards model (region-specific treatment effects only)
# theta: S x 1 vector with region-specific treatment effects
# mu0: prior mean vector for theta
# Sig0: prior covariance matrix for theta
# eta.vec: S x 1 vector of gamma shape hyperparameters for region-specific baseline hazards
# phi.vec: S x 1 vector of gamma rate hyperparameters for region-specific baseline hazards
# y.mat: S x 2 matrix of event times for all region-by-treatment combinations
# nu.mat: S x 2 matrix of statuses (1 if failed, 0 if censored) for all region-by-trt combinations
# z.mat: S x 2 matrix of trt indicators (1 if trt, 0 otherwise) for all region-by-trt combinations
# eta.star.vec: piece of posterior density (S x 1 vector)
log.ph.dens2 <- function( theta, mu0, Sig0, eta.vec, phi.vec, y.mat, nu.mat, z.mat,
                          eta.star.vec){
  
  # Calculate piece of posterior density (not on log scale)
  phi.star.vec <- phi.vec + rowSums(y.mat * exp(z.mat * theta))
  
  # Log of full conditional of theta_l
  log.val <- sum( -eta.star.vec * log(phi.star.vec) ) + sum( z.mat * nu.mat * theta ) -
    .5 * t(theta - mu0) %*% solve(Sig0) %*% (theta - mu0)
  
  return(log.val)
  
}



### Function that fits Cox proportional hazards models
###   PH1 model: model used to estimate global treatment effect
###   PH2 model: model used to estimate region-specific treatment effects
# Returns posterior summaries and acceptance ratios for random walk Metropolis-Hastings
# dat: condensed dataset where each region has two observations (one for each treatment group)
# mu0.ph1: prior mean for PH1 model (scalar)
# mu0.ph2: prior mean vector (S x 1) for PH2 model
# Sig0.ph1: prior variance for PH1 model (scalar)
# Sig0.ph2: prior covariance matrix (S x S) for PH2 model
# eta0.ph1: shape hyperparameter for prior on baseline hazard for PH1 model (scalar)
# eta0.ph2: shape hyperparameter vector (S x 1) for priors on baseline hazards for PH2 model
# phi0.ph1: rate hyperparameter for prior on baseline hazard for PH1 model (scalar)
# phi0.ph2: rate hyperparameter vector (S x 1) for priors on baseline hazards for PH2 model
# sd.MH1: tuning parameter in MH algorithm for PH1 model (stand. dev. of proposal distribution)
# sd.MH2: tuning parameter in MH algorithm for PH2 model (stand. dev. of proposal distribution)
# n.draws: total number of samples (before burn in and thinning)
# burn.in: burn in to remove from samples beginning of samples
# thin: amount by which to thin
# gamma0: value of interest in null hypothesis to be tested
ph.models <- function( dat, mu0.ph1, mu0.ph2, Sig0.ph1, Sig0.ph2, eta0.ph1, eta0.ph2,
                       phi0.ph1, phi0.ph2, sd.MH1, sd.MH2, n.draws, burn.in, thin, gamma0 ){
  
  
  ### Rearrange data into S x 2 matrices with rows corresponding to region and columns
  ###   corresponding to treatment group (first column: control; second column: treatment)
  y.mat <- matrix(dat$y.iq, nrow = 4, ncol = 2, byrow = TRUE)
  nu.mat <- matrix(dat$nu.iq, nrow = 4, ncol = 2, byrow = TRUE)
  z.mat <- matrix(dat$trt.group, nrow = 4, ncol = 2, byrow = TRUE)
  
  
  ### Calculate unchanging pieces of posterior densities of treatment effects (not on log scale)
  eta.star <- sum(nu.mat) + eta0.ph1          # used for global treatment effect
  eta.star.vec <- eta0.ph2 + rowSums(nu.mat)  # used for region-specific treatment effects
  
  
  ### Save prior hyperparameters as matrices
  mu0.ph1 <- matrix(mu0.ph1, ncol = 1)
  mu0.ph2 <- matrix(mu0.ph2, ncol = 1)
  Sig0.ph1 <- as.matrix(Sig0.ph1)
  Sig0.ph2 <- as.matrix(Sig0.ph2)
  
  
  ### Vectors/matrices to store sampled treatment effects and MH acceptance indicators
  gamma.ph1 <- numeric(n.draws)                      # store global treatment effects
  gamma.ph2 <- matrix(0, nrow = n.draws, ncol = S)   # store region-specific treatment effects
  acc.glob <- numeric(n.draws - 1)      # store acceptance indicators for global effect (PH1)
  acc.reg <- numeric(n.draws - 1)       # store acceptance indicators for region effects (PH2)
  
  
  ### Algorithm
  
  # Sample from posterior distributions of treatment effects (global and region-specific)
  # using random walk Metropolis-Hastings algorithm
  for(m in 2:n.draws){
    
    # Sample global treatment effect (for PH1 model)
    gamma.glob.prop <- gamma.ph1[m-1] + rnorm(1, 0, sd.MH1)    # new proposal for global effect
    log.ratio1 <- log.ph.dens1( theta = gamma.glob.prop, mu0 = mu0.ph1, Sig0 = Sig0.ph1,
                                eta0 = eta0.ph1, phi0 = phi0.ph1, y.mat = y.mat,
                                nu.mat = nu.mat, z.mat = z.mat, eta.star = eta.star ) -
      log.ph.dens1( theta = gamma.ph1[m-1], mu0 = mu0.ph1, Sig0 = Sig0.ph1, eta0 = eta0.ph1,
                    phi0 = phi0.ph1, y.mat = y.mat, nu.mat = nu.mat, z.mat = z.mat,
                    eta.star = eta.star )
    
    # Accept or reject global treatment effect proposal
    if( runif(1) < exp(log.ratio1) ){
      gamma.ph1[m] <- gamma.glob.prop
      acc.glob[m-1] <- 1
    } else{
      gamma.ph1[m] <- gamma.ph1[m-1] 
    }
    
    
    # Sample region-specific treatment effects (multivariate random walk MH) (for PH2 model)
    gamma.reg.props <- gamma.ph2[m-1,] + rnorm(S, 0, sd.MH2[1])    # new proposals for all regions
    log.ratio2 <- log.ph.dens2( theta = gamma.reg.props, mu0 = mu0.ph2, Sig0 = Sig0.ph2,
                                eta.vec = eta0.ph2, phi.vec = phi0.ph2, y.mat = y.mat,
                                nu.mat = nu.mat, z.mat = z.mat, eta.star.vec = eta.star.vec ) -
      log.ph.dens2( theta = gamma.ph2[m-1,], mu0 = mu0.ph2, Sig0 = Sig0.ph2,
                    eta.vec = eta0.ph2, phi.vec = phi0.ph2, y.mat = y.mat,
                    nu.mat = nu.mat, z.mat = z.mat, eta.star.vec = eta.star.vec )
    
    # Accept or reject region-specific treatment effects proposal
    if( runif(1) < exp(log.ratio2) ){
      gamma.ph2[m,] <- gamma.reg.props
      acc.reg[m-1] <- 1
    } else{
      gamma.ph2[m,] <- gamma.ph2[m-1,]
    }
    
  }
  
  # Acceptance rates
  acc.rate.ph1 <- mean(acc.glob)    # acceptance rate for global treatment effect
  acc.rate.ph2 <- mean(acc.reg)     # acceptance rate for region-specific treatment effects
  
  
  ### Posterior summaries
  
  # Subset posterior samples to account for burn-in and thinning
  samps.keep <- seq( burn.in + 1, n.draws, by = thin )
  gamma.samp.ph1 <- gamma.ph1[ samps.keep ]
  gamma.samp.ph2 <- gamma.ph2[ samps.keep, ]
  
  
  # Posterior means
  mean.ph1 <- mean(gamma.samp.ph1)        # mean of global treatment effect
  mean.ph2 <- colMeans(gamma.samp.ph2)    # mean of region-specific treatment effects
  
  
  # Posterior standard deviations
  sd.ph1 <- sd(gamma.samp.ph1)            # stand. dev. of global treatment effect
  sd.ph2 <- apply(gamma.samp.ph2, 2, sd)  # stand. dev. of region-specific treatment effects
  
  
  # Posterior probabilities that treatment effects are less than gamma0
  less.gamma0.ph1 <- mean(gamma.samp.ph1 < gamma0)       # Pr(gamma_G < gamma0|D)
  less.gamma0.ph2 <- colMeans(gamma.samp.ph2 < gamma0)   # Pr(gamma_i < gamma0|D), i=1,...,S
  
  
  # Save posterior summaries and MH acceptance rates in list
  ph.mods.summaries <- list()
  ph.mods.summaries$ph1.mean <- mean.ph1
  ph.mods.summaries$ph2.mean <- mean.ph2
  ph.mods.summaries$ph1.sd <- sd.ph1
  ph.mods.summaries$ph2.sd <- sd.ph2
  ph.mods.summaries$ph1.less.gamma0 <- less.gamma0.ph1
  ph.mods.summaries$ph2.less.gamma0 <- less.gamma0.ph2
  ph.mods.summaries$acceptance.rate.ph1 <- acc.rate.ph1
  ph.mods.summaries$acceptance.rate.ph2 <- acc.rate.ph2
  
  
  # Return list
  return(ph.mods.summaries)
  
}
