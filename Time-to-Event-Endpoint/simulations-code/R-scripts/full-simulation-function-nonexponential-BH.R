############################################################################################
# FULL SIMULATION FUNCTION
#
# Compare BMA approach to two Cox proportional hazards models (CPHM)--one with a common
# treatment effect and one with a region-by-treatment interaction (no covariates)--and
# a Bayesian hierarchical model (BHM) used by the FDA
############################################################################################


### Source R scripts
source("../R-scripts/bma-asymptotic-function.R")
source("../R-scripts/proportional-hazards-models.R")


### Simulation function
# Return the proportion of simulations that each method rejects the null hypothesis
# S: number of regions
# T0: max number of distinct regions to consider
# n.i: S x 1 vector of regional sample size
# reg.names: S x 1 vector of region names in alphabetical order
# K: number of time intervals in piecewise exponential model
# max.time: maximum follow up time (from beginning of study)
# max.enroll.time: maximum enrollment time
# dropout.rate: dropout rate
# pw.times: vector of left bounds of intervals, each with its own underlying baseline hazard
# pw.haz: vector of underlying baseline hazards corresponding to each defined interval
# reg.te: S x 1 vector of region-specific treatment effects
# num.covs: number of covariates (p)
# cov.effects: p x 1 vector of true covariate effects
# cov.names: p x 1 vector of labels for covariate names
# cov.props: p x 1 vector of success probabilities for binary covariates (0 if no binary covs)
# cov.means: p x 1 vector of means for normal covariates (0 if no continuous covs)
# cov.sds: p x 1 vector of standard deviations for normal covariates (0 if no continuous covs)
# mod.priors: L x 1 vector of prior model probabilities, where L is number of models
# mu0: mean vector for prior on treatment and covariate effects
# Sig0: covariance matrix for prior on treatment and covariate effects
# eta0: S x K matrix of shape hyperparameters for prior on baseline hazards
# phi0: S x K matrix of rate hyperparameters for prior on baseline hazards
# gamma0: value in hypotheses for posterior probability Pr(gamma > gamma0|D)
# n.draws: number of posterior samples to draw when evaluating consistency
# epsilon.star: vector of possible minimal clinically important regional differences of trtmnt
#      effects - used for consistency measures (used with BMA)
# beta.star: probability cutoff for global inconsistency for which two regions are considered
#      to be clinically different (used with BMA)
# pi0: value for which we want to calculate probabilities with local consistency (used with BMA)
# nsims: number of datasets to simulate
# sd.CPHM1: tuning parameter used in Metropolis-Hastings algorithm for first CPHM
# sd.CPHM2: tuning parameter used in Metropolis-Hastings algorithm for second CPHM
# print.iters: whether or not to print iteration number in R consolue
run.sims <- function( S, T0, n.i, reg.names, K, max.time, max.enroll.time, dropout.rate, pw.times,
                      pw.haz, reg.te, num.covs, cov.effects, cov.names, cov.props, cov.means,
                      cov.sds, mod.priors, mu0, Sig0, eta0, phi0, gamma0, n.draws, epsilon.star,
                      beta.star, pi0, nsims, sd.CPHM1, sd.CPHM2, print.iters = FALSE ){
  
  ### Matrices to hold simulation results
  
  # BMA: store PMPs for all models
  num.mods <- num_models(S, T0)                                # number of models for BMA
  BMA.pmp.mat <- matrix( 0, nrow = nsims, ncol = num.mods)     # posterior model probabilities
  
  # BMA: store means, sd's, and probabilities for global and region-specific treatment effects
  BMA.gte.mean <- matrix( 0, nrow = nsims, ncol = 1)           # global treatment effect mean
  BMA.gte.sd <- matrix( 0, nrow = nsims, ncol = 1)             # global treatment effect sd
  BMA.gte.less.gamma0 <- matrix( 0, nrow = nsims, ncol = 1)    # Pr(gamma_G < gamma0|D)
  BMA.rte.mean <- matrix( 0, nrow = nsims, ncol = S)           # regional treatment effect means
  BMA.rte.sd <- matrix( 0, nrow = nsims, ncol = S)             # regional treatment effect sd's
  BMA.rte.less.gamma0 <- matrix( 0, nrow = nsims, ncol = S)    # Pr(gamma_i < gamma0|D), i=1,...,S
  
  # BMA: store means and sd's for covariate effects
  BMA.cov.mean <- matrix( 0, nrow = nsims, ncol = num.covs)    # covariate effect means
  BMA.cov.sd <- matrix( 0, nrow = nsims, ncol = num.covs)      # covariate effect sd's
  
  # BMA: store local consistency probabilities (PMDA and leave-one-out)
  num.eps <- length(epsilon.star)
  loc.consis.PMDA.TE <- matrix(0, nrow = nsims, ncol = S)
  loc.consis.PMDA.RR <- matrix(0, nrow = nsims, ncol = S)
  loc.consis.loo <- list()    # Store results in separate matrix for each epsilon.star value
  for(q in 1:num.eps){
    loc.consis.loo[[q]] <- matrix(0, nrow = nsims, ncol = S)
  }
  
  # BMA: store pairwise consistency/inconsistency probabilities and
  # epsilon.star-level global inconsistency probabilities
  prws.consis.combined <- matrix( 0, nrow = num.eps * S + num.eps, ncol = S )
  prws.inconsis.combined <- matrix( 0, nrow = num.eps * S + num.eps, ncol = S )
  glob.consis.prob <- matrix( 0, nrow = nsims, ncol = num.eps)
  
  # CPHM: store means, sd's, and probabilities for global and region-specific treatment effects
  CPHM.gte.mean <- matrix( 0, nrow = nsims, ncol = 1)          # global trt effect mean
  CPHM.gte.sd <- matrix( 0, nrow = nsims, ncol = 1)            # global trt effect sd
  CPHM.gte.less.gamma0 <- matrix( 0, nrow = nsims, ncol = 1)   # Pr(gamma_G < gamma0|D)
  CPHM.rte.mean <- matrix( 0, nrow = nsims, ncol = S)          # region-spec. trt effect mean
  CPHM.rte.sd <- matrix( 0, nrow = nsims, ncol = S)            # region-spec. trt effect sd
  CPHM.rte.less.gamma0 <- matrix( 0, nrow = nsims, ncol = S)   # Pr(gamma_i < gamma0|D), i=1,...,S
  
  # BHM: store means, sd's, and probabilities for global and region-specific treatment effects
  BHM.gte.mean <- matrix( 0, nrow = nsims, ncol = 1)           # global trt effect mean
  BHM.gte.sd <- matrix( 0, nrow = nsims, ncol = 1)             # global trt effect sd
  BHM.gte.less.gamma0 <- matrix( 0, nrow = nsims, ncol = 1)    # Pr(gamma_G < gamma0|D)
  BHM.rte.mean <- matrix( 0, nrow = nsims, ncol = S)           # region-spec. trt effect mean
  BHM.rte.sd <- matrix( 0, nrow = nsims, ncol = S)             # region-spec. trt effect sd
  BHM.rte.less.gamma0 <- matrix( 0, nrow = nsims, ncol = S)    # Pr(gamma_i < gamma0|D), i=1,...,S
  
  # Matrix to store computation times (in minutes) for each model
  sim.time <- matrix(0, nrow = nsims, ncol = 3)
  colnames(sim.time) <- c("BMA", "CPHMs", "BHM")
  
  
  ### Prior hyperparameters and algorithm parameters for CPHMs
  
  # Compile matrix of model region classifications (for BMA)
  mod.mat <- compile_modMat(S, T0)
  
  # Define piecewise hazard function used to generate each region's data
  haz <- function(t, x, betas, lb.interval, haz.interval){
    #haz.interval[findInterval(t, lb.interval)] * exp( cbind(x) %*% cbind(betas) )
    haz.interval[findInterval(t, lb.interval)] * exp( as.matrix(x) %*% cbind(betas) )
  }
  
  # Prior elicitation
  mu0.CPHM1 <- 0                  # prior mean vector for CPHM1
  mu0.CPHM2 <- numeric(S)         # prior mean vector for CPHM2
  Sig0.CPHM1 <- as.matrix(10000)  # prior covariance matrix for CPHM1
  Sig0.CPHM2 <- 10000 * diag(S)   # prior covariance matrix for CPHM2
  eta0.CPHM1 <- .01               # shape hyperparameter for prior on baseline hazard for CPHM1
  eta0.CPHM2 <- rep(.01, S)       # shape hyperparameter for prior on baseline hazard for CPHM2
  phi0.CPHM1 <- .01               # rate hyperparameters for prior on baseline hazard for CPHM1
  phi0.CPHM2 <- rep(.01, S)       # rate hyperparameters for prior on baseline hazard for CPHM2
  
  # Input/tuning parameters for random walk Metropolis Hastings algorithms
  n.draws.CPHM <- 55000     # total number of samples (before burn in and thinning)
  burn.in.CPHM <- 5000      # burn-in to remove from samples beginning of samples
  thin.CPHM <- 5            # amount by which to thin
  
  
  ### Algorithm parameters for BHM
  
  # List of parameters to sample
  parms.BHM <- c("mu", "mui", "tau2")      # parameters to sample
  
  # Input/tuning parameters for JAGS model
  n.draws.BHM <- 100000     # total number of samples including burn-in
  burn.in.BHM <- 1000       # burn-in to remove from samples beginning of samples
  thin.BHM <- 5             # amount by which to thin
  
  
  ### Run simulation
  for(j in 1:nsims){
    
    ## Generate dataset
    
    # Simulate design matrix
    sub.id.index <- cbind( c(1, cumsum(n.i[-S]) + 1),  # first column is first subject ID in region i
                           c(cumsum(n.i[-S]), N) )     # last column is last subject ID in region i
    id.vec <- 1:N                    # vector of subject IDs
    trt.vec <- NULL                  # vector of treatment assignments
    reg.vec <- rep(1:S, n.i)         # vector of region labels
    for(i in 1:S){
      trt.vec <- c( trt.vec, rep(c(0,1), c(ceiling(n.i[i]/2), floor(n.i[i]/2))) )
    }
    W.mat <- matrix(0, nrow = N, ncol = S)
    for(i in 1:S){
      W.mat[,i] <- ifelse( trt.vec == 1 & reg.vec == i, 1, 0 )
    }
    
    # Simulate covariates (if any)
    beta.covs <- matrix( 0, nrow = N, ncol = num.covs )
    if(num.covs != 0){
      names(cov.effects) <- cov.names
      colnames(beta.covs) <- cov.names
    }
    if(num.bin != 0){
      for(p in 1:num.bin){
        beta.covs[,p] <- rbinom(N, 1, cov.props[p])
      }
    }
    if(num.con != 0){
      for(p in 1:num.con){
        beta.covs[,num.bin+p] <- rnorm(N, cov.means[p], cov.sds[p])
      }
    }
    covs <- data.frame( trt = W.mat, beta.covs = beta.covs )
    colnames(covs) <- c(paste("trt", 1:S, sep = ""), colnames(beta.covs))
    coefs.vec <- c(trt = reg.te, cov.effects)
    
    # Simulate time from enrollment to theoretical event (up to max.time)
    dat1 <- simsurv( x = covs, hazard = haz, betas = coefs.vec, lb.interval = pw.times,
                     haz.interval = pw.haz, maxt = max.time )
    TTE.from.enrollment <- dat1$eventtime
    event.status <- dat1$status
    
    # Simulate enrollment time (uniformly between 0 and max.enroll.time)
    enroll.time <- runif(N, 0, max.enroll.time)
    
    # Simulate time from enrollment to stochastic censorship time (theoretical dropout time)
    dropout.time <- rexp(N, dropout.rate)
    
    # Calculate event and censorship times from start of study
    event.time <- enroll.time + TTE.from.enrollment
    censored.time <- enroll.time + dropout.time
    
    # Determine observed time and event indicator
    obs.time <- ifelse( event.time < censored.time, event.time, censored.time )
    obs.time <- ifelse( obs.time > max.time, max.time, obs.time )
    event.status <- ifelse( event.time < censored.time & event.time < max.time, 1, 0 )
    
    # Compile dataset
    dat.full <- data.frame( id = dat1$id, eventtime = obs.time, status = event.status,
                            region = reg.vec, trt = trt.vec, enroll.time = enroll.time )
    dat.full$followup.time <- dat.full$eventtime - dat.full$enroll.time
    
    # Create interval bounds such that approximately same number of events are in each interval
    # Identify interval in which each patient failed or was censored
    int.cuts <- c( 0, quantile( dat.full$eventtime[which(dat.full$status == 1)],
                                probs = seq(1/K, 1, length = K)[-K] ) )
    dat.full$interval <- apply( cbind(dat.full$eventtime), 1, function(x){ sum(x > int.cuts) } )
    
    
    
    ## BMA approach using Laplace approximation
    start.BMA.time <- Sys.time()
    bma.sim.results <- bma.asymptotic( dat = dat.full, W.mat = W.mat, int.cuts = int.cuts,
                                       mu0 = mu0, Sig0 = Sig0, eta0 = eta0, phi0 = phi0,
                                       mod.mat = mod.mat, mod.priors = mod.priors, gamma0 = gamma0,
                                       n.draws = n.draws, epsilon.star = epsilon.star,
                                       beta.star = beta.star, pi0 = pi0 )
    end.BMA.time <- Sys.time()
    
    # Save posterior summaries
    sim.time[j,1] <- difftime(end.BMA.time,          # computation time (min) for BMA approach
                              start.BMA.time, units = "mins")
    BMA.gte.mean[j,] <- bma.sim.results$gte.mean     # posterior mean of global treatment effect
    BMA.gte.sd[j,] <- bma.sim.results$gte.sd         # posterior sd of global treatment effect
    BMA.gte.less.gamma0[j,] <- bma.sim.results$gte.less.gamma0    # Pr(gamma_G < gamma0|D)
    BMA.rte.mean[j,] <- bma.sim.results$rte.means    # posterior means of regional trtmt effects
    BMA.rte.sd[j,] <- bma.sim.results$rte.sds        # posterior sd's of regional trtmt effects
    BMA.rte.less.gamma0[j,] <- bma.sim.results$rte.less.gamma0   # Pr(gamma_i < gamma0|D), i=1,...,S
    BMA.cov.mean[j,] <- bma.sim.results$cov.means    # posterior mean of covariate effects
    BMA.cov.sd[j,] <- bma.sim.results$cov.sds        # posterior sd's of covariate effects
    
    # Save PMPs and consistency measures
    BMA.pmp.mat[j,] <- bma.sim.results$PMPs          # posterior model probabilities
    loc.consis.PMDA.TE[j,] <- bma.sim.results$local.consistency.ratio.TE   # Pr(gamma_i/gamma_G > pi0|D)
    loc.consis.PMDA.RR[j,] <- bma.sim.results$local.consistency.ratio.RR   # Pr(RR_i/RR_G > pi0|D)
    # Pr(|gamma_i - gamma_(-i)| < -log(epsilon.star)|D) for all values of epsilon.star
    for(q in 1:num.eps){
      loc.consis.loo[[q]][j,] <- bma.sim.results$local.consistency.loo[q,]
    }
    # Pr(|gamma_i - gamma_j| < -log(epsilon.star)|D) - store results for all epsilon_star in one matrix
    prws.consis.combined.j <- matrix( 0, nrow = num.eps * S + num.eps, ncol = S )
    prws.inconsis.combined.j <- matrix( 0, nrow = num.eps * S + num.eps, ncol = S )
    for(q in 1:num.eps){
      # Pairwise consistency probabilities for q^th value of epsilon.star
      pc.mat.q <- bma.sim.results$pairwise.consistency[[q]]
      prws.consis.combined.j[(q-1)*S+q,1] <- epsilon.star[q]
      prws.consis.combined.j[((q-1)*S+q+1):((q-1)*S+q+S),] <- pc.mat.q
      
      # Pairwise inconsistency probabilities for q^th value of epsilon.star
      pic.mat.q <- bma.sim.results$pairwise.inconsistency[[q]]
      prws.inconsis.combined.j[(q-1)*S+q,1] <- epsilon.star[q]
      prws.inconsis.combined.j[((q-1)*S+q+1):((q-1)*S+q+S),] <- pic.mat.q
    }
    prws.consis.combined <- prws.consis.combined + prws.consis.combined.j
    prws.inconsis.combined <- prws.inconsis.combined + prws.inconsis.combined.j
    glob.consis.prob[j,] <- bma.sim.results$global.consistency     # global consistency probability
    
    
    
    ## Cox proportional hazards models (CPHMs)
    
    # Sufficient statistics
    nu.iq.vec <- numeric(2*S)
    y.iq.vec <- numeric(2*S)
    for(i in 1:S){
      nu.iq.vec[2*i-1] <- sum( dat.full$status[ dat.full$region == i & dat.full$trt == 0 ] )
      nu.iq.vec[2*i] <- sum( dat.full$status[ dat.full$region == i & dat.full$trt == 1 ] )
      y.iq.vec[2*i-1] <- sum( dat.full$eventtime[ dat.full$region == i & dat.full$trt == 0 ] )
      y.iq.vec[2*i] <- sum( dat.full$eventtime[ dat.full$region == i & dat.full$trt == 1 ] )
    }
    
    # Condensed dataset with sufficient statistics
    # Two observations per region (one for each treatment)
    dat.cond <- data.frame( region = rep(1:S, each = 2),
                            trt.group = rep(0:1, S),
                            nu.iq = nu.iq.vec,
                            y.iq = y.iq.vec )
    
    # Obtain posterior summaries of global and region-specific treatment effects
    start.CPHM.time <- Sys.time()
    CPHM.results <- ph.models( dat = dat.cond, mu0.ph1 = mu0.CPHM1, mu0.ph2 = mu0.CPHM2,
                               Sig0.ph1 = Sig0.CPHM1, Sig0.ph2 = Sig0.CPHM2, eta0.ph1 = eta0.CPHM1,
                               eta0.ph2 = eta0.CPHM2, phi0.ph1 = phi0.CPHM1, phi0.ph2 = phi0.CPHM2,
                               sd.MH1 = sd.CPHM1, sd.MH2 = sd.CPHM2, n.draws = n.draws.CPHM,
                               burn.in = burn.in.CPHM, thin = thin.CPHM, gamma0 = gamma0 )
    end.CPHM.time <- Sys.time()
    sim.time[j,2] <- difftime(end.CPHM.time,                  # computation time (min) for CPHMs
                              start.CPHM.time, units = "mins")
    
    # Save posterior summaries
    CPHM.gte.mean[j,] <- CPHM.results$ph1.mean                # posterior mean of global treatment effect
    CPHM.rte.mean[j,] <- CPHM.results$ph2.mean                # posterior means of regional trtmt effects
    CPHM.gte.sd[j,] <- CPHM.results$ph1.sd                    # posterior sd of global treatment effect
    CPHM.rte.sd[j,] <- CPHM.results$ph2.sd                    # posterior sd's of regional trtmt effects
    CPHM.gte.less.gamma0[j,] <- CPHM.results$ph1.less.gamma0  # Pr(gamma_G < gamma0|D)
    CPHM.rte.less.gamma0[j,] <- CPHM.results$ph2.less.gamma0  # Pr(gamma_i < gamma0|D), i=1,...,S
    
    
    
    ## Bayesian hierarchical model (BHM) - equivalent model used by FDA
    
    # Compute observed subgroup log hazard ratios for each region
    log.hr.crude <- numeric(S)
    for(i in 1:S){
      dat.i <- dat.full[ which(dat.full$region == i), ]
      log.hr.crude[i] <- log( sum(dat.i$status[which(dat.i$trt == 1)]) /
                                sum(dat.i$eventtime[which(dat.i$trt == 1)]) /
                                ( sum(dat.i$status[which(dat.i$trt == 0)]) /
                                    sum(dat.i$eventtime[which(dat.i$trt == 0)]) ) )
    }
    
    # Estimate variances used in likelihood as (1/nu_i0 + 1/nu_i1), where nu_i0 and nu_i1
    # are the number of events in the control and treatment group, respectively, for region i
    # (same estimation rule used by FDA)
    sig2.i <- tapply( dat.cond$nu.iq, dat.cond$region, function(x) sum(1/x) )
    
    # Normal hierarchical model with observed subgroup log hazard ratios as response
    data.BHM <- list( y = log.hr.crude, s2i = sig2.i, region = 1:S )
    start.BHM.time <- Sys.time()
    jags.BHM <- jags.model( file = "../R-source/jagsBHM.txt", data = data.BHM,
                            n.adapt = burn.in.BHM, n.chains = 1, quiet = TRUE )
    invisible(capture.output(             # suppress in-function text from displaying in console
      samps.BHM <- coda.samples( jags.BHM, n.iter = n.draws.BHM, thin = thin.BHM,
                                 variable.names = parms.BHM )
    ))
    end.BHM.time <- Sys.time()
    sim.time[j,3] <- difftime(end.BHM.time,                 # computation time (min) for BHM
                              start.BHM.time, units = "mins")
    
    # Save posterior summaries
    BHM.gte.mean[j,] <- mean(samps.BHM[[1]][,1])             # posterior mean of global treatment effect
    BHM.rte.mean[j,] <- colMeans(samps.BHM[[1]][,2:(S+1)])   # posterior means of random trtmt effects
    BHM.gte.sd[j,] <- sd(samps.BHM[[1]][,1])                 # posterior sd of global treatment effect
    BHM.rte.sd[j,] <- apply(samps.BHM[[1]][,2:(S+1)], 2, sd) # posterior sd's of random trtmt effects
    BHM.gte.less.gamma0[j,] <- mean(samps.BHM[[1]][,1] < gamma0)  # Pr(gamma_G < gamma0|D)
    BHM.rte.less.gamma0[j,] <- colMeans(samps.BHM[[1]][,2:(S+1)] < gamma0)  # Pr(gamma_i < gamma0|D)
    
    
    ## Print iterations
    if(print.iters == TRUE){
      print(j)
    }
    
  }
  
  
  ### Proportion of datasets that reject the null
  reject.proportions <- matrix( 0, nrow = 3, ncol = S+1 )
  rownames(reject.proportions) <- c("CPHM", "BMA", "BHM")
  colnames.vec <- numeric(S+1)
  
  ## Global treatment effect
  CPHM.reject.glob <- ifelse( CPHM.gte.less.gamma0 > .975, 1, 0 )
  BMA.reject.glob <- ifelse( BMA.gte.less.gamma0 > .975, 1, 0 )
  BHM.reject.glob <- ifelse( BHM.gte.less.gamma0 > .975, 1, 0 )
  reject.proportions[1,1] <- mean(CPHM.reject.glob)
  reject.proportions[2,1] <- mean(BMA.reject.glob)
  reject.proportions[3,1] <- mean(BHM.reject.glob)
  
  ## Regional treatment effects
  reject.proportions[1, 2:(S+1)] <- colMeans(CPHM.rte.less.gamma0 > .975)
  reject.proportions[2, 2:(S+1)] <- colMeans(BMA.rte.less.gamma0 > .975)
  reject.proportions[3, 2:(S+1)] <- colMeans(BHM.rte.less.gamma0 > .975)
  colnames(reject.proportions) <- c("Global", reg.names)
  
  
  ### Posterior means of regional treatment effects
  rte.post.means <- rbind( colMeans(CPHM.rte.mean), colMeans(BMA.rte.mean),
                           colMeans(BHM.rte.mean) )
  colnames(rte.post.means) <- reg.names
  rownames(rte.post.means) <- c("CPHM", "BMA", "BHM")
  
  
  ### Bias and MSE for regional treatment effects
  CPHM.diff.rte.mat <- matrix( 0, nrow = nsims, ncol = S )
  BMA.diff.rte.mat <- matrix( 0, nrow = nsims, ncol = S )
  BHM.diff.rte.mat <- matrix( 0, nrow = nsims, ncol = S )
  for(i in 1:S){
    CPHM.diff.rte.mat[,i] <- CPHM.rte.mean[,i] - reg.te[i]   # (estimate - true effect) for CPHM
    BMA.diff.rte.mat[,i] <- BMA.rte.mean[,i] - reg.te[i]     # (estimate - true effect) for BMA
    BHM.diff.rte.mat[,i] <- BHM.rte.mean[,i] - reg.te[i]     # (estimate - true effect) for BHM
  }
  reg.trtmt.bias <- rbind( colMeans(CPHM.diff.rte.mat), colMeans(BMA.diff.rte.mat),
                           colMeans(BHM.diff.rte.mat) )
  colnames(reg.trtmt.bias) <- reg.names
  rownames(reg.trtmt.bias) <- c("CPHM", "BMA", "BHM")
  reg.trtmt.mse <- rbind( colMeans(CPHM.diff.rte.mat^2), colMeans(BMA.diff.rte.mat^2),
                          colMeans(BHM.diff.rte.mat^2) )
  colnames(reg.trtmt.mse) <- reg.names
  rownames(reg.trtmt.mse) <- c("CPHM", "BMA", "BHM")
  
  
  ### Average results of pairwise consistency and global consistency matrices
  prws.cons.avg <- prws.consis.combined / nsims
  prws.incons.avg <- prws.inconsis.combined / nsims
  glob.cons.avg <- cbind( epsilon.star, colMeans(glob.consis.prob) )
  colnames(prws.cons.avg) <- reg.names
  colnames(prws.incons.avg) <- reg.names
  colnames(glob.cons.avg) <- c( "epsilon_star", "GlobConsisProb")
  
  
  ### Median results of LOO local consistency
  loc.consis.loo.meds <- matrix( 0, nrow = S, ncol = num.eps )
  for(q in 1:num.eps){
    loc.consis.loo.meds[,q] <- apply( loc.consis.loo[[q]], 2, median )
  }
  rownames(loc.consis.loo.meds) <- reg.names
  colnames(loc.consis.loo.meds) <- epsilon.star
  
  
  ### Rename columns
  colnames(BMA.pmp.mat) <- 1:num.mods
  colnames(BMA.rte.mean) <- reg.names
  colnames(BMA.gte.mean) <- "GlobTrtEffect"
  colnames(loc.consis.PMDA.TE) <- reg.names
  colnames(loc.consis.PMDA.RR) <- reg.names
  
  
  results.list <- list( Rejection.Rate = reject.proportions,
                        Reg.Trtmt.Effect.Post.Means = rte.post.means,
                        Reg.Trtmt.Effect.Bias = reg.trtmt.bias,
                        Reg.Trtmt.Effect.MSE = reg.trtmt.mse,
                        Regional.Trtmt.Effects = BMA.rte.mean,
                        Global.Trtmt.Effect = BMA.gte.mean,
                        Model.PMP.Values = BMA.pmp.mat,
                        Local.Consistency.PMDA.TE = loc.consis.PMDA.TE,
                        Local.Consistency.PMDA.RR = loc.consis.PMDA.RR,
                        Local.Consistency.LOO = loc.consis.loo.meds,
                        Pairwise.Consistency = prws.cons.avg,
                        Pairwise.Inconsistency = prws.incons.avg,
                        Global.Consistency = glob.cons.avg,
                        Computation.Time.Minutes = sim.time )
  return(results.list)
  
}

