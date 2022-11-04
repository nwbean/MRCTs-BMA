############################################################################################
# FULL SIMULATION FUNCTION
#
# Compare the BMA survival-only model approach and to two Cox proportional
# hazards models (CPHM)--one with a common treatment effect and one with a
# region-by-treatment interaction (no covariates)--and fit both approaches to
# survival-only data and joint data
############################################################################################


### Source R scripts
source("../R-scripts/data-simulation-functions.R")         # functions to simulate all data
source("../R-scripts/bma-survival-function.R")             # function for BMA (survival only)
source("../R-scripts/proportional-hazards-models.R")       # function for CPHMs


### Load libraries
library(MASS)


### Simulation function
# Return the proportion of simulations that each method rejects the null hypothesis
# S: number of regions
# T0: max number of distinct regions to consider
# n.i: S x 1 vector of regional sample size
# reg.names: S x 1 vector of region names in alphabetical order
# pw.hazs: Q x 1 vector of underlying baseline hazards corresponding to each defined interval
# pw.times: Q x 1 vector of left bounds of intervals, each with its own underlying baseline hazard
# reg.te.surv: S x 1 vector of region-specific treatment effects for survival model
# sigma.e: likelihood standard deviation of longitudinal submodel (i.e., measurement error)
# b.mean: 1- or 2-dimensional vector of mean(s) of subject-specific random intercepts (and slopes)
# b.sd: standard deviation of subject-specific random intercepts (only if REs = "int")
# b.cov: 2 x 2 covariance matrix for subject-specific random effects (only if REs = "int-slope")
# alpha.vec: 1- or 2-dimensional vector of true association parameter(s) between two submodels
# vst.times: K x 1 vector of visit times when longitudinal measures are observed (K: number of visits)
# vst.means.0: K x 1 vector of means of longitudinal measures for each visit (control group)
# vst.means.1: K x 1 vector of means of longitudinal measures for each visit (trtmt group)
# max.time: maximum follow up time (from beginning of study)
# dropout.rate: subject dropout rate
# num.covs.long: number of covariates for long. model (p.long)
# num.covs.surv: number of covariates for survival model (p.surv)
# num.bin.long: number of binary covariates for long. model
# num.bin.surv: number of binary covariates for survival model
# num.con.long: number of continuous covariates for long. model
# num.con.surv: number of continuous covariates for survival model
# cov.long.effects: p.long x 1 vector of true covariate effects for long. model
# cov.surv.effects: p.surv x 1 vector of true covariate effects
# cov.long.names: p.long x 1 vector of labels for covariate names for long. model
# cov.surv.names: p.surv x 1 vector of labels for covariate names
# cov.long.props: p.long x 1 vector of success probs. for long. binary covariates (0 if no binary covs)
# cov.surv.props: p.surv x 1 vector of success probs. for survival binary covariates (0 if no binary covs)
# cov.long.means: p.long x 1 vector of means for long. normal covariates (0 if no continuous covs)
# cov.surv.means: p.surv x 1 vector of means for survival normal covariates (0 if no continuous covs)
# cov.long.sds: p.long x 1 vector of stand. devs. for long. normal covariates (0 if no continuous covs)
# cov.surv.sds: p.surv x 1 vector of stand. devs. for survival normal covariates (0 if no continuous covs)
# Q: number of time intervals in survival model (estimate Q*S constant baseline hazards)
# lin.spline.knots: locations of knots for linear splines (optional)
# a0.surv: pre-specified value for prior model probabilities (survival)
# a0.long: pre-specified value for prior model probabilities (longitudinal)
# mod.priors: L x 1 vector of prior model probabilities, where L is number of models (optional)
# mu0.trt.s: prior mean (scalar) for fixed survival treatment effects
# tau0.trt.s: prior precision (scalar) for fixed survival treatment effects
# mu0.def: default prior mean (scalar) for fixed effects
# tau0.def: default prior precision (scalar) for fixed effects
# eta0.err: prior shape (scalar) for error precision in longitudinal likelihood
# phi0.err: prior rate (scalar) for error precision in longitudinal likelihood
# eta0.rand: prior shape (scalar) for precision of subject-specific random intercepts
# phi0.rand: prior rate (scalar) for precision of subject-specific random intercepts
# mu0.alpha: prior mean (scalar) for association parameter connected to random intercepts
# tau0.alpha: prior precision (scalar) for association parameter connected to random intercepts
# max.iter: maximum number of iterations in optimization algorithm
# optim.toler.params: tolerance for when convergence of parameters has been satisfied in
#      optimization algorithm
# optim.toler.b: tolerance for when convergence of random effects has been satisfied in
#      optimization algorithm
# gamma0.surv: value in hypotheses for posterior probability Pr(gamma > gamma0|D) (longitudianal)
# gamma0.long: value in hypotheses for posterior probability Pr(gamma > gamma0|D) (survival)
# n.draws: number of posterior samples to draw when evaluating consistency
# epsilon.star: vector of possible minimal clinically important regional differences of trtmnt
#      effects - used for consistency measures (used with BMA)
# beta.star: probability cutoff for global inconsistency for which two regions are considered
#      to be clinically different (used with BMA)
# pi0: value for which we want to calculate probabilities with local consistency (used with BMA)
# nsims: number of datasets to simulate
# print.iters: whether or not to print iteration number in R console
run.sims <- function(
  # Data simulation inputs
  S, T0, n.i, reg.names, pw.hazs, pw.times, reg.te.surv, sigma.e, b.mean, b.sd = NULL,
  b.cov = NULL, alpha.vec, vst.times, vst.means.0, vst.means.1, max.time, dropout.rate,
  # Covariate details
  num.covs.long = 0, num.covs.surv = 0, num.bin.long = 0, num.bin.surv = 0,
  num.con.long = 0, num.con.surv = 0, cov.long.effects = NULL, cov.surv.effects = NULL,
  cov.long.names = NULL, cov.surv.names = NULL, cov.long.props = 0, cov.surv.props = 0,
  cov.long.means = 0, cov.surv.means = 0, cov.long.sds = 0, cov.surv.sds = 0,
  # BMA-JM inputs
  Q, lin.spline.knots = NULL, a0.surv = NULL, a0.long = NULL, mod.priors = NULL, mu.Y, Sig.Y,
  mu.X, Sig.X, mu.alpha, Sig.alpha, eta.iq, phi.iq, eta.tau, phi.tau, nu.G, C0.inv.G, REs,
  max.iter = 10, optim.toler.params = .005, optim.toler.b = .05, gamma0.surv, gamma0.long,
  n.draws, epsilon.star, beta.star, pi0,
  # General simulation inputs
  nsims, print.iters = FALSE ){
  
  
  ### Compile matrix of model region classifications (for both BMA approaches)
  mod.mat <- compile_modMat(S, T0)
  
  
  ### Matrices to hold simulation results
  
  # BMA-S: store PMPs for all models (joint data)
  BMA.S.num.mods <- num_models(S, T0)                     # number of models for BMA
  BMA.S.pmp.mat.joint <- matrix(0, nrow = nsims,
                                ncol = BMA.S.num.mods)    # posterior model probabilities
  
  # BMA-S: store means, sd's, and probabilities for global and region-specific treatment effects (joint data)
  BMA.S.gte.mean.joint <- matrix(0, nrow = nsims, ncol = 1)           # global treatment effect mean
  BMA.S.gte.sd.joint <- matrix(0, nrow = nsims, ncol = 1)             # global treatment effect sd
  BMA.S.gte.less.gamma0.joint <- matrix(0, nrow = nsims, ncol = 1)    # Pr(gamma_G < gamma0|D)
  BMA.S.rte.mean.joint <- matrix(0, nrow = nsims, ncol = S)           # regional treatment effect means
  BMA.S.rte.sd.joint <- matrix(0, nrow = nsims, ncol = S)             # regional treatment effect sd's
  BMA.S.rte.less.gamma0.joint <- matrix(0, nrow = nsims, ncol = S)    # Pr(gamma_i < gamma0|D), i=1,...,S
  
  # BMA-S: store PMPs for all models (survival data)
  BMA.S.pmp.mat.surv <- matrix(0, nrow = nsims,
                               ncol = BMA.S.num.mods)     # posterior model probabilities
  
  # BMA-S: store means, sd's, and probabilities for global and region-specific treatment effects (surv. data)
  BMA.S.gte.mean.surv <- matrix(0, nrow = nsims, ncol = 1)           # global treatment effect mean
  BMA.S.gte.sd.surv <- matrix(0, nrow = nsims, ncol = 1)             # global treatment effect sd
  BMA.S.gte.less.gamma0.surv <- matrix(0, nrow = nsims, ncol = 1)    # Pr(gamma_G < gamma0|D)
  BMA.S.rte.mean.surv <- matrix(0, nrow = nsims, ncol = S)           # regional treatment effect means
  BMA.S.rte.sd.surv <- matrix(0, nrow = nsims, ncol = S)             # regional treatment effect sd's
  BMA.S.rte.less.gamma0.surv <- matrix(0, nrow = nsims, ncol = S)    # Pr(gamma_i < gamma0|D), i=1,...,S
  
  # CPHM: store means, sd's, and probabilities for global and region-specific treatment effects (joint data)
  CPHM.gte.mean.joint <- matrix(0, nrow = nsims, ncol = 1)          # global trt effect mean
  CPHM.gte.sd.joint <- matrix(0, nrow = nsims, ncol = 1)            # global trt effect sd
  CPHM.gte.less.gamma0.joint <- matrix(0, nrow = nsims, ncol = 1)   # Pr(gamma_G < gamma0|D)
  CPHM.rte.mean.joint <- matrix(0, nrow = nsims, ncol = S)          # region-spec. trt effect mean
  CPHM.rte.sd.joint <- matrix(0, nrow = nsims, ncol = S)            # region-spec. trt effect sd
  CPHM.rte.less.gamma0.joint <- matrix(0, nrow = nsims, ncol = S)   # Pr(gamma_i < gamma0|D), i=1,...,S
  
  # CPHM: store means, sd's, and probabilities for global and region-specific treatment effects (surv. data)
  CPHM.gte.mean.surv <- matrix(0, nrow = nsims, ncol = 1)          # global trt effect mean
  CPHM.gte.sd.surv <- matrix(0, nrow = nsims, ncol = 1)            # global trt effect sd
  CPHM.gte.less.gamma0.surv <- matrix(0, nrow = nsims, ncol = 1)   # Pr(gamma_G < gamma0|D)
  CPHM.rte.mean.surv <- matrix(0, nrow = nsims, ncol = S)          # region-spec. trt effect mean
  CPHM.rte.sd.surv <- matrix(0, nrow = nsims, ncol = S)            # region-spec. trt effect sd
  CPHM.rte.less.gamma0.surv <- matrix(0, nrow = nsims, ncol = S)   # Pr(gamma_i < gamma0|D), i=1,...,S
  
  # Matrix to store computation times (in minutes) for each model
  sim.time <- matrix(0, nrow = nsims, ncol = 4)
  colnames(sim.time) <- c("BMA-S.joint", "CPHMs-joint", "BMA-S-surv", "CPHMs-surv")
  
  
  ### Construct pieces of data generation that do not change across iterations
  
  # Construct design matrix
  N <- sum(n.i)                                      # total number of subjects
  sub.id.index <- cbind( c(1, cumsum(n.i[-S]) + 1),  # first column is first subject ID in region i
                         c(cumsum(n.i[-S]), N) )     # last column is last subject ID in region i
  trt.vec <- NULL                  # vector of treatment assignments
  reg.vec <- rep(1:S, n.i)         # vector of region labels
  for(i in 1:S){
    trt.vec <- c( trt.vec, rep(c(0,1), c(ceiling(n.i[i]/2), floor(n.i[i]/2))) )
  }
  int.mat <- matrix(0, nrow = N, ncol = S)     # N x S matrix to store region indicators for intercepts
  W.mat <- matrix(0, nrow = N, ncol = S)       # N x S matrix to store region-by-trt indicators
  for(i in 1:S){
    int.mat[,i] <- ifelse( reg.vec == i, 1, 0 )
    W.mat[,i] <- ifelse( trt.vec == 1 & reg.vec == i, 1, 0 )
  }
  
  # Simulate covariates (if any) for longitudinal model
  beta.long.covs <- matrix( 0, nrow = N, ncol = num.covs.long )
  if(num.covs.long != 0){
    names(cov.long.effects) <- cov.long.names
    colnames(beta.long.covs) <- cov.long.names
  }
  if(num.bin.long != 0){
    for(p in 1:num.bin.long){
      beta.long.covs[,p] <- rbinom(N, 1, cov.long.props[p])
    }
  }
  if(num.con.long != 0){
    for(p in 1:num.con.long){
      beta.long.covs[,num.bin.long+p] <- rnorm(N, cov.long.means[p], cov.long.sds[p])
    }
  }
  W.long <- data.frame( int = int.mat, trt = W.mat, covs.long = beta.long.covs )
  
  # Simulate covariates (if any) for survival model
  beta.surv.covs <- matrix( 0, nrow = N, ncol = num.covs.surv )
  if(num.covs.surv != 0){
    names(cov.surv.effects) <- cov.surv.names
    colnames(beta.surv.covs) <- cov.surv.names
  }
  if(num.bin.surv != 0){
    for(p in 1:num.bin.surv){
      beta.surv.covs[,p] <- rbinom(N, 1, cov.surv.props[p])
    }
  }
  if(num.con.surv != 0){
    for(p in 1:num.con.surv){
      beta.Y.covs[,num.bin.surv+p] <- rnorm(N, cov.surv.means[p], cov.surv.sds[p])
    }
  }
  W.surv <- data.frame( trt = W.mat, covs.surv = beta.surv.covs )
  coefs.surv.vec <- c(trt = reg.te.surv, cov = cov.surv.effects)
  
  
  ### Prior model probabilities for BMA-S (survival model only)
  mod.priors.BMA.S <- construct_modPriors(S, T0, alpha = a0.surv)   # prior model probabilities
  
  
  ### Prior hyperparameters and algorithm parameters for CPHMs
  
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
  sd.CPHM1 <- .2            # used in MH algorithm for CPH1 model
  sd.CPHM2 <- .13           # used in MH algorithm for CPH2 model
  
  
  ### Run simulation
  for(m in 1:nsims){
    
    ###  Joint data with random effects
    
    ## Generate dataset
    
    # Simulate subject-specific random effects
    if(REs == "int"){
      b.mat <- cbind( rnorm(N, b.mean, b.sd) )
    } else if(REs == "int-slope"){
      b.mat <- mvrnorm(N, b.mean, b.cov)
    }
    
    # Simulate both survival and longitudinal datasets
    joint.dat <- sim.joint.data( N = N, reg.vec = reg.vec, trt.vec = trt.vec,
                                 base.haz.vec = pw.hazs, pw.times = pw.times,
                                 W.surv = W.surv, W.long = W.long,
                                 coefs.surv.vec = coefs.surv.vec, sigma.e = sigma.e,
                                 b.mat = b.mat, alpha.vec = alpha.vec, vst.times = vst.times,
                                 vst.means.0 = vst.means.0, vst.means.1 = vst.means.1,
                                 dropout.rate = dropout.rate, max.time = max.time )
    surv.dat <- joint.dat$surv.dat          # save survival dataset
    long.dat <- joint.dat$long.dat          # save longitudinal dataset
    
    
    
    ## BMA-S approach using Laplace approximation (survival model only) - joint data
    
    # Create interval bounds such that approximately same number of events are in each interval
    # Identify interval in which each patient failed or was censored
    int.cuts.BMA.S <- c( 0, quantile( surv.dat$eventtime[which(surv.dat$status == 1)],
                                      probs = seq(1/Q, 1, length = Q)[-Q] ) )
    surv.dat$interval <- apply( cbind(surv.dat$eventtime), 1,
                                function(x){ sum(x > int.cuts.BMA.S) } )
    
    # Fit survival-only model using BMA algorithm
    start.BMA.S.time <- Sys.time()
    bma.S.sim.results <- bma.survival( dat = surv.dat, W.mat = W.mat, int.cuts = int.cuts.BMA.S,
                                       mu0 = mu.Y, Sig0 = Sig.Y, eta0 = eta.iq, phi0 = phi.iq,
                                       mod.mat = mod.mat, mod.priors = mod.priors.BMA.S,
                                       gamma0 = gamma0.surv, n.draws = n.draws )
    end.BMA.S.time <- Sys.time()
    
    # Save posterior summaries
    sim.time[m,1] <- difftime(end.BMA.S.time,          # computation time (min) for BMA-S approach
                              start.BMA.S.time, units = "mins")
    BMA.S.gte.mean.joint[m,] <- bma.S.sim.results$gte.mean   # posterior mean of global treatment effect
    BMA.S.gte.sd.joint[m,] <- bma.S.sim.results$gte.sd       # posterior sd of global treatment effect
    BMA.S.gte.less.gamma0.joint[m,] <- bma.S.sim.results$gte.less.gamma0  # Pr(gamma_G < gamma0|D)
    BMA.S.rte.mean.joint[m,] <- bma.S.sim.results$rte.means  # posterior means of regional trtmt effects
    BMA.S.rte.sd.joint[m,] <- bma.S.sim.results$rte.sds      # posterior sd's of regional trtmt effects
    BMA.S.rte.less.gamma0.joint[m,] <- bma.S.sim.results$rte.less.gamma0 # Pr(gamma_i < gamma0|D), i=1,...,S
    BMA.S.pmp.mat.joint[m,] <- bma.S.sim.results$PMPs        # posterior model probabilities
    
    
    
    ## Cox proportional hazards models (CPHMs) (survival model only) - joint data
    
    # Sufficient statistics
    nu.iq.vec <- numeric(2*S)
    y.iq.vec <- numeric(2*S)
    for(i in 1:S){
      nu.iq.vec[2*i-1] <- sum( surv.dat$status[ surv.dat$region == i & surv.dat$trt == 0 ] )
      nu.iq.vec[2*i] <- sum( surv.dat$status[ surv.dat$region == i & surv.dat$trt == 1 ] )
      y.iq.vec[2*i-1] <- sum( surv.dat$eventtime[ surv.dat$region == i & surv.dat$trt == 0 ] )
      y.iq.vec[2*i] <- sum( surv.dat$eventtime[ surv.dat$region == i & surv.dat$trt == 1 ] )
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
                               burn.in = burn.in.CPHM, thin = thin.CPHM, gamma0 = gamma0.surv )
    end.CPHM.time <- Sys.time()
    sim.time[m,2] <- difftime(end.CPHM.time,                  # computation time (min) for CPHMs
                              start.CPHM.time, units = "mins")
    
    # Save posterior summaries
    CPHM.gte.mean.joint[m,] <- CPHM.results$ph1.mean                # posterior mean of global treatment effect
    CPHM.rte.mean.joint[m,] <- CPHM.results$ph2.mean                # posterior means of regional trtmt effects
    CPHM.gte.sd.joint[m,] <- CPHM.results$ph1.sd                    # posterior sd of global treatment effect
    CPHM.rte.sd.joint[m,] <- CPHM.results$ph2.sd                    # posterior sd's of regional trtmt effects
    CPHM.gte.less.gamma0.joint[m,] <- CPHM.results$ph1.less.gamma0  # Pr(gamma_G < gamma0|D)
    CPHM.rte.less.gamma0.joint[m,] <- CPHM.results$ph2.less.gamma0  # Pr(gamma_i < gamma0|D), i=1,...,S
    
    
    
    ###  Survival data without random effects
    
    ## Generate dataset
    
    # Simulate both survival and longitudinal datasets (alpha.vec = 0)
    joint.dat <- sim.joint.data( N = N, reg.vec = reg.vec, trt.vec = trt.vec,
                                 base.haz.vec = pw.hazs, pw.times = pw.times,
                                 W.surv = W.surv, W.long = W.long,
                                 coefs.surv.vec = coefs.surv.vec, sigma.e = sigma.e,
                                 b.mat = b.mat, alpha.vec = 0, vst.times = vst.times,
                                 vst.means.0 = vst.means.0, vst.means.1 = vst.means.1,
                                 dropout.rate = dropout.rate, max.time = max.time )
    surv.dat <- joint.dat$surv.dat          # save survival dataset
    
    
    
    ## BMA-S approach using Laplace approximation (survival model only) - survival data
    
    # Create interval bounds such that approximately same number of events are in each interval
    # Identify interval in which each patient failed or was censored
    int.cuts.BMA.S <- c( 0, quantile( surv.dat$eventtime[which(surv.dat$status == 1)],
                                      probs = seq(1/Q, 1, length = Q)[-Q] ) )
    surv.dat$interval <- apply( cbind(surv.dat$eventtime), 1,
                                function(x){ sum(x > int.cuts.BMA.S) } )
    
    # Fit survival-only model using BMA algorithm
    start.BMA.S.time <- Sys.time()
    bma.S.sim.results <- bma.survival( dat = surv.dat, W.mat = W.mat, int.cuts = int.cuts.BMA.S,
                                       mu0 = mu.Y, Sig0 = Sig.Y, eta0 = eta.iq, phi0 = phi.iq,
                                       mod.mat = mod.mat, mod.priors = mod.priors.BMA.S,
                                       gamma0 = gamma0.surv, n.draws = n.draws )
    end.BMA.S.time <- Sys.time()
    
    # Save posterior summaries
    sim.time[m,3] <- difftime(end.BMA.S.time,          # computation time (min) for BMA-S approach
                              start.BMA.S.time, units = "mins")
    BMA.S.gte.mean.surv[m,] <- bma.S.sim.results$gte.mean   # posterior mean of global treatment effect
    BMA.S.gte.sd.surv[m,] <- bma.S.sim.results$gte.sd       # posterior sd of global treatment effect
    BMA.S.gte.less.gamma0.surv[m,] <- bma.S.sim.results$gte.less.gamma0  # Pr(gamma_G < gamma0|D)
    BMA.S.rte.mean.surv[m,] <- bma.S.sim.results$rte.means  # posterior means of regional trtmt effects
    BMA.S.rte.sd.surv[m,] <- bma.S.sim.results$rte.sds      # posterior sd's of regional trtmt effects
    BMA.S.rte.less.gamma0.surv[m,] <- bma.S.sim.results$rte.less.gamma0 # Pr(gamma_i < gamma0|D), i=1,...,S
    BMA.S.pmp.mat.surv[m,] <- bma.S.sim.results$PMPs        # posterior model probabilities
    
    
    
    ## Cox proportional hazards models (CPHMs) (survival model only) - survival data
    
    # Sufficient statistics
    nu.iq.vec <- numeric(2*S)
    y.iq.vec <- numeric(2*S)
    for(i in 1:S){
      nu.iq.vec[2*i-1] <- sum( surv.dat$status[ surv.dat$region == i & surv.dat$trt == 0 ] )
      nu.iq.vec[2*i] <- sum( surv.dat$status[ surv.dat$region == i & surv.dat$trt == 1 ] )
      y.iq.vec[2*i-1] <- sum( surv.dat$eventtime[ surv.dat$region == i & surv.dat$trt == 0 ] )
      y.iq.vec[2*i] <- sum( surv.dat$eventtime[ surv.dat$region == i & surv.dat$trt == 1 ] )
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
                               burn.in = burn.in.CPHM, thin = thin.CPHM, gamma0 = gamma0.surv )
    end.CPHM.time <- Sys.time()
    sim.time[m,4] <- difftime(end.CPHM.time,                  # computation time (min) for CPHMs
                              start.CPHM.time, units = "mins")
    
    # Save posterior summaries
    CPHM.gte.mean.surv[m,] <- CPHM.results$ph1.mean                # posterior mean of global treatment effect
    CPHM.rte.mean.surv[m,] <- CPHM.results$ph2.mean                # posterior means of regional trtmt effects
    CPHM.gte.sd.surv[m,] <- CPHM.results$ph1.sd                    # posterior sd of global treatment effect
    CPHM.rte.sd.surv[m,] <- CPHM.results$ph2.sd                    # posterior sd's of regional trtmt effects
    CPHM.gte.less.gamma0.surv[m,] <- CPHM.results$ph1.less.gamma0  # Pr(gamma_G < gamma0|D)
    CPHM.rte.less.gamma0.surv[m,] <- CPHM.results$ph2.less.gamma0  # Pr(gamma_i < gamma0|D), i=1,...,S
    
    
    
    ### Print iterations
    if(print.iters == TRUE){
      print(m)
    }
    
  }
  
  
  ### Proportion of datasets that reject the null (surv.)
  reject.proportions.surv <- matrix( 0, nrow = 4, ncol = S+1 )
  rownames(reject.proportions.surv) <- c("CPHM-surv", "CPHM-joint", "BMA-S-surv", "BMA-S-joint")
  colnames.vec <- numeric(S+1)
  
  ## Global treatment effect (surv.)
  CPHM.reject.glob.joint <- ifelse( CPHM.gte.less.gamma0.joint > .975, 1, 0 )
  CPHM.reject.glob.surv <- ifelse( CPHM.gte.less.gamma0.surv > .975, 1, 0 )
  BMA.S.reject.glob.joint <- ifelse( BMA.S.gte.less.gamma0.joint > .975, 1, 0 )
  BMA.S.reject.glob.surv <- ifelse( BMA.S.gte.less.gamma0.surv > .975, 1, 0 )
  reject.proportions.surv[1,1] <- mean(CPHM.reject.glob.surv)
  reject.proportions.surv[2,1] <- mean(CPHM.reject.glob.joint)
  reject.proportions.surv[3,1] <- mean(BMA.S.reject.glob.surv)
  reject.proportions.surv[4,1] <- mean(BMA.S.reject.glob.joint)
  
  ## Regional treatment effects (surv.)
  reject.proportions.surv[1, 2:(S+1)] <- colMeans(CPHM.rte.less.gamma0.surv > .975)
  reject.proportions.surv[2, 2:(S+1)] <- colMeans(CPHM.rte.less.gamma0.joint > .975)
  reject.proportions.surv[3, 2:(S+1)] <- colMeans(BMA.S.rte.less.gamma0.surv > .975)
  reject.proportions.surv[4, 2:(S+1)] <- colMeans(BMA.S.rte.less.gamma0.joint > .975)
  colnames(reject.proportions.surv) <- c("Global", reg.names)
  
  
  ### Posterior means of regional treatment effects (surv.)
  rte.post.means.surv <- rbind( colMeans(CPHM.rte.mean.surv), colMeans(CPHM.rte.mean.joint),
                                colMeans(BMA.S.rte.mean.surv), colMeans(BMA.S.rte.mean.joint) )
  colnames(rte.post.means.surv) <- reg.names
  rownames(rte.post.means.surv) <- c("CPHM-surv", "CPHM-joint", "BMA-S-surv", "BMA-S-joint")
  
  
  ### Posterior standard deviations of regional treatment effects (surv.)
  rte.post.sds.surv <- rbind( colMeans(CPHM.rte.sd.surv), colMeans(CPHM.rte.sd.joint),
                                colMeans(BMA.S.rte.sd.surv), colMeans(BMA.S.rte.sd.joint) )
  colnames(rte.post.sds.surv) <- reg.names
  rownames(rte.post.sds.surv) <- c("CPHM-surv", "CPHM-joint", "BMA-S-surv", "BMA-S-joint")
  
  
  ### Bias and MSE for regional treatment effects (surv.)
  CPHM.diff.rte.mat.joint <- matrix( 0, nrow = nsims, ncol = S )
  CPHM.diff.rte.mat.surv <- matrix( 0, nrow = nsims, ncol = S )
  BMA.S.diff.rte.mat.joint <- matrix( 0, nrow = nsims, ncol = S )
  BMA.S.diff.rte.mat.surv <- matrix( 0, nrow = nsims, ncol = S )
  for(i in 1:S){
    CPHM.diff.rte.mat.joint[,i] <- CPHM.rte.mean.joint[,i] -
      reg.te.surv[i]   # (estimate - true effect) for CPHM - joint   
    CPHM.diff.rte.mat.surv[,i] <- CPHM.rte.mean.surv[,i] -
      reg.te.surv[i]   # (estimate - true effect) for CPHM - surv
    BMA.S.diff.rte.mat.joint[,i] <- BMA.S.rte.mean.joint[,i] -
      reg.te.surv[i] # (estimate - true effect) for BMA-S - joint
    BMA.S.diff.rte.mat.surv[,i] <- BMA.S.rte.mean.surv[,i] -
      reg.te.surv[i] # (estimate - true effect) for BMA-S - surv
  }
  reg.trtmt.bias.surv <- rbind( colMeans(CPHM.diff.rte.mat.surv), colMeans(CPHM.diff.rte.mat.joint),
                                colMeans(BMA.S.diff.rte.mat.surv), colMeans(BMA.S.diff.rte.mat.joint) )
  colnames(reg.trtmt.bias.surv) <- reg.names
  rownames(reg.trtmt.bias.surv) <- c("CPHM-surv", "CPHM-joint", "BMA-S-surv", "BMA-S-joint")
  reg.trtmt.mse.surv <- rbind( colMeans(CPHM.diff.rte.mat.surv^2), colMeans(CPHM.diff.rte.mat.joint^2),
                               colMeans(BMA.S.diff.rte.mat.surv^2), colMeans(BMA.S.diff.rte.mat.joint^2) )
  colnames(reg.trtmt.mse.surv) <- reg.names
  rownames(reg.trtmt.mse.surv) <- c("CPHM-surv", "CPHM-joint", "BMA-S-surv", "BMA-S-joint")
  
  
  ### Return results
  results.list <- list( Rejection.Rate.Surv = reject.proportions.surv,
                        Reg.Trtmt.Effect.Post.Means.Surv = rte.post.means.surv,
                        Reg.Trtmt.Effect.Post.SDs.Surv = rte.post.sds.surv,
                        Reg.Trtmt.Effect.Bias.Surv = reg.trtmt.bias.surv,
                        Reg.Trtmt.Effect.MSE.Surv = reg.trtmt.mse.surv,
                        Computation.Time.Minutes = sim.time )
  return(results.list)
  
}

