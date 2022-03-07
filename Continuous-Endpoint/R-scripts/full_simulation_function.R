############################################################################################
# FULL SIMULATION FUNCTION
#
# Compare BMA to Bayesian hierarchical model (BHM) with treatment as a random effect and
# fixed effect linear regression model (FELM) with common treatment effect and region as
# a covariate
#
# Save all PMPs from BMA to use when determining a guideline for assessing consistency
############################################################################################



### Simulation function
# Return the proportion of simulations that each method rejects the null hypothesis
run.sims <- function(N, S, T0, reg.names, reg.allctn, trtmt.allctn, trtmt.means, cntrl.means,
                     sd.Y, num.bin, num.con, cov.props, cov.means, cov.sds, modPriors, m0,
                     Sig0, delta0, nu0, gamma0, epsilon.star, beta.star, pi0, n.draws,
                     num.sims, print.iters = FALSE){
  
  ### Matrices to hold simulation results
  # Store hypothesis test results
  BMA.results <- matrix( 0, nrow = num.sims, ncol = S+1 )     # store probabilities Pr(gamma > gamma0|D)
  FELM.results <- matrix( 0, nrow = num.sims, ncol = S+1 )    # store lower bound of conf. ints for trtmt effects
  BHM.results <- matrix( 0, nrow = num.sims, ncol = S+1 )     # store lower bound of cred. ints for trtmt effects
  
  # Store PMP for both the correct model and the most likely model, and PMPs for all models
  BMA.mod.results <- matrix(0, nrow = num.sims, ncol = 3)
  colnames(BMA.mod.results) <- c("Correct.PMP", "Most.Likely.Model", "MLM.PMP")
  num.mods <- num_models(S, T0)
  BMA.PMP.vals <- matrix(0, nrow = num.sims, ncol = num.mods)
  
  # Store point estimates of regional treatment effects for each method
  BMA.rte.estimates <- matrix( 0, nrow = num.sims, ncol = S )
  FELM.rte.estimates <- matrix( 0, nrow = num.sims, ncol = S )
  BHM.rte.estimates <- matrix( 0, nrow = num.sims, ncol = S )
  
  # Store point estimates of regional intercepts for each method
  BMA.ri.estimates <- matrix( 0, nrow = num.sims, ncol = S )
  FELM.ri.estimates <- matrix( 0, nrow = num.sims, ncol = S )
  BHM.ri.estimates <- matrix( 0, nrow = num.sims, ncol = S )
  
  # Store point estimates of global treatment effect for BMA
  BMA.gte.estimates <- matrix(0, nrow = num.sims, ncol = 1)
  
  # Store local consistency probabilities (MHLW and leave-one-out)
  num.eps <- length(epsilon.star)
  loc.consis.MHLW <- matrix(0, nrow = num.sims, ncol = S)
  loc.consis.loo <- list()    # Store results in separate matrix for each epsilon.star value
  for(q in 1:num.eps){
    loc.consis.loo[[q]] <- matrix(0, nrow = num.sims, ncol = S)
  }
  
  # Store pairwise consistency/inconsistency probabilities and
  # epsilon.star-level global inconsistency probabilities
  prws.consis.combined <- matrix( 0, nrow = num.eps * S + num.eps, ncol = S )
  prws.inconsis.combined <- matrix( 0, nrow = num.eps * S + num.eps, ncol = S )
  glob.consis.prob <- matrix( 0, nrow = num.sims, ncol = num.eps)
  glob.inconsis.prob <- matrix( 0, nrow = num.sims, ncol = num.eps)
  
  
  ### Compile matrix of model region classifications (for BMA)
  mod.mat <- compile_modMat(S, T0)
  
  
  ### Run simulation
  for(j in 1:num.sims){
    
    ## Generate dataset
    single.sim.dat <- generate_data(N, S, reg.names, reg.allctn, trtmt.allctn, trtmt.means, cntrl.means,
                                    sd.Y, num.bin, num.con, cov.props, cov.means, cov.sds)
    Y <- single.sim.dat$Y
    region.labs <- single.sim.dat$RegionLabels
    trtmt.indctr <- single.sim.dat$TreatmentIndicator
    X <- single.sim.dat$Covariates
    
    
    ## BMA
    bma.sim.results <- bma_single(Y, region.labs, trtmt.indctr, X, T0, modPriors, m0, Sig0, delta0, nu0,
                                  gamma0, epsilon.star, beta.star, pi0, n.draws)
    BMA.results[j,1] <- bma.sim.results$GlobalTrtmtStats[3,1]
    for(i in 1:S){
      BMA.results[j,i+1] <- bma.sim.results$RegionalTrtmtStats[i,3]
    }
    BMA.gte.estimates[j,1] <- bma.sim.results$GlobalTrtmtStats[1,1]     # global treatment effect
    BMA.rte.estimates[j,] <- bma.sim.results$RegionalTrtmtStats[,1]     # posterior means of regional trtmt effects
    BMA.ri.estimates[j,] <- bma.sim.results$RegionalIntStats[,1]        # posterior means of regional intercepts
    BMA.PMP.vals[j,] <- bma.sim.results$pmp                             # PMP values for all models
    loc.consis.MHLW[j,] <- bma.sim.results$RegionalTrtmtStats[,4]       # Pr(gamma_i/gamma > pi0|D)
    # Pr(|gamma_i - gamma_(-i)| < epsilon.star|D) for all values of epsilon.star
    for(q in 1:num.eps){
      loc.consis.loo[[q]][j,] <- bma.sim.results$LOOLocalConsistency[,q]
    }
    #loc.consis.loo <- loc.consis.loo + bma.sim.results$LOOLocalConsistency
    # Pr(|gamma_i - gamma_j| < epsilon.star|D) - store results for all epsilon_star in one matrix
    prws.consis.combined.j <- matrix( 0, nrow = num.eps * S + num.eps, ncol = S )
    prws.inconsis.combined.j <- matrix( 0, nrow = num.eps * S + num.eps, ncol = S )
    for(q in 1:num.eps){
      # Pairwise consistency probabilities for q^th value of epsilon.star
      pc.mat.q <- bma.sim.results$PairwiseConsistency[[q]]
      prws.consis.combined.j[(q-1)*S+q,1] <- epsilon.star[q]
      prws.consis.combined.j[((q-1)*S+q+1):((q-1)*S+q+S),] <- pc.mat.q
      
      # Pairwise inconsistency probabilities for q^th value of epsilon.star
      pic.mat.q <- bma.sim.results$PairwiseInconsistency[[q]]
      prws.inconsis.combined.j[(q-1)*S+q,1] <- epsilon.star[q]
      prws.inconsis.combined.j[((q-1)*S+q+1):((q-1)*S+q+S),] <- pic.mat.q
    }
    prws.consis.combined <- prws.consis.combined + prws.consis.combined.j
    prws.inconsis.combined <- prws.inconsis.combined + prws.inconsis.combined.j
    glob.consis.prob[j,] <- bma.sim.results$GlobalConsistencyProb        # global consistency probability
    glob.inconsis.prob[j,] <- bma.sim.results$GlobalInconsistencyProb    # global inconsistency probability
    
    
    ## Prepare input data for all JAGS models
    unique.region <- sort( unique(region.labs) )
    region.number <- numeric(N)
    for(i in 1:length(unique.region)){
      region.number <- ifelse( region.labs == unique.region[i], i, region.number )
    }
    data.jags <- list(Y = Y, region.number = region.number, trtmt.indctr = trtmt.indctr)
    
    
    ## Fixed effects linear regression model for global test using JAGS ("rjags" package)
    # Fixed effects for intercept, regions, and treatment
    seed.inits1 <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 1*num.sims + j)
    jagsFELMglob.sim <- jags.model(file = "../R-source/jagsFELM_glob.txt", data = data.jags,
                                   inits = seed.inits1, n.adapt = 1000, n.chains = 4, quiet = TRUE)
    FELM.glob.vars <- numeric(S+2)
    FELM.glob.vars[1] <- "beta0"
    FELM.glob.vars[S+2] <- "gamma"
    for(i in 1:S){
      FELM.glob.vars[i+1] <- paste("betaReg[", i, "]", sep = "")
    }
    invisible(capture.output(             # suppress in-function text from displaying in console
      jagsFELMglob.samps <- jags.samples(jagsFELMglob.sim, n.iter = 2500, thin = 1, quiet = TRUE,
                                         variable.names = FELM.glob.vars)
    ))
    FELM.results[j,1] <- quantile( jagsFELMglob.samps$gamma, .025 )
    
    
    ## Fixed effects linear regression model for regional tests using JAGS ("rjags" package)
    # Fixed effects for intercept, regions, treatment, and region*treatment interactions
    seed.inits2 <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2*num.sims + j)
    jagsFELMreg.sim <- jags.model(file = "../R-source/jagsFELM_reg.txt", data = data.jags,
                                  inits = seed.inits2, n.adapt = 1000, n.chains = 4, quiet = TRUE)
    FELM.reg.vars <- numeric(2*S + 2)
    FELM.reg.vars[1] <- "beta0"
    FELM.reg.vars[S+2] <- "gamma"
    for(i in 1:S){
      FELM.reg.vars[i+1] <- paste("betaReg[", i, "]", sep = "")
      FELM.reg.vars[i+S+2] <- paste("gammaReg[", i, "]", sep = "")
    }
    invisible(capture.output(             # suppress in-function text from displaying in console
      jagsFELMreg.samps <- jags.samples(jagsFELMreg.sim, n.iter = 2500, thin = 1, quiet = TRUE,
                                        variable.names = FELM.reg.vars)
    ))
    for(i in 1:S){
      FELM.results[j,i+1] <- quantile( jagsFELMreg.samps$gamma + jagsFELMreg.samps[[i+S+2]], .025 )
      FELM.ri.estimates[j,i] <- mean( jagsFELMreg.samps$beta0 + jagsFELMreg.samps[[i+1]] )
      FELM.rte.estimates[j,i] <- mean( jagsFELMreg.samps$gamma + jagsFELMreg.samps[[i+S+2]] )
    }
    
    
    ## Bayesian hierarchical model using JAGS ("rjags" package)
    # Fixed intercept and fixed treatment slope
    # Random intercept and random treatment slope for each region
    seed.inits3 <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 3*num.sims + j)
    jagsBHM.sim <- jags.model(file = "../R-source/jagsBHM.txt", data = data.jags, inits = seed.inits3,
                              n.adapt = 1000, n.chains = 4, quiet = TRUE)
    BHM.vars <- numeric(2*S + 2)
    BHM.vars[1] <- "beta0"
    BHM.vars[S+2] <- "beta1"
    for(i in 1:S){
      BHM.vars[i+1] <- paste("rand0[", i, "]", sep = "")
      BHM.vars[i+S+2] <- paste("rand1[", i, "]", sep = "")
    }
    invisible(capture.output(             # suppress in-function text from displaying in console
      jagsBHM.samps <- jags.samples(jagsBHM.sim, n.iter = 2500, thin = 1, quiet = TRUE,
                                    variable.names = BHM.vars)
    ))
    BHM.results[j,1] <- quantile( jagsBHM.samps$beta1, .025 )
    for(i in 1:S){
      BHM.results[j,i+1] <- quantile( jagsBHM.samps$beta1 + jagsBHM.samps[[i+S+2]], .025 )
      BHM.ri.estimates[j,i] <- mean( jagsBHM.samps$beta0 + jagsBHM.samps[[i+2]] )
      BHM.rte.estimates[j,i] <- mean( jagsBHM.samps$beta1 + jagsBHM.samps[[i+S+2]] )
    }
    
    if(print.iters == TRUE){
      print(j)
    }
    
  }
  
  
  ### Proportion of datasets that reject the null
  reject.proportions <- matrix( 0, nrow = 3, ncol = S+1 )
  rownames(reject.proportions) <- c("BMA", "FELM", "BHM")
  colnames.vec <- numeric(S+1)
  colnames.vec[1] <- "Global"
  
  ## Global treatment effect
  BMA.reject.glob <- ifelse( BMA.results[,1] > .975, 1, 0 )
  FELM.reject.glob <- ifelse( FELM.results[,1] > 0, 1, 0 )
  BHM.reject.glob <- ifelse( BHM.results[,1] > 0, 1, 0 )
  reject.proportions[1,1] <- mean(BMA.reject.glob)
  reject.proportions[2,1] <- mean(FELM.reject.glob)
  reject.proportions[3,1] <- mean(BHM.reject.glob)
  
  ## Regional treatment effects
  for(i in 1:S){
    colnames.vec[i+1] <- reg.names[i]
    BMA.reject.reg <- ifelse( BMA.results[,i+1] > .975, 1, 0 )
    FELM.reject.reg <- ifelse( FELM.results[,i+1] > 0, 1, 0 )
    BHM.reject.reg <- ifelse( BHM.results[,i+1] > 0, 1, 0 )
    reject.proportions[1,i+1] <- mean(BMA.reject.reg)
    reject.proportions[2,i+1] <- mean(FELM.reject.reg)
    reject.proportions[3,i+1] <- mean(BHM.reject.reg)
  }
  colnames(reject.proportions) <- colnames.vec
  
  
  ### Average point estimates of regional treatment effects
  avg.rte.estimates <- rbind( colMeans(BMA.rte.estimates),
                              colMeans(FELM.rte.estimates),
                              colMeans(BHM.rte.estimates) )
  colnames(avg.rte.estimates) <- reg.names
  rownames(avg.rte.estimates) <- c("BMA", "FELM", "BHM")
  
  
  ### Average point estimates of regional intercepts
  avg.ri.estimates <- rbind( colMeans(BMA.ri.estimates),
                             colMeans(FELM.ri.estimates),
                             colMeans(BHM.ri.estimates) )
  colnames(avg.ri.estimates) <- reg.names
  rownames(avg.ri.estimates) <- c("BMA", "FELM", "BHM")
  
  
  ### Bias and MSE for regional treatment effects
  BMA.bias.rte.mat <- matrix( 0, nrow = num.sims, ncol = S )
  FELM.bias.rte.mat <- matrix( 0, nrow = num.sims, ncol = S )
  BHM.bias.rte.mat <- matrix( 0, nrow = num.sims, ncol = S )
  for(i in 1:S){
    true.reg.te <- trtmt.means[i] - cntrl.means[i]                  # true treatment effect for region i
    BMA.bias.rte.mat[,i] <- BMA.rte.estimates[,i] - true.reg.te     # (estimate - true effect) for BMA
    FELM.bias.rte.mat[,i] <- FELM.rte.estimates[,i] - true.reg.te   # (estimate - true effect) for FELM
    BHM.bias.rte.mat[,i] <- BHM.rte.estimates[,i] - true.reg.te     # (estimate - true effect) for BHM
  }
  reg.trtmt.bias <- rbind( colMeans(BMA.bias.rte.mat),
                           colMeans(FELM.bias.rte.mat),
                           colMeans(BHM.bias.rte.mat) )
  colnames(reg.trtmt.bias) <- reg.names
  rownames(reg.trtmt.bias) <- c("BMA", "FELM", "BHM")
  reg.trtmt.mse <- rbind( colMeans(BMA.bias.rte.mat^2),
                          colMeans(FELM.bias.rte.mat^2),
                          colMeans(BHM.bias.rte.mat^2) )
  colnames(reg.trtmt.mse) <- reg.names
  rownames(reg.trtmt.mse) <- c("BMA", "FELM", "BHM")
  
  
  ### Bias and MSE for regional intercepts
  BMA.bias.ri.mat <- matrix( 0, nrow = num.sims, ncol = S )
  FELM.bias.ri.mat <- matrix( 0, nrow = num.sims, ncol = S )
  BHM.bias.ri.mat <- matrix( 0, nrow = num.sims, ncol = S )
  for(i in 1:S){
    true.reg.int <- cntrl.means[i]                                 # true intercept for region i
    BMA.bias.ri.mat[,i] <- BMA.ri.estimates[,i] - true.reg.int     # (estimate - true intercept) for BMA
    FELM.bias.ri.mat[,i] <- FELM.ri.estimates[,i] - true.reg.int   # (estimate - true intercept) for FELM
    BHM.bias.ri.mat[,i] <- BHM.ri.estimates[,i] - true.reg.int     # (estimate - true intercept) for BHM
  }
  reg.int.bias <- rbind( colMeans(BMA.bias.ri.mat),
                         colMeans(FELM.bias.ri.mat),
                         colMeans(BHM.bias.ri.mat) )
  colnames(reg.int.bias) <- reg.names
  rownames(reg.int.bias) <- c("BMA", "FELM", "BHM")
  reg.int.mse <- rbind( colMeans(BMA.bias.ri.mat^2),
                        colMeans(FELM.bias.ri.mat^2),
                        colMeans(BHM.bias.ri.mat^2) )
  colnames(reg.int.mse) <- reg.names
  rownames(reg.int.mse) <- c("BMA", "FELM", "BHM")
  
  
  ### Average posterior model probability for all models
  BMA.avg.PMP <- cbind( 1:num.mods, colMeans(BMA.PMP.vals) )
  colnames(BMA.avg.PMP) <- c("Model", "AveragePMP")
  
  
  ### Average results of pairwise consistency and global consistency matrices
  prws.cons.avg <- prws.consis.combined / num.sims
  prws.incons.avg <- prws.inconsis.combined / num.sims
  glob.cons.avg <- cbind( epsilon.star, colMeans(glob.consis.prob) )
  glob.incons.avg <- cbind( epsilon.star, colMeans(glob.inconsis.prob) )
  colnames(prws.cons.avg) <- reg.names
  colnames(prws.incons.avg) <- reg.names
  colnames(glob.cons.avg) <- c( "epsilon_star", "GlobConsisProb")
  colnames(glob.incons.avg) <- c( "epsilon_star", "GlobInconsisProb")
  
  
  ### Median results of LOO local consistency
  loc.consis.loo.meds <- matrix( 0, nrow = S, ncol = num.eps )
  for(q in 1:num.eps){
    loc.consis.loo.meds[,q] <- apply( loc.consis.loo[[q]], 2, median )
  }
  rownames(loc.consis.loo.meds) <- reg.names
  colnames(loc.consis.loo.meds) <- epsilon.star
  
  
  ### Rename columns
  colnames(BMA.PMP.vals) <- 1:num.mods
  colnames(BMA.rte.estimates) <- reg.names
  colnames(BMA.gte.estimates) <- "GlobTrtEffect"
  colnames(loc.consis.MHLW) <- reg.names
  
  
  results.list <- list( Rejection.Rate = reject.proportions,
                        Avg.Trtmt.Effect.Estimates = avg.rte.estimates,
                        Avg.Reg.Intercept.Estimates = avg.ri.estimates,
                        Avg.PMP = BMA.avg.PMP,
                        Reg.Trtmt.Effect.Bias = reg.trtmt.bias,
                        Reg.Trtmt.Effect.MSE = reg.trtmt.mse,
                        Reg.Intercept.Bias = reg.int.bias,
                        Reg.Intercept.MSE = reg.int.mse,
                        
                        Regional.Intercepts = BMA.ri.estimates,
                        Regional.Trtmt.Effects = BMA.rte.estimates,
                        Global.Trtmt.Effect = BMA.gte.estimates,
                        Model.PMP.Values = BMA.PMP.vals,
                        Local.Consistency.MHLW = loc.consis.MHLW,
                        Local.Consistency.LOO = loc.consis.loo.meds,
                        Pairwise.Consistency = prws.cons.avg,
                        Pairwise.Inconsistency = prws.incons.avg,
                        Global.Inconsistency = glob.incons.avg,
                        Global.Consistency = glob.cons.avg)
  return(results.list)
  
}

