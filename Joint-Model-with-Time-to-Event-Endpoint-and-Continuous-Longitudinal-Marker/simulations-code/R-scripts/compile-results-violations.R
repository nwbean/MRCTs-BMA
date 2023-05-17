##################################################################################################
# COMPILE RESULTS (E.G., AVERAGE POSTERIOR MEANS) FROM SIMULATION SCENARIOS IN WHICH MODEL
# ASSUMPTIONS ARE VIOLATED
##################################################################################################



### Set working directory
setwd("/Users/nathanbean/Documents/Dissertation/Project 3/Project3_Simulations/cluster-results")


### Simulation details

# Simulation version
sim.v1 <- "equal-samp-alpha0p5"
sim.v2 <- "equal-samp-log-gamma-re"
sim.v3 <- "equal-samp-nonPH"
sim.v4 <- "equal-samp-inf-cens"
S <- 4      # number of regions

# Total number of result files from simulation
tot.num.results <- 2500                              # use for each scenario
num.scenarios <- 5                                   # number of scenarios "original" versions
num.per.scen <- tot.num.results / num.scenarios      # number of result files per scenario



### Combine results within each scenario for rejection rate, bias, and MSE

## Rejection rate results (v1)
rr.file.1 <- paste("./", sim.v1, "/rejection_rates/rr_", sim.v1, "-", 1, ".csv", sep = "")
rr.scen.mat.1 <- read.csv(rr.file.1, header = TRUE)
rr.scen.mat.1[1:4, 2:(S+2)] <- 0         # matrix to store rejection rates for scenario

## Rejection rate results (v2)
rr.file.2 <- paste("./", sim.v2, "/rejection_rates/rr_", sim.v2, "-", 1, ".csv", sep = "")
rr.scen.mat.2 <- read.csv(rr.file.2, header = TRUE)
rr.scen.mat.2[1:4, 2:(S+2)] <- 0         # matrix to store rejection rates for scenario

## Rejection rate results (v3)
rr.file.3 <- paste("./", sim.v3, "/rejection_rates/rr_", sim.v3, "-", 1, ".csv", sep = "")
rr.scen.mat.3 <- read.csv(rr.file.3, header = TRUE)
rr.scen.mat.3[1:4, 2:(S+2)] <- 0         # matrix to store rejection rates for scenario

## Rejection rate results (v4)
rr.file.4 <- paste("./", sim.v4, "/rejection_rates/rr_", sim.v4, "-", 1, ".csv", sep = "")
rr.scen.mat.4 <- read.csv(rr.file.4, header = TRUE)
rr.scen.mat.4[1:4, 2:(S+2)] <- 0         # matrix to store rejection rates for scenario

## MSE (regional treatment effects) results (v1)
mse.rte.file.1 <- paste("./", sim.v1, "/mse/rte_mse_", sim.v1, "-", 1, ".csv", sep = "")
mse.rte.scen.mat.1 <- read.csv(mse.rte.file.1, header = TRUE)
mse.rte.scen.mat.1[1:4, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario

## MSE (regional treatment effects) results (v2)
mse.rte.file.2 <- paste("./", sim.v2, "/mse/rte_mse_", sim.v2, "-", 1, ".csv", sep = "")
mse.rte.scen.mat.2 <- read.csv(mse.rte.file.2, header = TRUE)
mse.rte.scen.mat.2[1:4, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario

## MSE (regional treatment effects) results (v3)
mse.rte.file.3 <- paste("./", sim.v3, "/mse/rte_mse_", sim.v3, "-", 1, ".csv", sep = "")
mse.rte.scen.mat.3 <- read.csv(mse.rte.file.3, header = TRUE)
mse.rte.scen.mat.3[1:4, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario

## MSE (regional treatment effects) results (v4)
mse.rte.file.4 <- paste("./", sim.v4, "/mse/rte_mse_", sim.v4, "-", 1, ".csv", sep = "")
mse.rte.scen.mat.4 <- read.csv(mse.rte.file.4, header = TRUE)
mse.rte.scen.mat.4[1:4, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario

## Bias (regional treatment effects) results (v1)
bias.file.1 <- paste("./", sim.v1, "/bias/rte_bias_", sim.v1, "-", 1, ".csv", sep = "")
bias.scen.mat.1 <- read.csv(bias.file.1, header = TRUE)
bias.scen.mat.1[1:4, 2:(S+1)] <- 0      # matrix to store rte bias for scenario

## Bias (regional treatment effects) results (v2)
bias.file.2 <- paste("./", sim.v2, "/bias/rte_bias_", sim.v2, "-", 1, ".csv", sep = "")
bias.scen.mat.2 <- read.csv(bias.file.2, header = TRUE)
bias.scen.mat.2[1:4, 2:(S+1)] <- 0      # matrix to store rte bias for scenario

## Bias (regional treatment effects) results (v3)
bias.file.3 <- paste("./", sim.v3, "/bias/rte_bias_", sim.v3, "-", 1, ".csv", sep = "")
bias.scen.mat.3 <- read.csv(bias.file.3, header = TRUE)
bias.scen.mat.3[1:4, 2:(S+1)] <- 0      # matrix to store rte bias for scenario

## Bias (regional treatment effects) results (v4)
bias.file.4 <- paste("./", sim.v4, "/bias/rte_bias_", sim.v4, "-", 1, ".csv", sep = "")
bias.scen.mat.4 <- read.csv(bias.file.4, header = TRUE)
bias.scen.mat.4[1:4, 2:(S+1)] <- 0      # matrix to store rte bias for scenario

## Posterior means
means.scen.mat.1 <- NULL
means.scen.mat.2 <- NULL
means.scen.mat.3 <- NULL
means.scen.mat.4 <- NULL
# ## Posterior means results (v1)
# means.file.1 <- paste("./", sim.v1, "/mean/mean_", sim.v1, "-", 1, ".csv", sep = "")
# means.scen.mat.1 <- read.csv(means.file.1, header = TRUE)
# means.scen.mat.1[1:4, 2:(S+1)] <- 0      # matrix to store posterior means for scenario
# 
# ## Posterior means results (v2)
# means.file.2 <- paste("./", sim.v2, "/mean/mean_", sim.v2, "-", 1, ".csv", sep = "")
# means.scen.mat.2 <- read.csv(means.file.2, header = TRUE)
# means.scen.mat.2[1:4, 2:(S+1)] <- 0      # matrix to store posterior means for scenario
# 
# ## Posterior means results (v3)
# means.file.3 <- paste("./", sim.v3, "/mean/mean_", sim.v3, "-", 1, ".csv", sep = "")
# means.scen.mat.3 <- read.csv(means.file.3, header = TRUE)
# means.scen.mat.3[1:4, 2:(S+1)] <- 0      # matrix to store posterior means for scenario
# 
# ## BPosterior means results (v4)
# means.file.4 <- paste("./", sim.v4, "/mean/mean_", sim.v4, "-", 1, ".csv", sep = "")
# mean.scen.mat.4 <- read.csv(means.file.4, header = TRUE)
# mean.scen.mat.4[1:4, 2:(S+1)] <- 0      # matrix to store posterior means for scenario

rr.list.1 <- list()              # list to store rr.scen.mat for each scenario (v1)
rr.list.2 <- list()              # list to store rr.scen.mat for each scenario (v2)
rr.list.3 <- list()              # list to store rr.scen.mat for each scenario (v3)
rr.list.4 <- list()              # list to store rr.scen.mat for each scenario (v4)
mse.rte.list.1 <- list()         # list to store mse.rte.scen.mat for each scenario (v1)
mse.rte.list.2 <- list()         # list to store mse.rte.scen.mat for each scenario (v2)
mse.rte.list.3 <- list()         # list to store mse.rte.scen.mat for each scenario (v3)
mse.rte.list.4 <- list()         # list to store mse.rte.scen.mat for each scenario (v4)
bias.list.1 <- list()            # list to store bias.scen.mat for each scenario (v1)
bias.list.2 <- list()            # list to store bias.scen.mat for each scenario (v2)
bias.list.3 <- list()            # list to store bias.scen.mat for each scenario (v3)
bias.list.4 <- list()            # list to store bias.scen.mat for each scenario (v4)
means.list.1 <- list()           # list to store means.scen.mat for each scenario (v1)
means.list.2 <- list()           # list to store means.scen.mat for each scenario (v2)
means.list.3 <- list()           # list to store means.scen.mat for each scenario (v3)
means.list.4 <- list()           # list to store means.scen.mat for each scenario (v4)

which.scen <- 1                  # scenario indicator
for(j in 1:tot.num.results){
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen) (v1)
  rr.file.1 <- paste("./", sim.v1, "/rejection_rates/rr_", sim.v1, "-", j, ".csv", sep = "")
  rr.part.mat.1 <- read.csv(rr.file.1, header = TRUE)
  rr.scen.mat.1[1:4, 2:(S+2)] <- rr.scen.mat.1[1:4, 2:(S+2)] + rr.part.mat.1[1:4, 2:(S+2)]
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen) (v2)
  rr.file.2 <- paste("./", sim.v2, "/rejection_rates/rr_", sim.v2, "-", j, ".csv", sep = "")
  rr.part.mat.2 <- read.csv(rr.file.2, header = TRUE)
  rr.scen.mat.2[1:4, 2:(S+2)] <- rr.scen.mat.2[1:4, 2:(S+2)] + rr.part.mat.2[1:4, 2:(S+2)]
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen) (v3)
  rr.file.3 <- paste("./", sim.v3, "/rejection_rates/rr_", sim.v3, "-", j, ".csv", sep = "")
  rr.part.mat.3 <- read.csv(rr.file.3, header = TRUE)
  rr.scen.mat.3[1:4, 2:(S+2)] <- rr.scen.mat.3[1:4, 2:(S+2)] + rr.part.mat.3[1:4, 2:(S+2)]
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen) (v4)
  rr.file.4 <- paste("./", sim.v4, "/rejection_rates/rr_", sim.v4, "-", j, ".csv", sep = "")
  rr.part.mat.4 <- read.csv(rr.file.4, header = TRUE)
  rr.scen.mat.4[1:4, 2:(S+2)] <- rr.scen.mat.4[1:4, 2:(S+2)] + rr.part.mat.4[1:4, 2:(S+2)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen) (v1)
  mse.rte.file.1 <- paste("./", sim.v1, "/mse/rte_mse_", sim.v1, "-", j, ".csv", sep = "")
  mse.rte.part.mat.1 <- read.csv(mse.rte.file.1, header = TRUE)
  mse.rte.scen.mat.1[1:4, 2:(S+1)] <- mse.rte.scen.mat.1[1:4, 2:(S+1)] +
    mse.rte.part.mat.1[1:4, 2:(S+1)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen) (v2)
  mse.rte.file.2 <- paste("./", sim.v2, "/mse/rte_mse_", sim.v2, "-", j, ".csv", sep = "")
  mse.rte.part.mat.2 <- read.csv(mse.rte.file.2, header = TRUE)
  mse.rte.scen.mat.2[1:4, 2:(S+1)] <- mse.rte.scen.mat.2[1:4, 2:(S+1)] +
    mse.rte.part.mat.2[1:4, 2:(S+1)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen) (v3)
  mse.rte.file.3 <- paste("./", sim.v3, "/mse/rte_mse_", sim.v3, "-", j, ".csv", sep = "")
  mse.rte.part.mat.3 <- read.csv(mse.rte.file.3, header = TRUE)
  mse.rte.scen.mat.3[1:4, 2:(S+1)] <- mse.rte.scen.mat.3[1:4, 2:(S+1)] +
    mse.rte.part.mat.3[1:4, 2:(S+1)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen) (v4)
  mse.rte.file.4 <- paste("./", sim.v4, "/mse/rte_mse_", sim.v4, "-", j, ".csv", sep = "")
  mse.rte.part.mat.4 <- read.csv(mse.rte.file.4, header = TRUE)
  mse.rte.scen.mat.4[1:4, 2:(S+1)] <- mse.rte.scen.mat.4[1:4, 2:(S+1)] +
    mse.rte.part.mat.4[1:4, 2:(S+1)]
  
  # Matrix of bias (rte) with part of simulation results (1 out of num.per.scen) (v1)
  bias.file.1 <- paste("./", sim.v1, "/bias/rte_bias_", sim.v1, "-", j, ".csv", sep = "")
  bias.part.mat.1 <- read.csv(bias.file.1, header = TRUE)
  bias.scen.mat.1[1:4, 2:(S+1)] <- bias.scen.mat.1[1:4, 2:(S+1)] +
    bias.part.mat.1[1:4, 2:(S+1)]
  
  # Matrix of bias (rte) with part of simulation results (1 out of num.per.scen) (v2)
  bias.file.2 <- paste("./", sim.v2, "/bias/rte_bias_", sim.v2, "-", j, ".csv", sep = "")
  bias.part.mat.2 <- read.csv(bias.file.2, header = TRUE)
  bias.scen.mat.2[1:4, 2:(S+1)] <- bias.scen.mat.2[1:4, 2:(S+1)] +
    bias.part.mat.2[1:4, 2:(S+1)]
  
  # Matrix of bias (rte) with part of simulation results (1 out of num.per.scen) (v3)
  bias.file.3 <- paste("./", sim.v3, "/bias/rte_bias_", sim.v3, "-", j, ".csv", sep = "")
  bias.part.mat.3 <- read.csv(bias.file.3, header = TRUE)
  bias.scen.mat.3[1:4, 2:(S+1)] <- bias.scen.mat.3[1:4, 2:(S+1)] +
    bias.part.mat.3[1:4, 2:(S+1)]
  
  # Matrix of bias (rte) with part of simulation results (1 out of num.per.scen) (v4)
  bias.file.4 <- paste("./", sim.v4, "/bias/rte_bias_", sim.v4, "-", j, ".csv", sep = "")
  bias.part.mat.4 <- read.csv(bias.file.4, header = TRUE)
  bias.scen.mat.4[1:4, 2:(S+1)] <- bias.scen.mat.4[1:4, 2:(S+1)] +
    bias.part.mat.4[1:4, 2:(S+1)]
  
  # Matrix of posterior means with part of simulation results (1 out of num.per.scen) (v1)
  means.file.1 <- paste("./", sim.v1, "/mean/mean_", sim.v1, "-", j, ".csv", sep = "")
  means.part.mat.1 <- read.csv(means.file.1, header = TRUE)
  means.scen.mat.1 <- rbind(means.scen.mat.1, means.part.mat.1[2,])
  
  # Matrix of posterior means with part of simulation results (1 out of num.per.scen) (v2)
  means.file.2 <- paste("./", sim.v2, "/mean/mean_", sim.v2, "-", j, ".csv", sep = "")
  means.part.mat.2 <- read.csv(means.file.2, header = TRUE)
  means.scen.mat.2 <- rbind(means.scen.mat.2, means.part.mat.2[2,])
  
  # Matrix of posterior means with part of simulation results (1 out of num.per.scen) (v3)
  means.file.3 <- paste("./", sim.v3, "/mean/mean_", sim.v3, "-", j, ".csv", sep = "")
  means.part.mat.3 <- read.csv(means.file.3, header = TRUE)
  means.scen.mat.3 <- rbind(means.scen.mat.3, means.part.mat.3[2,])
  
  # Matrix of posterior means with part of simulation results (1 out of num.per.scen) (v4)
  means.file.4 <- paste("./", sim.v4, "/mean/mean_", sim.v4, "-", j, ".csv", sep = "")
  means.part.mat.4 <- read.csv(means.file.4, header = TRUE)
  means.scen.mat.4 <- rbind(means.scen.mat.4, means.part.mat.4[2,])
  
  # After iterating through num.per.scen files, store results in lists and change scenario indicator
  if(j %% num.per.scen == 0){
    
    # Average results
    rr.scen.mat.1[1:4, 2:(S+2)] <- rr.scen.mat.1[1:4, 2:(S+2)] / num.per.scen
    rr.scen.mat.2[1:4, 2:(S+2)] <- rr.scen.mat.2[1:4, 2:(S+2)] / num.per.scen
    rr.scen.mat.3[1:4, 2:(S+2)] <- rr.scen.mat.3[1:4, 2:(S+2)] / num.per.scen
    rr.scen.mat.4[1:4, 2:(S+2)] <- rr.scen.mat.4[1:4, 2:(S+2)] / num.per.scen
    mse.rte.scen.mat.1[1:4, 2:(S+1)] <- mse.rte.scen.mat.1[1:4, 2:(S+1)] / num.per.scen
    mse.rte.scen.mat.2[1:4, 2:(S+1)] <- mse.rte.scen.mat.2[1:4, 2:(S+1)] / num.per.scen
    mse.rte.scen.mat.3[1:4, 2:(S+1)] <- mse.rte.scen.mat.3[1:4, 2:(S+1)] / num.per.scen
    mse.rte.scen.mat.4[1:4, 2:(S+1)] <- mse.rte.scen.mat.4[1:4, 2:(S+1)] / num.per.scen
    bias.scen.mat.1[1:4, 2:(S+1)] <- bias.scen.mat.1[1:4, 2:(S+1)] / num.per.scen
    bias.scen.mat.2[1:4, 2:(S+1)] <- bias.scen.mat.2[1:4, 2:(S+1)] / num.per.scen
    bias.scen.mat.3[1:4, 2:(S+1)] <- bias.scen.mat.3[1:4, 2:(S+1)] / num.per.scen
    bias.scen.mat.4[1:4, 2:(S+1)] <- bias.scen.mat.4[1:4, 2:(S+1)] / num.per.scen
    
    # Store results in list
    scen.lab <- paste("scenario_", which.scen, sep = "")
    rr.list.1[[scen.lab]] <- rr.scen.mat.1             # save scenario-specific matrix of rr's
    rr.list.2[[scen.lab]] <- rr.scen.mat.2             # save scenario-specific matrix of rr's
    rr.list.3[[scen.lab]] <- rr.scen.mat.3             # save scenario-specific matrix of rr's
    rr.list.4[[scen.lab]] <- rr.scen.mat.4             # save scenario-specific matrix of rr's
    mse.rte.list.1[[scen.lab]] <- mse.rte.scen.mat.1   # save scenario-specific matrix of mse.rte's
    mse.rte.list.2[[scen.lab]] <- mse.rte.scen.mat.2   # save scenario-specific matrix of mse.rte's
    mse.rte.list.3[[scen.lab]] <- mse.rte.scen.mat.3   # save scenario-specific matrix of mse.rte's
    mse.rte.list.4[[scen.lab]] <- mse.rte.scen.mat.4   # save scenario-specific matrix of mse.rte's
    bias.list.1[[scen.lab]] <- bias.scen.mat.1         # save scenario-specific matrix of biases
    bias.list.2[[scen.lab]] <- bias.scen.mat.2         # save scenario-specific matrix of biases
    bias.list.3[[scen.lab]] <- bias.scen.mat.3         # save scenario-specific matrix of biases
    bias.list.4[[scen.lab]] <- bias.scen.mat.4         # save scenario-specific matrix of biases
    means.list.1[[scen.lab]] <- means.scen.mat.1       # save scenario-specific matrix of means
    means.list.2[[scen.lab]] <- means.scen.mat.2       # save scenario-specific matrix of means
    means.list.3[[scen.lab]] <- means.scen.mat.3       # save scenario-specific matrix of means
    means.list.4[[scen.lab]] <- means.scen.mat.4       # save scenario-specific matrix of means
    
    # Reset matrices and prepare for next scenario
    which.scen <- which.scen + 1                     # increase scenario indicator
    rr.scen.mat.1[1:4, 2:(S+2)] <- 0
    rr.scen.mat.2[1:4, 2:(S+2)] <- 0
    rr.scen.mat.3[1:4, 2:(S+2)] <- 0
    rr.scen.mat.4[1:4, 2:(S+2)] <- 0
    mse.rte.scen.mat.1[1:4, 2:(S+1)] <- 0
    mse.rte.scen.mat.2[1:4, 2:(S+1)] <- 0
    mse.rte.scen.mat.3[1:4, 2:(S+1)] <- 0
    mse.rte.scen.mat.4[1:4, 2:(S+1)] <- 0
    bias.scen.mat.1[1:4, 2:(S+1)] <- 0
    bias.scen.mat.2[1:4, 2:(S+1)] <- 0
    bias.scen.mat.3[1:4, 2:(S+1)] <- 0
    bias.scen.mat.4[1:4, 2:(S+1)] <- 0
    means.scen.mat.1 <- NULL
    means.scen.mat.2 <- NULL
    means.scen.mat.3 <- NULL
    means.scen.mat.4 <- NULL
    
  }
  
}



### Operating characteristics

## Rejection rates

# Global rejection rates
grr.v1 <- c( rr.list.1[[1]][2,2], rr.list.1[[2]][2,2], rr.list.1[[3]][2,2],
             rr.list.1[[4]][2,2], rr.list.1[[5]][2,2] )
grr.v2 <- c( rr.list.2[[1]][2,2], rr.list.2[[2]][2,2], rr.list.2[[3]][2,2],
             rr.list.2[[4]][2,2], rr.list.2[[5]][2,2] )
grr.v3 <- c( rr.list.3[[1]][2,2], rr.list.3[[2]][2,2], rr.list.3[[3]][2,2],
             rr.list.3[[4]][2,2], rr.list.3[[5]][2,2] )
grr.v4 <- c( rr.list.4[[1]][2,2], rr.list.4[[2]][2,2], rr.list.4[[3]][2,2],
             rr.list.4[[4]][2,2], rr.list.4[[5]][2,2] )

# True positive rates
tpr.v1 <- c( rowMeans(rr.list.1[[1]][2,3:6]), rowMeans(rr.list.1[[2]][2,4:6]),
             rowMeans(rr.list.1[[3]][2,5:6]), rr.list.1[[4]][2,6], NULL )
tpr.v2 <- c( rowMeans(rr.list.2[[1]][2,3:6]), rowMeans(rr.list.2[[2]][2,4:6]),
             rowMeans(rr.list.2[[3]][2,5:6]), rr.list.2[[4]][2,6], NULL )
tpr.v3 <- c( rowMeans(rr.list.3[[1]][2,3:6]), rowMeans(rr.list.3[[2]][2,4:6]),
             rowMeans(rr.list.3[[3]][2,5:6]), rr.list.3[[4]][2,6], NULL )
tpr.v4 <- c( rowMeans(rr.list.4[[1]][2,3:6]), rowMeans(rr.list.4[[2]][2,4:6]),
             rowMeans(rr.list.4[[3]][2,5:6]), rr.list.4[[4]][2,6], NULL )

# False positive rates
fpr.v1 <- c( NULL, rr.list.1[[2]][2,3], rowMeans(rr.list.1[[3]][2,3:4]),
             rowMeans(rr.list.1[[4]][2,3:5]), rowMeans(rr.list.1[[5]][2,3:6]) )
fpr.v2 <- c( NULL, rr.list.2[[2]][2,3], rowMeans(rr.list.2[[3]][2,3:4]),
             rowMeans(rr.list.2[[4]][2,3:5]), rowMeans(rr.list.2[[5]][2,3:6]) )
fpr.v3 <- c( NULL, rr.list.3[[2]][2,3], rowMeans(rr.list.3[[3]][2,3:4]),
             rowMeans(rr.list.3[[4]][2,3:5]), rowMeans(rr.list.3[[5]][2,3:6]) )
fpr.v4 <- c( NULL, rr.list.4[[2]][2,3], rowMeans(rr.list.4[[3]][2,3:4]),
             rowMeans(rr.list.4[[4]][2,3:5]), rowMeans(rr.list.4[[5]][2,3:6]) )


## Relative MSE for regional treatment effects with v1 as reference

rel.MSE.mat.alt <- matrix(0, nrow = S, ncol = 4)
rel.MSE.mat.null <- matrix(0, nrow = S, ncol = 4)

# Scenario 1: 0 nulls
# alternative regions
rel.MSE.mat.alt[1,] <- c(rowMeans(mse.rte.list.1[[1]][2,2:5]), rowMeans(mse.rte.list.2[[1]][2,2:5]),
                         rowMeans(mse.rte.list.3[[1]][2,2:5]), rowMeans(mse.rte.list.4[[1]][2,2:5])) /
  rowMeans(mse.rte.list.1[[1]][2,2:5])

# Scenario 2: 1 null
# alternative regions
rel.MSE.mat.alt[2,] <- c(rowMeans(mse.rte.list.1[[2]][2,3:5]), rowMeans(mse.rte.list.2[[2]][2,3:5]),
                         rowMeans(mse.rte.list.3[[2]][2,3:5]), rowMeans(mse.rte.list.4[[2]][2,3:5])) /
  rowMeans(mse.rte.list.1[[2]][2,3:5])
# null regions
rel.MSE.mat.null[1,] <- c(mse.rte.list.1[[2]][2,2], mse.rte.list.2[[2]][2,2],
                          mse.rte.list.3[[2]][2,2], mse.rte.list.4[[2]][2,2]) /
  mse.rte.list.1[[2]][2,2]

# Scenario 3: 2 nulls
# alternative regions
rel.MSE.mat.alt[3,] <- c(rowMeans(mse.rte.list.1[[3]][2,4:5]), rowMeans(mse.rte.list.2[[3]][2,4:5]),
                         rowMeans(mse.rte.list.3[[3]][2,4:5]), rowMeans(mse.rte.list.4[[3]][2,4:5])) /
  rowMeans(mse.rte.list.1[[3]][2,4:5])
# null regions
rel.MSE.mat.null[2,] <- c(rowMeans(mse.rte.list.1[[3]][2,2:3]), rowMeans(mse.rte.list.2[[3]][2,2:3]),
                          rowMeans(mse.rte.list.3[[3]][2,2:3]), rowMeans(mse.rte.list.4[[3]][2,2:3])) /
  rowMeans(mse.rte.list.1[[3]][2,2:3])

# Scenario 4: 3 nulls
# alternative regions
rel.MSE.mat.alt[4,] <- c(mse.rte.list.1[[4]][2,5], mse.rte.list.2[[4]][2,5],
                         mse.rte.list.3[[4]][2,5], mse.rte.list.4[[4]][2,5]) /
  mse.rte.list.1[[4]][2,5]
# null regions
rel.MSE.mat.null[3,] <- c(rowMeans(mse.rte.list.1[[4]][2,2:4]), rowMeans(mse.rte.list.2[[4]][2,2:4]),
                          rowMeans(mse.rte.list.3[[4]][2,2:4]), rowMeans(mse.rte.list.4[[4]][2,2:4])) /
  rowMeans(mse.rte.list.1[[4]][2,2:4])

# Scenario 5: 4 nulls
# null regions
rel.MSE.mat.null[4,] <- c(rowMeans(mse.rte.list.1[[5]][2,2:5]), rowMeans(mse.rte.list.2[[5]][2,2:5]),
                          rowMeans(mse.rte.list.3[[5]][2,2:5]), rowMeans(mse.rte.list.4[[5]][2,2:5])) /
  rowMeans(mse.rte.list.1[[5]][2,2:5])



## Bias for region-specific treatment effects

bias.mat.alt <- matrix(0, nrow = S, ncol = 4)
bias.mat.null <- matrix(0, nrow = S, ncol = 4)

# Scenario 1: 0 nulls
# alternative regions
bias.mat.alt[1,] <- c(rowMeans(bias.list.1[[1]][2,2:5]), rowMeans(bias.list.2[[1]][2,2:5]),
                         rowMeans(bias.list.3[[1]][2,2:5]), rowMeans(bias.list.4[[1]][2,2:5]))

# Scenario 2: 1 null
# alternative regions
bias.mat.alt[2,] <- c(rowMeans(bias.list.1[[2]][2,3:5]), rowMeans(bias.list.2[[2]][2,3:5]),
                         rowMeans(bias.list.3[[2]][2,3:5]), rowMeans(bias.list.4[[2]][2,3:5]))
# null regions
bias.mat.null[1,] <- c(bias.list.1[[2]][2,2], bias.list.2[[2]][2,2],
                       bias.list.3[[2]][2,2], bias.list.4[[2]][2,2])

# Scenario 3: 2 nulls
# alternative regions
bias.mat.alt[3,] <- c(rowMeans(bias.list.1[[3]][2,4:5]), rowMeans(bias.list.2[[3]][2,4:5]),
                         rowMeans(bias.list.3[[3]][2,4:5]), rowMeans(bias.list.4[[3]][2,4:5]))
# null regions
bias.mat.null[2,] <- c(rowMeans(bias.list.1[[3]][2,2:3]), rowMeans(bias.list.2[[3]][2,2:3]),
                          rowMeans(bias.list.3[[3]][2,2:3]), rowMeans(bias.list.4[[3]][2,2:3]))

# Scenario 4: 3 nulls
# alternative regions
bias.mat.alt[4,] <- c(bias.list.1[[4]][2,5], bias.list.2[[4]][2,5],
                      bias.list.3[[4]][2,5], bias.list.4[[4]][2,5])
# null regions
bias.mat.null[3,] <- c(rowMeans(bias.list.1[[4]][2,2:4]), rowMeans(bias.list.2[[4]][2,2:4]),
                          rowMeans(bias.list.3[[4]][2,2:4]), rowMeans(bias.list.4[[4]][2,2:4]))

# Scenario 5: 4 nulls
# null regions
bias.mat.null[4,] <- c(rowMeans(bias.list.1[[5]][2,2:5]), rowMeans(bias.list.2[[5]][2,2:5]),
                          rowMeans(bias.list.3[[5]][2,2:5]), rowMeans(bias.list.4[[5]][2,2:5]))



## Mean region-specific treatment effects

mean.mat.alt <- matrix(0, nrow = S, ncol = 4)
mean.mat.null <- matrix(0, nrow = S, ncol = 4)

# Scenario 1: 0 nulls
# alternative regions
mean.mat.alt[1,] <- c(mean(as.matrix(means.list.1[[1]][,1:4])), mean(as.matrix(means.list.2[[1]][,1:4])),
                      mean(as.matrix(means.list.3[[1]][,1:4])), mean(as.matrix(means.list.4[[1]][,1:4])))

# Scenario 2: 1 null
# alternative regions
mean.mat.alt[2,] <- c(mean(as.matrix(means.list.1[[2]][,2:4])), mean(as.matrix(means.list.2[[2]][,2:4])),
                      mean(as.matrix(means.list.3[[2]][,2:4])), mean(as.matrix(means.list.4[[2]][,2:4])))
# null regions
mean.mat.null[1,] <- c(mean(as.matrix(means.list.1[[2]][,1])), mean(as.matrix(means.list.2[[2]][,1])),
                       mean(as.matrix(means.list.3[[2]][,1])), mean(as.matrix(means.list.4[[2]][,1])))

# Scenario 3: 2 nulls
# alternative regions
mean.mat.alt[3,] <- c(mean(as.matrix(means.list.1[[3]][,3:4])), mean(as.matrix(means.list.2[[3]][,3:4])),
                      mean(as.matrix(means.list.3[[3]][,3:4])), mean(as.matrix(means.list.4[[3]][,3:4])))
# null regions
mean.mat.null[2,] <- c(mean(as.matrix(means.list.1[[3]][,1:2])), mean(as.matrix(means.list.2[[3]][,1:2])),
                       mean(as.matrix(means.list.3[[3]][,1:2])), mean(as.matrix(means.list.4[[3]][,1:2])))

# Scenario 4: 3 nulls
# alternative regions
mean.mat.alt[4,] <- c(mean(as.matrix(means.list.1[[4]][,4])), mean(as.matrix(means.list.2[[4]][,4])),
                      mean(as.matrix(means.list.3[[4]][,4])), mean(as.matrix(means.list.4[[4]][,4])))
# null regions
mean.mat.null[3,] <- c(mean(as.matrix(means.list.1[[4]][,1:3])), mean(as.matrix(means.list.2[[4]][,1:3])),
                       mean(as.matrix(means.list.3[[4]][,1:3])), mean(as.matrix(means.list.4[[4]][,1:3])))

# Scenario 5: 4 nulls
# null regions
mean.mat.null[4,] <- c(mean(as.matrix(means.list.1[[5]][,1:4])), mean(as.matrix(means.list.2[[5]][,1:4])),
                       mean(as.matrix(means.list.3[[5]][,1:4])), mean(as.matrix(means.list.4[[5]][,1:4])))



### Bias at times t=0 and t=60 for non-proportional hazards case
te.t0 <- log(.868)                            # treatment effect at time t=0 (alt regions)
te.t60.alt <- log(.868) - abs(log(.868))*2    # treatment effect at time t=60 (alt regions)
te.t60.null <- 0 - abs(log(.868))*2           # treatment effect at time t=60 (null regions)
mean.mat.alt[,3] - te.t0           # bias of alternative regions at t=0
mean.mat.null[,3] - 0              # bias of null regions at t=0
mean.mat.alt[,3] - te.t60.alt      # bias of alternative regions at t=60
mean.mat.null[,3] - te.t60.null    # bias of null regions at t=60



### Plot rejection rates and relative MSE with scenario (number null regions) as x-axis

# Create PDF file
#five.plot.file <- "/Users/nathanbean/Downloads/violations.pdf"
#pdf(five.plot.file, height = 6.46, width = 13)

# Adjust margins and layout
par( mar = c(4.1, 4.1, 1, 1) )
layout( mat=matrix(c(1,1,2,2,3,3,4,5,5,6,6,7), nrow = 2, byrow = TRUE) )

# Global rejection rate plot
grr.mat <- matrix( rbind(grr.v1, grr.v2, grr.v3), byrow = FALSE, nrow = 1 )
barplot( c(grr.mat), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0), 5) ),
         col = rep(c("gray50", "gray70", "gray90"), 6) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16), labels = 0:4, tick = FALSE, line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "Global Rejection Rate", side = 2, line = 2.8, cex = 1.1 )
abline( h = .025, lty = 3, lwd = 2 )
text(.5, .95, "A", cex = 2)
legend( "topright", bty = "n", legend = c("v1", "v2", "v3"),
        fill = c("gray50", "gray70", "gray90"), cex = 1.1 )
box()

# Relative MSE plot for alternative regions
MSE.alt.lim <- 1.8
#MSE.alt.lim <- 1.4    # use for "equal-samp-alpha1"
barplot( c(t(rel.MSE.mat.alt[,-c(1,4)])), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim),
         yaxt = "n", space = rep(c(.5,0), 4), col = rep(c("gray60", "gray90"), 4) )
axis( 1, at = c(1.5, 4, 6.5, 9), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
text(.5, .95*MSE.alt.lim, "B", cex = 2)
legend( "topright", bty = "n", fill = c("gray60", "gray90"), cex = 1.1,
        legend = c("v2 (alt. regions)", "v3 (alt. regions)") )
box()

# Relative MSE plot for null regions
MSE.null.lim <- 1.8
barplot( c(t(rel.MSE.mat.null[,-c(1,4)])), xlab = "", ylab = "", ylim = c(0,MSE.null.lim),
         yaxt = "n", space = rep(c(.5,0), 4), col = rep(c("gray60", "gray90"), 4) )
axis( 1, at = c(1.5, 4, 6.5, 9), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.null.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
abline( h = 1, lty = 3, lwd = 2 )
text(.5, .95*MSE.null.lim, "C", cex = 2)
legend( "topright", bty = "n", fill = c("gray60", "gray90"), cex = 1.1,
        legend = c("v2 (null regions)", "v3 (null regions)") )
box()

# Empty plot
plot.new()

# True positive rate plot
TPR.lim <- .6
tpr.mat <- matrix( rbind(tpr.v1, tpr.v2, tpr.v3), byrow = FALSE, nrow = 1 )
barplot( c(tpr.mat), xlab = "", ylab = "", ylim = c(0,TPR.lim),
         yaxt = "n", space = c( rep(c(.5,0,0), 4) ),
         col = rep(c("gray50", "gray70", "gray90"), 4) )
axis( 1, at = c(2, 5.5, 9, 12.5), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, TPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "True Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(.5, .95*TPR.lim, "D", cex = 2)
legend( "topright", bty = "n", legend = c("v1", "v2", "v3"),
        fill = c("gray50", "gray70", "gray90"), cex = 1.1 )
box()

# False positive rate plot
FPR.lim <- 0.16
fpr.mat <- matrix( rbind(fpr.v1, fpr.v2, fpr.v3), byrow = FALSE, nrow = 1 )
barplot( c(fpr.mat), xlab = "", ylab = "", ylim = c(0,FPR.lim),
         yaxt = "n", space = c( rep(c(.5,0,0), 4) ),
         col = rep(c("gray50", "gray70", "gray90"), 4) )
axis( 1, at = c(2, 5.5, 9, 12.5), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, FPR.lim, by = 0.02), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.8, cex = 1.1 )
mtext( "False Positive Rate", side = 2, line = 3.1, cex = 1.1 )
text(.5, .95*FPR.lim, "E", cex = 2)
abline( h = .025, lty = 3, lwd = 2 )
legend( "topright", bty = "n", legend = c("v1", "v2", "v3"),
        fill = c("gray50", "gray70", "gray90"), cex = 1.1 )
box()

# Empty plot
plot.new()

layout( mat=matrix(1, nrow = 1, ncol = 1) )

# Save PDF
#dev.off()



### Plot rejection rates and relative MSE with scenario (number null regions) as x-axis for the
### case in which the subject-specific random intercepts are not normally distributed

# Create PDF file
five.plot.file <- "/Users/nathanbean/Downloads/SA-log-gamma-re.pdf"
pdf(five.plot.file, height = 6.46, width = 13)

# Adjust margins and layout
par( mar = c(4.1, 4.1, 1, 1) )
layout( mat=matrix(c(1,1,2,2,3,3,4,5,5,6,6,7), nrow = 2, byrow = TRUE) )

# Global rejection rate plot
grr.mat <- matrix( rbind(grr.v1, grr.v2), byrow = FALSE, nrow = 1 )
barplot( c(grr.mat), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0), 5) ),
         col = rep(c("gray60", "gray90"), 6) )
axis( 1, at = c(1.5, 4, 6.5, 9, 11.5), labels = 0:4, tick = FALSE, line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "Global Rejection Rate", side = 2, line = 2.8, cex = 1.1 )
abline( h = .025, lty = 3, lwd = 2 )
text(.5, .95, "A", cex = 2)
legend( "topright", bty = "n", legend = c("normal", "log-gamma"),
        fill = c("gray60", "gray90"), cex = 1.1 )
box()

# Relative MSE plot for alternative regions
MSE.alt.lim <- 1.2    # use for "equal-samp-alpha1"
barplot( c(t(rel.MSE.mat.alt[,-c(1,3,4)])), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim),
         yaxt = "n", space = rep(.5, 4), col = rep(c("gray80"), 4) )
axis( 1, at = c(1, 2.5, 4, 5.5), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
text(.5, .95*MSE.alt.lim, "B", cex = 2)
legend( "topright", bty = "n", fill = c("gray80"), cex = 1.1,
        legend = c("log-gamma (alt. regions)") )
box()

# Relative MSE plot for null regions
MSE.null.lim <- 1.2
barplot( c(t(rel.MSE.mat.null[,-c(1,3,4)])), xlab = "", ylab = "", ylim = c(0,MSE.null.lim),
         yaxt = "n", space = rep(.5, 4), col = rep(c("gray80"), 4) )
axis( 1, at = c(1, 2.5, 4, 5.5), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.null.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
abline( h = 1, lty = 3, lwd = 2 )
text(.5, .95*MSE.null.lim, "C", cex = 2)
legend( "topright", bty = "n", fill = c("gray80"), cex = 1.1,
        legend = c("log-gamma (null regions)") )
box()

# Empty plot
plot.new()

# True positive rate plot
TPR.lim <- .6
tpr.mat <- matrix( rbind(tpr.v1, tpr.v2), byrow = FALSE, nrow = 1 )
barplot( c(tpr.mat), xlab = "", ylab = "", ylim = c(0,TPR.lim),
         yaxt = "n", space = c( rep(c(.5,0), 4) ),
         col = rep(c("gray60", "gray90"), 4) )
axis( 1, at = c(1.5, 4, 6.5, 9), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, TPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "True Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(.5, .95*TPR.lim, "D", cex = 2)
legend( "topright", bty = "n", legend = c("normal", "log-gamma"),
        fill = c("gray60", "gray90"), cex = 1.1 )
box()

# False positive rate plot
FPR.lim <- 0.16
fpr.mat <- matrix( rbind(fpr.v1, fpr.v2), byrow = FALSE, nrow = 1 )
barplot( c(fpr.mat), xlab = "", ylab = "", ylim = c(0,FPR.lim),
         yaxt = "n", space = c( rep(c(.5,0), 4) ),
         col = rep(c("gray60", "gray90"), 4) )
axis( 1, at = c(1.5, 4, 6.5, 9), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, FPR.lim, by = 0.02), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.8, cex = 1.1 )
mtext( "False Positive Rate", side = 2, line = 3.1, cex = 1.1 )
text(.5, .95*FPR.lim, "E", cex = 2)
abline( h = .025, lty = 3, lwd = 2 )
legend( "topright", bty = "n", legend = c("normal", "log-gamma"),
        fill = c("gray60", "gray90"), cex = 1.1 )
box()

# Empty plot
plot.new()

layout( mat=matrix(1, nrow = 1, ncol = 1) )

# Save PDF
dev.off()



### Plot distributions of subject-specific random intercepts for general simulation scenarios
### (normal) and violations scenario (log-gamma)

# Create PDF file
re.plot.file <- "/Users/nathanbean/Downloads/rand-effects-distns.pdf"
pdf(re.plot.file, height = 6.46, width = 10)

# Adjust margins
par( mar = c(4.1, 4.1, 1, 1) )

# Log-gamma pdf
dloggamma <- function(x, a, b){
  b^a * exp(a * x) * exp( -b * exp(x)) / gamma(a)
}

# Distribution values
norm.mean <- 0            # mean of normal distribution
norm.sd <- 1.057          # standard deviation of normal distribution
lg.shape <- 1.314376      # shape parameter of log-gamma distribution
lg.rate <- 0.8581164      # rate parameter of log-gamma distribution

# Plot distributions
x.vals <- seq(-5, 5, by = .01)
lg.vals <- dloggamma(x.vals, lg.shape, lg.rate)   # density values of log-gamma distribution
plot(x.vals, dnorm(x.vals, norm.mean, norm.sd), type = "l", xlim = c(-5, 5), ylim = c(0, .5),
     lwd = 2, col = "blue2", xlab = expression(b[ij]), ylab = "Density")
lines(x.vals, lg.vals, col = "red3", lwd = 2, lty = 2)
abline(v = 0, lty = 3)
legend("topright", c("normal", "log-gamma"), lwd = 2, lty = c(1,2), col = c("blue2", "red3"),
       bty = "n")

# Save PDF
dev.off()

