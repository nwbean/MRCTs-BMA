##################################################################################################
# COMPILE RESULTS FROM SIMULATIONS
##################################################################################################



### Set working directory
setwd("/Users/nathanbean/Documents/Dissertation/Project 1/Project1_Simulations_v2/cluster-results")


### Simulation details

# Simulation version
sim.version <- "n1508"
#sim.version <- "n1508-beta20"
#sim.version <- "n7650"
#sim.version <- "n7650-beta20"
#sim.version <- "n754"
#sim.version <- "n754-beta20"
#sim.version <- "null-half"
#sim.version <- "alt-half"
#sim.version <- "diff-effects-1"
#sim.version <- "diff-effects-2"
#sim.version <- "diff-effects-3"
#sim.version <- "SA1"
#sim.version <- "n1508-alphaNeg10"
#sim.version <- "n1508-alphaNeg4"
#sim.version <- "n1508-alphaNeg2"
#sim.version <- "n1508-alpha2"
#sim.version <- "n1508-alpha4"
#sim.version <- "n1508-alpha10"
S <- 5      # number of regions

# Total number of result files from simulation
tot.num.results <- 60
#tot.num.results <- 1200                              # use for versions "n7650" and "7650-beta20"
#tot.num.results <- 20                                # use for version "diff-effects-1-3"
num.scenarios <- 6                                    # number of unique scenarios
#num.scenarios <- 1                                   # use for version "diff-effects-1-3"
num.per.scen <- tot.num.results / num.scenarios       # number of result files per scenario



### Combine results within each scenario for rejection rate, bias, and MSE

## Rejection rate results
rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", 1, ".csv", sep = "")
rr.scen.mat <- read.csv(rr.file, header = TRUE)
rr.scen.mat[1:3, 2:(S+2)] <- 0                         # matrix to store rejection rates for scenario

## Bias (regional treatment effects) results
bias.rte.file <- paste("./", sim.version, "/bias/rte_bias_", sim.version, "-", 1, ".csv", sep = "")
bias.rte.scen.mat <- read.csv(bias.rte.file, header = TRUE)
bias.rte.scen.mat[1:3, 2:(S+1)] <- 0                   # matrix to store rte biases for scenario

## MSE (regional treatment effects) results
mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", 1, ".csv", sep = "")
mse.rte.scen.mat <- read.csv(mse.rte.file, header = TRUE)
mse.rte.scen.mat[1:3, 2:(S+1)] <- 0                    # matrix to store rte MSEs for scenario

## Bias (regional intercepts) results
bias.ri.file <- paste("./", sim.version, "/bias/ri_bias_", sim.version, "-", 1, ".csv", sep = "")
bias.ri.scen.mat <- read.csv(bias.ri.file, header = TRUE)
bias.ri.scen.mat[1:3, 2:(S+1)] <- 0                    # matrix to store rte biases for scenario

## MSE (regional intercepts) results
mse.ri.file <- paste("./", sim.version, "/mse/ri_mse_", sim.version, "-", 1, ".csv", sep = "")
mse.ri.scen.mat <- read.csv(mse.ri.file, header = TRUE)
mse.ri.scen.mat[1:3, 2:(S+1)] <- 0                     # matrix to store ri MSEs for scenario


gte.scen.mat <- NULL
rte.scen.mat <- NULL
pmp.scen.mat <- NULL
lc.MHLW.scen.mat <- NULL
lc.file <- paste("./", sim.version, "/loc_consis/lc_loo_", sim.version, "-", 1, ".csv", sep = "")
lc.example.mat <- read.csv(lc.file, header = TRUE)
lc.loo.scen.mat <- matrix( 0, nrow = dim(lc.example.mat)[1], ncol = dim(lc.example.mat)[2] )
pc.file <- paste("./", sim.version, "/prws_consis/pc_", sim.version, "-", 1, ".csv", sep = "")
pc.example.mat <- read.csv(pc.file, header = TRUE)
pc.scen.mat <- matrix( 0, nrow = dim(pc.example.mat)[1], ncol = dim(pc.example.mat)[2] )

gte.list <- list()                # list to store gte.scen.mat for each scenario
rte.list <- list()                # list to store rte.scen.mat for each scenario
rr.list <- list()                 # list to store rr.scen.mat for each scenario
bias.rte.list <- list()           # list to store bias.rte.scen.mat for each scenario
mse.rte.list <- list()            # list to store mse.rte.scen.mat for each scenario
bias.ri.list <- list()            # list to store bias.rte.scen.mat for each scenario
mse.ri.list <- list()             # list to store mse.ri.scen.mat for each scenario
pmp.list <- list()                # list to store pmp.scen.mat for each scenario
lc.MHLW.list <- list()            # list to store lc.MHLW.scen.mat for each scenario
lc.loo.list <- list()             # list to store lc.loo.scen.mat (leave one out) for each scenario
pc.list <- list()                 # list to store pc.scen.mat for each scenario


which.scen <- 1                   # scenario indicator
for(j in 1:tot.num.results){
  
  # Matrix of global treatment effects with part of simulation results (1 out of num.per.scen)
  gte.file <- paste("./", sim.version, "/gte/gte_", sim.version, "-", j, ".csv", sep = "")
  gte.part.mat <- read.csv(gte.file, header = TRUE)
  gte.scen.mat <- rbind(gte.scen.mat, gte.part.mat)
  
  # Matrix of regional treatment effects with part of simulation results (1 out of num.per.scen)
  rte.file <- paste("./", sim.version, "/rte/rte_", sim.version, "-", j, ".csv", sep = "")
  rte.part.mat <- read.csv(rte.file, header = TRUE)
  rte.scen.mat <- rbind(rte.scen.mat, rte.part.mat)
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen)
  rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", j, ".csv", sep = "")
  rr.part.mat <- read.csv(rr.file, header = TRUE)
  rr.scen.mat[1:3, 2:(S+2)] <- rr.scen.mat[1:3, 2:(S+2)] + rr.part.mat[1:3, 2:(S+2)]
  
  # Matrix of bias (rte) with part of simulation results (1 out of num.per.scen)
  bias.rte.file <- paste("./", sim.version, "/bias/rte_bias_", sim.version, "-", j, ".csv", sep = "")
  bias.rte.part.mat <- read.csv(bias.rte.file, header = TRUE)
  bias.rte.scen.mat[1:3, 2:(S+1)] <- bias.rte.scen.mat[1:3, 2:(S+1)] + bias.rte.part.mat[1:3, 2:(S+1)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen)
  mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", j, ".csv", sep = "")
  mse.rte.part.mat <- read.csv(mse.rte.file, header = TRUE)
  mse.rte.scen.mat[1:3, 2:(S+1)] <- mse.rte.scen.mat[1:3, 2:(S+1)] + mse.rte.part.mat[1:3, 2:(S+1)]
  
  # Matrix of bias (ri) with part of simulation results (1 out of num.per.scen)
  bias.ri.file <- paste("./", sim.version, "/bias/ri_bias_", sim.version, "-", j, ".csv", sep = "")
  bias.ri.part.mat <- read.csv(bias.ri.file, header = TRUE)
  bias.ri.scen.mat[1:3, 2:(S+1)] <- bias.ri.scen.mat[1:3, 2:(S+1)] + bias.ri.part.mat[1:3, 2:(S+1)]
  
  # Matrix of MSE (ri) with part of simulation results (1 out of num.per.scen)
  mse.ri.file <- paste("./", sim.version, "/mse/ri_mse_", sim.version, "-", j, ".csv", sep = "")
  mse.ri.part.mat <- read.csv(mse.ri.file, header = TRUE)
  mse.ri.scen.mat[1:3, 2:(S+1)] <- mse.ri.scen.mat[1:3, 2:(S+1)] + mse.ri.part.mat[1:3, 2:(S+1)]
  
  # Matrix of PMPs with part of simulation results (1 out of num.per.scen)
  pmp.file <- paste("./", sim.version, "/pmp/pmp_", sim.version, "-", j, ".csv", sep = "")
  pmp.part.mat <- read.csv(pmp.file, header = TRUE)
  pmp.scen.mat <- rbind(pmp.scen.mat, pmp.part.mat)
  
  # Matrix of local consistency probabilities (MHLW) with part of simulation results (1 out of num.per.scen)
  lc.MHLW.file <- paste("./", sim.version, "/loc_consis/lc_MHLW_", sim.version, "-", j, ".csv", sep = "")
  lc.MHLW.part.mat <- read.csv(lc.MHLW.file, header = TRUE)
  lc.MHLW.scen.mat <- rbind(lc.MHLW.scen.mat, lc.MHLW.part.mat)
  
  # Matrix of local consistency probabilities (LOO) with part of simulation results (1 out of num.per.scen)
  lc.loo.file <- paste("./", sim.version, "/loc_consis/lc_loo_", sim.version, "-", j, ".csv", sep = "")
  lc.loo.part.mat <- read.csv(lc.loo.file, header = TRUE)
  lc.loo.scen.mat <- lc.loo.scen.mat + lc.loo.part.mat
  
  # Matrix of pairwise consistency probabilities with part of simulation results (1 out of num.per.scen)
  pc.file <- paste("./", sim.version, "/prws_consis/pc_", sim.version, "-", j, ".csv", sep = "")
  pc.part.mat <- read.csv(pc.file, header = TRUE)
  pc.scen.mat <- pc.scen.mat + pc.part.mat
  
  # After iterating through num.per.scen files, store results in lists and change scenario indicator
  if(j %% num.per.scen == 0){
    
    # Average results
    rr.scen.mat[1:3, 2:(S+2)] <- rr.scen.mat[1:3, 2:(S+2)] / num.per.scen
    bias.rte.scen.mat[1:3, 2:(S+1)] <- bias.rte.scen.mat[1:3, 2:(S+1)] / num.per.scen
    mse.rte.scen.mat[1:3, 2:(S+1)] <- mse.rte.scen.mat[1:3, 2:(S+1)] / num.per.scen
    bias.ri.scen.mat[1:3, 2:(S+1)] <- bias.ri.scen.mat[1:3, 2:(S+1)] / num.per.scen
    mse.ri.scen.mat[1:3, 2:(S+1)] <- mse.ri.scen.mat[1:3, 2:(S+1)] / num.per.scen
    
    # Store results in list
    scen.lab <- paste("scenario_", which.scen, sep = "")
    gte.list[[scen.lab]] <- gte.scen.mat               # save scenario-specific matrix of gte's
    rte.list[[scen.lab]] <- rte.scen.mat               # save scenario-specific matrix of rte's
    rr.list[[scen.lab]] <- rr.scen.mat                 # save scenario-specific matrix of rr's
    bias.rte.list[[scen.lab]] <- bias.rte.scen.mat     # save scenario-specific matrix of bias.rte's
    mse.rte.list[[scen.lab]] <- mse.rte.scen.mat       # save scenario-specific matrix of mse.rte's
    bias.ri.list[[scen.lab]] <- bias.ri.scen.mat       # save scenario-specific matrix of bias.ri's
    mse.ri.list[[scen.lab]] <- mse.ri.scen.mat         # save scenario-specific matrix of mse.ri's
    pmp.list[[scen.lab]] <- pmp.scen.mat               # save scenario-specific matrix of PMPs
    lc.MHLW.list[[scen.lab]] <- lc.MHLW.scen.mat       # save scenario-specific matrix of cp's
    lc.loo.list[[scen.lab]] <- lc.loo.scen.mat /       # save scenario-specific matrix of cp's (LOO)
      num.per.scen
    
    # Save scenario-specific matrix of pc's (for all values of epsilon.star)
    pc.list[[scen.lab]] <- pc.scen.mat / num.per.scen
    
    # Reset matrices and prepare for next scenario
    which.scen <- which.scen + 1                       # increase scenario indicator
    gte.scen.mat <- NULL
    rte.scen.mat <- NULL
    rr.scen.mat[1:3, 2:(S+2)] <- 0
    bias.rte.scen.mat[1:3, 2:(S+1)] <- 0
    mse.rte.scen.mat[1:3, 2:(S+1)] <- 0
    bias.ri.scen.mat[1:3, 2:(S+1)] <- 0
    mse.ri.scen.mat[1:3, 2:(S+1)] <- 0
    pmp.scen.mat <- NULL
    lc.MHLW.scen.mat <- NULL
    lc.loo.scen.mat <- matrix( 0, nrow = dim(lc.loo.scen.mat)[1], ncol = dim(lc.loo.scen.mat)[2] )
    pc.scen.mat <- matrix( 0, nrow = dim(pc.scen.mat)[1], ncol = dim(pc.scen.mat)[2] )
    
  }
  
}



### Combine rejection rates according to FPR and TPR for each scenario

# Scenario 1: 0 nulls
round( rr.list[[1]][,2], 3 )                      # Global rejection rates
round( rowMeans(rr.list[[1]][,3:7]), 3 )          # TPR averages

# Scenario 2: 1 null
round( rr.list[[2]][,2], 3 )                      # Global rejection rates
round( rowMeans(rr.list[[2]][,4:7]), 3 )          # TPR averages
round( rr.list[[2]][,3], 3 )                      # FPR averages

# Scenario 3: 2 nulls
round( rr.list[[3]][,2], 3 )                      # Global rejection rates
round( rowMeans(rr.list[[3]][,5:7]), 3 )          # TPR averages
round( rowMeans(rr.list[[3]][,3:4]), 3 )          # FPR averages

# Scenario 4: 3 nulls
round( rr.list[[4]][,2], 3 )                      # Global rejection rates
round( rowMeans(rr.list[[4]][,6:7]), 3 )          # TPR averages
round( rowMeans(rr.list[[4]][,3:5]), 3 )          # FPR averages

# Scenario 5: 4 nulls
round( rr.list[[5]][,2], 3 )                      # Global rejection rates
round( rr.list[[5]][,7], 3 )                      # TPR averages
round( rowMeans(rr.list[[5]][,3:6]), 3 )          # FPR averages

# Scenario 6: 5 nulls
round( rr.list[[6]][,2], 3 )                      # Global rejection rates
round( rowMeans(rr.list[[6]][,3:7]), 3 )          # FPR averages


# Global rejection rates
grr.FELM <- c( rr.list[[1]][2,2], rr.list[[2]][2,2], rr.list[[3]][2,2],
               rr.list[[4]][2,2], rr.list[[5]][2,2], rr.list[[6]][2,2] )
grr.BMA <- c( rr.list[[1]][1,2], rr.list[[2]][1,2], rr.list[[3]][1,2],
              rr.list[[4]][1,2], rr.list[[5]][1,2], rr.list[[6]][1,2] )
grr.BHM <- c( rr.list[[1]][3,2], rr.list[[2]][3,2], rr.list[[3]][3,2],
              rr.list[[4]][3,2], rr.list[[5]][3,2], rr.list[[6]][3,2] )

# True postive rates
tpr.FELM <- c( rowMeans(rr.list[[1]][2,3:7]), rowMeans(rr.list[[2]][2,4:7]),
               rowMeans(rr.list[[3]][2,5:7]), rowMeans(rr.list[[4]][2,6:7]),
               rr.list[[5]][2,7], NULL )
tpr.BMA <- c( rowMeans(rr.list[[1]][1,3:7]), rowMeans(rr.list[[2]][1,4:7]),
              rowMeans(rr.list[[3]][1,5:7]), rowMeans(rr.list[[4]][1,6:7]),
              rr.list[[5]][1,7], NULL )
tpr.BHM <- c( rowMeans(rr.list[[1]][3,3:7]), rowMeans(rr.list[[2]][3,4:7]),
              rowMeans(rr.list[[3]][3,5:7]), rowMeans(rr.list[[4]][3,6:7]),
              rr.list[[5]][3,7], NULL )

# False postive rates
fpr.FELM <- c( NULL, rr.list[[2]][2,3],
               rowMeans(rr.list[[3]][2,3:4]), rowMeans(rr.list[[4]][2,3:5]),
               rowMeans(rr.list[[5]][2,3:6]), rowMeans(rr.list[[6]][2,3:7]) )
fpr.BMA <- c( NULL, rr.list[[2]][1,3],
              rowMeans(rr.list[[3]][1,3:4]), rowMeans(rr.list[[4]][1,3:5]),
              rowMeans(rr.list[[5]][1,3:6]), rowMeans(rr.list[[6]][1,3:7]) )
fpr.BHM <- c( NULL, rr.list[[2]][3,3],
              rowMeans(rr.list[[3]][3,3:4]), rowMeans(rr.list[[4]][3,3:5]),
              rowMeans(rr.list[[5]][3,3:6]), rowMeans(rr.list[[6]][3,3:7]) )



### Combine regional treatment effect biases according to null and alternative regions for each scenario

# Scenario 1: 0 nulls
round( rowMeans(bias.rte.list[[1]][,2:6]), 7 )    # alternative regions

# Scenario 2: 1 null
round( rowMeans(bias.rte.list[[2]][,3:6]), 7 )    # alternative regions
round( bias.rte.list[[2]][,2], 6 )                # null regions

# Scenario 3: 2 nulls
round( rowMeans(bias.rte.list[[3]][,4:6]), 7 )    # alternative regions
round( rowMeans(bias.rte.list[[3]][,2:3]), 7 )    # null regions

# Scenario 4: 3 nulls
round( rowMeans(bias.rte.list[[4]][,5:6]), 6 )    # alternative regions
round( rowMeans(bias.rte.list[[4]][,2:4]), 7 )    # null regions

# Scenario 5: 4 nulls
round( bias.rte.list[[5]][,6], 7 )                # alternative regions
round( rowMeans(bias.rte.list[[5]][,2:5]), 7 )    # null regions

# Scenario 6: 5 nulls
round( rowMeans(bias.rte.list[[6]][,2:6]), 7 )    # null regions



## Relative MSE for regional treatment effects with FELM as reference

rel.MSE.mat.alt <- matrix(0, nrow = S, ncol = 3)
rel.MSE.mat.null <- matrix(0, nrow = S, ncol = 3)

# Scenario 1: 0 nulls
# alternative regions
rel.MSE.mat.alt[1,] <- rowMeans(mse.rte.list[[1]][,2:6]) / rowMeans(mse.rte.list[[1]][2,2:6])

# Scenario 2: 1 null
# alternative regions
rel.MSE.mat.alt[2,] <- rowMeans(mse.rte.list[[2]][,3:6]) / rowMeans(mse.rte.list[[2]][2,3:6])
# null regions
rel.MSE.mat.null[1,] <- mse.rte.list[[2]][,2] / mse.rte.list[[2]][2,2]

# Scenario 3: 2 nulls
# alternative regions
rel.MSE.mat.alt[3,] <- rowMeans(mse.rte.list[[3]][,4:6]) / rowMeans(mse.rte.list[[3]][2,4:6])
# null regions
rel.MSE.mat.null[2,] <- rowMeans(mse.rte.list[[3]][,2:3]) / rowMeans(mse.rte.list[[3]][2,2:3])

# Scenario 4: 3 nulls
# alternative regions
rel.MSE.mat.alt[4,] <- rowMeans(mse.rte.list[[4]][,5:6]) / rowMeans(mse.rte.list[[4]][2,5:6])
# null regions
rel.MSE.mat.null[3,] <- rowMeans(mse.rte.list[[4]][,2:4]) / rowMeans(mse.rte.list[[4]][2,2:4])

# Scenario 5: 4 nulls
# alternative regions
rel.MSE.mat.alt[5,] <- mse.rte.list[[5]][,6] / mse.rte.list[[5]][2,6]
# null regions
rel.MSE.mat.null[4,] <- rowMeans(mse.rte.list[[5]][,2:5]) / rowMeans(mse.rte.list[[5]][2,2:5])

# Scenario 6: 5 nulls
# null regions
rel.MSE.mat.null[5,] <- rowMeans(mse.rte.list[[6]][,2:6]) / rowMeans(mse.rte.list[[6]][2,2:6])



### Plot rejection rates and relative MSE with scenario (number null regions) as x-axis

par( mar = c(4.1, 4.1, 1, 1) )
layout( mat=matrix(c(1,1,2,2,3,3,4,5,5,6,6,7), nrow = 2, byrow = TRUE) )

# Global rejection rate plot
grr.mat <- matrix( rbind(grr.FELM, grr.BMA, grr.BHM), byrow = FALSE, nrow = 1 )
barplot( c(grr.mat), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0), 6) ),
         col = rep(c("gray50", "gray70", "gray90"), 6) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16, 19.5), labels = 0:5, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "Global Rejection Rate", side = 2, line = 2.8, cex = 1.1 )
abline( h = .025, lty = 3, lwd = 2 )
text(.5, .95, "A", cex = 2)
legend( "topright", bty = "n", legend = c("FELM", "BMA", "BHM"),
        fill = c("gray50","gray70", "gray90"), cex = 1.1 )
box()

# True positive rate plot
TPR.lim <- 1
tpr.mat <- matrix( rbind(tpr.FELM, tpr.BMA, tpr.BHM), byrow = FALSE, nrow = 1 )
barplot( c(tpr.mat), xlab = "", ylab = "", ylim = c(0,TPR.lim),
         yaxt = "n", space = c( rep(c(.5,0,0), 5) ),
         col = rep(c("gray50", "gray70", "gray90"), 6) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16), labels = 0:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, TPR.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "True Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(.5, .95*TPR.lim, "B", cex = 2)
legend( "topright", bty = "n", legend = c("FELM", "BMA", "BHM"),
        fill = c("gray50","gray70", "gray90"), cex = 1.1 )
box()

# False positive rate plot
if( sim.version == "n1508" )  FPR.lim <- 0.6
if( sim.version == "null-half" )  FPR.lim <- 0.8
if( sim.version == "alt-half" )  FPR.lim <- 0.5
fpr.mat <- matrix( rbind(fpr.FELM, fpr.BMA, fpr.BHM), byrow = FALSE, nrow = 1 )
barplot( c(fpr.mat), xlab = "", ylab = "", ylim = c(0,FPR.lim),
         yaxt = "n", space = c( rep(c(.5,0,0), 5) ),
         col = rep(c("gray50", "gray70", "gray90"), 6) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16), labels = 1:5, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
if( sim.version == "n1508" ){
  axis( 2, at = seq(0, FPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
} else if( sim.version == "alt-half" ){
  axis( 2, at = seq(0, FPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
} else{
  axis( 2, at = seq(0, FPR.lim, by = 0.2), las = 1, cex.axis = 1.2 )
}
mtext( "Number of Null Regions", side = 1, line = 2.8, cex = 1.1 )
mtext( "False Positive Rate", side = 2, line = 2.6, cex = 1.1 )
text(.5, .95*FPR.lim, "C", cex = 2)
abline( h = .025, lty = 3, lwd = 2 )
legend( "topright", bty = "n", legend = c("FELM", "BMA", "BHM"),
        fill = c("gray50","gray70", "gray90"), cex = 1.1 )
box()


# Empty plot
plot.new()


# Relative MSE plot for alternative regions
if( sim.version == "n1508" )  MSE.alt.lim <- 1.5
if( sim.version == "null-half" )  MSE.alt.lim <- 1.75
if( sim.version == "alt-half" )  MSE.alt.lim <- 1.2
barplot( c(t(rel.MSE.mat.alt[,-2])), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim),
         yaxt = "n", space = rep(c(.5,0), 5), col = rep(c("gray60", "gray90"), 5) )
axis( 1, at = c(1.5, 4, 6.5, 9, 11.5), labels = 0:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
if( sim.version == "n1508" ){
  axis( 2, at = seq(0, MSE.alt.lim, by = 0.25), las = 1, cex.axis = 1.2 )
  mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
} else if( sim.version == "null-half" ){
  axis( 2, at = seq(0, MSE.alt.lim, by = 0.25), las = 1, cex.axis = 1.2 )
  mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
} else if( sim.version == "alt-half" ){
  axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
  mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
} else{
  axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
  mtext( "Relative MSE", side = 2, line = 2.8, cex = 1.1 )
}
text(.5, .95*MSE.alt.lim, "D", cex = 2)
legend( "top", bty = "n", legend = c("BMA (alt. regions)", "BHM (alt. regions)"),
        fill = c("gray60", "gray90"), cex = 1.1 )
box()

# Relative MSE plot for null regions
if( sim.version == "n1508" )  MSE.null.lim <- 1.5
if( sim.version == "null-half" )  MSE.null.lim <- 1.2
if( sim.version == "alt-half" )  MSE.null.lim <- 1.75
barplot( c(t(rel.MSE.mat.null[,-2])), xlab = "", ylab = "", ylim = c(0,MSE.null.lim),
         yaxt = "n", space = rep(c(.5,0), 5), col = rep(c("gray60", "gray90"), 5) )
axis( 1, at = c(1.5, 4, 6.5, 9, 11.5), labels = 1:5, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
if( sim.version == "n1508" ){
  axis( 2, at = seq(0, MSE.null.lim, by = 0.25), las = 1, cex.axis = 1.2 )
  mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
} else if( sim.version == "null-half" ){
  axis( 2, at = seq(0, MSE.null.lim, by = 0.2), las = 1, cex.axis = 1.2 )
  mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
} else if( sim.version == "alt-half" ){
  axis( 2, at = seq(0, MSE.null.lim, by = 0.25), las = 1, cex.axis = 1.2 )
  mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
} else{
  axis( 2, at = seq(0, MSE.null.lim, by = 0.2), las = 1, cex.axis = 1.2 )
  mtext( "Relative MSE", side = 2, line = 2.8, cex = 1.1 )
}
abline( h = 1, lty = 3, lwd = 2 )
text(.5, .95*MSE.null.lim, "E", cex = 2)
legend( "topright", bty = "n", legend = c("BMA (null regions)", "BHM (null regions)"),
        fill = c("gray60", "gray90"), cex = 1.1 )
box()

# Empty plot
plot.new()

layout( mat=matrix(1, nrow = 1, ncol = 1) )





### Plot results for sim.version = "diff-effects-1-3"

# Save rejection rates and relative MSE for each scenario
if(sim.version == "diff-effects-1"){
  rr.mat1 <- as.matrix( rr.list[[1]][c(2,1,3),-1] )
  rel.MSE.mat.diff1 <- matrix( 0, nrow = 3, ncol = S )
  for(i in 1:S)  rel.MSE.mat.diff1[,i] <- mse.rte.list[[1]][,i+1] / mse.rte.list[[1]][2,i+1]
}
if(sim.version == "diff-effects-2"){
  rr.mat2 <- as.matrix( rr.list[[1]][c(2,1,3),-1] )
  rel.MSE.mat.diff2 <- matrix( 0, nrow = 3, ncol = S )
  for(i in 1:S)  rel.MSE.mat.diff2[,i] <- mse.rte.list[[1]][,i+1] / mse.rte.list[[1]][2,i+1]
}
if(sim.version == "diff-effects-3"){
  rr.mat3 <- as.matrix( rr.list[[1]][c(2,1,3),-1] )
  rel.MSE.mat.diff3 <- matrix( 0, nrow = 3, ncol = S )
  for(i in 1:S)  rel.MSE.mat.diff3[,i] <- mse.rte.list[[1]][,i+1] / mse.rte.list[[1]][2,i+1]
}


par( mar = c(2.6, 3.6, 1, 1) )
par( mfrow = c(3,2) )

# Rejection rates plot - sim.version = "diff-effects-1"
barplot( c(rr.mat1), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0), 6) ),
         col = rep(c("gray50", "gray70", "gray90"), 6) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16, 19.5), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Global", "Region 1", "Region 2", "Region 3", "Region 4", "Region 5" ) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16, 19.5), tick = FALSE, line = .3, cex.axis = .75,
      labels = c( "", expression(paste( "(", gamma[1], " = 0.017)" )),
                  expression(paste( "(", gamma[2], " = 0.026)" )),
                  expression(paste( "(", gamma[3], " = 0.034)" )),
                  expression(paste( "(", gamma[4], " = 0.043)" )),
                  expression(paste( "(", gamma[5], " = 0.051)" )) ) )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = .9 )
mtext( "Rejection Rate", side = 2, line = 2.3, cex = .7 )
#text(.1, .95, "A", cex = 1.2)
legend( x = 5.5, y = 1.01, bty = "n", legend = c("FELM", "BMA", "BHM"), horiz = TRUE,
        fill = c("gray50","gray70", "gray90"), cex = .9 )
box()

# Relative MSE plot - sim.version = "diff-effects-1"
barplot( c(rel.MSE.mat.diff1[-2,]), xlab = "", ylab = "", ylim = c(0,1.2),
         yaxt = "n", space = rep(c(.5,0), 5), col = rep(c("gray60", "gray90"), 5) )
axis( 1, at = c(1.5, 4, 6.5, 9, 11.5), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Region 1", "Region 2", "Region 3", "Region 4", "Region 5" ) )
axis( 1, at = c(1.5, 4, 6.5, 9, 11.5), tick = FALSE, line = .3, cex.axis = .75,
      labels = c( expression(paste( "(", gamma[1], " = 0.017)" )),
                  expression(paste( "(", gamma[2], " = 0.026)" )),
                  expression(paste( "(", gamma[3], " = 0.034)" )),
                  expression(paste( "(", gamma[4], " = 0.043)" )),
                  expression(paste( "(", gamma[5], " = 0.051)" )) ) )
axis( 2, at = seq(0, 1.2, by = 0.2), las = 1, cex.axis = .9 )
abline( h = 1, lty = 3, lwd = 1.5 )
mtext( "Relative MSE", side = 2, line = 2.3, cex = .7 )
#text(.3, .95*1.2, "B", cex = 1.2)
legend( x = .01, y = 1.22, bty = "n", legend = c("BMA", "BHM"), horiz = TRUE,
        fill = c("gray60", "gray90"), cex = .9 )
box()

# Rejection rates plot - sim.version = "diff-effects-2"
barplot( c(rr.mat2), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0), 6) ),
         col = rep(c("gray50", "gray70", "gray90"), 6) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16, 19.5), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Global", "Region 1", "Region 2", "Region 3", "Region 4", "Region 5" ) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16, 19.5), tick = FALSE, line = .3, cex.axis = .75,
      labels = c( "", expression(paste( "(", gamma[1], " = 0.017)" )),
                  expression(paste( "(", gamma[2], " = 0.017)" )),
                  expression(paste( "(", gamma[3], " = 0.017)" )),
                  expression(paste( "(", gamma[4], " = 0.034)" )),
                  expression(paste( "(", gamma[5], " = 0.034)" )) ) )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = .9 )
mtext( "Rejection Rate", side = 2, line = 2.3, cex = .7 )
#text(.1, .95, "C", cex = 1.2)
legend( x = 5.5, y = 1.01, bty = "n", legend = c("FELM", "BMA", "BHM"), horiz = TRUE,
        fill = c("gray50","gray70", "gray90"), cex = .9 )
box()

# Relative MSE plot - sim.version = "diff-effects-2"
barplot( c(rel.MSE.mat.diff2[-2,]), xlab = "", ylab = "", ylim = c(0,1.2),
         yaxt = "n", space = rep(c(.5,0), 5), col = rep(c("gray60", "gray90"), 5) )
axis( 1, at = c(1.5, 4, 6.5, 9, 11.5), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Region 1", "Region 2", "Region 3", "Region 4", "Region 5" ) )
axis( 1, at = c(1.5, 4, 6.5, 9, 11.5), tick = FALSE, line = .3, cex.axis = .75,
      labels = c( expression(paste( "(", gamma[1], " = 0.017)" )),
                  expression(paste( "(", gamma[2], " = 0.017)" )),
                  expression(paste( "(", gamma[3], " = 0.017)" )),
                  expression(paste( "(", gamma[4], " = 0.034)" )),
                  expression(paste( "(", gamma[5], " = 0.034)" )) ) )
axis( 2, at = seq(0, 1.2, by = 0.2), las = 1, cex.axis = .9 )
abline( h = 1, lty = 3, lwd = 1.5 )
mtext( "Relative MSE", side = 2, line = 2.3, cex = .7 )
#text(.3, .95*1.2, "D", cex = 1.2)
legend( x = .01, y = 1.22, bty = "n", legend = c("BMA", "BHM"), horiz = TRUE,
        fill = c("gray60", "gray90"), cex = .9 )
box()

# Rejection rates plot - sim.version = "diff-effects-3"
barplot( c(rr.mat3), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0), 6) ),
         col = rep(c("gray50", "gray70", "gray90"), 6) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16, 19.5), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Global", "Region 1", "Region 2", "Region 3", "Region 4", "Region 5" ) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16, 19.5), tick = FALSE, line = .3, cex.axis = .75,
      labels = c( "", expression(paste( "(", gamma[1], " = 0.017)" )),
                  expression(paste( "(", gamma[2], " = 0.034)" )),
                  expression(paste( "(", gamma[3], " = 0.034)" )),
                  expression(paste( "(", gamma[4], " = 0.051)" )),
                  expression(paste( "(", gamma[5], " = 0.051)" )) ) )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = .9 )
mtext( "Rejection Rate", side = 2, line = 2.3, cex = .7 )
#text(.1, .95, "E", cex = 1.2)
legend( x = 5.5, y = 1.03, bty = "n", legend = c("FELM", "BMA", "BHM"), horiz = TRUE,
        fill = c("gray50","gray70", "gray90"), cex = .9 )
box()

# Relative MSE plot - sim.version = "diff-effects-3"
barplot( c(rel.MSE.mat.diff3[-2,]), xlab = "", ylab = "", ylim = c(0,1.2),
         yaxt = "n", space = rep(c(.5,0), 5), col = rep(c("gray60", "gray90"), 5) )
axis( 1, at = c(1.5, 4, 6.5, 9, 11.5), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Region 1", "Region 2", "Region 3", "Region 4", "Region 5" ) )
axis( 1, at = c(1.5, 4, 6.5, 9, 11.5), tick = FALSE, line = .3, cex.axis = .75,
      labels = c( expression(paste( "(", gamma[1], " = 0.017)" )),
                  expression(paste( "(", gamma[2], " = 0.034)" )),
                  expression(paste( "(", gamma[3], " = 0.034)" )),
                  expression(paste( "(", gamma[4], " = 0.051)" )),
                  expression(paste( "(", gamma[5], " = 0.051)" )) ) )
axis( 2, at = seq(0, 1.2, by = 0.2), las = 1, cex.axis = .9 )
abline( h = 1, lty = 3, lwd = 1.5 )
mtext( "Relative MSE", side = 2, line = 2.3, cex = .7 )
#text(.3, .95*1.2, "F", cex = 1.2)
legend( x = .01, y = 1.22, bty = "n", legend = c("BMA", "BHM"), horiz = TRUE,
        fill = c("gray60", "gray90"), cex = .9 )
box()

par( mfrow = c(1,1) )





### See file "compile_global_consistency_results.R" to compile results for global consistency measure





### Compile results for local consistency measures

## Local consistency using MHLW ratio (not recommended approach)

# Median local consistency probabilities by scenario
round( apply( lc.MHLW.list[[1]], 2, median ), 3 )
round( apply( lc.MHLW.list[[2]], 2, median ), 3 )
round( apply( lc.MHLW.list[[3]], 2, median ), 3 )
round( apply( lc.MHLW.list[[4]], 2, median ), 3 )
round( apply( lc.MHLW.list[[5]], 2, median ), 3 )
round( apply( lc.MHLW.list[[6]], 2, median ), 3 )


## Local consistency using absolute difference from leave-one-out global effect

# Median local consistency probabilities by scenario
# (Technically the average of the median probabilities for each scenario,
# averaged across the result files per scenario)
round( lc.loo.list[[1]], 3 )
round( lc.loo.list[[2]], 3 )
round( lc.loo.list[[3]], 3 )
round( lc.loo.list[[4]], 3 )
round( lc.loo.list[[5]], 3 )
round( lc.loo.list[[6]], 3 )




### Pairwise consistency

# Pairwise consistency matrices by scenario for all values of epsilon.star
round( pc.list[[1]], 3 )
round( pc.list[[2]], 3 )
round( pc.list[[3]], 3 )
round( pc.list[[4]], 3 )
round( pc.list[[5]], 3 )
round( pc.list[[6]], 3 )
