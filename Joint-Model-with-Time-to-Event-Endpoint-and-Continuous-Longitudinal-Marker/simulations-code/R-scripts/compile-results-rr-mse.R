##################################################################################################
# COMPILE RESULTS (E.G., REJECTION RATES, MSE) FROM MOST SIMULATION SCENARIOS
# (EXCLUDES SENSITIVITY ANALYSES AND GLOBAL CONSISTENCY RESULTS)
##################################################################################################



### Set working directory
setwd("/Users/nathanbean/Documents/Dissertation/Project 3/Project3_Simulations/cluster-results")


### Simulation details

# Simulation version
#sim.version <- "equal-samp-alpha0"
#sim.version <- "equal-samp-alpha0p15"
sim.version <- "equal-samp-alpha0p5"
#sim.version <- "equal-samp-alpha1"
#sim.version <- "original-alpha0"
#sim.version <- "original-alpha0p15"
#sim.version <- "original-alpha0p5"
#sim.version <- "original-alpha1"
S <- 4      # number of regions

# Total number of result files from simulation
tot.num.results <- 2500                              # use for all "equal-samp" versions
#tot.num.results <- 1000                              # use for all "original" versions
num.scenarios <- 5                                   # number of scenarios "equal-samp" versions
#num.scenarios <- 2                                   # number of scenarios "original" versions
num.per.scen <- tot.num.results / num.scenarios      # number of result files per scenario



### Combine results within each scenario for rejection rate, bias, and MSE

## Rejection rate results
rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", 1, ".csv", sep = "")
rr.scen.mat <- read.csv(rr.file, header = TRUE)
rr.scen.mat[1:4, 2:(S+2)] <- 0           # matrix to store rejection rates for scenario

## MSE (regional treatment effects) results
mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", 1, ".csv", sep = "")
mse.rte.scen.mat <- read.csv(mse.rte.file, header = TRUE)
mse.rte.scen.mat[1:4, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario

## Bias (regional treatment effects) results
#bias.rte.file <- paste("./", sim.version, "/bias/rte_bias_", sim.version, "-", 1, ".csv", sep = "")
#bias.rte.scen.mat <- read.csv(bias.rte.file, header = TRUE)
#bias.rte.scen.mat[1:4, 2:(S+1)] <- 0     # matrix to store rte biases for scenario

## MSE (regional treatment effects) results
#var.rte.file <- paste("./", sim.version, "/sd/sd_", sim.version, "-", 1, ".csv", sep = "")
#var.rte.scen.mat <- read.csv(var.rte.file, header = TRUE)
#var.rte.scen.mat[1:4, 1:S] <- 0      # matrix to store rte SDs for scenario
#row.names(var.rte.scen.mat) <- c("CPHM", "BMA-JM", "BMA-S", "BHM")


#gte.scen.mat <- NULL
#rte.scen.mat <- NULL
#pmp.scen.mat <- NULL
#ap.scen.mat <- NULL
#pc.file <- paste("./", sim.version, "/prws_consis/pc_", sim.version, "-", 1, ".csv", sep = "")
#pc.example.mat <- read.csv(pc.file, header = TRUE)
#pc.scen.mat <- matrix( 0, nrow = dim(pc.example.mat)[1], ncol = dim(pc.example.mat)[2] )

#gte.list <- list()                # list to store gte.scen.mat for each scenario
#rte.list <- list()                # list to store rte.scen.mat for each scenario
rr.list <- list()                 # list to store rr.scen.mat for each scenario
mse.rte.list <- list()            # list to store mse.rte.scen.mat for each scenario
#bias.rte.list <- list()           # list to store bias.rte.scen.mat for each scenario
#var.rte.list <- list()            # list to store var.rte.scen.mat for each scenario
#pmp.list <- list()                # list to store pmp.scen.mat for each scenario
#ap.list <- list()                 # list to store ap.scen.mat for each scenario
#lc.PMDA.RR.list <- list()         # list to store lc.PMDA.RR.scen.mat for each scenario
#lc.PMDA.TE.list <- list()         # list to store lc.PMDA.TE.scen.mat for each scenario
#lc.loo.list <- list()             # list to store lc.loo.scen.mat (leave one out) for each scenario
#pc.list <- list()                 # list to store pc.scen.mat for each scenario


which.scen <- 1                   # scenario indicator
for(j in 1:tot.num.results){
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen)
  rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", j, ".csv", sep = "")
  rr.part.mat <- read.csv(rr.file, header = TRUE)
  rr.scen.mat[1:4, 2:(S+2)] <- rr.scen.mat[1:4, 2:(S+2)] + rr.part.mat[1:4, 2:(S+2)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen)
  mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", j, ".csv", sep = "")
  mse.rte.part.mat <- read.csv(mse.rte.file, header = TRUE)
  mse.rte.scen.mat[1:4, 2:(S+1)] <- mse.rte.scen.mat[1:4, 2:(S+1)] + mse.rte.part.mat[1:4, 2:(S+1)]
  
  # Matrix of bias (rte) with part of simulation results (1 out of num.per.scen)
  #bias.rte.file <- paste("./", sim.version, "/bias/rte_bias_", sim.version, "-", j, ".csv", sep = "")
  #bias.rte.part.mat <- read.csv(bias.rte.file, header = TRUE)
  #bias.rte.scen.mat[1:4, 2:(S+1)] <- bias.rte.scen.mat[1:4, 2:(S+1)] + bias.rte.part.mat[1:4, 2:(S+1)]
  
  # Matrix of variances (rte) with part of simulation results (1 out of num.per.scen)
  #var.rte.file <- paste("./", sim.version, "/sd/sd_", sim.version, "-", j, ".csv", sep = "")
  #var.rte.part.mat <- read.csv(var.rte.file, header = TRUE)^2
  #var.rte.scen.mat[1:4, 1:S] <- var.rte.scen.mat[1:4, 1:S] + var.rte.part.mat[1:4, 1:S]
  
  # Matrix of global treatment effects with part of simulation results (1 out of num.per.scen)
  #gte.file <- paste("./", sim.version, "/gte/gte_", sim.version, "-", j, ".csv", sep = "")
  #gte.part.mat <- read.csv(gte.file, header = TRUE)
  #gte.scen.mat <- rbind(gte.scen.mat, gte.part.mat)
  
  # Matrix of regional treatment effects with part of simulation results (1 out of num.per.scen)
  #rte.file <- paste("./", sim.version, "/rte/rte_", sim.version, "-", j, ".csv", sep = "")
  #rte.part.mat <- read.csv(rte.file, header = TRUE)
  #rte.scen.mat <- rbind(rte.scen.mat, rte.part.mat)
  
  # Matrix of PMPs with part of simulation results (1 out of num.per.scen)
  #pmp.file <- paste("./", sim.version, "/pmp/pmp_", sim.version, "-", j, ".csv", sep = "")
  #pmp.part.mat <- read.csv(pmp.file, header = TRUE)
  #pmp.scen.mat <- rbind(pmp.scen.mat, pmp.part.mat[-c(1:3),])
  
  # Matrix of association par. summary stats with part of simulation results (1 out of num.per.scen)
  #ap.file <- paste("./", sim.version, "/assos_par/assos_par1_", sim.version, "-", j, ".csv", sep = "")
  #ap.part.mat <- read.csv(ap.file, header = TRUE)
  #ap.scen.mat <- rbind(ap.scen.mat, ap.part.mat)
  
  # Matrix of local consistency probabilities (PMDA using risk reduction) with part of
  # simulation results (1 out of num.per.scen)
  #lc.PMDA.RR.file <- paste("./", sim.version, "/loc_consis/lc_PMDA_RR_", sim.version, "-", j, ".csv", sep = "")
  #lc.PMDA.RR.part.mat <- read.csv(lc.PMDA.RR.file, header = TRUE)
  #lc.PMDA.RR.scen.mat <- rbind(lc.PMDA.RR.scen.mat, lc.PMDA.RR.part.mat)
  
  # Matrix of local consistency probabilities (PMDA using treatment effect) with part of
  # simulation results (1 out of num.per.scen)
  #lc.PMDA.TE.file <- paste("./", sim.version, "/loc_consis/lc_PMDA_TE_", sim.version, "-", j, ".csv", sep = "")
  #lc.PMDA.TE.part.mat <- read.csv(lc.PMDA.TE.file, header = TRUE)
  #lc.PMDA.TE.scen.mat <- rbind(lc.PMDA.TE.scen.mat, lc.PMDA.TE.part.mat)
  
  # Matrix of local consistency probabilities (LOO) with part of simulation results (1 out of num.per.scen)
  #lc.loo.file <- paste("./", sim.version, "/loc_consis/lc_loo_", sim.version, "-", j, ".csv", sep = "")
  #lc.loo.part.mat <- read.csv(lc.loo.file, header = TRUE)
  #lc.loo.scen.mat <- lc.loo.scen.mat + lc.loo.part.mat
  
  # Matrix of pairwise consistency probabilities with part of simulation results (1 out of num.per.scen)
  #pc.file <- paste("./", sim.version, "/prws_consis/pc_", sim.version, "-", j, ".csv", sep = "")
  #pc.part.mat <- read.csv(pc.file, header = TRUE)
  #pc.scen.mat <- pc.scen.mat + pc.part.mat
  
  # After iterating through num.per.scen files, store results in lists and change scenario indicator
  if(j %% num.per.scen == 0){
    
    # Average results
    rr.scen.mat[1:4, 2:(S+2)] <- rr.scen.mat[1:4, 2:(S+2)] / num.per.scen
    mse.rte.scen.mat[1:4, 2:(S+1)] <- mse.rte.scen.mat[1:4, 2:(S+1)] / num.per.scen
    #bias.rte.scen.mat[1:4, 2:(S+1)] <- bias.rte.scen.mat[1:4, 2:(S+1)] / num.per.scen
    #var.rte.scen.mat[1:4, 1:S] <- var.rte.scen.mat[1:4, 1:S] / num.per.scen
    
    # Store results in list
    scen.lab <- paste("scenario_", which.scen, sep = "")
    rr.list[[scen.lab]] <- rr.scen.mat                 # save scenario-specific matrix of rr's
    mse.rte.list[[scen.lab]] <- mse.rte.scen.mat       # save scenario-specific matrix of mse.rte's
    #bias.rte.list[[scen.lab]] <- bias.rte.scen.mat     # save scenario-specific matrix of bias.rte's
    #var.rte.list[[scen.lab]] <- var.rte.scen.mat       # save scenario-specific matrix of var.rte's
    #gte.list[[scen.lab]] <- gte.scen.mat               # save scenario-specific matrix of gte's
    #rte.list[[scen.lab]] <- rte.scen.mat               # save scenario-specific matrix of rte's
    #pmp.list[[scen.lab]] <- pmp.scen.mat               # save scenario-specific matrix of PMPs
    #ap.list[[scen.lab]] <- ap.scen.mat                 # save scenario-specific matrix of assos. pars
    #lc.PMDA.RR.list[[scen.lab]] <- lc.PMDA.RR.scen.mat # save scenario-specific matrix of cp's (PMDA - RR)
    #lc.PMDA.TE.list[[scen.lab]] <- lc.PMDA.TE.scen.mat # save scenario-specific matrix of cp's (PMDA - TE)
    #lc.loo.list[[scen.lab]] <- lc.loo.scen.mat /       # save scenario-specific matrix of cp's (LOO)
    #  num.per.scen
    
    # Save scenario-specific matrix of pc's (for all values of epsilon.star)
    #pc.list[[scen.lab]] <- pc.scen.mat / num.per.scen
    
    # Reset matrices and prepare for next scenario
    which.scen <- which.scen + 1                       # increase scenario indicator
    rr.scen.mat[1:4, 2:(S+2)] <- 0
    mse.rte.scen.mat[1:4, 2:(S+1)] <- 0
    #bias.rte.scen.mat[1:4, 2:(S+1)] <- 0
    #var.rte.scen.mat[1:4, 1:S] <- 0
    #gte.scen.mat <- NULL
    #rte.scen.mat <- NULL
    #pmp.scen.mat <- NULL
    #ap.scen.mat <- NULL
    #lc.PMDA.RR.scen.mat <- NULL
    #lc.PMDA.TE.scen.mat <- NULL
    #lc.loo.scen.mat <- matrix( 0, nrow = dim(lc.loo.scen.mat)[1], ncol = dim(lc.loo.scen.mat)[2] )
    #pc.scen.mat <- matrix( 0, nrow = dim(pc.scen.mat)[1], ncol = dim(pc.scen.mat)[2] )
    
  }
  
}



##################################################################################################
# Results for the following versions:
#   "equal-samp-alpha0"
#   "equal-samp-alpha0p15"
#   "equal-samp-alpha0p5"
#   "equal-samp-alpha1"
##################################################################################################


### Assess PMPs
round( colMeans(pmp.list[[1]]), 4)
round( colMeans(pmp.list[[2]]), 4)
round( colMeans(pmp.list[[3]]), 4)
round( colMeans(pmp.list[[4]]), 4)
round( colMeans(pmp.list[[5]]), 4)



### Association parameter summaries
round( colMeans(ap.list[[1]]), 4)
round( colMeans(ap.list[[2]]), 4)
round( colMeans(ap.list[[3]]), 4)
round( colMeans(ap.list[[4]]), 4)
round( colMeans(ap.list[[5]]), 4)


### Combine rejection rates according to FPR and TPR for each scenario

# Scenario 1: 0 nulls
round( rr.list[[1]][,2], 3 )                      # Global rejection rates
round( rowMeans(rr.list[[1]][,3:6]), 3 )          # TPR averages

# Scenario 2: 1 null
round( rr.list[[2]][,2], 3 )                      # Global rejection rates
round( rowMeans(rr.list[[2]][,4:6]), 3 )          # TPR averages
round( rr.list[[2]][,3], 3 )                      # FPR averages

# Scenario 3: 2 nulls
round( rr.list[[3]][,2], 3 )                      # Global rejection rates
round( rowMeans(rr.list[[3]][,5:6]), 3 )          # TPR averages
round( rowMeans(rr.list[[3]][,3:4]), 3 )          # FPR averages

# Scenario 4: 3 nulls
round( rr.list[[4]][,2], 3 )                      # Global rejection rates
round( rr.list[[4]][,6], 3 )                      # TPR averages
round( rowMeans(rr.list[[4]][,3:5]), 3 )          # FPR averages

# Scenario 5: 4 nulls
round( rr.list[[5]][,2], 3 )                      # Global rejection rates
round( rowMeans(rr.list[[5]][,3:6]), 3 )          # FPR averages


# Global rejection rates
grr.CPHM <- c( rr.list[[1]][1,2], rr.list[[2]][1,2], rr.list[[3]][1,2],
               rr.list[[4]][1,2], rr.list[[5]][1,2] )
grr.BMA.JM <- c( rr.list[[1]][2,2], rr.list[[2]][2,2], rr.list[[3]][2,2],
                 rr.list[[4]][2,2], rr.list[[5]][2,2] )
grr.BMA.S <- c( rr.list[[1]][3,2], rr.list[[2]][3,2], rr.list[[3]][3,2],
                rr.list[[4]][3,2], rr.list[[5]][3,2] )

# True positive rates
tpr.CPHM <- c( rowMeans(rr.list[[1]][1,3:6]), rowMeans(rr.list[[2]][1,4:6]),
               rowMeans(rr.list[[3]][1,5:6]), rr.list[[4]][1,6], NULL )
tpr.BMA.JM <- c( rowMeans(rr.list[[1]][2,3:6]), rowMeans(rr.list[[2]][2,4:6]),
                 rowMeans(rr.list[[3]][2,5:6]), rr.list[[4]][2,6], NULL )
tpr.BMA.S <- c( rowMeans(rr.list[[1]][3,3:6]), rowMeans(rr.list[[2]][3,4:6]),
                rowMeans(rr.list[[3]][3,5:6]), rr.list[[4]][3,6], NULL )

# False positive rates
fpr.CPHM <- c( NULL, rr.list[[2]][1,3], rowMeans(rr.list[[3]][1,3:4]),
               rowMeans(rr.list[[4]][1,3:5]), rowMeans(rr.list[[5]][1,3:6]) )
fpr.BMA.JM <- c( NULL, rr.list[[2]][2,3], rowMeans(rr.list[[3]][2,3:4]),
                 rowMeans(rr.list[[4]][2,3:5]), rowMeans(rr.list[[5]][2,3:6]) )
fpr.BMA.S <- c( NULL, rr.list[[2]][3,3], rowMeans(rr.list[[3]][3,3:4]),
                rowMeans(rr.list[[4]][3,3:5]), rowMeans(rr.list[[5]][3,3:6]) )



### Posterior means of region-specific trtmt effects for null and alt. regions for each scenario

# Scenario 1: 0 nulls
round( mean(as.matrix(rte.list[[1]])), 7 )          # alternative regions

# Scenario 2: 1 null
round( mean(as.matrix(rte.list[[2]][,2:4])), 7 )    # alternative regions
round( mean(rte.list[[2]][,1]), 6 )                 # null regions

# Scenario 3: 2 nulls
round( mean(as.matrix(rte.list[[3]][,3:4])), 7 )    # alternative regions
round( mean(as.matrix(rte.list[[3]][,1:2])), 7 )    # null regions

# Scenario 4: 3 nulls
round( mean(rte.list[[4]][,4]), 6 )                 # alternative regions
round( mean(as.matrix(rte.list[[4]][,1:3])), 7 )    # null regions

# Scenario 5: 4 nulls
round( mean(as.matrix(rte.list[[5]])), 7 )          # null regions



### Posterior stand. devs of region-specific trtmt effects for null and alt. regions for each scenario

# Scenario 1: 0 nulls
round( sd(as.matrix(rte.list[[1]])), 7 )            # alternative regions

# Scenario 2: 1 null
round( sd(as.matrix(rte.list[[2]][,2:4])), 7 )      # alternative regions
round( sd(rte.list[[2]][,1]), 6 )                   # null regions

# Scenario 3: 2 nulls
round( sd(as.matrix(rte.list[[3]][,3:4])), 7 )      # alternative regions
round( sd(as.matrix(rte.list[[3]][,1:2])), 7 )      # null regions

# Scenario 4: 3 nulls
round( sd(rte.list[[4]][,4]), 6 )                   # alternative regions
round( sd(as.matrix(rte.list[[4]][,1:3])), 7 )      # null regions

# Scenario 5: 4 nulls
round( sd(as.matrix(rte.list[[5]])), 7 )            # null regions



### Combine regional treatment effect biases according to null and alternative regions for each scenario

# Scenario 1: 0 nulls
round( rowMeans(bias.rte.list[[1]][,2:5]), 7 )    # alternative regions

# Scenario 2: 1 null
round( rowMeans(bias.rte.list[[2]][,3:5]), 7 )    # alternative regions
round( bias.rte.list[[2]][,2], 6 )                # null regions

# Scenario 3: 2 nulls
round( rowMeans(bias.rte.list[[3]][,4:5]), 7 )    # alternative regions
round( rowMeans(bias.rte.list[[3]][,2:3]), 7 )    # null regions

# Scenario 4: 3 nulls
round( bias.rte.list[[4]][,5], 6 )                # alternative regions
round( rowMeans(bias.rte.list[[4]][,2:4]), 7 )    # null regions

# Scenario 5: 4 nulls
round( rowMeans(bias.rte.list[[5]][,2:5]), 7 )    # null regions



## Relative MSE for regional treatment effects with CPHM as reference

rel.MSE.mat.alt <- matrix(0, nrow = S, ncol = 3)
rel.MSE.mat.null <- matrix(0, nrow = S, ncol = 3)

# Scenario 1: 0 nulls
# alternative regions
rel.MSE.mat.alt[1,] <- rowMeans(mse.rte.list[[1]][1:3,2:5]) / rowMeans(mse.rte.list[[1]][1,2:5])

# Scenario 2: 1 null
# alternative regions
rel.MSE.mat.alt[2,] <- rowMeans(mse.rte.list[[2]][1:3,3:5]) / rowMeans(mse.rte.list[[2]][1,3:5])
# null regions
rel.MSE.mat.null[1,] <- mse.rte.list[[2]][1:3,2] / mse.rte.list[[2]][1,2]

# Scenario 3: 2 nulls
# alternative regions
rel.MSE.mat.alt[3,] <- rowMeans(mse.rte.list[[3]][1:3,4:5]) / rowMeans(mse.rte.list[[3]][1,4:5])
# null regions
rel.MSE.mat.null[2,] <- rowMeans(mse.rte.list[[3]][1:3,2:3]) / rowMeans(mse.rte.list[[3]][1,2:3])

# Scenario 4: 3 nulls
# alternative regions
rel.MSE.mat.alt[4,] <- mse.rte.list[[4]][1:3,5] / mse.rte.list[[4]][1,5]
# null regions
rel.MSE.mat.null[3,] <- rowMeans(mse.rte.list[[4]][1:3,2:4]) / rowMeans(mse.rte.list[[4]][1,2:4])

# Scenario 5: 4 nulls
# null regions
rel.MSE.mat.null[4,] <- rowMeans(mse.rte.list[[5]][1:3,2:5]) / rowMeans(mse.rte.list[[5]][1,2:5])



### Plot rejection rates and relative MSE with scenario (number null regions) as x-axis

# Create PDF file
five.plot.file <- paste( "/Users/nathanbean/Downloads/", sim.version, ".pdf", sep = "" )
pdf(five.plot.file, height = 6.46, width = 13)

# Adjust margins and layout
par( mar = c(4.1, 4.1, 1, 1) )
layout( mat=matrix(c(1,1,2,2,3,3,4,5,5,6,6,7), nrow = 2, byrow = TRUE) )

# Global rejection rate plot
grr.mat <- matrix( rbind(grr.CPHM, grr.BMA.JM, grr.BMA.S), byrow = FALSE, nrow = 1 )
barplot( c(grr.mat), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0), 5) ),
         col = rep(c("gray50", "gray70", "gray90"), 6) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16), labels = 0:4, tick = FALSE, line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "Global Rejection Rate", side = 2, line = 2.8, cex = 1.1 )
abline( h = .025, lty = 3, lwd = 2 )
text(.5, .95, "A", cex = 2)
legend( "topright", bty = "n", legend = c("CPHM", "BMA-JM", "BMA-S"),
        fill = c("gray50", "gray70", "gray90"), cex = 1.1 )
box()

# Relative MSE plot for alternative regions
MSE.alt.lim <- 1.2
#MSE.alt.lim <- 1.4    # use for "equal-samp-alpha1"
barplot( c(t(rel.MSE.mat.alt[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim),
         yaxt = "n", space = rep(c(.5,0), 4), col = rep(c("gray60", "gray90"), 4) )
axis( 1, at = c(1.5, 4, 6.5, 9), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
text(.5, .95*MSE.alt.lim, "B", cex = 2)
legend( "topright", bty = "n", fill = c("gray60", "gray90"), cex = 1.1,
        legend = c("BMA-JM (alt. regions)", "BMA-S (alt. regions)") )
box()

# Relative MSE plot for null regions
MSE.null.lim <- 1.2
barplot( c(t(rel.MSE.mat.null[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.null.lim),
         yaxt = "n", space = rep(c(.5,0), 4), col = rep(c("gray60", "gray90"), 4) )
axis( 1, at = c(1.5, 4, 6.5, 9), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.null.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
abline( h = 1, lty = 3, lwd = 2 )
text(.5, .95*MSE.null.lim, "C", cex = 2)
legend( "topright", bty = "n", fill = c("gray60", "gray90"), cex = 1.1,
        legend = c("BMA-JM (null regions)", "BMA-S (null regions)") )
box()

# Empty plot
plot.new()

# True positive rate plot
TPR.lim <- .6
tpr.mat <- matrix( rbind(tpr.CPHM, tpr.BMA.JM, tpr.BMA.S), byrow = FALSE, nrow = 1 )
barplot( c(tpr.mat), xlab = "", ylab = "", ylim = c(0,TPR.lim),
         yaxt = "n", space = c( rep(c(.5,0,0), 4) ),
         col = rep(c("gray50", "gray70", "gray90"), 4) )
axis( 1, at = c(2, 5.5, 9, 12.5), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, TPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "True Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(.5, .95*TPR.lim, "D", cex = 2)
legend( "topright", bty = "n", legend = c("CPHM", "BMA-JM", "BMA-S"),
        fill = c("gray50", "gray70", "gray90"), cex = 1.1 )
box()

# False positive rate plot
FPR.lim <- 0.16
fpr.mat <- matrix( rbind(fpr.CPHM, fpr.BMA.JM, fpr.BMA.S), byrow = FALSE, nrow = 1 )
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
legend( "topright", bty = "n", legend = c("CPHM", "BMA-JM", "BMA-S"),
        fill = c("gray50", "gray70", "gray90"), cex = 1.1 )
box()

# Empty plot
plot.new()

layout( mat=matrix(1, nrow = 1, ncol = 1) )

# Save PDF
dev.off()





##################################################################################################
# Results for the following versions:
#   "original-alpha0"
#   "original-alpha0p15"
#   "original-alpha0p5"
#   "original-alpha1"
##################################################################################################

### Plot results for sim.version = "original-alpha__"

# Save rejection rates and relative MSE for same treatment effect scenario
rr.mat1 <- as.matrix( rr.list[[1]][1:3,-1] )
rel.MSE.mat.diff1 <- matrix( 0, nrow = 3, ncol = S )
for(i in 1:S)  rel.MSE.mat.diff1[,i] <- mse.rte.list[[1]][1:3,i+1] / mse.rte.list[[1]][1,i+1]

# Save rejection rates and relative MSE for different treatment effects scenario
rr.mat2 <- as.matrix( rr.list[[2]][1:3,-1] )
rel.MSE.mat.diff2 <- matrix( 0, nrow = 3, ncol = S )
for(i in 1:S)  rel.MSE.mat.diff2[,i] <- mse.rte.list[[2]][1:3,i+1] / mse.rte.list[[2]][1,i+1]

# Create PDF file
four.plot.file <- paste( "/Users/nathanbean/Downloads/", sim.version, ".pdf", sep = "" )
pdf(four.plot.file, height = 6.93, width = 10.73)

par( mar = c(3.6, 3.6, 1, 1) )
par( mfrow = c(2,2) )

# Rejection rates plot - sim.version = "original-alpha__" (same treatment effect)
barplot( c(rr.mat1), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0), 5) ),
         col = rep(c("gray50", "gray70", "gray90"), 5) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Global", "Asia", "Europe", "N. America", "Rest of World" ) )
axis( 1, at = c(5.5, 9, 12.5, 16), tick = FALSE, line = .3, cex.axis = .8,
      labels = c("(n = 711,", "(n = 3296,", "(n = 2847,", "(n = 2486,") )
axis( 1, at = c(5.5, 9, 12.5, 16), tick = FALSE, line = 1.2, cex.axis = .8,
      labels = rep("HR = 0.868)", S) )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = .9 )
mtext( "Rejection Rate", side = 2, line = 2.3, cex = .9 )
text(.3, .95, "A", cex = 1.2)
legend( "topright", bty = "n", legend = c("CPHM", "BMA-JM", "BMA-S"),
        fill = c("gray50", "gray70", "gray90"), cex = .9 )
box()

# Relative MSE plot - sim.version = "original-alpha__" (same treatment effect)
barplot( c(rel.MSE.mat.diff1[-1,]), xlab = "", ylab = "", ylim = c(0,1.2),
         yaxt = "n", space = rep(c(.5,0), 4),
         col = rep(c("gray60", "gray90"), 4) )
axis( 1, at = c(1.5, 4, 6.5, 9), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Asia", "Europe", "North America", "Rest of World" ) )
axis( 1, at = c(1.5, 4, 6.5, 9), tick = FALSE, line = .3, cex.axis = .8,
      labels = c("(n = 711,", "(n = 3296,", "(n = 2847,", "(n = 2486,") )
axis( 1, at = c(1.5, 4, 6.5, 9), tick = FALSE, line = 1.2, cex.axis = .8,
      labels = rep("HR = 0.868)", S) )
axis( 2, at = seq(0, 1.2, by = 0.2), las = 1, cex.axis = .9 )
abline( h = 1, lty = 3, lwd = 1.5 )
mtext( "Relative MSE", side = 2, line = 2.3, cex = .9 )
text(.4, .95*1.2, "B", cex = 1.2)
legend( "topright", bty = "n", legend = c("BMA-JM", "BMA-S"),
        fill = c("gray60", "gray90"), cex = .9 )
box()

# Rejection rates plot - sim.version = "original-alpha__" (different treatment effect)
barplot( c(rr.mat2), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0), 5) ),
         col = rep(c("gray50", "gray70", "gray90"), 5) )
axis( 1, at = c(2, 5.5, 9, 12.5, 16), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Global", "Asia", "Europe", "N. America", "Rest of World" ) )
axis( 1, at = c(5.5, 9, 12.5, 16), tick = FALSE, line = .3, cex.axis = .8,
      labels = c("(n = 711,", "(n = 3296,", "(n = 2847,", "(n = 2486,") )
axis( 1, at = c(5.5, 9, 12.5, 16), tick = FALSE, line = 1.2, cex.axis = .8,
      labels = c("HR = 0.62)", "HR = 0.82)", "HR = 1.01)", "HR = 0.83)") )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = .9 )
mtext( "Rejection Rate", side = 2, line = 2.3, cex = .9 )
text(.3, .95, "C", cex = 1.2)
legend( "topright", bty = "n", legend = c("CPHM", "BMA-JM", "BMA-S"),
        fill = c("gray50", "gray70", "gray90"), cex = .9 )
box()

# Relative MSE plot - sim.version = "original-alpha__" (different treatment effect)
MSE.ylim.diff <- 1.8     # use for "original-alpha0p5"
#MSE.ylim.diff <- 2       # use for "original-alpha2"
barplot( c(rel.MSE.mat.diff2[-1,]), xlab = "", ylab = "", ylim = c(0, MSE.ylim.diff),
         yaxt = "n", space = rep(c(.5,0), 4),
         col = rep(c("gray60", "gray90"), 4) )
axis( 1, at = c(1.5, 4, 6.5, 9), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Asia", "Europe", "North America", "Rest of World" ) )
axis( 1, at = c(1.5, 4, 6.5, 9), tick = FALSE, line = .3, cex.axis = .8,
      labels = c("(n = 711,", "(n = 3296,", "(n = 2847,", "(n = 2486,") )
axis( 1, at = c(1.5, 4, 6.5, 9), tick = FALSE, line = 1.2, cex.axis = .8,
      labels = c("HR = 0.62)", "HR = 0.82)", "HR = 1.01)", "HR = 0.83)") )
axis( 2, at = seq(0, MSE.ylim.diff, by = 0.2), las = 1, cex.axis = .9 )
abline( h = 1, lty = 3, lwd = 1.5 )
mtext( "Relative MSE", side = 2, line = 2.3, cex = .9 )
text(.4, .95*MSE.ylim.diff, "D", cex = 1.2)
legend( "topright", bty = "n", legend = c("BMA-JM", "BMA-S"),
        fill = c("gray60", "gray90"), cex = .9 )
box()

par( mfrow = c(1,1) )

# Save PDF
dev.off()

