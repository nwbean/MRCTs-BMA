##################################################################################################
# COMPILE RESULTS (E.G., REJECTION RATES, MSE) FROM TWO SIMULATION SCENARIOS
# COMPARING DIFFERENT BHM PRIORS ON THE HIERARCHICAL STANDARD DEVIATION OR PRECISION
##################################################################################################



### Set working directory
setwd("/Users/nathanbean/Documents/Dissertation/Project 2/Project2_Simulations/cluster-results")


### Simulation details

# Simulation version
sim.version <- "n9340-equal-samp-BHM-priors"
#sim.version <- "n9340-original-BHM-priors"
S <- 4      # number of regions

# Total number of result files from simulation
tot.num.results <- 50                                # use for version "n9340-equal-samp"
#tot.num.results <- 40                                # use for version "n9340-original"
num.scenarios <- 5                                   # number of unique scenarios "n9340-equal-samp"
#num.scenarios <- 2                                   # number of unique scenarios "n9340-original"
num.per.scen <- tot.num.results / num.scenarios       # number of result files per scenario



### Combine results within each scenario for rejection rate, bias, and MSE

## Rejection rate results
rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", 1, ".csv", sep = "")
rr.scen.mat <- read.csv(rr.file, header = TRUE)
rr.scen.mat[1:5, 2:(S+2)] <- 0           # matrix to store rejection rates for scenario

## Bias (regional treatment effects) results
bias.rte.file <- paste("./", sim.version, "/bias/rte_bias_", sim.version, "-", 1, ".csv", sep = "")
bias.rte.scen.mat <- read.csv(bias.rte.file, header = TRUE)
bias.rte.scen.mat[1:5, 2:(S+1)] <- 0     # matrix to store rte biases for scenario

## MSE (regional treatment effects) results
mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", 1, ".csv", sep = "")
mse.rte.scen.mat <- read.csv(mse.rte.file, header = TRUE)
mse.rte.scen.mat[1:5, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario

## Hierarchical variance summary statistics
hier.var.file <- paste("./", sim.version, "/BHM-hier-var/ct_", sim.version, "-", 1, ".csv", sep = "")
hier.var.scen.mat <- read.csv(hier.var.file, header = TRUE)
hier.var.scen.mat[1:3, 1:4] <- 0      # matrix to store summary stats for scenario
rownames(hier.var.scen.mat) <- mse.rte.scen.mat[3:5,1]

## Gelman-Rubin multivariate potential scale reduction factor
mpsrf.file <- paste("./", sim.version, "/gelman-rubin/multi_", sim.version, "-", 1, ".csv", sep = "")
mpsrf.scen.mat <- read.csv(mpsrf.file, header = TRUE)
mpsrf.scen.mat[1:3, 1:6] <- 0      # matrix to store summary stats for scenario
rownames(mpsrf.scen.mat) <- mse.rte.scen.mat[3:5,1]

## Gelman-Rubin potential scale reduction factors - point estimates
pe.psrf.file <- paste("./", sim.version, "/gelman-rubin/pe_", sim.version, "-", 1, ".csv", sep = "")
pe.psrf.scen.mat <- read.csv(pe.psrf.file, header = TRUE)
pe.psrf.scen.mat[1:3, 1:(S + 2)] <- 0      # matrix to store summary stats for scenario
rownames(pe.psrf.scen.mat) <- mse.rte.scen.mat[3:5,1]

## Gelman-Rubin potential scale reduction factors - upper bound of 95% credible interval
uci.psrf.file <- paste("./", sim.version, "/gelman-rubin/uci_", sim.version, "-", 1, ".csv", sep = "")
uci.psrf.scen.mat <- read.csv(uci.psrf.file, header = TRUE)
uci.psrf.scen.mat[1:3, 1:(S + 2)] <- 0      # matrix to store summary stats for scenario
rownames(uci.psrf.scen.mat) <- mse.rte.scen.mat[3:5,1]

rr.list <- list()                 # list to store rr.scen.mat for each scenario
bias.rte.list <- list()           # list to store bias.rte.scen.mat for each scenario
mse.rte.list <- list()            # list to store mse.rte.scen.mat for each scenario
hier.var.list <- list()           # list to store hier.var.scen.mat for each scenario
mpsrf.list <- list()              # list to store mpsrf.scen.mat for each scenario
pe.psrf.list <- list()            # list to store pe.psrf.scen.mat for each scenario
uci.psrf.list <- list()           # list to store uci.psrf.scen.mat for each scenario


which.scen <- 1                   # scenario indicator
for(j in 1:tot.num.results){
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen)
  rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", j, ".csv", sep = "")
  rr.part.mat <- read.csv(rr.file, header = TRUE)
  rr.scen.mat[1:5, 2:(S+2)] <- rr.scen.mat[1:5, 2:(S+2)] + rr.part.mat[1:5, 2:(S+2)]
  
  # Matrix of bias (rte) with part of simulation results (1 out of num.per.scen)
  bias.rte.file <- paste("./", sim.version, "/bias/rte_bias_", sim.version, "-", j, ".csv", sep = "")
  bias.rte.part.mat <- read.csv(bias.rte.file, header = TRUE)
  bias.rte.scen.mat[1:5, 2:(S+1)] <- bias.rte.scen.mat[1:5, 2:(S+1)] + bias.rte.part.mat[1:5, 2:(S+1)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen)
  mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", j, ".csv", sep = "")
  mse.rte.part.mat <- read.csv(mse.rte.file, header = TRUE)
  mse.rte.scen.mat[1:5, 2:(S+1)] <- mse.rte.scen.mat[1:5, 2:(S+1)] + mse.rte.part.mat[1:5, 2:(S+1)]
  
  # Matrix of hier. var. summary stats with part of simulation results (1 out of num.per.scen)
  hier.var.file <- paste("./", sim.version, "/BHM-hier-var/ct_", sim.version, "-", j, ".csv", sep = "")
  hier.var.part.mat <- read.csv(hier.var.file, header = TRUE)
  hier.var.scen.mat <- hier.var.scen.mat + hier.var.part.mat
  
  # Matrix of multivariate PSRF with part of simulation results (1 out of num.per.scen)
  mpsrf.file <- paste("./", sim.version, "/gelman-rubin/multi_", sim.version, "-", j, ".csv", sep = "")
  mpsrf.part.mat <- read.csv(mpsrf.file, header = TRUE)
  mpsrf.scen.mat <- mpsrf.scen.mat + mpsrf.part.mat
  
  # Matrix of PSRF point estimates with part of simulation results (1 out of num.per.scen)
  pe.psrf.file <- paste("./", sim.version, "/gelman-rubin/pe_", sim.version, "-", j, ".csv", sep = "")
  pe.psrf.part.mat <- read.csv(pe.psrf.file, header = TRUE)
  pe.psrf.scen.mat <- pe.psrf.scen.mat + pe.psrf.part.mat
  
  # Matrix of PSRF upper cred. int. bounds with part of simulation results (1 out of num.per.scen)
  uci.psrf.file <- paste("./", sim.version, "/gelman-rubin/uci_", sim.version, "-", j, ".csv", sep = "")
  uci.psrf.part.mat <- read.csv(uci.psrf.file, header = TRUE)
  uci.psrf.scen.mat <- uci.psrf.scen.mat + uci.psrf.part.mat
  
  # After iterating through num.per.scen files, store results in lists and change scenario indicator
  if(j %% num.per.scen == 0){
    
    # Average results
    rr.scen.mat[1:5, 2:(S+2)] <- rr.scen.mat[1:5, 2:(S+2)] / num.per.scen
    bias.rte.scen.mat[1:5, 2:(S+1)] <- bias.rte.scen.mat[1:5, 2:(S+1)] / num.per.scen
    mse.rte.scen.mat[1:5, 2:(S+1)] <- mse.rte.scen.mat[1:5, 2:(S+1)] / num.per.scen
    hier.var.scen.mat <- hier.var.scen.mat / num.per.scen
    mpsrf.scen.mat <- mpsrf.scen.mat / num.per.scen
    pe.psrf.scen.mat <- pe.psrf.scen.mat / num.per.scen
    uci.psrf.scen.mat <- uci.psrf.scen.mat / num.per.scen
    
    # Store results in list
    scen.lab <- paste("scenario_", which.scen, sep = "")
    rr.list[[scen.lab]] <- rr.scen.mat                 # save scenario-specific matrix of rr's
    bias.rte.list[[scen.lab]] <- bias.rte.scen.mat     # save scenario-specific matrix of bias.rte's
    mse.rte.list[[scen.lab]] <- mse.rte.scen.mat       # save scenario-specific matrix of mse.rte's
    hier.var.list[[scen.lab]] <- hier.var.scen.mat     # save scenario-specific matrix of hier. var.
    mpsrf.list[[scen.lab]] <- mpsrf.scen.mat           # save scenario-specific matrix of multi PSRF
    pe.psrf.list[[scen.lab]] <- pe.psrf.scen.mat       # save scenario-specific matrix of PSRF p.e.
    uci.psrf.list[[scen.lab]] <- uci.psrf.scen.mat     # save scenario-specific matrix of PSRF u.c.i.
    
    # Reset matrices and prepare for next scenario
    which.scen <- which.scen + 1                       # increase scenario indicator
    rr.scen.mat[1:5, 2:(S+2)] <- 0
    bias.rte.scen.mat[1:5, 2:(S+1)] <- 0
    mse.rte.scen.mat[1:5, 2:(S+1)] <- 0
    hier.var.scen.mat[1:3,] <- 0
    mpsrf.scen.mat[1:3,] <- 0
    pe.psrf.scen.mat[1:3,] <- 0
    uci.psrf.scen.mat[1:3,] <- 0
    
  }
  
}



##################################################################################################
# Results for "n9340-equal-samp-BHM-priors"
##################################################################################################


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
grr.BMA <- c( rr.list[[1]][2,2], rr.list[[2]][2,2], rr.list[[3]][2,2],
              rr.list[[4]][2,2], rr.list[[5]][2,2] )
grr.BHM.Gam <- c( rr.list[[1]][3,2], rr.list[[2]][3,2], rr.list[[3]][3,2],
                  rr.list[[4]][3,2], rr.list[[5]][3,2] )
grr.BHM.Unif <- c( rr.list[[1]][4,2], rr.list[[2]][4,2], rr.list[[3]][4,2],
                   rr.list[[4]][4,2], rr.list[[5]][4,2] )
grr.BHM.HC <- c( rr.list[[1]][5,2], rr.list[[2]][5,2], rr.list[[3]][5,2],
                 rr.list[[4]][5,2], rr.list[[5]][5,2] )

# True positive rates
tpr.CPHM <- c( rowMeans(rr.list[[1]][1,3:6]), rowMeans(rr.list[[2]][1,4:6]),
               rowMeans(rr.list[[3]][1,5:6]), rr.list[[4]][1,6], NULL )
tpr.BMA <- c( rowMeans(rr.list[[1]][2,3:6]), rowMeans(rr.list[[2]][2,4:6]),
              rowMeans(rr.list[[3]][2,5:6]), rr.list[[4]][2,6], NULL )
tpr.BHM.Gam <- c( rowMeans(rr.list[[1]][3,3:6]), rowMeans(rr.list[[2]][3,4:6]),
                  rowMeans(rr.list[[3]][3,5:6]), rr.list[[4]][3,6], NULL )
tpr.BHM.Unif <- c( rowMeans(rr.list[[1]][4,3:6]), rowMeans(rr.list[[2]][4,4:6]),
                   rowMeans(rr.list[[3]][4,5:6]), rr.list[[4]][4,6], NULL )
tpr.BHM.HC <- c( rowMeans(rr.list[[1]][5,3:6]), rowMeans(rr.list[[2]][5,4:6]),
                 rowMeans(rr.list[[3]][5,5:6]), rr.list[[4]][5,6], NULL )

# False positive rates
fpr.CPHM <- c( NULL, rr.list[[2]][1,3], rowMeans(rr.list[[3]][1,3:4]),
               rowMeans(rr.list[[4]][1,3:5]), rowMeans(rr.list[[5]][1,3:6]) )
fpr.BMA <- c( NULL, rr.list[[2]][2,3], rowMeans(rr.list[[3]][2,3:4]),
              rowMeans(rr.list[[4]][2,3:5]), rowMeans(rr.list[[5]][2,3:6]) )
fpr.BHM.Gam <- c( NULL, rr.list[[2]][3,3], rowMeans(rr.list[[3]][3,3:4]),
                  rowMeans(rr.list[[4]][3,3:5]), rowMeans(rr.list[[5]][3,3:6]) )
fpr.BHM.Unif <- c( NULL, rr.list[[2]][4,3], rowMeans(rr.list[[3]][4,3:4]),
                   rowMeans(rr.list[[4]][4,3:5]), rowMeans(rr.list[[5]][4,3:6]) )
fpr.BHM.HC <- c( NULL, rr.list[[2]][5,3], rowMeans(rr.list[[3]][5,3:4]),
                 rowMeans(rr.list[[4]][5,3:5]), rowMeans(rr.list[[5]][5,3:6]) )



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

rel.MSE.mat.alt <- matrix(0, nrow = S, ncol = 5)
rel.MSE.mat.null <- matrix(0, nrow = S, ncol = 5)

# Scenario 1: 0 nulls
# alternative regions
rel.MSE.mat.alt[1,] <- rowMeans(mse.rte.list[[1]][,2:5]) / rowMeans(mse.rte.list[[1]][1,2:5])

# Scenario 2: 1 null
# alternative regions
rel.MSE.mat.alt[2,] <- rowMeans(mse.rte.list[[2]][,3:5]) / rowMeans(mse.rte.list[[2]][1,3:5])
# null regions
rel.MSE.mat.null[1,] <- mse.rte.list[[2]][,2] / mse.rte.list[[2]][1,2]

# Scenario 3: 2 nulls
# alternative regions
rel.MSE.mat.alt[3,] <- rowMeans(mse.rte.list[[3]][,4:5]) / rowMeans(mse.rte.list[[3]][1,4:5])
# null regions
rel.MSE.mat.null[2,] <- rowMeans(mse.rte.list[[3]][,2:3]) / rowMeans(mse.rte.list[[3]][1,2:3])

# Scenario 4: 3 nulls
# alternative regions
rel.MSE.mat.alt[4,] <- mse.rte.list[[4]][,5] / mse.rte.list[[4]][1,5]
# null regions
rel.MSE.mat.null[3,] <- rowMeans(mse.rte.list[[4]][,2:4]) / rowMeans(mse.rte.list[[4]][1,2:4])

# Scenario 5: 4 nulls
# null regions
rel.MSE.mat.null[4,] <- rowMeans(mse.rte.list[[5]][,2:5]) / rowMeans(mse.rte.list[[5]][1,2:5])



### Plot rejection rates and relative MSE with scenario (number null regions) as x-axis

# Create PDF file
five.plot.file <- paste( "/Users/nathanbean/Downloads/", sim.version, ".pdf", sep = "" )
pdf(five.plot.file, height = 6.46, width = 13)

# Adjust margins and layout
par( mar = c(4.1, 4.1, 1, 1) )
layout( mat=matrix(c(1,1,2,2,3,3,4,5,5,6,6,7), nrow = 2, byrow = TRUE) )
colors.BHM.priors <- c("mediumorchid3", "cornflowerblue", "red3", "indianred1", "firebrick4")

# Global rejection rate plot
grr.mat <- matrix( rbind(grr.CPHM, grr.BMA, grr.BHM.Gam, grr.BHM.Unif, grr.BHM.HC),
                   byrow = FALSE, nrow = 1 )
barplot( c(grr.mat), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0,0,0), 5) ), col = rep(colors.BHM.priors, 6) )
axis( 1, at = c(3, 8.5, 14, 19.5, 25), labels = 0:4, tick = FALSE, line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "Global Rejection Rate", side = 2, line = 2.8, cex = 1.1 )
abline( h = .025, lty = 3, lwd = 2 )
text(.5, .95, "A", cex = 2)
legend( "topright", bty = "n", legend = c("CPHM", "BMA", "BHM-G", "BHM-U", "BHM-HC"),
        fill = colors.BHM.priors, cex = 1.1 )
box()

# Relative MSE plot for alternative regions
MSE.alt.lim <- 1.2
barplot( c(t(rel.MSE.mat.alt[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim),
         yaxt = "n", space = rep(c(.5,0,0,0), 4), col = rep(colors.BHM.priors[-1], 4) )
axis( 1, at = c(2.5, 7, 11.5, 16), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
text(.5, .95*MSE.alt.lim, "B", cex = 2)
legend( "topright", bty = "n", legend = c("BMA (alt. regions)", "BHM-G (alt. regions)",
                                          "BHM-U (alt. regions)", "BHM-HC (alt. regions)"),
        fill = colors.BHM.priors[-1], ncol = 2, cex = 1.1 )
box()

# Relative MSE plot for null regions
MSE.null.lim <- 1.2
barplot( c(t(rel.MSE.mat.null[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.null.lim),
         yaxt = "n", space = rep(c(.5,0,0,0), 4), col = rep(colors.BHM.priors[-1], 4) )
axis( 1, at = c(2.5, 7, 11.5, 16), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.null.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
abline( h = 1, lty = 3, lwd = 2 )
text(.5, .95*MSE.null.lim, "C", cex = 2)
legend( "topright", bty = "n", legend = c("BMA (null regions)", "BHM-G (null regions)",
                                          "BHM-U (null regions)", "BHM-HC (null regions)"),
        fill = colors.BHM.priors[-1], ncol = 2, cex = 1.1 )
box()

# Empty plot
plot.new()

# True positive rate plot
TPR.lim <- .6
tpr.mat <- matrix( rbind(tpr.CPHM, tpr.BMA, tpr.BHM.Gam, tpr.BHM.Unif, tpr.BHM.HC),
                   byrow = FALSE, nrow = 1 )
barplot( c(tpr.mat), xlab = "", ylab = "", ylim = c(0,TPR.lim),
         yaxt = "n", space = c( rep(c(.5,0,0,0,0), 4) ), col = rep(colors.BHM.priors, 4) )
axis( 1, at = c(3, 8.5, 14, 19.5), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, TPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "True Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(.5, .95*TPR.lim, "D", cex = 2)
legend( "topright", bty = "n", legend = c("CPHM", "BMA", "BHM-G", "BHM-U", "BHM-HC"),
        fill = colors.BHM.priors, cex = 1.1 )
box()

# False positive rate plot
FPR.lim <- 0.14
fpr.mat <- matrix( rbind(fpr.CPHM, fpr.BMA, fpr.BHM.Gam, fpr.BHM.Unif, fpr.BHM.HC),
                   byrow = FALSE, nrow = 1 )
barplot( c(fpr.mat), xlab = "", ylab = "", ylim = c(0,FPR.lim),
         yaxt = "n", space = c( rep(c(.5,0,0,0,0), 4) ), col = rep(colors.BHM.priors, 4) )
axis( 1, at = c(3, 8.5, 14, 19.5), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, FPR.lim, by = 0.02), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.8, cex = 1.1 )
mtext( "False Positive Rate", side = 2, line = 3.1, cex = 1.1 )
text(.5, .95*FPR.lim, "E", cex = 2)
abline( h = .025, lty = 3, lwd = 2 )
legend( "topright", bty = "n", legend = c("CPHM", "BMA", "BHM-G", "BHM-U", "BHM-HC"),
        fill = colors.BHM.priors, cex = 1.1 )
box()

# Empty plot
plot.new()

layout( mat=matrix(1, nrow = 1, ncol = 1) )

# Save PDF
dev.off()





##################################################################################################
# Results for version "n9340-original"
##################################################################################################

### Plot results for sim.version = "n9340-original"

# Save rejection rates and relative MSE for same treatment effect scenario
rr.mat1 <- as.matrix( rr.list[[1]][,-1] )
rel.MSE.mat.diff1 <- matrix( 0, nrow = 5, ncol = S )
for(i in 1:S)  rel.MSE.mat.diff1[,i] <- mse.rte.list[[1]][,i+1] / mse.rte.list[[1]][1,i+1]

# Save rejection rates and relative MSE for different treatment effects scenario
rr.mat2 <- as.matrix( rr.list[[2]][,-1] )
rel.MSE.mat.diff2 <- matrix( 0, nrow = 5, ncol = S )
for(i in 1:S)  rel.MSE.mat.diff2[,i] <- mse.rte.list[[2]][,i+1] / mse.rte.list[[2]][1,i+1]

# Create PDF file
four.plot.file <- paste( "/Users/nathanbean/Downloads/", sim.version, ".pdf", sep = "" )
pdf(four.plot.file, height = 6.93, width = 10.73)

par( mar = c(3.6, 3.6, 1, 1) )
par( mfrow = c(2,2) )
colors.BHM.priors <- c("mediumorchid3", "cornflowerblue", "red3", "indianred1", "firebrick4")

# Rejection rates plot - sim.version = "n9340-original" (same treatment effect)
barplot( c(rr.mat1), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0,0,0), 5) ),
         col = rep(colors.BHM.priors, 5) )
axis( 1, at = c(3, 8.5, 14, 19.5, 25), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Global", "Asia", "Europe", "N. America", "Rest of World" ) )
axis( 1, at = c(8.5, 14, 19.5, 25), tick = FALSE, line = .3, cex.axis = .8,
      labels = c("(n = 711,", "(n = 3296,", "(n = 2847,", "(n = 2486,") )
axis( 1, at = c(8.5, 14, 19.5, 25), tick = FALSE, line = 1.2, cex.axis = .8,
      labels = rep("HR = 0.868)", S) )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = .9 )
mtext( "Rejection Rate", side = 2, line = 2.3, cex = .9 )
text(.3, .95, "A", cex = 1.2)
legend( "topright", bty = "n", legend = c("CPHM", "BMA", "BHM-G", "BHM-U", "BHM-HC"),
        fill = colors.BHM.priors, cex = .9 )
box()

# Relative MSE plot - sim.version = "n9340-original" (same treatment effect)
barplot( c(rel.MSE.mat.diff1[-1,]), xlab = "", ylab = "", ylim = c(0,1.2),
         yaxt = "n", space = rep(c(.5,0,0,0), 4), col = rep(colors.BHM.priors[-1], 4) )
axis( 1, at = c(2.5, 7, 11.5, 16), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Asia", "Europe", "North America", "Rest of World" ) )
axis( 1, at = c(2.5, 7, 11.5, 16), tick = FALSE, line = .3, cex.axis = .8,
      labels = c("(n = 711,", "(n = 3296,", "(n = 2847,", "(n = 2486,") )
axis( 1, at = c(2.5, 7, 11.5, 16), tick = FALSE, line = 1.2, cex.axis = .8,
      labels = rep("HR = 0.868)", S) )
axis( 2, at = seq(0, 1.2, by = 0.2), las = 1, cex.axis = .9 )
abline( h = 1, lty = 3, lwd = 1.5 )
mtext( "Relative MSE", side = 2, line = 2.3, cex = .9 )
text(.4, .95*1.2, "B", cex = 1.2)
legend( "topright", bty = "n", legend = c("BMA", "BHM-G", "BHM-U", "BHM-HC"), ncol = 2,
        fill = colors.BHM.priors[-1], cex = .9 )
box()

# Rejection rates plot - sim.version = "n9340-original" (different treatment effect)
barplot( c(rr.mat2), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0,0,0), 5) ),
         col = rep(colors.BHM.priors, 5) )
axis( 1, at = c(3, 8.5, 14, 19.5, 25), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Global", "Asia", "Europe", "N. America", "Rest of World" ) )
axis( 1, at = c(8.5, 14, 19.5, 25), tick = FALSE, line = .3, cex.axis = .8,
      labels = c("(n = 711,", "(n = 3296,", "(n = 2847,", "(n = 2486,") )
axis( 1, at = c(8.5, 14, 19.5, 25), tick = FALSE, line = 1.2, cex.axis = .8,
      labels = c("HR = 0.62)", "HR = 0.82)", "HR = 1.01)", "HR = 0.83)") )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = .9 )
mtext( "Rejection Rate", side = 2, line = 2.3, cex = .9 )
text(.3, .95, "C", cex = 1.2)
legend( "topright", bty = "n", legend = c("CPHM", "BMA", "BHM-G", "BHM-U", "BHM-HC"),
        fill = colors.BHM.priors, cex = .9 )
box()

# Relative MSE plot - sim.version = "n9340-original" (different treatment effect)
barplot( c(rel.MSE.mat.diff2[-1,]), xlab = "", ylab = "", ylim = c(0,1.6),
         yaxt = "n", space = rep(c(.5,0,0,0), 4), col = rep(colors.BHM.priors[-1], 4) )
axis( 1, at = c(2.5, 7, 11.5, 16), tick = FALSE, line = -.7, cex.axis = .9,
      labels = c( "Asia", "Europe", "North America", "Rest of World" ) )
axis( 1, at = c(2.5, 7, 11.5, 16), tick = FALSE, line = .3, cex.axis = .8,
      labels = c("(n = 711,", "(n = 3296,", "(n = 2847,", "(n = 2486,") )
axis( 1, at = c(2.5, 7, 11.5, 16), tick = FALSE, line = 1.2, cex.axis = .8,
      labels = c("HR = 0.62)", "HR = 0.82)", "HR = 1.01)", "HR = 0.83)") )
axis( 2, at = seq(0, 1.6, by = 0.2), las = 1, cex.axis = .9 )
abline( h = 1, lty = 3, lwd = 1.5 )
mtext( "Relative MSE", side = 2, line = 2.3, cex = .9 )
text(.4, .95*1.6, "D", cex = 1.2)
legend( "topright", bty = "n", legend = c("BMA", "BHM-G", "BHM-U", "BHM-HC"),
        fill = colors.BHM.priors[-1], cex = .9 )
box()

par( mfrow = c(1,1) )

# Save PDF
dev.off()

