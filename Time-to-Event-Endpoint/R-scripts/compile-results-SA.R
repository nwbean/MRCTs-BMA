##################################################################################################
# COMPILE RESULTS FROM SIMULATIONS FOR SENSITIVITY ANALYSES
##################################################################################################



### Set working directory
setwd("/Users/nathanbean/Documents/Dissertation/Project 2/Project2_Simulations/cluster-results")


### Simulation details

# Simulation version
sim.version <- "n9340-equal-samp-SA-alpha"
#sim.version <- "n9340-equal-samp-SA-mu0-Sig0"
#sim.version <- "n9340-equal-samp-SA-K"
S <- 4      # number of regions

# Total number of result files from simulation
tot.num.results <- 100                               # use for all sensitivity analysis (SA) versions
num.scenarios <- 5                                   # number of unique scenarios for most versions
num.per.scen <- tot.num.results / num.scenarios       # number of result files per scenario



### Combine results within each scenario for rejection rate, bias, and MSE

## Determine number of rows in the saved csv files for rejection rate, bias, and MSE
if( sim.version == "n9340-equal-samp-SA-alpha" ){
  num.rows <- 11
} else if( sim.version == "n9340-equal-samp-SA-mu0-Sig0" ){
  num.rows <- 6
} else if( sim.version == "n9340-equal-samp-SA-K" ){
  num.rows <- 6
}

## Rejection rate results
rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", 1, ".csv", sep = "")
rr.scen.mat <- read.csv(rr.file, header = TRUE)
rr.scen.mat[1:num.rows, 2:(S+2)] <- 0           # matrix to store rejection rates for scenario

## Bias (regional treatment effects) results
bias.rte.file <- paste("./", sim.version, "/bias/rte_bias_", sim.version, "-", 1, ".csv", sep = "")
bias.rte.scen.mat <- read.csv(bias.rte.file, header = TRUE)
bias.rte.scen.mat[1:num.rows, 2:(S+1)] <- 0     # matrix to store rte biases for scenario

## MSE (regional treatment effects) results
mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", 1, ".csv", sep = "")
mse.rte.scen.mat <- read.csv(mse.rte.file, header = TRUE)
mse.rte.scen.mat[1:num.rows, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario

rr.list <- list()                 # list to store rr.scen.mat for each scenario
bias.rte.list <- list()           # list to store bias.rte.scen.mat for each scenario
mse.rte.list <- list()            # list to store mse.rte.scen.mat for each scenario

which.scen <- 1                   # scenario indicator
for(j in 1:tot.num.results){
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen)
  rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", j, ".csv", sep = "")
  rr.part.mat <- read.csv(rr.file, header = TRUE)
  rr.scen.mat[1:num.rows, 2:(S+2)] <- rr.scen.mat[1:num.rows, 2:(S+2)] +
    rr.part.mat[1:num.rows, 2:(S+2)]
  
  # Matrix of bias (rte) with part of simulation results (1 out of num.per.scen)
  bias.rte.file <- paste("./", sim.version, "/bias/rte_bias_", sim.version, "-", j, ".csv", sep = "")
  bias.rte.part.mat <- read.csv(bias.rte.file, header = TRUE)
  bias.rte.scen.mat[1:num.rows, 2:(S+1)] <- bias.rte.scen.mat[1:num.rows, 2:(S+1)] +
    bias.rte.part.mat[1:num.rows, 2:(S+1)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen)
  mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", j, ".csv", sep = "")
  mse.rte.part.mat <- read.csv(mse.rte.file, header = TRUE)
  mse.rte.scen.mat[1:num.rows, 2:(S+1)] <- mse.rte.scen.mat[1:num.rows, 2:(S+1)] +
    mse.rte.part.mat[1:num.rows, 2:(S+1)]
  
  # After iterating through num.per.scen files, store results in lists and change scenario indicator
  if(j %% num.per.scen == 0){
    
    # Average results
    rr.scen.mat[1:num.rows, 2:(S+2)] <- rr.scen.mat[1:num.rows, 2:(S+2)] / num.per.scen
    bias.rte.scen.mat[1:num.rows, 2:(S+1)] <- bias.rte.scen.mat[1:num.rows, 2:(S+1)] / num.per.scen
    mse.rte.scen.mat[1:num.rows, 2:(S+1)] <- mse.rte.scen.mat[1:num.rows, 2:(S+1)] / num.per.scen
    
    # Store results in list
    scen.lab <- paste("scenario_", which.scen, sep = "")
    rr.list[[scen.lab]] <- rr.scen.mat                 # save scenario-specific matrix of rr's
    bias.rte.list[[scen.lab]] <- bias.rte.scen.mat     # save scenario-specific matrix of bias.rte's
    mse.rte.list[[scen.lab]] <- mse.rte.scen.mat       # save scenario-specific matrix of mse.rte's
    
    # Reset matrices and prepare for next scenario
    which.scen <- which.scen + 1                       # increase scenario indicator
    rr.scen.mat[1:num.rows, 2:(S+2)] <- 0
    bias.rte.scen.mat[1:num.rows, 2:(S+1)] <- 0
    mse.rte.scen.mat[1:num.rows, 2:(S+1)] <- 0
    
  }
  
}


# Global rejection rates
grr.CPHM <- c( rr.list[[1]][1,2], rr.list[[2]][1,2], rr.list[[3]][1,2],
               rr.list[[4]][1,2], rr.list[[5]][1,2] )
grr.BHM <- c( rr.list[[1]][2,2], rr.list[[2]][2,2], rr.list[[3]][2,2],
              rr.list[[4]][2,2], rr.list[[5]][2,2] )
grr.BMA <- matrix(0, nrow = num.rows-2, ncol = num.scenarios)
for(i in 3:num.rows){
  grr.BMA[i-2,] <- c( rr.list[[1]][i,2], rr.list[[2]][i,2], rr.list[[3]][i,2],
                      rr.list[[4]][i,2], rr.list[[5]][i,2] )
}
rownames(grr.BMA) <- rr.list[[1]][3:num.rows,1]

# True positive rates
tpr.CPHM <- c( rowMeans(rr.list[[1]][1,3:6]), rowMeans(rr.list[[2]][1,4:6]),
               rowMeans(rr.list[[3]][1,5:6]), rr.list[[4]][1,6], NULL )
tpr.BHM <- c( rowMeans(rr.list[[1]][2,3:6]), rowMeans(rr.list[[2]][2,4:6]),
              rowMeans(rr.list[[3]][2,5:6]), rr.list[[4]][2,6], NULL )
tpr.BMA <- matrix(0, nrow = num.rows-2, ncol = num.scenarios-1)
for(i in 3:num.rows){
  tpr.BMA[i-2,] <- c( rowMeans(rr.list[[1]][i,3:6]), rowMeans(rr.list[[2]][i,4:6]),
                      rowMeans(rr.list[[3]][i,5:6]), rr.list[[4]][i,6], NULL )
}
rownames(tpr.BMA) <- rr.list[[1]][3:num.rows,1]

# False positive rates
fpr.CPHM <- c( NULL, rr.list[[2]][1,3], rowMeans(rr.list[[3]][1,3:4]),
               rowMeans(rr.list[[4]][1,3:5]), rowMeans(rr.list[[5]][1,3:6]) )
fpr.BHM <- c( NULL, rr.list[[2]][2,3], rowMeans(rr.list[[3]][2,3:4]),
              rowMeans(rr.list[[4]][2,3:5]), rowMeans(rr.list[[5]][2,3:6]) )
fpr.BMA <- matrix(0, nrow = num.rows-2, ncol = num.scenarios-1)
for(i in 3:num.rows){
  fpr.BMA[i-2,] <- c( NULL, rr.list[[2]][i,3], rowMeans(rr.list[[3]][i,3:4]),
                      rowMeans(rr.list[[4]][i,3:5]), rowMeans(rr.list[[5]][i,3:6]) )
}
rownames(fpr.BMA) <- rr.list[[1]][3:num.rows,1]



### Combine regional treatment effect biases according to null and alternative regions for each scenario

## Relative MSE for regional treatment effects with CPHM as reference

rel.MSE.mat.alt <- matrix(0, nrow = S, ncol = num.rows)
rel.MSE.mat.null <- matrix(0, nrow = S, ncol = num.rows)

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





##################################################################################################
# Plot for version "n9340-equal-samp-SA-alpha"
##################################################################################################

### Plot rejection rates and relative MSE with scenario (number null regions) as x-axis

# Names for legend
names.SA1 <- c( "CPHM", "BHM",
                expression(paste("BMA: ", alpha[0], " = -5", sep = "")),
                expression(paste("BMA: ", alpha[0], " = -2", sep = "")),
                expression(paste("BMA: ", alpha[0], " = -1", sep = "")),
                expression(paste("BMA: ", alpha[0], " = -0.5", sep = "")),
                expression(paste("BMA: ", alpha[0], " = 0", sep = "")),
                expression(paste("BMA: ", alpha[0], " = 0.5", sep = "")),
                expression(paste("BMA: ", alpha[0], " = 1", sep = "")),
                expression(paste("BMA: ", alpha[0], " = 2", sep = "")),
                expression(paste("BMA: ", alpha[0], " = 5", sep = "")) )

# Colors and textures of bars
colors.SA1 <- c("mediumorchid3", "red3", "darkblue", "blue", "deepskyblue", "steelblue3",
                "cornflowerblue", "steelblue3", "deepskyblue", "blue", "darkblue")
densities.SA1 <- c(-1,-1,30,30,30,30,-1,-1,-1,-1,-1)
angles.SA1 <- 45

# Create PDF file
SA1.plot.file <- paste( "/Users/nathanbean/Downloads/", sim.version, ".pdf", sep = "" )
pdf(SA1.plot.file, height = 12, width = 10)

# Adjust margins and layout
par( mar = c(4.1, 4.1, 1, 1) )
par( mfrow = c(3,2) )
layout( mat=matrix(c(1,1,2,3,4,5), nrow = 3, byrow = TRUE) )

# Global rejection rate plot
grr.mat <- matrix( rbind(grr.CPHM, grr.BHM, grr.BMA), byrow = FALSE, nrow = 1 )
barplot( c(grr.mat), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(1.5,0,0,0,0,0,0,0,0,0,0), 5) ),
         col = rep(colors.SA1, 6), density = rep(densities.SA1, 6), angle = angles.SA1 )
axis( 1, at = c(7, 19.5, 32, 44.5, 57), labels = 0:4, tick = FALSE, line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "Global Rejection Rate", side = 2, line = 2.8, cex = 1.1 )
abline( h = .025, lty = 3, lwd = 2 )
text(.1, .95, "A", cex = 2)
legend( "topright", bty = "n", legend = names.SA1, fill = colors.SA1,
        density = densities.SA1, angle = angles.SA1, cex = 1.1 )
box()

# Relative MSE plot for alternative regions
MSE.alt.lim <- 1.4
barplot( c(t(rel.MSE.mat.alt[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim),
         yaxt = "n", space = rep(c(1.5,0,0,0,0,0,0,0,0,0), 4),
         col = rep(colors.SA1[-1], 4), density = rep(densities.SA1[-1], 4), angle = angles.SA1 )
axis( 1, at = c(6.5, 18, 29.5, 41), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE (alt. regions)", side = 2, line = 2.8, cex = 1.1 )
text(1, .95*MSE.alt.lim, "B", cex = 2)
legend( "topright", bty = "n", legend = names.SA1[-1], fill = colors.SA1[-1],
        density = densities.SA1[-1], angle = angles.SA1, ncol = 3, cex = 1.1 )
box()

# Relative MSE plot for null regions
MSE.null.lim <- 1.4
barplot( c(t(rel.MSE.mat.null[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.null.lim),
         yaxt = "n", space = rep(c(1.5,0,0,0,0,0,0,0,0,0), 4), col = rep(colors.SA1[-1], 4),
         density = rep(densities.SA1[-1], 4), angle = angles.SA1 )
axis( 1, at = c(6.5, 18, 29.5, 41), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.null.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE (null regions)", side = 2, line = 2.8, cex = 1.1 )
abline( h = 1, lty = 3, lwd = 2 )
text(1, .95*MSE.null.lim, "C", cex = 2)
legend( "topright", bty = "n", legend = names.SA1[-1], fill = colors.SA1[-1],
        density = densities.SA1[-1], angle = angles.SA1, ncol = 3, cex = 1.1 )
box()

# True positive rate plot
TPR.lim <- .8
tpr.mat <- matrix( rbind(tpr.CPHM, tpr.BHM, tpr.BMA), byrow = FALSE, nrow = 1 )
barplot( c(tpr.mat), xlab = "", ylab = "", ylim = c(0,TPR.lim),
         yaxt = "n", space = c( rep(c(1.5,0,0,0,0,0,0,0,0,0,0), 4) ),
         col = rep(colors.SA1, 6), density = rep(densities.SA1, 6), angle = angles.SA1 )
axis( 1, at = c(7, 19.5, 32, 44.5), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, TPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "True Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(1, .95*TPR.lim, "D", cex = 2)
legend( "topright", bty = "n", legend = names.SA1, fill = colors.SA1,
        density = densities.SA1, angle = angles.SA1, cex = 1.1 )
box()

# False positive rate plot
FPR.lim <- 0.5
fpr.mat <- matrix( rbind(fpr.CPHM, fpr.BHM, fpr.BMA), byrow = FALSE, nrow = 1 )
barplot( c(fpr.mat), xlab = "", ylab = "", ylim = c(0,FPR.lim),
         yaxt = "n", space = c( rep(c(1.5,0,0,0,0,0,0,0,0,0,0), 4) ),
         col = rep(colors.SA1, 6), density = rep(densities.SA1, 6), angle = angles.SA1 )
axis( 1, at = c(7, 19.5, 32, 44.5), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, FPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.8, cex = 1.1 )
mtext( "False Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(1, .95*FPR.lim, "E", cex = 2)
abline( h = .025, lty = 3, lwd = 2 )
legend( "topright", bty = "n", legend = names.SA1, fill = colors.SA1,
        density = densities.SA1, angle = angles.SA1, cex = 1.1 )
box()

layout( mat=matrix(1, nrow = 1, ncol = 1) )

# Save PDF
dev.off()





##################################################################################################
# Plot for version "n9340-equal-samp-SA-mu0-Sig0"
##################################################################################################

### Plot rejection rates and relative MSE with scenario (number null regions) as x-axis

# Names for legend
names.SA2 <- c( "CPHM", "BHM",
                expression(paste("BMA: exp(", mu[0], ") = 0.7", sep = "")),
                expression(paste("BMA: exp(", mu[0], ") = 1.05", sep = "")),
                expression(paste("BMA: exp(", mu[0], ") = 1.3", sep = "")),
                expression(paste("BMA: exp(", mu[0], ") = 1.5", sep = "")) )

# Colors and textures of bars
colors.SA2 <- c("mediumorchid3", "red3", "cornflowerblue", "steelblue3", "deepskyblue", "darkblue")

# Create PDF file
SA2.plot.file <- paste( "/Users/nathanbean/Downloads/", sim.version, ".pdf", sep = "" )
pdf(SA2.plot.file, height = 12, width = 10)

# Adjust margins and layout
par( mar = c(4.1, 4.1, 1, 1) )
par( mfrow = c(3,2) )
layout( mat=matrix(c(1,1,2,3,4,5), nrow = 3, byrow = TRUE) )

# Global rejection rate plot
grr.mat <- matrix( rbind(grr.CPHM, grr.BHM, grr.BMA), byrow = FALSE, nrow = 1 )
barplot( c(grr.mat), xlab = "", ylab = "", ylim = c(0,1), yaxt = "n",
         space = c( rep(c(1.5,0,0,0,0,0), 5) ), col = rep(colors.SA2, 6) )
axis( 1, at = c(4.5, 12, 19.5, 27, 34.5), labels = 0:4, tick = FALSE, line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "Global Rejection Rate", side = 2, line = 2.8, cex = 1.1 )
abline( h = .025, lty = 3, lwd = 2 )
text(.6, .95, "A", cex = 2)
legend( "topright", bty = "n", legend = names.SA2, fill = colors.SA2, cex = 1.1 )
box()

# Relative MSE plot for alternative regions
MSE.alt.lim <- 1.2
barplot( c(t(rel.MSE.mat.alt[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim), yaxt = "n",
         space = rep(c(1.5,0,0,0,0), 4), col = rep(colors.SA2[-1], 4) )
axis( 1, at = c(4, 10.5, 17, 23.5), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE (alt. regions)", side = 2, line = 2.8, cex = 1.1 )
text(1.25, .95*MSE.alt.lim, "B", cex = 2)
legend( "topright", bty = "n", legend = names.SA2[-1], fill = colors.SA2[-1],
        ncol = 2, cex = 1.1 )
box()

# Relative MSE plot for null regions
MSE.null.lim <- 1.2
barplot( c(t(rel.MSE.mat.null[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.null.lim),
         yaxt = "n", space = rep(c(1.5,0,0,0,0), 4), col = rep(colors.SA2[-1], 4) )
axis( 1, at = c(4, 10.5, 17, 23.5), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.null.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE (null regions)", side = 2, line = 2.8, cex = 1.1 )
abline( h = 1, lty = 3, lwd = 2 )
text(1.25, .95*MSE.null.lim, "C", cex = 2)
legend( "topright", bty = "n", legend = names.SA2[-1], fill = colors.SA2[-1],
        ncol = 2, cex = 1.1 )
box()

# True positive rate plot
TPR.lim <- .6
tpr.mat <- matrix( rbind(tpr.CPHM, tpr.BHM, tpr.BMA), byrow = FALSE, nrow = 1 )
barplot( c(tpr.mat), xlab = "", ylab = "", ylim = c(0,TPR.lim), yaxt = "n",
         space = c( rep(c(1.5,0,0,0,0,0), 4) ), col = rep(colors.SA2, 6) )
axis( 1, at = c(4.5, 12, 19.5, 27), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, TPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "True Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(1.25, .95*TPR.lim, "D", cex = 2)
legend( "topright", bty = "n", legend = names.SA2, fill = colors.SA2, cex = 1.1 )
box()

# False positive rate plot
FPR.lim <- 0.16
fpr.mat <- matrix( rbind(fpr.CPHM, fpr.BHM, fpr.BMA), byrow = FALSE, nrow = 1 )
barplot( c(fpr.mat), xlab = "", ylab = "", ylim = c(0,FPR.lim), yaxt = "n",
         space = c( rep(c(1.5,0,0,0,0,0), 4) ), col = rep(colors.SA2, 6) )
axis( 1, at = c(4.5, 12, 19.5, 27), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, FPR.lim, by = 0.02), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.8, cex = 1.1 )
mtext( "False Positive Rate", side = 2, line = 3.1, cex = 1.1 )
text(1.25, .95*FPR.lim, "E", cex = 2)
abline( h = .025, lty = 3, lwd = 2 )
legend( "topright", bty = "n", legend = names.SA2, fill = colors.SA2, cex = 1.1 )
box()

layout( mat=matrix(1, nrow = 1, ncol = 1) )

# Save PDF
dev.off()





##################################################################################################
# Plot for version "n9340-equal-samp-SA-K"
##################################################################################################

### Plot rejection rates and relative MSE with scenario (number null regions) as x-axis

# Names for legend
names.SA3 <- c( "CPHM", "BHM", "BMA: K = 4", "BMA: K = 8", "BMA: K = 12", "BMA: K = 16" )

# Colors and textures of bars
colors.SA3 <- c("mediumorchid3", "red3", "cornflowerblue", "steelblue3", "deepskyblue", "darkblue")

# Create PDF file
SA3.plot.file <- paste( "/Users/nathanbean/Downloads/", sim.version, ".pdf", sep = "" )
pdf(SA3.plot.file, height = 12, width = 10)

# Adjust margins and layout
par( mar = c(4.1, 4.1, 1, 1) )
par( mfrow = c(3,2) )
layout( mat=matrix(c(1,1,2,3,4,5), nrow = 3, byrow = TRUE) )

# Global rejection rate plot
grr.mat <- matrix( rbind(grr.CPHM, grr.BHM, grr.BMA), byrow = FALSE, nrow = 1 )
barplot( c(grr.mat), xlab = "", ylab = "", ylim = c(0,1), yaxt = "n",
         space = c( rep(c(1.5,0,0,0,0,0), 5) ), col = rep(colors.SA3, 6) )
axis( 1, at = c(4.5, 12, 19.5, 27, 34.5), labels = 0:4, tick = FALSE, line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "Global Rejection Rate", side = 2, line = 2.8, cex = 1.1 )
abline( h = .025, lty = 3, lwd = 2 )
text(.6, .95, "A", cex = 2)
legend( "topright", bty = "n", legend = names.SA3, fill = colors.SA3, cex = 1.1 )
box()

# Relative MSE plot for alternative regions
MSE.alt.lim <- 1.2
barplot( c(t(rel.MSE.mat.alt[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim), yaxt = "n",
         space = rep(c(1.5,0,0,0,0), 4), col = rep(colors.SA3[-1], 4) )
axis( 1, at = c(4, 10.5, 17, 23.5), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE (alt. regions)", side = 2, line = 2.8, cex = 1.1 )
text(1.25, .95*MSE.alt.lim, "B", cex = 2)
legend( "topright", bty = "n", legend = names.SA3[-1], fill = colors.SA3[-1],
        ncol = 2, cex = 1.1 )
box()

# Relative MSE plot for null regions
MSE.null.lim <- 1.2
barplot( c(t(rel.MSE.mat.null[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.null.lim),
         yaxt = "n", space = rep(c(1.5,0,0,0,0), 4), col = rep(colors.SA3[-1], 4) )
axis( 1, at = c(4, 10.5, 17, 23.5), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.null.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE (null regions)", side = 2, line = 2.8, cex = 1.1 )
abline( h = 1, lty = 3, lwd = 2 )
text(1.25, .95*MSE.null.lim, "C", cex = 2)
legend( "topright", bty = "n", legend = names.SA3[-1], fill = colors.SA3[-1],
        ncol = 2, cex = 1.1 )
box()

# True positive rate plot
TPR.lim <- .5
tpr.mat <- matrix( rbind(tpr.CPHM, tpr.BHM, tpr.BMA), byrow = FALSE, nrow = 1 )
barplot( c(tpr.mat), xlab = "", ylab = "", ylim = c(0,TPR.lim), yaxt = "n",
         space = c( rep(c(1.5,0,0,0,0,0), 4) ), col = rep(colors.SA3, 6) )
axis( 1, at = c(4.5, 12, 19.5, 27), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, TPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "True Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(1.25, .95*TPR.lim, "D", cex = 2)
legend( "topright", bty = "n", legend = names.SA3, fill = colors.SA3, cex = 1.1 )
box()

# False positive rate plot
FPR.lim <- 0.14
fpr.mat <- matrix( rbind(fpr.CPHM, fpr.BHM, fpr.BMA), byrow = FALSE, nrow = 1 )
barplot( c(fpr.mat), xlab = "", ylab = "", ylim = c(0,FPR.lim), yaxt = "n",
         space = c( rep(c(1.5,0,0,0,0,0), 4) ), col = rep(colors.SA3, 6) )
axis( 1, at = c(4.5, 12, 19.5, 27), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, FPR.lim, by = 0.02), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.8, cex = 1.1 )
mtext( "False Positive Rate", side = 2, line = 3.1, cex = 1.1 )
text(1.25, .95*FPR.lim, "E", cex = 2)
abline( h = .025, lty = 3, lwd = 2 )
legend( "topright", bty = "n", legend = names.SA3, fill = colors.SA3, cex = 1.1 )
box()

layout( mat=matrix(1, nrow = 1, ncol = 1) )

# Save PDF
dev.off()
