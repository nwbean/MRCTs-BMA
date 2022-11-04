##################################################################################################
# COMPILE RESULTS FROM SIMULATIONS FOR SENSITIVITY ANALYSES
##################################################################################################



### Set working directory
setwd("/Users/nathanbean/Documents/Dissertation/Project 3/Project3_Simulations/cluster-results")


### Simulation details

# Simulation version
sim.version <- "SA-model-priors"
S <- 4      # number of regions

# Total number of result files from simulation
tot.num.results <- 2500                               # use for version "SA-model-priors"
num.scenarios <- 5                                    # number of unique scenarios
num.per.scen <- tot.num.results / num.scenarios       # number of result files per scenario



### Combine results within each scenario for rejection rate, bias, and MSE

## Determine number of rows in the saved csv files for rejection rate, bias, and MSE
num.rows <- 12

## Rejection rate results
rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", 1, ".csv", sep = "")
rr.scen.mat <- read.csv(rr.file, header = TRUE)
rr.scen.mat[1:num.rows, 2:(S+2)] <- 0           # matrix to store rejection rates for scenario

## MSE (regional treatment effects) results
mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", 1, ".csv", sep = "")
mse.rte.scen.mat <- read.csv(mse.rte.file, header = TRUE)
mse.rte.scen.mat[1:num.rows, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario

rr.list <- list()                 # list to store rr.scen.mat for each scenario
mse.rte.list <- list()            # list to store mse.rte.scen.mat for each scenario

which.scen <- 1                   # scenario indicator
for(j in 1:tot.num.results){
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen)
  rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", j, ".csv", sep = "")
  rr.part.mat <- read.csv(rr.file, header = TRUE)
  rr.scen.mat[1:num.rows, 2:(S+2)] <- rr.scen.mat[1:num.rows, 2:(S+2)] +
    rr.part.mat[1:num.rows, 2:(S+2)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen)
  mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", j, ".csv", sep = "")
  mse.rte.part.mat <- read.csv(mse.rte.file, header = TRUE)
  mse.rte.scen.mat[1:num.rows, 2:(S+1)] <- mse.rte.scen.mat[1:num.rows, 2:(S+1)] +
    mse.rte.part.mat[1:num.rows, 2:(S+1)]
  
  # After iterating through num.per.scen files, store results in lists and change scenario indicator
  if(j %% num.per.scen == 0){
    
    # Average results
    rr.scen.mat[1:num.rows, 2:(S+2)] <- rr.scen.mat[1:num.rows, 2:(S+2)] / num.per.scen
    mse.rte.scen.mat[1:num.rows, 2:(S+1)] <- mse.rte.scen.mat[1:num.rows, 2:(S+1)] / num.per.scen
    
    # Store results in list
    scen.lab <- paste("scenario_", which.scen, sep = "")
    rr.list[[scen.lab]] <- rr.scen.mat                 # save scenario-specific matrix of rr's
    mse.rte.list[[scen.lab]] <- mse.rte.scen.mat       # save scenario-specific matrix of mse.rte's
    
    # Reset matrices and prepare for next scenario
    which.scen <- which.scen + 1                       # increase scenario indicator
    rr.scen.mat[1:num.rows, 2:(S+2)] <- 0
    mse.rte.scen.mat[1:num.rows, 2:(S+1)] <- 0
    
  }
  
}


# Global rejection rates
grr.CPHM <- c( rr.list[[1]][1,2], rr.list[[2]][1,2], rr.list[[3]][1,2],
               rr.list[[4]][1,2], rr.list[[5]][1,2] )
grr.BMA.S <- c( rr.list[[1]][2,2], rr.list[[2]][2,2], rr.list[[3]][2,2],
                rr.list[[4]][2,2], rr.list[[5]][2,2] )
grr.BMA.JM <- matrix(0, nrow = num.rows-3, ncol = num.scenarios)
for(i in 4:num.rows){
  grr.BMA.JM[i-3,] <- c( rr.list[[1]][i,2], rr.list[[2]][i,2], rr.list[[3]][i,2],
                         rr.list[[4]][i,2], rr.list[[5]][i,2] )
}
rownames(grr.BMA.JM) <- rr.list[[1]][4:num.rows,1]

# True positive rates
tpr.CPHM <- c( rowMeans(rr.list[[1]][1,3:6]), rowMeans(rr.list[[2]][1,4:6]),
               rowMeans(rr.list[[3]][1,5:6]), rr.list[[4]][1,6], NULL )
tpr.BMA.S <- c( rowMeans(rr.list[[1]][2,3:6]), rowMeans(rr.list[[2]][2,4:6]),
                rowMeans(rr.list[[3]][2,5:6]), rr.list[[4]][2,6], NULL )
tpr.BMA.JM <- matrix(0, nrow = num.rows-3, ncol = num.scenarios-1)
for(i in 3:num.rows){
  tpr.BMA.JM[i-3,] <- c( rowMeans(rr.list[[1]][i,3:6]), rowMeans(rr.list[[2]][i,4:6]),
                         rowMeans(rr.list[[3]][i,5:6]), rr.list[[4]][i,6], NULL )
}
rownames(tpr.BMA.JM) <- rr.list[[1]][4:num.rows,1]

# False positive rates
fpr.CPHM <- c( NULL, rr.list[[2]][1,3], rowMeans(rr.list[[3]][1,3:4]),
               rowMeans(rr.list[[4]][1,3:5]), rowMeans(rr.list[[5]][1,3:6]) )
fpr.BMA.S <- c( NULL, rr.list[[2]][2,3], rowMeans(rr.list[[3]][2,3:4]),
                rowMeans(rr.list[[4]][2,3:5]), rowMeans(rr.list[[5]][2,3:6]) )
fpr.BMA.JM <- matrix(0, nrow = num.rows-3, ncol = num.scenarios-1)
for(i in 4:num.rows){
  fpr.BMA.JM[i-3,] <- c( NULL, rr.list[[2]][i,3], rowMeans(rr.list[[3]][i,3:4]),
                         rowMeans(rr.list[[4]][i,3:5]), rowMeans(rr.list[[5]][i,3:6]) )
}
rownames(fpr.BMA.JM) <- rr.list[[1]][4:num.rows,1]



### Combine regional treatment effect biases according to null and alternative regions for each scenario

## Relative MSE for regional treatment effects with CPHM as reference

rel.MSE.mat.alt <- matrix(0, nrow = S, ncol = num.rows-1)
rel.MSE.mat.null <- matrix(0, nrow = S, ncol = num.rows-1)

# Scenario 1: 0 nulls
# alternative regions
rel.MSE.mat.alt[1,] <- rowMeans(mse.rte.list[[1]][-3,2:5]) / rowMeans(mse.rte.list[[1]][1,2:5])

# Scenario 2: 1 null
# alternative regions
rel.MSE.mat.alt[2,] <- rowMeans(mse.rte.list[[2]][-3,3:5]) / rowMeans(mse.rte.list[[2]][1,3:5])
# null regions
rel.MSE.mat.null[1,] <- mse.rte.list[[2]][-3,2] / mse.rte.list[[2]][1,2]

# Scenario 3: 2 nulls
# alternative regions
rel.MSE.mat.alt[3,] <- rowMeans(mse.rte.list[[3]][-3,4:5]) / rowMeans(mse.rte.list[[3]][1,4:5])
# null regions
rel.MSE.mat.null[2,] <- rowMeans(mse.rte.list[[3]][-3,2:3]) / rowMeans(mse.rte.list[[3]][1,2:3])

# Scenario 4: 3 nulls
# alternative regions
rel.MSE.mat.alt[4,] <- mse.rte.list[[4]][-3,5] / mse.rte.list[[4]][1,5]
# null regions
rel.MSE.mat.null[3,] <- rowMeans(mse.rte.list[[4]][-3,2:4]) / rowMeans(mse.rte.list[[4]][1,2:4])

# Scenario 5: 4 nulls
# null regions
rel.MSE.mat.null[4,] <- rowMeans(mse.rte.list[[5]][-3,2:5]) / rowMeans(mse.rte.list[[5]][1,2:5])





##################################################################################################
# Plot for version "SA-mod-priors"
##################################################################################################

### Plot rejection rates and relative MSE with scenario (number null regions) as x-axis

# Names for legend
names.SA <- c( "CPHM", "BMA-S",
               expression(paste("BMA-JM: ", a[X], " = -1, ", a[Y], " = -1", sep = "")),
               expression(paste("BMA-JM: ", a[X], " = -1, ", a[Y], " = 0", sep = "")),
               expression(paste("BMA-JM: ", a[X], " = -1, ", a[Y], " = 1", sep = "")),
               expression(paste("BMA-JM: ", a[X], " = 0, ", a[Y], " = -1", sep = "")),
               expression(paste("BMA-JM: ", a[X], " = 0, ", a[Y], " = 0", sep = "")),
               expression(paste("BMA-JM: ", a[X], " = 0, ", a[Y], " = 1", sep = "")),
               expression(paste("BMA-JM: ", a[X], " = 1, ", a[Y], " = -1", sep = "")),
               expression(paste("BMA-JM: ", a[X], " = 1, ", a[Y], " = 0", sep = "")),
               expression(paste("BMA-JM: ", a[X], " = 1, ", a[Y], " = 1", sep = "")) )

# Colors and textures of bars
colors.SA <- c("red3", "mediumorchid3", "cornflowerblue", "cornflowerblue", "cornflowerblue",
               "blue", "blue", "blue", "darkblue", "darkblue", "darkblue")
densities.SA <- c(-1,-1,25,-1,50,25,-1,50,25,-1,50)
angles.SA <- 45

# Create PDF file
SA.plot.file <- paste( "/Users/nathanbean/Downloads/", sim.version, ".pdf", sep = "" )
pdf(SA.plot.file, height = 12, width = 10)

# Adjust margins and layout
par( mar = c(4.1, 4.1, 1, 1) )
par( mfrow = c(3,2) )
layout( mat=matrix(c(1,1,2,3,4,5), nrow = 3, byrow = TRUE) )

# Global rejection rate plot
grr.mat <- matrix( rbind(grr.CPHM, grr.BMA.S, grr.BMA.JM), byrow = FALSE, nrow = 1 )
barplot( c(grr.mat), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(1.5,0,0,0,0,0,0,0,0,0,0), 5) ),
         col = rep(colors.SA, 6), density = rep(densities.SA, 6), angle = angles.SA )
axis( 1, at = c(7, 19.5, 32, 44.5, 57), labels = 0:4, tick = FALSE, line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "Global Rejection Rate", side = 2, line = 2.8, cex = 1.1 )
abline( h = .025, lty = 3, lwd = 2 )
text(.1, .95, "A", cex = 2)
legend( "topright", bty = "n", legend = names.SA, fill = colors.SA,
        density = densities.SA, angle = angles.SA, cex = 1.1 )
box()

# Relative MSE plot for alternative regions
MSE.alt.lim <- 1.4
barplot( c(t(rel.MSE.mat.alt[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim),
         yaxt = "n", space = rep(c(1.5,0,0,0,0,0,0,0,0,0), 4),
         col = rep(colors.SA[-1], 4), density = rep(densities.SA[-1], 4), angle = angles.SA )
axis( 1, at = c(6.5, 18, 29.5, 41), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE (alt. regions)", side = 2, line = 2.8, cex = 1.1 )
text(1, .95*MSE.alt.lim, "B", cex = 2)
legend( "topright", bty = "n", legend = names.SA[-1], fill = colors.SA[-1],
        density = densities.SA[-1], angle = angles.SA, ncol = 2, cex = 1.1 )
box()

# Relative MSE plot for null regions
MSE.null.lim <- 1.4
barplot( c(t(rel.MSE.mat.null[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.null.lim),
         yaxt = "n", space = rep(c(1.5,0,0,0,0,0,0,0,0,0), 4), col = rep(colors.SA[-1], 4),
         density = rep(densities.SA[-1], 4), angle = angles.SA )
axis( 1, at = c(6.5, 18, 29.5, 41), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.null.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE (null regions)", side = 2, line = 2.8, cex = 1.1 )
abline( h = 1, lty = 3, lwd = 2 )
text(1, .95*MSE.null.lim, "C", cex = 2)
legend( "topright", bty = "n", legend = names.SA[-1], fill = colors.SA[-1],
        density = densities.SA[-1], angle = angles.SA, ncol = 2, cex = 1.1 )
box()

# True positive rate plot
TPR.lim <- .7
tpr.mat <- matrix( rbind(tpr.CPHM, tpr.BMA.S, tpr.BMA.JM), byrow = FALSE, nrow = 1 )
barplot( c(tpr.mat), xlab = "", ylab = "", ylim = c(0,TPR.lim),
         yaxt = "n", space = c( rep(c(1.5,0,0,0,0,0,0,0,0,0,0), 4) ),
         col = rep(colors.SA, 6), density = rep(densities.SA, 6), angle = angles.SA )
axis( 1, at = c(7, 19.5, 32, 44.5), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, TPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "True Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(1, .95*TPR.lim, "D", cex = 2)
legend( "topright", bty = "n", legend = names.SA, fill = colors.SA,
        density = densities.SA, angle = angles.SA, cex = 1.1 )
box()

# False positive rate plot
FPR.lim <- 0.3
fpr.mat <- matrix( rbind(fpr.CPHM, fpr.BMA.S, fpr.BMA.JM), byrow = FALSE, nrow = 1 )
barplot( c(fpr.mat), xlab = "", ylab = "", ylim = c(0,FPR.lim),
         yaxt = "n", space = c( rep(c(1.5,0,0,0,0,0,0,0,0,0,0), 4) ),
         col = rep(colors.SA, 6), density = rep(densities.SA, 6), angle = angles.SA )
axis( 1, at = c(7, 19.5, 32, 44.5), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, FPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.8, cex = 1.1 )
mtext( "False Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(1, .95*FPR.lim, "E", cex = 2)
abline( h = .025, lty = 3, lwd = 2 )
legend( "topright", bty = "n", legend = names.SA, fill = colors.SA,
        density = densities.SA, angle = angles.SA, cex = 1.1 )
box()

layout( mat=matrix(1, nrow = 1, ncol = 1) )

# Save PDF
dev.off()

