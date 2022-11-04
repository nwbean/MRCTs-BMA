##################################################################################################
# COMPILE RESULTS (E.G., REJECTION RATES, MSE) FOR SIMULATIONS IN WHICH ALTERNATIVE/NULL
# REGIONS ARE HALF THE SIZE OF NULL/ALTERNATIVE REGIONS
##################################################################################################



### Set working directory
setwd("/Users/nathanbean/Documents/Dissertation/Project 3/Project3_Simulations/cluster-results")


### Simulation details

# Simulation version
sim.version1 <- "equal-samp-alpha0p5"    # use for files 1-500 and 2001-2500
sim.version2 <- "equal-samp-alt-half"
#sim.version2 <- "equal-samp-null-half"
S <- 4      # number of regions

# Total number of result files from simulation
tot.num.results <- 2500
num.scenarios <- 5                                   # number of unique scenarios
num.per.scen <- tot.num.results / num.scenarios      # number of result files per scenario



### Combine results within each scenario for rejection rate, bias, and MSE

## Rejection rate results
rr.file <- paste("./", sim.version1, "/rejection_rates/rr_", sim.version1, "-", 1, ".csv", sep = "")
rr.scen.mat <- read.csv(rr.file, header = TRUE)
rr.scen.mat[1:4, 2:(S+2)] <- 0           # matrix to store rejection rates for scenario

## MSE (regional treatment effects) results
mse.rte.file <- paste("./", sim.version1, "/mse/rte_mse_", sim.version1, "-", 1, ".csv", sep = "")
mse.rte.scen.mat <- read.csv(mse.rte.file, header = TRUE)
mse.rte.scen.mat[1:4, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario


rr.list <- list()                 # list to store rr.scen.mat for each scenario
mse.rte.list <- list()            # list to store mse.rte.scen.mat for each scenario


which.scen <- 1                   # scenario indicator
for(j in 1:tot.num.results){
  
  if(which.scen %in% c(1,5)){
    sim.version <- sim.version1
  } else{
    sim.version <- sim.version2
  }
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen)
  rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", j, ".csv", sep = "")
  rr.part.mat <- read.csv(rr.file, header = TRUE)
  rr.scen.mat[1:4, 2:(S+2)] <- rr.scen.mat[1:4, 2:(S+2)] + rr.part.mat[1:4, 2:(S+2)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen)
  mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", j, ".csv", sep = "")
  mse.rte.part.mat <- read.csv(mse.rte.file, header = TRUE)
  mse.rte.scen.mat[1:4, 2:(S+1)] <- mse.rte.scen.mat[1:4, 2:(S+1)] + mse.rte.part.mat[1:4, 2:(S+1)]
  
  if(j %% num.per.scen == 0){
    
    # Average results
    rr.scen.mat[1:4, 2:(S+2)] <- rr.scen.mat[1:4, 2:(S+2)] / num.per.scen
    mse.rte.scen.mat[1:4, 2:(S+1)] <- mse.rte.scen.mat[1:4, 2:(S+1)] / num.per.scen
    
    # Store results in list
    scen.lab <- paste("scenario_", which.scen, sep = "")
    rr.list[[scen.lab]] <- rr.scen.mat                 # save scenario-specific matrix of rr's
    mse.rte.list[[scen.lab]] <- mse.rte.scen.mat       # save scenario-specific matrix of mse.rte's
    
    # Reset matrices and prepare for next scenario
    which.scen <- which.scen + 1                       # increase scenario indicator
    rr.scen.mat[1:4, 2:(S+2)] <- 0
    mse.rte.scen.mat[1:4, 2:(S+1)] <- 0
    
  }
  
}



##################################################################################################
# Results for the following versions:
#   "equal-samp-alt-half"
#   "equal-samp-null-half"
##################################################################################################


### Combine rejection rates according to FPR and TPR for each scenario

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



### Relative MSE for regional treatment effects with CPHM as reference

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
five.plot.file <- paste( "/Users/nathanbean/Downloads/", sim.version2, ".pdf", sep = "" )
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
barplot( c(t(rel.MSE.mat.alt[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim),
         yaxt = "n", space = rep(c(.5,0), 4), col = rep(c("gray60", "gray90"), 4) )
axis( 1, at = c(1.5, 4, 6.5, 9), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE", side = 2, line = 3, cex = 1.1 )
text(.5, .95*MSE.alt.lim, "B", cex = 2)
if(sim.version2 == "equal-samp-alt-half"){
  legend( "topright", bty = "n", fill = c("gray60", "gray90"), cex = 1.1,
          legend = c("BMA-JM (alt. regions)", "BMA-S (alt. regions)") )
} else{
  legend( "top", bty = "n", fill = c("gray60", "gray90"), cex = 1.1,
          legend = c("BMA-JM (alt. regions)", "BMA-S (alt. regions)") )
}
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
if(sim.version2 == "equal-samp-alt-half"){
  FPR.lim <- 0.12
} else{
  FPR.lim <- 0.18
}
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

