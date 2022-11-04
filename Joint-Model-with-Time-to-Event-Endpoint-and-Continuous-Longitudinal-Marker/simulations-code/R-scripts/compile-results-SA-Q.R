##################################################################################################
# COMPILE RESULTS (E.G., REJECTION RATES, MSE) FROM MOST SIMULATION SCENARIOS
# (EXCLUDES SENSITIVITY ANALYSES AND GLOBAL CONSISTENCY RESULTS)
##################################################################################################



### Set working directory
setwd("/Users/nathanbean/Documents/Dissertation/Project 3/Project3_Simulations/cluster-results")


### Simulation details

# Simulation version
sim.version <- "equal-samp-alpha0p5"
#sim.version <- "equal-samp-Q5"
#sim.version <- "equal-samp-Q12"
S <- 4      # number of regions

# Total number of result files from simulation
tot.num.results <- 2500                              # use for all "equal-samp" versions
num.scenarios <- 5                                   # number of scenarios "equal-samp" versions
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

rr.list <- list()                 # list to store rr.scen.mat for each scenario
mse.rte.list <- list()            # list to store mse.rte.scen.mat for each scenario


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
  
  # After iterating through num.per.scen files, store results in lists and change scenario indicator
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
#   "equal-samp-alpha0p5"
#   "equal-samp-Q5"
#   "equal-samp-Q12"
##################################################################################################


### Combine results from all three simulation studies with different values of Q
### NOTE: must run version "equal-samp-alpha0p5" first
if(sim.version == "equal-samp-alpha0p5"){
  
  rr.list.new <- rr.list
  mse.list.new <- mse.rte.list
  for(i in 1:(S+1)){
    
    # Rejection rates
    rr.mat <- rbind( rr.list[[i]][1,2:(S+2)],
                     rr.list[[i]][3,2:(S+2)],
                     rep(0, S + 1),
                     rr.list[[i]][2,2:(S+2)],
                     rep(0, S + 1) )
    rownames(rr.mat) <- c("CPHM", "BMA-S", "BMA-JM-Q5", "BMA-JM-Q8", "BMA-JM-Q12")
    rr.list.new[[i]] <- rr.mat
    
    # MSE
    mse.mat <- rbind( mse.rte.list[[i]][1,2:(S+1)],
                      mse.rte.list[[i]][3,2:(S+1)],
                      rep(0, S),
                      mse.rte.list[[i]][2,2:(S+1)],
                      rep(0, S) )
    rownames(mse.mat) <- c("CPHM", "BMA-S", "BMA-JM-Q5", "BMA-JM-Q8", "BMA-JM-Q12")
    mse.list.new[[i]] <- mse.mat
    
  }
  
} else if(sim.version == "equal-samp-Q5"){
  
  for(i in 1:(S+1)){
    
    # Rejection rates
    rr.list.new[[i]][3,] <- rr.list[[i]][2,2:(S+2)]
    
    # MSE
    mse.list.new[[i]][3,] <- mse.rte.list[[i]][2,2:(S+1)]
    
  }
  
} else if(sim.version == "equal-samp-Q12"){
  
  for(i in 1:(S+1)){
    
    # Rejection rates
    rr.list.new[[i]][5,] <- rr.list[[i]][2,2:(S+2)]
    
    # MSE
    mse.list.new[[i]][5,] <- mse.rte.list[[i]][2,2:(S+1)]
    
  }
  
}



### Combine rejection rates according to FPR and TPR for each scenario

# Global rejection rates
grr.CPHM <- c( rr.list.new[[1]][1,1], rr.list.new[[2]][1,1], rr.list.new[[3]][1,1],
               rr.list.new[[4]][1,1], rr.list.new[[5]][1,1] )
grr.BMA.S <- c( rr.list.new[[1]][2,1], rr.list.new[[2]][2,1], rr.list.new[[3]][2,1],
                rr.list.new[[4]][2,1], rr.list.new[[5]][2,1] )
grr.BMA.JM.Q5 <- c( rr.list.new[[1]][3,1], rr.list.new[[2]][3,1], rr.list.new[[3]][3,1],
                    rr.list.new[[4]][3,1], rr.list.new[[5]][3,1] )
grr.BMA.JM.Q8 <- c( rr.list.new[[1]][4,1], rr.list.new[[2]][4,1], rr.list.new[[3]][4,1],
                    rr.list.new[[4]][4,1], rr.list.new[[5]][4,1] )
grr.BMA.JM.Q12 <- c( rr.list.new[[1]][5,1], rr.list.new[[2]][5,1], rr.list.new[[3]][5,1],
                     rr.list.new[[4]][5,1], rr.list.new[[5]][5,1] )

# True positive rates
tpr.CPHM <- c( rowMeans(rr.list.new[[1]][1,2:5]), rowMeans(rr.list.new[[2]][1,3:5]),
               rowMeans(rr.list.new[[3]][1,4:5]), rr.list.new[[4]][1,5], NULL )
tpr.BMA.S <- c( rowMeans(rr.list.new[[1]][2,2:5]), rowMeans(rr.list.new[[2]][2,3:5]),
                 rowMeans(rr.list.new[[3]][2,4:5]), rr.list.new[[4]][2,5], NULL )
tpr.BMA.JM.Q5 <- c( rowMeans(rr.list.new[[1]][3,2:5]), rowMeans(rr.list.new[[2]][3,3:5]),
                    rowMeans(rr.list.new[[3]][3,4:5]), rr.list.new[[4]][3,5], NULL )
tpr.BMA.JM.Q8 <- c( rowMeans(rr.list.new[[1]][4,2:5]), rowMeans(rr.list.new[[2]][4,3:5]),
                    rowMeans(rr.list.new[[3]][4,4:5]), rr.list.new[[4]][4,5], NULL )
tpr.BMA.JM.Q12 <- c( rowMeans(rr.list.new[[1]][5,2:5]), rowMeans(rr.list.new[[2]][5,3:5]),
                     rowMeans(rr.list.new[[3]][5,4:5]), rr.list.new[[4]][5,5], NULL )

# False positive rates
fpr.CPHM <- c( NULL, rr.list.new[[2]][1,2], rowMeans(rr.list.new[[3]][1,2:3]),
               rowMeans(rr.list.new[[4]][1,2:4]), rowMeans(rr.list.new[[5]][1,2:5]) )
fpr.BMA.S <- c( NULL, rr.list.new[[2]][2,2], rowMeans(rr.list.new[[3]][2,2:3]),
                 rowMeans(rr.list.new[[4]][2,2:4]), rowMeans(rr.list.new[[5]][2,2:5]) )
fpr.BMA.JM.Q5 <- c( NULL, rr.list.new[[2]][3,2], rowMeans(rr.list.new[[3]][3,2:3]),
                    rowMeans(rr.list.new[[4]][3,2:4]), rowMeans(rr.list.new[[5]][3,2:5]) )
fpr.BMA.JM.Q8 <- c( NULL, rr.list.new[[2]][4,2], rowMeans(rr.list.new[[3]][4,2:3]),
                    rowMeans(rr.list.new[[4]][4,2:4]), rowMeans(rr.list.new[[5]][4,2:5]) )
fpr.BMA.JM.Q12 <- c( NULL, rr.list.new[[2]][5,2], rowMeans(rr.list.new[[3]][5,2:3]),
                     rowMeans(rr.list.new[[4]][5,2:4]), rowMeans(rr.list.new[[5]][5,2:5]) )



### Combine regional treatment effect MSE according to null and alternative regions for each scenario

## Relative MSE for regional treatment effects with CPHM as reference

rel.MSE.mat.alt <- matrix(0, nrow = S, ncol = 5)
rel.MSE.mat.null <- matrix(0, nrow = S, ncol = 5)

# Scenario 1: 0 nulls
# alternative regions
rel.MSE.mat.alt[1,] <- rowMeans(mse.list.new[[1]]) / rowMeans(mse.list.new[[1]][1,])

# Scenario 2: 1 null
# alternative regions
rel.MSE.mat.alt[2,] <- rowMeans(mse.list.new[[2]][,2:4]) / rowMeans(mse.list.new[[2]][1,2:4])
# null regions
rel.MSE.mat.null[1,] <- mse.list.new[[2]][,1] / mse.list.new[[2]][1,1]

# Scenario 3: 2 nulls
# alternative regions
rel.MSE.mat.alt[3,] <- rowMeans(mse.list.new[[3]][,3:4]) / rowMeans(mse.list.new[[3]][1,3:4])
# null regions
rel.MSE.mat.null[2,] <- rowMeans(mse.list.new[[3]][,1:2]) / rowMeans(mse.list.new[[3]][1,1:2])

# Scenario 4: 3 nulls
# alternative regions
rel.MSE.mat.alt[4,] <- mse.list.new[[4]][,4] / mse.list.new[[4]][1,4]
# null regions
rel.MSE.mat.null[3,] <- rowMeans(mse.list.new[[4]][,1:3]) / rowMeans(mse.list.new[[4]][1,1:3])

# Scenario 5: 4 nulls
# null regions
rel.MSE.mat.null[4,] <- rowMeans(mse.list.new[[5]]) / rowMeans(mse.list.new[[5]][1,])



### Plot rejection rates and relative MSE with scenario (number null regions) as x-axis

# Names for legend
names.SA <- c( "CPHM", "BMA-S", "BMA-JM: Q = 5", "BMA-JM: Q = 8", "BMA-JM: Q = 12" )

# Colors of bars
colors.SA <- c("red3", "mediumorchid3", "cornflowerblue", "blue", "darkblue")

# Create PDF file
five.plot.file <- paste( "/Users/nathanbean/Downloads/SA-Q.pdf", sep = "" )
pdf(five.plot.file, height = 12, width = 10)

# Adjust margins and layout
par( mar = c(4.1, 4.1, 1, 1) )
par( mfrow = c(3,2) )
layout( mat=matrix(c(1,1,2,3,4,5), nrow = 3, byrow = TRUE) )

# Global rejection rate plot
grr.mat <- matrix( rbind(grr.CPHM, grr.BMA.S, grr.BMA.JM.Q5, grr.BMA.JM.Q8, grr.BMA.JM.Q12),
                   byrow = FALSE, nrow = 1 )
barplot( c(grr.mat), xlab = "", ylab = "", ylim = c(0,1),
         yaxt = "n", space = c( rep(c(.5,0,0,0,0), 5) ), col = rep(colors.SA, 6) )
axis( 1, at = c(3, 8.5, 14, 19.5, 25), labels = 0:4, tick = FALSE, line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, 1, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "Global Rejection Rate", side = 2, line = 2.8, cex = 1.1 )
abline( h = .025, lty = 3, lwd = 2 )
text(.05, .95, "A", cex = 2)
legend( "topright", bty = "n", legend = names.SA, fill = colors.SA, cex = 1.1 )
box()

# Relative MSE plot for alternative regions
MSE.alt.lim <- 1.2
barplot( c(t(rel.MSE.mat.alt[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim),
         yaxt = "n", space = rep(c(.5,0,0,0), 4), col = rep(colors.SA[-1], 4) )
axis( 1, at = c(2.5, 7, 11.5, 16), labels = 0:3, tick = FALSE, line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.alt.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE (alt. regions)", side = 2, line = 2.8, cex = 1.1 )
text(.5, .95*MSE.alt.lim, "B", cex = 2)
legend( "topright", bty = "n", legend = names.SA[-1], fill = colors.SA[-1],
        ncol = 2, cex = 1.1 )
box()

# Relative MSE plot for null regions
MSE.null.lim <- 1.2
barplot( c(t(rel.MSE.mat.null[,-1])), xlab = "", ylab = "", ylim = c(0,MSE.null.lim),
         yaxt = "n", space = rep(c(.5,0,0,0), 4), col = rep(colors.SA[-1], 4) )
axis( 1, at = c(2.5, 7, 11.5, 16), labels = 1:4, tick = FALSE, line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.null.lim, by = 0.2), las = 1, cex.axis = 1.2 )
mtext( "Relative MSE (null regions)", side = 2, line = 2.8, cex = 1.1 )
abline( h = 1, lty = 3, lwd = 2 )
text(.5, .95*MSE.null.lim, "C", cex = 2)
legend( "topright", bty = "n", legend = names.SA[-1], fill = colors.SA[-1],
        ncol = 2, cex = 1.1 )
box()

# True positive rate plot
TPR.lim <- .7
tpr.mat <- matrix( rbind(tpr.CPHM, tpr.BMA.S, tpr.BMA.JM.Q5, tpr.BMA.JM.Q8, tpr.BMA.JM.Q12),
                   byrow = FALSE, nrow = 1 )
barplot( c(tpr.mat), xlab = "", ylab = "", ylim = c(0,TPR.lim), yaxt = "n",
         space = c( rep(c(.5,0,0,0,0), 4) ), col = rep(colors.SA, 4) )
axis( 1, at = c(3, 8.5, 14, 19.5), labels = 0:3, tick = FALSE, line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, TPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
mtext( "True Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(.5, .95*TPR.lim, "D", cex = 2)
legend( "topright", bty = "n", legend = names.SA, fill = colors.SA, cex = 1.1 )
box()

# False positive rate plot
FPR.lim <- 0.3
fpr.mat <- matrix( rbind(fpr.CPHM, fpr.BMA.S, fpr.BMA.JM.Q5, fpr.BMA.JM.Q8, fpr.BMA.JM.Q12),
                   byrow = FALSE, nrow = 1 )
barplot( c(fpr.mat), xlab = "", ylab = "", ylim = c(0,FPR.lim), yaxt = "n",
         space = c( rep(c(.5,0,0,0,0), 4) ), col = rep(colors.SA, 4) )
axis( 1, at = c(3, 8.5, 14, 19.5), labels = 1:4, tick = FALSE, line = -.3, cex.axis = 1.2 )
axis( 2, at = seq(0, FPR.lim, by = 0.1), las = 1, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.8, cex = 1.1 )
mtext( "False Positive Rate", side = 2, line = 2.8, cex = 1.1 )
text(.5, .95*FPR.lim, "E", cex = 2)
abline( h = .025, lty = 3, lwd = 2 )
legend( "topright", bty = "n", legend = names.SA, fill = colors.SA, cex = 1.1 )
box()

layout( mat=matrix(1, nrow = 1, ncol = 1) )

# Save PDF
dev.off()

