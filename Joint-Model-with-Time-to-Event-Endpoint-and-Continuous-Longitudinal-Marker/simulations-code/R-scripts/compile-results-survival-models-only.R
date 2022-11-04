##################################################################################################
# COMPILE RESULTS (E.G., REJECTION RATES, MSE) FROM SIMULATION SCENARIO
# COMPARING SURVIVAL MODELS ONLY
##################################################################################################



### Set working directory
setwd("/Users/nathanbean/Documents/Dissertation/Project 3/Project3_Simulations/cluster-results")


### Simulation details

# Simulation version
sim.version <- "survival-models-only"
S <- 4      # number of regions

# Total number of result files from simulation
tot.num.results <- 20
num.scenarios <- 5                                    # number of unique scenarios for most versions
num.per.scen <- tot.num.results / num.scenarios       # number of result files per scenario



### Combine results within each scenario for rejection rate, bias, and MSE

## Rejection rate results
rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", 1, ".csv", sep = "")
rr.scen.mat <- read.csv(rr.file, header = TRUE)
rr.scen.mat[1:4, 2:(S+2)] <- 0           # matrix to store rejection rates for scenario

## Bias (regional treatment effects) results
bias.rte.file <- paste("./", sim.version, "/bias/rte_bias_", sim.version, "-", 1, ".csv", sep = "")
bias.rte.scen.mat <- read.csv(bias.rte.file, header = TRUE)
bias.rte.scen.mat[1:4, 2:(S+1)] <- 0     # matrix to store rte biases for scenario

## MSE (regional treatment effects) results
mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", 1, ".csv", sep = "")
mse.rte.scen.mat <- read.csv(mse.rte.file, header = TRUE)
mse.rte.scen.mat[1:4, 2:(S+1)] <- 0     # matrix to store MSE biases for scenario

## MSE (regional treatment effects) results
var.rte.file <- paste("./", sim.version, "/sd/sd_", sim.version, "-", 1, ".csv", sep = "")
var.rte.scen.mat <- read.csv(var.rte.file, header = TRUE)
var.rte.scen.mat[1:4, 1:S] <- 0      # matrix to store rte SDs for scenario
row.names(var.rte.scen.mat) <- c("CPHM-surv", "CPHM-joint", "BMA-S-surv", "BMA-S-joint")


rr.list <- list()                 # list to store rr.scen.mat for each scenario
bias.rte.list <- list()           # list to store bias.rte.scen.mat for each scenario
mse.rte.list <- list()           # list to store mse.rte.scen.mat for each scenario
var.rte.list <- list()            # list to store var.rte.scen.mat for each scenario


which.scen <- 1                   # scenario indicator
for(j in 1:tot.num.results){
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen)
  rr.file <- paste("./", sim.version, "/rejection_rates/rr_", sim.version, "-", j, ".csv", sep = "")
  rr.part.mat <- read.csv(rr.file, header = TRUE)
  rr.scen.mat[1:4, 2:(S+2)] <- rr.scen.mat[1:4, 2:(S+2)] + rr.part.mat[1:4, 2:(S+2)]
  
  # Matrix of bias (rte) with part of simulation results (1 out of num.per.scen)
  bias.rte.file <- paste("./", sim.version, "/bias/rte_bias_", sim.version, "-", j, ".csv", sep = "")
  bias.rte.part.mat <- read.csv(bias.rte.file, header = TRUE)
  bias.rte.scen.mat[1:4, 2:(S+1)] <- bias.rte.scen.mat[1:4, 2:(S+1)] + bias.rte.part.mat[1:4, 2:(S+1)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen)
  mse.rte.file <- paste("./", sim.version, "/mse/rte_mse_", sim.version, "-", j, ".csv", sep = "")
  mse.rte.part.mat <- read.csv(mse.rte.file, header = TRUE)
  mse.rte.scen.mat[1:4, 2:(S+1)] <- mse.rte.scen.mat[1:4, 2:(S+1)] + mse.rte.part.mat[1:4, 2:(S+1)]
  
  # Matrix of variances (rte) with part of simulation results (1 out of num.per.scen)
  var.rte.file <- paste("./", sim.version, "/sd/sd_", sim.version, "-", j, ".csv", sep = "")
  var.rte.part.mat <- read.csv(var.rte.file, header = TRUE)^2
  var.rte.scen.mat[1:4, 1:S] <- var.rte.scen.mat[1:4, 1:S] + var.rte.part.mat[1:4, 1:S]
  
  # After iterating through num.per.scen files, store results in lists and change scenario indicator
  if(j %% num.per.scen == 0){
    
    # Average results
    rr.scen.mat[1:4, 2:(S+2)] <- rr.scen.mat[1:4, 2:(S+2)] / num.per.scen
    bias.rte.scen.mat[1:4, 2:(S+1)] <- bias.rte.scen.mat[1:4, 2:(S+1)] / num.per.scen
    var.rte.scen.mat[1:4, 1:S] <- var.rte.scen.mat[1:4, 1:S] / num.per.scen
    
    # Store results in list
    scen.lab <- paste("scenario_", which.scen, sep = "")
    rr.list[[scen.lab]] <- rr.scen.mat                 # save scenario-specific matrix of rr's
    bias.rte.list[[scen.lab]] <- bias.rte.scen.mat     # save scenario-specific matrix of bias.rte's
    mse.rte.list[[scen.lab]] <- mse.rte.scen.mat       # save scenario-specific matrix of mse.rte's
    var.rte.list[[scen.lab]] <- var.rte.scen.mat       # save scenario-specific matrix of var.rte's
    
    # Reset matrices and prepare for next scenario
    which.scen <- which.scen + 1                       # increase scenario indicator
    rr.scen.mat[1:4, 2:(S+2)] <- 0
    bias.rte.scen.mat[1:4, 2:(S+1)] <- 0
    mse.rte.scen.mat[1:4, 2:(S+1)] <- 0
    var.rte.scen.mat[1:4, 1:S] <- 0
    
  }
  
}


rr.list
bias.rte.list
mse.rte.list
var.rte.list 



### Combine rejection rates according to FPR and TPR for each scenario

# True positive rates
tpr.CPHM.surv <- c( rowMeans(rr.list[[1]][1,3:6]), rowMeans(rr.list[[2]][1,4:6]),
                    rowMeans(rr.list[[3]][1,5:6]), rr.list[[4]][1,6], NULL )
tpr.CPHM.joint <- c( rowMeans(rr.list[[1]][2,3:6]), rowMeans(rr.list[[2]][2,4:6]),
                     rowMeans(rr.list[[3]][2,5:6]), rr.list[[4]][2,6], NULL )
tpr.BMA.S.surv <- c( rowMeans(rr.list[[1]][3,3:6]), rowMeans(rr.list[[2]][3,4:6]),
                     rowMeans(rr.list[[3]][3,5:6]), rr.list[[4]][3,6], NULL )
tpr.BMA.S.joint <- c( rowMeans(rr.list[[1]][4,3:6]), rowMeans(rr.list[[2]][4,4:6]),
                      rowMeans(rr.list[[3]][4,5:6]), rr.list[[4]][4,6], NULL )

# False positive rates
fpr.CPHM.surv <- c( NULL, rr.list[[2]][1,3], rowMeans(rr.list[[3]][1,3:4]),
                    rowMeans(rr.list[[4]][1,3:5]), rowMeans(rr.list[[5]][1,3:6]) )
fpr.CPHM.joint <- c( NULL, rr.list[[2]][2,3], rowMeans(rr.list[[3]][2,3:4]),
                     rowMeans(rr.list[[4]][2,3:5]), rowMeans(rr.list[[5]][2,3:6]) )
fpr.BMA.S.surv <- c( NULL, rr.list[[2]][3,3], rowMeans(rr.list[[3]][3,3:4]),
                     rowMeans(rr.list[[4]][3,3:5]), rowMeans(rr.list[[5]][3,3:6]) )
fpr.BMA.S.joint <- c( NULL, rr.list[[2]][4,3], rowMeans(rr.list[[3]][4,3:4]),
                      rowMeans(rr.list[[4]][4,3:5]), rowMeans(rr.list[[5]][4,3:6]) )



## MSE for regional treatment effects

MSE.mat.alt <- matrix(0, nrow = S, ncol = 4)
MSE.mat.null <- matrix(0, nrow = S, ncol = 4)

# Scenario 1: 0 nulls
# alternative regions
MSE.mat.alt[1,] <- rowMeans(mse.rte.list[[1]][1:4,2:5])

# Scenario 2: 1 null
# alternative regions
MSE.mat.alt[2,] <- rowMeans(mse.rte.list[[2]][1:4,3:5])
# null regions
MSE.mat.null[1,] <- mse.rte.list[[2]][1:4,2]

# Scenario 3: 2 nulls
# alternative regions
MSE.mat.alt[3,] <- rowMeans(mse.rte.list[[3]][1:4,4:5])
# null regions
MSE.mat.null[2,] <- rowMeans(mse.rte.list[[3]][1:4,2:3])

# Scenario 4: 3 nulls
# alternative regions
MSE.mat.alt[4,] <- mse.rte.list[[4]][1:4,5]
# null regions
MSE.mat.null[3,] <- rowMeans(mse.rte.list[[4]][1:4,2:4])

# Scenario 5: 4 nulls
# null regions
MSE.mat.null[4,] <- rowMeans(mse.rte.list[[5]][1:4,2:5])



## Variance for regional treatment effects

var.mat.alt <- matrix(0, nrow = S, ncol = 4)
var.mat.null <- matrix(0, nrow = S, ncol = 4)

# Scenario 1: 0 nulls
# alternative regions
var.mat.alt[1,] <- rowMeans(var.rte.list[[1]][1:4,1:4])

# Scenario 2: 1 null
# alternative regions
var.mat.alt[2,] <- rowMeans(var.rte.list[[2]][1:4,2:4])
# null regions
var.mat.null[1,] <- var.rte.list[[2]][1:4,1]

# Scenario 3: 2 nulls
# alternative regions
var.mat.alt[3,] <- rowMeans(var.rte.list[[3]][1:4,3:4])
# null regions
var.mat.null[2,] <- rowMeans(var.rte.list[[3]][1:4,1:2])

# Scenario 4: 3 nulls
# alternative regions
var.mat.alt[4,] <- var.rte.list[[4]][1:4,4]
# null regions
var.mat.null[3,] <- rowMeans(var.rte.list[[4]][1:4,1:3])

# Scenario 5: 4 nulls
# null regions
var.mat.null[4,] <- rowMeans(var.rte.list[[5]][1:4,1:4])



### Plot rejection rates and relative MSE with scenario (number null regions) as x-axis

# Create PDF file
surv.mods.plot.file <- paste( "/Users/nathanbean/Downloads/", sim.version, ".pdf", sep = "" )
#pdf(surv.mods.plot.file, height = 6.93, width = 10.73)
pdf(surv.mods.plot.file, height = 7.62, width = 11.80)

# Adjust margins and layout
par( mar = c(4.1, 4.8, 1.2, 1.2) )
par( mfrow = c(2,2) )

# Names for legend
legend.names <- c(expression(paste("CPHM (", alpha, " = 0)", sep = "")),
                  expression(paste("CPHM (", alpha, " = 1)", sep = "")),
                  expression(paste("BMA-S (", alpha, " = 0)", sep = "")),
                  expression(paste("BMA-S (", alpha, " = 1)", sep = "")))

# MSE plot for alternative regions
MSE.alt.lim <- 0.06
barplot( c(t(MSE.mat.alt)), xlab = "", ylab = "", ylim = c(0,MSE.alt.lim),
         yaxt = "n", space = rep(c(.5,0,0,0), 4), density = c(-1,30,-1,30), angle = 45,
         col = rep(c("gray30", "gray30", "gray80", "gray80"), 4) )
axis( 1, at = c(2.5, 7, 11.5, 16), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.alt.lim, by = 0.01), las = 1, cex.axis = 1.2 )
mtext( "MSE (alternative regions)", side = 2, line = 3.7, cex = 1.1 )
text(.5, .95*MSE.alt.lim, "A", cex = 2)
legend( "topright", bty = "n", fill = c("gray30", "gray30", "gray80", "gray80"),
        density = c(-1,30,-1,30), angle = 45, cex = 1.1, ncol = 2, legend = legend.names )
box()

# MSE plot for null regions
MSE.null.lim <- 0.06
barplot( c(t(MSE.mat.null)), xlab = "", ylab = "", ylim = c(0,MSE.null.lim),
         yaxt = "n", space = rep(c(.5,0,0,0), 4), density = c(-1,30,-1,30), angle = 45,
         col = rep(c("gray30", "gray30", "gray80", "gray80"), 4) )
axis( 1, at = c(2.5, 7, 11.5, 16), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, MSE.null.lim, by = 0.01), las = 1, cex.axis = 1.2 )
mtext( "MSE (null regions)", side = 2, line = 3.7, cex = 1.1 )
abline( h = 1, lty = 3, lwd = 2 )
text(.5, .95*MSE.null.lim, "B", cex = 2)
legend( "topright", bty = "n", fill = c("gray30", "gray30", "gray80", "gray80"),
        density = c(-1,30,-1,30), angle = 45, cex = 1.1, ncol = 2, legend = legend.names )
box()

# Variance plot for alternative regions
var.alt.lim <- 0.015
barplot( c(t(var.mat.alt)), xlab = "", ylab = "", ylim = c(0,var.alt.lim), yaxt = "n",
         space = rep(c(.5,0,0,0), 4), density = c(-1,30,-1,30), angle = 45,
         col = rep(c("gray30", "gray30", "gray80", "gray80"), 4) )
axis( 1, at = c(2.5, 7, 11.5, 16), labels = 0:3, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
abline( h = 1, lty = 3, lwd = 2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, var.alt.lim, by = 0.003), las = 1, cex.axis = 1.2 )
mtext( "Variance (alternative regions)", side = 2, line = 3.7, cex = 1.1 )
text(.5, .95*var.alt.lim, "C", cex = 2)
legend( "topright", bty = "n", fill = c("gray30", "gray30", "gray80", "gray80"),
        density = c(-1,30,-1,30), angle = 45, cex = 1.1, ncol = 2, legend = legend.names )
box()

# Variance plot for null regions
var.null.lim <- 0.015
barplot( c(t(var.mat.null)), xlab = "", ylab = "", ylim = c(0,var.null.lim), yaxt = "n",
         space = rep(c(.5,0,0,0), 4), density = c(-1,30,-1,30), angle = 45,
         col = rep(c("gray30", "gray30", "gray80", "gray80"), 4)  )
axis( 1, at = c(2.5, 7, 11.5, 16), labels = 1:4, tick = FALSE,
      line = -.3, cex.axis = 1.2 )
mtext( "Number of Null Regions", side = 1, line = 2.3, cex = 1.1 )
axis( 2, at = seq(0, var.null.lim, by = 0.003), las = 1, cex.axis = 1.2 )
mtext( "Variance (null regions)", side = 2, line = 3.7, cex = 1.1 )
abline( h = 1, lty = 3, lwd = 2 )
text(.5, .95*var.null.lim, "D", cex = 2)
legend( "topright", bty = "n", fill = c("gray30", "gray30", "gray80", "gray80"),
        density = c(-1,30,-1,30), angle = 45, cex = 1.1, ncol = 2, legend = legend.names )
box()

# Save PDF
dev.off()

