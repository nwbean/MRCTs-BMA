##################################################################################################
# COMPILE RESULTS (E.G., REJECTION RATES, MSE) FROM SIMULATION SCENARIOS
# COMPARING VALUES OF T0 IN {2, 3, 4}
##################################################################################################



### Set working directory
setwd("/Users/nathanbean/Documents/Dissertation/Project 3/Project3_Simulations/cluster-results")


### Simulation details

# Simulation version
sim.T0.2 <- "original-alpha0p5"
sim.T0.3 <- "original-T0equal3"
sim.T0.4 <- "original-T0equal4"
S <- 4      # number of regions

# Total number of result files from simulation
tot.num.results.2 <- 1000                            # use for T0 = 2
tot.num.results.34 <- 4000                           # use for T0 = 3, 4
num.scenarios <- 2                                   # number of scenarios "original" versions
num.per.scen.2 <- tot.num.results.2 / num.scenarios  # number of result files per scenario (T0 = 2)
num.per.scen.34 <- tot.num.results.34 /
  num.scenarios             # number of result files per scenario (T0 = 3, 4)



### Combine results within each scenario for rejection rate, bias, and MSE

## Rejection rate results (T0 = 2)
rr.file.2 <- paste("./", sim.T0.2, "/rejection_rates/rr_", sim.T0.2, "-", 1, ".csv", sep = "")
rr.scen.mat.2 <- read.csv(rr.file.2, header = TRUE)
rr.scen.mat.2[1:4, 2:(S+2)] <- 0         # matrix to store rejection rates for scenario

## Rejection rate results (T0 = 3)
rr.file.3 <- paste("./", sim.T0.3, "/rejection_rates/rr_", sim.T0.3, "-", 1, ".csv", sep = "")
rr.scen.mat.3 <- read.csv(rr.file.3, header = TRUE)
rr.scen.mat.3[1:4, 2:(S+2)] <- 0         # matrix to store rejection rates for scenario

## Rejection rate results (T0 = 4)
rr.file.4 <- paste("./", sim.T0.4, "/rejection_rates/rr_", sim.T0.4, "-", 1, ".csv", sep = "")
rr.scen.mat.4 <- read.csv(rr.file.4, header = TRUE)
rr.scen.mat.4[1:4, 2:(S+2)] <- 0         # matrix to store rejection rates for scenario

## MSE (regional treatment effects) results (T0 = 2)
mse.rte.file.2 <- paste("./", sim.T0.2, "/mse/rte_mse_", sim.T0.2, "-", 1, ".csv", sep = "")
mse.rte.scen.mat.2 <- read.csv(mse.rte.file.2, header = TRUE)
mse.rte.scen.mat.2[1:4, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario

## MSE (regional treatment effects) results (T0 = 3)
mse.rte.file.3 <- paste("./", sim.T0.3, "/mse/rte_mse_", sim.T0.3, "-", 1, ".csv", sep = "")
mse.rte.scen.mat.3 <- read.csv(mse.rte.file.3, header = TRUE)
mse.rte.scen.mat.3[1:4, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario

## MSE (regional treatment effects) results (T0 = 4)
mse.rte.file.4 <- paste("./", sim.T0.4, "/mse/rte_mse_", sim.T0.4, "-", 1, ".csv", sep = "")
mse.rte.scen.mat.4 <- read.csv(mse.rte.file.4, header = TRUE)
mse.rte.scen.mat.4[1:4, 2:(S+1)] <- 0      # matrix to store rte MSEs for scenario

rr.list.2 <- list()              # list to store rr.scen.mat for each scenario (T0 = 2)
rr.list.3 <- list()              # list to store rr.scen.mat for each scenario (T0 = 3)
rr.list.4 <- list()              # list to store rr.scen.mat for each scenario (T0 = 4)
mse.rte.list.2 <- list()         # list to store mse.rte.scen.mat for each scenario (T0 = 2)
mse.rte.list.3 <- list()         # list to store mse.rte.scen.mat for each scenario (T0 = 3)
mse.rte.list.4 <- list()         # list to store mse.rte.scen.mat for each scenario (T0 = 4)


which.scen.2 <- 1                # scenario indicator (T0 = 2)
which.scen.34 <- 1               # scenario indicator (T0 = 3, 4)
for(j in 1:tot.num.results.34){
  
  if(j <= tot.num.results.2){
    
    # Matrix of rejection rates with part of simulation results (1 out of num.per.scen) (T0 = 2)
    rr.file.2 <- paste("./", sim.T0.2, "/rejection_rates/rr_", sim.T0.2, "-", j, ".csv", sep = "")
    rr.part.mat.2 <- read.csv(rr.file.2, header = TRUE)
    rr.scen.mat.2[1:4, 2:(S+2)] <- rr.scen.mat.2[1:4, 2:(S+2)] + rr.part.mat.2[1:4, 2:(S+2)]
    
    # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen) (T0 = 2)
    mse.rte.file.2 <- paste("./", sim.T0.2, "/mse/rte_mse_", sim.T0.2, "-", j, ".csv", sep = "")
    mse.rte.part.mat.2 <- read.csv(mse.rte.file.2, header = TRUE)
    mse.rte.scen.mat.2[1:4, 2:(S+1)] <- mse.rte.scen.mat.2[1:4, 2:(S+1)] +
      mse.rte.part.mat.2[1:4, 2:(S+1)]
    
    # After iterating through num.per.scen files, store results in lists and change scenario indicator
    if(j %% num.per.scen.2 == 0){
      
      # Average results
      rr.scen.mat.2[1:4, 2:(S+2)] <- rr.scen.mat.2[1:4, 2:(S+2)] / num.per.scen.2
      mse.rte.scen.mat.2[1:4, 2:(S+1)] <- mse.rte.scen.mat.2[1:4, 2:(S+1)] / num.per.scen.2
      
      # Store results in list
      scen.lab <- paste("scenario_", which.scen.2, sep = "")
      rr.list.2[[scen.lab]] <- rr.scen.mat.2             # save scenario-specific matrix of rr's
      mse.rte.list.2[[scen.lab]] <- mse.rte.scen.mat.2   # save scenario-specific matrix of mse.rte's
      
      # Reset matrices and prepare for next scenario
      which.scen.2 <- which.scen.2 + 1                   # increase scenario indicator
      rr.scen.mat.2[1:4, 2:(S+2)] <- 0
      mse.rte.scen.mat.2[1:4, 2:(S+1)] <- 0
      
    }
    
  }
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen) (T0 = 3)
  rr.file.3 <- paste("./", sim.T0.3, "/rejection_rates/rr_", sim.T0.3, "-", j, ".csv", sep = "")
  rr.part.mat.3 <- read.csv(rr.file.3, header = TRUE)
  rr.scen.mat.3[1:4, 2:(S+2)] <- rr.scen.mat.3[1:4, 2:(S+2)] + rr.part.mat.3[1:4, 2:(S+2)]
  
  # Matrix of rejection rates with part of simulation results (1 out of num.per.scen) (T0 = 4)
  rr.file.4 <- paste("./", sim.T0.4, "/rejection_rates/rr_", sim.T0.4, "-", j, ".csv", sep = "")
  rr.part.mat.4 <- read.csv(rr.file.4, header = TRUE)
  rr.scen.mat.4[1:4, 2:(S+2)] <- rr.scen.mat.4[1:4, 2:(S+2)] + rr.part.mat.4[1:4, 2:(S+2)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen) (T0 = 3)
  mse.rte.file.3 <- paste("./", sim.T0.3, "/mse/rte_mse_", sim.T0.3, "-", j, ".csv", sep = "")
  mse.rte.part.mat.3 <- read.csv(mse.rte.file.3, header = TRUE)
  mse.rte.scen.mat.3[1:4, 2:(S+1)] <- mse.rte.scen.mat.3[1:4, 2:(S+1)] +
    mse.rte.part.mat.3[1:4, 2:(S+1)]
  
  # Matrix of MSE (rte) with part of simulation results (1 out of num.per.scen) (T0 = 4)
  mse.rte.file.4 <- paste("./", sim.T0.4, "/mse/rte_mse_", sim.T0.4, "-", j, ".csv", sep = "")
  mse.rte.part.mat.4 <- read.csv(mse.rte.file.4, header = TRUE)
  mse.rte.scen.mat.4[1:4, 2:(S+1)] <- mse.rte.scen.mat.4[1:4, 2:(S+1)] +
    mse.rte.part.mat.4[1:4, 2:(S+1)]
  
  # After iterating through num.per.scen files, store results in lists and change scenario indicator
  if(j %% num.per.scen.34 == 0){
    
    # Average results
    rr.scen.mat.3[1:4, 2:(S+2)] <- rr.scen.mat.3[1:4, 2:(S+2)] / num.per.scen.34
    rr.scen.mat.4[1:4, 2:(S+2)] <- rr.scen.mat.4[1:4, 2:(S+2)] / num.per.scen.34
    mse.rte.scen.mat.3[1:4, 2:(S+1)] <- mse.rte.scen.mat.3[1:4, 2:(S+1)] / num.per.scen.34
    mse.rte.scen.mat.4[1:4, 2:(S+1)] <- mse.rte.scen.mat.4[1:4, 2:(S+1)] / num.per.scen.34
    
    # Store results in list
    scen.lab <- paste("scenario_", which.scen.34, sep = "")
    rr.list.3[[scen.lab]] <- rr.scen.mat.3             # save scenario-specific matrix of rr's
    rr.list.4[[scen.lab]] <- rr.scen.mat.4             # save scenario-specific matrix of rr's
    mse.rte.list.3[[scen.lab]] <- mse.rte.scen.mat.3   # save scenario-specific matrix of mse.rte's
    mse.rte.list.4[[scen.lab]] <- mse.rte.scen.mat.4   # save scenario-specific matrix of mse.rte's
    
    # Reset matrices and prepare for next scenario
    which.scen.34 <- which.scen.34 + 1                 # increase scenario indicator
    rr.scen.mat.3[1:4, 2:(S+2)] <- 0
    rr.scen.mat.4[1:4, 2:(S+2)] <- 0
    mse.rte.scen.mat.3[1:4, 2:(S+1)] <- 0
    mse.rte.scen.mat.4[1:4, 2:(S+1)] <- 0
    
  }
  
}


### Plot results to compare rejection rates and MSE for T0 = 2, 3, 4

# Save rejection rates and relative MSE for same treatment effect scenario
rr.mat1 <- as.matrix( rbind( rr.list.2[[1]][2,-1], rr.list.3[[1]][2,-1], rr.list.4[[1]][2,-1] ) )
rel.MSE.mat.diff1 <- matrix( 0, nrow = 3, ncol = S )
for(i in 1:S){
  rel.MSE.mat.diff1[,i] <- rbind(mse.rte.list.2[[1]][2,i+1], mse.rte.list.3[[1]][2,i+1],
                                 mse.rte.list.4[[1]][2,i+1]) / mse.rte.list.4[[1]][2,i+1]
}

# Save rejection rates and relative MSE for different treatment effects scenario
rr.mat2 <- as.matrix( rbind( rr.list.2[[2]][2,-1], rr.list.3[[2]][2,-1], rr.list.4[[2]][2,-1] ) )
rel.MSE.mat.diff2 <- matrix( 0, nrow = 3, ncol = S )
for(i in 1:S){
  rel.MSE.mat.diff2[,i] <- rbind(mse.rte.list.2[[2]][2,i+1], mse.rte.list.3[[2]][2,i+1],
                                 mse.rte.list.4[[2]][2,i+1]) / mse.rte.list.4[[2]][2,i+1]
}

# Create PDF file
four.plot.file <- paste( "/Users/nathanbean/Downloads/T0-comparison.pdf", sep = "" )
pdf(four.plot.file, height = 6.93, width = 10.73)

par( mar = c(3.6, 3.6, 1, 1) )
par( mfrow = c(2,2) )

# Rejection rates plot (same treatment effect)
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
legend( "topright", bty = "n", fill = c("gray50", "gray70", "gray90"),
        legend = c(expression(paste(D[0], " = 2")),
                   expression(paste(D[0], " = 3")),
                   expression(paste(D[0], " = 4"))), cex = .9 )
box()

# Relative MSE plot (same treatment effect)
barplot( c(rel.MSE.mat.diff1[-1,]), xlab = "", ylab = "", ylim = c(0,1.2),
         yaxt = "n", space = rep(c(.5,0), 4),
         col = rep(c("gray50", "gray80"), 4) )
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
legend( "topright", bty = "n", fill = c("gray50", "gray80"),
        legend = c(expression(paste(D[0], " = 2")), expression(paste(D[0], " = 3"))), cex = .9 )
box()

# Rejection rates plot (different treatment effect)
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
legend( "topright", bty = "n", fill = c("gray50", "gray70", "gray90"),
        legend = c(expression(paste(D[0], " = 2")),
                   expression(paste(D[0], " = 3")),
                   expression(paste(D[0], " = 4"))), cex = .9 )
box()

# Relative MSE plot (different treatment effect)
MSE.ylim.diff <- 1.2
barplot( c(rel.MSE.mat.diff2[-1,]), xlab = "", ylab = "", ylim = c(0, MSE.ylim.diff),
         yaxt = "n", space = rep(c(.5,0), 4),
         col = rep(c("gray50", "gray80"), 4) )
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
legend( "topright", bty = "n", fill = c("gray50", "gray80"),
        legend = c(expression(paste(D[0], " = 2")), expression(paste(D[0], " = 3"))), cex = .9 )
box()

par( mfrow = c(1,1) )

# Save PDF
dev.off()

