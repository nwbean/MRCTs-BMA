##################################################################################################
# COMPILE GLOBAL CONSISTENCY RESULTS FROM SIMULATIONS
#
# Execute code from either local computer or from cluster
##################################################################################################



### Set working directory and source Rcpp files
setwd("/Users/nathanbean/Documents/Dissertation/Project 2/Project2_Simulations/cluster-results")
library(Rcpp)
sourceCpp("../R-source/bma_functions_tte_endpoint.cpp")


### Simulation details

# Simulation version
sim.version <- "n9340-equal-samp"
#sim.version <- "n9340-equal-samp-beta0_8"
#sim.version <- "n9340-equal-samp-beta0_2"
#sim.version <- "n9340-original"
S <- 4      # number of regions

# Total number of result files from simulation
tot.num.results <- 50                                # use for version "n9340-equal-samp"
#tot.num.results <- 40                                # use for version "n9340-original"
num.scenarios <- 5                                   # number of unique scenarios
#num.scenarios <- 2                                   # use for version "n9340-original"
num.per.scen <- tot.num.results / num.scenarios       # number of result files per scenario



### Concatenate PMPs within each scenario for global consistency ratios

gc.scen.mat <- NULL
gc.list <- list()          # list to store gc.scen.mat for each scenario

which.scen <- 1             # scenario indicator
for(j in 1:tot.num.results){
  
  # Matrix of global consistency probabilities with fraction of simulation results (1 out of num.per.scen)
  gc.file <- paste("./", sim.version, "/glob_consis/gc_", sim.version, "-", j, ".csv", sep = "")
  gc.part.mat <- read.csv(gc.file, header = TRUE)
  gc.scen.mat <- cbind(gc.scen.mat, gc.part.mat[,2])
  
  # After iterating through num.per.scen files, store results in lists and change scenario indicator
  if(j %% num.per.scen == 0){
    
    # Store results for global consistency in list
    rownames(gc.scen.mat) <- gc.part.mat[,1]                # list epsilon.star values as row names
    scen.lab <- paste("scenario_", which.scen, sep = "")
    gc.list[[scen.lab]] <- rowMeans(gc.scen.mat)            # save scenario-specific matrix of results
    
    # Reset matrix and prepare for next scenario
    which.scen <- which.scen + 1                      # increase scenario indicator
    gc.scen.mat <- NULL

  }
    
}


### Compile results for global consistency probabilities

## Save results for different scenarios
if(sim.version == "n9340-equal-samp"){
  gc.n9340.1 <- gc.list[[1]]
  gc.n9340.2 <- gc.list[[2]]
  gc.n9340.3 <- gc.list[[3]]
  gc.n9340.4 <- gc.list[[4]]
  gc.n9340.5 <- gc.list[[5]]
}
if(sim.version == "n9340-equal-samp-beta0_8"){
  gc.n9340.beta0_8.1 <- gc.list[[1]]
  gc.n9340.beta0_8.2 <- gc.list[[2]]
  gc.n9340.beta0_8.3 <- gc.list[[3]]
  gc.n9340.beta0_8.4 <- gc.list[[4]]
  gc.n9340.beta0_8.5 <- gc.list[[5]]
}
if(sim.version == "n9340-equal-samp-beta0_2"){
  gc.n9340.beta0_2.1 <- gc.list[[1]]
  gc.n9340.beta0_2.2 <- gc.list[[2]]
  gc.n9340.beta0_2.3 <- gc.list[[3]]
  gc.n9340.beta0_2.4 <- gc.list[[4]]
  gc.n9340.beta0_2.5 <- gc.list[[5]]
}
if(sim.version == "n9340-original"){
  gc.n9340.original.1 <- gc.list[[1]]
  gc.n9340.original.2 <- gc.list[[2]]
}



###################################################################################
# Version: n9340-equal-samp (and versions with different beta* values)
###################################################################################

### Plot global consistency probabilities for varying values of epsilon.star

## COLOR SCALE

# Create PDF file
equal.samp.plot.file <- paste( "/Users/nathanbean/Downloads/glob-consis-equal-samp.pdf", sep = "" )
#pdf(equal.samp.plot.file, height = 4, width = 13)
pdf(equal.samp.plot.file, height = 10, width = 6.5)

epsilon.star <- as.numeric( names(gc.list[[1]]) )
par( mar = c(3.9, 4.1, 2.0, 1.2) )
#par( mfrow = c(1,3) )
par( mfrow = c(3,1) )

# Global consistency plot for beta = 0.2
plot( epsilon.star, gc.n9340.beta0_2.1, xlab = "", ylab = "", xlim = c(min(epsilon.star),1),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "midnightblue", cex.axis = 1.2 )
axis( 2, at = c(0, .2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(beta*"*", " = 0.2") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n9340.beta0_2.2, lwd = 3, col = "lightblue" )
lines( epsilon.star, gc.n9340.beta0_2.3, lwd = 3, col = "cornflowerblue" )
lines( epsilon.star, gc.n9340.beta0_2.4, lwd = 3, col = "midnightblue", lty = 2 )
lines( epsilon.star, gc.n9340.beta0_2.5, lwd = 3, col = "lightblue", lty = 2 )
legend( "bottomleft", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                        "3 Null Regions", "4 Null Regions"), bty = "n", lwd = 2,
        col = c("midnightblue", "lightblue", "cornflowerblue", "midnightblue",
                "lightblue"), lty = c(1,1,1,2,2), cex = 1.1 )

# Global consistency plot for beta = 0.5
plot( epsilon.star, gc.n9340.1, xlab = "", ylab = "", xlim = c(min(epsilon.star),1),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "midnightblue", cex.axis = 1.2 )
axis( 2, at = c(0, .2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(beta*"*", " = 0.5") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n9340.2, lwd = 3, col = "lightblue" )
lines( epsilon.star, gc.n9340.3, lwd = 3, col = "cornflowerblue" )
lines( epsilon.star, gc.n9340.4, lwd = 3, col = "midnightblue", lty = 2 )
lines( epsilon.star, gc.n9340.5, lwd = 3, col = "lightblue", lty = 2 )
legend( "bottomleft", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                        "3 Null Regions", "4 Null Regions"), bty = "n", lwd = 2,
        col = c("midnightblue", "lightblue", "cornflowerblue", "midnightblue",
                "lightblue"), lty = c(1,1,1,2,2), cex = 1.1 )

# Global consistency plot for beta = 0.8
plot( epsilon.star, gc.n9340.beta0_8.1, xlab = "", ylab = "", xlim = c(min(epsilon.star),1),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "midnightblue", cex.axis = 1.2 )
axis( 2, at = c(0, .2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(beta*"*", " = 0.8") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n9340.beta0_8.2, lwd = 3, col = "lightblue" )
lines( epsilon.star, gc.n9340.beta0_8.3, lwd = 3, col = "cornflowerblue" )
lines( epsilon.star, gc.n9340.beta0_8.4, lwd = 3, col = "midnightblue", lty = 2 )
lines( epsilon.star, gc.n9340.beta0_8.5, lwd = 3, col = "lightblue", lty = 2 )
legend( "bottomleft", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                        "3 Null Regions", "4 Null Regions"), bty = "n", lwd = 2,
        col = c("midnightblue", "lightblue", "cornflowerblue", "midnightblue",
                "lightblue"), lty = c(1,1,1,2,2), cex = 1.1 )

# Save PDF
dev.off()
par( mfrow = c(1,1) )



###################################################################################
# Version: n9340-original
###################################################################################

### Plot global consistency probabilities for varying values of epsilon.star

## COLOR SCALE

# Create PDF file
original.plot.file <- paste( "/Users/nathanbean/Downloads/glob-consis-original.pdf", sep = "" )
pdf(original.plot.file, height = 5, width = 6.5)

epsilon.star <- as.numeric( names(gc.list[[1]]) )
par( mar = c(3.9, 4.1, 2.0, 1.2) )

# Global consistency plot for beta = 0.5
plot( epsilon.star, gc.n9340.original.1, xlab = "", ylab = "", xlim = c(min(epsilon.star),1),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "midnightblue", cex.axis = 1.2 )
axis( 2, at = c(0, .2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(beta*"*", " = 0.5") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n9340.original.2, lwd = 3, col = "cornflowerblue" )
legend( "bottomleft", c("Same region-specific treatment hazard ratios",
                        "Different region-specific treatment hazard ratios"), bty = "n",
        lwd = 2, col = c("midnightblue", "cornflowerblue"), lty = c(1,1), cex = 1.1 )

# Save PDF
dev.off()

