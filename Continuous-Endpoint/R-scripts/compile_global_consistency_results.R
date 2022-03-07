##################################################################################################
# COMPILE GLOBAL CONSISTENCY RESULTS FROM SIMULATIONS
#
# Execute code from either local computer or from cluster
##################################################################################################



### Set working directory and source Rcpp files
setwd("/Users/nathanbean/Documents/Dissertation/Project 1/Project1_Simulations_v2/cluster-results")
library(Rcpp)
#sourceCpp("../R-source/bma_functions_rcpp_consistency.cpp")
sourceCpp("../R-source/bma_functions_rcpp_code.cpp")


### Simulation details

# Simulation version
sim.version <- "n1508"
#sim.version <- "n1508-beta20"
#sim.version <- "n754"
#sim.version <- "n754-beta20"
#sim.version <- "n7650"
#sim.version <- "n7650-beta20"
#sim.version <- "null-half"
#sim.version <- "alt-half"
#sim.version <- "diff-effects-1"
#sim.version <- "diff-effects-2"
#sim.version <- "diff-effects-3"
S <- 5      # number of regions

# Total number of result files from simulation
tot.num.results <- 60                               # use for most versions
#tot.num.results <- 1200                            # use for version "n7650" and "n7650-beta20"
#tot.num.results <- 20                              # use for version "diff-effects-1-3"
num.scenarios <- 6                                  # number of unique scenarios
#num.scenarios <- 1                                 # use for version "diff-effects-1-3"
num.per.scen <- tot.num.results / num.scenarios     # number of result files per scenario



### Concatenate PMPs within each scenario for global consistency ratios

gc.scen.mat <- NULL
gic.scen.mat <- NULL
gc.list <- list()          # list to store gc.scen.mat for each scenario
gic.list <- list()         # list to store gic.scen.mat for each scenario

which.scen <- 1             # scenario indicator
for(j in 1:tot.num.results){
  
  # Matrix of global consistency probabilites with fraction of simulation results (1 out of num.per.scen)
  gc.file <- paste("./", sim.version, "/glob_consis/gc_", sim.version, "-", j, ".csv", sep = "")
  gc.part.mat <- read.csv(gc.file, header = TRUE)
  gc.scen.mat <- cbind(gc.scen.mat, gc.part.mat[,2])
  
  # Matrix of global inconsistency probabilites with fraction of simulation results (1 out of num.per.scen)
  gic.file <- paste("./", sim.version, "/glob_inconsis/ginc_", sim.version, "-", j, ".csv", sep = "")
  gic.part.mat <- read.csv(gic.file, header = TRUE)
  gic.scen.mat <- cbind(gic.scen.mat, gic.part.mat[,2])
  
  # After iterating through num.per.scen files, store results in lists and change scenario indicator
  if(j %% num.per.scen == 0){
    
    # Store results for global consistency in list
    rownames(gc.scen.mat) <- gc.part.mat[,1]                # list epsilon.star values as row names
    scen.lab <- paste("scenario_", which.scen, sep = "")
    gc.list[[scen.lab]] <- rowMeans(gc.scen.mat)            # save scenario-specific matrix of results
    
    # Store results for global inconsistency in list
    rownames(gic.scen.mat) <- gic.part.mat[,1]              # list epsilon.star values as row names
    scen.lab <- paste("scenario_", which.scen, sep = "")
    gic.list[[scen.lab]] <- rowMeans(gic.scen.mat)          # save scenario-specific matrix of results
    
    # Reset matrix and prepare for next scenario
    which.scen <- which.scen + 1                      # increase scenario indicator
    gc.scen.mat <- NULL
    gic.scen.mat <- NULL

  }
    
}


### Compile results for global (in)consistency probabilities

## Save results for different scenarios
if(sim.version == "n1508"){
  gc.n1508.1 <- gc.list[[1]]
  gc.n1508.2 <- gc.list[[2]]
  gc.n1508.3 <- gc.list[[3]]
  gc.n1508.4 <- gc.list[[4]]
  gc.n1508.5 <- gc.list[[5]]
  gc.n1508.6 <- gc.list[[6]]
}
if(sim.version == "n1508-beta20"){
  gc.n1508.beta20.1 <- gc.list[[1]]
  gc.n1508.beta20.2 <- gc.list[[2]]
  gc.n1508.beta20.3 <- gc.list[[3]]
  gc.n1508.beta20.4 <- gc.list[[4]]
  gc.n1508.beta20.5 <- gc.list[[5]]
  gc.n1508.beta20.6 <- gc.list[[6]]
}
if(sim.version == "n754"){
  gc.n754.1 <- gc.list[[1]]
  gc.n754.2 <- gc.list[[2]]
  gc.n754.3 <- gc.list[[3]]
  gc.n754.4 <- gc.list[[4]]
  gc.n754.5 <- gc.list[[5]]
  gc.n754.6 <- gc.list[[6]]
}
if(sim.version == "n754-beta20"){
  gc.n754.beta20.1 <- gc.list[[1]]
  gc.n754.beta20.2 <- gc.list[[2]]
  gc.n754.beta20.3 <- gc.list[[3]]
  gc.n754.beta20.4 <- gc.list[[4]]
  gc.n754.beta20.5 <- gc.list[[5]]
  gc.n754.beta20.6 <- gc.list[[6]]
}
if(sim.version == "n7650"){
  gc.n7650.1 <- gc.list[[1]]
  gc.n7650.2 <- gc.list[[2]]
  gc.n7650.3 <- gc.list[[3]]
  gc.n7650.4 <- gc.list[[4]]
  gc.n7650.5 <- gc.list[[5]]
  gc.n7650.6 <- gc.list[[6]]
}
if(sim.version == "n7650-beta20"){
  gc.n7650.beta20.1 <- gc.list[[1]]
  gc.n7650.beta20.2 <- gc.list[[2]]
  gc.n7650.beta20.3 <- gc.list[[3]]
  gc.n7650.beta20.4 <- gc.list[[4]]
  gc.n7650.beta20.5 <- gc.list[[5]]
  gc.n7650.beta20.6 <- gc.list[[6]]
}


### Plot global consistency probabilities for varying values of epsilon.star - proposed method #1

## GRAY SCALE

epsilon.star <- as.numeric( names(gc.list[[1]]) )
par( mar = c(3.9, 4.1, 2.0, 1.2) )
par( mfrow = c(2,3) )

# Global consistency plot for n = 754, beta = 0.5
plot( epsilon.star, gc.n754.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "black", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 754,  ", beta*"*", " = 0.5") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n754.2, lwd = 3, col = "gray60" )
lines( epsilon.star, gc.n754.3, lwd = 3, col = "gray85" )
lines( epsilon.star, gc.n754.4, lwd = 3, col = "gray85", lty = 2 )
lines( epsilon.star, gc.n754.5, lwd = 3, col = "gray60", lty = 2 )
lines( epsilon.star, gc.n754.6, lwd = 3, col = "black", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                     "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("black", "gray60", "gray85", "gray85", "gray60", "black"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

# Global consistency plot for n = 1508, beta = 0.5
plot( epsilon.star, gc.n1508.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "black", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 1508,  ", beta*"*", " = 0.5") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n1508.2, lwd = 3, col = "gray60" )
lines( epsilon.star, gc.n1508.3, lwd = 3, col = "gray85" )
lines( epsilon.star, gc.n1508.4, lwd = 3, col = "gray85", lty = 2 )
lines( epsilon.star, gc.n1508.5, lwd = 3, col = "gray60", lty = 2 )
lines( epsilon.star, gc.n1508.6, lwd = 3, col = "black", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                         "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("black", "gray60", "gray85", "gray85", "gray60", "black"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

# Global consistency plot for n = 7650, beta = 0.5
plot( epsilon.star, gc.n7650.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "black", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 7650,  ", beta*"*", " = 0.5") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n7650.2, lwd = 3, col = "gray60" )
lines( epsilon.star, gc.n7650.3, lwd = 3, col = "gray85" )
lines( epsilon.star, gc.n7650.4, lwd = 3, col = "gray85", lty = 2 )
lines( epsilon.star, gc.n7650.5, lwd = 3, col = "gray60", lty = 2 )
lines( epsilon.star, gc.n7650.6, lwd = 3, col = "black", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                         "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("black", "gray60", "gray85", "gray85", "gray60", "black"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

# Global consistency plot for n = 754, beta = 0.2
plot( epsilon.star, gc.n754.beta20.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "black", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 754,  ", beta*"*", " = 0.2") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n754.beta20.2, lwd = 3, col = "gray60" )
lines( epsilon.star, gc.n754.beta20.3, lwd = 3, col = "gray85" )
lines( epsilon.star, gc.n754.beta20.4, lwd = 3, col = "gray85", lty = 2 )
lines( epsilon.star, gc.n754.beta20.5, lwd = 3, col = "gray60", lty = 2 )
lines( epsilon.star, gc.n754.beta20.6, lwd = 3, col = "black", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                         "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("black", "gray60", "gray85", "gray85", "gray60", "black"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

# Global consistency plot for n = 1508, beta = 0.2
plot( epsilon.star, gc.n1508.beta20.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "black", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 1508,  ", beta*"*", " = 0.2") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n1508.beta20.2, lwd = 3, col = "gray60" )
lines( epsilon.star, gc.n1508.beta20.3, lwd = 3, col = "gray85" )
lines( epsilon.star, gc.n1508.beta20.4, lwd = 3, col = "gray85", lty = 2 )
lines( epsilon.star, gc.n1508.beta20.5, lwd = 3, col = "gray60", lty = 2 )
lines( epsilon.star, gc.n1508.beta20.6, lwd = 3, col = "black", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                         "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("black", "gray60", "gray85", "gray85", "gray60", "black"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

# Global consistency plot for n = 7650, beta = 0.2
plot( epsilon.star, gc.n7650.beta20.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "black", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 7650,  ", beta*"*", " = 0.2") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n7650.beta20.2, lwd = 3, col = "gray60" )
lines( epsilon.star, gc.n7650.beta20.3, lwd = 3, col = "gray85" )
lines( epsilon.star, gc.n7650.beta20.4, lwd = 3, col = "gray85", lty = 2 )
lines( epsilon.star, gc.n7650.beta20.5, lwd = 3, col = "gray60", lty = 2 )
lines( epsilon.star, gc.n7650.beta20.6, lwd = 3, col = "black", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                         "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("black", "gray60", "gray85", "gray85", "gray60", "black"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

par( mfrow = c(1,1) )




## COLOR SCALE

epsilon.star <- as.numeric( names(gc.list[[1]]) )
par( mar = c(3.9, 4.1, 2.0, 1.2) )
par( mfrow = c(2,3) )

# Global consistency plot for n = 754, beta = 0.5
plot( epsilon.star, gc.n754.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "midnightblue", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 754,  ", beta*"*", " = 0.5") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n754.2, lwd = 3, col = "cornflowerblue" )
lines( epsilon.star, gc.n754.3, lwd = 3, col = "lightblue" )
lines( epsilon.star, gc.n754.4, lwd = 3, col = "lightblue", lty = 2 )
lines( epsilon.star, gc.n754.5, lwd = 3, col = "cornflowerblue", lty = 2 )
lines( epsilon.star, gc.n754.6, lwd = 3, col = "midnightblue", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                         "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("midnightblue", "cornflowerblue", "lightblue",
                "lightblue", "cornflowerblue", "midnightblue"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

# Global consistency plot for n = 1508, beta = 0.5
plot( epsilon.star, gc.n1508.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "midnightblue", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 1508,  ", beta*"*", " = 0.5") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n1508.2, lwd = 3, col = "cornflowerblue" )
lines( epsilon.star, gc.n1508.3, lwd = 3, col = "lightblue" )
lines( epsilon.star, gc.n1508.4, lwd = 3, col = "lightblue", lty = 2 )
lines( epsilon.star, gc.n1508.5, lwd = 3, col = "cornflowerblue", lty = 2 )
lines( epsilon.star, gc.n1508.6, lwd = 3, col = "midnightblue", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                         "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("midnightblue", "cornflowerblue", "lightblue",
                "lightblue", "cornflowerblue", "midnightblue"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

# Global consistency plot for n = 7650, beta = 0.5
plot( epsilon.star, gc.n7650.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "midnightblue", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 7650,  ", beta*"*", " = 0.5") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n7650.2, lwd = 3, col = "cornflowerblue" )
lines( epsilon.star, gc.n7650.3, lwd = 3, col = "lightblue" )
lines( epsilon.star, gc.n7650.4, lwd = 3, col = "lightblue", lty = 2 )
lines( epsilon.star, gc.n7650.5, lwd = 3, col = "cornflowerblue", lty = 2 )
lines( epsilon.star, gc.n7650.6, lwd = 3, col = "midnightblue", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                         "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("midnightblue", "cornflowerblue", "lightblue",
                "lightblue", "cornflowerblue", "midnightblue"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

# Global consistency plot for n = 754, beta = 0.2
plot( epsilon.star, gc.n754.beta20.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "midnightblue", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 754,  ", beta*"*", " = 0.2") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n754.beta20.2, lwd = 3, col = "cornflowerblue" )
lines( epsilon.star, gc.n754.beta20.3, lwd = 3, col = "lightblue" )
lines( epsilon.star, gc.n754.beta20.4, lwd = 3, col = "lightblue", lty = 2 )
lines( epsilon.star, gc.n754.beta20.5, lwd = 3, col = "cornflowerblue", lty = 2 )
lines( epsilon.star, gc.n754.beta20.6, lwd = 3, col = "midnightblue", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                         "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("midnightblue", "cornflowerblue", "lightblue",
                "lightblue", "cornflowerblue", "midnightblue"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

# Global consistency plot for n = 1508, beta = 0.2
plot( epsilon.star, gc.n1508.beta20.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "midnightblue", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 1508,  ", beta*"*", " = 0.2") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n1508.beta20.2, lwd = 3, col = "cornflowerblue" )
lines( epsilon.star, gc.n1508.beta20.3, lwd = 3, col = "lightblue" )
lines( epsilon.star, gc.n1508.beta20.4, lwd = 3, col = "lightblue", lty = 2 )
lines( epsilon.star, gc.n1508.beta20.5, lwd = 3, col = "cornflowerblue", lty = 2 )
lines( epsilon.star, gc.n1508.beta20.6, lwd = 3, col = "midnightblue", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                         "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("midnightblue", "cornflowerblue", "lightblue",
                "lightblue", "cornflowerblue", "midnightblue"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

# Global consistency plot for n = 7650, beta = 0.2
plot( epsilon.star, gc.n7650.beta20.1, xlab = "", ylab = "", xlim = c(0,max(epsilon.star)),
      ylim = c(0,1), xaxs = "i", yaxs = "i", yaxt = "n", type = "l", lwd = 3,
      main = "",  col = "midnightblue", cex.axis = 1.2 )
axis( 2, at = c(.2, .4, .6, .8, 1), labels = TRUE, las = 1, cex.axis = 1.2 )
title( expression( paste(italic(N), " = 7650,  ", beta*"*", " = 0.2") ), line = 0.7,
       cex.main = 1.5 )
mtext( expression(epsilon), side = 1, line = 2.3, cex = 1.5 )
mtext( "Probability", side = 2, line = 2.8, cex = 1 )
lines( epsilon.star, gc.n7650.beta20.2, lwd = 3, col = "cornflowerblue" )
lines( epsilon.star, gc.n7650.beta20.3, lwd = 3, col = "lightblue" )
lines( epsilon.star, gc.n7650.beta20.4, lwd = 3, col = "lightblue", lty = 2 )
lines( epsilon.star, gc.n7650.beta20.5, lwd = 3, col = "cornflowerblue", lty = 2 )
lines( epsilon.star, gc.n7650.beta20.6, lwd = 3, col = "midnightblue", lty = 2 )
legend( "bottomright", c("0 Null Regions", "1 Null Region", "2 Null Regions",
                         "3 Null Regions", "4 Null Regions", "5 Null Regions"), bty = "n", lwd = 2,
        col = c("midnightblue", "cornflowerblue", "lightblue",
                "lightblue", "cornflowerblue", "midnightblue"),
        lty = c(1,1,1,2,2,2), cex = 1.1 )

par( mfrow = c(1,1) )

