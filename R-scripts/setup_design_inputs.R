##################################################################################################
# SETUP DESIGN INPUTS
#
# GENERATE CSV FILE WITH FUNCTION INPUTS FOR SIMULATIONS
#
# FOR S REGIONS, THE GENERATED CSV FILE OF DESIGN INPUTS CONSIDERS S+1 SCENARIOS:
#   0 NULL REGIONS, 1 NULL REGION, . . . , S NULL REGIONS
##################################################################################################


# Load Rcpp code
library(Rcpp)
setwd("/Users/nathanbean/Documents/Dissertation/Project 1/Project1_Simulations_v2/cluster-scripts")
sourceCpp("../R-source/bma_functions_rcpp_code.cpp")



# Simulation details
sim.version <- "n1508"              # version of simulation
#sim.id.number <- ???                # seed for simulation version "n1508"
#sim.id.number <- 754                # seed for simulation version "n754"
#sim.id.number <- 7650               # seed for simulation version "n7650"
#sim.id.number <- 15081              # seed for simulation version "alt_half_n1508"
#sim.id.number <- 15082              # seed for simulation version "null_half_n1508"
#sim.id.number <- 15083              # seed for simulation version "diff_effects_1"
#sim.id.number <- 15084              # seed for simulation version "diff_effects_2"
#sim.id.number <- 15085              # seed for simulation version "diff_effects_3"



##############################################################################
# S = 5 regions (Asia Pacific, Eastern Europe, Other, US, Western Europe)
#
# Input designed to mimic GlaxoSmithKline COPD study
# (ClinicalTrials.gov Identifier: NCT02105974)
##############################################################################

S <- 5                                                  # number of regions
T0 <- 5                                                 # max number of distinct regions to consider
mod.mat <- compile_modMat(S, T0)                        # matrix with all possible models in model space
num.mods <- num_models(S, T0)                           # number of models in model space
reg.names <- c("Reg1", "Reg2", "Reg3", "Reg4", "Reg5")  # region names in alphabetical order
N <- 1508                                               # total number of subjects (v. "n1508")
#N <- 754                                                # total number of subjects (v. "n754")
#N <- 7650                                               # total number of subjects (v. "n7650")
sd.Y <- 0.205                                           # common st. dev. for both groups and all regions
reg.allctn <- rep(1/S, S)                               # region sample size allocation (must sum to 1)
reg.allctn.alt.half <- matrix(0, nrow = S+1, ncol = S)  # sample size allocation for alt. half of null
reg.allctn.null.half <- matrix(0, nrow = S+1, ncol = S) # sample size allocation for null half of alt.
reg.allctn.alt.half[c(1,S+1),] <- reg.allctn
reg.allctn.null.half[c(1,S+1),] <- reg.allctn
for(i in 1:(S-1)){
  reg.allctn.alt.half[i+1,] <- round( c( rep( 2/(S+i), i), rep( 1/(S+i), S-i) ), digits = 5 )
  reg.allctn.null.half[i+1,] <- round( c( rep( 1/(2*S-i), i), rep( 2/(2*S-i), S-i) ), digits = 5 )
}
reg.allctn.alt.half[,S] <- 1 - rowSums(reg.allctn.alt.half[,1:(S-1)])    # Force to sum to exactly 1.00
reg.allctn.null.half[,S] <- 1 - rowSums(reg.allctn.null.half[,1:(S-1)])  # Force to sum to exactly 1.00
for(i in 1:(S-1)){
  reg.allctn.alt.half <- rbind( 
    reg.allctn.alt.half,
    round( c( rep( 2/(S+i), i), rep( 1/(S+i), S-i) ), digits = 5 ) )
  reg.allctn.null.half <- rbind( 
    reg.allctn.null.half,
    round( c( rep( 1/(2*S-i), i), rep( 2/(2*S-i), S-i) ), digits = 5 ) )
}


reg.allctn.alt.half <- rbind( reg.allctn.alt.half, reg.allctn )
reg.allctn.null.half <- rbind( reg.allctn.null.half, reg.allctn )
cntrl.mean <- .082                                      # mean of control group for all regions
trtmt.mean <- .116                                      # mean of treatment groups
#trtmt.means <- c(.017, .026, .034, .043, .051) + .082   # treatment group means for "diff_effects_1"
#trtmt.means <- c(.017, .017, .017, .034, .034) + .082   # treatment group means for "diff_effects_2"
#trtmt.means <- c(.017, .034, .034, .051, .051 ) + .082  # treatment group means for "diff_effects_3"


# Determine number of rows of design inputs
#Number of unique scenarios, number  and random seeds for each scenario
num.scenarios <- S + 1                          # number of unique scenarios
#num.scenarios <- 1                              # number of unique scenarios for "diff_effects_1-3"
num.total.ds <- 10000                           # number of total datasets to simulate per scenario
# num.ds.row = 1000 for v. "n1508"
# num.ds.row = 1000 for v. "n754"
# num.ds.row = 50 for v. "n7650"
# num.ds.row = 500 for v. "diff_effects_1-3"
num.ds.row <- 1000                              # number of datasets to simulate per row of design inputs
num.repeat.rows <- ceiling(num.total.ds / num.ds.row)     # number of repeated rows per scenario
total.num.rows <- num.repeat.rows * num.scenarios         # total number of rows for design inputs
set.seed(sim.id.number)
scenario.seeds <- sample(1:(2^20), total.num.rows, replace = FALSE)


### Create matrix of design inputs
# Column 1: row indicator for unique design inputs
# Column 2: seed number
# Column 3: number of regions (S)
# Column 4: max number of distinct regions to consider (T0)
# Column 5: total number of subjects (N)
# Column 6: proportion of subjects allocated to treatment (1 minus this proportion allocated to control)
# Column 7: control mean
# Column 8: common standard deviation for both groups and all regions
# Columns 9-(8+S): region names
# Columns (9+S)-(8+2*S): sample size allocation for each region
# Columns (9+2*S)-(8+3*S): treatment group means for each region
inputs.mat1.1 <- matrix(0, nrow = total.num.rows, ncol = 8)
inputs.mat1.2 <- matrix(0, nrow = total.num.rows, ncol = S)
inputs.mat1.3 <- matrix(0, nrow = total.num.rows, ncol = S)
inputs.mat1.4 <- matrix(0, nrow = total.num.rows, ncol = S)
inputs.mat1.4.1 <- matrix(1, nrow = S, ncol = S)
inputs.mat1.4.1[lower.tri(inputs.mat1.4.1)] <- 2
inputs.mat1.4.2 <- rbind(inputs.mat1.4.1, rep(2, S))

inputs.mat1.1[,1] <- 1:total.num.rows
inputs.mat1.1[,2] <- scenario.seeds
inputs.mat1.1[,3] <- S
inputs.mat1.1[,4] <- T0
inputs.mat1.1[,5] <- N
inputs.mat1.1[,6] <- .5
inputs.mat1.1[,7] <- cntrl.mean
inputs.mat1.1[,8] <- sd.Y
for(i in 1:S){
  inputs.mat1.2[,i] <- reg.names[i]
  inputs.mat1.3[,i] <- reg.allctn[i]
  inputs.mat1.4.2[,i] <- ifelse(inputs.mat1.4.2[,i] == 1, trtmt.mean, inputs.mat1.4.2[,i] )
  inputs.mat1.4.2[,i] <- ifelse(inputs.mat1.4.2[,i] == 2, cntrl.mean, inputs.mat1.4.2[,i] )
}
if( sim.version == "alt_half_n1508" ){
  # Adjust regional sample size allocation so alternative regions are half of null regions
  for(j in 1:num.scenarios){
    for(k in ((j-1)*num.repeat.rows + 1):(j*num.repeat.rows))
    inputs.mat1.3[k,] <- reg.allctn.alt.half[j,]
  }
}
if( sim.version == "null_half_n1508" ){
  # Adjust regional sample size allocation so null regions are half of alternative regions
  for(j in 1:num.scenarios){
    for(k in ((j-1)*num.repeat.rows + 1):(j*num.repeat.rows))
      inputs.mat1.3[k,] <- reg.allctn.null.half[j,]
  }
}
rep.row <- function(x, n){
  matrix( rep(x, each = n), nrow=n )
}
for(i in 1:num.scenarios){
  first.row <- (i-1)*num.repeat.rows + 1
  last.row <- i*num.repeat.rows
  if( sim.version == "diff_effects_1" |
      sim.version == "diff_effects_2" |
      sim.version == "diff_effects_3" ){
    inputs.mat1.4[first.row:last.row,] <- rep.row( trtmt.means, num.repeat.rows )
  } else{
    inputs.mat1.4[first.row:last.row,] <- rep.row(inputs.mat1.4.2[i,], num.repeat.rows)
  }
}


inputs.mat1 <- data.frame(inputs.mat1.1, inputs.mat1.2, inputs.mat1.3, inputs.mat1.4)
colnames(inputs.mat1) <- c( "scenario", "seed", "S", "T0", "N", "trtmt.allctn", "cntrl.mean", "sd.Y",
                            "reg1.name", "reg2.name", "reg3.name", "reg4.name", "reg5.name",
                            "reg1.allctn", "reg2.allctn", "reg3.allctn", "reg4.allctn",
                            "reg5.allctn", "reg1.trtmt.mean", "reg2.trtmt.mean", "reg3.trtmt.mean",
                            "reg4.trtmt.mean", "reg5.trtmt.mean" )

# Save as .csv file
file.name <- paste("./../cluster-input/input_mat_", sim.version, ".csv", sep="")
write.csv(inputs.mat1, file.name, quote = FALSE, row.names = FALSE)


