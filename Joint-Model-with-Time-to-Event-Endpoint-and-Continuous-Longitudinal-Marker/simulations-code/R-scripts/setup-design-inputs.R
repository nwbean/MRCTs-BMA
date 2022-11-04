##################################################################################################
# SETUP DESIGN INPUTS
#
# GENERATE CSV FILE WITH FUNCTION INPUTS FOR SIMULATIONS
#
# FOR MOST CASES WITH S REGIONS, THE GENERATED CSV FILE OF DESIGN INPUTS
#   CONSIDERS S+1 SCENARIOS: 0 NULL REGIONS, 1 NULL REGION, . . . , S NULL REGIONS
##################################################################################################


# Load Rcpp code
library(Rcpp)
setwd("/Users/nathanbean/Documents/Dissertation/Project 3/Project3_Simulations/cluster-scripts")
sourceCpp("../R-source/bma-functions-joint-model.cpp")



# Simulation details
#sim.version <- "original_n9340_T0equal2"   # version of simulation
sim.version <- "equal_samp_n9340_T0equal2" # version of simulation
#sim.version <- "alt_half_n9340"     # version of simulation
#sim.version <- "null_half_n9340"    # version of simulation
#sim.version <- "survival-models-only" # version of simulation
#sim.version <- "full-joint-model"   # version of simulation
#sim.version <- "SA_n9340"           # version of simulation
#sim.id.number <- 93400              # seed for simulation version "original_n9340_T0equal2"
sim.id.number <- 93403              # seed for simulation version "equal_samp_n9340_T0equal2"
#sim.id.number <- 93404              # seed for simulation version "alt_half_n9340"
#sim.id.number <- 93405              # seed for simulation version "null_half_n9340"
#sim.id.number <- 93406              # seed for simulation version "survival-models-only"
#sim.id.number <- 93407              # seed for simulation version "full-joint-model"



#################################################################################
# S = 4 regions (e.g., Europe, North America, Asia, Rest of the World)
#
# Input designed to mimic phase 3 LEADER trial
# (ClinicalTrials.gov Identifier: NCT01179048)
#################################################################################

S <- 4                                                  # number of regions
#T0 <- 4                                                 # max number of distinct regions to consider
T0 <- 2                                                 # max number of distinct regions to consider
reg.names <- c("Reg1", "Reg2", "Reg3", "Reg4")          # region names in alphabetical order
#reg.names <- c("Asia", "Europe", "NorthAmerica", "RoW") # region names in alphabetical order (v. original)
N <- 9340                                               # total number of subjects (v. "n9340")
reg.allctn <- rep(1/S, S)                               # region sample size allocation (must sum to 1)
#reg.allctn <- c(711, 3296, 2847, 2486) / N              # region sample size allocation (must sum to 1, v. original)
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
hr.trtmt <- .868                    # hazard ratio of trtmt group with effect compared to control
hr.control <- 1                     # hazard ratio of trtmt group without effect compared to control
trtmt.hrs <- c(.62, .82, 1.01, .83)   # treatment hazard ratios for "original_n9430"
max.fu.time <- 60                   # maximum follow up time from beginning of study (months)
dropout.rate <- .0082               # dropout rate
sigma.e <- .886                     # sd (measurement error) in longitudinal model                      


# Determine number of rows of design inputs
#Number of unique scenarios, number  and random seeds for each scenario
#num.scenarios <- 2                              # number of unique scenarios for "original_n9340"
num.scenarios <- S + 1                          # number of unique scenarios
#num.scenarios <- 1                              # number of unique scenarios for "diff_effects_1-3"
num.total.ds <- 10000                           # number of total datasets to simulate per scenario
#num.total.ds <- 2000                           # number of total datasets for "n9340-equal-samp"
#num.total.ds <- 3000                           # number of total datasets for "n9340-equal-samp-T0equal2"
# Set num.ds.row = 20 for v. "original_n9340_T0equal2"
# Set num.ds.row = 20 for v. "equal_samp_n9340_T0equal2"
# Set num.ds.row = 20 for v. "alt_half_n9340"
# Set num.ds.row = 20 for v. "null_half_n9340"
# Set num.ds.row = 2500 for v. "survival-models-only"
# Set num.ds.row = 400 for v. "survival-models-only"
num.ds.row <- 20                     # number of datasets to simulate per row of design inputs
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
# Column 6: maximum follow up time from beginning of study (months)
# Column 7: dropout rate
# Column 8: likelihood standard deviation of longitudinal submodel (i.e., measurement error)
# Columns 9-(8+S): region names
# Columns (9+S)-(8+2*S): sample size allocation for each region
# Columns (9+2*S)-(8+3*S): hazard ratios of treatment group vs control group for each region
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
inputs.mat1.1[,6] <- max.fu.time
inputs.mat1.1[,7] <- dropout.rate
inputs.mat1.1[,8] <- sigma.e
for(i in 1:S){
  inputs.mat1.2[,i] <- reg.names[i]
  inputs.mat1.3[,i] <- reg.allctn[i]
  inputs.mat1.4.2[,i] <- ifelse(inputs.mat1.4.2[,i] == 1, hr.trtmt, inputs.mat1.4.2[,i] )
  inputs.mat1.4.2[,i] <- ifelse(inputs.mat1.4.2[,i] == 2, hr.control, inputs.mat1.4.2[,i] )
}
if( sim.version == "alt_half_n9340" ){
  # Adjust regional sample size allocation so alternative regions are half of null regions
  for(j in 1:num.scenarios){
    for(k in ((j-1)*num.repeat.rows + 1):(j*num.repeat.rows))
    inputs.mat1.3[k,] <- reg.allctn.alt.half[j,]
  }
}
if( sim.version == "null_half_n9340" ){
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
  if( sim.version %in% c("original_n9340", "original_n9340_T0equal2") ){
    if(i == 1){
      inputs.mat1.4[first.row:last.row,] <- rep.row( rep(hr.trtmt, 4), num.repeat.rows )
    } else{
      inputs.mat1.4[first.row:last.row,] <- rep.row( trtmt.hrs, num.repeat.rows )
    }
  } else if( sim.version == "diff_effects_1" |
      sim.version == "diff_effects_2" |
      sim.version == "diff_effects_3" ){
    inputs.mat1.4[first.row:last.row,] <- rep.row( trtmt.hrs, num.repeat.rows )
  } else{
    inputs.mat1.4[first.row:last.row,] <- rep.row(inputs.mat1.4.2[i,], num.repeat.rows)
  }
}


inputs.mat1 <- data.frame(inputs.mat1.1, inputs.mat1.2, inputs.mat1.3, inputs.mat1.4)
colnames(inputs.mat1) <- c( "scenario", "seed", "S", "T0", "N", "max.followup.time",
                            "dropout.rate", "sigma.e",
                            "reg1.name", "reg2.name", "reg3.name", "reg4.name",
                            "reg1.allctn", "reg2.allctn", "reg3.allctn", "reg4.allctn",
                            "reg1.trtmt.hr", "reg2.trtmt.hr", "reg3.trtmt.hr", "reg4.trtmt.hr" )


# Save as .csv file
file.name <- paste("./../cluster-input/input_mat_", sim.version, ".csv", sep="")
write.csv(inputs.mat1, file.name, quote = FALSE, row.names = FALSE)


