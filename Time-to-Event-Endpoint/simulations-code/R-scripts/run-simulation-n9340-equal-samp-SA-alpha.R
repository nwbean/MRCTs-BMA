#########################################################################################
# RUN SIMULATION
#
# Execute code from either local computer or from cluster
#########################################################################################


### Values to update for each simulation
### Simulation version and number of datasets to generate
sim.version <- "n9340-equal-samp-SA-alpha"
nsims <- 500


### Code for cluster
options(echo=TRUE)                             # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx <- as.numeric(args[1])


### Determine if code is being run on Mac or cluster
if(getwd() == "/Users/nathanbean/Documents/Dissertation/Project 2/Project2_Simulations/cluster-scripts" ||
   getwd() == "/Users/nathanbean"){
  node.idx <- 1
  root <- "/Users/nathanbean/Documents/Dissertation/Project 2/Project2_Simulations/cluster-scripts"
  setwd(root)
  library(simsurv)
} else{
  root <- ".";
  setwd(root)
  #install.packages("simsurv", repos="http://cran.rstudio.com/")
  library(simsurv, lib.loc = "../R-packages")
}



### Load libraries and source code
library(Rcpp)
library(rjags)
sourceCpp("../R-source/bma_functions_tte_endpoint.cpp")
source("../R-scripts/full-simulation-function-SA-alpha.R")



### Load inputs for simulation functions
input.file <- paste("../cluster-input/input_mat_SA_n9340.csv")
input.vals <- read.csv(input.file, header = TRUE)


## Scenario-specific simulation values
row.id <- which(input.vals$scenario == node.idx)          # extract row with inputs for given scenario
rand.seed <- input.vals[row.id, 2]                        # random seed


## Parameters for data generation
# Values mimic LEADER trial (ClinicalTrials.gov Identifier: NCT01179048)
S <- input.vals[row.id, 3]                                # number of regions
T0 <- input.vals[row.id, 4]                               # max number of distinct regions to consider
N <- input.vals[row.id, 5]                                # total number of subjects
reg.names <- input.vals[row.id, 10:(9+S)]                 # region names in alphabetical order
reg.allctn <- as.numeric(input.vals[row.id, (10+S):(9+2*S)]) # region sample size allocation (must sum to 1)
n.i <- floor(N * reg.allctn) +
  c(rep(1, N %% S), rep(0, S - (N %% S)))                 # regional sample size
max.time <- input.vals[row.id, 6]                         # maximum follow up time (from beginning of study)
max.enroll.time <- input.vals[row.id, 7]                  # maximum enrollment time
dropout.rate <- input.vals[row.id, 8]                     # dropout rate
base.haz <- input.vals[row.id, 9]                         # underlying constant baseline hazard
reg.hr <- as.numeric(input.vals[row.id, (10+2*S):(9+3*S)])   # region-specific hazard ratios
reg.te <- log(reg.hr)                                     # region-specific treatment effects
cov.effects <- NULL                                       # true covariate effects


## Covariate information
num.bin <- 0                                     # number of binary covariates
num.con <- 0                                     # number of normal covariates
cov.props <- 0                                   # success probabilities for binary covariates (0 if num.bin = 0)
cov.means <- 0                                   # means for normal covariates (0 if num.con = 0)
cov.sds <- 0                                     # standard deviations for normal covariates (0 if num.con = 0)
num.covs <- num.bin + num.con                    # number of covariates
cov.names <- paste( "Cov", 1:num.covs, sep = "") # labels for covariates


## Prior elicitation
K <- 8                                           # number of intervals for piecewise exponential model
hr.pred <- 1.3                                   # prediction of true hazard ratio
trt.pred <- log(hr.pred)                         # prediction of true treatment effect
cov.pred <- NULL                                 # predictions of true covariate effects
mu0 <- c( rep(trt.pred, S),
          rep(cov.pred, num.covs) )              # mean vector for prior on treatment and covariate effects
if(num.covs > 0){
  Sig0 <- diag( c( rep((10 * abs(trt.pred))^2, S),
                   (10 * abs(cov.pred))^2 ) )    # covariance matrix for prior on treatment and covariate effect
} else{
  Sig0 <- diag( rep((10 * abs(trt.pred))^2, S) ) # covariance matrix for prior on treatment and covariate effect
}
eta0 <- matrix( 0.01, nrow = S, ncol = K )       # shape hyperparameters for prior on baseline hazards
phi0 <- matrix( 0.01, nrow = S, ncol = K )       # rate hyperparameters for prior on baseline hazards
alpha0 <- c(-5, -2, -1, -.5, 0, .5, 1, 2, 5)     # pre-specified values for prior model probabilities


## Values for BMA
gamma0 <- 0                         # value for posterior prob. Pr(gamma > gamma0|D)
n.draws <- 10000                    # number of posterior samples to draw when evaluating consistency
pi0 <- 1/S                          # value for post. prob. Pr(gamma_i/gamma > pi0|D) - (Japanese MHLW consistency)
epsilon.star <- 1/hr.pred           # minimal clinically important regional difference, used in most measures of consistency
beta.star <- .5                     # probability cutoff for which consider two regions to be clinically different,
#                                   #   used in Pr(|gamma_i - gamma_j| > -log(epsilon.star)|D) >= beta.star
#                                   #   (used for global inconsistency)


## Tuning parameters for competitor Cox proportional hazards models
sd.CPHM1 <- .2       # used in Metropolis-Hastings algorithm for CPH1 model
sd.CPHM2 <- .13      # used in Metropolis-Hastings algorithm for CPH2 model



### Run simulations
set.seed(rand.seed)
start.sim.time <- Sys.time()
sim.results <- run.sims( S = S, T0 = T0, n.i = n.i, reg.names = reg.names, K = K, max.time = max.time,
                         max.enroll.time = max.enroll.time, dropout.rate = dropout.rate,
                         base.haz = base.haz, reg.te = reg.te, num.covs = num.covs,
                         cov.effects = cov.effects, cov.names = cov.names, cov.props = cov.props,
                         cov.means = cov.means, cov.sds = cov.sds, alpha0 = alpha0, mu0 = mu0,
                         Sig0 = Sig0, eta0 = eta0, phi0 = phi0, gamma0 = gamma0, n.draws = n.draws,
                         epsilon.star = epsilon.star, beta.star = beta.star, pi0 = pi0,
                         sd.CPHM1 = sd.CPHM1, sd.CPHM2 = sd.CPHM2, nsims = nsims, print.iters = TRUE )
end.sim.time <- Sys.time()
total.sim.time <- end.sim.time - start.sim.time
total.sim.time



### Save results (rejection rate, PMPs, estimates, MSE, bias, consistency measures) as .csv files
dir.create( paste("../cluster-results/", sim.version, sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/rejection_rates", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/bias", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/mse", sep=""), showWarnings = FALSE )
output.file.rr <- paste("../cluster-results/", sim.version, "/rejection_rates/rr_", sim.version, "-", row.id, ".csv", sep="")
output.file.bias.rte <- paste("../cluster-results/", sim.version, "/bias/rte_bias_", sim.version, "-", row.id, ".csv", sep="")
output.file.mse.rte <- paste("../cluster-results/", sim.version, "/mse/rte_mse_", sim.version, "-", row.id, ".csv", sep="")
write.csv(sim.results$Rejection.Rate, output.file.rr, quote = FALSE)
write.csv(sim.results$Reg.Trtmt.Effect.Bias, output.file.bias.rte, quote = FALSE)
write.csv(sim.results$Reg.Trtmt.Effect.MSE, output.file.mse.rte, quote = FALSE)

