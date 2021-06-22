#########################################################################################
# RUN SIMULATION
#
# Execute code from either local computer or from cluster
#########################################################################################


### Values to update for each simulation
### Simulation version and number of datasets to generate
sim.version <- "SA8"   # double control intercept, true mean difference
num.sims <- 1000


### Code for cluster
options(echo=TRUE)                             # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx <- as.numeric(args[1])



### Determine if code is being run on Mac or cluster
if(getwd() == "/Users/nathanbean/Documents/Dissertation/Project 1/Project1_Simulations_v2/cluster-scripts" ||
   getwd() == "/Users/nathanbean"){
  node.idx <- 1
  root <- "/Users/nathanbean/Documents/Dissertation/Project 1/Project1_Simulations_v2/cluster-scripts"
  setwd(root)
} else{
  root <- ".";
  setwd(root)
} 



### Load libraries and source Rcpp/R code
library(Rcpp)
library(rjags)
sourceCpp("../R-source/bma_functions_rcpp_code.cpp")
source("../R-scripts/full_simulation_function.R")



### Load inputs for simulation functions
input.file <- paste("../cluster-input/input_mat_equal_samp_n1508.csv")
input.vals <- read.csv(input.file, header = TRUE)


## Scenario-specific simulation values
row.id <- which(input.vals$scenario == node.idx)          # extract row with inputs for given scenario
rand.seed <- input.vals[row.id, 2]                        # random seed


## Parameters for data generation
S <- input.vals[row.id, 3]                                # number of regions
T0 <- input.vals[row.id, 4]                               # max number of distinct regions to consider
N <- input.vals[row.id, 5]                                # total number of subjects
trtmt.allctn <- c(input.vals[row.id, 6],
                  1 - input.vals[row.id, 6])              # treatment group allocation (must sum to 1)
cntrl.means <- rep(input.vals[row.id, 7], S)              # regional means for control groups
sd.Y <- input.vals[row.id, 8]                             # common st. dev. for both groups and all regions
reg.names <- numeric(S)                                   # region names in alphabetical order
reg.allctn <- numeric(S)                                  # region sample size allocation (must sum to 1)
trtmt.means <- numeric(S)                                 # regional means for treatment groups
for(i in 1:S){
  reg.names[i] <- as.character(input.vals[row.id, 8+i])
  reg.allctn[i] <- input.vals[row.id, 8+S+i]
  trtmt.means[i] <- input.vals[row.id, 8+2*S+i]
}


## Covariate information
num.bin <- 0                                    # number of binary covariates
num.con <- 0                                    # number of normal covariates
cov.props <- 0                                  # success probabilites for binary covariates (0 if num.bin = 0)
cov.means <- 0                                  # means for normal covariates (0 if num.con = 0)
cov.sds <- 0                                    # standard deviations for normal covariates (0 if num.con = 0)


## Prior specification
mod.prior.alpha <- 0                            # determines type of prior assumption for models
mean.diff <- .034                               # assumed mean difference (treatment mean minus control mean)
cntrl.int <- .164                               # assumed intercept for control group (approximately 20% increase from true value)
m0 <- construct_m0(S, num.bin, num.con, reg_int_mean = cntrl.int, reg_trt_mean = 0, cov_mean = 0)
int.mltplr <- 10                                # amount to multiply control intercept by in diagonal of Sig0
Sig0.int.var <- (cntrl.int * int.mltplr)^2      # set cntrl.int * int.mltplr as standard deviation
trtmt.mltplr <- 10                              # amount to multiply mean difference by in diagonal of Sig0
Sig0.var <- (mean.diff * trtmt.mltplr)^2        # set mean.diff * trtmt.mltplr as standard deviation
Sig0 <- construct_Sig0(S, num.bin, num.con, reg_int_diag = Sig0.int.var, reg_trt_diag = Sig0.var,
                       cov_diag = Sig0.var, reg_int_offdiag = 0, reg_trt_offdiag = 0,
                       cov_offdiag = 0, ri_rt_offdiag = 0, ri_cov_offdiag = 0, rt_cov_offdiag = 0)
delta0 <- .001
nu0 <- .001
modPriors <- construct_modPriors(S, T0, alpha = mod.prior.alpha)


## Values for BMA
n.draws <- 100000                   # number of draws from each posterior dist.
gamma0 <- 0                         # value for posterior prob. Pr(gamma > gamma0|D)
pi0 <- 1/S                          # value for post. prob. Pr(gamma_i/gamma > pi0|D) - (Japanese MHLW consistency)
epsilon.star <- 0.018               # minimal clinically important difference, used in most measures of consistency
beta.star <- .5                     # probability (1 - beta_star) for which consider two regions to be clinically different,
#                                   #   used in Pr(|gamma_i - gamma_j| > epsilon.star|D) >= (1-beta.star)
#                                   #   (used for global inconsistency)



### Run simulations
set.seed(rand.seed)
start.sim.time <- Sys.time()
sim.results <- run.sims(N, S, T0, reg.names, reg.allctn, trtmt.allctn, trtmt.means,
                        cntrl.means, sd.Y, num.bin, num.con, cov.props, cov.means,
                        cov.sds, modPriors, m0, Sig0, delta0, nu0, gamma0, epsilon.star,
                        beta.star, pi0, n.draws, num.sims, print.iters = TRUE)
end.sim.time <- Sys.time()
total.sim.time <- end.sim.time - start.sim.time
total.sim.time



### Print results in console
sim.details <- data.frame( Regions = reg.names,
                           Sample.Allctn = reg.allctn,
                           Cntrl.Means = cntrl.means,
                           Trtmt.Means = trtmt.means )
cat(
  cat("Simulation #", node.idx, sep = ""), "\n",
  "Total simulation time:", round(total.sim.time, 3), "\n"
)

sim.details

cat("Rejection rates (true positive rate / false positive rate):")
sim.results$Rejection.Rate

cat("Average PMP for each model:")
sim.results$Avg.PMP

cat("Average point estimates (posterior means) for regional treatment effects:")
sim.results$Avg.Trtmt.Effect.Estimates

cat("Bias of regional treatment effect estimates for each method:")
sim.results$Reg.Trtmt.Effect.Bias

cat("MSE of regional treatment effect estimates for each method:")
sim.results$Reg.Trtmt.Effect.MSE

cat("Average point estimates (posterior means) for regional intercepts:")
sim.results$Avg.Reg.Intercept.Estimates

cat("Bias of regional intercept estimates for each method:")
sim.results$Reg.Intercept.Bias

cat("MSE of regional intercept estimates for each method:")
sim.results$Reg.Intercept.MSE



### Save results (rejection rate, average model PMPs, bias, MSE) as .csv files
dir.create( paste("../cluster-results/", sim.version, sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/rejection_rates", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/pmp", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/avg_pmp", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/bias", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/mse", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/rte", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/gte", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/loc_consis", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/prws_consis", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/glob_consis", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/glob_inconsis", sep=""), showWarnings = FALSE )
output.file.rr <- paste("../cluster-results/", sim.version, "/rejection_rates/rr_", sim.version, "-", row.id, ".csv", sep="")
output.file.pmp <- paste("../cluster-results/", sim.version, "/pmp/pmp_", sim.version, "-", row.id, ".csv", sep="")
output.file.avg.pmp <- paste("../cluster-results/", sim.version, "/avg_pmp/avg_pmp_", sim.version, "-", row.id, ".csv", sep="")
output.file.bias.rte <- paste("../cluster-results/", sim.version, "/bias/rte_bias_", sim.version, "-", row.id, ".csv", sep="")
output.file.mse.rte <- paste("../cluster-results/", sim.version, "/mse/rte_mse_", sim.version, "-", row.id, ".csv", sep="")
output.file.bias.ri <- paste("../cluster-results/", sim.version, "/bias/ri_bias_", sim.version, "-", row.id, ".csv", sep="")
output.file.mse.ri <- paste("../cluster-results/", sim.version, "/mse/ri_mse_", sim.version, "-", row.id, ".csv", sep="")
output.file.rte <- paste("../cluster-results/", sim.version, "/rte/rte_", sim.version, "-", row.id, ".csv", sep="")
output.file.gte <- paste("../cluster-results/", sim.version, "/gte/gte_", sim.version, "-", row.id, ".csv", sep="")
output.file.lc.MHLW <- paste("../cluster-results/", sim.version, "/loc_consis/lc_MHLW_", sim.version, "-", row.id, ".csv", sep="")
output.file.lc.loo <- paste("../cluster-results/", sim.version, "/loc_consis/lc_loo_", sim.version, "-", row.id, ".csv", sep="")
output.file.pc <- paste("../cluster-results/", sim.version, "/prws_consis/pc_", sim.version, "-", row.id, ".csv", sep="")
output.file.gc <- paste("../cluster-results/", sim.version, "/glob_consis/gc_", sim.version, "-", row.id, ".csv", sep="")
output.file.ginc <- paste("../cluster-results/", sim.version, "/glob_inconsis/ginc_", sim.version, "-", row.id, ".csv", sep="")
write.csv(sim.results$Rejection.Rate, output.file.rr, quote = FALSE)
write.csv(sim.results$Model.PMP.Values, output.file.pmp, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Avg.PMP, output.file.avg.pmp, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Reg.Trtmt.Effect.Bias, output.file.bias.rte, quote = FALSE)
write.csv(sim.results$Reg.Trtmt.Effect.MSE, output.file.mse.rte, quote = FALSE)
write.csv(sim.results$Reg.Intercept.Bias, output.file.bias.ri, quote = FALSE)
write.csv(sim.results$Reg.Intercept.MSE, output.file.mse.ri, quote = FALSE)
write.csv(sim.results$Regional.Trtmt.Effects, output.file.rte, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Global.Trtmt.Effect, output.file.gte, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Local.Consistency.MHLW, output.file.lc.MHLW, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Local.Consistency.LOO, output.file.lc.loo, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Pairwise.Consistency, output.file.pc, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Global.Consistency, output.file.gc, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Global.Inconsistency, output.file.ginc, quote = FALSE, row.names = FALSE)

