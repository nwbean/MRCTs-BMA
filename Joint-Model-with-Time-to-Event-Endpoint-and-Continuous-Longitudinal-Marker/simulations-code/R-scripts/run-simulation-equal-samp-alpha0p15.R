#########################################################################################
# RUN SIMULATION
#
# Execute code from either local computer or from cluster
#########################################################################################


### Values to update for each simulation
### Simulation version and number of datasets to generate
sim.version <- "equal-samp-alpha0p15"
nsims <- 20


### Code for cluster
options(echo=TRUE)                             # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx <- as.numeric(args[1])


### Determine if code is being run on Mac or cluster
if(getwd() == "/Users/nathanbean/Documents/Dissertation/Project 3/Project3_Simulations/cluster-scripts" ||
   getwd() == "/Users/nathanbean"){
  node.idx <- 1
  root <- "/Users/nathanbean/Documents/Dissertation/Project 3/Project3_Simulations/cluster-scripts"
  setwd(root)
} else{
  root <- ".";
  setwd(root)
}



### Load libraries and source code
library(Rcpp)
library(rjags)
sourceCpp("../R-source/bma-functions-joint-model.cpp")
source("../R-scripts/full-simulation-function.R")



### Details for data generation

## Load inputs for simulation functions
input.file <- paste("../cluster-input/input_mat_equal_samp_n9340_T0equal2.csv")
input.vals <- read.csv(input.file, header = TRUE)


## Scenario-specific simulation values
row.id <- which(input.vals$scenario == node.idx)          # extract row with inputs for given scenario
rand.seed <- input.vals[row.id, 2]                        # random seed


## Parameters for data generation
# Values mimic LEADER trial (ClinicalTrials.gov Identifier: NCT01179048)
S <- input.vals[row.id, 3]                       # number of regions
T0 <- input.vals[row.id, 4]                      # max number of distinct regions to consider
N <- input.vals[row.id, 5]                       # total number of subjects
reg.names <- input.vals[row.id, 9:(8+S)]         # region names in alphabetical order
reg.allctn <- as.numeric(input.vals[row.id, (9+S):(8+2*S)]) # region sample size allocation (must sum to 1)
n.i <- floor(N * reg.allctn) +
  c(rep(1, N %% S), rep(0, S - (N %% S)))        # regional sample size
max.time <- input.vals[row.id, 6]                # maximum follow up time (months, from beginning of study)
dropout.rate <- input.vals[row.id, 7]            # dropout rate
reg.hr <- as.numeric(input.vals[row.id, (9+2*S):(8+3*S)])   # region-specific hazard ratios
reg.te.surv <- log(reg.hr)                       # region-specific treatment effects (survival model)
cov.long.effects <- NULL                         # true covariate effects for long. model
cov.surv.effects <- NULL                         # true covariate effects for survival model
sigma.e <- input.vals[row.id, 8]                 # sd (measurement error) in longitudinal model


## True association parameter connecting two submodels via shared subject-specific random intercepts
alpha.vec <- 0.15


## Underlying piecewise constant baseline hazards
pw.hazs <- seq(.0026, .0045, .0001)
pw.times <- seq(0, 60, length.out = length(pw.hazs) + 1)  # left bounds of intervals for piecewise baseline hazards
pw.times <- pw.times[ !(pw.times %in% max.time) ]         # remove maximum time (max. time isn't a left bound)        


## Covariate information
num.bin.long <- 0       # number of binary covariates for long. model
num.bin.surv <- 0       # number of binary covariates for survival model
num.con.long <- 0       # number of normal covariates for long. model
num.con.surv <- 0       # number of normal covariates for survival model
cov.long.props <- 0     # success probabilities for long. binary covariates (0 if num.bin.long = 0)
cov.surv.props <- 0     # success probabilities for survival binary covariates (0 if num.bin.surv = 0)
cov.long.means <- 0     # means for long. normal covariates (0 if num.con.long = 0)
cov.surv.means <- 0     # means for survival normal covariates (0 if num.con.surv = 0)
cov.long.sds <- 0       # standard deviations for long. normal covariates (0 if num.con.long = 0)
cov.surv.sds <- 0       # standard deviations for survival normal covariates (0 if num.con.surv = 0)
num.covs.long <- num.bin.long + num.con.long      # number of covariates for long. model
num.covs.surv <- num.bin.surv + num.con.surv      # number of covariates for survival model
if(num.covs.long > 0){
  cov.long.names <- paste0( "Cov.long.", 1:num.covs.long)  # names of covariates in long. model
} else{
  cov.long.names <- NULL
}
if(num.covs.surv > 0){
  cov.surv.names <- paste0( "Cov.surv.", 1:num.covs.surv)  # names of covariates in survival model
} else{
  cov.surv.names <- NULL
}


## Visit times (months) and visit-specific mean measurements for treatment and control group
vst.times <- c(0, 3, 6, 12, 18, 24, 30, 36, 42, 48, 54)
vst.means.0 <- c(8.7, 8.2, 8.1, 8.0, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 8.0)
vst.means.1 <- c(8.7, 7.2, 7.2, 7.3, 7.3, 7.4, 7.4, 7.5, 7.5, 7.5, 7.6)


## Distributional values of subject-specific random intercepts
REs <- "int"                     # type of random effects to include ("int" or "int-slope")
if(REs == "int"){
  b.mean <- 0                    # mean
  b.sd <- 1.057                  # standard deviation
} else if(REs == "int-slope"){
  b.mean <- c(0, 0)              # mean
  #b.sd <- c(1.057, .03)          # standard deviation
}



### Model specifications for longitudinal and survival submodels

## Number of time intervals (surv.) and location of knots for linear splines (long.)
Q <- 8                           # number of time intervals with separate baseline hazards
lin.spline.knots <- c(3, 18)     # locations of knots for linear splines


## Prior elicitation for parameters in survival submodel
hr.pred.surv <- 1.3                      # prediction of true hazard ratio (surv.)
trt.pred.surv <- log(hr.pred.surv)       # prediction of true treatment effect (surv.)
cov.pred.surv <- NULL                    # predictions of true covariate effects (surv.)
mu.Y <- c( rep(trt.pred.surv, S),        # mean vector for prior on trt and covariate effects (surv.)
          rep(cov.pred.surv, num.covs.surv) )
Sig.Y <- diag( (10 * abs(mu.Y))^2 )      # covariance matrix for prior on trt and covariate effects (surv.)
eta.iq <- matrix( 0.01, nrow = S, ncol = Q )   # shape hyperparameters for prior on baseline hazards (surv.)
phi.iq <- matrix( 0.01, nrow = S, ncol = Q )   # rate hyperparameters for prior on baseline hazards (surv.)
if(REs == "int"){
  mu.alpha <- as.matrix(0)               # prior mean for association parameter
  Sig.alpha <- as.matrix(1000)           # prior variance for association parameter
} else if(REs == "int-slope"){
  mu.alpha <- matrix(0, nrow = 2, ncol = 1)   # prior mean vector for association parameters
  Sig.alpha <- 1000 * diag(2)            # prior covariance matrix for association parameters
}


## Prior elicitation for parameters in longitudinal submodel
prior.int.mean.long <- 8.6               # prior mean for all regression effects (long.)
prior.reg.mean.long <- -0.25             # prior mean for all regression effects (long.)
num.effects.long <- 2*S + 1 + length(lin.spline.knots) +
  S*length(lin.spline.knots) + num.covs.long     # number of regression effects minus S intercepts (long.)
mu.X <- c( rep(prior.int.mean.long, S),  # mean vector for prior on regression effects (long.)
           rep(prior.reg.mean.long, num.effects.long) )
Sig.X <- diag( (10 * abs(mu.X))^2 )      # cov. matrix for prior on regression effects (long.)
eta.tau <- .001            # shape for gamma prior on precision of likelihood error (long.)
phi.tau <- .001            # rate for gamma prior on precision of likelihood error (long.)
if(REs == "int"){
  nu.G <- 1                # degrees of freedom for Wishart prior on G^-1 (precision matrix of REs) (long.)
  C0.G <- as.matrix(1)  # scale matrix in Wishart prior on G^-1
} else if(REs == "int-slope"){
  nu.G <- 2                # degrees of freedom for Wishart prior on G^-1 (precision matrix of REs) (long.)
  C0.G <- diag(2)          # scale matrix in Wishart prior on G^-1
}
C0.inv.G <- solve(C0.G)          # inverse of scale matrix in Wishart prior on G^-1


## Prior elicitation for models in model space
a0.surv <- 0                     # pre-specified value for prior model probabilities (surv.)
a0.long <- 0                     # pre-specified value for prior model probabilities (long.)
# prior model probabilities \propto exp(D.surv.l*a0.surv + D.long.l*a0.long),
#    where D.surv.l and D.long.l are the number of distinct trt effects for the l^th
#    survival and longitudinal submodels, respectively.
# NOTE: a vector of custom prior model probabilities can be read into the BMA function
#    using the function argument "mod.priors" In this case, the function arguments
#    "a0.surv" and "a0.long" should not be used.



### Values for BMA
max.iter <- 10       # maximum number of iterations in optimization algorithm for each model
optim.toler.params <- .005   # tolerance for parameter convergence in optimization algorithm
optim.toler.b <- .05         # tolerance for random effects convergence in optimization algorithm
gamma0.surv <- 0     # value for posterior prob. Pr(gamma_G > gamma0|D) (survival model)
gamma0.long <- 0     # value for posterior prob. Pr(gamma_G > gamma0|D) (longitudinal model)
n.draws <- 10000     # number of posterior samples to draw when evaluating consistency
pi0 <- 1/S           # value for post. prob. Pr(gamma_i/gamma > piO|D) - (Japanese MHLW consistency)
epsilon.star <- 1/hr.pred.surv
#                    # range of possible MCIRDS, used in most measures of consistency
beta.star <- .5      # probability cutoff for which consider two regions to be clinically different,
#                    # used in Pr(|gamma_i - gamma_j| > -log(epsilon.star)|D) >= beta.star
#                    # (used for global inconsistency)



### Run simulations
set.seed(rand.seed)
start.sim.time <- Sys.time()
sim.results <- run.sims(
  # Data simulation inputs
  S = S, T0 = T0, n.i = n.i, reg.names = reg.names, pw.hazs = pw.hazs,
  pw.times = pw.times, reg.te.surv = reg.te.surv, sigma.e = sigma.e,
  b.mean = b.mean, b.sd = b.sd, alpha.vec = alpha.vec, vst.times = vst.times,
  vst.means.0 = vst.means.0, vst.means.1 = vst.means.1, max.time = max.time,
  dropout.rate = dropout.rate,
  # Covariate details
  num.covs.long = num.covs.long, num.covs.surv = num.covs.surv,
  num.bin.long = num.bin.long, num.bin.surv = num.bin.surv,
  num.con.long = num.con.long, num.con.surv = num.con.surv,
  cov.long.effects = cov.long.effects, cov.surv.effects = cov.surv.effects,
  cov.long.names = cov.long.names, cov.surv.names = cov.surv.names,
  cov.long.props = cov.long.props, cov.surv.props = cov.surv.props,
  cov.long.means = cov.long.means, cov.surv.means = cov.surv.means,
  cov.long.sds = cov.long.sds, cov.surv.sds = cov.surv.sds,
  # BMA-JM inputs
  Q = Q, lin.spline.knots = lin.spline.knots, a0.surv = a0.surv, a0.long = a0.long,
  mu.Y = mu.Y, Sig.Y = Sig.Y, mu.X = mu.X, Sig.X = Sig.X, mu.alpha = mu.alpha,
  Sig.alpha = Sig.alpha, eta.iq = eta.iq, phi.iq = phi.iq, eta.tau = eta.tau,
  phi.tau = phi.tau, nu.G = nu.G, C0.inv.G = C0.inv.G, REs = REs,
  gamma0.surv = gamma0.surv, gamma0.long = gamma0.long, max.iter = max.iter,
  optim.toler.params = optim.toler.params, optim.toler.b = optim.toler.b,
  n.draws = n.draws, epsilon.star = epsilon.star, beta.star = beta.star, pi0 = pi0,
  # General simulation inputs
  nsims = nsims, print.iters = TRUE )
end.sim.time <- Sys.time()
total.sim.time <- end.sim.time - start.sim.time
total.sim.time



### Save results (rejection rate, PMPs, estimates, MSE, bias, consistency measures) as .csv files
dir.create( paste("../cluster-results/", sim.version, sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/rejection_rates", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/pmp", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/bias", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/mse", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/mean", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/sd", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/assos_par", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/loc_consis", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/prws_consis", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/prws_inconsis", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/glob_consis", sep=""), showWarnings = FALSE )
dir.create( paste("../cluster-results/", sim.version, "/comp_time", sep=""), showWarnings = FALSE )
output.file.rr <- paste("../cluster-results/", sim.version, "/rejection_rates/rr_", sim.version, "-", row.id, ".csv", sep="")
output.file.pmp <- paste("../cluster-results/", sim.version, "/pmp/pmp_", sim.version, "-", row.id, ".csv", sep="")
output.file.bias.rte <- paste("../cluster-results/", sim.version, "/bias/rte_bias_", sim.version, "-", row.id, ".csv", sep="")
output.file.mse.rte <- paste("../cluster-results/", sim.version, "/mse/rte_mse_", sim.version, "-", row.id, ".csv", sep="")
output.file.mean <- paste("../cluster-results/", sim.version, "/mean/mean_", sim.version, "-", row.id, ".csv", sep="")
output.file.sd <- paste("../cluster-results/", sim.version, "/sd/sd_", sim.version, "-", row.id, ".csv", sep="")
output.file.ap1 <- paste("../cluster-results/", sim.version, "/assos_par/assos_par1_", sim.version, "-", row.id, ".csv", sep="")
output.file.ap2 <- paste("../cluster-results/", sim.version, "/assos_par/assos_par2_", sim.version, "-", row.id, ".csv", sep="")
output.file.lc.PMDA.TE <- paste("../cluster-results/", sim.version, "/loc_consis/lc_PMDA_TE_", sim.version, "-", row.id, ".csv", sep="")
output.file.lc.PMDA.RR <- paste("../cluster-results/", sim.version, "/loc_consis/lc_PMDA_RR_", sim.version, "-", row.id, ".csv", sep="")
output.file.lc.loo <- paste("../cluster-results/", sim.version, "/loc_consis/lc_loo_", sim.version, "-", row.id, ".csv", sep="")
output.file.pc <- paste("../cluster-results/", sim.version, "/prws_consis/pc_", sim.version, "-", row.id, ".csv", sep="")
output.file.pi <- paste("../cluster-results/", sim.version, "/prws_inconsis/pi_", sim.version, "-", row.id, ".csv", sep="")
output.file.gc <- paste("../cluster-results/", sim.version, "/glob_consis/gc_", sim.version, "-", row.id, ".csv", sep="")
output.file.ct <- paste("../cluster-results/", sim.version, "/comp_time/ct_", sim.version, "-", row.id, ".csv", sep="")
write.csv(sim.results$Rejection.Rate.Surv, output.file.rr, quote = FALSE)
write.csv(sim.results$Joint.Model.PMP.Values, output.file.pmp, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Reg.Trtmt.Effect.Bias.Surv, output.file.bias.rte, quote = FALSE)
write.csv(sim.results$Reg.Trtmt.Effect.MSE.Surv, output.file.mse.rte, quote = FALSE)
write.csv(sim.results$Reg.Trtmt.Effect.Post.Means.Surv, output.file.mean, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Reg.Trtmt.Effect.Post.SDs.Surv, output.file.sd, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Rand.Int.Association.Param, output.file.ap1, quote = FALSE, row.names = FALSE)
#write.csv(sim.results$Rand.Slope.Association.Param, output.file.ap2, quote = FALSE, row.names = FALSE)
#write.csv(sim.results$Local.Consistency.PMDA.TE, output.file.lc.PMDA.TE, quote = FALSE, row.names = FALSE)
#write.csv(sim.results$Local.Consistency.PMDA.RR, output.file.lc.PMDA.RR, quote = FALSE, row.names = FALSE)
#write.csv(sim.results$Pairwise.Consistency.Surv, output.file.pc, quote = FALSE, row.names = FALSE)
#write.csv(sim.results$Pairwise.Inconsistency.Surv, output.file.pi, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Global.Consistency.Surv, output.file.gc, quote = FALSE, row.names = FALSE)
write.csv(sim.results$Computation.Time.Minutes, output.file.ct, quote = FALSE, row.names = FALSE)

