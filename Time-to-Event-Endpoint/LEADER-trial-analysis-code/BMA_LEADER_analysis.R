###############################################################################################
# Execute BMA approach for the LEADER trial data
#
# ClinicalTrials.gov Identifier: NCT01179048
###############################################################################################


### Set working directory
setwd("D:/Users/mse0nxb1/Documents/Project 2 - TTE Endpoint")


### Source R scripts with additional functions
source("bma_functions_tte_endpoint.R")
source("bma-asymptotic-function.R")


########## NOTE: CODE FOR DATA PREPARATION REMOVED ##########

### Format data into necessary form
# Dataset should have the following columns:
#   eventtime: observed time (event time or censored time)
#   status: event indicator (1 = event, 0 = censored)
#   trt: treatment group indicator (1 = treatment group, 0 = control group)
#   region: numeric values for regions (1 to number of regions, sorted numerically)


### Derive values from data
S <- length(unique(dat$region))        # number of regions
T0 <- S                                # max number of distinct regions to consider
N <- nrow(dat)                         # total number of subjects


### Construct design matrix
W.mat <- matrix(0, nrow = N, ncol = S)
for(i in 1:S){
  W.mat[,i] <- ifelse( dat$trt == 1 & dat$region == unique(dat$region)[i], 1, 0 )
}


### Compile matrix of models with distinct region classifications when S = T0 = 4
mod.mat <- rbind( c(1, 1, 1, 1),
                  c(1, 1, 1, 2),
                  c(1, 1, 2, 1),
                  c(1, 1, 2, 2),
                  c(1, 1, 2, 3),
                  c(1, 2, 1, 1),
                  c(1, 2, 1, 2),
                  c(1, 2, 1, 3),
                  c(1, 2, 2, 1),
                  c(1, 2, 2, 2),
                  c(1, 2, 2, 3),
                  c(1, 2, 3, 1),
                  c(1, 2, 3, 2),
                  c(1, 2, 3, 3),
                  c(1, 2, 3, 4) )


### Prior elicitation
K <- 8                                           # number of intervals for piecewise exponential model
hr.pred <- 1.3                                   # prediction of true hazard ratio
trt.pred <- log(hr.pred)                         # prediction of true treatment effect
mu0 <- rep(trt.pred, S)                          # mean vector for prior on treatment effects
Sig0 <- diag( rep((10 * abs(trt.pred))^2, S) )   # covariance matrix for prior on treatment effects
eta0 <- matrix( 0.01, nrow = S, ncol = K )       # shape hyperparameters for prior on baseline hazards
phi0 <- matrix( 0.01, nrow = S, ncol = K )       # rate hyperparameters for prior on baseline hazards
alpha0 <- 0                                      # pre-specified value for prior model probabilities
mod.priors <- construct_modPriors(S, T0, mod.mat,  # prior model probabilities \propto exp(D.l*alpha0)
                                  alpha = alpha0)


### Create interval bounds such that approximately same number of events are in each interval,
###   and identify interval in which each patient failed or was censored
int.cuts <- c( 0, quantile( dat$eventtime[which(dat$status == 1)],
                            probs = seq(1/K, 1, length = K)[-K] ) )
dat$interval <- apply( cbind(dat$eventtime), 1, function(x){ sum(x > int.cuts) } )


### Values for BMA
gamma0 <- 0          # value for posterior prob. Pr(gamma > gamma0|D)
n.draws <- 10000     # number of posterior samples to draw when evaluating consistency
pi0 <- 1/S           # value for post. prob. Pr(gamma_i/gamma > pi0|D) - (Japanese MHLW consistency)
epsilon.star <- c( seq(.5, .7, by = .1),
                   seq(.71, .76, by = .01), 1/1.3,
                   seq(.77, .86, by = .005), 1/1.15,
                   seq(.87, 1, by = .005) )
#                    # range of possible MCIRDs, used in most measures of consistency
beta.star <- .5      # probability cutoff for which consider two regions to be clinically different,
#                    #   used in Pr(|gamma_i - gamma_j| > -log(epsilon.star)|D) >= beta.star
#                    #   (used for global inconsistency)


### Run BMA approach
set.seed(801)
start.sim.time <- Sys.time()
BMA.results <- bma.asymptotic( dat = dat, W.mat = W.mat, int.cuts = int.cuts,
                               mu0 = mu0, Sig0 = Sig0, eta0 = eta0, phi0 = phi0,
                               mod.mat = mod.mat, mod.priors = mod.priors, gamma0 = gamma0,
                               n.draws = n.draws, epsilon.star = epsilon.star,
                               beta.star = beta.star, pi0 = pi0 )
end.sim.time <- Sys.time()
total.sim.time <- end.sim.time - start.sim.time
total.sim.time


### Regression parameter summaries

# Global treatment effect summaries
gte.mean <- as.numeric(BMA.results$gte.mean)                  # posterior mean of global treatment effect
gte.sd <- as.numeric(BMA.results$gte.sd)                      # posterior sd of global treatment effect
gte.less.gamma0 <- as.numeric(BMA.results$gte.less.gamma0)    # Pr(gamma_G < gamma0|D)

# Region-specific treatment effect summaries
region.names <- as.character(unique(dat$region.name))
rte.means <- BMA.results$rte.means                # posterior means of regional trtmt effects
rownames(rte.means) <- "PosteriorMean"
colnames(rte.means) <- region.names
rte.sds <- BMA.results$rte.sds                    # posterior sd's of regional trtmt effects
rownames(rte.sds) <- "PosteriorStandDev"
colnames(rte.sds) <- region.names
rte.less.gamma0 <- BMA.results$rte.less.gamma0    # Pr(gamma_i < gamma0|D), i=1,...,S
rownames(rte.less.gamma0) <- "PosteriorProbLessThanGamma0"
colnames(rte.less.gamma0) <- region.names

# 95% central credible intervals for global and region-specific treatment effects
cred.ints <- BMA.results$credible.intervals
rownames(cred.ints) <- c( "Global", region.names )


### Consistency results

# Posterior model probabilities for all models
PMPs <- cbind(mod.mat, BMA.results$PMPs)
rownames(PMPs) <- paste("Model", 1:nrow(mod.mat), sep = "_")
colnames(PMPs) <- c( region.names, "PMP" )

# epsilon.star-level global consistency probability
global.consistency <- BMA.results$global.consistency
rownames(global.consistency) <- paste("epsilon.star", epsilon.star, sep = "=")
colnames(global.consistency) <- "GlobalConsistency"

# Region-specific local consistency probabilities (variations of Japanese PMDA approach)
loc.consis.PMDA.TE <- BMA.results$local.consistency.ratio.TE   # Pr(gamma_i/gamma_G > pi0|D)
rownames(loc.consis.PMDA.TE) <- "PMDA_LocalConsistency_TreatmentEffect"
colnames(loc.consis.PMDA.TE) <- region.names
loc.consis.PMDA.RR <- BMA.results$local.consistency.ratio.RR   # Pr(RR_i/RR_G > pi0|D)
rownames(loc.consis.PMDA.RR) <- "PMDA_LocalConsistency_RiskReduction"
colnames(loc.consis.PMDA.RR) <- region.names

# Region-specific epsilon.star-level local consistency probabilities
# Pr(|gamma_i - gamma_(-i)| < -log(epsilon.star)|D) for all values of epsilon.star
local.consistency.loo <- BMA.results$local.consistency.loo
rownames(local.consistency.loo) <- paste("epsilon.star", epsilon.star, sep = "=")
colnames(local.consistency.loo) <- region.names

# epsilon.star-level pairwise consistency and inconsistency probabilities
# Pr(|gamma_i - gamma_j| < -log(epsilon.star)|D) for all values of epsilon.star
# Each element of list corresponds to matrix of probabilities for a different epsilon.star
# Rows/columns of each epsilon.star-specific matrix correspond to regions (alphabetically)
pairwise.consistency <- BMA.results$pairwise.consistency
names(pairwise.consistency) <- paste("epsilon.star", epsilon.star, sep = "=")
pairwise.inconsistency <- BMA.results$pairwise.inconsistency
names(pairwise.inconsistency) <- paste("epsilon.star", epsilon.star, sep = "=")


### Print all results
gte.mean
gte.sd
gte.less.gamma0
rte.means
rte.sds
rte.less.gamma0
cred.ints
PMPs
global.consistency
loc.consis.PMDA.TE
loc.consis.PMDA.RR
local.consistency.loo
pairwise.consistency
pairwise.inconsistency

