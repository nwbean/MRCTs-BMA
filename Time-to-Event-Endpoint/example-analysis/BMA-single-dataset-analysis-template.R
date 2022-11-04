###############################################################################################
# Execute BMA approach for one dataset (template code)
#
# Note: users should replace all occurrences of "NULL" below with appropriate values
###############################################################################################


### Set working directory
wd.file <- "name-of-folder-with-dataset-and-code"
setwd(wd.file)


### Load libraries and source Rcpp/R code
library(Rcpp)
sourceCpp("bma_functions_tte_endpoint.cpp")
source("bma-asymptotic-function.R")


### Read in dataset
dat <- read.csv("file-path.csv", header = TRUE)
# Dataset should have the following columns:
#   eventtime: observed time (event time or censored time)
#   status: event indicator (1 = event, 0 = censored)
#   trt: treatment group indicator (1 = treatment group, 0 = control group)
#   region: numeric values for regions (1 to number of regions, sorted numerically)
# Note: additional columns with optional covariates can be included



###############################################################################################
# NOTE: REPLACE ALL INSTANCES OF "NULL" BELOW WITH APPROPRIATE VALUES
###############################################################################################


### Derive values from data
S <- length(unique(dat$region))     # number of regions
T0 <- S                             # max number of distinct regions to consider
N <- nrow(dat)                      # total number of subjects
X.mat <- NULL                       # N x p matrix of p optional covariates (excluding treatment)


### Construct design matrix
W.trt.mat <- matrix(0, nrow = N, ncol = S)
for(i in 1:S){
  W.trt.mat[,i] <- ifelse( dat$trt == 1 & dat$region == unique(dat$region)[i], 1, 0 )
}
W.mat <- cbind(W.trt.mat, X.mat)


### Prior elicitation
K <- NULL                                        # number of intervals for piecewise exponential model
hr.pred <- NULL                                  # prediction of true hazard ratio
trt.pred <- log(hr.pred)                         # prediction of true treatment effect
cov.preds <- NULL                                # predictions of true covariate effects (NULL if no covariates)
mu0 <- c( rep(trt.pred, S),
          rep(cov.preds, length(cov.preds)) )    # mean vector for prior on treatment and covariate effects
if(length(cov.preds) > 0){
  Sig0 <- diag( c( rep((10 * abs(trt.pred))^2, S),
                   (10 * abs(cov.preds))^2 ) )   # covariance matrix for prior on treatment and covariate effects
} else{
  Sig0 <- diag( rep((10 * abs(trt.pred))^2, S) ) # covariance matrix for prior on treatment effects
}
eta0 <- matrix( NULL, nrow = S, ncol = K )       # shape hyperparameters for prior on baseline hazards
phi0 <- matrix( NULL, nrow = S, ncol = K )       # rate hyperparameters for prior on baseline hazards
alpha0 <- NULL                                   # pre-specified value for prior model probabilities
mod.priors <- construct_modPriors(S, T0, alpha = alpha0)   # prior model probabilities \propto exp(T0*alpha0)


### Values for BMA
mod.mat <- compile_modMat(S, T0)    # matrix of model-specific region classifications
gamma0 <- 0                 # value for posterior prob. Pr(gamma > gamma0|D)
n.draws <- 10000            # number of posterior samples to draw when evaluating consistency
pi0 <- 1/S                  # value for post. prob. Pr(gamma_i/gamma > pi0|D) - (Japanese MHLW consistency)
epsilon.star <- NULL        # value (or range) of possible MCIRD(s), used in most measures of consistency
beta.star <- .5             # probability cutoff for which consider two regions to be clinically different,
#                           #   used in Pr(|gamma_i - gamma_j| > -log(epsilon.star)|D) >= beta.star
#                           #   (used for global inconsistency)


### Create interval bounds such that approximately same number of events are in each interval
### Identify interval in which each patient failed or was censored
int.cuts <- c( 0, quantile( dat$eventtime[which(dat$status == 1)],
                            probs = seq(1/K, 1, length = K)[-K] ) )
dat$interval <- apply( cbind(dat$eventtime), 1, function(x){ sum(x > int.cuts) } )


### Run BMA approach
set.seed(NULL)
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
rte.means <- BMA.results$rte.means                # posterior means of regional trtmt effects
rownames(rte.means) <- "PosteriorMean"
colnames(rte.means) <- unique(dat$region)
rte.sds <- BMA.results$rte.sds                    # posterior sd's of regional trtmt effects
rownames(rte.sds) <- "PosteriorStandDev"
colnames(rte.sds) <- unique(dat$region)
rte.less.gamma0 <- BMA.results$rte.less.gamma0    # Pr(gamma_i < gamma0|D), i=1,...,S
rownames(rte.less.gamma0) <- "PosteriorProbLessThanGamma0"
colnames(rte.less.gamma0) <- unique(dat$region)

# 95% central credible intervals for global and region-specific treatment effects
cred.ints <- BMA.results$credible.intervals
rownames(cred.ints) <- c( "Global", unique(dat$region) )

# Covariate effect summaries (NULL if p = 0)
cov.means <- BMA.results$cov.means                # posterior mean of covariate effects
rownames(cov.means) <- "PosteriorMean"
colnames(cov.means) <- colnames(X.mat)
cov.sds <- BMA.results$cov.sds                    # posterior sd's of covariate effects
rownames(cov.sds) <- "PosteriorStandDev"
colnames(cov.sds) <- colnames(X.mat)



### Consistency results

# Posterior model probabilities for all models
PMPs <- cbind(mod.mat, BMA.results$PMPs)
rownames(PMPs) <- paste("Model", 1:nrow(mod.mat), sep = "_")
colnames(PMPs) <- c( unique(dat$region), "PMP" )

# epsilon.star-level global consistency probability
global.consistency <- BMA.results$global.consistency
rownames(global.consistency) <- paste("epsilon.star", epsilon.star, sep = "=")
colnames(global.consistency) <- "GlobalConsistency"

# Region-specific local consistency probabilities (variations of Japanese PMDA approach)
loc.consis.PMDA.TE <- BMA.results$local.consistency.ratio.TE   # Pr(gamma_i/gamma_G > pi0|D)
rownames(loc.consis.PMDA.TE) <- "PMDA_LocalConsistency_TreatmentEffect"
colnames(loc.consis.PMDA.TE) <- unique(dat$region)
loc.consis.PMDA.RR <- BMA.results$local.consistency.ratio.RR   # Pr(RR_i/RR_G > pi0|D)
rownames(loc.consis.PMDA.RR) <- "PMDA_LocalConsistency_RiskReduction"
colnames(loc.consis.PMDA.RR) <- unique(dat$region)

# Region-specific epsilon.star-level local consistency probabilities
# Pr(|gamma_i - gamma_(-i)| < -log(epsilon.star)|D) for all values of epsilon.star
local.consistency.loo <- BMA.results$local.consistency.loo
rownames(local.consistency.loo) <- paste("epsilon.star", epsilon.star, sep = "=")
colnames(local.consistency.loo) <- unique(dat$region)

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
cov.means
cov.sds
PMPs
global.consistency
loc.consis.PMDA.TE
loc.consis.PMDA.RR
local.consistency.loo
pairwise.consistency
pairwise.inconsistency

