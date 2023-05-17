###############################################################################################
# Execute BMA joint modeling approach for one dataset (template code)
#
# Note: users should replace all occurrences of "NULL" below with appropriate values
###############################################################################################


### Set working directory
wd.file <- "name-of-folder-with-dataset-and-code"
setwd(wd.file)


### Load libraries and source Rcpp/R code
library(Rcpp)
sourceCpp("bma-functions-joint-model.cpp")                 # Rcpp code with extra functions
source("../R-scripts/data-simulation-functions.R")         # functions to simulate all data
source("../R-scripts/bma-joint-model-function.R")          # function for BMA (joint model)


### Read in datasets for the survival and longitudinal submodels

# Survival dataset
surv.dat <- read.csv("file-path-of-survival-data.csv", header = TRUE)
# Dataset should have the following columns:
#   subjid: subject ID
#   eventtime: observed time (event time or censored time)
#   status: event indicator (1 = event, 0 = censored)
#   trt: treatment group indicator (1 = treatment group, 0 = control group)
#   region: numeric values for regions (1 to number of regions, sorted numerically)
# Note: additional columns with optional covariates can be included

# Longitudinal dataset
long.dat <- read.csv("file-path-of-longitudinal-data.csv", header = TRUE)
# Dataset should have the following columns:
#   subjid: subject ID
#   X.ijk: observed value of longitudinal marker at time k for subject j from region i
#   time: time when X.ijk was observed (e.g., visit time)
#   trt: treatment group indicator (1 = treatment group, 0 = control group)
#   region: numeric values for regions (1 to number of regions, sorted numerically)
# Note: additional columns with optional covariates can be included



###############################################################################################
# NOTE: REPLACE ALL INSTANCES OF "NULL" BELOW WITH APPROPRIATE VALUES
###############################################################################################


### Compile matrix of model region classifications (for both survival and longitudinal submodels)
S <- length(unique(surv.dat$region))  # number of regions
T0 <- NULL                            # max number of distinct regions to consider
mod.mat <- compile_modMat(S, T0)      # region classifications


### Names of additional covariates (if any) included in the survival and longitudinal submodels
cov.surv.names <- NULL
cov.long.names <- NULL


### Create interval bounds for the baseline hazard such that approximately same number of events
### are in each interval
Q <- NULL      # number of time intervals with separate baseline hazards
int.cuts <- c( 0, quantile( surv.dat$eventtime[which(surv.dat$status == 1)],
                            probs = seq(1/Q, 1, length = Q)[-Q] ) )


### Location of knots for linear splines in the longitudinal submodel
lin.spline.knots <- NULL         # locations of knots for linear splines
REs <- "int"                     # type of random effects to include ("int" for random intercepts only)


### Prior elicitation for parameters in survival submodel
hr.pred.surv <- NULL                     # prediction of true hazard ratio (surv.)
trt.pred.surv <- log(hr.pred.surv)       # prediction of true treatment effect (surv.)
cov.pred.surv <- NULL                    # predictions of true covariate effects (surv.)
mu.Y <- rep(trt.pred.surv, S)            # mean vector for prior on trt and covariate effects (surv.)
Sig.Y <- diag( (10 * abs(mu.Y))^2 )      # covariance matrix for prior on trt and covariate effects (surv.)
eta.iq <- matrix( NULL, nrow = S, ncol = Q )  # shape hyperparameters for prior on baseline hazards (surv.)
phi.iq <- matrix( NULL, nrow = S, ncol = Q )  # rate hyperparameters for prior on baseline hazards (surv.)
mu.alpha <- as.matrix(NULL)              # prior mean for association parameter
Sig.alpha <- as.matrix(NULL)             # prior variance for association parameter


### Prior elicitation for parameters in longitudinal submodel
prior.int.mean.long <- NULL              # prior mean for all regression effects (long.)
prior.reg.mean.long <- NULL              # prior mean for all regression effects (long.)
num.effects.long <- 2*S + 1 + length(lin.spline.knots) +
  S*length(lin.spline.knots)             # number of regression effects minus S intercepts (long.)
mu.X <- c( rep(prior.int.mean.long, S),  # mean vector for prior on regression effects (long.)
           rep(prior.reg.mean.long, num.effects.long) )
Sig.X <- diag( (10 * abs(mu.X))^2 )      # cov. matrix for prior on regression effects (long.)
eta.tau <- NULL            # shape for gamma prior on precision of likelihood error (long.)
phi.tau <- NULL            # rate for gamma prior on precision of likelihood error (long.)
nu.G <- NULL               # degrees of freedom for Wishart prior on G^-1 (precision matrix of REs) (long.)
C0.G <- as.matrix(NULL)    # scale matrix in Wishart prior on G^-1
C0.inv.G <- solve(C0.G)    # inverse of scale matrix in Wishart prior on G^-1


### Prior elicitation for models in model space
a0.surv <- NULL                  # pre-specified value for prior model probabilities (surv.)
a0.long <- NULL                  # pre-specified value for prior model probabilities (long.)
# prior model probabilities \propto exp(D.surv.l*a0.surv + D.long.l*a0.long),
#    where D.surv.l and D.long.l are the number of distinct trt effects for the l^th
#    survival and longitudinal submodels, respectively.
# NOTE: a vector of custom prior model probabilities can be read into the BMA function
#    using the function argument "mod.priors" In this case, the function arguments
#    "a0.surv" and "a0.long" should not be used.
mod.priors <- NULL


### Values for BMA
max.iter <- 10       # maximum number of iterations in optimization algorithm for each model
optim.toler.params <- .005   # tolerance for parameter convergence in optimization algorithm
optim.toler.b <- .05         # tolerance for random effects convergence in optimization algorithm
gamma0.surv <- 0     # value for posterior prob. Pr(gamma_G > gamma0|D) (survival model)
gamma0.long <- 0     # value for posterior prob. Pr(gamma_G > gamma0|D) (longitudinal model)
n.draws <- 10000     # number of posterior samples to draw when evaluating consistency
pi0 <- 1/S           # value for post. prob. Pr(gamma_i/gamma > piO|D) - (Japanese MHLW consistency)
epsilon.star <- NULL
#                    # range of possible MCIRDS, used in most measures of consistency
beta.star <- .5      # probability cutoff for which consider two regions to be clinically different,
#                    # used in Pr(|gamma_i - gamma_j| > -log(epsilon.star)|D) >= beta.star
#                    # (used for global inconsistency)


### Create interval bounds such that approximately same number of events are in each interval
int.cuts <- c( 0, quantile( surv.dat$eventtime[which(surv.dat$status == 1)],
                            probs = seq(1/Q, 1, length = Q)[-Q] ) )


### Run BMA approach
set.seed(NULL)
start.time <- Sys.time()
BMA.results <- bma.joint.model( surv.dat = surv.dat, cov.surv.names = cov.surv.names, int.cuts = int.cuts,
                                long.dat = long.dat, lin.spline.knots = lin.spline.knots,
                                cov.long.names = cov.long.names, mod.mat = mod.mat, a0.surv = a0.surv,
                                a0.long = a0.long, mod.priors = mod.priors, mu.Y = mu.Y, Sig.Y = Sig.Y,
                                mu.X = mu.X, Sig.X = Sig.X, mu.alpha = mu.alpha, Sig.alpha = Sig.alpha,
                                eta.iq = eta.iq, phi.iq = phi.iq, eta.tau = eta.tau, phi.tau = phi.tau,
                                nu.G = nu.G, C0.inv.G = C0.inv.G, REs = REs, max.iter = max.iter,
                                optim.toler.params = optim.toler.params, optim.toler.b = optim.toler.b,
                                gamma0.surv = gamma0.surv, gamma0.long = gamma0.long, n.draws = n.draws,
                                epsilon.star = epsilon.star, beta.star = beta.star, pi0 = pi0 )
end.time <- Sys.time()
total.time <- end.time - start.time
total.time



### Regression parameter summaries - survival submodel

# Global treatment effect summaries
gte.surv.mean <- as.numeric(BMA.results$gte.surv.mean)        # posterior mean of global treatment effect
gte.surv.sd <- as.numeric(BMA.results$gte.surv.sd)            # posterior sd of global treatment effect
gte.surv.less.gamma0 <- as.numeric(BMA.results$gte.surv.less.gamma0)    # Pr(gamma_{Y,G} < gamma0|D)

# Region-specific treatment effect summaries
rte.surv.means <- BMA.results$rte.surv.means      # posterior means of region trtmt effects
rownames(rte.surv.means) <- "PosteriorMean"
colnames(rte.surv.means) <- unique(surv.dat$region)
rte.surv.sds <- BMA.results$rte.surv.sds          # posterior sd's of region trtmt effects
rownames(rte.surv.sds) <- "PosteriorStandDev"
colnames(rte.surv.sds) <- unique(surv.dat$region)
rte.surv.less.gamma0 <- BMA.results$rte.surv.less.gamma0    # Pr(gamma_{Y,i} < gamma0|D), i=1,...,S
rownames(rte.surv.less.gamma0) <- "PosteriorProbLessThanGamma0"
colnames(rte.surv.less.gamma0) <- unique(surv.dat$region)

# 95% central credible intervals for global and region-specific treatment effects
cred.ints.surv <- BMA.results$credible.intervals.surv
rownames(cred.ints.surv) <- c( "Global", unique(surv.dat$region) )

# Covariate effect summaries (NULL if p = 0)
cov.surv.means <- BMA.results$cov.surv.means      # posterior mean of covariate effects
rownames(cov.surv.means) <- "PosteriorMean"
cov.surv.sds <- BMA.results$cov.surv.sds          # posterior sd's of covariate effects
rownames(cov.surv.sds) <- "PosteriorStandDev"



### Regression parameter summaries - longitudinal submodel

# Global treatment effect summaries
gte.long <- rbind(BMA.results$gte.long.mean,         # posterior mean of global treatment effect
                  BMA.results$gte.long.sd,           # posterior sd of global treatment effect
                  BMA.results$gte.long.probs[1,],    # Pr(gamma_{X,G} < gamma0|D)
                  BMA.results$gte.long.probs[2,])    # Pr(gamma_{X,G} > gamma0|D)
rownames(gte.long) <- c("PosteriorMean", "PosteriorSD", "ProbLessThanGamma0", "ProbGrtrThanGamma0")

# Region-specific treatment effect summaries
rte.long.means <- BMA.results$rte.long.means      # posterior means of region trtmt effects
rte.long.sds <- BMA.results$rte.long.sds          # posterior sd's of region trtmt effects
rte.long.less.gamma0 <- BMA.results$rte.long.less.gamma0    # Pr(gamma_{X,i} < gamma0|D), i=1,...,S
rte.long.grtr.gamma0 <- BMA.results$rte.long.grtr.gamma0    # Pr(gamma_{X,i} > gamma0|D), i=1,...,S

# Covariate effect summaries (NULL if p = 0)
cov.long.means <- BMA.results$cov.long.means      # posterior mean of covariate effects
rownames(cov.long.means) <- "PosteriorMean"
cov.long.sds <- BMA.results$cov.long.sds          # posterior sd's of covariate effects
rownames(cov.long.sds) <- "PosteriorStandDev"



### Association parameter summary
alpha.mean <- BMA.results$alpha.mean
alpha.sd <- BMA.results$alpha.sd
cred.int.alpha <- BMA.results$credible.interval.alpha



### Posterior model probabilities

# Posterior model probabilities for all models - first number in model name is the longitudinal
# submodel, the second number is the survival submodel
PMPs <- data.frame(
  long.mod = rep(1:nrow(mod.mat), each = nrow(mod.mat)),
  surv.mod = rep(1:nrow(mod.mat), nrow(mod.mat)),
  long.reg = mod.mat[rep(1:nrow(mod.mat), each = nrow(mod.mat)),],
  surv.reg = mod.mat[rep(1:nrow(mod.mat), nrow(mod.mat)),],
  PMP = BMA.results$PMPs)
rownames(PMPs) <- paste("Model", rep(1:nrow(mod.mat), each = nrow(mod.mat)),
                        rep(1:nrow(mod.mat), nrow(mod.mat)), sep = "_")
mar.PMPs.surv <- data.frame(
  surv.mod = 1:nrow(mod.mat),
  surv.reg = mod.mat,
  PMP = tapply(PMPs$PMP, PMPs$surv.mod, sum))    # marginal PMPs for survival submodels
mar.PMPs.long <- data.frame(
  long.mod = 1:nrow(mod.mat),
  long.reg = mod.mat,
  tapply(PMPs$PMP, PMPs$long.mod, sum))    # marginal PMPs for longitudianl submodels



### Consistency results (survival treatment effects)

# epsilon.star-level global consistency probability
global.consistency.surv <- BMA.results$global.consistency.surv
rownames(global.consistency.surv) <- paste("epsilon.star", epsilon.star, sep = "=")
colnames(global.consistency.surv) <- "GlobalConsistency"

# Region-specific local consistency probabilities (variations of Japanese PMDA approach)
loc.consis.PMDA.TE <- BMA.results$local.consistency.ratio.TE   # Pr(gamma_{Y,i}/gamma_{Y,G} > pi0|D)
rownames(loc.consis.PMDA.TE) <- "PMDA_LocalConsistency_TreatmentEffect"
colnames(loc.consis.PMDA.TE) <- unique(surv.dat$region)
loc.consis.PMDA.RR <- BMA.results$local.consistency.ratio.RR   # Pr(RR_i/RR_G > pi0|D)
rownames(loc.consis.PMDA.RR) <- "PMDA_LocalConsistency_RiskReduction"
colnames(loc.consis.PMDA.RR) <- unique(surv.dat$region)

# epsilon.star-level pairwise consistency and inconsistency probabilities
# Pr(|gamma_{Y,i} - gamma_{Y,j}| < -log(epsilon.star)|D) for all values of epsilon.star
# Each element of list corresponds to matrix of probabilities for a different epsilon.star
# Rows/columns of each epsilon.star-specific matrix correspond to regions (alphabetically)
pairwise.consistency.surv <- BMA.results$pairwise.consistency.surv
names(pairwise.consistency.surv) <- paste("epsilon.star", epsilon.star, sep = "=")
pairwise.inconsistency.surv <- BMA.results$pairwise.inconsistency.surv
names(pairwise.inconsistency.surv) <- paste("epsilon.star", epsilon.star, sep = "=")


### Print all results
gte.surv.mean
gte.surv.sd
gte.surv.less.gamma0
rte.surv.means
rte.surv.sds
rte.surv.less.gamma0
cred.ints.surv
cov.surv.means
cov.surv.sds
gte.long
rte.long.means
rte.long.sds
rte.long.less.gamma0
rte.long.grtr.gamma0
alpha.mean
alpha.sd
cred.int.alpha
PMPs
mar.PMPs.surv
mar.PMPs.long
global.consistency.surv
loc.consis.PMDA.TE
loc.consis.PMDA.RR
pairwise.consistency.surv
pairwise.inconsistency.surv

