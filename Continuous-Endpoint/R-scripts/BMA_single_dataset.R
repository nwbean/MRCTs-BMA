#########################################################################################
# Execute BMA approach for one generated dataset
#
# Note: users can replace generated dataset with real data and
#       change values of elicited priors
#########################################################################################


### Load libraries and source Rcpp/R code
library(Rcpp)
sourceCpp("../R-source/bma_functions_rcpp_code.cpp")


### Input values to generate data

S <- 5                                      # number of regions
T0 <- 5                                     # max number of distinct regions to consider
N <- 1508                                   # total number of subjects
trtmt.allctn <- c(.5, .5)                   # treatment group allocation (must sum to 1)
cntrl.means <- rep(.082, S)                 # regional means for control groups
sd.Y <- 0.205                               # common st. dev. for both groups and all regions
reg.names <- c("Reg1", "Reg2", "Reg3",
               "Reg4", "Reg5")              # region names in alphabetical order
reg.allctn <- rep(.2, S)                    # region sample size allocation (must sum to 1)
trtmt.means <- c(.082, rep(.116, S-1))      # regional means for treatment groups
num.bin <- 0                                # number of binary covariates
num.con <- 0                                # number of normal covariates
cov.props <- 0                              # success probabilities for binary covariates (0 if num.bin = 0)
cov.means <- 0                              # means for normal covariates (0 if num.con = 0)
cov.sds <- 0                                # standard deviations for normal covariates (0 if num.con = 0)


### Generate data
set.seed(1)
dat <- generate_data(N, S, reg.names, reg.allctn, trtmt.allctn, trtmt.means, cntrl.means,
                     sd.Y, num.bin, num.con, cov.props, cov.means, cov.sds)
Y <- dat$Y                               # numeric vector of length N with primary outcome for each subject
region.labs <- dat$RegionLabels          # character vector of length N with region labels for each subject
trtmt.indctr <- dat$TreatmentIndicator   # integer vector of length N with 0 (control) and 1 (treatment) values
X <- dat$Covariates                      # N x p matrix of covariates (number of columns equals 0 if no covariates)


### Prior elicitation
mod.prior.alpha <- 0                            # determines type of prior assumption for models
mean.diff <- .04                                # assumed mean difference (treatment mean minus control mean)
cntrl.int <- .1                                 # assumed intercept for control group (approximately 20% increase from true value)
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
gamma0 <- 0                         # value in hypotheses for posterior prob. Pr(gamma > gamma0|D)
pi0 <- 1/S                          # value for post. prob. Pr(gamma_i/gamma > pi0|D) for ratio local consistency method
epsilon.star <- .018                # minimal clinically important difference, used in most measures of consistency
beta.star <- .5                     # probability (1 - beta_star) for which consider two regions to be clinically different,
#                                   #   used in Pr(|gamma_i - gamma_j| > epsilon.star|D) >= (1-beta.star)
#                                   #   (used for global inconsistency)



##### If using real data, replace the following NULL values in the commented-out code below
##### with the indicated values/variables

# ### Read in data
# Y <- NULL               # numeric vector of length N with primary outcome for each subject
# region.labs <- NULL     # character vector of length N with region labels for each subject
# trtmt.indctr <- NULL    # integer vector of length N with 0 (control) and 1 (treatment) values
# X <- NULL               # N x p matrix of covariates. If p = 0, set X = matrix(0, nrow = N, ncol = 0)
# S <- length( unique(region.labs) )  # number of regions
# T0 <- NULL                          # maximum number of distinct treatment effects to consider for each model (1 <= T0 <= S)
# 
# 
# ### Prior elicitation
# m0 <- NULL                          # (2*S + p) x 1 vector column with mean hyperparameters for normal prior on regression parameters
# Sig0 <- NULL                        # (2*S + p) x (2*S + p) covariance matrix hyperparameter for normal prior on regression parameters
# delta0 <- NULL                      # shape hyperparameter for gamma prior on precision parameter tau
# nu0 <- NULL                         # rate hyperparameter for gamma prior on precision parameter tau (prior mean of delta0/nu0)
# modPriors <- NULL                   # L-dimensional numeric vector of prior model probabilities, where L is number of models
# 
# 
# ## Values for BMA
# n.draws <- 100000                   # number of draws from each posterior dist.
# gamma0 <- 0                         # value in hypotheses for posterior prob. Pr(gamma > gamma0|D)
# pi0 <- 1/S                          # value for post. prob. Pr(gamma_i/gamma > pi0|D) for ratio local consistency method
# epsilon.star <- NULL                # minimal clinically important regional difference, used in most measures of consistency
# beta.star <- .5                     # probability (1 - beta_star) for which consider two regions to be clinically different,
# #                                   #   used in Pr(|gamma_i - gamma_j| > epsilon.star|D) >= (1-beta.star)
# #                                   #   (used for global inconsistency)


### Run BMA approach
start.sim.time <- Sys.time()
bma.results <- bma_single(Y, region.labs, trtmt.indctr, X, T0, modPriors, m0, Sig0, delta0,
                          nu0, gamma0, epsilon.star, beta.star, pi0, n.draws)
end.sim.time <- Sys.time()
total.sim.time <- end.sim.time - start.sim.time
total.sim.time


### Regression parameter summaries
GlobalTrtmtStats <- bma.results$GlobalTrtmtStats            # global treatment effect
#                                                           # PostProbGamma0 = Pr(gamma_G > gamma0|D)
RegionalTrtmtStats <- bma.results$RegionalTrtmtStats[,1:3]  # region-specific treatment effect
#                                                           # PostProbGamma0 = Pr(gamma_i > gamma0|D)
RegionalIntStats <- bma.results$RegionalIntStats            # region-specific intercepts
CovariateStats <- bma.results$CovariateStats                # covariate effects (NULL if p = 0)


### Consistency results

# PMPs for all models
PMPs <- bma.results$pmp

# epsilon.star-level global consistency probability
GlobConsProb <- bma.results$GlobalConsistencyProb

# Region-specific local consistency probabilities Pr(gamma_i/gamma > pi0|D) (Japanese MHLW approach)
LocConsRatio <- bma.results$RegionalTrtmtStats[,4]

# Region-specific epsilon.star-level local consistency probabilities Pr(|gamma_i - gamma_(-i)| < epsilon.star|D)
LocConsAbsDiff <- bma.results$LOOLocalConsistency

# epsilon.star-level pairwise consistency probabilities Pr(|gamma_i - gamma_j| < epsilon.star|D)
prws.cons <- bma.results$PairwiseConsistency

