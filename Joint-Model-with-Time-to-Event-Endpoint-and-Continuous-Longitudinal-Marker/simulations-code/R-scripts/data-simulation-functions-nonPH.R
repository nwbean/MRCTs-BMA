############################################################################################
# FUNCTIONS TO SIMULATE SURVIVAL AND LONGITUDINAL DATA
############################################################################################


### Survival function conditional on surviving beyond time point t_{q-1} - calculate probability
### that individual i survives past time t
#      t.val: time point of interest (t)
#      t.qm1: time point t_{q-1} (i.e., beginning of qth interval)
#      lambda.q: baseline hazard for qth interval
#      e.lp: exp(lp) for a given linear predictor lp (portion that does not depend on time)
#      tv.effect: time-varying fixed effect
S.int.q.fun <- function(t.val, t.qm1, lambda.q, e.lp, tv.effect, trt){
  
  if(trt == 0){
    surv.prob <- exp( -lambda.q * e.lp * tv.effect^-1 *
                        (exp(t.val * tv.effect) - exp(t.qm1 * tv.effect)) )
  } else if(trt == 1){
    surv.prob <- exp( -lambda.q * e.lp * (t.val - t.qm1) )
  }
  return(surv.prob)
  
}


### Inverse of survival function conditional on surviving beyond time point t_{q-1} - calculate
### time t corresponding to a given survival probability
#      u: given survival probability
#      t.qm1: time point t_{q-1} (i.e., beginning of qth interval)
#      lambda.q: baseline hazard for qth interval
#      e.lp: exp(lp) for a given linear predictor lp (portion that does not depend on time)
#      tv.effect: time-varying fixed effect
S.int.q.inv.fun <- function(u, t.qm1, lambda.q, e.lp, tv.effect, trt){
  
  if(trt == 0){
    t.val <- tv.effect^-1 * log( -log(u) * tv.effect /
                                   (lambda.q * e.lp) + exp(t.qm1 * tv.effect) )
  } else if(trt == 1){
    t.val <- -log(u) / (lambda.q * e.lp) + t.qm1
  }
  return(t.val)
  
}

### Function to sample N event times from piecewise exponential model with longitudinal trajectory
#      N: number of subjects to sample
#      lambda.vec: constant baseline hazards for each interval
#      t.knots: knots where new time intervals begin
#      e.lp.vec: N x 1 vector of exp(lp) values for linear predictors lp (portions that do not
#                depend on time)
#      tv.effect: time-varying fixed effect
sim.events <- function(N, lambda.vec, t.knots, e.lp.vec, tv.effect, trt.vec){
  
  # Empty vectors/matrices to store results
  Q <- length(t.knots) - 1          # number of intervals
  time.vec <- rep(NA, N)            # store event times
  
  # Loop through subjects
  for(i in 1:N){
    
    # Loop through time intervals
    for(q in 1:Q){
      
      # Calculate probability that individual i survives past time t_{q+1}
      # given they survived past t_{q}
      if(trt.vec[i] == 0){
        pi.iq <- exp( -lambda.vec[q] * e.lp.vec[i] * tv.effect^-1 *
                        (exp(t.knots[q+1] * tv.effect) - exp(t.knots[q] * tv.effect)) )
      } else if(trt.vec[i] == 1){
        pi.iq <- exp( -lambda.vec[q] * e.lp.vec[i] * (t.knots[q+1] - t.knots[q]) )
      }
      if(t.knots[q+1] == Inf)  pi.iq <- 0
      
      # Draw indicator for whether individual i survived past time t_{q+1} (i.e., no event in interval q)
      # (1: no event in interval q; 0: event in interval q)
      surv.q <- rbinom(1, 1, pi.iq)       # random Bernoulli draw with probability pi.iq
      
      # Sample event time if individual i has an event in interval q
      if(surv.q == 0){
        
        # Sample random uniform draw between S(t_q|t > t_{q-1}) and S(t_{q-1}|t > t_{q-1})
        # and identify corresponding time value
        S.t.q <- S.int.q.fun(t.knots[q+1], t.knots[q], lambda.vec[q], e.lp.vec[i], tv.effect,
                             trt.vec[i])
        S.t.qm1 <- S.int.q.fun(t.knots[q], t.knots[q], lambda.vec[q], e.lp.vec[i], tv.effect,
                               trt.vec[i])
        u <- runif(1, min = S.t.q, max = S.t.qm1)
        time.vec[i] <- S.int.q.inv.fun(u, t.knots[q], lambda.vec[q], e.lp.vec[i], tv.effect,
                                       trt.vec[i])
        
        # Break the loop that iterates through time intervals
        break
        
      }
      
    }
    
  }
  
  return(time.vec)
}



### Function to compile survival and longitudinal datasets
#      N: number of subjects to sample
#      reg.vec: vector of region indicators (1 through S, where S is number of regions)
#      trt.vec: vector of treatment indicators (1: treatment; 0: placebo)
#      base.haz.vec: Q x 1 vector of constant baseline hazards for each interval
#      pw.times: Q x 1 vector of knots where new time intervals begin
#      W.surv: N x (S + p.surv) matrix of region-by-trt indicators and covariates (survival)
#      W.long: N x (2*S + p.long) matrix of region and reg-by-trt indicators and covariates (longitudinal)
#      coefs.surv.vec: (S + p.surv) x 1 vector of regression coefficients (survival)
#      sigma.e: likelihood standard deviation of longitudinal submodel (i.e., measurement error)
#      b.mat: N x 1 (or N x 2) matrix of individual-specific random intercepts (and slopes)
#      alpha.vec: association parameter(s) (multiplied by random effects in survival submodel)
#      vst.times: K x 1 vector of visit times when longitudinal measures are observed
#      vst.means.0: K x 1 vector of means of longitudinal measures for each visit (control group)
#      vst.means.1: K x 1 vector of means of longitudinal measures for each visit (trtmt group)
#      dropout.rate: dropout rate
#      max.time: maximum follow up time
sim.joint.data <- function( N, reg.vec, trt.vec, base.haz.vec, pw.times, W.surv, W.long,
                            coefs.surv.vec, sigma.e, b.mat, alpha.vec, vst.times, vst.means.0,
                            vst.means.1, dropout.rate, max.time ){
  
  ## Simulate survival data
  
  # N x 1 vector of exp(lp) where lp is the linear predictor (portion that does not depend on time)
  e.lp.vec <- exp( b.mat %*% cbind(alpha.vec) + as.matrix(W.surv) %*% cbind(coefs.surv.vec) )
  tv.effect <- abs(log(.868)) * 2 / max.time     # fixed effect for time varying effect
  
  # Simulate time from enrollment to event (months)
  event.time <- sim.events(N = N, lambda.vec = base.haz.vec, t.knots = c(pw.times, Inf),
                           e.lp.vec = e.lp.vec, tv.effect = tv.effect, trt.vec = trt.vec)
  
  # Simulate time from enrollment to stochastic censorship time (theoretical dropout time)
  censored.time <- rexp(N, dropout.rate)
  
  # Determine observed time and event indicator
  obs.time <- ifelse( event.time < censored.time, event.time, censored.time )
  obs.time <- ifelse( obs.time > max.time, max.time, obs.time )
  event.status <- ifelse( event.time < censored.time & event.time < max.time, 1, 0 )
  
  # Compile dataset
  if(ncol(W.surv) > S){       # check if covariates are included for survival model
    surv.dat <- data.frame( subjid = 1:N,
                            eventtime = obs.time,
                            status = event.status,
                            trt = trt.vec,
                            region = reg.vec,
                            W.surv[,(S+1):ncol(W.surv)] )
  } else{
    surv.dat <- data.frame( subjid = 1:N,
                            eventtime = obs.time,
                            status = event.status,
                            trt = trt.vec,
                            region = reg.vec )
  }
  
  
  ## Simulate longitudinal data
  
  # Number of visits per subject
  n.vsts <- sapply(surv.dat$subjid,                   # number of visits per subject
                   function(x){ sum( vst.times <= surv.dat$eventtime[x] ) } )
  vst.vec <- unlist( sapply(n.vsts, function(x){ vst.times[1:x] }) )   # vector of visits
  
  # Convert longitudinal data to long format
  e.ijk.vec <- rnorm(sum(n.vsts), 0, sigma.e)         # sample random errors
  if(ncol(b.mat) == 1){
    b.mat.long <- as.matrix(rep(b.mat, n.vsts))       # long format for random effects
  } else if(ncol(b.mat) == 2){
    b.mat.long <- cbind(rep(b.mat[,1], n.vsts),       # long format for random effects
                        rep(b.mat[,2], n.vsts))
  }
  trt.vec.long <- rep(trt.vec, n.vsts)                # long format for treatment vector
  reg.vec.long <- rep(reg.vec, n.vsts)                # long format for region indicator vector
  
  # Sample longitudinal responses for baseline visit
  X.ijk <- rep(0, sum(n.vsts))
  which.obs.bsln.0 <- which(vst.vec == 0 & trt.vec.long == 0)
  which.obs.bsln.1 <- which(vst.vec == 0 & trt.vec.long == 1)
  X.ijk[which.obs.bsln.0] <- b.mat.long[which.obs.bsln.0,1] + vst.means.0[1] +
    e.ijk.vec[which.obs.bsln.0]        # baseline values for control group
  X.ijk[which.obs.bsln.1] <- b.mat.long[which.obs.bsln.1,1] + vst.means.1[1] +
    e.ijk.vec[which.obs.bsln.1]        # baseline values for treatment group
  
  # Sample longitudinal responses for all other visits
  if(ncol(b.mat) == 1){     # random intercepts only
    
    for(k in 2:length(vst.times)){
      # Observations in control group for visit k
      which.obs.0 <- which(vst.vec == vst.times[k] & trt.vec.long == 0)
      # Observations in treatment group for visit k
      which.obs.1 <- which(vst.vec == vst.times[k] & trt.vec.long == 1)
      # Sample observed longitudinal responses for control group for visit k
      X.ijk[which.obs.0] <- b.mat.long[which.obs.0,1] + vst.means.0[k] + e.ijk.vec[which.obs.0]
      # Sample observed longitudinal responses for treatment group for visit k
      X.ijk[which.obs.1] <- b.mat.long[which.obs.1,1] + vst.means.1[k] + e.ijk.vec[which.obs.1]
    }
    
  } else{     # random intercepts and slopes
    
    for(k in 2:length(vst.times)){
      # Observations in control group for visit k
      which.obs.0 <- which(vst.vec == vst.times[k] & trt.vec.long == 0)
      # Observations in treatment group for visit k
      which.obs.1 <- which(vst.vec == vst.times[k] & trt.vec.long == 1)
      # Sample observed longitudinal responses for control group for visit k
      X.ijk[which.obs.0] <- b.mat.long[which.obs.0,1] + b.mat.long[which.obs.0,2] * vst.times[k] +
        vst.means.0[k] + e.ijk.vec[which.obs.0]
      # Sample observed longitudinal responses for treatment group for visit k
      X.ijk[which.obs.1] <- b.mat.long[which.obs.1,1] + b.mat.long[which.obs.1,2] * vst.times[k] +
        vst.means.1[k] + e.ijk.vec[which.obs.1]
    }
    
  }
  
  # Create data frame to store longitudinal data in long format
  if(ncol(W.long) > 2*S){       # check if covariates are included for longitudinal model
    covs.long.long <- W.long[rep(1:N, n.vsts), (2*S + 1):ncol(W.long)]  # long format for covariates
    long.dat <- data.frame( subjid = rep(1:N, n.vsts),
                            X.ijk = X.ijk,
                            time = vst.vec,
                            trt = trt.vec.long,
                            region = reg.vec.long,
                            covs.long.long )
  } else{
    long.dat <- data.frame( subjid = rep(1:N, n.vsts),
                            X.ijk = X.ijk,
                            time = vst.vec,
                            trt = trt.vec.long,
                            region = reg.vec.long )
  }
  
  ## Return both datasets saved as a list
  data.list <- list(surv.dat = surv.dat, long.dat = long.dat)
  return(data.list)
  
}



### Plot longitudinal trajectories (X.ijk)
# obs.means.0 <- tapply(long.dat$X.ijk[which(long.dat$trt == 0)],
#                       long.dat$time[which(long.dat$trt == 0)], mean)
# obs.means.1 <- tapply(long.dat$X.ijk[which(long.dat$trt == 1)],
#                       long.dat$time[which(long.dat$trt == 1)], mean)
# plot(vst.times, obs.means.0, type = "l", col = "blue",
#      ylim = c( min(c(obs.means.0, obs.means.1)) - .2, max(c(obs.means.0, obs.means.1)) + .2 ))
# lines(vst.times, obs.means.1, col = "red")



### Kaplan-Meier plot from simulated dataset (red for x_i = 0, blue for x_i = 1)
# library(survival)
# plot( survfit(Surv(eventtime, status) ~ trt, data = surv.dat), xlab = "Time",
#       ylab = "Overall survival probability", col = c("red", "blue") )


## Survival function - calculate survival probability at time t
## (used to plot survival curve)
##    t: time point of interest
##    lambda.vec: constant baseline hazards for each interval
##    t.knots: boundaries of time intervals (include 0 and Inf)
##    e.lp: exp(lp) for a given linear predictor lp minus alpha * (alpha_X,2 + b_ij2)
# S.fun <- function(t.val, lambda.vec, t.knots, e.lp){
# 
#   # Identify which interval contains t
#   int.q <- max( which(t.val > t.knots) )
# 
#   # Calculate survival probability
#   log.surv.prob <- 0
#   for(q in 1:int.q){
# 
#     if(q == int.q){
#       log.surv.prob <- log.surv.prob - lambda.vec[q] * e.lp * (t.val - t.knots[q])
#     }
#     if(q < int.q){
#       log.surv.prob <- log.surv.prob - lambda.vec[q] * e.lp * (t.knots[q+1] - t.knots[q])
#     }
# 
#   }
# 
#   return(exp(log.surv.prob))
# 
# }


## Survival curves for x_i = 0 and x_i = 1
# t.vals <- seq(.1, 60, by = 0.1)
# e.lp.0 <- mean( e.lp.vec[which(trt.vec == 0)] )
# e.lp.1 <- mean( e.lp.vec[which(trt.vec == 1)] )
# surv.probs.0 <- numeric(length(t.vals))
# surv.probs.1 <- numeric(length(t.vals))
# for(i in 1:length(t.vals)){
#   surv.probs.0[i] <- S.fun(t.vals[i], lambda.vec = base.haz.vec, t.knots = pw.times, e.lp.0)
#   surv.probs.1[i] <- S.fun(t.vals[i], lambda.vec = base.haz.vec, t.knots = pw.times, e.lp.1)
# }
# lines(t.vals, surv.probs.0, col = "red", lty = 2)   # survival curve for subjects with x_i = 0
# lines(t.vals, surv.probs.1, col = "blue", lty = 2)  # survival curve for subjects with x_i = 1
# legend( "bottom", c("K-M curve when trt = 0", "K-M curve when trt = 1",
#                     "survival probability when trt = 0", "survival probability when trt = 1"),
#         col = c("red", "blue", "red", "blue"), lty = c(1,1,2,2), bty = "n" )


