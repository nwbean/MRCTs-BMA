#####################################################################################################
# R Functions Used in BMA Algorithm
#####################################################################################################


### Function that determines the number of models for specified S, T0
num_models <- function(S, T0){
  
  # Calculate number of models
  nModels <- 0
  for(p in 1:T0){
    
    t = 0
    for(j in 0:p){
      t = t + (-1)^(p-j) * j^S * choose(p,j)
    }
    
    nModels <- t / factorial(p) + nModels
    
  }
  return(nModels)
  
}



### Function to construct model priors
# User input must include the following:
#     S: number of regions
#     T0: maximum number of regions to consider in a given models
#     mod_mat: matrix of models with distinct region classifications
#     alpha: tuning parameter, where p(M_l) is proportional to exp{T_l^alpha} for the l^th model
#            and T_l is the number of distinct treatment effects for model M_l
#            If alpha = 0 (default), set each model prior proportional to 1 (i.e., uniform);
#            If alpha > 0, priors favor models with heterogenious treatment effects;
#            If alpha < 0, priors favor models with homorogenious treatment effects;
#     See Biostatistics "Bayesian adaptive basket trial design using model averaging" (Psioda et al., 2019)
construct_modPriors <- function(S, T0, mod_mat, alpha = 0){
  
  # Check function inputs
  if(S <= 0)  Rcpp::stop("S must be greater than 0");
  if(T0 <= 0)  Rcpp::stop("T0 must be greater than 0");
  if(S < T0)  Rcpp::stop("T0 must be less than or equal to S");
  
  # Determine number of models
  numMods <- num_models(S, T0)
  
  # Initialize model priors
  mod_priors <- numeric(numMods)
  num_distinct <- 0
  for(l in 1:numMods){
    num_distinct <- length(unique(mod_mat[l,]))
    mod_priors[l] <- exp(num_distinct * alpha)
  }
  mod_priors <- mod_priors / sum(mod_priors)
  
  return(mod_priors)
  
}



### Function to construct mu_l (updated location vector for regression priors)
update_mu_l <- function(mu0, modMat_l){
  
  # Extract necessary values
  S <- length(modMat_l)                      # number of regions
  unq_parms_l <- sort(unique(modMat_l))
  T_l <- length(unq_parms_l)                 # number of distinct region labels
  p_cov <- length(mu0) - S                   # number of covariates (if any)

# Construct mu_l - assume mean hyperparameters for region-specific trtmt effects are same across regions
mu_l <- numeric(T_l + p_cov)
mu_l_temp1 <- mu0[1:T_l]

if(p_cov > 0){
  
  mu_l_temp2 <- mu0[(S+1):length(mu0)]
  
  for(i in 1:length(mu_l_temp1)){
    mu_l[i] <- mu_l_temp1[i]
  }
  for(i in 1:length(mu_l_temp2)){
    mu_l[ length(mu_l_temp1) + i ] <- mu_l_temp2[i]
  }
  
} else{
  
  for(i in 1:length(mu_l_temp1)){
    mu_l[i] <- mu_l_temp1[i]
  }
  
}

return(mu_l)

}



### Function to construct Sig_l (updated dispersion matrix for regression priors)
update_Sig_l <- function(Sig0, modMat_l){
  
  # Extract necessary values
  S <- length(modMat_l)                    # number of regions
  unq_parms_l <- sort(unique(modMat_l))
  T_l <- length(unq_parms_l)               # number of distinct region labels
  p_cov <- ncol(Sig0) - S                  # number of covariates (if any)
  
  # Construct Sig_l - assume mean hyperparameters for regional trtmt effects are same across regions
  Sig_l <- matrix(0, nrow = T_l + p_cov, ncol = T_l + p_cov)
  Sig_l11 <- Sig0[ 1:T_l, 1:T_l ]
  
  if(p_cov > 0){
    
    Sig_l12 <- Sig0[ 1:T_l, (S+1):ncol(Sig0) ]
    Sig_l22 <- Sig0[ (S+1):ncol(Sig0), (S+1):ncol(Sig0) ]
    Sig_l <- cbind( rbind(Sig_l11, Sig_l12), rbind(Sig_l12.t(), Sig_l22) )
    
  } else{
    
    Sig_l <- Sig_l11
    
  }
  
  return(as.matrix(Sig_l))
  
}



### Function to construct W.l (updated design matrix)
update_W_l <- function(W, modMat.l){
  
  # Extract necessary values
  N <- nrow(W)                               # number of subjects
  S <- length(modMat.l)                      # number of regions
  unq.parms.l <- sort(unique(modMat.l))
  T.l <- length(unq.parms.l)                 # number of distinct region labels
  
  # Construct W.l (updated design matrix)
  W.lt <- matrix(0, nrow = N, ncol = T.l)    # columns indicate distinct regional trtmt assigments
  
  for(t in 1:T.l){
    for(i in 1:S){
      for(k in 1:N){
        if(modMat.l[i] == unq.parms.l[t] & W[k,i] == 1){
          W.lt[k,t] <- 1
        }
      }
    }
  }
  
  return(W.lt)
  
}



### Function to perform Bayesian model averaging
bma_fun <- function(values_vec, post_mod_probs){
  bma_result <- sum( values_vec * post_mod_probs )
  return(bma_result)
}



### Function to calculate pairwise consistency/inconsistency probabilities for model M_l, l=1,...,L
# S: number of regions
# modMat_l: S x 1 vector of region classification labels for model M_l
# n_draws: number of samples to draw from posterior distributions
# rte_means_l: T.l x 1 vector of posterior means of distinct region-specific treatment effects
# rte_sds_l: T.l x 1 vector of posterior sd`s of distinct region-specific treatment effects
# epsilon_star: vector of possible minimal clinically important differences of trtmnt effects
pairwise_consistency_l <- function(S, modMat_l, n_draws, rte_means_l, rte_sds_l, epsilon_star){
  
  # Determine number of distinct treatment effects for model M_l
  T.l <- length(rte_means_l)
  
  # For model l, create list of matrices storing pairwise consistency/inconsistency probs.
  #    (separate matrix for all q epsilon.star values)
  prob.prws.consis.l <- list(length(epsilon_star))
  prob.prws.inconsis.l <- list(length(epsilon_star))
  for(q in 1:length(epsilon_star)){
    prob.prws.consis.l[[q]] <- matrix(1, nrow = S, ncol = S)
    prob.prws.inconsis.l[[q]] <- matrix(1, nrow = S, ncol = S)
  }
  
  # Save values for each pairwise comparison that shares the t^th and w^th distinct treatment
  # effects under model M_l
  for(t in 1:T.l){
    
    # Sample n.draws distinct region-specific treatment effects gamma_(l,t)
    samp.l.t <- rnorm( n_draws, rte_means_l[t], rte_sds_l[t] )
    
    for(w in 1:T.l){
      
      # Sample n.draws distinct region-specific treatment effects gamma_(l,w)
      samp.l.w <- rnorm( n_draws, rte_means_l[w], rte_sds_l[w] )
      
      # Calculate |gamma_(l,t) - gamma_(l,w)| for regions that share the t^th and
      # w^th distinct treatment effects
      prws.abs.diff <- abs(samp.l.t - samp.l.w)
      
      # Identify which regions share i^th and w^th distinct treatment effects
      for(i in 1:(S-1)){
        for(j in (i+1):S){
          
          if(modMat_l[i] == t){
            if(modMat_l[j] == w){
              
              # Calculate Pr(|gamma_(l,t) - gamma_(l,w)| < -log(epsilon.star)|D,M_l)
              # for all epsilon.star values
              for(q in 1:length(epsilon_star)){
                
                # extract q^th matrix for model l to store pw consistency/inconsistency probabilities
                pwc.prob.mat.l.q <- prob.prws.consis.l[[q]]
                pwic.prob.mat.l.q <- prob.prws.inconsis.l[[q]]
                # vector to store indicators for pw consis. prob.
                less.epsilon.yes <- ifelse(prws.abs.diff < -log(epsilon_star[q]), 1, 0)
                
                # Calculate pairwise consistency probabilities between regions i and j for model l
                pwc.prob.mat.l.q[i,j] <- mean(less.epsilon.yes)
                pwc.prob.mat.l.q[j,i] <- mean(less.epsilon.yes)
                prob.prws.consis.l[[q]] <- pwc.prob.mat.l.q
                
                # Calculate pairwise inconsistency probabilities between regions i and j for model l
                pwic.prob.mat.l.q[i,j] <- 1 - mean(less.epsilon.yes)
                pwic.prob.mat.l.q[j,i] <- 1 - mean(less.epsilon.yes)
                prob.prws.inconsis.l[[q]] <- pwic.prob.mat.l.q
                
              }
              
            }
          }
          
        }
      }
      
    }
  }
  
  # Return S x S matrices with pairwise comparisons as list for model M_l
  prob.prws.comparisons.l <- list(2)
  prob.prws.comparisons.l[[1]] <- prob.prws.consis.l
  prob.prws.comparisons.l[[2]] <- prob.prws.inconsis.l
  return(prob.prws.comparisons.l)
  
}



### Function to calculate pairwise consistency and inconsistency probabilities
# S: number of regions
# pmp: L x 1 vector of posterior model probabilities
# prws_consis_list: list with all pairwise consistency matrices for each epsilon.star
#    value and each model
# prws_inconsis_list: list with all pairwise consistency matrices for each epsilon.star
#    value and each model
# epsilon_star: vector of possible minimal clinically important differences of trtmnt effects
pairwise_consistency <- function(S, pmp, prws_consis_list, prws_inconsis_list, epsilon_star){
  
  # Determine number of models
  numMods <- length(pmp)
  
  # List to store BMA pairwise consistency results (one matrix for each epsilon_star value)
  bma.prws.consis <- list(length(epsilon_star))
  bma.prws.inconsis <- list(length(epsilon_star))
  for(q in 1:length(epsilon_star)){
    
    bma.prws.consis.q <- matrix(0, nrow = S, ncol = S)
    bma.prws.inconsis.q <- matrix(0, nrow = S, ncol = S)
    for(i in 1:S){
      for(j in 1:S){
        
        # Extract Pr(|gamma_i - gamma_j| < -log(epsilon_star)|D, M_l) and
        # Pr(|gamma_i - gamma_j| > -log(epsilon_star)|D, M_l) for all L models
        pwc.probs.ij <- numeric(numMods)
        pwic.probs.ij <- numeric(numMods)
        for(l in 1:numMods){
          prws_consis_list.l <- prws_consis_list[[l]]
          prws_inconsis_list.l <- prws_inconsis_list[[l]]
          prws.consis.l.q <- prws_consis_list.l[[q]]
          prws.inconsis.l.q <- prws_inconsis_list.l[[q]]
          pwc.probs.ij[l] <- prws.consis.l.q[i,j]
          pwic.probs.ij[l] <- prws.inconsis.l.q[i,j]
        }
        
        # Posterior Pr(|gamma_i - gamma_j| < -log(epsilon_star)|D) and
        # posterior Pr(|gamma_i - gamma_j| > -log(epsilon_star)|D) 
        bma.prws.consis.q[i,j] <- rbind(pwc.probs.ij) %*% cbind(pmp)
        bma.prws.inconsis.q[i,j] <- rbind(pwic.probs.ij) %*% cbind(pmp)
        
      }
    }
    
    # Save matrix of pairwise consistency probabilities for q^th value of epsilon_star
    bma.prws.consis[[q]] <- bma.prws.consis.q
    bma.prws.inconsis[[q]] <- bma.prws.inconsis.q
    
  }
  
  # Return S x S matrices with pairwise comparisons as list for model M_l
  bma.prws.comps <- list(2)
  bma.prws.comps[[1]] <- bma.prws.consis
  bma.prws.comps[[2]] <- bma.prws.inconsis
  return(bma.prws.comps)
  
}



### Function to calculate epsilon.star-level global consistency/inconsistency probabilities
# S: number of regions
# modMat: L x S matrix of region classification labels for all L models
# pmp: L x 1 vector of posterior model probabilities
# bma_prws_consis: list of pairwise consistency matrices for all epsilon.star values
# bma_prws_inconsis: list of pairwise inconsistency matrices for all epsilon.star values
# epsilon_star: vector of possible minimal clinically important differences of trtmnt effects
# beta_star: probability cutoff for global inconsistency for which two regions are considered
#      to be clinically different
global_consistency <- function(S, modMat, pmp, bma_prws_consis, bma_prws_inconsis, epsilon_star, beta_star){
  
  # Determine number of models
  numMods <- length(pmp)
  
  # Determine if regional treatment effects are greater than a clinically meaningful regional
  # difference of one another. If so, determine which models allow these regions to differ and
  # store results in which.mods.keep (1 if model should be considered, 0 otherwise)
  which.mods.keep <- matrix(0, nrow = numMods, ncol = length(epsilon_star))
  for(q in 1:length(epsilon_star)){
    
    pwic.mat.q <- bma_prws_inconsis[[q]]  # Extract matrix of pwc probs. for q^th epsilon_star value
    for(l in 1:numMods){
      
      i.loop.continue <- 1                # Continue with i for-loop if equal to 1; break if 0
      for(i in 1:S){
        
        for(j in 1:S){
          
          # If Pr(|gamma_i - gamma_j| > -log(epsilon_star) | D) < beta_star, consider
          # Model l in calculation for global inconsistency only if model forces regions i and j
          # to differ
          
          if(pwic.mat.q[i,j] >= beta_star){
            if(modMat[l,i] != modMat[l,j]){
              
              # Consider Model l when adding PMPs for global inconsistency probability
              i.loop.continue <- 0
              which.mods.keep[l,q] <- 1
              break
              
            }
          }
          
        }
        
        # If j for-loop is stopped early and model M_l is considered, then break i for-loop
        if(i.loop.continue == 0)  break
        
      }
      
    }
    
  }
  
  # For q^th value of epsilon_star, calculate epsilon_star-level global inconsistency probability
  # for models which require any two regions to differ by more than epsilon_star, with
  # probability beta_star. The epsilon_star-level global consistency probability is
  # defined as 1 - epsilon_star-level global inconsistency probability.
  glob.inconsis.prob <- t(which.mods.keep) %*% cbind(pmp)
  glob.consis.prob <- 1.0 - glob.inconsis.prob
  
  # Return S x S vector with pairwise comparisons as list for model M_l
  return(glob.consis.prob)
  
}

