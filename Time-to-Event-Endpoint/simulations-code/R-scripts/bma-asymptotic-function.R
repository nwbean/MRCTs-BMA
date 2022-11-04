############################################################################################
# BMA Approach for Time-to-Event Endpoints Using Asymptotic Estimation of Posterior
# Distributions of Region-Specific Treatment Effects
#
# Piecewise Exponential Model
############################################################################################


### Functions for marginal posterior distributions of regression effects

## Approximate h(theta_l) function, where h(theta_l) = log( p(theta_l|D, M_l) )
##   i.e., the log of the full conditional distribution of theta_l, l = 1,...,L
# x: place holder for (T_l + p)-dimensional vector theta_l
# S: number of regions
# mu.l: prior mean vector for theta_l
# Sig.l: prior covariance matrix for theta_l
# phi.0: (S x K) matrix of phi.ik values, i=1,...,S, k=1,...,K
# eta.tilde: (S x K) matrix with values of eta.tilde.ik
# phi.tilde.summand.list: list with all S (n.i x K) matrices of summands for phi.tilde.ik
# dat.all.regs: full dataset for all subjects from all regions
# W.l: (N x (T_l + p)) design matrix
h.theta.l <- function( x, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde, phi.tilde.summand.list,
                       dat.all.regs, W.l ){
  
  # Verify that x is saved as a (T_l + p) x 1 matrix
  x <- cbind(x)
  
  # Calculate phi.tilde.ik for all i and k
  # Also calculate log of exp( sum_j^n.i delta_ijk * nu_ij * w_ijl' theta_l ) for all i
  K <- ncol(phi.0)
  phi.tilde <- matrix( 0, nrow = nrow(phi.0), ncol = ncol(phi.0) )
  log.exp.piece <- matrix( 0, nrow = nrow(phi.0), ncol = ncol(phi.0) )
  for(i in 1:S){
    
    # phi.tilde.ik
    reg.i <- sort( unique(dat.all.regs$region) )[i]
    W.il <- as.matrix( W.l[ which(dat.all.regs$region == reg.i), ] )
    phi.tilde[i,] <- phi.0[i,] + colSums( phi.tilde.summand.list[[i]] * as.numeric(exp(W.il %*% x)) )
    
    # log of exp( sum_j^n.i delta_ijk * nu_ij * w_ijl' theta_l )
    dat.i <- dat.all.regs[ which(dat.all.regs$region == reg.i), ]
    for(k in 1:K){
      log.exp.piece[i,k] <- sum( dat.i$status[ which(dat.i$interval == k) ] *
                                   as.numeric(W.il[ which(dat.i$interval == k), ] %*% x) )
    }
    
  }
  
  # Log of full conditional of theta_l
  log.val <- sum( eta.0 * log(phi.0) - lgamma(eta.0) + lgamma(eta.tilde) -
                    eta.tilde * log(phi.tilde) + log.exp.piece ) -
    .5*length(x)*log(2*pi) - .5*log(det(Sig.l)) - .5*t(x - mu.l) %*% solve(Sig.l) %*% (x - mu.l)
  
  return(log.val)
  
}


## Negative of h(theta_l) function
# See function "h.theta.l" for function input descriptions
neg.h.theta.l <- function( x, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde, phi.tilde.summand.list,
                           dat.all.regs, W.l ){
  val <- h.theta.l( x, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde, phi.tilde.summand.list,
                    dat.all.regs, W.l )
  return( -1 * val )
}


## Gradient of h(theta_l) function
# See function "h.theta.l" for function input descriptions
h.p.theta.l <- function( x, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde, phi.tilde.summand.list,
                         dat.all.regs, W.l ){
  
  # Verify that x is saved as a (T_l + p) x 1 matrix
  x <- cbind(x)
  
  # First piece
  p1 <- -solve(Sig.l) %*% (x - mu.l)
  
  # Second piece
  p2.i <- matrix(0, nrow = length(x), ncol = S)
  for(i in 1:S){
    reg.i <- sort( unique(dat.all.regs$region) )[i]
    c.hat.i <- phi.tilde.summand.list[[i]]
    W.il <- as.matrix( W.l[ which(dat.all.regs$region == reg.i), ] )
    p2.i[,i] <- t( rowSums( t( as.numeric( colSums( as.numeric(exp(W.il %*% x)) * c.hat.i ) +
                                             phi.0[i,] )^-1 * as.numeric(eta.tilde[i,]) *
                                 t( t(W.il) %*% ( as.numeric(exp(W.il %*% x)) * c.hat.i ) ) ) ) )
  }
  p2 <- rowSums(p2.i)
  
  # Third piece
  p3 <- colSums( W.l[ which(dat.all.regs$status == 1), ] )
  
  # Gradient of h(gamma_l) function
  grad <- p1 - p2 + p3
  return(grad)
  
}


## Hessian of h(theta_l) function
# See function "h.theta.l" for function input descriptions
h.pp.theta.l <- function( x, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde, phi.tilde.summand.list,
                          dat.all.regs, W.l ){
  
  # Verify that x is saved as a (T_l + p) x 1 matrix
  x <- cbind(x)
  
  # First piece
  p1 <- -solve(Sig.l)
  
  # Second piece
  p2 <- matrix(0, nrow = ncol(W.l), ncol = ncol(W.l))
  K <- ncol(phi.0)
  for(i in 1:S){
    reg.i <- sort( unique(dat.all.regs$region) )[i]
    W.il <- as.matrix( W.l[ which(dat.all.regs$region == reg.i), ] )
    for(k in 1:K){
      c.hat.ik <- phi.tilde.summand.list[[i]][,k]
      p2 <- p2 + eta.tilde[i,k] *
        ( (t( c.hat.ik * as.numeric(exp(W.il %*% x)) * W.il ) %*% W.il) *
            (phi.0[i,k] + sum(exp(W.il %*% x) * c.hat.ik) )^-1 -
            cbind( colSums( c.hat.ik * as.numeric(exp(W.il %*% x)) * W.il ) ) %*%
            rbind( colSums( c.hat.ik * as.numeric(exp(W.il %*% x)) * W.il ) ) *
            (phi.0[i,k] + sum( exp(W.il %*% x) * c.hat.ik ) )^-2 )
    }
  }
  
  # Hessian of h(gamma_(l,t)) function
  hess <- p1 - p2
  return(hess)
  
}


## Approximated log of marginal distribution of data given model M_l
# theta.mode: posterior mode of pr(theta|D,M_l)
# See function "h.theta.l" for function remaining input descriptions
log.p.D.Ml <- function( theta.mode, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde,
                        phi.tilde.summand.list, dat.all.regs, W.l ){
  
  # Negative inverse of Hessian matrix
  psi.mat <- -solve( h.pp.theta.l( theta.mode, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde,
                                   phi.tilde.summand.list, dat.all.regs, W.l ) )
  
  # h(theta_l) evaluated at posterior mode of theta
  log.post.theta <- h.theta.l( theta.mode, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde,
                               phi.tilde.summand.list, dat.all.regs, W.l )
  
  # Approximated log of marginal likelihood
  log.val <- .5 * length(theta.mode) * log(2*pi) + .5 * log(det(psi.mat)) + log.post.theta
  return(log.val)
  
}



### Function to run BMA algorithm using Laplace approximation for a single dataset
# Returns posterior summaries, approximated PMPs, and consistency measures
# dat: full dataset for all regions combined
# W.mat: design matrix with region-by-treatment indicators and optional covariates
# int.cuts: K x 1 vector of lower interval boundaries for all K intervals
# mu0: prior mean vector for regression effects (region-specific treatment effects and covariates)
# Sig0: prior covariance matrix for regression effects
# eta0: S x K matrix of shape hyperparameters for prior on baseline hazard
# phi0: S x K matrix of rate hyperparameters for prior on baseline hazard
# mod.mat: L x S matrix of possible region groupings, each row corresponding to one of L models
# mod.priors: L x 1 vector of prior model probabilities
# gamma0: value of interest in null hypothesis to be tested
# n.draws: number of posterior samples to draw when evaluating consistency
# epsilon.star: vector of possible minimal clinically important regional differences of trtmnt
#      effects - used for consistency measures
# beta.star: probability cutoff for global inconsistency for which two regions are considered
#      to be clinically different
# pi0: value for which we want to calculate probabilities with local consistency
bma.asymptotic <- function( dat, W.mat, int.cuts, mu0, Sig0, eta0, phi0, mod.mat,
                            mod.priors, gamma0, n.draws, epsilon.star, beta.star, pi0 ){

  
  # Derive additional values based on function inputs
  S <- length( unique(dat$region) )     # number of regions
  N <- nrow(dat)                        # total sample size
  n.i <- table(dat$region)              # regional sample size
  K <- length(int.cuts)                 # number of intervals
  num.covs <- ncol(W.mat) - S           # number of covariates
  
  
  # Calculate eta.tilde.ik for all i and k (value remains constant across models and sample parameters).
  # Also calculate ind.sum.mat, where ind.sum.mat_ik = sum_(j=1)^n.i (d_ijk * nu_ij)
  # (i.e., the number of patients who received treatment from region i who failed in interval k)
  eta.tilde <- matrix( 0, nrow = S, ncol = K )
  for(i in 1:S){
    for(k in 1:K){
      d.ijk.times.nu.ij <- ifelse( dat$region == i & dat$status == 1 &
                                     dat$interval == k, 1, 0 )
      eta.tilde[i,k] <- sum( d.ijk.times.nu.ij ) + eta0[i,k]
    }
  }
  
  
  # Calculate delta_ijk * (y_ij - m_{k-1}) + sum_(g=k+1)^K( delta_ijg * (m_k - m_{k-1}) )
  # for all subjects in region i for each time interval k. Save (n.i x K) matrix for each
  # region in a list. These summands are later used to construct phi.tilde.ik.
  phi.tilde.summand.list <- list()
  for(i in 1:S){
    
    dat.i <- dat[ which(dat$region == i), ]    # data for only region i
    
    # matrix to store summand for all subjects in region i for each time interval k, k=1,...,K
    summand.mat.i <- matrix( 0, nrow = nrow(dat.i), ncol = K )
    for(k in 1:K){
      summand.val <- ifelse( dat.i$interval == k, dat.i$eventtime - int.cuts[k], 0 )
      summand.val <- ifelse( dat.i$interval > k, int.cuts[k+1] - int.cuts[k], summand.val )
      summand.mat.i[,k] <- summand.val
    }
    phi.tilde.summand.list[[i]] <- summand.mat.i
    
  }
  
  
  # Vector to store log marginal likelihoods for each model
  L <- nrow(mod.mat)        # number of models in model space
  log.p.D <- numeric(L)
  
  
  # Matrices to hold model-specific posterior summary statistics for each region
  rte.mean.l.mat <- matrix( 0, nrow = L, ncol = S )
  rte.sd.l.mat <- matrix( 0, nrow = L, ncol = S )
  rte.gamma0.prob.l.mat <- matrix( 0, nrow = L, ncol = S )
  
  
  # Matrices to hold model-specific posterior summary statistics for global treatment effect
  gte.mean.l <- numeric(L)
  gte.sd.l <- numeric(L)
  gte.less.gamma0.l <- numeric(L)
  
  
  # Matrices to hold model-specific posterior summary statistics for covariate effects
  beta.mean.l.mat <- matrix( 0, nrow = L, ncol = num.covs )
  beta.sd.l.mat <- matrix( 0, nrow = L, ncol = num.covs )
  
  
  # Matrices to hold consistency results
  loc.consis.loo.l <- list()    # list to store matrices (one for each epsilon.star) for LOO loc. cons.
  for(q in 1:length(epsilon.star)){
    loc.consis.loo.l[[q]] <- matrix( 0, nrow = L, ncol = S )
  }
  pmda.TE.loc.cons.l <- matrix( 0, nrow = L, ncol = S )   # PMDA local consistency for all models (TE)
  pmda.RR.loc.cons.l <- matrix( 0, nrow = L, ncol = S )   # PMDA local consistency for all models (RR)
  prws.consis.list <- list()    # list to store pairwise consistency results
  prws.inconsis.list <- list()  # list to store pairwise inconsistency results
  
  
  # Loop through models
  for(l in 1:L){
    
    # Details for model M_l
    mod.l <- mod.mat[l,]               # region groupings for model l
    T.l <- length( unique(mod.l) )     # number of distinct treatment effects T_l
    n.l <- numeric(T.l)                # vector to store combined sample sizes
    num.regs.l <- numeric(T.l)         # vector to store number of regions in each group
    rte.means.l <- numeric(T.l)        # vector to store distinct gamma_(l,t) means
    rte.var.l <- numeric(T.l)          # vector to store distinct gamma_(l,t) variances
    log.ml.cntrbtn <- numeric(T.l)     # contribution of each region to marginal likelihood
    
    # Update dimensions of mean vector and covariance matrix hyperparameters in prior
    # distributions, and update design matrix W.l
    mu0.l <- update_mu_l( mu0, mod.mat[l,] )
    Sig0.l <- update_Sig_l( Sig0, mod.mat[l,] )
    W.l <- update_W_l( W.mat, mod.mat[l,] )
    
    # Calculate posterior mode of p(theta_l|D, M_l)
    # Nelder-Mead method (default) is robust and nearly as fast as L-BFGS-B method
    if( ncol(W.l) > 1 ){
      theta.mode.l <- optim( par = numeric(ncol(W.l)), fn = neg.h.theta.l,
                             method = "Nelder-Mead", S = S, mu.l = mu0.l, Sig.l = Sig0.l,
                             phi.0 = phi0, eta.0 = eta0, eta.tilde = eta.tilde,
                             phi.tilde.summand.list = phi.tilde.summand.list,
                             dat.all.regs = dat, W.l = W.l )$par
    } else{
      theta.mode.l <- optim( par = numeric(ncol(W.l)), fn = neg.h.theta.l, method = "Brent",
                             lower = -100, upper = 100, S = S, mu.l = mu0.l, Sig.l = Sig0.l,
                             phi.0 = phi0, eta.0 = eta0, eta.tilde = eta.tilde,
                             phi.tilde.summand.list = phi.tilde.summand.list,
                             dat.all.regs = dat, W.l = W.l )$par
    }
    
    # Calculate posterior covariance matrix of theta_l|D,M_l
    theta.cov.l <- -h.pp.theta.l( x = theta.mode.l, S = S, mu.l = mu0.l, Sig.l = Sig0.l,
                                  phi.0 = phi0, eta.0 = eta0, eta.tilde = eta.tilde,
                                  phi.tilde.summand.list = phi.tilde.summand.list,
                                  dat.all.regs = dat, W.l = W.l )^(-1)
    
    # Calculate log of marginal likelihood for model M_l
    log.p.D[l] <- log.p.D.Ml( theta.mode = theta.mode.l, S = S, mu.l = mu0.l, Sig.l = Sig0.l,
                              phi.0 = phi0, eta.0 = eta0, eta.tilde = eta.tilde,
                              phi.tilde.summand.list = phi.tilde.summand.list,
                              dat.all.regs = dat, W.l = W.l )
    # Calculate combined sample size, number of events, and sum of event times (y)
    # for regions that share the t^th distinct trtmt effect
    which.regs.l <- list()        # list to store groups of region labels
    for(t in 1:T.l){
      
      which.regs.Tl <- NULL     # turn into vector to store which regions are the same
      for(i in 1:S){
        if(mod.l[i] == t){
          
          # Combined regional values based on groupings for model M_l
          n.l[t] = n.l[t] + n.i[i]
          num.regs.l[t] <- num.regs.l[t] + 1
          which.regs.Tl <- c( which.regs.Tl, i )
          
        }
        
        # Store labels of regions that are grouped together
        which.regs.l[[t]] <- which.regs.Tl
        
      }
    }
    
    
    # Calculate summary stats for posterior distribution of gamma_(l,t) and store stats
    # with regions that share the treatment effect gamma_(l,t), t=1,...,T_l, l=1,...,L
    for(t in 1:T.l){
      
      # Save approximated posterior mean and variance for regional treatment effects
      rte.means.l[t] <- theta.mode.l[t]
      rte.var.l[t] <- theta.cov.l[t,t]
      
      # Calculate Pr(gamma_(l,t) < gamma0|D, M_l)
      rte.prob.less.gamma0 <- pnorm( gamma0, rte.means.l[t], sqrt(rte.var.l[t]) )
      
      # Store posterior values for each region
      for(k in 1:S){
        if(mod.l[k] == t){
          
          # Posterior summary stats
          rte.mean.l.mat[l,k] <- rte.means.l[t]
          rte.sd.l.mat[l,k] <- sqrt(rte.var.l[t])
          rte.gamma0.prob.l.mat[l,k] <- rte.prob.less.gamma0
          
        }
      }
      
    }
    
    
    # Calculate global treatment effect summary stats for model M_l
    gte.mean.l[l] <- sum( rte.means.l * n.l/N )
    gte.sd.l[l] <- sqrt( sum( (n.l/N)^2 * rte.var.l ) )
    gte.less.gamma0.l[l] <- pnorm( gamma0, gte.mean.l[l], gte.sd.l[l] )
    
    
    # Store summary stats for covariates for model M_l
    if(num.covs > 0){
      beta.mean.l.mat[l,] <- theta.mode.l[ (T.l+1):(T.l+num.covs) ]
      beta.sd.l.mat[l,] <- sqrt( diag(theta.cov.l)[ (T.l+1):(T.l+num.covs) ] )
    }
    
    
    # Obtain Pr(gamma_(l,t)/gamma > pi0|D,M_l) for model M_l (PMDA local consistency with
    # treatment effects) and Pr( (1 - exp{gamma_(l,t)})/(1 - exp{gamma}) > pi0|D,M_l) for model M_l
    # (PMDA local consistency with risk reduction)
    for(t in 1:T.l){
      
      # Calculate Pr(gamma_(l,t)/gamma > pi0|D,M_l) for a given value of t
      pmda.TE.loc.cons.lt <- mean( rnorm(10000, rte.means.l[t], sqrt(rte.var.l[t])) /
                                     rnorm(10000, gte.mean.l[l], gte.sd.l[l]) > pi0 )
      
      # Calculate Pr( (1 - exp{gamma_(l,t)})/(1 - exp{gamma}) > pi0|D,M_l) for a given value of t
      pmda.RR.loc.cons.lt <- mean( (1 - exp( rnorm(10000, rte.means.l[t], sqrt(rte.var.l[t])) )) /
                                     (1 - exp( rnorm(10000, gte.mean.l[l], gte.sd.l[l]) )) > pi0 )
      
      # Store posterior values for each region
      for(k in 1:S){
        if(mod.l[k] == t){
          
          # Posterior summary stats
          pmda.TE.loc.cons.l[l,k] <- pmda.TE.loc.cons.lt
          pmda.RR.loc.cons.l[l,k] <- pmda.RR.loc.cons.lt
          
        }
      }
      
    }
    
    
    # Calculate leave-one-out global treatment effect
    gte.means.loo.l <- numeric(S)
    gte.sd.loo.l <- numeric(S)
    for(i in 1:S){
      t.loo <- mod.l[i]       # classification of region to leave out under M_l
      num.regs.loo <- sum( mod.l == t.loo )
      if( num.regs.loo == 1 ){
        rte.means.loo <- rte.means.l[-t.loo]
        rte.vars.loo <- rte.var.l[-t.loo]
        n.l.loo <- n.l[-t.loo]
        N.loo <- sum(n.l.loo)
        gte.means.loo.l[i] <- sum( rte.means.loo * n.l.loo/N.loo )
        gte.sd.loo.l[i] <- sqrt( sum( rte.vars.loo * (n.l.loo/N.loo)^2 ) )
      } else{
        
        # Update data to account for removed i^th region
        dat.loo <- dat[ which(dat$region != i), ]                  # update dat
        N.loo <- N - n.i[i]                                        # update N
        n.l.loo <- n.l
        n.l.loo[t.loo] <- n.l[t.loo] - n.i[i]                      # update n.i
        W.l.loo <- cbind(W.l[ which(dat$region != i), ])           # update W.l
        eta.tilde.loo <- cbind(eta.tilde[-i,])                     # update eta.tilde.loo
        phi.tilde.summand.list.loo <- phi.tilde.summand.list[-i]   # update phi.tilde.summand.list
        
        # Calculate posterior mode of p(theta_l|D, M_l) without the i^th region
        if( ncol(W.l.loo) > 1 ){
          theta.l.loo <- optim( par = numeric(ncol(W.l.loo)), fn = neg.h.theta.l,
                                method = "Nelder-Mead", S = S-1, mu.l = mu0.l, Sig.l = Sig0.l,
                                phi.0 = cbind(phi0[-i,]), eta.0 = cbind(eta0[-i,]), eta.tilde = eta.tilde.loo,
                                phi.tilde.summand.list = phi.tilde.summand.list.loo,
                                dat.all.regs = dat.loo, W.l = W.l.loo )$par
        } else{
          theta.l.loo <- optim( par = numeric(ncol(W.l.loo)), fn = neg.h.theta.l, method = "Brent",
                                lower = -100, upper = 100, S = S-1, mu.l = mu0.l, Sig.l = Sig0.l,
                                phi.0 = cbind(phi0[-i,]), eta.0 = cbind(eta0[-i,]), eta.tilde = eta.tilde.loo,
                                phi.tilde.summand.list = phi.tilde.summand.list.loo,
                                dat.all.regs = dat.loo, W.l = W.l.loo )$par
        }
        
        # Calculate posterior covariance matrix of theta_l|D,M_l
        theta.cov.l.loo <- -h.pp.theta.l( x = theta.l.loo, S = S-1, mu.l = mu0.l, Sig.l = Sig0.l,
                                          phi.0 = phi0, eta.0 = eta0, eta.tilde = eta.tilde.loo,
                                          phi.tilde.summand.list = phi.tilde.summand.list.loo,
                                          dat.all.regs = dat.loo, W.l = W.l.loo )^(-1)
        
        rte.means.loo <- theta.l.loo[1:T.l]
        rte.vars.loo <- diag(theta.cov.l.loo)[1:T.l]
        gte.means.loo.l[i] <- sum( rte.means.loo * n.l.loo/N.loo )
        gte.sd.loo.l[i] <- sqrt( sum( rte.vars.loo * (n.l.loo/N.loo)^2 ) )
      }
    }
    
    
    # Obtain Pr(epsilon.star < exp{gamma_(l,t) - gamma_(-i)} < 1/epsilon.star|D,M_l) conditional
    # on M_l, l=1,...,L
    # Equivalent to Pr(|gamma_(l,t) - gamma_(-i)| < -log(epsilon.star)|D,M_l)
      
    # Save values for each region under model M_l
    consis.diff.loo <- matrix(0, nrow = n.draws, ncol = S)
    for(i in 1:S){
      
      # Identify distinct treatment effect number for i^th region
      t.i <- mod.l[i]
      
      # Obtain samples of gamma_(l,t) and gamma_(-i), where i is in the set Omega_{l,t}
      gamma.ti.samp <- rnorm(n.draws, rte.means.l[t.i], sqrt(rte.var.l[t.i]))
      gamma.glob.no.i <- rnorm(n.draws, gte.means.loo.l[i], gte.sd.loo.l[i])
      
      # Calculate difference |gamma_(l,t) - gamma_(-i)| for each region
      consis.diff.loo[,i] <- abs( gamma.ti.samp - gamma.glob.no.i )
      
    }
    
    # Calculate probability that difference |gamma_(l,t) - gamma_(-i)| is less than
    # -log(epsilon.star) for each value of epsilon.star
    for(q in 1:length(epsilon.star)){
      loc.consis.loo.l[[q]][l,] <- colMeans( consis.diff.loo < -log(epsilon.star[q]) )
    }
    
    
    # Obtain pairwise consistency probabilities Pr(|gamma_i - gamma_j| < -log(epsilon_star)|D,M_l)
    # conditional on M_l, l=1,...,L, for all pairwise comparisons between regions and all
    # specified values of epsilon_star. Also obtain pairwise inconsistency probabilities
    # Pr(|gamma_i - gamma_j| > -log(epsilon_star)|D,M_l).
    
    # For model l, create list of matrices storing pairwise consistency/inconsistency
    #    probabilities (separate matrix for all q epsilon_star values)
    prws.comps.list <- pairwise_consistency_l(S = S, modMat_l = mod.mat[l,], n_draws = n.draws,
                                              rte_means_l = rte.means.l, rte_sds_l = sqrt(rte.var.l),
                                              epsilon_star = epsilon.star)
    prws.consis.list[[l]] <- prws.comps.list[[1]]
    prws.inconsis.list[[l]] <- prws.comps.list[[2]]
    
  }
  
  
  # Calculate posterior model probability for each model
  max.log.p.D <- max(log.p.D)
  ml.prop <- exp( log.p.D - max.log.p.D )
  pmp.vec <- ( ml.prop * mod.priors ) / sum( ml.prop * mod.priors )
  
  
  # Obtain n.samp samples from posterior distributions of global and region-specific treatment
  # effects (mixture of normal distributions), and obtain 95% central credible intervals
  cred.ints <- matrix(0, nrow = S + 1, ncol = 2)      # store bounds of 95% central credible ints.
  rownames(cred.ints) <- c("Global", paste("Reg", 1:S, sep = ""))
  colnames(cred.ints) <- c("LowBound", "UpBound")
  n.samp <- 100000       # number of draws to sample from posterior distributions
  which.mods.glob <- sample(1:L, n.samp, pmp.vec, replace = TRUE)   # sample model indicators
  glob.samp <- rnorm(n.samp, gte.mean.l[which.mods.glob], gte.sd.l[which.mods.glob])
  cred.ints[1,] <- quantile(glob.samp, c(.025, .975))
  for(i in 1:S){
    which.mods.reg <- sample(1:L, n.samp, pmp.vec, replace = TRUE)   # sample model indicators
    reg.i.samp <- rnorm(n.samp, rte.mean.l.mat[which.mods.reg, i], rte.sd.l.mat[which.mods.reg, i])
    cred.ints[i+1,] <- quantile(reg.i.samp, c(.025, .975))
  }
  
  
  # Calculate epsilon.star-level local consistency (leave one out) for all values of epsilon.star
  loc.consis.loo <- matrix(0, nrow = length(epsilon.star), ncol = S)
  row.names(loc.consis.loo) <- epsilon.star
  colnames(loc.consis.loo) <- 1:S
  for(q in 1:length(epsilon.star)){
    loc.consis.loo[q,] <- t(pmp.vec) %*% loc.consis.loo.l[[q]]
  }
  
  
  # Calculate epsilon.star-level pairwise consistency for all values of epsilon.star
  bma.prws.comps <- pairwise_consistency( S = S, pmp = pmp.vec,
                                          prws_consis_list = prws.consis.list,
                                          prws_inconsis_list = prws.inconsis.list,
                                          epsilon_star = epsilon.star )
  bma.prws.consis <- bma.prws.comps[[1]]
  bma.prws.inconsis <- bma.prws.comps[[2]]
  
  
  # Calculate epsilon.star-level global consistency/inconsistency probabilities
  # The probability for epsilon.star-level global inconsistency is calculated as the sum of the PMPs
  # for all models which do not allow regional treatment effects to differ by more than epsilon.star
  # with beta.star probability. If Pr(|gamma_i - gamma_j| > -log(epsilon.star)| D) >= beta.star
  # for i != j, then we consider all models in which regions i and j differ. Identify these models
  # for all pairwise comparisons, and then add up the PMPs for all of these models.
  glob.consis.vec <- global_consistency( S = S, modMat = mod.mat, pmp = pmp.vec,
                                         bma_prws_consis = bma.prws.consis,
                                         bma_prws_inconsis = bma.prws.inconsis,
                                         epsilon_star = epsilon.star,
                                         beta_star = beta.star )
  
  
  # Save posterior summaries, PMPs, and consistency measures in list
  BMA.results <- list()
  BMA.results$PMPs <- pmp.vec
  BMA.results$gte.mean <- rbind(pmp.vec) %*% cbind(gte.mean.l)
  BMA.results$gte.sd <- rbind(pmp.vec) %*% cbind(gte.sd.l)
  BMA.results$gte.less.gamma0 <- rbind(pmp.vec) %*% cbind(gte.less.gamma0.l)
  BMA.results$rte.means <- rbind(pmp.vec) %*% rte.mean.l.mat
  BMA.results$rte.sds <- rbind(pmp.vec) %*% rte.sd.l.mat
  BMA.results$rte.less.gamma0 <- rbind(pmp.vec) %*% rte.gamma0.prob.l.mat
  BMA.results$credible.intervals <- cred.ints
  BMA.results$cov.means <- rbind(pmp.vec) %*% beta.mean.l.mat
  BMA.results$cov.sds <- rbind(pmp.vec) %*% beta.sd.l.mat
  BMA.results$local.consistency.ratio.TE <- pmp.vec %*% pmda.TE.loc.cons.l
  BMA.results$local.consistency.ratio.RR <- pmp.vec %*% pmda.RR.loc.cons.l
  BMA.results$local.consistency.loo <- loc.consis.loo
  BMA.results$global.consistency <- glob.consis.vec
  BMA.results$pairwise.consistency <- bma.prws.consis
  BMA.results$pairwise.inconsistency <- bma.prws.inconsis
  
  
  # Return list
  return(BMA.results)

}

