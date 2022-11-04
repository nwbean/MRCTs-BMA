############################################################################################
# BMA Approach for Joint Models with a Time-to-Event Endpoint and a Longitudinal
# Continuous Marker Using Asymptotic Estimation of Posterior Distributions of
# Region-Specific Treatment Effects
#
# Allow for multiple prior model probabilities to be elicited
#
# Piecewise Exponential Baseline Hazards
# Shared Random Effects
############################################################################################


### Load libraries
library(MASS)
library(lme4)
library(survival)


### Function for h(theta_{Y,l}), where h(theta_{Y,l}) = log( p(theta_{Y,l}|xi, b, D, M_l) )
###   i.e., the log of the full conditional distribution of theta_{Y,l}, l = 1,...,L
# x: place holder for (D_{Y,l} + p_Y)-dimensional vector theta_{Y,l}
# S: number of regions
# mu.Yl: prior mean vector for theta_{Y,l}
# Sig.Yl: prior covariance matrix for theta_{Y,l}
# eta.iq: (S x Q) matrix of shape hyperparameters for gamma prior on piecewise constant baseline hazards
# phi.iq: (S x Q) matrix of rate hyperparameters for gamma prior on piecewise constant baseline hazards
# eta.tilde: (S x Q) matrix with values of eta.tilde.iq
# phi.tilde.summand.list: list with all S (n.i x Q) matrices of summands for phi.tilde.iq
# surv.dat: full survival dataset for all subjects from all regions
# W.Yl: (N x (D_{Y,l} + p_Y)) design matrix for survival submodel
# b.mat: N x 1 (or 2) matrix of subj-specific random effects (rand. intercept and optionally rand. slope)
# alpha.vec: 1 (or 2) x 1 vector of association parameters corresponding to shared random effects
h.theta.Yl <- function( x, S, mu.Yl, Sig.Yl, eta.iq, phi.iq, eta.tilde, phi.tilde.summand.list,
                        surv.dat, W.Yl, b.mat, alpha.vec ){

  # Verify that x is saved as a (D_{Y,l} + p_Y) x 1 matrix
  x <- cbind(x)
  alpha.vec <- cbind(alpha.vec)
  W.Yl <- cbind(W.Yl)
  mu.Yl <- as.matrix(mu.Yl)
  Sig.Yl <- as.matrix(Sig.Yl)

  # Calculate phi.tilde.iq for all i and q
  # Also calculate log of exp( sum_j^n.i delta_ijq * nu_ij * w_{Y,lij}' theta_{Y,l} ) for all i
  Q <- ncol(phi.iq)
  phi.tilde <- matrix( 0, nrow = nrow(phi.iq), ncol = ncol(phi.iq) )
  log.exp.piece <- matrix( 0, nrow = nrow(phi.iq), ncol = ncol(phi.iq) )
  for(i in 1:S){

    # phi.tilde.iq
    reg.i <- sort( unique(surv.dat$region) )[i]
    b.i <- as.matrix( b.mat[ which(surv.dat$region == reg.i), ] )   # random effects for reg. i
    W.Yli <- as.matrix( W.Yl[ which(surv.dat$region == reg.i), ] )  # surv. design matrix for reg. i
    exp.lp.i <- as.numeric(exp( b.i %*% alpha.vec + W.Yli %*% x ))  # exp. of linear predictor
    phi.tilde[i,] <- phi.iq[i,] + colSums( phi.tilde.summand.list[[i]] * exp.lp.i )

    # log of exp( sum_j^n.i delta_ijq * nu_ij * w_{Y,lij}' theta_{Y,l} )
    dat.i <- surv.dat[ which(surv.dat$region == reg.i), ]
    for(q in 1:Q){
      log.exp.piece[i,q] <- sum( dat.i$status[ which(dat.i$interval == q) ] *
                                   as.numeric(W.Yli[ which(dat.i$interval == q), ] %*% x) )
    }

  }

  # Log of full conditional of theta_{Y,l}
  log.val <- sum( eta.iq * log(phi.iq) - lgamma(eta.iq) + lgamma(eta.tilde) -
                    eta.tilde * log(phi.tilde) + log.exp.piece ) -
    .5*length(x)*log(2*pi) - .5*log(det(Sig.Yl)) - .5*t(x - mu.Yl) %*% solve(Sig.Yl) %*% (x - mu.Yl)

  return(log.val)

}


### Negative of h(theta_{Y,l}) function
# See function "h.theta.Yl" for function input descriptions
neg.h.theta.Yl <- function( x, S, mu.Yl, Sig.Yl, eta.iq, phi.iq, eta.tilde, phi.tilde.summand.list,
                            surv.dat, W.Yl, b.mat, alpha.vec ){
  val <- h.theta.Yl( x, S, mu.Yl, Sig.Yl, eta.iq, phi.iq, eta.tilde, phi.tilde.summand.list,
                     surv.dat, W.Yl, b.mat, alpha.vec )
  return( -1 * val )
}


### Function for h(alpha), where h(alpha) = log( p(alpha|xi, b, D, M_l) )
###   i.e., the log of the full conditional distribution of alpha
# x: place holder for 1 (or 2) x 1 vector alpha
# S: number of regions
# mu.a: 1 (or 2) x 1 prior mean vector for alpha
# Sig.a: 1 x 1 (or 2 x 2) prior covariance matrix for alpha
# eta.iq: (S x Q) matrix of shape hyperparameters for gamma prior on piecewise constant baseline hazards
# phi.iq: (S x Q) matrix of rate hyperparameters for gamma prior on piecewise constant baseline hazards
# eta.tilde: (S x Q) matrix with values of eta.tilde.iq
# phi.tilde.summand.list: list with all S (n.i x Q) matrices of summands for phi.tilde.iq
# surv.dat: full survival dataset for all subjects from all regions
# W.Yl: (N x (D_{Y,l} + p_Y)) design matrix for survival submodel
# theta.Yl: (D_{Y,l} + p_Y) x 1 vector of regression effects for survival submodel
# b.mat: N x 1 (or 2) matrix of subj-specific random effects (rand. intercept and optionally rand. slope)
h.alpha <- function( x, S, mu.a, Sig.a, eta.iq, phi.iq, eta.tilde, phi.tilde.summand.list,
                     surv.dat, W.Yl, theta.Yl, b.mat ){

  # Verify that x is saved as a 1 x 1 (or 2 x 1) matrix
  x <- cbind(x)
  theta.Yl <- cbind(theta.Yl)

  # Calculate phi.tilde.iq for all i and q
  # Also calculate log of exp( sum_j^n.i delta_ijq * nu_ij * b_ij' alpha ) for all i
  Q <- ncol(phi.iq)
  phi.tilde <- matrix( 0, nrow = nrow(phi.iq), ncol = ncol(phi.iq) )
  log.exp.piece <- matrix( 0, nrow = nrow(phi.iq), ncol = ncol(phi.iq) )
  for(i in 1:S){

    # phi.tilde.iq
    reg.i <- sort( unique(surv.dat$region) )[i]
    b.i <- as.matrix( b.mat[ which(surv.dat$region == reg.i), ] )   # random effects for reg. i
    W.Yli <- as.matrix( W.Yl[ which(surv.dat$region == reg.i), ] )  # surv. design matrix for reg. i
    exp.lp.i <- as.numeric(exp( b.i %*% x + W.Yli %*% theta.Yl ))   # exp. of linear predictor
    phi.tilde[i,] <- phi.iq[i,] + colSums( phi.tilde.summand.list[[i]] * exp.lp.i )

    # log of exp( sum_j^n.i delta_ijq * nu_ij * alpha' b_ij )
    dat.i <- surv.dat[ which(surv.dat$region == reg.i), ]
    for(q in 1:Q){
      log.exp.piece[i,q] <- sum( dat.i$status[ which(dat.i$interval == q) ] *
                                   as.numeric(b.i[ which(dat.i$interval == q), ] %*% x) )
    }

  }

  # Log of full conditional of alpha
  log.val <- sum( eta.iq * log(phi.iq) - lgamma(eta.iq) + lgamma(eta.tilde) -
                    eta.tilde * log(phi.tilde) + log.exp.piece ) -
    .5*length(x)*log(2*pi) - .5*log(det(Sig.a)) - .5*t(x - mu.a) %*% solve(Sig.a) %*% (x - mu.a)

  #log.val <- sum( log.exp.piece )
  return(log.val)

}


### Negative of h(alpha) function
# See function "h.alpha" for function input descriptions
neg.h.alpha <- function( x, S, mu.a, Sig.a, eta.iq, phi.iq, eta.tilde, phi.tilde.summand.list,
                         surv.dat, W.Yl, theta.Yl, b.mat ){
  val <- h.alpha( x, S, mu.a, Sig.a, eta.iq, phi.iq, eta.tilde, phi.tilde.summand.list,
                  surv.dat, W.Yl, theta.Yl, b.mat )
  return( -1 * val )
}


# ### Function for h(alpha, theta_{Y,l}), where h(alpha, theta_{Y,l}) = log( p(alpha, theta_{Y,l}|xi, b, D, M_l) )
# ###   i.e., the log of the full conditional distribution of alpha and theta_{Y,l}, l = 1,...,L
# # x: place holder for (r + D_{Y,l} + p_Y)-dimensional vector c(alpha, theta_{Y,l})
# # S: number of regions
# # mu.a: prior mean vector for alpha
# # Sig.a: prior covariance matrix for alpha
# # mu.Yl: prior mean vector for theta_{Y,l}
# # Sig.Yl: prior covariance matrix for theta_{Y,l}
# # eta.iq: (S x Q) matrix of shape hyperparameters for gamma prior on piecewise constant baseline hazards
# # phi.iq: (S x Q) matrix of rate hyperparameters for gamma prior on piecewise constant baseline hazards
# # eta.tilde: (S x Q) matrix with values of eta.tilde.iq
# # phi.tilde.summand.list: list with all S (n.i x Q) matrices of summands for phi.tilde.iq
# # surv.dat: full survival dataset for all subjects from all regions
# # W.Yl: (N x (D_{Y,l} + p_Y)) design matrix for survival submodel
# # b.mat: N x 1 (or 2) matrix of subj-specific random effects (rand. intercept and optionally rand. slope)
# h.alpha.theta.Yl <- function( x, S, mu.a, Sig.a, mu.Yl, Sig.Yl, eta.iq, phi.iq, eta.tilde,
#                               phi.tilde.summand.list, surv.dat, W.Yl, b.mat ){
#   
#   # Verify that x is saved as a (D_{Y,l} + p_Y) x 1 matrix
#   alpha.vec <- cbind(x[1:length(mu.a)])
#   theta.Yl.vec <- cbind(x[-c(1:length(mu.a))])
#   W.Yl <- cbind(W.Yl)
#   mu.a <- as.matrix(mu.a)
#   Sig.a <- as.matrix(Sig.a)
#   mu.Yl <- as.matrix(mu.Yl)
#   Sig.Yl <- as.matrix(Sig.Yl)
#   
#   # Calculate phi.tilde.iq for all i and q
#   # Also calculate log of exp( sum_j^n.i delta_ijq * nu_ij * w_{Y,lij}' theta_{Y,l} ) for all i
#   Q <- ncol(phi.iq)
#   phi.tilde <- matrix( 0, nrow = nrow(phi.iq), ncol = ncol(phi.iq) )
#   log.exp.piece <- matrix( 0, nrow = nrow(phi.iq), ncol = ncol(phi.iq) )
#   for(i in 1:S){
#     
#     # phi.tilde.iq
#     reg.i <- sort( unique(surv.dat$region) )[i]
#     b.i <- as.matrix( b.mat[ which(surv.dat$region == reg.i), ] )   # random effects for reg. i
#     W.Yli <- as.matrix( W.Yl[ which(surv.dat$region == reg.i), ] )  # surv. design matrix for reg. i
#     exp.lp.i <- as.numeric(exp( b.i %*% alpha.vec + W.Yli %*% theta.Yl.vec ))  # exp. of linear predictor
#     phi.tilde[i,] <- phi.iq[i,] + colSums( phi.tilde.summand.list[[i]] * exp.lp.i )
#     
#     # log of exp( sum_j^n.i delta_ijq * nu_ij * w_{Y,lij}' theta_{Y,l} )
#     dat.i <- surv.dat[ which(surv.dat$region == reg.i), ]
#     for(q in 1:Q){
#       log.exp.piece[i,q] <- sum( dat.i$status[ which(dat.i$interval == q) ] *
#                                    as.numeric(W.Yli[ which(dat.i$interval == q), ] %*% theta.Yl.vec +
#                                                 b.i[ which(dat.i$interval == q), ] %*% alpha.vec) )
#     }
#     
#   }
#   
#   # Log of full conditional of c(alpha, theta_{Y,l})
#   log.val <- sum( eta.iq * log(phi.iq) - lgamma(eta.iq) + lgamma(eta.tilde) -
#                     eta.tilde * log(phi.tilde) + log.exp.piece ) -
#     ( .5*length(theta.Yl.vec)*log(2*pi) - .5*log(det(Sig.Yl)) -
#         .5*t(theta.Yl.vec - mu.Yl) %*% solve(Sig.Yl) %*% (theta.Yl.vec - mu.Yl) ) -
#     ( .5*length(alpha.vec)*log(2*pi) - .5*log(det(Sig.a)) -
#         .5*t(alpha.vec - mu.a) %*% solve(Sig.a) %*% (alpha.vec - mu.a) )
#   
#   return(log.val)
#   
# }
# 
# 
# ### Negative of h(alpha, theta_{Y,l}) function
# # See function "h.alpha.theta.Yl" for function input descriptions
# neg.h.alpha.theta.Yl <- function( x, S, mu.a, Sig.a, mu.Yl, Sig.Yl, eta.iq, phi.iq, eta.tilde,
#                                   phi.tilde.summand.list, surv.dat, W.Yl, b.mat ){
#   val <- h.alpha.theta.Yl( x, S, mu.a, Sig.a, mu.Yl, Sig.Yl, eta.iq, phi.iq, eta.tilde,
#                            phi.tilde.summand.list, surv.dat, W.Yl, b.mat )
#   return( -1 * val )
# }


### Function for h(b_ij), where h(b_ij) = log( p(b_ij|xi, b_(-ij), D, M_l) )
###   i.e., the log of the full conditional distribution of alpha
# x: place holder for 1 (or 2) x 1 vector for b_ij
# reg.i: region number (1 - S) for subject of interest
# subj.j: subject number within region (1 - n.i) for subject of interest
# G.inv: 1 x 1 (or 2 x 2) inverted prior covariance matrix for b_ij
# eta.i: (1 x S) vector of shape hyperparameters for gamma prior on piecewise constant bh's for region i
# phi.i: (1 x S) vector of rate hyperparameters for gamma prior on piecewise constant bh's for region i
# eta.tilde.i: (1 x Q) matrix with values of eta.tilde.iq for region i
# phi.tilde.summand.i: (n.i x Q) matrix of summands for phi.tilde.iq for region i
# obs.times.i: observed times (events or censoring) for all subjects in region i
# surv.dat.ij: row in survival dataset corresponding to subject j from region i
# W.Yli: (n.i x (D_{Y,l} + p_Y)) design matrix for region i for survival submodel
# X.ij: K_ij x 1 vector of observed longitudinal responses for subject j from region i
# W.Xlij: design matrix for subject j in region i for longitudinal submodel
# z.t.mat.ij: K_ij x 1 (or 2) matrix with first column of 1s and (optionally) second column of
#      observed times X.ij (if random slope is included)
# theta.Yl: (D_{Y,l} + p_Y) x 1 vector of regression effects for survival submodel
# theta.Xl: vector of regression effects for longitudinal submodel
# b.i: n.i x 1 (or 2) matrix of subj-specific random effects for region i
# alpha.vec: 1 (or 2) x 1 vector of association parameters corresponding to shared random effects
# tau.val: scalar precision of error term in longitudinal likelihood
h.b.ij <- function( x, reg.i, subj.j, G.inv, eta.i, phi.i, eta.tilde.i, phi.tilde.summand.i,
                    obs.times.i, surv.dat.ij, W.Yli, X.ij, W.Xlij, z.t.mat.ij, theta.Yl, theta.Xl,
                    b.i, alpha.vec, tau.val ){
  
  # Verify that x is saved as a 1 x 1 (or 1 x 2)
  x <- rbind(x)
  eta.i <- rbind(eta.i)
  phi.i <- rbind(phi.i)
  theta.Yl <- cbind(theta.Yl)
  theta.Xl <- cbind(theta.Xl)
  alpha.vec <- cbind(alpha.vec)
  z.t.mat.ij <- cbind(z.t.mat.ij)
  b.i <- cbind(b.i)
  G.inv <- as.matrix(G.inv)
  
  # Calculate phi.tilde.iq for the given i of interest and for all q
  b.i[subj.j,] <- x              # update n.i^th row of matrix of random effects for reg. i
  Q <- ncol(phi.iq)
  exp.lp.i <- as.numeric(exp( b.i %*% alpha.vec + W.Yli %*% theta.Yl ))   # exp. of linear predictor
  phi.tilde.i <- phi.iq[reg.i,] + colSums( phi.tilde.summand.i * exp.lp.i )    # 1 x q matrix
  
  # Calculate log of exp( delta_ijq * nu_ij * b_ij' alpha )
  log.exp.piece <- surv.dat.ij$status * as.numeric(x %*% alpha.vec)
  
  # Log of longitudinal likelihood for subject j from region i
  lklhd.mean <- W.Xlij %*% theta.Xl + z.t.mat.ij %*% matrix(x, ncol = 1)
  log.long.lklhd.ij <- -.5*length(X.ij)*log(2*pi) + .5*length(X.ij)*log(tau.val) -
    .5 * tau.val * t(X.ij - lklhd.mean) %*% (X.ij - lklhd.mean)
  
  # Log of full conditional of b_ij
  log.val <- sum( eta.i * log(phi.i) - lgamma(eta.i) + lgamma(eta.tilde.i) -
                    eta.tilde.i * log(phi.tilde.i) ) + log.exp.piece +
    log.long.lklhd.ij - .5*length(x)*log(2*pi) + .5*log(det(G.inv)) - .5 * x %*% G.inv %*% t(x)
  
  return(log.val)
  
}
# h.b <- function( x, G.inv, eta.iq, phi.iq, eta.tilde, phi.tilde.summand.list, surv.dat, W.Yl, X.ijk,
#                  W.Xl, long.dat.subjid, z.t.mat, theta.Yl, theta.Xl, alpha.vec, tau.val ){
#   
#   # Verify that x is saved as an N x 1 (or N x 2) matrix
#   x <- cbind(x)
#   theta.Yl <- cbind(theta.Yl)
#   theta.Xl <- cbind(theta.Xl)
#   X.ijk <- cbind(X.ijk)
#   alpha.vec <- cbind(alpha.vec)
#   z.t.mat <- cbind(z.t.mat)
#   G.inv <- as.matrix(G.inv)
#   
#   # Calculate phi.tilde.iq for all i and q
#   # Also calculate log of exp( sum_j^n.i delta_ijq * nu_ij * b_ij' alpha ) for all i
#   Q <- ncol(phi.iq)
#   phi.tilde <- matrix( 0, nrow = nrow(phi.iq), ncol = ncol(phi.iq) )
#   log.exp.piece <- matrix( 0, nrow = nrow(phi.iq), ncol = ncol(phi.iq) )
#   for(i in 1:S){
#     
#     # phi.tilde.iq
#     reg.i <- sort( unique(surv.dat$region) )[i]
#     x.i <- cbind( x[ which(surv.dat$region == reg.i), ] )   # random effects for reg. i
#     W.Yli <- as.matrix( W.Yl[ which(surv.dat$region == reg.i), ] )  # surv. design matrix for reg. i
#     exp.lp.i <- as.numeric(exp( x.i %*% alpha.vec + W.Yli %*% theta.Yl ))   # exp. of linear predictor
#     phi.tilde[i,] <- phi.iq[i,] + colSums( phi.tilde.summand.list[[i]] * exp.lp.i )
#     
#     # log of exp( sum_j^n.i delta_ijq * nu_ij * alpha' b_ij )
#     dat.i <- surv.dat[ which(surv.dat$region == reg.i), ]
#     for(q in 1:Q){
#       log.exp.piece[i,q] <- sum( dat.i$status[ which(dat.i$interval == q) ] *
#                                    as.numeric(x[ which(dat.i$interval == q), ] %*% alpha.vec) )
#     }
#     
#   }
#   
#   # Log of longitudinal likelihood for all subject
#   if(ncol(x) == 1){
#     x.long.mat.new <- cbind( rep(x[,1], table(long.dat.subjid)) )
#   } else if(ncol(x) == 2){
#     x.long.mat.new <- cbind( rep(x[,1], table(long.dat.subjid)),
#                              rep(x[,2], table(long.dat.subjid)) )
#   }
#   lklhd.mean <- W.Xl %*% theta.Xl + rowSums(x.long.mat.new * z.t.mat)
#   log.long.lklhd.ij <- -.5*length(X.ijk)*log(2*pi) + .5*length(X.ijk)*log(tau.val) -
#     .5 * tau.val * t(X.ijk - lklhd.mean) %*% (X.ijk - lklhd.mean)
#   
#   # Log of prior distribution for all subjects
#   if(ncol(x) == 1){
#     log.prior <- -.5*length(x)*log(2*pi) + .5*log(det(G.inv)) - .5 * t(x) %*% x %*% G.inv
#   } else if(ncol(x) > 1){
#     log.prior <- sum( apply(x, 1, function(x)
#       -.5*length(x)*log(2*pi) + .5*log(det(G.inv)) - .5 * rbind(x) %*% G.inv %*% cbind(x)) )
#   }
#   
#   # Log of full conditional of b_ij
#   log.val <- sum( eta.iq * log(phi.iq) - lgamma(eta.iq) + lgamma(eta.tilde) -
#                     eta.tilde * log(phi.tilde) + log.exp.piece ) +
#     log.long.lklhd.ij + log.prior
#   
#   return(log.val)
#   
# }


### Negative of h(b.ij) function
# See function "h.b.ij" for function input descriptions
neg.h.b.ij <- function( x, reg.i, subj.j, G.inv, eta.i, phi.i, eta.tilde.i, phi.tilde.summand.i,
                        obs.times.i, surv.dat.ij, W.Yli, X.ij, W.Xlij, z.t.mat.ij, theta.Yl, theta.Xl,
                        b.i, alpha.vec, tau.val ){
  val <- h.b.ij( x, reg.i, subj.j, G.inv, eta.i, phi.i, eta.tilde.i, phi.tilde.summand.i,
                 obs.times.i, surv.dat.ij, W.Yli, X.ij, W.Xlij, z.t.mat.ij, theta.Yl, theta.Xl,
                 b.i, alpha.vec, tau.val )
  return( -1 * val )
}
# neg.h.b <- function( x, G.inv, eta.iq, phi.iq, eta.tilde, phi.tilde.summand.list, surv.dat, W.Yl, X.ijk,
#                      W.Xl, long.dat.subjid, z.t.mat, theta.Yl, theta.Xl, alpha.vec, tau.val ){
#   val <- h.b( x, G.inv, eta.iq, phi.iq, eta.tilde, phi.tilde.summand.list, surv.dat, W.Yl, X.ijk,
#               W.Xl, long.dat.subjid, z.t.mat, theta.Yl, theta.Xl, alpha.vec, tau.val )
#   return( -1 * val )
# }


### Hessian of h(theta_{Y,l}) function
# See function "h.theta.Yl" for function input descriptions
h.pp.theta.Yl <- function( x, S, Sig.Yl, phi.iq, eta.tilde, phi.tilde.summand.list,
                           surv.dat, W.Yl, b.mat, alpha.vec ){
  
  # Verify that x is saved as a (D_{Y,l} + p_Y) x 1 matrix
  x <- cbind(x)
  alpha.vec <- cbind(alpha.vec)
  W.Yl <- cbind(W.Yl)
  Sig.Yl <- as.matrix(Sig.Yl)
  
  # First piece
  p1 <- -solve(Sig.Yl)
  
  # Second piece
  p2 <- matrix(0, nrow = ncol(W.Yl), ncol = ncol(W.Yl))
  Q <- ncol(phi.iq)
  for(i in 1:S){
    
    # Calculate the exponent of the linear predictor
    reg.i <- sort( unique(surv.dat$region) )[i]
    b.i <- as.matrix( b.mat[ which(surv.dat$region == reg.i), ] )    # random effects for reg. i
    W.Yli <- as.matrix( W.Yl[ which(surv.dat$region == reg.i), ] )
    exp.lp <- as.numeric(exp(b.i %*% alpha.vec + W.Yli %*% x))
    
    # Calculate second piece of Hessian matrix
    for(q in 1:Q){
      c.hat.iq <- phi.tilde.summand.list[[i]][,q]
      p2 <- p2 + eta.tilde[i,q] *
        ( (t( c.hat.iq * exp.lp * W.Yli ) %*% W.Yli) *
            (phi.iq[i,q] + sum(exp.lp * c.hat.iq) )^-1 -
            cbind( colSums( c.hat.iq * exp.lp * W.Yli ) ) %*%
            rbind( colSums( c.hat.iq * exp.lp * W.Yli ) ) *
            (phi.iq[i,q] + sum( exp.lp * c.hat.iq ) )^-2 )
    }
  }
  
  # Hessian of h(theta_{Y,l}) function
  hess <- p1 - p2
  return(hess)
  
}


### Hessian of h(alpha) function
# See function "h.alpha" for function input descriptions
h.pp.alpha <- function( x, S, Sig.a, phi.iq, eta.tilde, phi.tilde.summand.list,
                        surv.dat, W.Yl, theta.Yl, b.mat ){
  
  # Verify that x is saved as a 1 x 1 (or 2 x 1) matrix
  x <- cbind(x)
  theta.Yl <- cbind(theta.Yl)
  Sig.a <- as.matrix(Sig.a)
  
  # First piece
  p1 <- -solve(Sig.a)
  
  # Second piece
  p2 <- matrix(0, nrow = ncol(b.mat), ncol = ncol(b.mat))
  Q <- ncol(phi.iq)
  for(i in 1:S){
    
    # Calculate the exponent of the linear predictor
    reg.i <- sort( unique(surv.dat$region) )[i]
    b.i <- as.matrix( b.mat[ which(surv.dat$region == reg.i), ] )    # random effects for reg. i
    W.Yli <- as.matrix( W.Yl[ which(surv.dat$region == reg.i), ] )
    exp.lp <- as.numeric(exp(b.i %*% x + W.Yli %*% theta.Yl))
    
    # Calculate second piece of Hessian matrix
    for(q in 1:Q){
      c.hat.iq <- phi.tilde.summand.list[[i]][,q]
      p2 <- p2 + eta.tilde[i,q] *
        ( (t( c.hat.iq * exp.lp * b.i ) %*% b.i) *
            (phi.iq[i,q] + sum(exp.lp * c.hat.iq) )^-1 -
            cbind( colSums( c.hat.iq * exp.lp * b.i ) ) %*%
            rbind( colSums( c.hat.iq * exp.lp * b.i ) ) *
            (phi.iq[i,q] + sum( exp.lp * c.hat.iq ) )^-2 )
    }
  }
  
  # Hessian of h(alpha) function
  hess <- p1 - p2
  return(hess)
  
}


### Function for piece of Hessian of h(theta_{Y,l}) (or h(alpha)) after taking derivative with
### respect to both theta_{Y,l} and alpha
# See function "h.theta.Yl" or "h.alpha" for function input descriptions
h.pp.theta.Yl.alpha <- function( S, phi.iq, eta.tilde, phi.tilde.summand.list,
                                 surv.dat, W.Yl, theta.Yl, alpha.vec, b.mat ){
  
  # Verify that inputs are saved as matrices
  alpha.vec <- cbind(alpha.vec)
  theta.Yl <- cbind(theta.Yl)
  
  # Piece of Hessian matrix
  hess <- matrix(0, nrow = length(theta.Yl), ncol = length(alpha.vec))
  Q <- ncol(phi.iq)
  for(i in 1:S){
    
    # Calculate the exponent of the linear predictor
    reg.i <- sort( unique(surv.dat$region) )[i]
    b.i <- as.matrix( b.mat[ which(surv.dat$region == reg.i), ] )    # random effects for reg. i
    W.Yli <- as.matrix( W.Yl[ which(surv.dat$region == reg.i), ] )
    exp.lp <- as.numeric(exp(b.i %*% alpha.vec + W.Yli %*% theta.Yl))
    
    # Calculate piece of Hessian matrix
    for(q in 1:Q){
      c.hat.iq <- phi.tilde.summand.list[[i]][,q]
      hess <- hess - eta.tilde[i,q] *
        ( (t( c.hat.iq * exp.lp * W.Yli ) %*% b.i) *
            (phi.iq[i,q] + sum(exp.lp * c.hat.iq) )^-1 -
            cbind( colSums( c.hat.iq * exp.lp * W.Yli ) ) %*%
            rbind( colSums( c.hat.iq * exp.lp * b.i ) ) *
            (phi.iq[i,q] + sum( exp.lp * c.hat.iq ) )^-2 )
    }
  }
  return(hess)
  
}


### Hessian of h(G^-1) function - either 1 x 1 or 3 x 3 depending on if only random intercepts
### or both random intercepts and slopes are used (depends on number of unique parameters in G^-1)
# G.inv: 1 x 1 (or 2 x 2) inverted prior covariance matrix for b_ij
# nu.G: degrees of freedom for Wishart prior on inverse covariance matrix of random effects
# N: total number of subjects
h.pp.G.inv <- function( G.inv, nu.G, N ){
  
  # Verify that G.inv is saved as a matrix
  G.inv <- as.matrix(G.inv)
  
  # Calculate (nu0 - r + N - 1)
  nu.G.piece <- nu.G - nrow(G.inv) + N - 1
  
  # Identify if only random intercepts or both random intercepts and slopes are inlcuded
  if(nrow(G.inv) == 1){
    
    # Extract elements of G^-1
    g.inv.11 <- as.numeric(G.inv)
    
    # Hessian piece for G^-1 = g.inv.11
    g.inv.11 <- as.numeric(G.inv)
    hess <- as.matrix( -nu.G.piece / (2 * g.inv.11^2) )
    
  } else if(nrow(G.inv) == 2){
    
    # Extract elements of G^-1
    g.inv.11 <- as.numeric(G.inv[1,1])
    g.inv.12 <- as.numeric(G.inv[1,2])
    g.inv.22 <- as.numeric(G.inv[2,2])
    
    # Calculate various pieces of Hessian matrix
    g.inv.11.h <- -g.inv.22^2 * nu.G.piece / (2 * (g.inv.11 * g.inv.22 - 2 * g.inv.12)^2)
    g.inv.22.h <- -g.inv.11^2 * nu.G.piece / (2 * (g.inv.11 * g.inv.22 - 2 * g.inv.12)^2)
    g.inv.12.h <- -2 * nu.G.piece / (g.inv.11 * g.inv.22 - 2 * g.inv.12)^2
    g.inv.11.22.h <- -g.inv.12 * nu.G.piece / (g.inv.11 * g.inv.22 - 2 * g.inv.12)^2
    g.inv.11.12.h <- g.inv.22 * nu.G.piece / (g.inv.11 * g.inv.22 - 2 * g.inv.12)^2
    g.inv.22.12.h <- g.inv.11 * nu.G.piece / (g.inv.11 * g.inv.22 - 2 * g.inv.12)^2
    
    # Compile Hessian matrix
    hess <- matrix(0, nrow = 3, ncol = 3)
    hess[1,1] <- g.inv.11.h
    hess[2,2] <- g.inv.22.h
    hess[3,3] <- g.inv.12.h
    hess[1,2] <- hess[2,1] <- g.inv.11.22.h
    hess[1,3] <- hess[3,1] <- g.inv.11.12.h
    hess[2,3] <- hess[3,2] <- g.inv.22.12.h
    
  }
  
  return(hess)
  
}


### Hessian of h(b_ij) function
# See function "h.b.ij" for function input descriptions
# h.pp.b.ij <- function( x, subj.j, G.inv, phi.i, eta.tilde.i, phi.tilde.summand.i, W.Yli,
#                        z.t.mat.ij, b.i, theta.Yl, alpha.vec, tau.val ){
#   
#   # Verify that x is saved as a 1 x 1 (or 1 x 2)
#   x <- rbind(x)
#   phi.i <- rbind(phi.i)
#   theta.Yl <- cbind(theta.Yl)
#   alpha.vec <- cbind(alpha.vec)
#   z.t.mat.ij <- cbind(z.t.mat.ij)
#   b.i <- cbind(b.i)
#   G.inv <- as.matrix(G.inv)
#   W.Yli <- as.matrix(W.Yli)
#   
#   # First piece
#   p1 <- -G.inv
#   
#   # Second piece
#   p2 <- tau.val * t(z.t.mat.ij) %*% z.t.mat.ij
#   
#   # Third piece
#   p3 <- 0
#   b.i[subj.j,] <- x       # replace row j of region i
#   exp.lp <- as.numeric(exp(b.i %*% alpha.vec + W.Yli %*% theta.Yl))  # exponent of the linear predictor
#   Q <- ncol(phi.iq)
#   for(q in 1:Q){
#     c.hat.iq <- phi.tilde.summand.i[,q]
#     c.hat.g.hat.sum <- sum(c.hat.iq * exp.lp)
#     p3 <- p3 + eta.tilde[i,q] *
#       ( c.hat.g.hat.sum * (phi.iq[i,q] + c.hat.g.hat.sum)^-1 -
#           c.hat.g.hat.sum^2 * (phi.iq[i,q] + c.hat.g.hat.sum)^-2 )
#   }
#   p3 <- p3 * alpha.vec %*% t(alpha.vec)
#   
#   # Hessian of h(b_ij) function
#   hess <- p1 - p2 - p3
#   return(hess)
#   
# }


### Function for approximated log marginal likelihood for model M_l, l = 1,...,L
# surv.dat: full survival dataset for all subjects from all regions
# long.dat: full longitudinal dataset for all subjects from all regions
# W.Yl: (N x (D_{Y,l} + p_Y)) design matrix for survival submodel
# W.Xl: design matrix for longitudinal submodel
# eta.iq: (S x Q) matrix of shape hyperparameters for gamma prior on piecewise constant baseline hazards
# phi.iq: (S x Q) matrix of rate hyperparameters for gamma prior on piecewise constant baseline hazards
# eta.tilde: (S x Q) matrix with values of eta.tilde.iq
# phi.tilde.summand.list: list with all S (n.i x Q) matrices of summands for phi.tilde.iq
# theta.Yl: (D_{Y,l} + p_Y) x 1 vector of regression effects for survival submodel
# theta.Xl: vector of regression effects for longitudinal submodel
# alpha.vec: 1 (or 2) x 1 vector of association parameters corresponding to shared random effects
# tau.val: scalar precision of error term in longitudinal likelihood
# b.mat: N x 1 (or N x 2) matrix of subj-specific REs (rand. intercept and optionally rand. slope)
# G.inv: 1 x 1 (or 2 x 2) inverted prior covariance matrix for b_ij
# mu.Yl: prior mean vector for theta_{Y,l}
# Sig.Yl: prior covariance matrix for theta_{Y,l}
# mu.Xl: prior mean vector for theta_{X,l}
# Sig.Xl: prior covariance matrix for theta_{X,l}
# mu.alpha: prior mean vector (either 1 x 1 or 2 x 1) for association parameter(s)
# Sig.alpha: prior covariance matrix (either 1 x 1 or 2 x 2) for association parameter(s)
# eta.tau: scalar shape hyperparameters for prior on precision of longitudinal likelihood error
# phi.tau: scalar rate hyperparameters for prior on precision of longitudinal likelihood error
# nu.G: degrees of freedom for Wishart prior on inverse covariance matrix of random effects
# C0.inv.G: inverse scale matrix for Wishart prior on inverse covariance matrix of random effects
approx.lml <- function( surv.dat, long.dat, W.Yl, W.Xl, eta.iq, phi.iq, eta.tilde,
                        phi.tilde.summand.list, theta.Yl, theta.Xl, alpha.vec, tau.val,
                        b.mat, G.inv, mu.Yl, Sig.Yl, mu.Xl, Sig.Xl, mu.alpha, Sig.alpha,
                        eta.tau, phi.tau, nu.G, C0.inv.G ){
  
  # Verify that input vectors/matrices are saved as matrices
  W.Yl <- cbind(W.Yl)
  W.Xl <- cbind(W.Xl)
  theta.Yl <- cbind(theta.Yl)
  theta.Xl <- cbind(theta.Xl)
  alpha.vec <- cbind(alpha.vec)
  mu.Yl <- as.matrix(mu.Yl)
  Sig.Yl <- as.matrix(Sig.Yl)
  mu.Xl <- as.matrix(mu.Xl)
  Sig.Xl <- as.matrix(Sig.Xl)
  mu.alpha <- as.matrix(mu.alpha)
  Sig.alpha <- as.matrix(Sig.alpha)
  C0.inv.G <- as.matrix(C0.inv.G)
  
  # Calculate phi.tilde.iq for all i and q
  # Also calculate log of exp( sum_j^n.i delta_ijq * nu_ij * w_{Y,lij}' theta_{Y,l} ) for all i
  S <- length(unique(surv.dat$region))
  Q <- ncol(phi.iq)
  phi.tilde <- matrix( 0, nrow = nrow(phi.iq), ncol = ncol(phi.iq) )
  log.exp.piece <- matrix( 0, nrow = nrow(phi.iq), ncol = ncol(phi.iq) )
  for(i in 1:S){
    
    # phi.tilde.iq
    reg.i <- sort( unique(surv.dat$region) )[i]
    b.i <- as.matrix( b.mat[ which(surv.dat$region == reg.i), ] )   # random effects for reg. i
    W.Yli <- as.matrix( W.Yl[ which(surv.dat$region == reg.i), ] )  # surv. design matrix for reg. i
    exp.lp.i <- as.numeric(exp( b.i %*% alpha.vec + W.Yli %*% theta.Yl ))  # exp. of linear predictor
    phi.tilde[i,] <- phi.iq[i,] + colSums( phi.tilde.summand.list[[i]] * exp.lp.i )
    
    # log of exp( sum_j^n.i delta_ijq * nu_ij * [alpha' b_ij + w_{Y,lij}' theta_{Y,l}] )
    dat.i <- surv.dat[ which(surv.dat$region == reg.i), ]
    for(q in 1:Q){
      log.exp.piece[i,q] <- sum( dat.i$status[ which(dat.i$interval == q) ] *
                                   (as.numeric(b.i[ which(dat.i$interval == q), ] %*% alpha.vec) +
                                      as.numeric(W.Yli[ which(dat.i$interval == q), ] %*% theta.Yl)) )
    }
    
  }
  
  # Log of survival likelihood for all subject after integrating out baseline hazards
  # (maximized at mode of parameters and REs)
  surv.lklhd <- sum( eta.iq * log(phi.iq) - lgamma(eta.iq) + lgamma(eta.tilde) -
                       eta.tilde * log(phi.tilde) + log.exp.piece )
  
  # Log of longitudinal likelihood for all subject (maximized at mode of parameters and REs)
  if(ncol(b.mat) == 1){
    z.t.mat <- cbind( rep(1, nrow(long.dat)) )
    b.long.mat.new <- cbind( rep(b.mat[,1], table(long.dat$subjid)) )
  } else if(ncol(b.mat) == 2){
    z.t.mat <- cbind( rep(1, nrow(long.dat)), long.dat$X.ijk )
    b.long.mat.new <- cbind( rep(b.mat[,1], table(long.dat$subjid)),
                             rep(b.mat[,2], table(long.dat$subjid)) )
  }
  long.mean <- W.Xl %*% theta.Xl + rowSums(b.long.mat.new * z.t.mat)
  long.lklhd <- -.5*length(long.dat$X.ijk)*log(2*pi) + .5*length(long.dat$X.ijk)*log(tau.val) -
    .5 * tau.val * t(long.dat$X.ijk - long.mean) %*% (long.dat$X.ijk - long.mean)
  
  # Maximized log-likelihood (after baseline hazards integrated out)
  piece1 <- surv.lklhd + long.lklhd
  
  # Term involving number of unknown fixed effects
  num.param <- sum( c(length(theta.Yl), length(theta.Xl), length(alpha.vec), length(tau.val),
                      nrow(G.inv) + nrow(G.inv) * (nrow(G.inv) - 1)/2) )
  piece2 <- num.param/2 * log(2*pi)
  
  # Identify element locations for different fixed effects in Hessian matrix
  dims.alpha <- 1:length(alpha.vec)                                               # alpha
  dims.theta.Yl <- (max(dims.alpha)+1):(max(dims.alpha)+length(theta.Yl))         # theta.Yl
  dims.theta.Xl <- (max(dims.theta.Yl)+1):(max(dims.theta.Yl)+length(theta.Xl))   # theta.Xl
  dims.tau <- max(dims.theta.Xl) + 1                                              # tau
  dims.G.inv <- (dims.tau+1):(dims.tau+nrow(G.inv)+nrow(G.inv)*(nrow(G.inv)-1)/2) # elements of G.inv
  
  # Calculate pieces of Hessian matrix
  alpha.h <- h.pp.alpha( x = alpha.vec, S = S, Sig.a = Sig.alpha, phi.iq = phi.iq,
                         eta.tilde = eta.tilde, phi.tilde.summand.list = phi.tilde.summand.list,
                         surv.dat = surv.dat, W.Yl = W.Yl, theta.Yl = theta.Yl, b.mat = b.mat )
  theta.Yl.h <- h.pp.theta.Yl( x = theta.Yl, S = S, Sig.Yl = Sig.Yl, phi.iq = phi.iq,
                               eta.tilde = eta.tilde, phi.tilde.summand.list = phi.tilde.summand.list,
                               surv.dat = surv.dat, W.Yl = W.Yl, b.mat = b.mat, alpha.vec = alpha.vec )
  alpha.theta.Yl.h <- h.pp.theta.Yl.alpha( S = S, phi.iq = phi.iq, eta.tilde = eta.tilde,
                                           phi.tilde.summand.list = phi.tilde.summand.list,
                                           surv.dat = surv.dat, W.Yl = W.Yl, theta.Yl = theta.Yl,
                                           alpha.vec = alpha.vec, b.mat = b.mat )
  theta.Xl.h <- -as.numeric(tau.val) * ( t(W.Xl) %*% W.Xl + solve(Sig.Xl) )
  tau.h <- -as.numeric(tau.val)^-2 * (.5 * nrow(long.dat) + .5 * eta.tau - 1 )
  G.inv.h <- h.pp.G.inv( G.inv = G.inv, nu.G = nu.G, N = nrow(surv.dat) )
  
  # Compile Hessian matrix evaluated at maximized values
  hess.mat <- matrix(0, nrow = num.param, ncol = num.param)
  hess.mat[dims.alpha, dims.alpha] <- alpha.h
  hess.mat[dims.theta.Yl, dims.theta.Yl] <- theta.Yl.h
  hess.mat[dims.theta.Yl, dims.alpha] <- alpha.theta.Yl.h
  hess.mat[dims.alpha, dims.theta.Yl] <- t(alpha.theta.Yl.h)
  hess.mat[dims.theta.Xl, dims.theta.Xl] <- theta.Xl.h
  hess.mat[dims.tau, dims.tau] <- tau.h
  hess.mat[dims.G.inv, dims.G.inv] <- G.inv.h
  
  # Negative inverse Hessian of fixed effects evaluated at maximized values
  piece3 <- -.5 * log(det(-hess.mat))
  
  # Maximized prior distributions (log scale)
  theta.Yl.pr <- -.5*length(theta.Yl)*log(2*pi) - .5*log(det(Sig.Yl)) -
    .5*t(theta.Yl - mu.Yl) %*% solve(Sig.Yl) %*% (theta.Yl - mu.Yl)
  theta.Xl.pr <- -.5*length(theta.Xl)*log(2*pi) - .5*log(det(Sig.Xl)) -
    .5*t(theta.Xl - mu.Xl) %*% solve(Sig.Xl) %*% (theta.Xl - mu.Xl)
  alpha.pr <- -.5*length(alpha.vec)*log(2*pi) - .5*log(det(Sig.alpha)) -
    .5*t(alpha.vec - mu.alpha) %*% solve(Sig.alpha) %*% (alpha.vec - mu.alpha)
  if(ncol(b.mat) == 1){
    b.pr <- -.5*length(b.mat)*log(2*pi) + .5*log(det(G.inv)) - .5 * t(b.mat) %*% b.mat %*% G.inv
  } else if(ncol(b.mat) > 1){
    b.pr <- sum( apply(x, 1, function(x)
      -.5*length(x)*log(2*pi) + .5*log(det(G.inv)) - .5 * rbind(x) %*% G.inv %*% cbind(x)) )
  }
  tau.pr <- dgamma(tau.val, shape = eta.tau, rate = phi.tau, log = TRUE)
  G.inv.pr <- .5 * (nu.G - nrow(G.inv) - 1) * log(det(G.inv)) - .5 * sum(diag(C0.inv.G %*% G.inv)) -
    .5 * nu.G * nrow(G.inv) * log(2) + .5 * nu.G * log(det(C0.inv.G)) -
    .25 * nrow(G.inv) * (nrow(G.inv) - 1) * log(pi) - sum( log(gamma(.5*(nu.G + 1 - 1:nrow(G.inv)))) )
  piece4 <- theta.Yl.pr + theta.Xl.pr + alpha.pr + b.pr + tau.pr + G.inv.pr
  
  
  # Approximated log marginal likelihood
  approx.log.ml <- piece1 + piece2 + piece3 + piece4
  
  return(approx.log.ml)
  
}


### Function to create survival design matrix for l^th survival partition
# surv.dat: full dataset for all regions combined for survival submodel
#    Dataset should have the following columns (with same names, case sensitive):
#      - subjid: subject ID
#      - eventtime: observed time (event time or censored time)
#      - status: event indicator (1 = event, 0 = censored)
#      - trt: treatment group indicator (1 = treatment group, 0 = control group)
#      - region: numeric values for regions (1 to number of regions)
#      - OPTIONAL COVARIATES CAN BE INCLUDED; NAMES SHOULD BE SPECIFIED IN SEPARATE VECTOR
# surv.mod.l: S-dimensional vector corresponding to l^th partition of regions into sets
# cov.surv.names: vector of names of optional baseline covariates to include in surv. model
update.W.Y <- function( surv.dat, surv.mod.l, cov.surv.names = NULL ){
  
  # Create frame for design matrix for survival submodel
  D.Y.l <- length(unique(surv.mod.l))
  p.Y <- length(cov.surv.names)
  W.Yl.mat <- matrix( 0, nrow = nrow(surv.dat), ncol = (D.Y.l + p.Y) )
  
  # Add region-specific treatment effects
  trt.names.l <- NULL        # store names of treatment effects
  for(i in 1:D.Y.l){
    which.regs.i <- which(surv.mod.l == i)       # identify regions in i^th set
    W.Yl.mat[,i] <- ifelse(surv.dat$trt == 1 & surv.dat$region %in% which.regs.i, 1, 0)
    trt.names.l <- c( trt.names.l, paste0("trt.", paste(which.regs.i, collapse = "")) )
  }
  
  # Add covariates (if any)
  if(p.Y > 0){
    for(p in 1:p.Y){
      W.Yl.mat[,S+p] <- surv.dat[, which(colnames(surv.dat) == cov.surv.names[p]) ]
    }
  }
  
  # Add column names to design matrix
  #    Columns 1:D.Y.l - reg-specific trt effects
  #    Columns (D.Y.l+1):(D.Y.l+p.X) - optional covariates
  colnames(W.Yl.mat) <- c(trt.names.l, cov.surv.names)
  
  return(W.Yl.mat)
  
}


### Function to create longitudinal design matrix for l^th longitudinal partition
# long.dat: full dataset for all regions combined for longitudinal submodel
#    Dataset should have the following columns (with same names, case sensitive):
#      - subjid: subject ID
#      - X.ijk: observed responses of longitudinal marker
#      - time: time points when longitudinal responses were observed
#      - trt: treatment group indicator (1 = treatment group, O = control group)
#      - region: numeric values for regions (1 to number of regions, sorted numerically)
#      - OPTIONAL COVARIATES CAN BE INCLUDED; NAMES SHOULD BE SPECIFIED IN SEPARATE VECTOR
# long.mod.l: S-dimensional vector corresponding to l^th partition of regions into sets
# lin.spline.knots: vector of time points of knots used for optional linear splines
# cov.long.names: vector of names of optional baseline covariates to include in long. model
update.W.X <- function( long.dat, long.mod.l, lin.spline.knots = NULL, cov.long.names = NULL ){
  
  ### Construct design matrices for longitudinal submodel
  S <- length(long.mod.l)
  D.X.l <- length(unique(long.mod.l))
  num.knots <- length(lin.spline.knots)
  p.X <- length(cov.long.names)
  W.Xl.mat <- matrix( 0, nrow = nrow(long.dat),
                      ncol = ((2 + num.knots)*D.X.l + S + num.knots + p.X + 1) )
  
  # Add region-specific intercepts
  int.names <- paste0("int.", 1:S)
  for(i in 1:S){
    W.Xl.mat[,i] <- ifelse(long.dat$region == i, 1, 0)
  }
  
  # Add region-specific treatment effects and region-specific intercepts
  trt.names.l <- NULL        # store names of treatment effects
  for(i in 1:D.X.l){
    which.regs.i <- which(long.mod.l == i)       # identify regions in i^th set
    W.Xl.mat[,S+1+num.knots+i] <- ifelse(long.dat$trt == 1 &
                                           long.dat$region %in% which.regs.i, 1, 0)
    trt.names.l <- c( trt.names.l, paste0("trt.", paste(which.regs.i, collapse = "")) )
  }
  
  # Add time column and D.X.l columns for time-by-trt interactions
  W.Xl.mat[,S+1] <- long.dat$time    # time
  time.trt.names.l <- NULL          # store names of time-by-trt interactions
  for(i in 1:D.X.l){
    which.regs.i <- which(long.mod.l == i)       # identify regions in i^th set
    W.Xl.mat[,(S+num.knots+1+D.X.l+i)] <-
      ifelse(long.dat$trt == 1 & long.dat$region %in% which.regs.i,
             long.dat$time, 0)                   # time-by-trt interaction
    time.trt.names.l <- c( time.trt.names.l,
                           paste0("time:trt.", paste(which.regs.i, collapse = "")) )
  }
  
  # Add columns for knots for linear splines (if any)
  if(num.knots > 0){
    
    # Names of knots (knot.1, knot.2, ... )
    knots.names <- paste0("knot.", 1:num.knots)
    knots.trt.names.l <- NULL    # store names of knots-by-trt interactions
    
    # Calculate (time - knot) if time > knot
    for(p in 1:num.knots){
      spline.p.val <- ifelse( long.dat$time - lin.spline.knots[p] > 0,
                              long.dat$time - lin.spline.knots[p], 0 )
      W.Xl.mat[,S+1+p] <- spline.p.val
      for(i in 1:D.X.l){
        which.regs.i <- which(long.mod.l == i)     # identify regions in i^th set
        W.Xl.mat[,((1+p)*D.X.l+S+num.knots+1+i)] <-
          ifelse(long.dat$trt == 1 & long.dat$region %in% which.regs.i,
                 spline.p.val, 0)                  # knot-by-trt interaction
        knots.trt.names.l <- c( knots.trt.names.l,
                                paste0(knots.names[p], ":trt.",
                                       paste(which.regs.i, collapse = "")) )
      }
    }
    
  } else{
    knots.names <- NULL
    knots.trt.names.l <- NULL
  }
  
  # Add covariates (if any)
  if(p.X > 0){
    for(p in 1:p.X){
      W.Xl.mat[,((2+num.knots)*D.X.l+S+num.knots+1+p)] <-
        long.dat[, which(colnames(long.dat) == cov.long.names[p]) ]
    }
  }
  
  # Add column names to longitudinal design matrix
  #    Columns 1:S - reg-specific intercepts
  #    Column (S+1) - time
  #    Columns (S+2):(S+num.knots+1) - optional spline knots
  #    Columns (S+num.knots+2):(S+num.knots+D.X.l+1) - reg-specific trt effects
  #    Columns (S+num.knots+D.X.l+2):(S+num.knots+2*D.X.l+1) - time-by-trt
  #    Columns (S+num.knots+2*D.X.l+2):(S+num.knots+(2+num.knots)*D.X.l+1) - knots-by-trt
  #    Columns (S+num.knots+(2+num.knots)*D.X.l+2):(S+num.knots+(2+num.knots)*D.X.l+1+p.X) -
  #         optional covariates
  colnames(W.Xl.mat) <- c(int.names, "time", knots.names, trt.names.l,
                          time.trt.names.l, knots.trt.names.l, cov.long.names)
  
  return(W.Xl.mat)
  
}


### Function to create prior mean vector and covariance matrix for l^th survival partition
# mu.Y: (S + p.Y) x 1 prior mean vector for survival treatment and covariate effects
# Sig.Y: (S + p.Y) x (S + p.Y) prior covariance matrix for survival treatment and covariate effects
# surv.mod.l: S-dimensional vector corresponding to l^th partition of regions into sets
# cov.surv.names: vector of names of optional baseline covariates to include in surv. model
update.mu.Sig.Y <- function( mu.Y, Sig.Y, surv.mod.l, cov.surv.names = NULL ){
  
  # Update mu.Y and Sig.Y for survival submodel l
  # (Assume all regions have same prior means and variances)
  S <- length(surv.mod.l)
  D.Y.l <- length(unique(surv.mod.l))
  if(is.null(cov.surv.names)){
    mu.Yl <- mu.Y[ 1:D.Y.l ]
    Sig.Yl <- Sig.Y[ 1:D.Y.l, 1:D.Y.l ]
  } else{
    p.Y <- length(cov.surv.names)
    mu.Yl <- mu.Y[ c(1:D.Y.l, (S+1):(S+p.Y)) ]
    Sig.Yl <- Sig.Y[ c(1:D.Y.l, (S+1):(S+p.Y)), c(1:D.Y.l, (S+1):(S+p.Y)) ]
  }
  
  # Return list with both mu.Yl and Sig.Yl
  mu.Sig.Yl <- list( mu.Yl = mu.Yl, Sig.Yl = Sig.Yl )
  return(mu.Sig.Yl)
  
}


### Function to create prior mean vector and covariance matrix for l^th longitudinal partition
# mu.X: prior mean vector for longitudinal treatment and covariate effects
# Sig.X: prior covariance matrix for longitudinal treatment and covariate effects
# long.mod.l: S-dimensional vector corresponding to l^th partition of regions into sets
# lin.spline.knots: vector of time points of knots used for optional linear splines
# cov.long.names: vector of names of optional baseline covariates to include in long. model
update.mu.Sig.X <- function( mu.X, Sig.X, long.mod.l, lin.spline.knots = NULL,
                             cov.long.names = NULL ){
  
  # Identify which indices to keep
  S <- length(long.mod.l)
  D.X.l <- length(unique(long.mod.l))
  num.knots <- length(lin.spline.knots)
  p.X <- length(cov.long.names)
  inx.keep <- 1:(S+1+num.knots)    # indices for intercepts, time, and knots (if any)
  inx.keep <- c(inx.keep,
                (S+2+num.knots):(S+1+num.knots+D.X.l),      # indices for trt effects
                (2*S+2+num.knots):(2*S+1+num.knots+D.X.l))  # indices for time-by-trt interactions
  if(num.knots > 0){
    for(p in 1:num.knots){
      inx.keep <- c(inx.keep, ((2+p)*S+2+num.knots):((2+p)*S+1+num.knots+D.X.l))
    }
  }
  if(p.X > 0){
    inx.keep <- c(inx.keep, ((3+num.knots)*S+2+num.knots):((3+num.knots)*S+1+num.knots+p.X))
  }
  
  # Update mu.X and Sig.X for survival submodel l
  # (Assume all regions have same prior means and variances)
  mu.Xl <- mu.X[ inx.keep ]
  Sig.Xl <- Sig.X[ inx.keep, inx.keep ]
  
  # Return list with both mu.Xl and Sig.Xl
  mu.Sig.Xl <- list( mu.Xl = mu.Xl, Sig.Xl = Sig.Xl )
  return(mu.Sig.Xl)
  
}


### Function to calculate log density of multivariate normal distribution
# x: observed vector
# mu: mean vector
# Sig: covariance matrix
# dmvn.log <- function( x, mu, Sig ){
#   
#   x <- as.matrix(x)
#   mu <- as.matrix(mu)
#   Sig <- as.matrix(Sig)
#   log.val <- -.5*length(x)*log(2*pi) - .5*log(det(Sig)) -
#     .5*t(x - mu) %*% solve(Sig) %*% (x - mu)
#   return(log.val)
#   
# }


### Function to run BMA algorithm to fit joint models for a single survival/longitudinal dataset
# Returns posterior summaries, approximated PMPs, and consistency measures
# surv.dat: full dataset for all regions combined for survival submodel
#    Dataset should have the following columns (with same names, case sensitive):
#      - subjid: subject ID
#      - eventtime: observed time (event time or censored time)
#      - status: event indicator (1 = event, 0 = censored)
#      - trt: treatment group indicator (1 = treatment group, 0 = control group)
#      - region: numeric values for regions (1 to number of regions)
#      - OPTIONAL COVARIATES CAN BE INCLUDED; NAMES SHOULD BE SPECIFIED IN SEPARATE VECTOR
# cov.surv.names: vector of names of optional baseline covariates to include in surv. model
# int.cuts: Q x 1 vector of lower interval boundaries for all Q intervals
# long.dat: full dataset for all regions combined for longitudinal submodel
#    Dataset should have the following columns (with same names, case sensitive):
#      - subjid: subject ID
#      - X.ijk: observed responses of longitudinal marker
#      - time: time points when longitudinal responses were observed
#      - trt: treatment group indicator (1 = treatment group, O = control group)
#      - region: numeric values for regions (1 to number of regions, sorted numerically)
#      - OPTIONAL COVARIATES CAN BE INCLUDED; NAMES SHOULD BE SPECIFIED IN SEPARATE VECTOR
# lin.spline.knots: vector of time points of knots used for optional linear splines
# cov.long.names: vector of names of optional baseline covariates to include in long. model
# mod.mat: L x S matrix of possible region groupings, each row corresponding to one of L models
# a0.surv: weight for number of distinct survival trt effects in prior model probs.
# a0.long: weight for number of distinct longitudinal trt effects in prior model probs.
# mod.priors: L x 1 vector of prior model probs. (alternative to specifying a0.surv and a0.long)
# mu.Y: (S + p.Y) x 1 prior mean vector for survival treatment and covariate effects
# Sig.Y: (S + p.Y) x (S + p.Y) prior covariance matrix for survival treatment and covariate effects
# mu.X: prior mean vector for longitudinal treatment and covariate effects
# Sig.X: prior covariance matrix for longitudinal treatment and covariate effects
# mu.alpha: prior mean vector (either 1 x 1 or 2 x 1) for association parameter(s)
# Sig.alpha: prior covariance matrix (either 1 x 1 or 2 x 2) for association parameter(s)
# eta.iq: S x Q matrix of shape hyperparameters for prior on baseline hazard
# phi.iq: S x Q matrix of rate hyperparameters for prior on baseline hazard
# eta.tau: scalar shape hyperparameters for prior on precision of longitudinal likelihood error
# phi.tau: scalar rate hyperparameters for prior on precision of longitudinal likelihood error
# nu.G: degrees of freedom for Wishart prior on inverse covariance matrix of random effects
# C0.inv.G: inverse scale matrix for Wishart prior on inverse covariance matrix of random effects
# REs: indicator of which type of subject-specific random effects to include in the models
#      - "int": random intercepts only
#      - "int-slope": random intercepts and random slopes
# gamma0.surv: value of interest in null hypothesis to be tested (surv. trt effect)
# gamma0.long: value of interest in null hypothesis to be tested (long. trt effect)
# max.iter: maximum number of iterations in optimization algorithm
# optim.toler.params: tolerance for when convergence of parameters has been satisfied in
#       optimization algorithm
# optim.toler.b: tolerance for when convergence of random effects has been satisfied in
#       optimization algorithm
# n.draws: number of posterior samples to draw when evaluating consistency
# epsilon.star: vector of possible minimal clinically important regional differences of trtmnt
#       effects - used for consistency measures
# beta.star: probability cutoff for global inconsistency for which two regions are considered
#      to be clinically different
# pi0: value for which we want to calculate probabilities with local consistency
bma.joint.model <- function( surv.dat, cov.surv.names = NULL, int.cuts, long.dat,
                             lin.spline.knots = NULL, cov.long.names = NULL, mod.mat,
                             a0.surv = NULL, a0.long = NULL, mod.priors = NULL,
                             mu.Y, Sig.Y, mu.X, Sig.X, mu.alpha, Sig.alpha, eta.iq, phi.iq,
                             eta.tau, phi.tau, nu.G, C0.inv.G, REs = "int", gamma0.surv,
                             gamma0.long, max.iter = 10, optim.toler.params = .005,
                             optim.toler.b = .05, n.draws, epsilon.star, beta.star, pi0 ){
  
  
  ### Derive additional values based on function inputs
  S <- length( unique(surv.dat$region) )     # number of regions
  N <- nrow(surv.dat)                        # total sample size
  n.i <- table(surv.dat$region)              # regional sample sizes
  N.obs <- nrow(long.dat)                    # total number of longitudinal observations
  num.obs.i <- table(long.dat$region)        # number of longitudinal observations per region
  Q <- length(int.cuts)                      # number of intervals
  num.knots <- length(lin.spline.knots)        # number of knots for linear splines in long. submodel
  p.Y <- length(cov.surv.names)              # number of covariates for survival submodel
  p.X <- length(cov.long.names)              # number of covariates for longitudinal submodel
  
  
  ### Identify interval in which each patient failed or was censored
  surv.dat$interval <- apply( cbind(surv.dat$eventtime), 1, function(x){ sum(x > int.cuts) } )
  
  
  ### Construct design matrices for survival and longitudinal submodels
  W.Y.mat <- update.W.Y(surv.dat = surv.dat, surv.mod.l = 1:S)
  W.X.mat <- update.W.X(long.dat = long.dat, long.mod.l = 1:S, lin.spline.knots = lin.spline.knots)
  
  
  ### Identify number of association parameters depending on if only random intercepts are included
  ### in the model or both random intercepts and slopes
  if(REs == "int"){
    num.alpha <- 1
  } else if(REs == "int-slope"){
    num.alpha <- 2
  }
  
  
  ### Calculate eta.tilde.iq for all i and q (value remains constant across models and sample parameters).
  ### Also calculate ind.sum.mat, where ind.sum.mat_iq = sum_(j=1)^n.i (d_ijq * nu_ij)
  ### (i.e., the number of patients who received treatment from region i who failed in interval q)
  eta.tilde <- matrix( 0, nrow = S, ncol = Q )
  for(i in 1:S){
    for(q in 1:Q){
      d.ijq.times.nu.ij <- ifelse( surv.dat$region == i & surv.dat$status == 1 &
                                     surv.dat$interval == q, 1, 0 )
      eta.tilde[i,q] <- sum( d.ijq.times.nu.ij ) + eta.iq[i,q]
    }
  }
  
  
  ### Calculate delta_ijq * (y_ij - m_{q-1}) + sum_(g=q+1)^Q( delta_ijg * (m_q - m_{q-1}) )
  ### for all subjects in region i for each time interval q. Save (n.i x Q) matrix for each
  ### region in a list. These summands are later used to construct phi.tilde.iq.
  phi.tilde.summand.list <- list()
  for(i in 1:S){
    
    surv.dat.i <- surv.dat[ which(surv.dat$region == i), ]    # data for only region i
    
    # Matrix to store summand for all subjects in region i for each time interval q, q=1,...,Q
    summand.mat.i <- matrix( 0, nrow = nrow(surv.dat.i), ncol = Q )
    for(q in 1:Q){
      summand.val <- ifelse( surv.dat.i$interval == q, surv.dat.i$eventtime - int.cuts[q], 0 )
      summand.val <- ifelse( surv.dat.i$interval > q, int.cuts[q+1] - int.cuts[q], summand.val )
      summand.mat.i[,q] <- summand.val
    }
    phi.tilde.summand.list[[i]] <- summand.mat.i
    
  }
  
  
  ### Calculate degrees of freedom of posterior distribution of G^-1
  nu.tilde.G <- nu.G + N
  
  
  ### Details about model space
  L.surv <- nrow(mod.mat)         # number of models in survival-only model space
  L.long <- nrow(mod.mat)         # number of models in longitudinal-only model space
  L.all <- L.surv * L.long        # number of models in joint model space
  mods.fit <- data.frame( mod.num = 1:L.all,                        # model number
                          long.mod = rep(1:L.long, each = L.long),  # labels for long. submodels
                          surv.mod = rep(1:L.surv, L.surv) )        # labels for surv. submodels
  mods.fit$D.surv.l <- apply(mod.mat[mods.fit[,2],], 1,
                             function(x) length(unique(x)))  # numb. distinct trt effects for surv. submodels
  mods.fit$D.long.l <- apply(mod.mat[mods.fit[,3],], 1,
                             function(x) length(unique(x)))  # numb. distinct trt effects for long. submodels
  mods.fit$log.mar.lklhd <- 0     # place holder for log marginal likelihood
  
  
  ### Matrices/vectors to hold model-specific posterior summary statistics
  
  # Region-specific treatment effects (survival)
  rte.surv.mean.l.mat <- matrix( 0, nrow = L.all, ncol = S )
  rte.surv.sd.l.mat <- matrix( 0, nrow = L.all, ncol = S )
  rte.surv.less.gamma0.prob.l.mat <- matrix( 0, nrow = L.all, ncol = S )
  
  # Global treatment effect (survival)
  gte.surv.mean.l <- numeric(L.all)
  gte.surv.sd.l <- numeric(L.all)
  gte.surv.less.gamma0.l <- numeric(L.all)
  
  # Lists of matrices to store summary statistics for region-specific treatment effects
  # (including interactions with time and splines)
  rte.long.mean.l.list <- list(2 + num.knots)
  rte.long.sd.l.list <- list(2 + num.knots)
  rte.long.less.gamma0.l.list <- list(2 + num.knots)
  rte.long.grtr.gamma0.l.list <- list(2 + num.knots)
  for(i in 1:(2 + num.knots)){
    rte.long.mean.l.list[[i]] <- matrix( 0, nrow = L.all, ncol = S )
    rte.long.sd.l.list[[i]] <- matrix( 0, nrow = L.all, ncol = S )
    rte.long.less.gamma0.l.list[[i]] <- matrix( 0, nrow = L.all, ncol = S )
    rte.long.grtr.gamma0.l.list[[i]] <- matrix( 0, nrow = L.all, ncol = S )
    rte.long.names <- colnames(W.X.mat)[((i)*S+num.knots+2):((1+i)*S+num.knots+1)]
    colnames(rte.long.mean.l.list[[i]]) <- rte.long.names
    colnames(rte.long.sd.l.list[[i]]) <- rte.long.names
    colnames(rte.long.less.gamma0.l.list[[i]]) <- rte.long.names
    colnames(rte.long.grtr.gamma0.l.list[[i]]) <- rte.long.names
  }
  
  # Global treatment effects (including interactions with time and splines) (longitudinal)
  gte.long.mean.l <- matrix(0, nrow = L.all, ncol = 2 + num.knots)
  gte.long.sd.l <- matrix(0, nrow = L.all, ncol = 2 + num.knots)
  gte.long.less.gamma0.l <- matrix(0, nrow = L.all, ncol = 2 + num.knots)
  gte.long.grtr.gamma0.l <- matrix(0, nrow = L.all, ncol = 2 + num.knots)
  if(num.knots > 0){
    trt.long.names <- c("trt", "time:trt", paste0("knot.", 1:num.knots, ":trt"))
  } else{
    trt.long.names <- c("trt", "time:trt")
  }
  colnames(gte.long.mean.l) <- trt.long.names
  colnames(gte.long.sd.l) <- trt.long.names
  colnames(gte.long.less.gamma0.l) <- trt.long.names
  colnames(gte.long.grtr.gamma0.l) <- trt.long.names
  
  # Intercepts and time effects (long.)
  intrcpt.time.long.means.l.mat <- matrix( 0, nrow = L.all, ncol = S + 1 + num.knots )
  intrcpt.time.long.sds.l.mat <- matrix( 0, nrow = L.all, ncol = S + 1 + num.knots )
  
  # Covariate effects (surv. and long.)
  covs.surv.means.l.mat <- matrix( 0, nrow = L.all, ncol = p.Y )
  covs.surv.sds.l.mat <- matrix( 0, nrow = L.all, ncol = p.Y) 
  covs.long.means.l.mat <- matrix( 0, nrow = L.all, ncol = p.X )
  covs.long.sds.l.mat <- matrix( 0, nrow = L.all, ncol = p.X )
  
  # Association parameter
  alpha.mean.l <- matrix(0, nrow = L.all, ncol = num.alpha)
  alpha.sd.l <- matrix(0, nrow = L.all, ncol = num.alpha)
  
  # Precision of error in longitudinal likelihood
  tau.means <- numeric(L.all)
  tau.sds <- numeric(L.all)
  
  # Standard deviation(s) (and possibly correlation) of random effects
  if(REs == "int"){
    rand.eff.sds <- numeric(L.all)
  } else if(REs == "int-slope"){
    rand.eff.sds <- matrix(0, nrow = L.all, ncol = 3)
    colnames(rand.eff.sds) <- c("IntSD", "SlopeSD", "Corr")
  }
  
  # Matrices to hold consistency results
  pmda.TE.loc.cons.l <- matrix( 0, nrow = L.all, ncol = S )   # PMDA local consistency for all models (TE)
  pmda.RR.loc.cons.l <- matrix( 0, nrow = L.all, ncol = S )   # PMDA local consistency for all models (RR)
  prws.consis.list <- list()    # list to store pairwise consistency results
  prws.inconsis.list <- list()  # list to store pairwise inconsistency results
  
  # Lists/vectors to store other quantities specific to each model
  theta.Yl.means <- list()        # list to store optimized modes of theta_{Y,l} for each model
  optim.times <- numeric(L.all)   # amount of time in optimization step for each model
  optim.iters <- numeric(L.all)   # number of iterations in optimization step for each model
  
  
  ### Create matrix that is a function of time depending on random effects included
  if(REs == "int"){
    G.inv.dim <- 1                    # number of rows/columns of G^-1
    z.time.mat <- matrix(1, nrow = nrow(W.X.mat), ncol = 1)
    RE.form <- "(1|subjid)"           # random effect piece of lmer model
  } else if(REs == "int-slope"){
    G.inv.dim <- 2                    # number of rows/columns of G^-1
    z.time.mat <- cbind( 1, long.dat$time )
    RE.form <- "(time|subjid)"        # random effect piece of lmer model
  }
  
  
  ### List of subject IDs for each region
  reg.subj.ids <- list()   # identify which subjects belong to which region
  for(i in 1:S){
    reg.subj.ids[[i]] <- surv.dat$subjid[which(surv.dat$region == i)]
  }
  
  
  ### Loop through all models and obtain posterior summaries
  prev.long.mod <- 0         # used to know if W.Xl should be updated and a new lmer model fit
  for(l in 1:L.all){
    
    # Identify survival partition and longitudinal partition
    surv.mod.l <- mod.mat[ mods.fit$surv.mod[l], ]
    long.mod.l <- mod.mat[ mods.fit$long.mod[l], ]
    D.Yl <- length(unique(surv.mod.l))    # number of distinct trt effects for surv. submodel l
    D.Xl <- length(unique(long.mod.l))    # number of distinct trt effects for long. submodel l
    n.surv.l <- numeric(D.Yl)             # vector to store combined sample sizes
    n.long.l <- numeric(D.Xl)             # vector to store combined number of long. observations
    num.regs.surv.l <- numeric(D.Yl)      # vector to store number of regions in each group (surv.)
    num.regs.long.l <- numeric(D.Xl)      # vector to store number of regions in each group (long.)
    
    # Update W.Yl, mu.Yl, and Sig.Yl
    W.Yl <- update.W.Y( surv.dat = surv.dat, surv.mod.l = surv.mod.l )
    mu.Sig.Yl <- update.mu.Sig.Y( mu.Y = mu.Y, Sig.Y = Sig.Y, surv.mod.l = surv.mod.l )
    mu.Yl <- mu.Sig.Yl$mu.Yl
    Sig.Yl <- mu.Sig.Yl$Sig.Yl
    
    # Set initial value(s) of alpha.new = 0 if the first model is being fit for the first time.
    # Otherwise, keep values of alpha.new from previous model as initial values
    if(l == 1){
      alpha.new <- as.matrix( rep(0, num.alpha) )
    }
    
    # Fit new survival model to use as initial values of theta_{Y,l} if current survival
    # partition is being included in joint model for first time (true while long.mod.l = 1).
    # Otherwise, use optimized parameter estimates from previous model
    if(mods.fit$long.mod[l] == 1){
      
      # Initial values of theta_{Y,l} from survival model
      formula.Yl <- paste0("Surv(eventtime, status) ~ as.factor(region) + ",
                           paste(colnames(W.Yl), collapse = " + "))
      surv.mod.init <- coxph(as.formula(formula.Yl), data = cbind(surv.dat, W.Yl))
      theta.Yl.new <- summary(surv.mod.init)$coefficients[S:(S + D.Yl + p.Y - 1),1]
      
    } else{
      
      # Initial values of theta_{Y,l} from former joint model when this current survival
      # partition was first fit (while long.mod.l = 1))
      theta.Yl.new <- theta.Yl.means[[ mods.fit$surv.mod[l] ]]
      
    }
    
    # Update W.Xl, mu.Xl, and Sig.Xl if longitudinal model is different from previous iteration,
    # and fit new longitudinal model to use as initial values of theta_{X,l}.
    # Otherwise, use optimized parameter estimates from previous model
    if(mods.fit$long.mod[l] != prev.long.mod){
      
      # Update design matrix and prior mean vector, and prior covariance matrix
      prev.long.mod <- mods.fit$long.mod[l]
      W.Xl <- update.W.X( long.dat = long.dat, long.mod.l = long.mod.l,
                          lin.spline.knots = lin.spline.knots )
      mu.Sig.Xl <- update.mu.Sig.X( mu.X = mu.X, Sig.X = Sig.X, long.mod.l = long.mod.l,
                                    lin.spline.knots = lin.spline.knots )
      mu.Xl <- mu.Sig.Xl$mu.Xl
      Sig.Xl <- mu.Sig.Xl$Sig.Xl
      
      # Initial values of theta_{X,l}, tau, G^-1, and random effects from longitudinal model (lmer)
      formula.Xl <- paste0("X.ijk ~ -1 + ",
                           paste(c(colnames(W.Xl), RE.form), collapse = " + "))
      long.mod.init <- lmer(as.formula(formula.Xl), data = cbind(long.dat, W.Xl))
      theta.Xl.new <- summary(long.mod.init)$coefficients[,1]
      tau.new <- summary(long.mod.init)$sigma^-2
      b.mat.new <- as.matrix(ranef(long.mod.init)$subjid)
      if(REs == "int"){
        G.inv.new <- matrix( summary(long.mod.init)$varcor$subjid, nrow = 1 )^-1
        b.long.mat.new <- cbind( rep(b.mat.new[,1], table(long.dat$subjid)) )
      } else if(REs == "int-slope"){
        G.inv.new <- matrix( summary(long.mod.init)$varcor$subjid, nrow = 2 )^-1
        b.long.mat.new <- cbind( rep(b.mat.new[,1], table(long.dat$subjid)),
                                 rep(b.mat.new[,2], table(long.dat$subjid)) )
      }
      
    }
    
    # Combine initial values of alpha and theta_{Y,l}
    #alpha.theta.Yl.new <- c(alpha.new, theta.Yl.new)
    
    # Optimization algorithm
    stop.optimization <- 0  # indicator if optimization algorithm should continue (0: yes, 1: no)
    which.iter <- 0         # reset iteration counter
    max.diff.params <- 100  # reset maximum difference in parameters estimates between iterations
    max.diff.b <- 100       # reset maximum difference in rand. effct. estimates between iterations
    start.optim.time <- Sys.time()
    while(stop.optimization == 0){
      
      # Update iteration
      which.iter <- which.iter + 1
      
      # Determine if this iteration should be the last (do not update REs on last iteration)
      # (i.e., (max.diff.params < optim.toler.params & max.diff.b < optim.toler.b) or
      # which.iter == max.iter)
      if( which.iter == max.iter | (max.diff.params < optim.toler.params &
                                    max.diff.b < optim.toler.b) ){
        stop.optimization <- 1
      }
      
      # Calculate mode of alpha (association parameter(s))
      alpha.cur <- alpha.new       # current value of alpha
      alpha.new <- optim( par = alpha.cur, fn = neg.h.alpha, method = "BFGS", S = S, mu.a = mu.alpha,
                          Sig.a = Sig.alpha, eta.iq = eta.iq, phi.iq = phi.iq, eta.tilde = eta.tilde,
                          phi.tilde.summand.list = phi.tilde.summand.list, surv.dat = surv.dat,
                          W.Yl = W.Yl, theta.Yl = theta.Yl.new, b.mat = b.mat.new )$par
      alpha.diff <- max( abs( alpha.new - alpha.cur ) )   # difference between new and old modes

      # Calculate mode of theta_{Y,l} (regression effects for survival submodel)
      theta.Yl.cur <- theta.Yl.new       # current value of theta.Yl
      theta.Yl.new <- optim( par = theta.Yl.cur, fn = neg.h.theta.Yl,
                             method = "BFGS", S = S, mu.Yl = mu.Yl, Sig.Yl = Sig.Yl,
                             eta.iq = eta.iq, phi.iq = phi.iq, eta.tilde = eta.tilde,
                             phi.tilde.summand.list = phi.tilde.summand.list, surv.dat = surv.dat,
                             W.Yl = W.Yl, b.mat = b.mat.new, alpha.vec = alpha.new )$par
      theta.Yl.diff <- max( abs( theta.Yl.new - theta.Yl.cur ) )   # diff. between new and old modes
      
      # # Calculate mode of alpha (association parameter(s)) and theta_{Y,l} (regression effects
      # # for survival submodel)
      # alpha.theta.Yl.cur <- alpha.theta.Yl.new       # current value of alpha and theta_{Y,l}
      # alpha.theta.Yl.new <- optim( par = alpha.theta.Yl.cur, fn = neg.h.alpha.theta.Yl,
      #                              method = "BFGS", S = S, mu.a = mu.alpha, Sig.a = Sig.alpha,
      #                              mu.Yl = mu.Yl, Sig.Yl = Sig.Yl, eta.iq = eta.iq,
      #                              phi.iq = phi.iq, eta.tilde = eta.tilde,
      #                              phi.tilde.summand.list = phi.tilde.summand.list,
      #                              surv.dat = surv.dat, W.Yl = W.Yl, b.mat = b.mat.new )$par
      # alpha.theta.Yl.diff <- max( abs( alpha.theta.Yl.new - alpha.theta.Yl.cur ) )
      # 
      # # Save updated values of alpha and theta_{Y,l}
      # alpha.new <- alpha.theta.Yl.new[1:length(mu.alpha)]          # updated values of alpha
      # theta.Yl.new <- alpha.theta.Yl.new[-c(1:length(mu.alpha))]   # updated values of theta_{Y,l}
      
      # Calculate mode of tau (precision of error in longitudinal likelihood)
      tau.cur <- tau.new
      eta.tilde.tau <- .5 * (nrow(long.dat) + ncol(W.Xl) + eta.tau)
      phi.sum.tau <- sum( (long.dat$X.ijk - W.Xl %*% theta.Xl.new -
                             rowSums(z.time.mat * b.long.mat.new) )^2 )
      phi.tilde.tau <- .5 * ( phi.tau + phi.sum.tau +
                                t(theta.Xl.new - mu.Xl) %*% solve(Sig.Xl) %*% (theta.Xl.new - mu.Xl) )
      tau.new <- (eta.tilde.tau-1)/phi.tilde.tau    # mode of posterior gamma distribution
      tau.diff <- abs(tau.cur - tau.new)   # difference between new and old mode
      
      # Calculate mode of theta_{X,full} (regression effects for longitudinal submodel)
      theta.Xl.cur <- theta.Xl.new
      Sig.tilde.Xl <- solve( t(W.Xl) %*% W.Xl + solve(Sig.Xl) )
      mu.tilde.Xl <- Sig.tilde.Xl %*%
        ( t(W.Xl) %*% (long.dat$X.ijk - rowSums(z.time.mat * b.long.mat.new)) +
            solve(Sig.Xl) %*% mu.Xl )
      theta.Xl.new <- mu.tilde.Xl
      theta.Xl.diff <- max( abs(theta.Xl.cur - theta.Xl.new) )   # diff. between new and old modes
      
      # Calculate mode of G^-1 (inverse of covariance matrix of random effects)
      G.inv.cur <- G.inv.new
      C0.tilde.G <- solve( C0.inv.G + t(b.mat.new) %*% b.mat.new )
      G.inv.new <- as.matrix((nu.tilde.G - G.inv.dim - 1) * C0.tilde.G)  # mode of post. Wishart distn.
      G.inv.diff <- max( abs(G.inv.cur - G.inv.new) )   # difference between new and old mode
      
      # Calculate mode of subject-specific random effects (only if stop.optimization == 0)
      if(stop.optimization == 0){
        
        b.mat.cur <- b.mat.new       # current value of random effects
        for(j in 1:N){
          
          # Identify region and patient number for subject ij
          which.subjid <- surv.dat$subjid[j]     # subject of interest (subject ij)
          subj.reg <- surv.dat$region[j]         # subject's region (region i)
          subj.num <- which(reg.subj.ids[[subj.reg]] == which.subjid)  # which patient in reg. i (1 - n.i)
          b.ij.current <- b.mat.cur[j,]
          
          # Subset data and random effects
          long.rows.ij <- which(long.dat$subjid == which.subjid)      # rows in long.dat for subject ij
          obs.times.i <- surv.dat$eventtime[which(surv.dat$region == subj.reg)]  # TTE times for region i
          surv.dat.ij <- surv.dat[ which(surv.dat$subjid == which.subjid), ]   # survival data for subject ij
          W.Yli <- W.Yl[ which(surv.dat$region == subj.reg), ]   # surv. design matrix for region i
          X.ij <- long.dat$X.ijk[long.rows.ij]                   # observed long. responses for subject ij
          W.Xlij <- W.Xl[ long.rows.ij, ]                        # long. design matrix for subject ij
          z.t.mat.ij <- z.time.mat[long.rows.ij,]
          b.i <- b.mat.cur[ which(surv.dat$region == subj.reg), ]    # random effects for region i
          
          # Find posterior mode of b_ij
          b.ij.new <- optim( par = b.ij.current, fn = neg.h.b.ij, method = "BFGS",
                             reg.i = subj.reg, subj.j = subj.num,
                             G.inv = G.inv.new, eta.i = eta.iq[subj.reg,], phi.i = phi.iq[subj.reg,],
                             eta.tilde.i = eta.tilde[subj.reg,],
                             phi.tilde.summand.i = phi.tilde.summand.list[[subj.reg]],
                             obs.times.i = obs.times.i, surv.dat.ij = surv.dat.ij, W.Yli = W.Yli,
                             X.ij = X.ij, W.Xlij = W.Xlij, z.t.mat.ij = z.t.mat.ij, theta.Yl = theta.Yl.new,
                             theta.Xl = theta.Xl.new, b.i = b.i, alpha.vec = alpha.new, tau.val = tau.new )$par
          
          # Save new mode of random effect(s) for subject j from region i
          b.mat.new[j,] <- b.ij.new
          
        }
        if(REs == "int"){
          b.long.mat.new <- cbind( rep(b.mat.new[,1], table(long.dat$subjid)) )
        } else if(REs == "int-slope"){
          b.long.mat.new <- cbind( rep(b.mat.new[,1], table(long.dat$subjid)),
                                   rep(b.mat.new[,2], table(long.dat$subjid)) )
        }
        max.diff.b <- max( abs( b.mat.new - b.mat.cur ) )   # difference between new and old modes
        
      }
      
      # Calculate maximum difference among all optimized parameters
      max.diff.params <- max( c(alpha.diff, theta.Yl.diff, tau.diff, theta.Xl.diff, G.inv.diff) )
      #max.diff.params <- max( c(alpha.theta.Yl.diff, tau.diff, theta.Xl.diff, G.inv.diff) )
      
    }
    end.optim.time <- Sys.time()
    optim.times[l] <- end.optim.time - start.optim.time
    optim.iters[l] <- which.iter - 1   # number of iterations required for convergence
    
    # Posterior mean and standard deviation for association parameter(s) (alpha)
    # (approximated with Laplace approximation)
    alpha.mean.l[l,] <- alpha.new
    alpha.sd.l[l,] <- diag(-h.pp.alpha( x = alpha.new, S = S, Sig.a = Sig.alpha,
                                        phi.iq = phi.iq, eta.tilde = eta.tilde,
                                        phi.tilde.summand.list = phi.tilde.summand.list,
                                        surv.dat = surv.dat, W.Yl = W.Yl,
                                        theta.Yl = theta.Yl.new, b.mat = b.mat.new ))^-.5
    
    # Posterior means and standard deviations for region-specific treatment effects (surv.)
    # (theta_{Y,l} approximated with Laplace approximation)
    theta.Yl.means[[l]] <- theta.Yl.new    # save theta.Yl means for to later use as initial values
    rte.surv.means.l <- theta.Yl.new[1:D.Yl]
    theta.Yl.sds <- diag(-h.pp.theta.Yl( x = theta.Yl.new, S = S, Sig.Yl = Sig.Yl,
                                         phi.iq = phi.iq, eta.tilde = eta.tilde,
                                         phi.tilde.summand.list = phi.tilde.summand.list,
                                         surv.dat = surv.dat, W.Yl = W.Yl, b.mat = b.mat.new,
                                         alpha.vec = alpha.new ))^-.5
    rte.surv.sds.l <- theta.Yl.sds[1:D.Yl]
    rte.l.d.surv.less.gamma0 <- pnorm(gamma0.surv, rte.surv.means.l, rte.surv.sds.l)
    
    # Posterior means and standard deviations for covariate effects (survival model)
    if(p.Y > 0){
      covs.surv.means.l.mat[l,] <- theta.Yl.new[ (D.Yl+1):(D.Yl+p.Y) ]
      covs.surv.sds.l.mat[l,] <- theta.Yl.sds[ (D.Yl+1):(D.Yl+p.Y) ]
    }
    
    # Posterior means and standard deviations for region-specific treatment effects
    # (including interactions with time and splines) (long.)
    theta.Xl.sds <- diag(as.numeric(tau.new)^-1 * Sig.tilde.Xl)^.5
    rte.long.mean.l <- matrix(0, nrow = 2 + num.knots, ncol = D.Xl)
    rte.long.sd.l <- matrix(0, nrow = 2 + num.knots, ncol = D.Xl)
    rte.long.less.gamma0.l <- matrix(0, nrow = 2 + num.knots, ncol = D.Xl)
    rte.long.grtr.gamma0.l <- matrix(0, nrow = 2 + num.knots, ncol = D.Xl)
    for(k in 1:(2 + num.knots)){
      rte.rows.l <- ((k-1)*D.Xl+S+num.knots+2):((k)*D.Xl+S+num.knots+1)
      rte.long.mean.l[k,] <- theta.Xl.new[ rte.rows.l ]
      rte.long.sd.l[k,] <- theta.Xl.sds[ rte.rows.l ]
      rte.long.less.gamma0.l[k,] <- pnorm(gamma0.long, theta.Xl.new[ rte.rows.l ],
                                          theta.Xl.sds[ rte.rows.l ])
      rte.long.grtr.gamma0.l[k,] <- 1 - rte.long.less.gamma0.l[k,]
    }
    
    # Posterior means and standard deviations for intercepts and time effects (long.)
    intrcpt.time.long.means.l.mat <- theta.Xl.new[ 1:(S+num.knots+1) ]
    intrcpt.time.long.sds.l.mat <- theta.Xl.sds[ 1:(S+num.knots+1) ]
    
    # Posterior means and standard deviations for covariate effects (longitudinal model)
    if(p.X > 0){
      covs.Xl.rows <- (S+num.knots+(2+num.knots)*D.Xl+2):(S+num.knots+(2+num.knots)*D.Xl+1+p.X)
      covs.long.means.l.mat[l,] <- theta.Xl.new[ covs.Xl.rows ]
      covs.long.sds.l.mat[l,] <- theta.Xl.sds[ covs.Xl.rows ]
    }
    
    # Posterior mean and standard deviation for precision of long. likelihood (tau)
    tau.means[l] <- eta.tilde.tau/phi.tilde.tau          # mean of posterior gamma distribution
    tau.sds[l] <- sqrt(eta.tilde.tau/phi.tilde.tau^2)    # SD of posterior gamma distribution
    
    # Posterior standard deviations (and potentially correlation) of random effects (G matrix)
    if(REs == "int"){
      rand.eff.sds[l] <- sqrt( G.inv.new^-1 )
    } else if(REs == "int-slope"){
      G.mat.new <- solve(G.inv.new)
      rand.eff.sds[l, 1:2] <- sqrt( diag(G.mat.new) )
      rand.eff.sds[l, 3] <- G.mat.new[1,2]
    }
    
    # Store stats with regions that share the treatment effect, d=1,...,D_{Y,l} (surv.)
    for(d in 1:D.Yl){
      
      # Store posterior values for each region
      for(k in 1:S){
        if(surv.mod.l[k] == d){
          
          # Posterior summary stats - region-specific treatment effects
          rte.surv.mean.l.mat[l,k] <- rte.surv.means.l[d]
          rte.surv.sd.l.mat[l,k] <- rte.surv.sds.l[d]
          rte.surv.less.gamma0.prob.l.mat[l,k] <- rte.l.d.surv.less.gamma0[d]
          
        }
      }
      
    }
    
    # Store stats with regions that share the treatment effect, d=1,...,D_{X,l} (long.)
    for(d in 1:D.Xl){
      
      # Store posterior values for each region
      for(k in 1:S){
        if(long.mod.l[k] == d){
          
          # Posterior summary stats - region-specific treatment effects
          for(j in 1:(2 + num.knots)){
            rte.long.mean.l.list[[j]][l,k] <- rte.long.mean.l[j,d]
            rte.long.sd.l.list[[j]][l,k] <- rte.long.sd.l[j,d]
            rte.long.less.gamma0.l.list[[j]][l,k] <- rte.long.less.gamma0.l[j,d]
            rte.long.grtr.gamma0.l.list[[j]][l,k] <- rte.long.grtr.gamma0.l[j,d]
          }
          
        }
      }
      
    }
    
    # Calculate combined sample size for regions that share the d^th distinct trtmt effect (surv.)
    which.regs.surv.l <- list()      # list to store groups of region labels
    for(d in 1:D.Yl){
      
      which.regs.surv.Dl <- NULL     # turn into vector to store which regions are the same
      for(i in 1:S){
        if(surv.mod.l[i] == d){
          
          # Combined regional values based on groupings for model M_l
          n.surv.l[d] <- n.surv.l[d] + n.i[i]
          num.regs.surv.l[d] <- num.regs.surv.l[d] + 1
          which.regs.surv.Dl <- c( which.regs.surv.Dl, i )
          
        }
        
        # Store labels of regions that are grouped together
        which.regs.surv.l[[d]] <- which.regs.surv.Dl
        
      }
    }
    
    # Calculate combined sample size for regions that share the d^th distinct trtmt effect (long.)
    which.regs.long.l <- list()      # list to store groups of region labels
    for(d in 1:D.Xl){
      
      which.regs.long.Dl <- NULL     # turn into vector to store which regions are the same
      for(i in 1:S){
        if(long.mod.l[i] == d){
          
          # Combined regional values based on groupings for model M_l
          n.long.l[d] <- n.long.l[d] + num.obs.i[i]
          num.regs.long.l[d] <- num.regs.long.l[d] + 1
          which.regs.long.Dl <- c( which.regs.long.Dl, i )
          
        }
        
        # Store labels of regions that are grouped together
        which.regs.long.l[[d]] <- which.regs.long.Dl
        
      }
    }
    
    # Calculate global treatment effect, summary stats for model M_l (surv.)
    gte.surv.mean.l[l] <- sum( rte.surv.means.l * n.surv.l/N )
    gte.surv.sd.l[l] <- sqrt( sum( (n.surv.l/N)^2 * rte.surv.sds.l^2 ) )
    gte.surv.less.gamma0.l[l] <- pnorm( gamma0.surv, gte.surv.mean.l[l], gte.surv.sd.l[l] )
    
    # Calculate global treatment effects, summary stats for model M_l (long.)
    for(k in 1:(2 + num.knots)){
      gte.long.mean.l[l,k] <- sum( rte.long.mean.l[k,] * n.long.l/N.obs )
      gte.long.sd.l[l,k] <- sum( rte.long.sd.l[k,] * n.long.l/N.obs )
      gte.long.less.gamma0.l[l,k] <- sum( rte.long.less.gamma0.l[k,] * n.long.l/N.obs )
      gte.long.grtr.gamma0.l[l,k] <- sum( rte.long.grtr.gamma0.l[k,] * n.long.l/N.obs )
    }
    
    # Calculate approximate log marginal likelihood for model M_l
    mods.fit[l,6] <- approx.lml( surv.dat = surv.dat, long.dat = long.dat, W.Yl = W.Yl, W.Xl = W.Xl,
                                 eta.iq = eta.iq, phi.iq = phi.iq, eta.tilde = eta.tilde,
                                 phi.tilde.summand.list = phi.tilde.summand.list,
                                 theta.Yl = theta.Yl.new, theta.Xl = theta.Xl.new,
                                 alpha.vec = alpha.new, tau.val = tau.new, b.mat = b.mat.new,
                                 G.inv = G.inv.new, mu.Yl = mu.Yl, Sig.Yl = Sig.Yl, mu.Xl = mu.Xl,
                                 Sig.Xl = Sig.Xl, mu.alpha = mu.alpha, Sig.alpha = Sig.alpha,
                                 eta.tau = eta.tau, phi.tau = phi.tau, nu.G = nu.G,
                                 C0.inv.G = C0.inv.G )
    
    # Obtain Pr(gamma_(l,d)/gamma_G > pi0|D,M_l) for model M_l (PMDA local consistency with
    # treatment effects) and Pr(1 - exp{gamma_(l,d)}) /(1 - exp{gamma_G}) > pi0|D,M_l) for model M_l
    # (PMDA local consistency with risk reduction)
    for(d in 1:D.Yl){
      
      # Calculate Pr(gamma_(l,d)/gamma_G > pi)|D,M_l) for a given value of d
      pmda.TE.loc.cons.ld <- mean( rnorm(n.draws, rte.surv.means.l[d], rte.surv.sds.l[d])/
                                     rnorm(n.draws, gte.surv.mean.l[l], gte.surv.sd.l[l]) > pi0)
      
      # Calculate Pr( (1 - exp{gamma_(l,d)})/(1 - exp{gamma_G}) > pi0|D,M_l) for a given value of d
      pmda.RR.loc.cons.ld <- mean( (1 - exp( rnorm(n.draws, rte.surv.means.l[d],
                                                   rte.surv.sds.l[d]) )) /
                                     (1 - exp( rnorm(n.draws, gte.surv.mean.l[l],
                                                     gte.surv.sd.l[l]) )) > pi0 )
      
      # Store posterior values for each region
      for(k in 1:S){
        if(surv.mod.l[k] == d){
          
          # Posterior summary stats
          pmda.TE.loc.cons.l[l,k] <- pmda.TE.loc.cons.ld
          pmda.RR.loc.cons.l[l,k] <- pmda.RR.loc.cons.ld
          
        }
      }
      
    }
    
    # Obtain pairwise consistency probabilities Pr(|gamma_i - gamma_j| < -log(epsilon_star)|D,M_l)
    # conditional on M_l, l=1,...,L, for all pairwise comparisons between regions and all
    # specified values of epsilon_star. Also obtain pairwise inconsistency probabilities
    # Pr(|gamma_i - gamma_j| > -log(epsilon_star)|D,M_l).
    
    # For model l, create list of matrices storing pairwise consistency/inconsistency
    #    probabilities (separate matrix for all q epsilon_star values)
    prws.comps.list <- pairwise_consistency_surv_l(S = S, modMat_l = surv.mod.l,
                                                   n_draws = n.draws,
                                                   rte_means_l = rte.surv.means.l,
                                                   rte_sds_l = rte.surv.sds.l,
                                                   epsilon_star = epsilon.star)
    prws.consis.list[[l]] <- prws.comps.list[[1]]
    prws.inconsis.list[[l]] <- prws.comps.list[[2]]
    
  }
  
  
  # Calculate posterior model probability for each model
  # Use prior model probabilities \propto exp(D.Yl * alpha0.surv + D.Xl*alpha0.long)
  pmp.mat <- matrix(0, nrow = L.all, ncol = length(a0.surv))
  for(q in 1:length(a0.surv)){
    mod.priors <- exp( mods.fit[,4] * a0.surv[q] + mods.fit[,5] * a0.long[q] ) /
      sum( exp( mods.fit[,4] * a0.surv[q] + mods.fit[,5] * a0.long[q] ) )
    max.log.p.D <- max(mods.fit[,6])
    ml.prop <- exp( mods.fit[,6] - max.log.p.D)
    pmp.vec <- ( ml.prop * mod.priors ) / sum( ml.prop * mod.priors )
    pmp.mat[,q] <- pmp.vec
  }
  
  
  # Save posterior summaries, PMPs, and consistency measures in list
  BMA.results <- list()
  BMA.results$gte.surv.less.gamma0 <- t(pmp.mat) %*% cbind(gte.surv.less.gamma0.l)
  BMA.results$rte.surv.means <- t(pmp.mat) %*% rte.surv.mean.l.mat
  BMA.results$rte.surv.less.gamma0 <- t(pmp.mat) %*% rte.surv.less.gamma0.prob.l.mat
  
  
  # Return list
  return(BMA.results)
  
}
