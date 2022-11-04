// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace std;
using namespace arma;

#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>
#include <complex>
#include <string>



// ### Function to calculate factorial
int factorial(int n){
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}



// ### Function that determines the number of models for specified S, T0
// [[Rcpp::export]]
int num_models(int S, int T0){
  
  // Check function inputs
  if(S <= 0)  Rcpp::stop("S must be greater than 0");
  if(T0 <= 0)  Rcpp::stop("T0 must be greater than 0");
  if(S < T0)  Rcpp::stop("T0 must be less than or equal to S");
  
  
  int nModels = 0;
  double t;
  
  // Calculate number of models
  for(int p = 1; p <= T0; p++){
    
    t = 0;
    for(int j = 0; j < (p+1); j++){
      t = t + std::pow(-1,p-j) * std::pow(j,S) * R::choose(p,j);
    }
    
    nModels = t / (double)factorial(p) + nModels;
    
  }
  return nModels;
  
}




// ### Function that creates a matrix containing lables for all distinct models given values of S, T0
// [[Rcpp::export]]
void recursive_fill(int i, int sMin, int sMax, int & S, int & T0, int & r, Rcpp::IntegerVector & m,
                    Rcpp::IntegerMatrix & models, arma::ivec & nDistinctParms){
  
  for (int k = sMin; k <= sMax; k++)
  {
    m[i] = k;
    if(i+1 == S)
    {
      if(max(m) <= T0)
      {
        models.row(r) = m;
        nDistinctParms[r] = max(m);
        r++;
      }
    }
    else
    {
      int i2 = i + 1;
      int sMin2 = 1;
      int sMax2 = 0;
      for(int s = 0; s <= i; s++){
        if(m[s] > sMax2){
          sMax2 = m[s];
        }
      }
      sMax2++;
      
      recursive_fill(i2, sMin2, sMax2, S, T0, r, m, models, nDistinctParms);	
    }
  }
}



// ### Function that outputs matrix containing lables for all distinct models given a value of S, T0
// [[Rcpp::export]]
Rcpp::IntegerMatrix compile_modMat(int S, int T0){
  
  // Check function inputs
  if(S <= 0)  Rcpp::stop("S must be greater than 0");
  if(T0 <= 0)  Rcpp::stop("T0 must be greater than 0");
  if(S < T0)  Rcpp::stop("T0 must be less than or equal to S");
  
  
  int regionInd = 0;        // region indicator: indicates column of modMat
  int rowNum = 0;           // row number: indicates row of modMat (corresponds to a unique model)
  int sMin = 1;             // starting value used in 'recursive_fill': keep initial value set at 1
  int sMax = 1;             // starting value used in 'recursive_fill': keep initial value set at 1
  int numMods = num_models(S, T0);                     // number of unique models for set S, T0
  Rcpp::IntegerVector currentMod(S);                   // stores distinct param labels for a given model
  Rcpp::IntegerMatrix modMat(numMods, S);              // stores distinct param labels for all models
  arma::ivec nDistinctParms = zeros<ivec>(numMods);    // stores number of distinct params for each model
  
  recursive_fill(regionInd, sMin, sMax, S, T0, rowNum, currentMod, modMat, nDistinctParms);
  
  return modMat;
  
}



// ### Function to create design matrix
// # region_vec: N-dimensional vector of region labels
// # trtmt_vec: N-dimensional vector of treatment indicators (1 = treatment, 0 = control)
// # X_mat: N x p covariate matrix
// [[Rcpp::export]]
Rcpp::NumericMatrix design_mat(Rcpp::CharacterVector & region_vec, Rcpp::IntegerVector & trtmt_vec){
  
  // Number of regions
  Rcpp::CharacterVector unq_reg = sort_unique(region_vec);
  int S = unq_reg.size();
  int N = region_vec.size();
  
  // Re-label region names using numbers
  // ("1" corresponds to first region alphabetically, and "S" to the last)
  Rcpp::IntegerVector reg_vec_num(N);
  for(int i = 1; i <= S; i++){
    for(int k = 0; k < N; k++){
      if(region_vec[k] == unq_reg[i-1]){
        reg_vec_num[k] = i;
      }
    }
  }
  
  // Create region indicator columns (V) and region*trtmt indicator columns (W) for design matrix
  Rcpp::NumericMatrix V(N, S);        // N x S matrix of zeros with region indicator columns
  Rcpp::NumericMatrix W(N, S);        // N x S matrix of zeros with region*trtmt indicator columns
  for(int i = 1; i <= S; i++){
    for(int k = 0; k < N; k++){
      if(reg_vec_num[k] == i){
        V(k,i-1) = 1;
        if(trtmt_vec[k] == 1){
          W(k,i-1) = 1;
        }
      }
    }
  }
  
  // Create design matrix W_mat
  Rcpp::NumericMatrix Z_mat = Rcpp::cbind(V, W);
  
  return Z_mat;
  
}



// ### Function to construct W_l (updated design matrix)
// [[Rcpp::export]]
arma::mat update_W_l(Rcpp::NumericMatrix & W, Rcpp::IntegerVector & modMat_l){
  
  // Extract necessary values
  int N = W.nrow();                                         // number of subjects
  int S = modMat_l.size();                                  // number of regions
  Rcpp::IntegerVector unq_parms_l = sort_unique(modMat_l);
  int T_l = unq_parms_l.size();                             // number of distinct region labels
  int p_cov = W.ncol() - S;                                 // number of covariates (if any)
  
  // Construct W_l (updated design matrix)
  arma::mat W_lt(N, T_l, fill::zeros) ;      // columns will indicate distinct regional trtmt assigments
  
  for(int t = 0; t < T_l; t++){
    for(int i = 0; i < S; i++){
      for(int k = 0; k < N; k++){
        
        if(modMat_l[i] == unq_parms_l[t] && W(k,i) == 1){
          W_lt(k,t) = 1;
        }
        
      }
    }
  }
  
  // Constuct W_l with W_lt and covarate matrix X_mat
  arma::mat W_l(N, T_l + p_cov);
  //arma::mat W_l1 = W( Rcpp::Range(0,N-1), Rcpp::Range(0,S-1) );               // first S cols of W
  
  if(p_cov > 0){
    
    arma::mat X_mat = W( Rcpp::Range(0,N-1), Rcpp::Range(S,W.ncol()-1) );   // covariates matrix
    W_l = arma::join_rows(W_lt, X_mat);
    
  } else{
    
    W_l = W_lt;
    
  }
  
  return W_l;
  
}



// ### Function to construct mu_l (updated location vector for regression priors)
// [[Rcpp::export]]
arma::colvec update_mu_l(Rcpp::NumericVector & mu0, Rcpp::IntegerVector & modMat_l){
  
  // Extract necessary values
  int S = modMat_l.size();                                  // number of regions
  Rcpp::IntegerVector unq_parms_l = sort_unique(modMat_l);
  int T_l = unq_parms_l.size();                             // number of distinct region labels
  int p_cov = mu0.size() - S;                               // number of covariates (if any)
  
  // Construct mu_l - assume mean hyperparameters for regional trtmt effects are same across regions
  arma::colvec mu_l(T_l + p_cov);
  Rcpp::NumericVector mu_l_temp1 = mu0[ Rcpp::Range(0, T_l-1) ];
  
  if(p_cov > 0){
    
    Rcpp::NumericVector mu_l_temp2 = mu0[ Rcpp::Range(S, mu0.size()-1) ];
    
    for(int i = 0; i < mu_l_temp1.size(); i++){
      mu_l(i) = mu_l_temp1(i);
    }
    for(int i = 0; i < mu_l_temp2.size(); i++){
      mu_l( mu_l_temp1.size() + i ) = mu_l_temp2(i);
    }
    
  } else{
    
    for(int i = 0; i < mu_l_temp1.size(); i++){
      mu_l(i) = mu_l_temp1(i);
    }
    
  }
  
  return( mu_l );
  
}



// ### Function to construct Sig_l (updated dispersion matrix for regression priors)
// [[Rcpp::export]]
arma::mat update_Sig_l(Rcpp::NumericMatrix & Sig0, Rcpp::IntegerVector & modMat_l){
  
  // Extract necessary values
  int S = modMat_l.size();                                  // number of regions
  Rcpp::IntegerVector unq_parms_l = sort_unique(modMat_l);
  int T_l = unq_parms_l.size();                             // number of distinct region labels
  int p_cov = Sig0.ncol() - S;                              // number of covariates (if any)
  
  // Construct Sig_l - assume mean hyperparameters for regional trtmt effects are same across regions
  arma::mat Sig_l(T_l + p_cov, T_l + p_cov);
  arma::mat Sig_l11 = Sig0( Rcpp::Range(0,T_l-1), Rcpp::Range(0,T_l-1) );
  
  if(p_cov > 0){
    
    arma::mat Sig_l12 = Sig0( Rcpp::Range(0,T_l-1), Rcpp::Range(S,Sig0.ncol()-1) );
    arma::mat Sig_l22 = Sig0( Rcpp::Range(S,Sig0.ncol()-1), Rcpp::Range(S,Sig0.ncol()-1) );
    Sig_l = arma::join_cols( arma::join_rows(Sig_l11, Sig_l12), arma::join_rows(Sig_l12.t(), Sig_l22) );
    
  } else{
    
    Sig_l = Sig_l11;
    
  }
  
  return( Sig_l );
  
}



// ### Function to calculate global treatment effect
// [[Rcpp::export]]
arma::mat global_effect(arma::colvec & weights_vec, arma::mat & samples_mat){
  arma::mat glob_samp = samples_mat * weights_vec;
  return glob_samp;
}



// ### Function to perform Bayesian model averaging
// [[Rcpp::export]]
double bma_fun(arma::colvec values_vec, arma::colvec post_mod_probs){
  double bma_result = sum( values_vec % post_mod_probs );
  return bma_result;
}



// ### Function to construct m0 hyperparameter vector for normal regression priors
// # User input must include the following:
//     S: number of regions
//     numBin: number of binary covariates (each covariate is iid Bernoulli)
//     numCon: number of continuous covariates (each covariate is iid normal)
//     reg_int_mean: value of m0 corresponding to regional intercepts (assumed to be the same)
//     reg_trt_mean: value of m0 corresponding to regional treatment effects (assumed to be the same)
//     cov_mean: value of m0 corresponding to covariate effects (assumed to be the same)
// [[Rcpp::export]]
arma::colvec construct_m0(int S, int numBin, int numCon, double reg_int_mean, double reg_trt_mean,
                          double cov_mean){
  
  // # Submatrix corresponding to regional intercepts
  arma::colvec m0(2*S + numBin + numCon);
  for(int i = 0; i < S; i++){
    m0[i] = reg_int_mean;
  }
  for(int i = S; i < 2*S; i++){
    m0[i] = reg_trt_mean;
  }
  for(int i = 2*S; i < 2*S + numBin + numCon; i++){
    m0[i] = cov_mean;
  }
  
  return m0;
  
}



// ### Function to construct Sig0 hyperparameter matrix for normal regression priors
// # User input must include the following:
//     S: number of regions
//     numBin: number of binary covariates (each covariate is iid Bernoulli)
//     numCon: number of continuous covariates (each covariate is iid normal)
//     reg_int_diag: diagonal elements of Sig0 corresponding to regional intercepts (assumed to be the same)
//     reg_trt_diag: diagonal elements of Sig0 corresponding to regional treatment effects (assumed to be the same)
//     cov_diag: diagonal elements of Sig0 corresponding to covariate effects (assumed to be the same)
//     reg_int_offdiag: off-diag elements corresponding to regional intercepts (assumed to be the same)
//     reg_trt_offdiag: off-diag elements corresponding to regional treatment effects (assumed to be the same)
//     cov_offdiag: off-diag elements corresponding to covariate effects (assumed to be the same)
//     ri_rt_offdiag: off-diag elements corresponding to reg. ints and reg. trtmt effects (assumed to be the same)
//     ri_cov_offdiag: off-diag elements corresponding to reg. ints and covariates (assumed to be the same)
//     rt_cov_offdiag: off-diag elements corresponding to reg. trtmt effects and covariates (assumed to be the same)
// [[Rcpp::export]]
arma::mat construct_Sig0(int S, int numBin, int numCon, double reg_int_diag, double reg_trt_diag,
                         double cov_diag, double reg_int_offdiag, double reg_trt_offdiag,
                         double cov_offdiag, double ri_rt_offdiag, double ri_cov_offdiag,
                         double rt_cov_offdiag){
  
  // ## Construct submatrices
  
  // # Submatrix corresponding to regional intercepts
  arma::mat Sig0_ri(S, S, fill::ones);
  Sig0_ri = Sig0_ri * reg_int_offdiag;
  Sig0_ri.diag() += reg_int_diag - reg_int_offdiag;
  
  // # Submatrix corresponding to regional treatment effects
  arma::mat Sig0_rt(S, S, fill::ones);
  Sig0_rt = Sig0_rt * reg_trt_offdiag;
  Sig0_rt.diag() += reg_trt_diag - reg_trt_offdiag;
  
  // # Submatrix corresponding to covariate effects
  arma::mat Sig0_cov(numBin + numCon, numBin + numCon, fill::ones);
  Sig0_cov = Sig0_cov * cov_offdiag;
  Sig0_cov.diag() += cov_diag - cov_offdiag;
  
  // # Submatrix corresponding to covariance between regional ints and regional trtmt effects
  arma::mat Sig0_ri_rt(S, S, fill::ones);
  Sig0_ri_rt = Sig0_ri_rt * ri_rt_offdiag;
  
  // # Submatrix corresponding to covariance between regional ints and covariate effects
  arma::mat Sig0_ri_cov(S, numBin + numCon, fill::ones);
  Sig0_ri_cov = Sig0_ri_cov * ri_cov_offdiag;
  
  // # Submatrix corresponding to covariance between regional trtmt effects and covariate effects
  arma::mat Sig0_rt_cov(S, numBin + numCon, fill::ones);
  Sig0_rt_cov = Sig0_rt_cov * rt_cov_offdiag;
  
  
  // ## Combine submatrices to construct Sig0
  arma::mat Sig0_1 = arma::join_rows(Sig0_ri, arma::join_rows(Sig0_ri_rt, Sig0_ri_cov));
  arma::mat Sig0_2 = arma::join_rows(Sig0_ri_rt.t(), arma::join_rows(Sig0_rt, Sig0_rt_cov));
  arma::mat Sig0_3 = arma::join_rows(Sig0_ri_cov.t(), arma::join_rows(Sig0_rt_cov.t(), Sig0_cov));
  arma::mat Sig0 = arma::join_cols(Sig0_1, arma::join_cols(Sig0_2, Sig0_3));
  
  return Sig0;
  
}



// ### Function to construct model priors
// # User input must include the following:
//     S: number of regions
//     T0: maximum number of regions to consider in a given models
//     alpha: tuning parameter, where p(M_l) is proportional to exp{T_l^alpha} for the l^th model
//            and T_l is the number of distinct treatment effects for model M_l
//            If alpha = 0 (default), set each model prior proportional to 1 (i.e., uniform);
//            If alpha > 0, priors favor models with heterogenious treatment effects;
//            If alpha < 0, priors favor models with homorogenious treatment effects;
//     See Biostatistics "Bayesian adaptive basket trial design using model averaging" (Psioda et al., 2019)
// [[Rcpp::export]]
Rcpp::NumericVector construct_modPriors(int S, int T0, double alpha = 0 ){
  
  // Check function inputs
  if(S <= 0)  Rcpp::stop("S must be greater than 0");
  if(T0 <= 0)  Rcpp::stop("T0 must be greater than 0");
  if(S < T0)  Rcpp::stop("T0 must be less than or equal to S");
  
  // # Determine number of models
  int numMods = num_models(S, T0);
  
  // # Initialize model priors
  Rcpp::NumericVector mod_priors(numMods);
  Rcpp::IntegerMatrix mod_mat = compile_modMat(S, T0);
  int num_distinct = 0;
  for(int l = 0; l < numMods; l++){
    num_distinct = unique(mod_mat.row(l)).size();
    mod_priors[l] = exp(num_distinct * alpha);
  }
  mod_priors = mod_priors / sum(mod_priors);
  
  return mod_priors;
  
}



// ### Function to calculate pairwise consis./inconsis. probs. for model M_l, l=1,...,L (surv.)
// # User input must include the following:
//     S: number of regions
//     modMat_l: S x 1 vector of region classification labels for model M_l
//     n_draws: number of samples to draw from posterior distributions
//     rte_means_l: T_l x 1 vector of posterior means of distinct region-specific treatment effects
//     rte_sds_l: T_l x 1 vector of posterior sd's of distinct region-specific treatment effects
//     epsilon_star: vector of possible minimal clinically important differences of trtmnt effects
// [[Rcpp::export]]
Rcpp::List pairwise_consistency_surv_l(int S, Rcpp::IntegerVector modMat_l, int n_draws,
                                       Rcpp::NumericVector rte_means_l,
                                       Rcpp::NumericVector rte_sds_l,
                                       Rcpp::NumericVector epsilon_star){
  
  // Check function inputs
  if(S <= 0)  Rcpp::stop("S must be greater than 0");
  
  // # Determine number of distinct treatment effects for model M_l
  int T_l = rte_means_l.size();
  
  // For model l, create list of matrices storing pairwise consistency/inconsistency probs.
  //    (separate matrix for all q epsilon_star values)
  Rcpp::List prob_prws_consis_l(epsilon_star.size());
  Rcpp::List prob_prws_inconsis_l(epsilon_star.size());
  for(int q = 0; q < epsilon_star.size(); q++){
    prob_prws_consis_l[q] = arma::mat(S, S, fill::ones);
    prob_prws_inconsis_l[q] = arma::mat(S, S, fill::ones);
  }
  
  // Save values for each pairwise comparison that shares the t^th and w^th distinct treatment
  // effects under model M_l
  for(int t = 0; t < T_l; t++){
    
    // Sample n_draws distinct region-specific treatment effects gamma_(l,t)
    Rcpp::NumericVector samp_l_t(n_draws);
    for(int k = 0; k < n_draws; k++){
        samp_l_t(k) = R::rnorm( rte_means_l(t), rte_sds_l(t) );
    }
    
    for(int w = 0; w < T_l; w++){
      
      // Sample n_draws distinct region-specific treatment effects gamma_(l,w)
      Rcpp::NumericVector samp_l_w(n_draws);
      for(int k = 0; k < n_draws; k++){
        samp_l_w(k) = R::rnorm( rte_means_l(w), rte_sds_l(w) );
      }
      
      // Calculate |gamma_(l,t) - gamma_(l,w)| for regions that share the t^th and
      // w^th distinct treatment effects
      arma::vec prws_abs_diff = abs(samp_l_t - samp_l_w);
      
      // Identify which regions share i^th and w^th distinct treatment effects
      for(int i = 0; i < S; i++){
        for(int j = i+1; j < S; j++){
          
          if(modMat_l(i) == t+1){
            if(modMat_l(j) == w+1){
              
              // Calculate Pr(|gamma_(l,t) - gamma_(l,w)| < -log(epsilon_star)|D,M_l)
              // for all epsilon_star values
              for(int q = 0; q < epsilon_star.size(); q++){
                
                // extract q^th matrix for model l to store pw consistency/inconsistency probabilities
                arma::mat pwc_prob_mat_l_q = prob_prws_consis_l[q];
                arma::mat pwic_prob_mat_l_q = prob_prws_inconsis_l[q];
                // vector to store indicators for pw consis. prob.
                arma::colvec less_epsilon_yes(n_draws, fill::zeros);
                for(int k = 0; k < n_draws; k++){
                  if(prws_abs_diff[k] < -log(epsilon_star[q])){
                    less_epsilon_yes[k] = 1.0;
                  }
                }
                
                // Calculate pairwise consistency probabilities between regions i and j for model l
                pwc_prob_mat_l_q(i,j) = mean(less_epsilon_yes);
                pwc_prob_mat_l_q(j,i) = mean(less_epsilon_yes);
                prob_prws_consis_l[q] = pwc_prob_mat_l_q;
                
                // Calculate pairwise inconsistency probabilities between regions i and j for model l
                pwic_prob_mat_l_q(i,j) = 1 - mean(less_epsilon_yes);
                pwic_prob_mat_l_q(j,i) = 1 - mean(less_epsilon_yes);
                prob_prws_inconsis_l[q] = pwic_prob_mat_l_q;
                
              }
              
            }
          }
          
        }
      }
      
    }
  }
  
  // Return S x S matrices with pairwise comparisons as list for model M_l
  Rcpp::List prob_prws_comparisons_l(2);
  prob_prws_comparisons_l(0) = prob_prws_consis_l;
  prob_prws_comparisons_l(1) = prob_prws_inconsis_l;
  return prob_prws_comparisons_l;
  
}



// ### Function to calculate pairwise consistency and inconsistency probabilities (surv. and long.)
// # User input must include the following:
//     S: number of regions
//     pmp: L x 1 vector of posterior model probabilities
//     prws_consis_list: list with all pairwise consistency matrices for each epsilon_star
//          value and each model
//     prws_inconsis_list: list with all pairwise consistency matrices for each epsilon_star
//          value and each model
//     epsilon_star: vector of possible minimal clinically important differences of trtmnt effects
// [[Rcpp::export]]
Rcpp::List pairwise_consistency(int S, Rcpp::NumericVector pmp,
                                Rcpp::List prws_consis_list, Rcpp::List prws_inconsis_list,
                                Rcpp::NumericVector epsilon_star){
  
  // Check function inputs
  if(S <= 0)  Rcpp::stop("S must be greater than 0");
  
  // # Determine number of models
  int numMods = pmp.size();
  
  // List to store BMA pairwise consistency results (one matrix for each epsilon_star value)
  Rcpp::List bma_prws_consis(epsilon_star.size());
  Rcpp::List bma_prws_inconsis(epsilon_star.size());
  for(int q = 0; q < epsilon_star.size(); q++){
    
    Rcpp::NumericMatrix bma_prws_consis_q(S,S);
    Rcpp::NumericMatrix bma_prws_inconsis_q(S,S);
    for(int i = 0; i < S; i++){
      for(int j = 0; j < S; j++){
        
        // Extract Pr(|gamma_i - gamma_j| < -log(epsilon_star)|D, M_l) and
        // Pr(|gamma_i - gamma_j| > -log(epsilon_star)|D, M_l) for all L models
        arma::colvec pwc_probs_ij(numMods, fill::zeros);
        arma::colvec pwic_probs_ij(numMods, fill::zeros);
        for(int l = 0; l < numMods; l++){
          Rcpp::List prws_consis_list_l = prws_consis_list[l];
          Rcpp::List prws_inconsis_list_l = prws_inconsis_list[l];
          arma::mat prws_consis_l_q = prws_consis_list_l[q];
          arma::mat prws_inconsis_l_q = prws_inconsis_list_l[q];
          pwc_probs_ij[l] = prws_consis_l_q(i,j);
          pwic_probs_ij[l] = prws_inconsis_l_q(i,j);
        }
        
        // Posterior Pr(|gamma_i - gamma_j| < -log(epsilon_star)|D) and
        // posterior Pr(|gamma_i - gamma_j| > -log(epsilon_star)|D) 
        bma_prws_consis_q(i,j) = bma_fun(pwc_probs_ij, pmp);
        bma_prws_inconsis_q(i,j) = bma_fun(pwic_probs_ij, pmp);
        
      }
    }
    
    // Save matrix of pairwise consistency probabilities for q^th value of epsilon_star
    bma_prws_consis[q] = bma_prws_consis_q;
    bma_prws_inconsis[q] = bma_prws_inconsis_q;
    
  }
  
  // Return S x S matrices with pairwise comparisons as list for model M_l
  Rcpp::List bma_prws_comps(2);
  bma_prws_comps(0) = bma_prws_consis;
  bma_prws_comps(1) = bma_prws_inconsis;
  return bma_prws_comps;
  
}



// ### Function to calculate epsilon_star-level global consis./inconsis. probs. (surv. and long.)
// # User input must include the following:
//     S: number of regions
//     modMat: L x S matrix of region classification labels for all L models
//     pmp: L x 1 vector of posterior model probabilities
//     bma_prws_consis: list of pairwise consistency matrices for all epsilon_star values
//     bma_prws_inconsis: list of pairwise inconsistency matrices for all epsilon_star values
//     epsilon_star: vector of possible minimal clinically important differences of trtmnt effects
//     beta_star: probability cutoff for global inconsistency for which two regions are
//        considered to be clinically different
// [[Rcpp::export]]
arma::vec global_consistency(int S, Rcpp::IntegerMatrix modMat, arma::vec pmp,
                             Rcpp::List bma_prws_consis, Rcpp::List bma_prws_inconsis,
                             Rcpp::NumericVector epsilon_star, double beta_star){
  
  // # Determine number of models
  int numMods = pmp.size();
  
  // Determine if regional treatment effects are greater than a clinically meaningful regional
  // difference of one another. If so, determine which models allow these regions to differ and
  // store results in which_mods_keep (1 if model should be considered, 0 otherwise)
  arma::mat which_mods_keep(numMods, epsilon_star.size(), fill::zeros);
  for(int q = 0; q < epsilon_star.size(); q++){
    
    arma::mat pwic_mat_q = bma_prws_inconsis[q];    // Extract matrix of pwc probs. for q^th epsilon_star value
    for(int l = 0; l < numMods; l++){
      
      int i_loop_continue = 1;     // Continue with i for-loop if equal to 1; break if 0
      for(int i = 0; i < S; i++){
        
        for(int j = i+1; j < S; j++){
          
          // If Pr(|gamma_i - gamma_j| > -log(epsilon_star) | D) < beta_star, consider
          // Model l in calculation for global inconsistency only if model forces regions i and j
          // to differ
          
          if(pwic_mat_q(i,j) >= beta_star){
            if(modMat(l,i) != modMat(l,j)){
              
              // Consider Model l when adding PMPs for global inconsistency probability
              i_loop_continue = 0;
              which_mods_keep(l,q) = 1;
              break;
              
            }
          }
          
        }
        
        // If j for-loop is stopped early and model M_l is considered, then break i for-loop
        if(i_loop_continue == 0)  break;
        
      }
      
    }
    
  }
  
  // For q^th value of epsilon_star, calculate epsilon_star-level global inconsistency probability
  // for models which require any two regions to differ by more than epsilon_star, with
  // probability beta_star. The epsilon_star-level global consistency probability is
  // defined as 1 - epsilon_star-level global inconsistency probability.
  arma::vec glob_inconsis_prob = which_mods_keep.t() * pmp;
  arma::vec glob_consis_prob = 1.0 - glob_inconsis_prob;
  
  // Return S x S vector with pairwise comparisons as list for model M_l
  return glob_consis_prob;
  
}

