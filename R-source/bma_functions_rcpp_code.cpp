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
Rcpp::NumericMatrix design_mat(Rcpp::CharacterVector & region_vec, Rcpp::IntegerVector & trtmt_vec,
                               Rcpp::NumericMatrix & X_mat){
  
  // Number of regions
  Rcpp::CharacterVector unq_reg = sort_unique(region_vec);
  int S = unq_reg.size();
  int N = region_vec.size();
  
  // Re-label region names using numbers
  // ("1" correpsonds to first region alphabetically, and "S" to the last)
  Rcpp::IntegerVector reg_vec_num(N);
  for(int i = 1; i <= S; i++){
    for(int k = 0; k < N; k++){
      if(region_vec[k] == unq_reg[i-1]){
        reg_vec_num[k] = i;
      }
    }
  }
  
  // Create region indicator columns (V) and region*trtmt indicator columns (Z) for design matrix
  Rcpp::NumericMatrix V(N, S);        // N x S matrix of zeros with region indicator columns
  Rcpp::NumericMatrix Z(N, S);        // N x S matrix of zeros with region*trtmt indicator columns
  for(int i = 1; i <= S; i++){
    for(int k = 0; k < N; k++){
      if(reg_vec_num[k] == i){
        V(k,i-1) = 1;
        if(trtmt_vec[k] == 1){
          Z(k,i-1) = 1;
        }
      }
    }
  }
  
  // Create design matrix W_mat
  Rcpp::NumericMatrix W_mat = Rcpp::cbind(V, Z, X_mat);
  
  return W_mat;
  
}



// ### Function to construct W_l (updated design matrix) and m.l and Sig.l (updated hyperparameters
// ### for regression priors)
// [[Rcpp::export]]
arma::mat update_W_l(Rcpp::NumericMatrix & W, Rcpp::IntegerVector & modMat_l){
  
  // Extract necessary values
  int N = W.nrow();                                         // number of subjects
  int S = modMat_l.size();                                  // number of regions
  Rcpp::IntegerVector unq_parms_l = sort_unique(modMat_l);
  int T_l = unq_parms_l.size();                             // number of distinct region labels
  int p_cov = W.ncol() - 2*S;                               // number of covariates (if any)
  
  // Construct W_l (updated design matrix)
  arma::mat W_lt(N, T_l, fill::zeros) ;      // columns will indicate distinct regional trtmt assigments
  
  for(int t = 0; t < T_l; t++){
    for(int i = 0; i < S; i++){
      for(int k = 0; k < N; k++){
        
        if(modMat_l[i] == unq_parms_l[t] && W(k,(S+i)) == 1){
          W_lt(k,t) = 1;
        }
        
      }
    }
  }
  
  // Constuct W_l with first S columns of W (region indicators), W_lt, and covarate matrix X_mat
  arma::mat W_l(N, S + T_l + p_cov);
  arma::mat W_l1 = W( Rcpp::Range(0,N-1), Rcpp::Range(0,S-1) );               // first S cols of W
  
  if(p_cov > 0){
    
    arma::mat X_mat = W( Rcpp::Range(0,N-1), Rcpp::Range(2*S,W.ncol()-1) );   // covariates matrix
    W_l = arma::join_rows( arma::join_rows(W_l1, W_lt), X_mat);
    
  } else{
    
    W_l = arma::join_rows(W_l1, W_lt);
    
  }
  
  return W_l;
  
}



// ### Function to construct m_l (updated location vector for regression priors)
// [[Rcpp::export]]
arma::colvec update_m_l(Rcpp::NumericVector & m0, Rcpp::IntegerVector & modMat_l){
  
  // Extract necessary values
  int S = modMat_l.size();                                  // number of regions
  Rcpp::IntegerVector unq_parms_l = sort_unique(modMat_l);
  int T_l = unq_parms_l.size();                             // number of distinct region labels
  int p_cov = m0.size() - 2*S;                              // number of covariates (if any)
  
  // Construct m_l - assume mean hyperparameters for regional trtmt effects are same across regions
  arma::colvec m_l(S + T_l + p_cov);
  Rcpp::NumericVector m_l_temp1 = m0[ Rcpp::Range(0, S+T_l-1) ];
  
  if(p_cov > 0){
    
    Rcpp::NumericVector m_l_temp2 = m0[ Rcpp::Range(2*S, m0.size()-1) ];
    
    for(int i = 0; i < m_l_temp1.size(); i++){
      m_l(i) = m_l_temp1(i);
    }
    for(int i = 0; i < m_l_temp2.size(); i++){
      m_l( m_l_temp1.size() + i ) = m_l_temp2(i);
    }
    
  } else{
    
    for(int i = 0; i < m_l_temp1.size(); i++){
      m_l(i) = m_l_temp1(i);
    }
    
  }
  
  return( m_l );
  
}



// ### Function to construct Sig_l (updated dispersion matrix for regression priors)
// [[Rcpp::export]]
arma::mat update_Sig_l(Rcpp::NumericMatrix & Sig0, Rcpp::IntegerVector & modMat_l){
  
  // Extract necessary values
  int S = modMat_l.size();                                  // number of regions
  Rcpp::IntegerVector unq_parms_l = sort_unique(modMat_l);
  int T_l = unq_parms_l.size();                             // number of distinct region labels
  int p_cov = Sig0.ncol() - 2*S;                              // number of covariates (if any)
  
  // Construct Sig_l - assume mean hyperparameters for regional trtmt effects are same across regions
  arma::mat Sig_l(S + T_l + p_cov, S + T_l + p_cov);
  arma::mat Sig_l11 = Sig0( Rcpp::Range(0,S+T_l-1), Rcpp::Range(0,S+T_l-1) );
  
  if(p_cov > 0){
    
    arma::mat Sig_l12 = Sig0( Rcpp::Range(0,S+T_l-1), Rcpp::Range(2*S,Sig0.ncol()-1) );
    arma::mat Sig_l22 = Sig0( Rcpp::Range(2*S,Sig0.ncol()-1), Rcpp::Range(2*S,Sig0.ncol()-1) );
    Sig_l = arma::join_cols( arma::join_rows(Sig_l11, Sig_l12), arma::join_rows(Sig_l12.t(), Sig_l22) );
    
  } else{
    
    Sig_l = Sig_l11;
    
  }
  
  return( Sig_l );
  
}



// ### Function to calculate posterior model probabilities given the data
// [[Rcpp::export]]
arma::colvec calc_pmp(arma::colvec & Y, Rcpp::NumericMatrix & W, Rcpp::NumericVector & m0,
                      Rcpp::NumericMatrix & Sig0, double & delta0, double & nu0,
                      Rcpp::IntegerMatrix & modMat, arma::colvec & modPriors){
  
  // Extract necessary values
  int numMods = modMat.nrow();
  int N = W.nrow();
  
  // Marginal likelihoods for the data given model M_l, l=1,...,L
  arma::colvec l_probDatMl_prop(numMods, fill::zeros);
  for(int l = 1; l <= numMods; l++){
    
    // Obtain W_l, m_l, and Sig_l (updated design matrix and hyperparameters for priors)
    Rcpp::IntegerVector modMat_row_l = modMat.row(l-1);
    arma::mat W_l = update_W_l(W, modMat_row_l);
    arma::colvec m_l = update_m_l(m0, modMat_row_l);
    arma::mat Sig_l = update_Sig_l(Sig0, modMat_row_l);
    
    // Proportional log marginal likelihood for the data given model M_l
    arma::mat I_N = eye(N,N);
    arma::mat I_WSW = I_N + W_l * Sig_l * W_l.t();
    arma::colvec Y_Wm = Y - W_l * m_l;
    double log_doub = log( 1 + 1/nu0 * as_scalar( Y_Wm.t() * I_WSW.i() * Y_Wm ) );
    l_probDatMl_prop[l-1] = -.5*std::real(log_det( I_WSW )) - ( (N+delta0)/2.0 )*log_doub;
    
  }
  
  // Proportional marginal likelihood for the data given model M_l
  double max_log_prob = max(l_probDatMl_prop);
  arma::colvec probDatMl_prop = exp( l_probDatMl_prop - max_log_prob );
  
  // Marginal likelihood for the data given model M_l
  arma::colvec probDatMl = probDatMl_prop / sum(probDatMl_prop);
  
  // Posterior probability for model M_l given the data
  arma::colvec postProbMl = ( probDatMl % modPriors ) / as_scalar( probDatMl.t() * modPriors );
  return postProbMl;
  
}



// ### Function to derive joint posterior distribution of regression parameters given model M.l
// [[Rcpp::export]]
void joint_post_dist_l(arma::colvec & Y, Rcpp::NumericMatrix & W, Rcpp::NumericVector & m0,
                       Rcpp::NumericMatrix & Sig0, double & delta0, double & nu0, Rcpp::IntegerVector & modMat_l,
                       double & df_l, arma::colvec & theta_tilde_l, arma::mat & disp_l){
  
  // Extract necessary values
  int N = W.nrow();
  
  // Obtain W_l, m_l, and Sig_l (updated design matrix and hyperparameters for priors)
  arma::mat W_l = update_W_l(W, modMat_l);
  arma::colvec m_l = update_m_l(m0, modMat_l);
  arma::mat Sig_l = update_Sig_l(Sig0, modMat_l);
  
  // Joint posterior distribution for regression parameters given model M_l
  arma::mat Wl_t_Wl = W_l.t() * W_l;
  arma::colvec theta_hat_l = Wl_t_Wl.i() * W_l.t() * Y;
  arma::mat Lam_l = ( Wl_t_Wl + Sig_l.i() ).i() * Sig_l.i();
  // Degrees of freedom for univariate t distribution for M_l
  df_l = N + delta0;
  // Location vector for univariate t distribution for M_l
  arma::mat I_Wl = eye(W_l.n_cols, W_l.n_cols);
  arma::mat I_N = eye(N,N);
  theta_tilde_l = Lam_l * m_l + (I_Wl - Lam_l) * theta_hat_l;
  // Dispersion matrix for univariate t distribution for M_l
  arma::colvec theta_m = theta_hat_l - m_l;
  double s2_tilde_l = as_scalar( Y.t()*(I_N - W_l*Wl_t_Wl.i()*W_l.t())*Y + theta_m.t()*Lam_l.t()*Wl_t_Wl*theta_m + nu0 ) / (N+delta0);
  disp_l = s2_tilde_l * ( Wl_t_Wl + Sig_l.i() ).i();
  
  // Joint posterior df, location vector, and dispersion matrix given M_l are updated using pointers
  
}



// ### Function that samples from a general univariate t distribution
// [[Rcpp::export]]
arma::vec get_t_samp(int n_draws, double df_val, double m_val, double disp_val){
  
  // Define vectors to store samples
  arma::vec stud_t_samp(n_draws, fill::zeros);     // vector to store draws from Student's t distribution
  arma::vec t_samp(n_draws, fill::zeros);          // vector to store draws from general t distribution
  
  for(int i = 0; i < n_draws; i++){
    // Sample n_draws from Student's t distribution
    stud_t_samp[i] = R::rt(df_val);
    
    // Transform Student's t distribution sample to obtain sample from general univariate t distribution
    t_samp[i] = stud_t_samp[i] * sqrt(disp_val) + m_val;
  }
  
  return t_samp;
  
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



// ### Function to execute BMA with one dataset
// # User input must include the following:
//     Y: response vector
//     region_labs: N-dimensional vector of region labels
//     trtmt_indctr: N-dimensional vector of treatment indicators (1 = treatment, 0 = control)
//     X: N x p covariate matrix
//     T0: maximum number of distinct treatment effects to consider in model space
//     modPriors: vector of priors for each model (e.g., uniform, favor homogeneity)
//     m0: (2S + p)-dimensional mean vector for normal prior for regression parameters
//     Sig0: (2S + p) x (2S + p) covariance matrix for normal prior for regression parameters
//     delta0: shape hyperparameter for gamma prior for tau
//     nu0: rate hyperparameter for gamma prior for tau (expected value = delta0 / nu0)
//     gamma0: value in hypotheses for which we want to calculate probabilities
//     epsilon_star: vector of possible minimal clinically important differences of trtmnt effects - used for consistency measures
//     beta_star: probability in (1 - beta_star) for global inconsistency for which consider two regions to be clinically different
//     pi0: value for which we want to calculate probabilities with local consistency
//     n_draws: number of draws to sample from each posterior distribution
// [[Rcpp::export]]
Rcpp::List bma_single(arma::colvec & Y, Rcpp::CharacterVector & region_labs, Rcpp::IntegerVector & trtmt_indctr,
                      Rcpp::NumericMatrix & X, int T0, arma::colvec & modPriors, Rcpp::NumericVector & m0,
                      Rcpp::NumericMatrix & Sig0, double delta0, double nu0, double gamma0,
                      Rcpp::NumericVector epsilon_star, double beta_star, double pi0, int n_draws){
  
  // Check function inputs
  if( Y.size() != region_labs.size() || Y.size() != trtmt_indctr.size() ){
    Rcpp::stop("Y, region_labs, and trtmt_indctr must be of the same length");
  }
  if( X.nrow() > 0 && X.nrow() != Y.size() ){
    Rcpp::stop("The number of rows in X must equal the length of Y");
  }
  if(T0 <= 0)  Rcpp::stop("T0 must be greater than 0");
  double modPriors_sum = round( sum(modPriors) * 100000000.0 ) / 100000000.0;   // round to 8th decimal
  if( modPriors_sum != 1 )  Rcpp::stop("Model priors must sum to 1");
  if(delta0 <= 0)  Rcpp::stop("delta0 must be greater than 0");
  if(nu0 <= 0)  Rcpp::stop("nu0 must be greater than 0");
  for(int i = 0; i < epsilon_star.size(); i++){
    if(epsilon_star(i) < 0)  Rcpp::stop("epsilon_star must be non-negative");
  }
  if(beta_star < 0 || beta_star > 1)  Rcpp::stop("beta_star must be between 0 and 1");
  if(pi0 < 0 || pi0 > 1)  Rcpp::stop("pi0 must be between 0 and 1");
  if(n_draws <= 0)  Rcpp::stop("n_draws must be greater than 0");
  
  
  
  // ## Step 0: Initial values
  
  // # Construct N x (2S + p) design matrix
  //     First set of S columns correpsond to region (1 if correct region, 0 otherwise)
  //         - order of columns correspond to alphabetical ordering of regions listed in "region_labs"
  //     Second set of S columns correspond to regional treatment assignment
  //         (1 if correct region and in treatment group, 0 otherwise)
  //         - order of columns correspond to alphabetical ordering of regions listed in "region_labs"
  //     Last p columns correpsond to covariates, if any
  Rcpp::NumericMatrix W = design_mat(region_labs, trtmt_indctr, X);
  if( W.ncol() != m0.size() ){
    Rcpp::stop("m0 must be (2*S + p + 1) x 1, where S and p are the number of regions and covariates");
  }
  if( W.ncol() != Sig0.ncol() || W.ncol() != Sig0.nrow() ){
    Rcpp::stop("Sig0 must be (2*S + p + 1) x (2*S + p + 1), where S and p are the number of regions and covariates");
  }
  
  
  // Regions
  Rcpp::CharacterVector unq_reg = sort_unique(region_labs);   // region labels
  int S = unq_reg.size();                                     // number of regions
  Rcpp::IntegerVector n_i(S);
  for(int i = 0; i < S; i++){
    n_i[i] = sum(W.column(i));                                // vector of regional sample sizes
  }
  if(S < T0)  Rcpp::stop("T0 must be less than or equal to S");
  if( num_models(S, T0) != modPriors.size() ){
    Rcpp::stop("The number of model priors must equal the number of possible models (determined by S and T0)");
  }
  
  
  // Covariates
  int p_cov = X.ncol();                                       // number of covariate effects
  Rcpp::CharacterVector cov_names;
  if(p_cov > 0){
    cov_names = colnames(X);                                  // covariate names
  }
  
  
  // Model prior specification
  int numMods = num_models(S, T0);                            // L, number of models in model space
  
  
  // Subset data to remove observations from i^th region for leave-one-out computations
  Rcpp::List W_list(S);                      // list of size S to store altered W matrices
  Rcpp::List Y_list(S);                      // list of size S to store altered Y vectors
  
  for(int i = 0; i < S; i++){
    
    // Identify observations in i^th region and remove rows from Y and W
    arma::colvec Y_loo(Y.n_elem - n_i[i], fill::zeros);      // vector to store Y without i^th region (leave-one-out)
    Rcpp::NumericMatrix W_loo(W.nrow() - n_i[i], W.ncol());  // matrix to store W without i^th region (leave-one-out)
    int reg_i_count = 0;
    for(int j = 0; j < W.nrow(); j++){
      
      if( W(j,i) == 0 ){
        Y_loo[reg_i_count] = Y[j];
        W_loo.row(reg_i_count) = W.row(j);
        reg_i_count = reg_i_count + 1;
      }
      
    }
    
    // Identify i^th and (i+S)^th columns of W_loo (now contains only 0's) and remove
    Rcpp::NumericMatrix W_loo2(W_loo.nrow(), W.ncol() - 2);
    int col_count = 0;
    for(int k = 0; k < W.ncol(); k++){
      
      if( k != i && k != i+S ){
        W_loo2.column(col_count) = W_loo.column(k);
        col_count = col_count + 1;
      }
      
    }
    
    Y_list[i] = Y_loo;
    W_list[i] = W_loo2;
    
  }
  
  
  
  // ## Step 1: Construct model matrix
  Rcpp::IntegerMatrix modMat = compile_modMat(S, T0);         // L x S matrix
  
  
  
  // ## Step 2: Calculate posterior model probabilities (L-dimensional vector)
  arma::colvec pmp = calc_pmp(Y, W, m0, Sig0, delta0, nu0, modMat, modPriors);
  
  
  
  // ## Step 3: Calculate joint posterior distribution for each model M_l, and then sample n_draws from
  // ##         marginal posterior distributions of gamma_(l,t) conditional on M_l, t=1,...,T_l, l=1,...,L.
  // ##         Sample from posterior distribution of global treatment effect and posterior distributions
  // ##         of other regression parameters (regional intercepts and covariate effects) conditional on M_l.
  
  
  // Matrices to store posterior summaries of regional treatment effects for each model
  arma::mat post_mean_rte(numMods, S, fill::zeros);         // posterior means of gamma_(l,t)|D,M_l
  arma::mat post_sd_rte(numMods, S, fill::zeros);           // posterior standard dev. of gamma_(l,t)|D,M_l
  arma::mat prob_grtr_gamma0_rte(numMods, S, fill::zeros);  // Pr(gamma_(l,t) > gamma0|D,M_l)
  arma::mat prob_grtr_pi0_consis(numMods, S, fill::zeros);  // Pr(gamma_(l,t)/gamma > pi0|D,M_l)
  Rcpp::List prob_less_eps_consis_loo(epsilon_star.size()); // list of size epsilon_star.size() for LOO local consis. probs.
  for(int q = 0; q < epsilon_star.size(); q++){
    arma::mat prob_loo_q(numMods, S, fill::zeros);          // save empty L x S matrix in prob_less_eps_consis_loo list
    prob_less_eps_consis_loo[q] = prob_loo_q;               //    for each value of epsilon_star
  }
  Rcpp::List prws_consis_list(numMods);                     // list of size numMods for pairwise consistency probabilities
  Rcpp::List prws_inconsis_list(numMods);                   // list of size numMods for pairwise consistency probabilities
  
  
  // Vectors to store posterior summaries of global treatment effect for each model
  arma::colvec post_mean_glob(numMods, fill::zeros);        // posterior mean of gamma|D,M_l
  arma::colvec post_sd_glob(numMods, fill::zeros);          // posterior standard dev. of gamma|D,M_l
  arma::colvec prob_grtr_gamma0_glob(numMods, fill::zeros); // Pr(gamma > gamma0|D,M_l)
  
  
  // Matrices to store posterior summaries of regional intercepts and covariate effects for each model
  arma::mat post_mean_ri(numMods, S, fill::zeros);          // posterior means of mu_i
  arma::mat post_sd_ri(numMods, S, fill::zeros);            // posterior standard dev. of mu_i
  arma::mat post_mean_cov(numMods, p_cov, fill::zeros);     // posterior means of beta_p
  arma::mat post_sd_cov(numMods, p_cov, fill::zeros);       // posterior standard dev. of beta_p
  
  
  // Loop through all possible models
  for(int l = 0; l < numMods; l++){
    
    Rcpp::IntegerVector modMat_l = modMat.row(l);      // store distinct trtmt effect labels for model M_l
    Rcpp::IntegerVector unq_parms_l = sort_unique(modMat_l);
    int T_l = unq_parms_l.size();                      // number of unique treatment effects for model M_l
    arma::mat samp_l_t(n_draws, T_l, fill::zeros);     // matrix to store posterior draws of regional trtmt effects
    arma::mat samp_l_mui(n_draws, S, fill::zeros);     // matrix to store posterior draws of regional intercepts
    arma::mat samp_l_b(n_draws, p_cov, fill::zeros);   // matrix to store posterior draws of covariate effects
    arma::colvec n_l(T_l, fill::zeros);                // vector to store combined sample sizes
    arma::mat consis_ratio(n_draws, T_l, fill::zeros); // matrix to store gamma_(l,t)/gamma used for local consistency
    arma::mat consis_diff_loo(n_draws, S, fill::zeros); // matrix to store |gamma_(l,t) - gamma_(-i)| used for local consistency
    
    
    // # Step 3.1: Calculate joint posterior distribution for model M_l using all regions.
    // #           Save degrees of freedom, location vector, and dispersion matrix of joint
    // #           posterior distribution
    double df_l = 0.0;
    arma::colvec theta_tilde_l(W.ncol() - S + T_l, fill::zeros);
    arma::mat disp_l(W.ncol() - S + T_l, W.ncol() - S + T_l, fill::zeros);
    joint_post_dist_l(Y, W, m0, Sig0, delta0, nu0, modMat_l, df_l, theta_tilde_l, disp_l);
    
    
    // Loop through distinct treatment effect parameters gamma_(l,t) for model M_l
    for(int t = 0; t < T_l; t++){
      
      // Save values for t distribution for gamma_(l,t)
      double theta_l_t = theta_tilde_l(S+t);
      double disp_l_t = disp_l(S+t,S+t);
      
      
      // # Step 3.2: Sample n_draws from marginal posterior distribution of gamma_(l,t) conditional
      // #           on M_l, t=1,...,T_l, l=1,...,L
      samp_l_t.col(t) = get_t_samp(n_draws, df_l, theta_l_t, disp_l_t);
      
      
      // Calculate summary stats for posterior distribution of gamma_(l,t) and store stats with regions
      // that share the treatment effect gamma_(l,t), t=1,...,T_l, l=1,...,L
      for(int i = 0; i < S; i++){
        
        // Save values for all regions that share the t^th distinct treatment effect for model M_l
        if(modMat(l,i) == t+1){
          
          // Calculate combined sample size for regions that share the t^th distinct treatment effect
          n_l[t] = n_l[t] + n_i[i];
          
          // Calculate posterior mean for regional treatment effects
          post_mean_rte(l,i) = mean(samp_l_t.col(t));
          
          // Calculate posterior standard deviation for regional treatment effects
          post_sd_rte(l,i) = stddev(samp_l_t.col(t));
          
          // Calculate probability that gamma_(l,t) is greater than gamma0 for regional treatment effects
          arma::vec grtr_gamma0_yes_rte(n_draws, fill::zeros);
          for(int k = 0; k < n_draws; k++){
            if(samp_l_t(k,t) > gamma0){
              grtr_gamma0_yes_rte[k] = 1.0;
            }
          }
          prob_grtr_gamma0_rte(l,i) = mean(grtr_gamma0_yes_rte);
          
        }
        
      }
      
    }
    
    
    // # Step 3.3: Obtain sample of n_draws for global treatment effect conditional on M_l, l=1,...,L
    
    double N_doub = W.nrow();
    arma::colvec samp_size_weights = n_l / N_doub;
    arma::mat global_samp_l = global_effect(samp_size_weights, samp_l_t);
    
    // Calculate posterior mean and standard deviation of global treatment effect for model M_l
    post_mean_glob[l] = mean(global_samp_l.col(0));
    post_sd_glob[l] = stddev(global_samp_l.col(0));
    
    // Calculate Pr(gamma > gamma0|D, M_l) for model M_l
    arma::vec grtr_gamma0_yes_glob(n_draws, fill::zeros);
    for(int k = 0; k < n_draws; k++){
      if(global_samp_l(k,0) > gamma0){
        grtr_gamma0_yes_glob[k] = 1.0;
      }
    }
    prob_grtr_gamma0_glob[l] = mean(grtr_gamma0_yes_glob);
    
    
    // # Step 3.4: Obtain samples of n_draws for other regression parameters conditional on M_l, l=1,...,L
    
    // Step 3.4.1: Sample n_draws from marginal posterior distribution for regional intercepts mu_i,
    //             i=1,...,S, conditional on M_l, l=1,...,L
    for(int i = 0; i < S; i++){
      
      // Save values for t distribution for regional intercept mu_i conditional on M_l
      double theta_l_mui = theta_tilde_l(i);
      double disp_l_mui = disp_l(i,i);
      
      // Obtain sample of n_draws from t distribution for regional intercept mu_i conditional on M_l
      samp_l_mui.col(i) = get_t_samp(n_draws, df_l, theta_l_mui, disp_l_mui);
      
      // Calculate posterior mean and standard deviation of regional intercept mu_i for model M_l
      post_mean_ri(l,i) = mean(samp_l_mui.col(i));
      post_sd_ri(l,i) = stddev(samp_l_mui.col(i));
      
    }
    
    
    // Step 3.4.2: Sample n_draws from marginal posterior distribution for covariate effects
    //             conditional on M_l, l=1,...,L
    if(p_cov > 0){
      
      for(int p = 0; p < p_cov; p++){
        
        // Save values for t distribution for covariate effect beta_p conditional on M_l
        double theta_l_betap = theta_tilde_l(S + T_l + p);
        double disp_l_betap = disp_l(S + T_l + p, S + T_l + p);
        
        // Obtain sample of n_draws from t distribution for covariate effect beta_p conditional on M_l
        samp_l_b.col(p) = get_t_samp(n_draws, df_l, theta_l_betap, disp_l_betap);
        
        // Calculate posterior mean and standard deviation of covariate effect beta_p for model M_l
        post_mean_cov(l,p) = mean(samp_l_b.col(p));
        post_sd_cov(l,p) = stddev(samp_l_b.col(p));
        
      }
      
    }
    
    
    // # Step 3.5: Obtain Pr(gamma_(l,t)/gamma > pi0|D,M_l) conditional on M_l, l=1,...,L
    
    // Save values for all regions that share the t^th distinct treatment effect for model M_l
    for(int t = 0; t < T_l; t++){
      
      // Calculate ratio gamma_(l,t)/gamma for regions that share the t^th distinct treatment effect
      consis_ratio.col(t) = samp_l_t.col(t) / global_samp_l;
      
      // Calculate Pr(gamma_(l,t)/gamma > pi0|D,M_l) for each region
      for(int i = 0; i < S; i++){
        
        if(modMat(l,i) == t+1){
          
          // Calculate probability gamma_(l,t)/gamma that ratio is greater than pi0
          arma::vec grtr_pi0_yes_consis(n_draws, fill::zeros);
          for(int k = 0; k < n_draws; k++){
            if(consis_ratio(k,t) > pi0){
              grtr_pi0_yes_consis[k] = 1.0;
            }
          }
          prob_grtr_pi0_consis(l,i) = mean(grtr_pi0_yes_consis);
          
        }
        
      }
      
    }
    
    
    // # Step 3.6: Calculate leave-one-out global treatment effects for model M_l without
    // #           i^th region. Save degrees of freedom, location vector, and dispersion matrix
    // #           of joint posterior distribution
    
    // Create matrix to store samples of leave-one-out global treatment effect without i^th region, i=1,...,S
    arma::mat global_samp_l_loo(n_draws, S, fill::zeros);
    
    // Loop through S regions
    for(int i = 0; i < S; i++){
      
      // Update Y, W, modMat_l, and T_l to account for removed i^th region
      arma::colvec Y_loo = Y_list[i];                   // Y without i^th region (leave-one-out)
      Rcpp::NumericMatrix W_loo = W_list[i];            // W without i^th region (leave-one-out)
      Rcpp::IntegerVector modMat_l_loo = modMat_l;      // store distinct trtmt effect labels for model M_l
      modMat_l_loo.erase(i);
      Rcpp::IntegerVector unq_parms_l_loo = sort_unique(modMat_l_loo);
      int T_l_loo = unq_parms_l_loo.size();
      
      // Update m0 to account for removed i^th region
      Rcpp::NumericVector m0_loo = m0;                  // delete i^th and (S+i)^th elements of m0
      m0_loo.erase(S+i);
      m0_loo.erase(i);
      
      // Update Sig0 to account for removed i^th region
      Rcpp::NumericMatrix Sig0_loo(Sig0.nrow() - 2, Sig0.ncol() - 2);
      int row_j_loo = 0;
      int col_k_loo = 0;
      for(int row_j = 0; row_j < Sig0.nrow(); row_j++){
        
        if(row_j != i && row_j != i+S){
          for(int col_k = 0; col_k < Sig0.ncol(); col_k++){
            
            if(col_k != i && col_k != i+S){
              Sig0_loo(row_j_loo,col_k_loo) = Sig0(row_j,col_k);
              col_k_loo = col_k_loo + 1;
            }
            
          }
          col_k_loo = 0;
          row_j_loo = row_j_loo + 1;
        }
        
      }
      
      // Joint posterior distribution without i^th region
      double df_l_loo = 0.0;
      arma::colvec theta_tilde_l_loo(W_loo.ncol() - (S-1) + T_l_loo, fill::zeros);
      arma::mat disp_l_loo(W_loo.ncol() - (S-1) + T_l_loo, W_loo.ncol() - (S-1) + T_l_loo, fill::zeros);
      joint_post_dist_l(Y_loo, W_loo, m0_loo, Sig0_loo, delta0, nu0, modMat_l_loo, df_l_loo, theta_tilde_l_loo, disp_l_loo);
      
      
      // # Step 3.6.1: Sample n_draws from marginal posterior distribution of gamma_(l,t) conditional
      // #             on M_l, t=1,...,T_l_loo, l=1,...,L
      
      arma::mat samp_l_t_loo(n_draws, T_l_loo, fill::zeros);    // matrix to store posterior draws of regional trtmt effects
      
      // Loop through distinct treatment effect parameters gamma_(l,t) for model M_l (without i^th region)
      for(int t = 0; t < T_l_loo; t++){
        
        // Save values for t distribution for gamma_(l,t)
        double theta_l_t_loo = theta_tilde_l_loo(S-1+t);
        double disp_l_t_loo = disp_l_loo(S-1+t,S-1+t);
        
        // Sample n_draws from marginal posterior distribution of gamma_(l,t)
        samp_l_t_loo.col(t) = get_t_samp(n_draws, df_l_loo, theta_l_t_loo, disp_l_t_loo);
        
      }
      
      
      // # Step 3.6.2: Obtain sample of n_draws for leave-one-out global treatment effect missing
      // #             i^th region conditional on M_l, l=1,...,L
      
      // Calculate combined sample size for regions that share the t^th distinct treatment effect
      // gamma_(l,t), t=1,...,T_l_loo, l=1,...,L. Exclude i^th region.
      arma::colvec n_l_loo = n_l;               // vector to store combined sample sizes
      int t_i = modMat_l[i];
      n_l_loo[t_i-1] = n_l_loo[t_i-1] - n_i[i];
      if(n_l_loo[t_i-1] == 0){
        n_l_loo.shed_row(t_i-1);
      }
      
      double N_loo = W_loo.nrow();
      arma::colvec samp_size_weights_loo = n_l_loo / N_loo;
      global_samp_l_loo.col(i) = global_effect(samp_size_weights_loo, samp_l_t_loo);
      
    }
    
    
    // # Step 3.7: Obtain Pr(|gamma_(l,t) - gamma_(-i)| < epsilon_star|D,M_l) conditional on M_l, l=1,...,L
    
    // Save values for each region under model M_l
    for(int i = 0; i < S; i++){
      
      // Identify distinct treatment effect number for i^th region
      int t_i = modMat_l[i];
      
      // Calculate difference |gamma_(l,t) - gamma_(-i)| for each region
      consis_diff_loo.col(i) = abs( samp_l_t.col(t_i-1) - global_samp_l_loo.col(i) );
      
    }
    
    // Calculate probability that differnece |gamma_(l,t) - gamma_(-i)| is less than
    // epsilon_star for each value of epsilon_star
    for(int q = 0; q < epsilon_star.size(); q++){
      
      // Extract L x S matrix corresponding to q^th value of epsilon_star to store loc. cons. probs.
      arma::mat prob_less_eps_consis_loo_q = prob_less_eps_consis_loo[q];
      for(int i = 0; i < S; i++){
        
        arma::vec less_eps_yes_consis_loo(n_draws, fill::zeros);
        for(int k = 0; k < n_draws; k++){
          if(consis_diff_loo(k,i) < epsilon_star[q]){
            less_eps_yes_consis_loo[k] = 1.0;
          }
        }
        prob_less_eps_consis_loo_q(l,i) = mean(less_eps_yes_consis_loo);
        
      }
      
      // Store L x S matrix for q^th value of epsilon_star into list
      prob_less_eps_consis_loo[q] = prob_less_eps_consis_loo_q;
      
    }
    
    
    // # Step 3.8: Obtain pairwise consistency probababilities
    // #           Pr(|gamma_i - gamma_j| < epsilon_star|D,M_l) conditional on M_l, l=1,...,L,
    // #           for all pairwise comparisons between regions and all specified values of epsilon_star.
    // #           Also obtain pairwise inconsistency probabilities
    // #           Pr(|gamma_i - gamma_j| > epsilon_star|D,M_l).
    
    // For model l, create list of matrices storing pairwise consistency/inconsistency probs.
    //    (separate matrix for all q epsilon_star values)
    Rcpp::List prob_prws_consis_l(epsilon_star.size());
    Rcpp::List prob_prws_inconsis_l(epsilon_star.size());
    for(int q = 0; q < epsilon_star.size(); q++){
      prob_prws_consis_l[q] = arma::mat(S, S, fill::ones);
      prob_prws_inconsis_l[q] = arma::mat(S, S, fill::ones);
    }
    
    // Save values for each pairwise comparison that shares the t^th and w^th distinct treatment effects under model M_l
    for(int t = 0; t < T_l; t++){
      for(int w = 0; w < T_l; w++){
        
        // Calculate |gamma_(l,t) - gamma_(l,w)| for regions that share the t^th and w^th distinct treatment effects
        arma::vec prws_abs_diff = abs(samp_l_t.col(t) - samp_l_t.col(w));    // |gamma_(l,t) - gamma_(l,w)|
        
        // Identify which regions share i^th and w^th distinct treatment effects
        for(int i = 0; i < S; i++){
          for(int j = i+1; j < S; j++){
            
            if(modMat(l,i) == t+1){
              if(modMat(l,j) == w+1){
                
                // Calculate Pr(|gamma_(l,t) - gamma_(l,w)| < epsilon_star|D,M_l) for all epsilon_star values
                for(int q = 0; q < epsilon_star.size(); q++){
                  
                  arma::mat pwc_prob_mat_l_q = prob_prws_consis_l[q];    // extract q^th matrix for model l to stare pw consis. probs.
                  arma::mat pwic_prob_mat_l_q = prob_prws_inconsis_l[q]; // extract q^th matrix for model l to stare pw inconsis. probs.
                  arma::colvec less_epsilon_yes(n_draws, fill::zeros);   // vector to store indicators for pw consis. prob.
                  for(int k = 0; k < n_draws; k++){
                    if(prws_abs_diff[k] < epsilon_star[q]){
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
    
    // Store S x S matrices with pairwise comparisons into lists
    prws_consis_list[l] = prob_prws_consis_l;
    prws_inconsis_list[l] = prob_prws_inconsis_l;
    
  }
  
  
  
  // ## Step 4: Calculate posterior means, standard deviations, and probabilities using BMA
  
  
  // # Step 4.1: Calculate posterior means, standard deviations, and probabilities for regional trtmt effects
  
  // Matrix to store BMA results for all S regions
  //    column 1: posterior mean
  //    column 2: posterior standard deviation
  //    column 3: posterior Pr(gamma_i > gamma0|D)
  //    column 4: posterior Pr(gamma_i/gamma > pi0|D)
  Rcpp::NumericMatrix bma_vals_rte(S,4);
  Rcpp::colnames(bma_vals_rte) = Rcpp::CharacterVector::create("PostMean", "PostSD", "PostProbGamma0", "MHLWProbLC");
  Rcpp::rownames(bma_vals_rte) = unq_reg;
  
  for(int i = 0; i < S; i++){
    
    // Posterior mean for regional treatment effects
    bma_vals_rte(i,0) = bma_fun(post_mean_rte.col(i), pmp);
    
    // Posterior standard deviation for regional treatment effects
    bma_vals_rte(i,1) = bma_fun(post_sd_rte.col(i), pmp);
    
    // Posterior Pr(gamma_i > gamma0|D) for regional treatment effects
    bma_vals_rte(i,2) = bma_fun(prob_grtr_gamma0_rte.col(i), pmp);
    
    // Posterior Pr(gamma_i/gamma > pi0|D) for local consistency assessment
    bma_vals_rte(i,3) = bma_fun(prob_grtr_pi0_consis.col(i), pmp);
    
  }
  
  // Matrix to store BMA LOO local consistency results Pr(|gamma_i - gamma_(-i)| < epsilon_star|D)
  //    (one row for each region, one column for each epsilon_star value)
  Rcpp::NumericMatrix bma_LOO_loc_cons(S,epsilon_star.size());
  Rcpp::colnames(bma_LOO_loc_cons) = epsilon_star;
  Rcpp::rownames(bma_LOO_loc_cons) = unq_reg;
  for(int q = 0; q < epsilon_star.size(); q++){
    
    arma::mat prob_less_eps_consis_loo_q = prob_less_eps_consis_loo[q];
    for(int i = 0; i < S; i++){
      
      // Posterior Pr(|gamma_i - gamma_(-i)| < epsilon_star|D) for leave-one-out local consistency assessment
      bma_LOO_loc_cons(i,q) = bma_fun(prob_less_eps_consis_loo_q.col(i), pmp);
      
    }
    
  }
  
  // List to store BMA pairwise consistency results (one matrix for each epsilon_star value)
  Rcpp::List bma_prws_consis(epsilon_star.size());
  Rcpp::List bma_prws_inconsis(epsilon_star.size());
  for(int q = 0; q < epsilon_star.size(); q++){
    
    Rcpp::NumericMatrix bma_prws_consis_q(S,S);
    Rcpp::NumericMatrix bma_prws_inconsis_q(S,S);
    for(int i = 0; i < S; i++){
      for(int j = 0; j < S; j++){
        
        // Extract Pr(|gamma_i - gamma_j| < epsilon_star|D, M_l) and
        // Pr(|gamma_i - gamma_j| > epsilon_star|D, M_l) for all L models
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
        
        // Posterior Pr(|gamma_i - gamma_j| < epsilon_star|D) and
        // posterior Pr(|gamma_i - gamma_j| > epsilon_star|D) 
        bma_prws_consis_q(i,j) = bma_fun(pwc_probs_ij, pmp);
        bma_prws_inconsis_q(i,j) = bma_fun(pwic_probs_ij, pmp);
        
      }
    }
    
    // Save matrix of pairwise consistency probabilities for q^th value of epsilon_star
    bma_prws_consis[q] = bma_prws_consis_q;
    bma_prws_inconsis[q] = bma_prws_inconsis_q;
    
  }
  
  
  // # Step 4.2: Calculate posterior means, standard deviations, and probabilities for global trtmt effect
  
  Rcpp::NumericMatrix bma_vals_glob(3,1);
  Rcpp::rownames(bma_vals_glob) = Rcpp::CharacterVector::create("PostMean", "PostSD", "PostProbGamma0");
  Rcpp::colnames(bma_vals_glob) = Rcpp::CharacterVector::create("Global");
  
  // Posterior mean for global treatment effect
  bma_vals_glob(0,0) = bma_fun(post_mean_glob, pmp);
  
  // Posterior standard deviation for global treatment effect
  bma_vals_glob(1,0) = bma_fun(post_sd_glob, pmp);
  
  // Posterior probability Pr(gamma > gamma0|D) for global treatment effect
  bma_vals_glob(2,0) = bma_fun(prob_grtr_gamma0_glob, pmp);
  
  
  // # Step 4.3: Calculate posterior means, standard deviations, and probabilities for regional intercepts
  
  // Matrix to store BMA results for all S regions
  //    column 1: posterior mean
  //    column 2: posterior standard deviation
  Rcpp::NumericMatrix bma_vals_ri(S,2);
  Rcpp::colnames(bma_vals_ri) = Rcpp::CharacterVector::create("PostMean", "PostSD");
  Rcpp::rownames(bma_vals_ri) = unq_reg;
  
  for(int i = 0; i < S; i++){
    
    // Posterior mean for regional intercepts
    bma_vals_ri(i,0) = bma_fun(post_mean_ri.col(i), pmp);
    
    // Posterior standard deviation for regional intercepts
    bma_vals_ri(i,1) = bma_fun(post_sd_ri.col(i), pmp);
    
  }
  
  
  // # Step 4.4: Calculate posterior means, standard deviations, and probabilities for covariate effects
  
  // Matrix to store BMA results for all p_cov regions
  //    column 1: posterior mean
  //    column 2: posterior standard deviation 
  Rcpp::NumericMatrix bma_vals_cov(p_cov,2);
  Rcpp::colnames(bma_vals_cov) = Rcpp::CharacterVector::create("PostMean", "PostSD");
  Rcpp::rownames(bma_vals_cov) = cov_names;
  
  if(p_cov > 0){
    
    for(int p = 0; p < p_cov; p++){
      
      // Posterior mean for covariate effects
      bma_vals_cov(p,0) = bma_fun(post_mean_cov.col(p), pmp);
      
      // Posterior standard deviation for covariate effects
      bma_vals_cov(p,1) = bma_fun(post_sd_cov.col(p), pmp);
      
    }
    
  }
  
  
  // # Step 4.5: Calculate epsilon_star-level global consistency/inconsistency probabilities
  //    The probability for epsilon-star global inconsistency is calculated as the sum of the PMPs for
  //    all models which do not allow regional treatment effects to differ by more than epsilon_star with
  //    (1 - beta_star) probability. If Pr(|gamma_i - gamma_j| > epsilon_star | D) >= 1 - beta_star
  //    for i != j, then we consider all models in which regions i and j differ. Identify these models
  //    for all pairwise comparisons, and then add up the PMPs for all of these models.
  
  
  // Determine if regional treatment effects are greater than a clinically meaninful difference
  // of one another. If so, determine which models allow these regions to differ and
  // store results in which_mods_keep (1 if model should be considered, 0 otherwise)
  arma::mat which_mods_keep(numMods, epsilon_star.size(), fill::zeros);
  for(int q = 0; q < epsilon_star.size(); q++){
    
    arma::mat pwic_mat_q = bma_prws_inconsis[q];    // Extract matrix of pwc probs. for q^th epsilon_star value
    for(int l = 0; l < numMods; l++){
      
      int i_loop_continue = 1;     // Continue with i for-loop if equal to 1; break if 0
      for(int i = 0; i < S; i++){
        
        for(int j = i+1; j < S; j++){
          
          // If Pr(|gamma_i - gamma_j| > epsilon_star | D) < 1 - beta_star, consider Model l in
          // calculation for global inconsistency only if model forces regions i and j to differ
          
          if(pwic_mat_q(i,j) >= 1 - beta_star){
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
  
  // For q^th value of epsilong_star, calculate epsilon_star-level global inconsistency probability
  // for models which require any two regions to differ by more than epsilon_star, with
  // probability (1 - beta_star). The epsilon_star-level global consistency probability is
  // defined as 1 - epsilon_star-level global inconsistency probability.
  arma::vec glob_inconsis_prob = which_mods_keep.t() * pmp;
  arma::vec glob_consis_prob = 1.0 - glob_inconsis_prob;
  
  
  // ## Return list with the following:
  //       (1) vector of posterior model probabilities
  //       (2) matrix of summary stats for regional treatment effects using BMA 
  //       (3) vector of summary stats for global treatment effect using BMA
  //       (4) matrix of summary stats for regional intercepts using BMA 
  //       (5) matrix of summary stats for covariate effects (if there are covariates) using BMA
  //       (6) matrix of leave-one-out local consistency probabilities
  //       (7) list of matrices with pairwise consistency probabilities
  //       (8) list of matrices with pairwise inconsistency probabilities
  //       (9) matrix with epsilon_star global consistency probability
  //       (10) matrix with epsilon_star global inconsistency probability
  Rcpp::List single_sim_results;
  single_sim_results["pmp"] = pmp;
  single_sim_results["RegionalTrtmtStats"] = bma_vals_rte;
  single_sim_results["GlobalTrtmtStats"] = bma_vals_glob;
  single_sim_results["RegionalIntStats"] = bma_vals_ri;
  if(p_cov > 0){
    single_sim_results["CovariateStats"] = bma_vals_cov;
  }
  single_sim_results["LOOLocalConsistency"] = bma_LOO_loc_cons;
  single_sim_results["PairwiseConsistency"] = bma_prws_consis;
  single_sim_results["PairwiseInconsistency"] = bma_prws_inconsis;
  single_sim_results["GlobalConsistencyProb"] = glob_consis_prob;
  single_sim_results["GlobalInconsistencyProb"] = glob_inconsis_prob;
  
  return single_sim_results;
  
}



// ### Function to simulate one dataset
// # User input must include the following:
//     N: sample size
//     S: number of regions
//     reg_names: S-dimensional character vector with region names
//     reg_allctn: S-dimensional vector of total sample size allocation to regions (must sum to 1)
//     trtmt_allctn: 2-dimensional vector of treatment allocations within each region
//                   (treatment first, then control; must sum to 1)
//     trtmt_means: S-dimensional vector of regional treatment means
//     cntrl_means: S-dimensional vector of regional control means
//     sd_Y: common standard deviation of the response Y for both treatment groups in all regions
//     numBin: number of binary covariates (each covariate is iid Bernoulli)
//     numCon: number of continuous covariates (each covariate is iid normal)
//     cov_props: numBin-dimensional vector of success proportions for binary covariates
//     cov_means: numCon-dimensional vector of means for continuous covariates
//     cov_sds: numCon-dimensional vector of standard deviations for continuous covariates
// [[Rcpp::export]]
Rcpp::List generate_data(int N, int S, Rcpp::CharacterVector reg_names, Rcpp::NumericVector & reg_allctn,
                         Rcpp::NumericVector & trtmt_allctn, Rcpp::NumericVector & trtmt_means,
                         Rcpp::NumericVector & cntrl_means, double sd_Y, int numBin, int numCon,
                         Rcpp::NumericVector & cov_props, Rcpp::NumericVector & cov_means,
                         Rcpp::NumericVector & cov_sds){
  
  // Check function inputs
  if(N <= 0)  Rcpp::stop("N must be greater than 0");
  if(S <= 0)  Rcpp::stop("S must be greater than 0");
  Rcpp::CharacterVector sorted_regions = clone(reg_names);
  sorted_regions = sorted_regions.sort(false);
  for(int s = 0; s < S; s++){
    if(reg_names[s] != sorted_regions[s]){
      Rcpp::stop("Must list region names alphabetically. Verify that reg_allctn, trtmt_means, cntrl_means correspond to correct order");
    }
  }
  if(reg_names.size() != S)  Rcpp::stop("The length of reg_names must equal S");
  if(reg_allctn.size() != S)  Rcpp::stop("The length of reg_allctn must equal S");
  if(sum(reg_allctn) != 1)  Rcpp::stop("reg_allctn must sum to 1");
  if(trtmt_allctn.size() != 2)  Rcpp::stop("The length of trtmt_allctn must equal 2");
  if(sum(trtmt_allctn) != 1)  Rcpp::stop("trtmt_allctn must sum to 1");
  if(trtmt_means.size() != S)  Rcpp::stop("The length of trtmt_means must equal S");
  if(cntrl_means.size() != S)  Rcpp::stop("The length of cntrl_means must equal S");
  if(sd_Y <= 0)  Rcpp::stop("sd_Y must be greater than 0");
  if(numBin < 0)  Rcpp::stop("numBin must be non-negative");
  if(numCon < 0)  Rcpp::stop("numCon must be non-negative");
  if(numBin > 0){
    if(cov_props.size() != numBin)  Rcpp::stop("The length of cov_props must equal numBin");
    for(int p = 0; p < cov_props.size(); p++){
      if(cov_props[p] < 0 || cov_props[p] > 1){
        Rcpp::stop("All entries of cov_props must be between 0 and 1, inclusively");
      }
    }
  }
  if(numCon > 0){
    if(cov_means.size() != numCon)  Rcpp::stop("The length of cov_means must equal numCon");
    if(cov_sds.size() != numCon)  Rcpp::stop("The length of cov_sds must equal numCon");
    for(int p = 0; p < cov_sds.size(); p++){
      if(cov_sds[p] <= 0)  Rcpp::stop("All entries of cov_sds must be greater than 0");
    }
  }
  
  
  // ## Determine regional sample sizes n_i and ensure they sum to N
  Rcpp::NumericVector n_i_approx = ceiling(N * reg_allctn);      // approximate regional sample sizes
  
  // Reduce regional sample sizes by 1 beginning with first region until sample sizes sum to N
  int i = 0;
  while( sum(n_i_approx) > N ){
    n_i_approx[i] = n_i_approx[i] - 1;
    i++;
  }
  Rcpp::NumericVector n_i = n_i_approx;       // S-dimensional vector of regional sample sizes (sum to N)
  
  
  // ## Determine region-by-treatment-group sample sizes n_ij and ensure they sum to n_i for each region
  Rcpp::NumericMatrix n_ij_approx(S,2);     // S x 2 matrix with region-by-trtmt_group sample sizes
  //    Rows correspond to region
  //    First column: treatment group; second column: control group
  for(int i = 0; i < S; i++){
    
    n_ij_approx.row(i) = floor(n_i[i] * trtmt_allctn);
    
    // Increase treatment group sample size by 1 if two treatment group sizes do not add up to n_i
    if( sum(n_ij_approx.row(i)) < n_i[i] ){
      n_ij_approx(i,0) = n_ij_approx(i,0) + 1;
    }
    
  }
  Rcpp::NumericMatrix n_ij = n_ij_approx;    // S x 2 matrix with region-by-trtmt_group sample sizes (rows sum to n_i)
  
  
  // ## Simulate response vector Y and record region labels and treatment indicators
  Rcpp::NumericVector Y(N);                 // N-dimensional vector to store normal outcomes
  Rcpp::CharacterVector region_labs(N);     // N-dimensional vector to store region labels
  Rcpp::IntegerVector trtmt_indctr(N);      // N-dimensional vector to store treatment indicators
  
  int n_ij_cumsum = 0;
  for(int i = 0; i < S; i++){
    
    for(int j = 0; j < 2; j++){
      
      for(int k = 0; k < n_ij(i,j); k++){
        
        if(j == 0){
          Y[n_ij_cumsum + k] = R::rnorm( trtmt_means[i], sd_Y );
          trtmt_indctr[n_ij_cumsum + k] = 1;
        }
        else{
          Y[n_ij_cumsum + k] = R::rnorm( cntrl_means[i], sd_Y );
          trtmt_indctr[n_ij_cumsum + k] = 0;
        }
        region_labs[n_ij_cumsum + k] = reg_names[i];
        
      }
      
      // Update the number of subjects in preceding regions and treatment groups
      n_ij_cumsum = n_ij_cumsum + n_ij(i,j);
      
    }
    
  }
  
  
  // ## Simulate numBin binary covariates and numCon continuous covariates
  int numCovs = numBin + numCon;
  Rcpp::NumericMatrix X(N, numCovs);
  for(int k = 0; k < N; k++){
    
    for(int p = 0; p < numBin; p++){
      // Draw p^th binary covariate from Bernoulli(cov_props[p])
      X(k,p) = R::rbinom(1, cov_props[p]);
    }
    
    for(int p = 0; p < numCon; p++){
      // Draw p^th continuous covariate from N(cov_means[p], cov_sds[p]^2)
      X(k,p+numBin) = R::rnorm(cov_means[p], cov_sds[p]);
    }
  }
  
  // Add column names for covariate matrix
  Rcpp::CharacterVector X_col_names(numBin + numCon);
  for(int p = 1; p <= numBin + numCon; p++){
    string p_char = to_string(p);
    string cov_string = "Cov";
    string col_name = cov_string + p_char;
    X_col_names[p-1] = col_name;
  }
  Rcpp::colnames(X) = X_col_names;
  
  
  // ## Return list with the following:
  //       (1) Y: N-dimensional response vector
  //       (2) region_labs: N-dimensional character vector of region labels (S unique values)
  //       (3) trtmt_indctr: N-dimensional vector of treatment indicators (1: treatment; 0: control)
  //       (4) X: N x p matrix of covariates
  Rcpp::List data_list;
  data_list["Y"] = Y;
  data_list["RegionLabels"] = region_labs;
  data_list["TreatmentIndicator"] = trtmt_indctr;
  data_list["Covariates"] = X;
  
  return data_list;
  
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
//     alpha: tuning parameter, where p(M_l) is proportional to T_l^alpha for the l^th model
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
    mod_priors[l] = pow(num_distinct, alpha);
  }
  mod_priors = mod_priors / sum(mod_priors);
  
  return mod_priors;
  
}

