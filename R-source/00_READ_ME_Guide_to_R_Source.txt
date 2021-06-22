GUIDE TO R SOURCE:


CPP files with Rcpp code and TXT files with JAGS scripts



bma_functions_rcpp_code.cpp
   - Rcpp code with main functions to generate data and run simulations. Returns global
     treatment effect estimates, region-specific treatment effects, rejection rates, bias
     and MSE (region-specific intercepts and region-specific treatment effects), PMPs,
     local constancy probabilities (Japanese MHLW approach and proposed leave-one-out absolute
     difference approach), and epsilon-level pairwise consistency probabilities



jagsBHM.txt
   - JAGS script for Bayesian hierarchical model



jagsFELM_glob.txt
   - JAGS script for first fixed-effects linear model used to calculate global treatment
     effect


jagsFELM_reg.txt
   - JAGS script for second fixed-effects linear model used to calculate region-specific
     treatment effects

