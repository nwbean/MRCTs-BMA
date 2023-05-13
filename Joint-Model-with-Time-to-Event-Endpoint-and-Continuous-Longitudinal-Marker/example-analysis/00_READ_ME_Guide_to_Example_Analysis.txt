GUIDE TO EXAMPLE ANALYSIS:


Code and simulated data that can be used to run the BMA algorithm (serves as an example).


bma-functions-joint-model.cpp - Rcpp code with functions to calculate different forms
   of consistency on time-to-event endpoint



bma-joint-model-function.R - main function that runs proposed BMA approach algorithm using
   Laplace approximations



BMA-single-dataset-analysis-example.R - runs proposed BMA approach using simulated datasets
   example_survival_data.csv and example_longitudinal_data. Additional instructions are
   included as comments in the file.



BMA-single-dataset-analysis-template.R - runs proposed BMA approach for any single data
   example (separate datasets for survival and longitudinal data--file contains instructions
   for how these datasets should be set up). Users can change values of elicited priors and
   model specifications by following the instructions that are included as comments in the
   file.



example_survival_data.csv - simulated survival dataset similar to the datasets created for
   the first scenario in simulation study "original-alpha0p5"
   - N = 9340
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - underlying treatment hazard ratios of .868 for all regions



example_longitudinal_data.csv - simulated longitudinal dataset corresponding to the survival
   dataset saved as "example_survival_dataset.csv"

