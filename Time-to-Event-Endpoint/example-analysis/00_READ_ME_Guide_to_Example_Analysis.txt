GUIDE TO EXAMPLE ANALYSIS:


Code and simulated data that can be used to run the BMA algorithm (serves as an example).


bma_functions_tte_endpoint.cpp - Rcpp code with functions to calculate different forms
   of consistency



bma-asymptotic-function.R - main function that runs proposed BMA approach algorithm using
   Laplace approximation



BMA-single-dataset-analysis-example.R - runs proposed BMA approach using example_data.csv.
   Additional instructions are included as comments in the file.



BMA-single-dataset-analysis-template.R - runs proposed BMA approach for any single dataset.
   Users can change values of elicited priors by following the instructions that are included
   as comments in the file.



example_data.csv - simulated dataset similar to the datasets created for the simulation
   study "n9340-original"
   - N = 9340
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - underlying treatment hazard ratios of .868 for all regions

