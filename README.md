# MRCTs-BMA

This repository contains all files and code used in the simulation studies that evaluate the Bayesian model averaging (BMA) approach for multi-regional clinical trials (MRCTs) with either a continuous or a time-to-event endpoint.


Continuous-Endpoint
   - Contains all files and code used in the simulation studies that evaluate the BMA approach
     for MRCTs with a continuous endpoint, as presented in the paper "Bayesian multi-regional
     clinical trials using model averaging" (DOI: kxab027) in _Biostatistics_. These files can
     be used to reproduce all simulation results.


Time-to-Event-Endpoint
   - Contains all files and code used in the simulation studies that evaluate the BMA approach
     for MRCTs with a time-to-event endpoint.


Joint Model with Time-to-Event-Endpoint and Continuous Longitudinal Marker
   - Contains all files and code used in the simulation studies that evaluate the BMA approach
     for MRCTs with both a time-to-event endpoint and a continuous longitudinal marker.


GUIDE TO SIMULATION CODE IN FOLDERS "CONTINUOUS-ENDPOINT", "TIME-TO-EVENT-ENDPOINT", and "JOINT MODEL WITH TIME-TO-EVENT-ENDPOINT AND CONTINUOUS LONGITUDINAL MARKER" (MORE-DETAILED GUIDES PROVIDED IN EACH FOLDER AND EACH SUBFOLDER):


In addition to the following folders found in this GitHub repository, blank folders should also be created titled "cluster-err" and "cluster-out". Error files and out files from the cluster will be saved here when running the simulations.


cluster-input
   - CSV files with input values used for simulation studies


cluster-log
   - Folders for each simulation study containing R logs with simulation details and
     output


cluster-results
   - Folders for each simulation study containing results for the global
     treatment effect estimates, region-specific treatment effects, rejection rates, bias
     and MSE (region-specific intercepts and region-specific treatment effects), PMPs,
     local constancy probabilities (Japanese MHLW approach and proposed leave-one-out absolute
     difference approach), and epsilon-level pairwise consistency probabilities
   - NOTE: results cannot be uploaded to GitHub due to size. All results can be reproduced
     using the files provided.


cluster-scripts
   - Batch files used to begin simulations on cluster


R-scripts
   - R scripts used in simulation process, including scripts to create CSV file with
     input values, run simulations, compile results, and assess consistency


R-source
   - CPP files with Rcpp code and TXT files with JAGS scripts
