GUIDE TO SIMULATION CODE (MORE-DETAILED GUIDES PROVIDED IN EACH FOLDER):



cluster-err
   - Error files from cluster when running simulation studies (files deleted after
     problems are corrected and simulations rerun - this folder is now empty)



cluster-input
   - CSV files with input values used for simulation studies



cluster-log
   - Folders for each simulation study containing R logs with simulation details and
     output



cluster-out
   - Out files from cluster when running simulation studies (files deleted - this folder
     is now empty)



cluster-results
   - Folders for each simulation study containing results for the global
     treatment effect estimates (survival and longitudinal), region-specific
     treatment effects (survival and longitudinal), rejection rates (survival),
     bias and MSE of region-specific treatment effects (survival), the association
     parameter, PMPs, local constancy probabilities (two variations of the Japanese
     MHLW approach and the proposed leave-one-out absolute difference approach)
     (survival), and epsilon-level pairwise consistency probabilities (survival)
   - NOTE: results cannot be uploaded to GitHub due to size. All results can be reproduced
     using the files provided.


cluster-scripts
   - Batch files used to begin simulations on cluster


R-scripts
   - R scripts used in simulation process, including scripts to create CSV file with
     input values, run simulations, and compile results


R-source
   - CPP file with Rcpp code
