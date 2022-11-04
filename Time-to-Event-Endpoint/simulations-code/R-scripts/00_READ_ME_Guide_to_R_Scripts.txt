GUIDE TO R SCRIPTS:


R scripts used in simulation process, including scripts to create CSV file with
input values, run simulations, and compile results.


bma-asymptotic-function.R - main function that runs proposed BMA approach algorithm using
     Laplace approximation



compile-global-consistency-results.R - compiles results for assessing global consistency
     with proposed epsilon-level global consistency probability from simulations, and creates
     plots for differing values of epsilon



compile-results-BHM-priors.R - compiles results for global treatment effect estimates, region-
     specific treatment effects, rejection rates, and bias and MSE of region-specific treatment
     effects for the following simulation studies:
        - n9340-equal-samp-BHM-priors
        - n9340-original-BHM-priors
     Also plots bar graphs to compare the BMA approach to the CPHMs and three BHMs (different
     priors on the hierarchical precision) with respect to the global rejection rate, the MSE
     for region-specific treatment effects, the true positive rate, and the false positive rate.



compile-results-rr-mse.R - compiles results for global treatment effect estimates, region-
     specific treatment effects, rejection rates, bias and MSE of region-specific treatment
     effects, PMPs, local constancy probabilities (two variations of the Japanese MHLW approach
     and the proposed leave-one-out absolute difference approach), and epsilon-level pairwise
     consistency probabilities for the following simulation studies:
        - n9340-equal-samp
        - n9340-original
        - n9340-null-half
        - n9340-alt-half
        - n9340-nonexponential-BH
     Also plots bar graphs to compare the BMA approach to the CPHMs and the BHM with respect to
     the global rejection rate, the MSE for region-specific treatment effects, the true positive
     rate, and the false positive rate.



compile-results-SA.R - compiles results for rejection rates, bias, and MSE (region-specific
     treatment effects) for the following simulation studies:
        - n9340-equal-samp-SA-alpha
        - n9340-equal-samp-SA-mu0-Sig0
        - n9340-equal-samp-SA-K
     Also plots bar graphs to compare the BMA approach to the CPHMs and the BHM with respect to
     the global rejection rate, the MSE for region-specific treatment effects, the true positive
     rate, and the false positive rate.



full-simulation-function.R - function that runs entire simulation for the following studies:
   - n9340-equal-samp
   - n9340-equal-samp-beta0_2
   - n9340-equal-samp-beta0_8
   - n9340-original
   - n9340-null-half
   - n9340-alt-half



full-simulation-function-BHM-priors.R - function that runs entire simulation for the
following studies:
   - n9340-equal-samp-BHM-priors
   - n9340-original-BHM-priors



full-simulation-function-nonexponential-BH.R - function that runs entire simulation for the study
     "n9340-nonexponential-BH"



full-simulation-function-SA-alpha0.R - function that runs entire simulation for the study
     "n9340-equal-samp-SA-alpha0"



full-simulation-function-SA-K.R - function that runs entire simulation for the study
     "n9340-equal-samp-SA-K"



full-simulation-function-SA-mu0-Sig0.R - function that runs entire simulation for the study
     "n9340-equal-samp-SA-mu0-Sig0"



install_packages.R - script that installs packages contained in folder "R-packages"



proportional-hazards-models.R - function that runs both CPHMs



run-simulation-n9340-alt-half.R - runs simulation with the following details:
   - alternative sample sizes half of null
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



run-simulation-n9340-equal-samp.R - runs simulation with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



run-simulation-n9340-equal-samp-beta0_2.R - runs simulation with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .2



run-simulation-n9340-equal-samp-beta0_8.R - runs simulation with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .8



run-simulation-n9340-equal-samp-BHM-priors.R - runs simulation with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5
   - compares BHMs with different priors on hierarchical precision



run-simulation-n9340-equal-samp-SA-alpha.R - runs simulation with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: values in {-5, -2, -1, -.5, 0, .5, 1, 2, 5}
   - beta.star: .5



run-simulation-n9340-equal-samp-SA-K.R - runs simulation with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: values in {4, 8, 12, 16}
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



run-simulation-n9340-equal-samp-SA-mu0-Sig0.R - runs simulation with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: values in {log(.7), log(1.05), log(1.3), log(1.5)} for all elements
   - Sig0: diagonals of (10 x mu0)^2
   - alpha0: 0
   - beta.star: .5



run-simulation-n9340-nonexponential-BH.R - runs simulation with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - baseline hazards of .02, .035, .055 over intervals (0, 2], (2, 3.75], (3.75, inf) 
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



run-simulation-n9340-null-half.R - runs simulation with the following details:
   - null sample sizes half of alternative
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



run-simulation-n9340-original.R - runs simulation with the following details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



run-simulation-n9340-original-BHM-priors.R - runs simulation with the following details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5
   - compares BHMs with different priors on hierarchical precision



setup-design-inputs.R - creates CSV files with input values for simulations and saves
     file to the "cluster-input" folder

