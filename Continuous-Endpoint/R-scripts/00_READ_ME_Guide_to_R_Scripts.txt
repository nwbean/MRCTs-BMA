GUIDE TO R SCRIPTS:


R scripts used in simulation process, including scripts to create CSV file with
input values, run simulations, and compile results.



BMA_single_dataset.R - runs proposed BMA approach for a single simulated dataset. Users can
     replace simulated dataset with real data and change values of elicited priors



compile_global_consistency_results.R - compiles results for assessing global consistency
     with proposed epsilon-level global consistency probability from simulations, and creates
     plots featured in main paper for differing values of epsilon



compile_results.R - compiles all results for all simulation versions and provides global
     treatment effect estimates, region-specific treatment effects, rejection rates, bias
     and MSE (region-specific intercepts and region-specific treatment effects), PMPs,
     local constancy probabilities (Japanese MHLW approach and proposed leave-one-out absolute
     difference approach), and epsilon-level pairwise consistency probabilities. Also plots bar
     graphs to compare the BMA approach to the FELM and BHM with respect to the global rejection
     rate, the true positive rate, the false positive rate, and the MSE. To compile results for
     global consistency, see "compile_global_consistency_results.R" script



full_simulation_function.R - function that runs entire simulation and returns the
     global treatment effect estimates, region-specific treatment effects, rejection
     rates, bias and MSE (region-specific intercepts and region-specific treatment
     effects), PMPs, local constancy probabilities (MHLW and proposed leave-one-out
     absolute difference), and epsilon-level pairwise consistency probabilities.



run-simulation-alt-half.R - runs simulation with the following details:
   - alternative sample sizes half of null
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-diff-effects-1.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - underlying region-specific treatment effects: .017, .026, .034, .043, .051
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-diff-effects-2.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - underlying region-specific treatment effects: .017, .017, .017, .034, .034
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-diff-effects-3.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - underlying region-specific treatment effects: .017, .034, .034, .051, .051
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-n1508.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-n1508-alpha2.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 2
   - beta.star: .5



run-simulation-n1508-alpha4.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 4
   - beta.star: .5



run-simulation-n1508-alpha10.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 10
   - beta.star: .5



run-simulation-n1508-alphaNeg2.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: -2
   - beta.star: .5



run-simulation-n1508-alphaNeg4.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: -4
   - beta.star: .5



run-simulation-n1508-alphaNeg10.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: -10
   - beta.star: .5



run-simulation-n1508-beta20.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .2



run-simulation-n754.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 754 (half of N = 1508)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-n754-beta20.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 754 (half of N = 1508)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .2



run-simulation-n7650.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 7650 (region sample sizes calculated to detect 90% power)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-n7650-beta20.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 7650 (region sample sizes calculated to detect 90% power)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .2



run-simulation-SA1.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .082 for regional intercepts and .034 for regional treatment effects
   - Sigma0: diagonals of (10 x .082)^2 for regional intercepts and
     (10 x .034)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-SA2.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .082 for regional intercepts and .017 for regional treatment effects
   - Sigma0: diagonals of (10 x .082)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-SA3.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .082 for regional intercepts and .068 for regional treatment effects
   - Sigma0: diagonals of (10 x .082)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-SA4.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .041 for regional intercepts and .017 for regional treatment effects
   - Sigma0: diagonals of (10 x .041)^2 for regional intercepts and
     (10 x .017)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-SA5.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .041 for regional intercepts and .034 for regional treatment effects
   - Sigma0: diagonals of (10 x .041)^2 for regional intercepts and
     (10 x .034)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-SA6.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .041 for regional intercepts and .068 for regional treatment effects
   - Sigma0: diagonals of (10 x .041)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-SA7.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .164 for regional intercepts and .017 for regional treatment effects
   - Sigma0: diagonals of (10 x .164)^2 for regional intercepts and
     (10 x .017)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-SA8.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .164 for regional intercepts and .034 for regional treatment effects
   - Sigma0: diagonals of (10 x .164)^2 for regional intercepts and
     (10 x .034)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-SA9.R - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .164 for regional intercepts and .068 for regional treatment effects
   - Sigma0: diagonals of (10 x .164)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



setup_design_inputs.R - creates CSV file with input values for simulations and saves
     file to the "cluster-input" folder

