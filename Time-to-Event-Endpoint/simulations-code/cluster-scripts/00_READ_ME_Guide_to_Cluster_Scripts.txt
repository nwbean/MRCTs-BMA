GUIDE TO CLUSTER SCRIPTS:


Batch files used to begin simulations on cluster



batch-n9340-alt-half.sh - runs simulation "n9340-alt-half" with the following details:
   - alternative sample sizes half of null
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



batch-n9340-equal-samp.sh - runs simulation "n9340-equal-samp" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



batch-n9340-equal-samp-beta0_2.sh - runs simulation "n9340-equal-samp-beta0_2"
with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .2



batch-n9340-equal-samp-beta0_8.sh - runs simulation "n9340-equal-samp-beta0_8"
with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .8



batch-n9340-equal-samp-BHM-priors.sh - runs simulation "n9340-equal-samp-BHM-priors"
with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5
   - compares BHMs with different priors on hierarchical precision



batch-n9340-equal-samp-SA-alpha.sh - runs simulation "n9340-equal-samp-SA-alpha"
with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: values in {-5, -2, -1, -.5, 0, .5, 1, 2, 5}
   - beta.star: .5



batch-n9340-equal-samp-SA-K.sh - runs simulation "n9340-equal-samp-SA-K" with the
following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: values in {4, 8, 12, 16}
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



batch-n9340-equal-samp-SA-mu0-Sig0.sh - runs simulation "n9340-equal-samp-SA-mu0-Sig0"
with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: values in {log(.7), log(1.05), log(1.3), log(1.5)} for all elements
   - Sig0: diagonals of (10 x mu0)^2
   - alpha0: 0
   - beta.star: .5



batch-n9340-nonexponential-BH.sh - runs simulation "n9340-nonexponential-BH" with the
following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - baseline hazards of .02, .035, .055 over intervals (0, 2], (2, 3.75], (3.75, inf) 
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



batch-n9340-null-half.sh - runs simulation "n9340-null-half" with the following details:
   - null sample sizes half of alternative
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



batch-n9340-original.sh - runs simulation "n9340-original" with the following details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



batch-n9340-original-BHM-priors.sh - runs simulation "n9340-original-BHM-priors"
with the following details:
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
