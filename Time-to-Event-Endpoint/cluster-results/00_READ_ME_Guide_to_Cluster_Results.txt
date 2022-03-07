GUIDE TO CLUSTER RESULTS:


Folders for the following simulation studies containing results for the global
     treatment effect estimates, region-specific treatment effects, rejection rates, bias
     and MSE (region-specific intercepts and region-specific treatment effects), PMPs,
     local constancy probabilities (two variations of the Japanese MHLW approach and the
     proposed leave-one-out absolute difference approach), and epsilon-level pairwise
     consistency probabilities



n9340-alt-half
   - alternative sample sizes half of null
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



n9340-equal-samp
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



n9340-equal-samp-beta0_2
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .2



n9340-equal-samp-beta0_8
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .8



n9340-equal-samp-SA-alpha
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: values in {-5, -2, -1, -.5, 0, .5, 1, 2, 5}
   - beta.star: .5



n9340-equal-samp-SA-K
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: values in {4, 8, 12, 16}
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



n9340-equal-samp-SA-mu0-Sig0
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: values in {log(.7), log(1.05), log(1.3), log(1.5)} for all elements
   - Sig0: diagonals of (10 x mu0)^2
   - alpha0: 0
   - beta.star: .5



n9340-nonexponential-BH
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - baseline hazards of .02, .035, .055 over intervals (0, 2], (2, 3.75], (3.75, inf) 
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



n9340-null-half
   - null sample sizes half of alternative
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5



n9340-original
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - K: 8
   - mu0: log(1.3) for all elements
   - Sig0: diagonals of (10 x log(1.3))^2
   - alpha0: 0
   - beta.star: .5
