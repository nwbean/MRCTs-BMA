GUIDE TO CLUSTER RESULTS:


Folders for the following simulation studies containing results for the global
     treatment effect estimates (survival and longitudinal), region-specific
     treatment effects (survival and longitudinal), rejection rates (survival),
     bias and MSE of region-specific treatment effects (survival), the association
     parameter, PMPs, local constancy probabilities (two variations of the Japanese
     MHLW approach and the proposed leave-one-out absolute difference approach)
     (survival), and epsilon-level pairwise consistency probabilities (survival)



batch-equal-samp-alpha0.sh - runs simulation "equal-samp-alpha0" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.0



batch-equal-samp-alpha0p15.sh - runs simulation "equal-samp-alpha0p15" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.15



batch-equal-samp-alpha0p5.sh - runs simulation "equal-samp-alpha0p5" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



batch-equal-samp-alpha1.sh - runs simulation "equal-samp-alpha1" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 1.0



batch-equal-samp-alt-half.sh - runs simulation "equal-samp-alt-half" with the following details:
   - alternative sample sizes half of null
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



batch-equal-samp-null-half.sh - runs simulation "equal-samp-null-half" with the following details:
   - null sample sizes half of alternative
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



batch-equal-samp-Q5.sh - runs simulation "equal-samp-Q5" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 5
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



batch-equal-samp-Q12.sh - runs simulation "equal-samp-Q12" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 12
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



batch-original-alpha0.sh - runs simulation "original-alpha0" with the following details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.0



batch-original-alpha0p15.sh - runs simulation "original-alpha0p15" with the following details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.15



batch-original-alpha0p5.sh - runs simulation "original-alpha0p5" with the following details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



batch-original-alpha1.sh - runs simulation "original-alpha1" with the following details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 1.0



batch-SA-model-priors.sh - runs simulation "SA-model-priors" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: values in {-1, 0, 1}
   - a_Y: values in {-1, 0, 1}
   - alpha: 0.5


batch-survival-models-only.sh - runs simulation "survival-models-only" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8 (survival-only BMA approach)
   - a_0: 0 (survival-only BMA approach)
   - alpha: 0.5
