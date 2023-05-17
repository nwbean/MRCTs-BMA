GUIDE TO R SCRIPTS:


R scripts used in simulation process, including scripts to create CSV file with
input values, run simulations, and compile results.


bma-joint-model-function.R - main function that runs the proposed joint BMA approach algorithm
     using Laplace approximations (used for all simulations except "SA-model-priors")



bma-joint-model-function-model-priors.R - main function that runs the proposed joint BMA
     approach algorithm using Laplace approximations (used for simulation "SA-model-priors")



bma-survival-function.R - main function that runs the survival-only BMA approach algorithm
     using Laplace approximation



compile-results-alt-null-half.R - plots bar graphs to compare the BMA joint modeling approach to
     the CPHMs and the survival-only BMA approach with respect to the global rejection rates, MSE
     of region-specific treatment effects, and true/false positive rates for the region-specific
     treatment effects on the time-to-event endpoint for the following simulation studies:
        - equal-samp-alt-half
        - equal-samp-null-half



compile-results-rr-mse.R - plots bar graphs to compare the BMA joint modeling approach to the
     CPHMs and the survival-only BMA approach with respect to the global rejection rates, MSE
     of region-specific treatment effects, and true/false positive rates for the region-specific
     treatment effects on the time-to-event endpoint for the following simulation studies:
        - equal-samp-alpha0
        - equal-samp-alpha0p15
        - equal-samp-alpha0p5
        - equal-samp-alpha1
        - original-alpha0
        - original-alpha0p15
        - original-alpha0p5
        - original-alpha1



compile-results-SA-model-priors.R - plots bar graphs to compare the BMA joint modeling approach
     to the CPHMs and the survival-only BMA approach with respect to the global rejection rates,
     MSE of region-specific treatment effects, and true/false positive rates for the region-specific
     treatment effects on the time-to-event endpoint for the following simulation study:
        - SA-model-priors



compile-results-SA-Q.R - plots bar graphs to compare the BMA joint modeling approach to the
     CPHMs and the survival-only BMA approach with respect to the global rejection rates, MSE
     of region-specific treatment effects, and true/false positive rates for the region-specific
     treatment effects on the time-to-event endpoint for the following simulation studies:
        - equal-samp-Q5
        - equal-samp-Q12



compile-results-survival-models-only.R - plots bar graphs to compare the simulated joint data
     versus simulated survival-only data for the CPHMs and the survival-only BMA approach with
     respect to the global rejection rates, MSE of region-specific treatment effects, and
     true/false positive rates for the region-specific treatment effects on the time-to-event
     endpoint for the following simulation study:
        - survival-models-only



compile-results-T0-comparison.R - plots bar graphs to compare the BMA joint modeling approach
     for different values of the maximum distinct treatment effects per submodel (T0) with
     respect to the global rejection rates, MSE of region-specific treatment effects, and
     true/false positive rates for the region-specific treatment effects on the time-to-event
     endpoint for the following simulation study:
        - original-alpha0p5
        - original-T0equal3
        - original-T0equal4



compile-results-violations.R - summarize bias, MSE, and posterior means for scenarios that
     violate model assumptions and compare these operating characteristics to the scenario
     in which all underlying model assumptions are satisfied. The results are summarized for
     the following simulation study:
        - equal-samp-alpha0p5
        - equal-samp-log-gamma-re
        - equal-samp-nonPH
        - equal-samp-inf-cens



data-simulation-functions.R - functions to simulate both survival and longitudinal data



data-simulation-functions-inf-cens.R - functions to simulate both survival and longitudinal
     data for simulation study with informative censoring



data-simulation-functions-nonPH.R - functions to simulate both survival and longitudinal
     data for simulation study with non-proportional hazards


full-simulation-function.R - function that runs entire simulation for the following studies:
   - equal-samp-alpha0
   - equal-samp-alpha0p15
   - equal-samp-alpha0p5
   - equal-samp-alpha1
   - equal-samp-alt-half
   - equal-samp-null-half
   - equal-samp-Q5
   - equal-samp-Q12
   - original-alpha0
   - original-alpha0p15
   - original-alpha0p5
   - original-alpha1



full-simulation-function-inf-cens.R - function that runs entire simulation for the following
     study:
   - equal-samp-inf-cens



full-simulation-function-log-gamma-re.R - function that runs entire simulation for the
     following study:
   - equal-samp-log-gamma-re



full-simulation-function-nonPH.R - function that runs entire simulation for the following
     study:
   - equal-samp-nonPH



full-simulation-function-SA-model-priors.R - function that runs entire simulation for the study
     "SA-model-priors"



full-simulation-function-survival-models-only.R - function that runs entire simulation for the
     study "survival-models-only"



proportional-hazards-models.R - function that runs both CPHMs



run-simulation-equal-samp-alpha0.R - runs simulation "equal-samp-alpha0" with the following
     details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.0



run-simulation-equal-samp-alpha0p15.R - runs simulation "equal-samp-alpha0p15" with the following
     details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.15



run-simulation-equal-samp-alpha0p5.R - runs simulation "equal-samp-alpha0p5" with the following
     details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



run-simulation-equal-samp-alpha1.R - runs simulation "equal-samp-alpha1" with the following
     details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 1.0



run-simulation-equal-samp-alt-half.R - runs simulation "equal-samp-alt-half" with the following
     details:
   - alternative sample sizes half of null
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



run-simulation-equal-samp-inf-cens.R - runs simulation "equal-samp-inf-cens" with the following
     details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5
   - informative censoring mechanism



run-simulation-equal-samp-log-gamma-re.R - runs simulation "equal-samp-log-gamma-re" with the
     following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5
   - subject-specific random intercepts drawn from log-gamma distribution



run-simulation-equal-samp-nonPH.R - runs simulation "equal-samp-nonPH" with the following
     details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5
   - non-proportional hazards



run-simulation-equal-samp-null-half.R - runs simulation "equal-samp-null-half" with the following
     details:
   - null sample sizes half of alternative
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



run-simulation-equal-samp-Q5.R - runs simulation "equal-samp-Q5" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 5
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



run-simulation-equal-samp-Q12.R - runs simulation "equal-samp-Q12" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 12
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



run-simulation-original-alpha0.R - runs simulation "original-alpha0" with the following details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.0



run-simulation-original-alpha0p15.R - runs simulation "original-alpha0p15" with the following
     details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.15



run-simulation-original-alpha0p5.R - runs simulation "original-alpha0p5" with the following
     details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5



run-simulation-original-alpha1.R - runs simulation "original-alpha1" with the following details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 1.0



run-simulation-original-T0equal3.R - runs simulation "original-T0equal3" with the following
     details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5
   - T0: 3



run-simulation-original-T0equal4.R - runs simulation "original-T0equal4" with the following
     details:
   - sample sizes of 711, 3296, 2847, 2486 for Regions 1-4
   - case 1: treatment hazard ratios of .868 for all regions
   - case 2: treatment hazard ratios of .62, .82, 1.01, .83 for Regions 1-4
   - N = 9340
   - Q: 8
   - a_X: 0
   - a_Y: 0
   - alpha: 0.5
   - T0: 4



run-simulation-SA-model-priors.R - runs simulation "SA-model-priors" with the following details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8
   - a_X: values in {-1, 0, 1}
   - a_Y: values in {-1, 0, 1}
   - alpha: 0.5



run-simulation-survival-models-only.R - runs simulation "survival-models-only" with the following
     details:
   - equal regional sample sizes
   - treatment hazard ratio of .868 for alternative regions
   - N = 9340
   - Q: 8 (survival-only BMA approach)
   - a_0: 0 (survival-only BMA approach)
   - alpha: 0.5



setup-design-inputs.R - creates CSV files with input values for simulations and saves
     file to the "cluster-input" folder

