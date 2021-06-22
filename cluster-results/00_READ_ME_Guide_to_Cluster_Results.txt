GUIDE TO CLUSTER RESULTS:


NOTE: results cannot be uploaded to GitHub due to size. All results can be reproduced
     using the files provided.


Running the simulation code in the R scripts creates folders containing results for the global
     treatment effect estimates, region-specific treatment effects, rejection rates, bias
     and MSE (region-specific intercepts and region-specific treatment effects), PMPs,
     local constancy probabilities (Japanese MHLW approach and proposed leave-one-out absolute
     difference approach), and epsilon-level pairwise consistency probabilities



alt-half
   - alternative sample sizes half of null
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



diff-effects-1
   - equal regional sample sizes
   - N = 1508
   - underlying region-specific treatment effects: .017, .026, .034, .043, .051
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



diff-effects-2
   - equal regional sample sizes
   - N = 1508
   - underlying region-specific treatment effects: .017, .017, .017, .034, .034
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



diff-effects-3
   - equal regional sample sizes
   - N = 1508
   - underlying region-specific treatment effects: .017, .034, .034, .051, .051
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



n1508
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



n1508 - original results
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5
   - contains results for epsilon-level global and pairwise consistency using many epsilon values



n1508-alpha2
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 2
   - beta.star: .5



n1508-alpha4
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 4
   - beta.star: .5



n1508-alpha10
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 10
   - beta.star: .5



n1508-alphaNeg2
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: -2
   - beta.star: .5



n1508-alphaNeg4
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: -4
   - beta.star: .5



n1508-alphaNeg10
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: -10
   - beta.star: .5



n1508-beta20
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .2



n754
   - equal regional sample sizes
   - N = 754 (half of N = 1508)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



n754-beta20
   - equal regional sample sizes
   - N = 754 (half of N = 1508)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .2



n7650
   - equal regional sample sizes
   - N = 7650 (region sample sizes calculated to detect 90% power)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



run-simulation-n7650-beta20
   - equal regional sample sizes
   - N = 7650 (region sample sizes calculated to detect 90% power)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .2



SA1
   - equal regional sample sizes
   - N = 1508
   - m0: .082 for regional intercepts and .034 for regional treatment effects
   - Sigma0: diagonals of (10 x .082)^2 for regional intercepts and
     (10 x .034)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



SA2
   - equal regional sample sizes
   - N = 1508
   - m0: .082 for regional intercepts and .017 for regional treatment effects
   - Sigma0: diagonals of (10 x .082)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



SA3
   - equal regional sample sizes
   - N = 1508
   - m0: .082 for regional intercepts and .068 for regional treatment effects
   - Sigma0: diagonals of (10 x .082)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



SA4
   - equal regional sample sizes
   - N = 1508
   - m0: .041 for regional intercepts and .017 for regional treatment effects
   - Sigma0: diagonals of (10 x .041)^2 for regional intercepts and
     (10 x .017)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



SA5
   - equal regional sample sizes
   - N = 1508
   - m0: .041 for regional intercepts and .034 for regional treatment effects
   - Sigma0: diagonals of (10 x .041)^2 for regional intercepts and
     (10 x .034)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



SA6
   - equal regional sample sizes
   - N = 1508
   - m0: .041 for regional intercepts and .068 for regional treatment effects
   - Sigma0: diagonals of (10 x .041)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



SA7
   - equal regional sample sizes
   - N = 1508
   - m0: .164 for regional intercepts and .017 for regional treatment effects
   - Sigma0: diagonals of (10 x .164)^2 for regional intercepts and
     (10 x .017)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



SA8
   - equal regional sample sizes
   - N = 1508
   - m0: .164 for regional intercepts and .034 for regional treatment effects
   - Sigma0: diagonals of (10 x .164)^2 for regional intercepts and
     (10 x .034)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



SA9
   - equal regional sample sizes
   - N = 1508
   - m0: .164 for regional intercepts and .068 for regional treatment effects
   - Sigma0: diagonals of (10 x .164)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5

