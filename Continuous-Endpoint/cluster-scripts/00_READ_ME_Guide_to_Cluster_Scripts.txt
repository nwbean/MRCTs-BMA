GUIDE TO CLUSTER SCRIPTS:


Batch files used to begin simulations on cluster



batch-simulation-alt-half.sh - runs simulation with the following details:
   - alternative sample sizes half of null
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-diff-effects-1.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - underlying region-specific treatment effects: .017, .026, .034, .043, .051
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-diff-effects-2.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - underlying region-specific treatment effects: .017, .017, .017, .034, .034
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-diff-effects-3.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - underlying region-specific treatment effects: .017, .034, .034, .051, .051
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-n1508.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-n1508-alpha2.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 2
   - beta.star: .5



batch-simulation-n1508-alpha4.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 4
   - beta.star: .5



batch-simulation-n1508-alpha10.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 10
   - beta.star: .5



batch-simulation-n1508-alphaNeg2.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: -2
   - beta.star: .5



batch-simulation-n1508-alphaNeg4.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: -4
   - beta.star: .5



batch-simulation-n1508-alphaNeg10.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: -10
   - beta.star: .5



batch-simulation-n1508-beta20.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .2



batch-simulation-n754.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 754 (half of N = 1508)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-n754-beta20.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 754 (half of N = 1508)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .2



batch-simulation-n7650.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 7650 (region sample sizes calculated to detect 90% power)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-n7650-beta20.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 7650 (region sample sizes calculated to detect 90% power)
   - m0: .10 for regional intercepts and .04 for regional treatment effects
   - Sigma0: diagonals of (10 x .10)^2 for regional intercepts and
     (10 x .04)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .2



batch-simulation-SA1.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .082 for regional intercepts and .034 for regional treatment effects
   - Sigma0: diagonals of (10 x .082)^2 for regional intercepts and
     (10 x .034)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-SA2.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .082 for regional intercepts and .017 for regional treatment effects
   - Sigma0: diagonals of (10 x .082)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-SA3.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .082 for regional intercepts and .068 for regional treatment effects
   - Sigma0: diagonals of (10 x .082)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-SA4.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .041 for regional intercepts and .017 for regional treatment effects
   - Sigma0: diagonals of (10 x .041)^2 for regional intercepts and
     (10 x .017)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-SA5.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .041 for regional intercepts and .034 for regional treatment effects
   - Sigma0: diagonals of (10 x .041)^2 for regional intercepts and
     (10 x .034)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-SA6.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .041 for regional intercepts and .068 for regional treatment effects
   - Sigma0: diagonals of (10 x .041)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-SA7.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .164 for regional intercepts and .017 for regional treatment effects
   - Sigma0: diagonals of (10 x .164)^2 for regional intercepts and
     (10 x .017)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-SA8.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .164 for regional intercepts and .034 for regional treatment effects
   - Sigma0: diagonals of (10 x .164)^2 for regional intercepts and
     (10 x .034)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5



batch-simulation-SA9.sh - runs simulation with the following details:
   - equal regional sample sizes
   - N = 1508
   - m0: .164 for regional intercepts and .068 for regional treatment effects
   - Sigma0: diagonals of (10 x .164)^2 for regional intercepts and
     (10 x .068)^2 for regional treatment effects
   - alpha0: 0
   - beta.star: .5

