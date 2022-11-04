******************************************************************************************************;
* Reproduce Primary Analysis and EDA Post Hoc Analysis
******************************************************************************************************;


/* Read in datasets */
*****CODE FOR DATA PREPARATION REMOVED*****;


/* Planned analysis using Cox's proportional hazards model – treatment only */
proc phreg data=ttedata;
    class TRT(ref="0");
    model EVENTTIME*STATUS(0) = TRT / alpha=.05 risklimits=both;
run;


/* Planned analysis using Cox's proportional hazards model - treatment and region with interaction */
proc phreg data=ttedata;
    class TRT(ref="0") REGION (ref="Rest of the World");
    model EVENTTIME*STATUS(0) = TRT REGION TRT*REGION / alpha=.05 risklimits=both;
    hazardratio TRT / diff=ref;
run;


/* FDA post hoc Bayesian hierarchical model on the observed region-specific log hazard ratios */
data bhmdat;
    input region $ yi sig2i;
    datalines;
    Asia -0.47446 0.06869369
    Europe -0.20469 0.008799172
    NoAmer 0.00965 0.009346611
    RoW -0.18237 0.01135161
;
run;


proc mcmc data=bhmdat seed=801 nbi=50000 nmc=10000 thin=10;
    parm mu tau;
    prior mu ~ normal(0, var=16);
    prior tau ~ gamma(shape=.001, iscale=.001);
    random mui ~ normal(mu, prec=tau) subject=region monitor=(mui);
    model yi ~ n(mui, var=sig2i);
run;
