# semicmprsk
Separable pathway effects of semi-competing risks via multi-state models

## semicompeting.R
Main functions
#### semicr.causal(A, T, dT, R, dR, X=NULL, weights=NULL, a1, a0, asm='markov', conf.int=NULL, nboot=0, ...)
#### Compare counterfactual cumulative incidences
    Inputs
    A: treatment
    T: time to terminal event
    dT: censoring indicator of terminal event
    R: time to intermediate event
    dR: censoring indicator of intermediate event
    X: covariates matrix
    weights: weights of each unit
    a1: the first hypothetical treatment vector under comparison (as treated group)
    a0: the second hypothetical treatment vector under comparison (as control group)
    asm: 'markov' for Markovness, 'semimarkov' for semi-Markovness
    conf.int: significance level (Null for no confidence interval)
    nboot: number of bootstrap resamplings to calculate confidence intervals (0 or 1 for analytic form)
#### semicr.sensitivity(A, T, dT, R, dR, X=NULL, weights=NULL, a1, a0, sens=seq(0,1,0.2), ...)
#### Senitivity analysis for pointwise treatment effect
    Inputs
    A: treatment
    T: time to terminal event
    dT: censoring indicator of terminal event
    R: time to intermediate event
    dR: censoring indicator of intermediate event
    X: covariates matrix
    weights: weights of each unit
    a1: the first hypothetical treatment vector under comparison (as treated group)
    a0: the second hypothetical treatment vector under comparison (as control group)
    sens: sequence of sensitivity parameters, between 0 and 1
#### testpath(A,T,dT,R,dR,X=NULL,weights=NULL)
#### Test separable pathway effects
    Outputs
    p1: p-value for the effect on path 0 -> 1
    p2: p-value for the effect on path 0 -> 2
    p3: p-value for the effect on path 2 -> 3 (the first for Markov, the second for semi-Markov)
    p0: p-value for overall effect (intention-to-treat)
    p23: p-value for the effect on path 0 -> 3, considering 1 and 3 as competing events
    p01: p-value for the effect on path 0 -> 1, considering 1 and 3 as competing events

## simulation1.R
Simulation (Part 1): generate data, draw estimated cumulative incidences, calculate bias

## simulation2.R
Simulation (Part 2): generate data, test pathway effects, evaluate confidence intervals

## sensitivity.R
Sensitivity analysis, where the dismissible treatment decomposition assumption is violated

## data.R
Real data application
