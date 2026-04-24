* mus229bayes.do  for Stata 17

capture log close

********** OVERVIEW OF mus229bayes.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 29  BAYESIAN METHODS: BASICS
* 29.2: INTRODUCTION
* 29.4: AN IID EXAMPLE
* 29.6: A LINEAR REGRESSION EXAMPLE
* 29.7: MODIFYING THE MH ALGORITHM
* 29.8: RE MODEL
* 29.9: BAYESIAN MODEL SELECTION
* 29.10: BAYESIAN PREDICTION
* 29.11: PROBIT EXAMPLE

* To run you need files
*   mus229acs.dta
* in your directory

* No Stata community-contributed commands are used

********** SETUP **********

clear all
set linesize 80
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* File mus229acs.dta is authors' extract from 
* U.S. American Community Survey 2010
* 25-65 yearolds working >= 35 hours per week

*********** 29.2: INTRODUCTORY EXAMPLE

* Read in earnings - schooling data
qui use mus229acs
describe earnings lnearnings age education
keep if _n <= 100
summarize earnings lnearnings age education

* Bayesian linear regression with uninformative prior
bayes, rseed(10101): regress lnearnings education age

* ML linear regression (same as OLS with i.i.d. errors)
regress lnearnings education age, noheader

* Plot of the density of the 10,000 draws for each parameter
quietly bayes, rseed(10101): regress lnearnings education age
bayesgraph kdensity _all, combine 

sum

* Not included - how to interpret sig2 ~ inv gamma(.01,.01)
* Interpret as invsig2 ~ gamma(.01,.01)
preserve
clear
set obs 10000
gen invsig2 = rgamma(.01,100)
gen sig2 = 1/invsig2
sum, d
restore

*********** 29.4: AN IID EXAMPLE

* Generate a sample of 50 observations on y
clear
qui set obs 50 
set seed 10101
gen y = rnormal(10,10)
summarize

* The MLE is the sample mean
mean y

* Bayesian posterior for mu with normal y and N(5,4) prior for mu
bayesmh y, likelihood(normal(100)) prior({y:_cons}, normal(5,4)}) ///
   rseed(10101) saving(mcmcdraws_iid, replace)

* Bayesian hypothesis test: Pr[mu > 10]
bayestest interval {y:_cons}, lower(10)

* Bayesian statistics for transformation of parameter mu
bayesstats summary ({y:_cons}^2)

* Diagnostic plots for MH posterior draws
bayesgraph diagnostics {y:_cons}, scale(1.1) 

capture matrix drop pstart
capture matrix drop pmeans
capture matrix drop psds
capture matrix drop p_all

* Manually obtain posterior means and st. devs for four different chains
set seed 10101
forvalues i=1/4 {
    local start = rnormal(10,6^2)
    quietly bayesmh y, likelihood(normal(100)) prior({y:_cons}, normal(5,4)) ///
    initial({y:_cons} `start') 
    matrix pstart =  (nullmat(pstart) \ `start')
    matrix pmeans = (nullmat(pmeans) \ e(mean))
    matrix psds = (nullmat(psds) \ e(sd))
}
matrix p_all = pstart,pmeans,psds
matrix list p_all, title("Start value, post. means and st. devs. for 4 chains")

* Compute within and between variation across the 4 MCMC runs
mata
    pmeans = st_matrix("pmeans")
    psd = st_matrix("psds")
    W = (psd'psd)/4
    meanpmeans = mean(pmeans)
    B = (pmeans-meanpmeans*J(4,1,1))'(pmeans-meanpmeans*J(4,1,1))/(4-1)
    // (1) Between, (2) within, (3) Between/total, (4) Total / Within 
    (B, W, B/(B+W), (B+W)/W)
end

bayesmh y, nchains(4) likelihood(normal(100)) rseed(10101)         ///
    prior({y:_cons}, normal(5,4)) initsummary
	
* Four different chains using nchains() option
bayesmh y, nchains(4) likelihood(normal(100)) rseed(10101)      ///
    prior({y:_cons}, normal(5,4)) init1({y:_cons} 24)           ///
	init2({y:_cons} 7) init3({y:_cons} 11) init4({y:_cons} -23) ///
	initsummary nomodelsummary  
	
* Gelman-Rubin convergence diagnostic
bayesstats grubin

capture matrix drop pmeans
capture matrix drop psds
capture matrix drop p_all

* Compare posterior for two different priors
quietly bayesmh y, likelihood(normal(100)) prior({y:_cons}, normal(4,8)) ///
      rseed(10101)
matrix pmeans = e(mean)
matrix psds =  e(sd)
qui bayesmh y, likelihood(normal(100)) prior({y:_cons}, chi2(4)) ///
    rseed(10101)  
matrix pmeans = pmeans \ e(mean)
matrix psds = psds \ e(sd)
matrix p_all = pmeans,psds
matrix list p_all, title("Post. means and st. devs. for 2 different priors")

* Summarize the unique retained draws 
use mcmcdraws_iid, clear
summarize

* Expand to get the 10,000 MH draws, including repeated draws
expand _frequency
sort _index
gen s = _n
summarize eq1_p1

* Not included - get 95% credible interval - for some reason lower bound differs
centile eq1_p1, centile(2.5, 97.5)

* Graph the first 50 draws of mu
quietly tsset s
tsline eq1_p1 if s < 50, scale(1.5) ytitle("Parameter mu")      ///
   xtitle("MCMC posterior draw number") saving(graph1.gph, replace) 

* Plot prior, likelihood and posterior densities
graph twoway (function prior=(1/(sqrt(8*_pi))*exp(-((x-5)^2)/8)), lstyle(p3) range(0 16)) ///
    (kdensity eq1_p1, ytitle("Densities") xtitle("mu") clstyle(p1))                    ///
    (function likelihood=(1/(sqrt(4.752*_pi))*exp(-((x-10.7)^2)/4.752)), lstyle(p2) range(0 16)), ///
    scale(1.5) legend(pos(11) col(1) ring(0) lab(1 "prior")         ///
    lab(2 "posterior") lab(3 "likelihood")) saving(graph2.gph, replace)
graph combine graph1.gph graph2.gph, iscale(1.0) ysize(2.5) xsize(6) rows(1)

* Not used - plot the acf  
* ac eq1_p1, scale(1.5) ytitle("Autocorrelations of mu") ///
*   xtitle("Lag of draw") saving(graph2.gph, replace)

*********** 29.6: A LINEAR REGRESSION EXAMPLE

* MLE for the regression (same as OLS with i.i.d. errors)
qui use mus229acs, clear 
quietly keep if _n <= 100
regress lnearnings education age

* Bayesian posterior with informative priors: Normal for b, inv gamma for s2 
bayesmh lnearnings education age, likelihood(normal({var}))  ///
   prior({lnearnings:education}, normal(0.06,0.0001))        ///
   prior({lnearnings:age}, normal(0.02,0.0001))              ///
   prior({lnearnings:_cons}, normal(10,100)) prior({var}, igamma(1,0.5)) ///
   rseed(10101) saving(mcmcdraws_fullregress, replace) 
estimates store fullregress

* MCMC statistics for all parameters 
bayesstats ess

* Diagnostic plots for MH posterior draws of beta_education
bayesgraph diagnostics {lnearnings:education}

* Trace plot for all four parameters
bayesgraph trace _all, combine

* Not included: matrix plot of parameter draws against each other 
bayesgraph matrix _all
* Bayesian posterior with informative priors: normal for b, inv gamma for s2 
* Change the prior on the intercept

* Bayesian posterior with prior for intercept tighter and centered on zero
qui bayesmh lnearnings education age, likelihood(normal({var}))    ///
    rseed(10101) prior({lnearnings:education}, normal(0.06,0.01))  ///
    prior({lnearnings:age}, normal(0.02,0.01))                     ///
    prior({lnearnings:_cons}, normal(0,1)) prior({var}, igamma(1,0.5))  
bayesstats summary

/*
qui bayesmh lnearnings education educcopy age,           /// 
    rseed(10101) likelihood(normal({var})) burnin(20000) ///
    prior({lnearnings:education}, normal(0.06,0.0001))   ///
    prior({lnearnings:educcopy}, normal(0.010,0.0001))   ///
    prior({lnearnings:age}, normal(0.02,0.01))           ///
    prior({lnearnings:_cons}, normal(10,100)) prior({var}, igamma(1,0.5))  
bayesstats summary
*/

* Bayesian posterior with same regressor appearing twice and informative prior
generate educcopy = education   
qui bayesmh lnearnings education educcopy,             /// 
    rseed(10101) likelihood(normal({var}))             ///
    prior({lnearnings:education}, normal(0.06,0.0001)) ///
    prior({lnearnings:educcopy}, normal(0.10,0.0001))  ///
    prior({lnearnings:_cons}, normal(10,100)) prior({var}, igamma(1,0.5))  
bayesstats summary

// Not included
bayesstats ess

* Posterior density for first and second half of draws for all parameters
bayesgraph kdensity _all, show(both) combine

// Command and output not included - multiple chains 
qui bayesmh lnearnings education educcopy, nchains(5)    /// 
    rseed(10101) likelihood(normal({var}))               ///
    prior({lnearnings:education}, normal(0.06,0.0001))   ///
    prior({lnearnings:educcopy}, normal(0.10,0.0001))    ///
    prior({lnearnings:_cons}, normal(10,100)) prior({var}, igamma(1,0.5))  

* Rc statistic from running five chains on preceding bayesmh command
bayesstats grubin

// Command and output not included - multiple chains and 20000 burnin
qui bayesmh lnearnings education educcopy, nchains(5) burnin(20000) /// 
    rseed(10101) likelihood(normal({var}))                  ///
    prior({lnearnings:education}, normal(0.06,0.0001))      ///
    prior({lnearnings:educcopy}, normal(0.10,0.0001))       ///
    prior({lnearnings:_cons}, normal(10,100)) prior({var}, igamma(1,0.5))  
bayesstats grubin
  
// Not included - increase prior variances to one
bayesmh lnearnings education educcopy,  burnin(20000) /// 
    rseed(10101) likelihood(normal({var}))            ///
    prior({lnearnings:education}, normal(0.06,1))     ///
    prior({lnearnings:educcopy}, normal(0.10,1))      ///
    prior({lnearnings:_cons}, normal(10,100)) prior({var}, igamma(1,0.5))  
bayesgraph kdensity _all, show(both) combine
  
*********** 29.7: MODIFYING THE MH ALGORITHM

* MH with blocking: var in separate block 
qui bayesmh lnearnings education age, likelihood(normal({var})) ///
    prior({lnearnings:education}, normal(0.06,0.0001)) ///
    prior({lnearnings:age}, normal(0.02,0.0001))       ///
    prior({lnearnings:_cons}, normal(10,100))           ///
    prior({var}, igamma(1,0.5)) block({var}) rseed(10101)
bayesstats summary
bayesstats ess
di "Overall acceptance rate = " e(arate)  

* Hybrid MH with Gibbs sampling subcomponent 
qui bayesmh lnearnings education age, likelihood(normal({var})) ///
    prior({lnearnings:education}, normal(0.06,0.0001)) ///
    prior({lnearnings:age}, normal(0.02,0.0001))       ///
    prior({lnearnings:_cons}, normal(10,100))           ///
    prior({var}, igamma(1,0.5)) block({var}, gibbs) rseed(10101)
bayesstats summary
bayesstats ess
di "Overall acceptance rate = " e(arate) 

*********** 29.8: HIERARCHICAL PRIOR FOR RANDOM EFFECTS

* Here regress earnings on inercept and age
* The random effect varies with the level of education

* Mixed-effects estimation of random-effects model
qui use mus229acs, clear
quietly keep if _n <= 100
mixed lnearnings age || education: , nolog

* Predict the random intercepts = intercept + RE 
predict u0, reffects
by education, sort: generate tolist = (_n==1)
generate randomint = _b[_cons] + u0
list education u0 randomint if tolist, clean

// Not included 
bayes: mixed lnearnings age || education:
  
. * Bayesian hierarchical model using bayesmh with multilevel specification
bayesmh lnearnings age V[education], rseed(10101)     ///
     likelihood(normal({var_0})) showreffects          ///
     prior({lnearnings:_cons}, normal(10,100))         ///
     prior({lnearnings:age}, normal(0,0.01))           ///
     prior({var_0}, igamma(0.0001, 0.0001))            ///
     prior({var_V}, igamma(0.0001, 0.0001))            ///
     block({lnearnings:} {var_0} {var_V} {V}, gibbs split)

// Not included
bayesgraph kdensity _all, show(both) combine

// Output not included but text included
* Bayesian hierarchical model random effects using split in block
fvset base none education     // Needed as constant omitted from regression
bayesmh lnearnings age i.education, noconstant rseed(10101)                 ///
    likelihood(normal({var_0})) burnin(20000)                               ///
    prior({lnearnings:i.education}, normal({lnearnings:_cons},{var_ed}))    ///
    block({lnearnings:i.education}, split)                                  ///
    prior({lnearnings:_cons}, normal(10,100)) block({lnearnings:_cons}, gibbs) ///
    prior({var_0}, igamma(0.0001, 0.0001)) block({var_0}, gibbs)             ///
    prior({lnearnings:age}, normal(0,0.01)) block({lnearnings:age}, gibbs)   ///
    prior({var_ed}, igamma(0.0001, 0.0001)) block({var_ed}, gibbs)           

// Not included
bayesgraph kdensity _all, show(both) combine

// Not included 
* Bayesian hierarchical model random effects without gibbs and split
fvset base none education // Need this as constant omitted from regression
bayesmh lnearnings age i.education,          ///
    likelihood(normal({var_0})) noconstant   /// 
    prior({lnearnings:i.education}, normal({lnearnings:_cons},{var_ed})) ///
    prior({lnearnings:_cons}, normal(10,100)) ///
    prior({var_0}, igamma(0.0001, 0.0001))   ///
    prior({lnearnings:age}, normal(0,0.01))  ///
    prior({var_ed}, igamma(0.0001, 0.0001))  ///
    mcmcsize(10000) burnin(10000) nomodelsummary rseed(10101) 
bayesgraph kdensity _all, show(both) combine
bayesgraph trace _all, combine

*********** 29.9: BAYESIAN MODEL SELECTION
 
* Bayesian posterior for small model without regressor education
qui bayesmh lnearnings age, likelihood(normal({var}))                   ///
    prior({lnearnings:age}, normal(0.02,0.0001))                        ///
    prior({lnearnings:_cons}, normal(10,100)) prior({var}, igamma(1,0.5)) ///
    saving(mcmcdraws_smallregress, replace) rseed(10101)
estimates store smallregress
di "Marginal likelihood for smaller model "  e(lml_lm) 

* Bayes factor with first-listed model the base model
bayesstats ic fullregress smallregress, bayesfactor

* Posterior odds Bayes factor with first-listed model the base model
bayestest model fullregress smallregress, prior(0.8 0.2)

* ML analysis
quietly regress lnearnings education age
scalar llfull = e(ll)
quietly regress lnearnings age
di "LL full = " llfull " LL small = " e(ll) "  Diff = " llfull-e(ll)

*********** 29.10:  BAYESIAN PREDICTION

* Bayesian posterior with informative priors: Normal for b, inv gamma for s2 
qui bayesmh lnearnings education age, likelihood(normal({var}))     ///
    rseed(10101) prior({lnearnings:education}, normal(0.06,0.0001)) ///
    prior({lnearnings:age}, normal(0.02,0.0001))                    ///
    prior({lnearnings:_cons}, normal(10,100))                       ///
    prior({var}, igamma(1,0.5)) saving(mcmcdraws_fullregress, replace) 

* Create 2 new observations (numbers 101 and 102) for prediction 
set obs 102
quietly replace education = 12 if _n == 101
quietly replace age = 40 if _n == 101
quietly replace education = 16 if _n == 102
quietly replace age = 40 if _n == 102

* Obtain 10,000 MCMC predictions for each of the 2 individuals
bayespredict {_ysim1} if _n > 100, saving(mcmcpredict, replace) rseed(10101)

* Summarize the MCMC predictions dataset
preserve
use mcmcpredict, clear
summarize
restore

* Summarize the MCMC predictions
bayesstats summary {_ysim} using mcmcpredict

* Directly obtain the mean of the posterior predictive distribution.
bayespredict plnearnings if _n > 100, mean rseed(10101)
list plnearnings if _n > 100, clean

*********** 29.10:  PROBIT EXAMPLE

* Create dependent variable for probit example
qui use mus229acs, clear 
quietly keep if _n <= 100
generate dhighearns = earnings > 75000
quietly save musrearnings, replace
summarize dearnings age education

* MLE for the probit regression
probit dhighearns education age, nolog

// Not included - use simpler bayes: command
bayes, rseed(10101): probit dhighearns education age

* Bayesian posterior for probit regression with flat priors for beta 
bayesmh dhighearns education age,  rseed(10101) likelihood(probit) ///
    prior({dhighearns:}, flat) saving(mcmcdraws_probit, replace)           

/* Not included - probit with tighter normal priors
bayesmh dhighearns education age, likelihood(probit)   ///
    prior({dhighearns: education age}, normal(0,0.1))  ///
    prior({dhighearns: _cons}, normal(0,100)) rseed(10101) 
*/	

* Program to compute AME = phi(x'b)*b_j to be used by bayespredict
program ameprog
    version 17.0
    args sum mu
    // Define locals with shorter names for convenience 
    local touse $BAYESPR_touse
    local theta $BAYESPR_theta
    local j $BAYESPR_passthruopts
    // Obtain the current MCMC value of the jth coefficient
    tempname betaj
    scalar `betaj' = `theta'[1, `j']
    // Compute phi(xb) and store it in temporary variable tmpv 
    tempvar tmpv
    generate double `tmpv' = normalden(invnormal(`mu')) if `touse'
    summarize `tmpv' if `touse', meanonly
    // Store the final result, AME, in temporary scalar sum 
    scalar `sum' = r(mean)*`betaj'
end

* Program that computes AMEs in each of 10,000 (the default) MCMC draws
bayespredict (ame1:@ameprog {_mu1}, passthruopts(1)) ///
    (ame2:@ameprog {_mu1}, passthruopts(2)) ///
    (ame3:@ameprog {_mu1}, passthruopts(3)), saving(amepred, replace)
	
* Summary statistics for the AMEs from 10,000 MCMC draws
bayesstats summary {ame1} {ame2} {ame3} using amepred

* Compare Bayesian AME to AME from ML estimation
qui probit dhighearns education age
margins, dydx(*)

*** Erase files created by this program and not used elsewhere in the book

erase mcmcdraws_iid.dta
erase musrearnings.dta
erase mcmcdraws_probit.dta
erase mcmcdraws_fullregress.dta
erase mcmcdraws_smallregress.dta
erase graph1.gph
erase graph2.gph
erase mcmcpredict.dta
erase mcmcpredict.ster
erase amepred.dta
erase amepred.ster

********** END
