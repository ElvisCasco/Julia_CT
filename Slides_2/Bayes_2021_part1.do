* Bayes_2021_part1.do May 2021 for Stata version 17 9also works in 16)
* based on mus229bayes.do  May 2021 for Stata version 17

capture log close

log using Bayes_2021_part1.txt, text replace

********** OVERVIEW OF Bayes_2021_part1.do **********

* Stata program based on Chapter 29 of
* A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2021)
* Stata Press

* To run you need files
*   mus229acs.dta
* in your directory

* No Stata user-written commands are used

********** SETUP **********

set more off
* version 17
clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* File mus229acs.dta is authors' extract from 
* U.S. American Community Survey 2010
* 25-65 yearolds working >= 35 hours per week

*********** 29.2: INTRODUCTORY EXAMPLE

* Read in and summarize earnings - schooling data
qui use mus229acs.dta, clear 
describe earnings lnearnings age education
qui keep if _n <= 100
summarize earnings lnearnings age education

* ML linear regression (same as OLS with iid errors)
regress lnearnings education age, noheader

* Bayesian linear regression with uninformative prior and Stata defaults
bayes, rseed(10101): regress lnearnings education age

* Diagnostic plots for MH posterior draws of beta_education
bayesgraph diagnostics {lnearnings:education}

* Check convergence using multiple chains
bayes, rseed(10101) nchains(5): regress lnearnings education age

* Give Gelman-Rubin Rc statistic for each parameter
bayesstats grubin 

* bayemh comamnd that gives same results as earlier (needs block to do this)
bayesmh lnearnings education age, likelihood(normal({sigma2}))  ///
   prior({lnearnings:education}, normal(0,10000))        ///
   prior({lnearnings:age}, normal(0,10000))              ///
   prior({lnearnings:_cons},normal(0,10000))             ///
   prior({sigma2},igamma(0.01,0.01)) rseed(10101)        ///
   block({lnearnings: education age _cons}) block({sigma2}) 

* bayesmh example with informative priors	   
bayesmh lnearnings education age, likelihood(normal({var}))  ///
   prior({lnearnings:education}, normal(0.06,0.0001))        ///
   prior({lnearnings:age}, normal(0.02,0.0001))              ///
   prior({lnearnings:_cons},normal(10,100))                  ///
   prior({var},igamma(1,0.5)) rseed(10101)

********** END
