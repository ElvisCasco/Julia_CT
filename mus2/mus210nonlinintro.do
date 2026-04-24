* mus210nonlinintro.do  for Stata 17

capture log close

********** OVERVIEW OF mus210nonlinintro.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 10: INTRODUCTION TO NONLINEAR REGRESSION
* 10.2: BINARY OUTCOME MODELS
* 10.3: PROBIT MODEL
* 10.4: MES AND COEFFICIENT INTERPRETATION
* 10.5: LOGIT MODEL
* 10.6: NONLINEAR LEAST SQUARES
* 10.7: OTHER NONLINEAR ESTIMATORS

* To run you need files
*   mus210mepsdocvisyoung.dta   
* in your directory

* No community-contributed commands are used

* To speed up program reduce reps in vce(bootstrap, reps(400))

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus210mepsdocvisyoung.dta from 2002 Medical Expenditure Panel Survey (MEPS)
* U.S. data on office-based physician visits by persons aged 25-64 years
* Same as used by Deb, Munkin and Trivedi (2006) 
* Excludes those receiving public insurance (Medicare and Medicaid)
* Restricted to those working in the private sector but not self-employed.

********** 10.2: BINARY OUTCOME MODELS

* Read in dataset, select one year of data, and describe key variables
qui use mus210mepsdocvisyoung
qui keep if year02 == 1
generate visit = docvis > 0 
label variable visit "= 1 if doctor visit"
describe visit private chronic female income

* Summmary of key variables
summarize visit private chronic female income

* The following constructs Figure 1
graph twoway function y=normal(0.5+x), range(-2.5 2.5)            ///
    title("Probit model with Pr(y=1) = {&Phi}(0.5+x)")            ///
    ytitle("Probability that y = 1") xtitle("Scalar regressor x") ///
    saving(graph1.gph, replace)
graph twoway function y=normalden(0.5+x), range(-2.5 2.5)         ///
    title("ME in probit model") xtitle("Scalar regressor x") ///
    ytitle("ME dPr(y=1)/dx") saving(graph2.gph, replace)
graph combine graph1.gph graph2.gph, ysize(2.5) xsize(6)          ///
    xcommon iscale(1.2) rows(1) 

********** 10.3: PROBIT MODEL

* Probit regression (command probit) with robust standard errors
probit visit private chronic female income, vce(robust)

* Probit regression (command probit) with default standard errors
probit visit private chronic female income, nolog

/*
* Not included - clustered standard errors example: cluster on age
probit visit private chronic female income, vce(cluster age) nolog
* Not included - bootstrap standard errors
probit visit private chronic female income, vce(bootstrap, reps(400)) nolog
probit visit private chronic female income, vce(bootstrap, cluster(age) reps(400)) nolog
*/

* Predicted probabilities from probit
qui probit visit private chronic female income, vce(robust)
predict phatprobit, pr
summarize visit phatprobit

********** 10.4 MARGINAL EFFECTS

* AMEs using calculus method
qui probit visit private chronic female income, vce(robust)
margins, dydx(*) 

* AMEs using finite-difference method
qui probit visit i.private i.chronic i.female income, vce(robust)
margins, dydx(*) 

* AMEs with interacted regressors
probit visit i.private##c.income i.chronic i.female, vce(robust) nolog noheader
margins, dydx(*) 

* OLS regression (command regress) for linear probability model
regress visit private chronic female income, vce(robust) noheader
predict phatols, xb
summarize visit phatols 

********** 10.5 LOGIT MODEL

* Logit regression (command logit)
logit visit private chronic female income, vce(robust)

* AMEs from logit regression
margins, dydx(*)

* Not included: Predicted probabilities from logit
predict phatlogit, pr
sum visit phatlogit
correlate phatprobit phatlogit

********** 10.6 NONLINEAR LEAST SQUARES

* Nonlinear least-squares regression (command nl) for probit model
nl (visit = normal({xb: private chronic female income}+{b0})), vce(robust)

********** 12.7 OTHER NONLINEAR ESTIMATORS 

* THIS DOES NOT SEEM TO HAVE BEEN INCLUDED

* Maximum likeihood estimation (command mlexp) for probit model
gen int q = 2*visit - 1
mlexp (ln(normal(q*{xb: private chronic female income _cons}))), vce(robust)

* Generalized linear model estimation (command glm) for probit model
glm visit private chronic female income, family(binomial) link(probit) vce(robust) 

// Not included - plot logit versus probit functions
preserve
clear
set obs 100
gen z = (_n - 50)/20
gen probit = normal(z)
gen logitrescaled = exp(1.7*z) / (1+exp(1.7*z))
twoway (line logitrescaled z) (line probit z)
restore

*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph

********** END



