* ct_ivgmm.do  based on mus06p1iv.do  May 2013 for Stata version 12.0

log using ct_ivgmm.txt, text replace 

********** OVERVIEW OF ct_ivgmm.do **********

* STATA Program 
* For A. Colin Cameron and Pravin K. Trivedi "Lectures in Microeconometrics" 
* IV, 2SLS and linear GMM

* To run you need file
*   mus06data.dta    
* in your directory

* Stata add-on condivreg is used

********** SETUP **********

set more off
version 12.0
set mem 10m
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* The original data is from the Medical Expenditure Panel Survey (MEPS)
* for individuals over the age of 65 years

********** INSTRUMENTAL VARIABLES EXAMPLE

* Read data, define global x2list (exogenous regressors), and summarize data
use mus06data.dta
global x2list totchr age female blhisp linc
keep if linc != .
describe ldrugexp hi_empunion $x2list
summarize ldrugexp hi_empunion $x2list 

* Define the exogenous regressors using global x2list
global x2list totchr age female blhisp linc

* OLS 
regress ldrugexp hi_empunion $x2list, vce(robust)

* Two available instruments for hi_empunion
describe ssiratio multlc
summarize ssiratio multlc 
correlate hi_empunion ssiratio multlc 

* IV estimator with ssiratio as single instrument for hi_empunion
ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust)

* 2SLS estimator with ssiratio and multlc as instruments for hi_empunion
ivregress 2sls ldrugexp (hi_empunion = ssiratio multlc) $x2list, vce(robust)

* Robust version of Hausman test using augmented regression
quietly ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust)
estat endogenous

* Robust Hausman-Wu test of endogeneity implemented manually
quietly regress hi_empunion ssiratio $x2list
quietly predict v1hat, resid
quietly regress ldrugexp hi_empunion v1hat $x2list, vce(robust)
test v1hat 

* Weak instrument tests - just-identified model
quietly ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust)
estat firststage, forcenonrobust all  

* Conditional test and confidence intervals when weak instruments 
condivreg ldrugexp (hi_empunion = ssiratio) $x2list, lm ar 2sls test(0)

* GMM estimator with ssiratio and multlc as instruments for hi_empunion
ivregress gmm ldrugexp (hi_empunion = ssiratio multlc) $x2list, vce(robust)

* Test of overidentifying restrictions following ivregress gmm
quietly ivregress gmm ldrugexp (hi_empunion = ssiratio multlc) $x2list, wmatrix(robust) 
estat overid

* Compare estimators 
quietly regress ldrugexp hi_empunion $x2list, vce(robust)
estimates store OLS
quietly ivregress 2sls ldrugexp (hi_empunion = ssiratio multlc) $x2list, vce(robust)
estimates store IV
quietly ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust)
estimates store TWOSLS
quietly ivregress gmm ldrugexp (hi_empunion = ssiratio multlc) $x2list, vce(robust)
estimates store GMM
estimates table OLS IV TWOSLS GMM, b(%9.4f) se(%9.3f) stats(N r2 F)

* CLOSE OUTPUT           
* log close
* clear
* exit

