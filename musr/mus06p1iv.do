* mus06p1iv.do  Oct 2009 for Stata version 11

cap log close

********** OVERVIEW OF mus06p1iv.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Chapter 6
* 6.3: INSTRUMENTAL VARIABLES EXAMPLE
* 6.4: WEAK INSTRUMENTS
* 6.5: BETTER INFERENCE WITH WEAK INSTRUMENTS
* 6.6: 3SLS SYSTEMS ESTIMATION

* To run you need files
*   mus06data.dta    
* in your directory

* Stata user-written commands
*   condivreg
*   ivreg2
*   jive
* are used

********** SETUP **********

set more off
version 11
clear all
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* The original data is from MEPS over 65 similar to chapter 3

********** 6.3: INSTRUMENTAL VARIABLES EXAMPLE

* Read data, define global x2list, and summarize data
use mus06data.dta
global x2list totchr age female blhisp linc 
summarize ldrugexp hi_empunion $x2list

// Following not included in book
summarize

* Summarize available instruments 
summarize ssiratio lowincome multlc firmsz if linc!=.

* IV estimation of a just-identified model with single endog regressor
ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust) first

* Compare 5 estimators and variance estimates for overidentified models
global ivmodel "ldrugexp (hi_empunion = ssiratio multlc) $x2list"
quietly ivregress 2sls $ivmodel, vce(robust)
estimates store TwoSLS
quietly ivregress gmm  $ivmodel, wmatrix(robust) 
estimates store GMM_het
quietly ivregress gmm  $ivmodel, wmatrix(robust) igmm
estimates store GMM_igmm
quietly ivregress gmm  $ivmodel, wmatrix(cluster age) 
estimates store GMM_clu
quietly ivregress 2sls  $ivmodel
estimates store TwoSLS_def
estimates table TwoSLS GMM_het GMM_igmm GMM_clu TwoSLS_def, b(%9.5f) se  

* Obtain OLS estimates to compare with preceding IV estimates
regress ldrugexp hi_empunion $x2list, vce(robust) 

* Robust Durbin-Wu-Hausman test of endogeneity implemented by estat endogenous
ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust)
estat endogenous

* Robust Durbin-Wu-Hausman test of endogeneity implemented manually
quietly regress hi_empunion ssiratio $x2list
quietly predict v1hat, resid
quietly regress ldrugexp hi_empunion v1hat $x2list, vce(robust)
test v1hat 

* Test of overidentifying restrictions following ivregress gmm
quietly ivregress gmm ldrugexp (hi_empunion = ssiratio multlc) ///
  $x2list, wmatrix(robust) 
estat overid

* Test of overidentifying restrictions following ivregress gmm
ivregress gmm ldrugexp (hi_empunion = ssiratio lowincome multlc firmsz) ///
  $x2list, wmatrix(robust) 
estat overid

* Regression with a dummy variable regressor
treatreg ldrugexp $x2list, treat(hi_empunion = ssiratio $x2list)

********** 6.4: WEAK INSTRUMENTS

* Correlations of endogenous regressor with instruments
correlate hi_empunion ssiratio lowincome multlc firmsz if linc!=.

* Weak instrument tests - just-identified model
quietly ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust)
estat firststage, forcenonrobust all  

* Weak instrument tests - two or more overidentifying restrictions
quietly ivregress gmm ldrugexp (hi_empunion = ssiratio lowincome multlc firmsz) ///
   $x2list, vce(robust)
estat firststage, forcenonrobust

* Compare 4 just-identified model estimates with different instruments
quietly regress ldrugexp hi_empunion $x2list, vce(robust)
estimates store OLS0
quietly ivregress 2sls ldrugexp (hi_empunion=ssiratio) $x2list, vce(robust)
estimates store IV_INST1
quietly estat firststage, forcenonrobust
scalar me1 = r(mineig)
quietly ivregress 2sls ldrugexp (hi_empunion=lowincome) $x2list, vce(robust)
estimates store IV_INST2
quietly estat firststage, forcenonrobust
scalar me2 = r(mineig)
quietly ivregress 2sls ldrugexp (hi_empunion=multlc) $x2list, vce(robust) 
estimates store IV_INST3
quietly estat firststage, forcenonrobust
scalar me3 = r(mineig)
quietly ivregress 2sls ldrugexp (hi_empunion=firmsz) $x2list, vce(robust)
estimates store IV_INST4
quietly estat firststage, forcenonrobust
scalar me4 = r(mineig)
estimates table OLS0 IV_INST1 IV_INST2 IV_INST3 IV_INST4, b(%8.4f) se  
display "Minimum eigenvalues are:     " me1 _s(2) me2 _s(2) me3 _s(2) me4

********** 6.5: BETTER INFERENCE WITH WEAK INSTRUMENTS

* Conditional test and confidence intervals when weak instruments 
condivreg ldrugexp (hi_empunion = ssiratio) $x2list, lm ar 2sls test(0)

* Variants of IV Estimators: 2SLS, LIML, JIVE, GMM_het, GMM-het using IVREG2
global ivmodel "ldrugexp (hi_empunion = ssiratio lowincome multlc firmsz) $x2list"
quietly ivregress 2sls $ivmodel, vce(robust)
estimates store TWOSLS
quietly ivregress liml $ivmodel, vce(robust)
estimates store LIML
quietly jive $ivmodel, robust
estimates store JIVE
quietly ivregress gmm $ivmodel, wmatrix(robust) 
estimates store GMM_het
quietly ivreg2 $ivmodel, gmm robust
estimates store IVREG2
estimates table TWOSLS LIML JIVE GMM_het IVREG2, b(%7.4f) se 

********** 6.6: 3SLS SYSTEMS ESTIMATION

* 3SLS estimation requires errors to be homoskedastic
reg3 (ldrugexp hi_empunion totchr age female blhisp linc) ///
  (hi_empunion ldrugexp totchr female blhisp ssiratio)


********** CLOSE OUTPUT **************
