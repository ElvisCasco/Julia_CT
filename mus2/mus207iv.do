* mus207iv.do  for Stata 17
  
capture log close

********** OVERVIEW OF mus207iv.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press 

* Chapter 7
* 7.2: SIMULTANEOUS EQUATIONS MODEL
* 7.3: IV ESTIMATION
* 7.4: IV EXAMPLE
* 7.5: WEAK INSTRUMENTS
* 7.6: DIAGNOSTICS AND TESTS FOR WEAK INSTRUMENTS
* 7.7: INFERENCE WITH WEAK INSTRUMENTS 
* 7.8: FINITE SAMPLE INFERENCE WITH WEAK INSTRUMENTS
* 7.8: THREE-STAGE LEAST-SQUARES SYSTEMS ESTIMATION

* To run you need files
*   mus207mepspresdrugs.dta    
* in your directory

* community-contributed commands
*   condivreg
*   ivreg2
*   jive
*   rivtest
*   weakiv (needs avar)
*   avar   (not directly used but needed for weakiv)
*   weakivtest
* are used

********** SETUP **********

clear all
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* The data is an extract from the U.S. Medical Expenditure Panel Survey 
* of individuals 65 years and older
* This file is a slightly reduced version of mus06data.dta used in 
* the first and revised editions of this book. 
* Reduction is due to rejection of ssiratio observations outside [0,1] interval
* The empirical estimates change little.  

********** CHAPTER 7.2: SIMULTANEOUS EQUATIONS MODEL

* Generate data for simultaneous equations model with b12=1; b21=c12=c21=1
qui set obs 1000
set seed 10101
matrix C = (1, .7, 1)                            // Variances 1, covariance 0.7
drawnorm u1 u2, n(1000) corr(C) cstorage(lower)  // Bivariate normal (u1, u2)
matrix C = (1, .3, 1)
drawnorm x1 x2, n(1000) corr(C) cstorage(lower)  // Bivariate normal (x1, x2)
generate y1 =(1/(1-.25))*(x1 + 1*(x2+u2) + u1)   // Reduced form for y1
generate y2 = 0.25*y1 + 1*x2 + u2                // Generating y2 given y1
summarize y1 y2 x1 x2 u1 u2
correlate y1 y2 x1 x2 u1 u2

* OLS is inconsistent 
regress y1 y2 x1, vce(robust) noheader

* IV with valid instrument (here x2 for y2) are consistent 
ivregress 2sls y1 (y2=x2) x1, vce(robust) noheader

* Not included - systems estimation of both structural equations
reg3 (y1 y2 x1) (y2 y1 x2)

********** CHAPTER 7.4: IV EXAMPLE

* Read data, define global x2list, and summarize data
clear all
qui use mus207mepspresdrugs
global x2list totchr age female blhisp linc 
summarize ldrugexp hi_empunion $x2list

* Summarize available instruments 
summarize ssiratio lowincome multlc firmsz if linc!=.

* IV estimation of a just-identified model with single endog regressor
ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust) first

* Compare five estimators and variance estimates for overidentified models
global ivmodel "ldrugexp (hi_empunion = ssiratio multlc) $x2list"
qui ivregress 2sls $ivmodel, vce(robust)
estimates store TwoSLS
qui ivregress gmm  $ivmodel, wmatrix(robust) 
estimates store GMM_het
qui ivregress gmm  $ivmodel, wmatrix(robust) igmm
estimates store GMM_igmm
qui ivregress gmm  $ivmodel, wmatrix(cluster age) 
estimates store GMM_clu
qui ivregress 2sls  $ivmodel
estimates store TwoSLS_def
estimates table TwoSLS GMM_het GMM_igmm GMM_clu TwoSLS_def, b(%9.5f) se  

* Obtain OLS estimates to compare with preceding IV estimates
regress ldrugexp hi_empunion $x2list, vce(robust) 

* Robust DWH test of endogeneity implemented by estat endogenous
ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust)
estat endogenous

* Robust DWH test of endogeneity implemented manually
qui regress hi_empunion ssiratio $x2list
qui predict v1hat, resid
qui regress ldrugexp hi_empunion $x2list v1hat, vce(robust)
test v1hat 

* Control function estimator adds first-stage residual as regressors
regress ldrugexp hi_empunion $x2list v1hat, vce(robust) noheader

* Test of overidentifying restrictions following ivregress gmm
qui ivregress gmm ldrugexp (hi_empunion = ssiratio multlc) $x2list, ///
    wmatrix(robust) 
estat overid

* Test of overidentifying restrictions following ivregress gmm
ivregress gmm ldrugexp (hi_empunion = ssiratio lowincome multlc firmsz) ///
    $x2list, wmatrix(robust) 
estat overid

ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust)

// Not included 
* IV estimation using the eregress command in a just-identified model
eregress ldrugexp $x2list, endogenous(hi_empunion = $x2list ssiratio) vce(robust)
ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust)
ivregress liml ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust)
* Stata 15 has replaced the treatreg command with etregress

* Regression with a binary endogenous regressor
etregress ldrugexp $x2list, treat(hi_empunion = ssiratio $x2list) vce(robust)

********** CHAPTER 7.5: WEAK INSTRUMENTS

* Example is based on Chapter 8 of T. Lancaster
* Introduction to Modern Bayesian Econometrics (Blackwell, 2007)

* Program for weak instruments simulation 
clear all
global numsims 1000 
program weakivsim, rclass 
    version 17
    drop _all
    set obs $numobs                      // Set sample size
    generate z = rnormal() 
    generate v = rnormal()
    generate u = rnormal() + $ecor*v     // corr(u,v)=ecor/sqrt(1+ecor^2)
    generate y2 = 0 + $pi*z + v          // Set instrument strength
    generate y1 = 0 + 2*$pi*z + (u+2*v)  // beta = 2
    regress y1 y2
    return scalar bols = _b[y2] 
    return scalar seols = _se[y2]
    return scalar tols = (_b[y2]-2)/_se[y2]
    ivregress 2sls y1 (y2=z) 
    return scalar biv = _b[y2]
    return scalar seiv = _se[y2] 
    return scalar tiv = (_b[y2]-2)/_se[y2]
    regress y2 z                         // First-stage regression 
    return scalar Fiv = e(F)
end 

* Simulation 1: beta=2, large sample, weak IV, independent errors 
global numobs 1000   // Large sample 
global pi 0.1        // Weak IV
global ecor 0.0      // correlation(u,v)=0
simulate bols=r(bols) seols=r(seols) tols=r(tols) biv=r(biv) seiv=r(seiv) ///
    tiv=r(tiv) Fiv=r(Fiv), seed(10101) reps($numsims) nolegend nodots: weakivsim
mean bols seols biv seiv Fiv
display "Concentration parameter = " $numobs*$pi*1*$pi  // As E[z^2]=1

* Simulation 2: beta=2, small sample, weak IV, independent errors 
global numobs 100    // Small sample 
global pi 0.1        // Weak IV
global ecor 0.0      // correlation(u,v)=0
simulate bols=r(bols) seols=r(seols) tols=r(tols) biv=r(biv) seiv=r(seiv) ///
    tiv=r(tiv) Fiv=r(Fiv), seed(10101) reps($numsims) nolegend nodots: weakivsim
mean bols seols biv seiv Fiv
display "Concentration parameter = " $numobs*$pi*1*$pi  // As E[z^2]=1

* Simulation 3: beta=2, small sample, strong IV, independent errors 
global numobs 100   // Small sample 
global pi 0.5       // Strong IV
global ecor 0.0     // correlation(u,v)=00 
simulate bols=r(bols) seols=r(seols) tols=r(tols) biv=r(biv) seiv=r(seiv) ///
    tiv=r(tiv) Fiv=r(Fiv), seed(10101) reps($numsims) nolegend nodots: weakivsim
mean bols seols biv seiv Fiv 
display "Concentration parameter = " $numobs*$pi*1*$pi  // As E[z^2]=1

* Simulation 4: beta=2, large sample, weak IV, correlated errors 
global numobs 1000  // Large sample 
global pi 0.1       // Weak IV
global ecor 0.5775  // Implies correlation(u,v)=.5775/sqrt(1+.5775^2)=0.5
simulate bols=r(bols) seols=r(seols) tols=r(tols) biv=r(biv) seiv=r(seiv) ///
    tiv=r(tiv) Fiv=r(Fiv), seed(10101) reps($numsims) nolegend nodots: weakivsim
mean bols seols biv seiv Fiv 
display "Concentration parameter = " $numobs*$pi*1*$pi  // As E[z^2]=1

* Not included 
* Simulation 4a: beta=2, small sample, weak IV, correlated errors 
global numobs 100   // Small sample 
global pi 0.1       // Weak IV
global ecor 0.5775  // Implies correlation(u,v)=.5775/sqrt(1+.5775^2)=0.5
simulate bols=r(bols) seols=r(seols) tols=r(tols) biv=r(biv) seiv=r(seiv) ///
    tiv=r(tiv) Fiv=r(Fiv), seed(10101) reps($numsims) nolegend nodots: weakivsim
mean bols seols biv seiv Fiv
display "Concentration parameter = " $numobs*$pi*1*$pi  // As E[z^2]=1

* Simulation 5: beta=2, small sample, strong IV, correlated errors 
global numobs 100   // Small sample 
global pi 0.5       // Strong IV
global ecor 0.5775  // Implies correlation(u,v)=.5775/sqrt(1+.5775^2)=0.5
simulate bols=r(bols) seols=r(seols) tols=r(tols) biv=r(biv) seiv=r(seiv) ///
    tiv=r(tiv) Fiv=r(Fiv), seed(10101) reps($numsims) nolegend nodots: weakivsim
mean bols seols biv seiv Fiv
display "Concentration parameter = " $numobs*$pi*1*$pi  // As E[z^2]=1

/* Following is not used and takes a long time
* Cauchy program
clear all
program cauchysim, rclass 
    version 17
    drop _all
    set obs 1000000                      // Set sample size
    generate t = rt(1) 
    sum t
    return scalar mu = r(mean) 
    return scalar sd = r(sd)
end 
simulate mu=r(mu) sd=r(sd), seed(10101) reps(1000) nolegend nodots: cauchysim 
mean mu sd
*/

********** CHAPTER 7.6: DIAGNOSTICS AND TESTS FOR WEAK INSTRUMENTS

* Correlations of endogenous regressor with instruments
qui use mus207mepspresdrugs, clear
correlate hi_empunion ssiratio lowincome multlc firmsz if linc!=.
display "Concentration parameter = " $numobs*2*$pi/(2^2+2*2*$ecor+1)

* Weak instrument tests - just-identified with heteroskedastic-robust errors
global x2list totchr age female blhisp linc 
qui ivregress 2sls ldrugexp (hi_empunion = ssiratio) $x2list, vce(robust)
estat firststage, forcenonrobust all  

* Weak instrument tests - just-identified using weakivtest
qui ivregress 2sls ldrugexp (hi_empunion = ssiratio) ///
    $x2list, vce(robust)
weakivtest

* Weak instrument tests - overidentified with heteroskedastic-robust errors
qui ivregress 2sls ldrugexp (hi_empunion = ssiratio lowincome multlc firmsz) ///
    $x2list, vce(robust)
estat firststage, forcenonrobust

weakivtest

* Not included - compare weakivtest to estat first stage
qui ivregress gmm ldrugexp (hi_empunion = ssiratio) $x2list
estat firststage, forcenonrobust all
qui ivregress gmm ldrugexp (hi_empunion = ssiratio multlc) $x2list
weakivtest
qui ivregress gmm ldrugexp (hi_empunion = lowincome firmsz)  $x2list
weakivtest

* Not included - cluster–robust standard errors example
qui ivregress gmm ldrugexp (hi_empunion = ssiratio lowincome multlc firmsz) ///
    $x2list, vce(robust)
weakivtest
clear mata     // Prevents conflict with one of the subsequent commands

* ivreg2 for overidentified model with heteroskedastic-robust standard errors
ivreg2 ldrugexp (hi_empunion = ssiratio lowincome multlc firmsz ) $x2list, ///
    robust first

// Not included
* ivreg2 for just-identified model with heteroskedastic-robust standard errors
ivreg2 ldrugexp (hi_empunion = ssiratio) $x2list, robust sfirst

* Compare four just-identified model estimates with different instruments
qui regress ldrugexp hi_empunion $x2list, vce(robust)
estimates store OLS0
qui ivreg2 ldrugexp (hi_empunion=ssiratio) $x2list, robust
estimates store IV_INST1
scalar f1 = e(widstat)
qui ivreg2 ldrugexp (hi_empunion=lowincome) $x2list, robust
estimates store IV_INST2
scalar f2 = e(widstat)
qui ivreg2 ldrugexp (hi_empunion=multlc) $x2list, robust 
estimates store IV_INST3
scalar f3 = e(widstat)
qui ivreg2 ldrugexp (hi_empunion=firmsz) $x2list, robust
estimates store IV_INST4
scalar f4 = e(widstat)
estimates table OLS0 IV_INST1 IV_INST2 IV_INST3 IV_INST4, b(%8.4f) se  
display "Robust first-stage F:  " f1 _s(2) f2 _s(2) f3 _s(2) f4

* Montiel-Olea critical value for 10% relative bias with alpha = 0.05
qui ivregress 2sls ldrugexp (hi_empunion=firmsz) $x2list, vce(robust)
qui weakivtest
di "Effective F = " r(F_eff) " with critical value = " r(c_TSLS_10)
clear mata     

********** CHAPTER 7.7: BETTER INFERENCE WITH WEAK INSTRUMENTS

qui use mus207mepspresdrugs
global x2list totchr age female blhisp linc 

* Anderson-Rubin test for overidentified sample with robust standard errors
regress ldrugexp ssiratio firmsz $x2list, vce(robust)
test ssiratio firmsz

* Anderson-Rubin test for beta = -0.6 rather than beta = 0
qui generate yfordiffbeta = ldrugexp - (-0.6)*hi_empunion
regress yfordiffbeta ssiratio firmsz $x2list, vce(robust)
test ssiratio firmsz

// Not included
* Can also difference out exogenous regressors
* Anderson-Rubin test for overidentified eample with robust standard errors
qui regress ldrugexp $x2list
qui predict resldrugexp, resid
qui regress ssiratio $x2list
qui predict resssiratio, resid
qui regress firmsz $x2list
qui predict resfirmsz, resid
regress resldrugexp resssiratio resfirmsz, vce(robust)
* Anderson-Rubin test for beta = -0.6 rather than beta = 0
qui generate yfordiffbeta2 = resldrugexp - (-0.6)*hi_empunion
regress yfordiffbeta2 resssiratio resfirmsz, vce(robust)

* Condivreg: Weak instrument robust inference for i.i.d. errors
condivreg ldrugexp (hi_empunion=ssiratio firmsz) $x2list, 2sls lm ar test(0)

* rivtest: Weak instrument robust inference for non-i.i.d. errors
qui ivregress 2sls ldrugexp (hi_empunion=ssiratio firmsz) $x2list, vce(robust)
rivtest, ci null(0) 

* weakiv: Weak instrument robust inference for non-i.i.d. errors
qui ivregress 2sls ldrugexp (hi_empunion=ssiratio firmsz) $x2list, ///
    vce(robust)
weakiv, null(0)

********** CHAPTER 7.8: OTHER ESTIMATORS

* Variants of IV estimators: 2SLS, LIML, JIVE, GMM_het, GMM-het using IVREG2
global ivmodel "ldrugexp $x2list (hi_empunion = ssiratio lowincome multlc firmsz)"
qui ivregress 2sls $ivmodel, vce(robust)
estimates store TWOSLS
qui ivregress liml $ivmodel, vce(robust)
estimates store LIML
qui jive $ivmodel, robust
estimates store JIVE
qui ivregress gmm $ivmodel, wmatrix(robust) 
estimates store GMM_het
qui ivreg2 $ivmodel, gmm2s robust
estimates store IVREG2
estimates table TWOSLS LIML JIVE GMM_het IVREG2, b(%7.4f) se 

// Not included: eregress in an over-identified model equals LIML 
eregress ldrugexp $x2list, ///
    endogenous(hi_empunion = $x2list ssiratio lowincome multlc firmsz)
  
********** CHAPTER 7.9: 3SLS SYSTEMS ESTIMATION

* 3SLS estimation requires errors to be homoskedastic
reg3 (ldrugexp hi_empunion totchr age female blhisp linc) ///
    (hi_empunion ldrugexp totchr female blhisp ssiratio)

********** END
