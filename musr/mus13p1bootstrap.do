* mus13p1bootstrap.do  Oct 2009 for Stata version 11

cap log close

********** OVERVIEW OF mus13p1bootstrap.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Chapter 13
* 13.3: BOOTSTRAP PAIRS USING THE VCE(BOOTSTRAP) OPTION
* 13.4: BOOTSTRAP PAIRS USING THE BOOTSTRAP COMMAND
* 13.5: BOOTSTRAPS WITH ASYMPTOTIC REFINEMENT
* 13.6: BOOTSTRAP PAIRS USING BSAMPLE AND SIMULATE
* 13.7: ALTERNATIVE RESAMPLING SCHEMES

* To run you need files
*   mus10data.dta    
*   mus16data.dta 
* in your directory
* No Stata user-written commands are used

* To speed up program reduce reps in 
*     vce(boot, reps(400)) and vce(boot, reps(999))
* or in bootstrap reps(999)   

********** SETUP **********

set more off
version 11
clear all
set linesize 81
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* This program analyzes the data used in chapter 10
* 2002 Medical Expenditure Panel Survey (MEPS)
* U.S. individuals aged 25-64 years working in private sector 
* but not self-employed
* and not receiving public insurance (Medicare and Medicaid) 
* Data due to Deb, Munkin and Trivedi (2006)

****** 13.3: BOOTSTRAP PAIRS USING THE VCE(BOOTSTRAP) OPTION

* Sample is only the first 50 observations of chapter 10 data
use mus10data.dta
quietly keep if year02 == 1
quietly drop if _n > 50
quietly keep docvis chronic age 
quietly save bootdata, replace

* Option vce(bootstrap) to compute bootstrap standard errors
poisson docvis chronic, vce(boot, reps(400) seed(10101) nodots)

* Bootstrap standard errors for different reps and seeds
quietly poisson docvis chronic, vce(boot, reps(50) seed(10101))
estimates store boot50
quietly poisson docvis chronic, vce(boot, reps(50) seed(20202))
estimates store boot50diff
quietly poisson docvis chronic, vce(boot, reps(2000) seed(10101))
estimates store boot2000
quietly poisson docvis chronic, vce(robust)
estimates store robust
estimates table boot50 boot50diff boot2000 robust, b(%8.5f) se(%8.5f)

* Option vce(boot, cluster) to compute cluster-bootstrap standard errors
poisson docvis chronic, vce(boot, cluster(age) reps(400) seed(10101) nodots)

poisson docvis chronic, vce(cluster age)

* Bootstrap confidence intervals: normal-based, percentile, BC, and BCa
quietly poisson docvis chronic, vce(boot, reps(999) seed(10101) bca)
estat bootstrap, all

* List the average of the bootstraps
matrix list e(b_bs)

****** 13.4: BOOTSTRAP PAIRS USING THE BOOTSTRAP COMMAND

* Bootstrap command applied to Stata estimation command
bootstrap, reps(400) seed(10101) nodots noheader: poisson docvis chronic

* Bootstrap standard-error estimate of the standard error of a coeff estimate
bootstrap _b _se, reps(400) seed(10101) nodots: poisson docvis chronic

* Program to return b and robust estimate V of the VCE
program poissrobust, eclass
  version 11
  tempname b V
  poisson docvis chronic, vce(robust)
  matrix `b' = e(b)
  matrix `V' = e(V)
  ereturn post `b' `V'
end

// Following output not given in the book 
* Check preceding program by running once
poissrobust
ereturn display

* Bootstrap standard-error estimate of robust standard errors
bootstrap _b _se, reps(400) seed(10101) nodots nowarn: poissrobust

* Set up the selection model two-step estimator data of chapter 16
use mus16data.dta, clear
generate y = ambexp
generate dy = y > 0
generate lny = ln(y)
global xlist age female educ blhisp totchr ins 

* Program to return b for Heckman 2-step estimator of selection model 
program hecktwostep, eclass
  version 11
  tempname b V
  tempvar xb 
  capture drop invmills
  probit dy $xlist
  predict `xb', xb
  generate invmills = normalden(`xb')/normprob(`xb')
  regress lny $xlist invmills
  matrix `b' = e(b)
  ereturn post `b' 
end

// Following not included in book
* Check preceding program by running once
hecktwostep
ereturn display

* Bootstrap for Heckman two-step estimator using chapter 16 example
bootstrap _b, reps(400) seed(10101) nodots nowarn: hecktwostep

// Following not included in book
* Check results
heckman lny $xlist, select(dy = $xlist) twostep

* Program to return (b1-b2) for hausman test of endogeneity
program hausmantest, eclass
  version 11
  tempname b bols biv
  regress ldrugexp hi_empunion totchr age female blhisp linc, vce(robust) 
  matrix `bols' = e(b)
  ivregress 2sls ldrugexp (hi_empunion = ssiratio) totchr age female blhisp linc, vce(robust)  
  matrix `biv' = e(b)
  matrix `b' = `bols' - `biv'
  ereturn post `b' 
end

// Following not included in book
* Check preceding program by running once
use mus06data.dta, clear
hausmantest
ereturn display

* Bootstrap estimates for Hausman test using chapter 6 example
use mus06data.dta, clear
bootstrap _b, reps(400) seed(10101) nodots nowarn: hausmantest

* Perform Hausman test on the potentially endogenous regressor
test hi_empunion 

// Following not included in book
* Perform Hausman test on all regressors
test hi_empunion totchr age female blhisp linc

* Bootstrap estimate of the standard error of the coefficient of variation
use bootdata, clear
bootstrap coeffvar=(r(sd)/r(mean)), reps(400) seed(10101) nodots   ///
   nowarn saving(coeffofvar, replace): summarize docvis

****** 13.5: BOOTSTRAPS WITH ASYMPTOTIC REFINEMENT

* Percentile-t for a single coefficient: Bootstrap the t statistic
use bootdata, clear
quietly poisson docvis chronic, vce(robust)
local theta = _b[chronic] 
local setheta = _se[chronic]
bootstrap tstar=((_b[chronic]-`theta')/_se[chronic]), seed(10101)        ///
  reps(999) nodots saving(percentilet, replace): poisson docvis chronic, ///
  vce(robust)

* Percentile-t p-value for symmetric two-sided Wald test of H0: theta = 0
use percentilet, clear
quietly count if abs(`theta'/`setheta') < abs(tstar)
display "p-value = " r(N)/_N

* Percentile-t critical values and confidence interval
_pctile tstar, p(2.5,97.5) 
scalar lb = `theta' + r(r1)*`setheta'
scalar ub = `theta' + r(r2)*`setheta'
display "2.5 and 97.5 percentiles of t* distn: " r(r1) ", " r(r2) _n ///
    "95 percent percentile-t confidence interval is  (" lb ","  ub ")"

use bootdata, clear
nbreg docvis chronic

******* 13.6: BOOTSTRAP PAIRS USING BSAMPLE AND SIMULATE

* Program to do one bootstrap replication
program onebootrep, rclass
  version 11  
  drop _all
  use bootdata
  bsample
  poisson docvis chronic, vce(robust)
  return scalar tstar = (_b[chronic]-$theta)/_se[chronic]
end

* Now do 999 bootstrap replications
use bootdata, clear
quietly poisson docvis chronic, vce(robust) 
global theta = _b[chronic]
global setheta = _se[chronic]
simulate tstar=r(tstar), seed(10101) reps(999) nodots  ///
  saving(percentilet2, replace): onebootrep

* Analyze the results to get the p-value
use percentilet2, clear
quietly count if abs($theta/$setheta) < abs(tstar)
display "p-value = " r(N)/_N  

* Program to do one bootstrap of B replications
program mybootstrap, rclass
  use bootdata, clear
  quietly poisson docvis chronic, vce(robust) 
  global theta = _b[chronic]
  global setheta = _se[chronic]
  simulate tstar=r(tstar), reps(999) nodots  ///
    saving(percentilet2, replace): onebootrep
  use percentilet2, clear
  quietly count if abs($theta/$setheta) < abs(tstar)
  return scalar pvalue =  r(N)/_N 
end 

* Check the program by running once
set seed 10101 
mybootstrap
display r(pvalue)

******* 13.7 ALTERNATIVE RESAMPLING SCHEMES (TO GET SE'S)

* Program to resample using bootstrap pairs
program bootpairs
  version 11 
  drop _all
  use bootdata
  bsample
  poisson docvis chronic
end

* Check the program by running once
bootpairs

* Bootstrap pairs for the parameters
simulate _b, seed(10101) reps(400) nodots: bootpairs
summarize

// Following not included in book
* Compare to bootstrap command
use bootdata, clear
bootstrap, reps(100) seed(10101): poisson docvis chronic

* Estimate the model with original actual data and save estimates
use bootdata
nbreg docvis chronic
predict muhat
global alpha = e(alpha)

* Program for parametric bootstrap generating from negative binomial
program bootparametric, eclass
  version 11
  capture drop nu dvhat
  generate nu = rgamma(1/$alpha,$alpha)
  generate dvhat = rpoisson(muhat*nu)
  nbreg dvhat chronic
end

// Following not included in book 
* Check the program by running once
set seed 10101
bootparametric

* Parametric bootstrap for the parameters
simulate _b, seed(10101) reps(400) nodots: bootparametric
summarize

* Program for residual bootstrap for OLS with iid errors
use bootdata, clear
quietly regress docvis chronic
predict uhat, resid
keep uhat
save residuals, replace
program bootresidual
  version 11 
  drop _all
  use residuals
  bsample                     
  merge using bootdata   
  regress docvis chronic
  predict xb
  generate ystar =  xb + uhat
  regress ystar chronic
end

// Following not included in book 
* Check the program by running once
bootresidual

* Residual bootstrap for the parameters
simulate _b, seed(10101) reps(400) nodots: bootresidual
summarize

* Wild bootstrap for OLS with iid errors
use bootdata, clear
program bootwild
  version 11 
  drop _all
  use bootdata 
  regress docvis chronic
  predict xb
  predict u, resid
  gen ustar = -0.618034*u
  replace ustar = 1.618034*u if runiform() > 0.723607
  gen ystar =  xb + ustar
  regress ystar chronic
end

// Following not included in book
* Check the program by running once
bootwild

* Wild bootstrap for the parameters
simulate _b, seed(10101) reps(400) nodots: bootwild
summarize

* Jackknife estimate of standard errors
use bootdata, replace
poisson docvis chronic, vce(jackknife, mse nodots) 

// Following not included in book
jacknife, mse nodots: poisson docvis chronic

********** CLOSE OUTPUT
