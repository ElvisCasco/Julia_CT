* ct_nonparboot.do  My 2013 for Stata version 12.0

********** OVERVIEW OF ct_nonparboot.do **********

log using ct_nonparboot.txt, text replace

* STATA Program by A. Colin Cameron
* Kernel density
* Nonparametric
* Bootstrap 

* To run you need files
*   nonparametric.dta  
*   bootstrap.dta    
* in your directory

*** THE BOOTSTRAPS TAKE TIME
*** TO SPEED UP COMMENT THEM OUT OR CHANGE NUMBER OF REPS TO 50

********** SETUP **********

set more off
version 10.0
set mem 10m
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* The original data are from the PSID Individual Level Final Release 1993 data
* From www.isr.umich.edu/src/psid  then choose Data Center 
* 4856 observations on 9 variables for Females 30 to 50 years 
* See mma09p1np.do for further description

******* OLS WITH DOCTOR VISITS DATA

* Read in data, select, describe and summarize key variables
use nonparametric.dta, clear
describe 

* Work with age 36 and nonmissing education data
keep if age == 36
drop if educatn == .
summarize

******* KERNEL DENSITY ESTIMATE

histogram lnhwage
histogram lnhwage, bin(30)
graph save ct_nonparametric1, replace
graph export ct_nonparametric1.wmf, replace
kdensity lnhwage
kdensity lnhwage, bw(0.21)
graph twoway (kdensity lnhwage, bw(0.21)) (kdensity lnhwage, bw(0.07) clstyle(p2)) (kdensity lnhwage, bw(0.63) clstyle(p3)), legend( label(1 "Default") label(2 "Half default") label(3 "Twice default") ) 
graph save ct_nonparametric2, replace
graph export ct_nonparametric2.wmf, replace
histogram lnhwage, kdensity  
kdensity lnhwage, normal

******* NONPARAMETRIC REGRESSION 

regress lnhwage educatn
* Kernel
lpoly lnhwage educatn, ci msize(medsmall)
graph save ct_nonparametric3, replace
graph export ct_nonparametric3.wmf, replace
* Local linear
lpoly lnhwage educatn, degree(1) ci
* Lowess
lowess lnhwage educatn
* Kernel for different bandwidths - default, halfdefault, twicedefault
graph twoway (lpoly lnhwage educatn, bw(1.5)) (lpoly lnhwage educatn, bw(0.75) clstyle(p2)) (lpoly lnhwage educatn, bw(3.0) clstyle(p3)), legend( label(1 "Default") label(2 "Half default") label(3 "Twice default") ) 
graph save ct_nonparametric4, replace
graph export ct_nonparametric4.wmf, replace
* Compare kernel, local linear, lowess with default bandwidths
graph twoway (lpoly lnhwage educ) (lpoly lnhwage educ, degree(1) clstyle(p2)) (lowess lnhwage educ, clstyle(p3)), legend( label(1 "Kernel") label(2 "Local linear") label(3 "lowess") ) 
graph save ct_nonparametric5, replace
graph export ct_nonparametric5.wmf, replace

****** BOOTSTRAP PAIRS USING THE VCE(BOOTSTRAP) OPTION

* Data description
use bootdata.dta, clear 
describe 
summarize

* Default and robust standard errors
poisson docvis chronic, nolog
poisson docvis chronic, nolog vce(robust)

* Compute bootstrap standard errors using option vce(bootstrap) to 
poisson docvis chronic, vce(boot, reps(400) seed(10101) nodots)

* Same using bootstrap prefix command
bootstrap, seed(10101) reps(400) nodots: poisson docvis chronic

* The jackknife
poisson docvis chronic, vce(jackknife)

* Bootstrap confidence intervals: normal-based, percentile, BC, and BCa
quietly poisson docvis chronic, vce(boot, reps(999) seed(10101) bca)
estat bootstrap, all

/* SKIP THE FOLLOWING AS TAKES TIME

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

*/

* MAY WANT TO SKIP THE FOLLOWING AS TAKES TIME

****** PERCENTILE-T BOOTSTRAP WITH ASYMPTOTIC REFINEMENT

* Percentile-t for a single coefficient: Bootstrap the t statistic
use bootdata.dta, clear 
quietly poisson docvis chronic, vce(robust)
local theta = _b[chronic] 
local setheta = _se[chronic]
bootstrap tstar=((_b[chronic]-`theta')/_se[chronic]), seed(10101)        ///
  reps(999) nodots saving(percentilet, replace): poisson docvis chronic, ///
  vce(robust)

* Simple plot of the distribution of the tstar
use percentilet, clear
summarize
centile tstar, c(2.5,50,97.5)
kdensity tstar, normal

* Fancier plot of the distribution of the tstar
use percentilet, clear
summarize
centile tstar, c(2.5,50,97.5)
kdensity tstar, generate(evalpoint densityest) xtitle("tstar from the bootstrap replications")
generate phistnorm = normalden(evalpoint)
label variable phistnorm "Standard normal"
label variable densityest "Bootstrap density"
label variable evalpoint "tstar"
graph twoway (scatter densityest evalpoint, connect(l) msize(tiny)) ///
  (scatter phistnorm evalpoint, connect(l) msize(tiny))
graph export ct_bootstrap5.wmf, replace

* Percentile-t p-value for symmetric two-sided Wald test of H0: theta = 0
use percentilet, clear
quietly count if abs(`theta'/`setheta') < abs(tstar)
display "p-value = " r(N)/(_N+1)

* Percentile-t critical values and confidence interval
_pctile tstar, p(2.5,97.5) 
scalar lb = `theta' + r(r1)*`setheta'
scalar ub = `theta' + r(r2)*`setheta'
display "2.5 and 97.5 percentiles of t* distn: " r(r1) ", " r(r2) _n ///
    "95 percent percentile-t confidence interval is:  (" lb   ","  ub ")"

********** CLOSE OUTPUT
* log close
* clear
* exit


