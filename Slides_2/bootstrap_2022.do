* bootstrap_2022.do  February 2022 for Stata version 16

log using bootstrap_2022.txt, text replace

********** OVERVIEW OF bootstrap2020.do **********

* STATA Program by A. Colin Cameron
* Bootstrap 

* To run you need files
*   bootstrap.dta    
* in your directory

* And this uses Stata addon
*     boottest

********** SETUP **********

clear all
set more off
version 16
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* The original data are from the PSID Individual Level Final Release 1993 data
* From www.isr.umich.edu/src/psid  then choose Data Center 
* 4856 observations on 9 variables for Females 30 to 50 years 
* See mma09p1np.do for further description

****** BOOTSTRAP PAIRS USING THE VCE(BOOTSTRAP) OPTION

* Data description
use bootdata.dta, clear 
describe 

* Summmarize and Poisson with robust se's 
summarize
poisson docvis chronic, nolog vce(robust)

* Default standard errors are way too small
poisson docvis chronic, nolog

* Compute bootstrap standard errors using option vce(bootstrap) to 
poisson docvis chronic, vce(boot, reps(400) seed(10101) nodots)

* Same using bootstrap prefix command
bootstrap, seed(10101) reps(400) nodots: poisson docvis chronic

* The jackknife
poisson docvis chronic, vce(jackknife)

* Bootstrap confidence intervals: normal-based, percentile, BC, and BCa
quietly poisson docvis chronic, vce(boot, reps(999) seed(10101) bca)
estat bootstrap, all

* The following as takes time

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

****** PERCENTILE-T BOOTSTRAP WITH ASYMPTOTIC REFINEMENT

* Percentile-t for a single coefficient: Bootstrap the t statistic
use bootdata.dta, clear 
quietly poisson docvis chronic, vce(robust)
local theta = _b[chronic] 
local setheta = _se[chronic]
bootstrap tstar=((_b[chronic]-`theta')/_se[chronic]), seed(10101)        ///
  reps(999) nodots saving(percentilet, replace): poisson docvis chronic, ///
  vce(robust)

* Simple plot of the distribution of the 999 tstar
use percentilet, clear
summarize
centile tstar, c(2.5,50,97.5)
kdensity tstar, student(48) xtitle("tstar") title("Density of tsar") scale(1.1)
* graph export bootstrap05fig.wmf, replace

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
* graph export ct_bootstrap1.wmf, replace

* Percentile-t critical values and confidence interval
centile  tstar, c(2.5, 97.5)

* Percentile-t p-value for symmetric two-sided Wald test of H0: theta = 0
use percentilet, clear
quietly count if abs(`theta'/`setheta') < abs(tstar)
display "p-value = " r(N)/(_N+1)

******* WILD CLUSTER BOOTSTRAP 

* Use asymptotically incorrect default standard errors
use bootdata.dta, clear 
poisson docvis chronic 
boottest chronic, seed(10101)

* Use asympototically correct heteroskedastic-robust standard errors
use bootdata.dta, clear 
poisson docvis chronic, vce(robust)
boottest chronic, seed(10101)

********** CLOSE OUTPUT
* log close
* clear
* exit


