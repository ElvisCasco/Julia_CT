* cluster2017part2.do  based on cluster2015part2.do  August 2017 for Stata 14

log using cluster2017part2.txt, text replace

********** OVERVIEW OF cluster2017part2.do **********

* Based on "Microeconometrics using Stata" 
* by A. Colin Cameron and Pravin K. Trivedi (2008)

* Chapter 9
*  CLUSTERED DATA - MOULTON TYPE

* To run you need files
*   mus08psidextract.dta
*   mus09vietnam_ex2.dta
* in your directory

* Stata user-written commands 
*  ivreg2
*  cgmreg.ado is at http://cameron.econ.ucdavis.edu/research/papers.html
*  but is commented out so no need to use
* are used

********** SETUP **********

clear all
set more off
version 14.0
set linesize 90
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus09vietnam_ex2.dta
* World Bank Vietnam Living Standards survey 1997-98. 
* Data from Cameron and Trivedi (2005, p.848)
* 5740 individuals

* mus08psidextract.dta
* The Panel Study of Income Dynamics
* Same as Stata website file psidextract.dta
* Data due to  Baltagi and Khanti-Akom (1990) 
* This is corrected version of data in Cornwell and Rupert (1988).
* 595 individuals for years 1976-82.

******* CLUSTERED DATA - MOULTON TYPE CLUSTERING

* Read in Vietnam clustered data and summarize
* Persons are in households are in communes

use mus09vietnam_ex2.dta, clear

* One observation is missing so drop it
list in 27733 / 27735
drop if missing(pharvis)

describe
summarize

* Some key variables
summarize pharvis lnhhexp illness commune

* xtset is currently on commune
xtset 

* We also want household identifier
* Use fact that lnhhexp is unique to each household in these data
quietly egen hh = group(lnhhexp)
bysort hh: generate person_in_hh = _n
sum hh person_in_hh
* Generate unique person in commune identifier
generate person = 100*hh + person_in_hh

* Summary statistics
sum pharvis lnhhexp illness AGE hh person_in_hh person

* OLS estimation with cluster-robust standard errors
* Cluster on household and then on commume
quietly regress pharvis lnhhexp illness
estimates store OLS_iid
quietly regress pharvis lnhhexp illness, vce(robust)
estimates store OLS_het
quietly regress pharvis lnhhexp illness, vce(cluster hh)
estimates store OLS_hh
quietly regress pharvis lnhhexp illness, vce(cluster commune) 
estimates store OLS_comm
estimates table OLS_iid OLS_het OLS_hh OLS_comm, b(%10.4f) se stats(r2 N)

* Cluster on hh
xtset hh person_in_hh
xtdescribe

* If instead cluster on commune and person
* xtset commune person

* Within-cluster correlation of pharvis
loneway pharvis hh

quietly xtreg pharvis, mle
display "Intra-class correlation for household: " e(rho)
quietly correlate pharvis L1.pharvis
display "Correlation for adjoining household:   " r(rho)

* OLS, RE and FE estimation with clustering on household and on village
quietly regress pharvis lnhhexp illness, vce(cluster hh)
estimates store OLS_hh
quietly xtreg pharvis lnhhexp illness, re
estimates store RE_hh
quietly xtreg pharvis lnhhexp illness, fe
estimates store FE_hh
quietly xtset commune
quietly regress pharvis lnhhexp illness, vce(cluster commune)
estimates store OLS_vill
quietly xtreg pharvis lnhhexp illness, re
estimates store RE_vill
quietly xtreg pharvis lnhhexp illness, fe
estimates store FE_vill
estimates table OLS_hh RE_hh FE_hh OLS_vill RE_vill FE_vill, b(%7.5f) se(%7.4f)

* Simplest case of random intercept same as xtreg, mle
mixed pharvis lnhhexp illness || hh:, nolog 
xtreg pharvis lnhhexp illness, mle

* Mixed model with random interecept and a random slope
mixed pharvis lnhhexp illness || hh: illness, nolog covariance(unstructured)

/* Comment out as takes time
* Twoway random effects with error e_g + e_h + e_igh, g=ILLDAYS h=hh
mixed pharvis lnhhexp illness || _all: R.ILLDAYS || hh: , mle
*/

******* TWOWAY CLUSTER ROBUST

* Twoway cluster-robust example - use addon cgmreg.ado

/*
* Twoway cluster-robust example - use addon cgmreg.ado
cgmreg pharvis lnhhexp illness, cluster(hh AGE)

Note: +/- means the corresponding matrix is added/subtracted

Calculating cov part for variables:  hh (+)
Calculating cov part for variables:  hh AGE (-)
Calculating cov part for variables:  AGE (+)
                                                 Number of obs     =    27765
                                                 Num clusvars      =    2
                                                 Num combinations  =    3
                                                 G(hh)             =    5740
                                                 G(AGE)            =    98
------------------------------------------------------------------------------
     pharvis |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
     lnhhexp |   .0247682   .0139228     1.78   0.075      -.00252    .0520563
     illness |    .624163   .0190994    32.68   0.000     .5867288    .6615971
       _cons |   .0590868   .0382208     1.55   0.122    -.0158245    .1339981
------------------------------------------------------------------------------
*/

* Twoway cluster-robust example - use add on ivreg2
ivreg2 pharvis lnhhexp illness, cluster(hh AGE)

************ PAIRS CLUSTER PERCENTILE t

* Reduce sample size to speed up bootstraps 
preserve
keep if commune <= 20

regress pharvis lnhhexp illness, vce(cluster commune)

* For cluster bootstrap can only xtset the cluster variable
xtset commune

* Following gives BC and BCa cluster bootstrap with asymptotic refinement
regress pharvis lnhhexp illness, ///
   vce(boot, cluster(commune) seed(10101)  reps(999) bca)
estat bootstrap, all

* Percentile-t for a single coefficient: Bootstrap the t statistic
quietly regress pharvis lnhhexp illness, vce(cluster commune)
local theta = _b[lnhhexp] 
local setheta = _se[lnhhexp]
bootstrap tstar=((_b[lnhhexp]-`theta')/_se[lnhhexp]), seed(10101) ///
  reps(999) saving(percentilet, replace):  ///
  regress pharvis lnhhexp illness, vce(cluster commune)

* Percentile-t p-value for symmetric two-sided Wald test of H0: theta = 0
use percentilet, clear
quietly count if abs(`theta'/`setheta') < abs(tstar)
display "p-value = " r(N)/_N

* Percentile-t critical values and confidence interval
_pctile tstar, p(2.5,97.5) 
scalar lb = `theta' + r(r1)*`setheta'
scalar ub = `theta' + r(r2)*`setheta'
display "2.5 and 97.5 percentiles of t* distn: " r(r1) ", " r(r2) _n ///
    "95 percent percentile-t confidence interval is:  (" lb   ","  ub ")"
 
restore

********** CLOSE OUTPUT
* log close
* clear
* exit

