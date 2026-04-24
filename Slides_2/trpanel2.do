* trpanel2.do  September 2013 for Stata 12.0
* Based on mus09p1panlin2.do

********** OVERVIEW OF trpanel2.do **********

log using trpanel2.txt, text replace

* Based on "Microeconometrics using Stata" 
* by A. Colin Cameron and Pravin K. Trivedi (2008)

* Chapter 9
* 9.2: PANEL INSTRUMENTAL VARIABLE ESTIMATORS
* 9.3: HAUSMAN TAYLOR ESTIMATOR
* 9.4: ARELLANO-BOND ESTIMATOR 
* 9.5: MIXED LINEAR MODELS
* 9.6: CLUSTERED DATA

* To run you need files
*   mus08psidextract.dta
*   abdata.dta
*   mus09vietnam_ex2.dta
* in your directory

* Stata user-written commands 
*  pvar, sgmm and helm
* are used

********** SETUP **********

set more off
version 12.0
clear all
set linesize 90
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus08psidextract.dta
* The Panel Study of Income Dynamics
* Same as Stata website file psidextract.dta
* Data due to  Baltagi and Khanti-Akom (1990) 
* This is corrected version of data in Cornwell and Rupert (1988).
* 595 individuals for years 1976-82.

* mus09vietnam_ex2.dta
* World Bank Vietnam Living Standards survey 1997-98. 
* Data from Cameron and Trivedi (2005, p.848)
* 5740 individuals

******* PANEL IV ESTIMATOR

* Panel IV example: FE with wks instrumented by external instrument ms
use mus08psidextract.dta, clear
xtivreg lwage exp exp2 (wks = ms), fe 

******* HAUSMAN-TAYLOR ESTIMATOR

* Hausman-Taylor example of Baltagi and Khanti-Akom (1990)
use mus08psidextract.dta, clear
xthtaylor lwage occ south smsa ind exp exp2 wks ms union fem blk ed,  ///
  endog(exp exp2 wks ms union ed)

// Hausman-Taylor with panel bootstrap SEs or with jackknife
// to get panel robust standard errors
* xthtaylor lwage occ south smsa ind exp exp2 wks ms union fem blk ed, ///
*   endog(exp exp2 wks ms union ed) vce(boot, reps(400) nodots seed(10101))
* xthtaylor lwage occ south smsa ind exp exp2 wks ms union fem blk ed, ///
*   endog(exp exp2 wks ms union ed) vce(jackknife)

******* ARELLANO-BOND ESTIMATOR

* 2SLS or one-step GMM for a pure time-series AR(2) panel model
use mus08psidextract.dta, clear
xtabond lwage, lags(2) vce(robust)

* Optimal or two-step GMM for a pure time-series AR(2) panel model
xtabond lwage, lags(2) twostep vce(robust)

* Reduce the number of instruments for a pure time-series AR(2) panel model
xtabond lwage, lags(2) vce(robust) maxldep(1)

* Optimal or two-step GMM for a dynamic panel model
xtabond lwage occ south smsa ind, lags(2) maxldep(3)     ///
  pre(wks,lag(1,2)) endogenous(ms,lag(0,2))              ///
  endogenous(union,lag(0,2)) twostep vce(robust) artests(3)

* Test whether error is serially correlated
estat abond

* Test of overidentifying restrictions (first estimate with no vce(robust))
quietly xtabond lwage occ south smsa ind, lags(2) maxldep(3) ///
  pre(wks,lag(1,2)) endogenous(ms,lag(0,2))              ///
  endogenous(union,lag(0,2)) twostep artests(3)
estat sargan

* Arellano/Bover or Blundell/Bond for a dynamic panel model
xtdpdsys lwage occ south smsa ind, lags(2) maxldep(3)    ///
  pre(wks,lag(1,2)) endogenous(ms,lag(0,2))              ///
  endogenous(union,lag(0,2)) twostep vce(robust) artests(3)

******* EXTENSION TO SHORT PANEL VAR MODEL WITH FIXED EFFECTS

* Use Arellano Bond data from Stata in file abdata.dta
* n is log employment;  k is log capital;  ys is log output;  w is log wage
* This requires user-written addons: pvar, sgmm and helm
* Written by Love
* Article: Love, I. and Ziccino, L. (2006), "Financial Development and 
* Dynamic Investment Behaviour: Evidence from Panel VAR," 
* Quarterly Review of Economics and Finance, 46, 190-210.
use abdata.dta, clear
summarize
xtset
* Do Helmert transformations
helm n 
helm k
helm ys
helm w
* Just two variables
pvar n k, lag(2) gmm impulse monte 50 decomp list_imp
graph export ct_panel24k.wmf, replace
* Now three variables
pvar n k ys, lag(2) gmm impulse monte 50 decomp list_imp


******* MIXED LINEAR MODELS - NOTE THAT STATA 13 HAS CHANGED THIS

* Random intercept model estimated using xtmixed
use mus08psidextract.dta, clear
xtmixed lwage exp exp2 wks ed || id:, mle

* Cluster robust standard errors after xtmixed using bootstrap
* xtset id
* bootstrap, reps(400) seed(10101) cluster(id) nodots: ///
*   xtmixed lwage exp exp2 wks ed || id:, mle

/* Above bootstrap commented out as takes a long time 
   It yields output
. * Cluster robust standard errors after xtmixed using bootstrap
. xtset id
       panel variable:  id (balanced)
. bootstrap, reps(400) seed(10101) cluster(id) nodots: ///
>   xtmixed lwage exp exp2 wks ed || id:, mle

Mixed-effects ML regression                     Number of obs      =      4165
Group variable: id                              Number of groups   =       595
                                                Obs per group: min =         7
                                                               avg =       7.0
                                                               max =         7
                                                Wald chi2(4)       =   2092.79
Log likelihood =  293.69563                     Prob > chi2        =    0.0000
                                    (Replications based on 595 clusters in id)
------------------------------------------------------------------------------
             |   Observed   Bootstrap                         Normal-based
       lwage |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
         exp |   .1079955   .0041447    26.06   0.000     .0998721    .1161189
        exp2 |  -.0005202   .0000831    -6.26   0.000    -.0006831   -.0003573
         wks |   .0008365   .0008458     0.99   0.323    -.0008212    .0024943
          ed |   .1378559   .0099856    13.81   0.000     .1182844    .1574273
       _cons |   2.989858   .1510383    19.80   0.000     2.693829    3.285888
------------------------------------------------------------------------------

------------------------------------------------------------------------------
                             |   Observed   Bootstrap         Normal-based
  Random-effects Parameters  |   Estimate   Std. Err.     [95% Conf. Interval]
-----------------------------+------------------------------------------------
id: Identity                 |
                   sd(_cons) |   .8509015   .0259641      .8015045    .9033428
-----------------------------+------------------------------------------------
                sd(Residual) |   .1536109     .00824      .1382808    .1706406
------------------------------------------------------------------------------
LR test vs. linear regression: chibar2(01) =  4576.13 Prob >= chibar2 = 0.0000
*/

* Random-slopes model estimated using xtmixed
xtmixed lwage exp exp2 wks ed || id: exp wks, covar(unstructured) mle

* Random-coefficients model estimated using xtrc
quietly set matsize 600
* The following works in Stata 10 but 11, 12, 13
* xtrc lwage exp wks, i(id)
* Error: variables dropped due to collinearity. xtrc requires that you be able to fit 
* OLS regression to each panel without dropping variables due to collinearity.
* Problem: some people have the same hrs every year
* So instead
xtrc lwage exp, i(id)
 
* List the estimated variance matrix
matrix list e(Sigma)

* Two-way random-effects model estimated using xtmixed
xtmixed lwage exp exp2 wks ed || _all: R.t || id: , mle

// Following not included in book
* gllamm
* gllamm lwage exp exp2 wks ed, i(id) nip(10) adapt
* xtmixed lwage exp exp2 wks ed || id:, mle

/***** UPDATE FOR STATA 13  

* Stata 13 has replaced xtmixed with mixed
* and mixed now does cluster robust standard errors so is better

* The above program is written for Stata 12 
* and will run under Stata 13 under version 12.0 control

* In Stata 13 run the following code

* Random intercept model estimated using xtmixed
use mus08psidextract.dta, clear
mixed lwage exp exp2 wks ed || id:, mle

* For cluster-robust clustering on id either of the following works
mixed lwage exp exp2 wks ed || id:, mle vce(robust)
mixed lwage exp exp2 wks ed || id:, mle vce(cluster id)

* Random-slopes model estimated using mixed
mixed lwage exp exp2 wks ed || id: exp wks, covar(unstructured) mle
mixed lwage exp exp2 wks ed || id: exp wks, covar(unstructured) mle vce(robust)

* Two-way random-effects model estimated using xtmixed
mixed lwage exp exp2 wks ed || _all: R.t || id: , mle
* Robust and cluster vce's no longer available

*/

******* 9.6: CLUSTERED DATA

* Read in Vietnam clustered data and summarize
use mus09vietnam_ex2.dta, clear
summarize pharvis lnhhexp illness commune

* OLS estimation with cluster-robust standard errors
quietly regress pharvis lnhhexp illness
estimates store OLS_iid
quietly regress pharvis lnhhexp illness, vce(robust)
estimates store OLS_het
quietly regress pharvis lnhhexp illness, vce(cluster lnhhexp)
estimates store OLS_hh
quietly regress pharvis lnhhexp illness, vce(cluster commune) 
estimates store OLS_vill
estimates table OLS_iid OLS_het OLS_hh OLS_vill, b(%10.4f) se stats(r2 N)

* Generate integer-valued household and person identifiers and xtset
quietly egen hh = group(lnhhexp)
sort hh
by hh: generate person = _n
xtset hh person

drop if missing(hh)

xtdescribe

* Within-cluster correlation of pharvis
* Best is to use loneway
loneway pharvis hh
* The following is very similar
quietly xtreg pharvis, mle
display "Intra-class correlation for household: " e(rho)
* THe following is very crude
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
estimates table OLS_hh RE_hh FE_hh OLS_vill RE_vill FE_vill, b(%7.4f) se

* Hierarchical linear model with household and village variance components
xtmixed pharvis lnhhexp illness || commune: || hh:, mle difficult

/*
* STATA 13: Note it gives variances and not standard deviations
xtmixed pharvis lnhhexp illness || commune: || hh:, mle difficult
* And for cluster-robust clustering on commune
xtmixed pharvis lnhhexp illness || commune: || hh:, mle difficult vce(robust)
xtmixed pharvis lnhhexp illness || commune: || hh:, mle difficult cluster(commune)
* The following does not work as not same or higher level than commune
* xtmixed pharvis lnhhexp illness || commune: || hh:, mle difficult cluster(hh)
*/

********** CLOSE OUTPUT
* log close
* clear
* exit

