* mus09p1panlin2.do Oct 2009 for Stata version 11

cap log close

********** OVERVIEW OF mus09p1panlin2.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Chapter 9
* 9.2: PANEL INSTRUMENTAL VARIABLE ESTIMATORS
* 9.3: HAUSMAN TAYLOR ESTIMATOR
* 9.4: ARELLANO-BOND ESTIMATOR 
* 9.5: MIXED LINEAR MODELS
* 9.6: CLUSTERED DATA

* To run you need files
*   mus08psidextract.dta
*   mus09vietnam_ex2.dta
* in your directory
* No Stata user-written commands are used

********** SETUP **********

set more off
version 11
clear all
set memory 30m
set linesize 90
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus08psidextract.dta
* PSID. Same as Stata website file psidextract.dta
* Data due to  Baltagi and Khanti-Akom (1990) 
* This is corrected version of data in Cornwell and Rupert (1988).
* 595 individuals for years 1976-82

* mus09vietnam_ex2.dta
* World Bank Vietnam Living Standards survey 1997-98. 
* Data from Cameron and Trivedi (2005, p.848)
* 5740 individuals

******* 9.2: PANEL IV ESTIMATOR

* Panel IV example: FE with wks instrumented by external instrument ms
use mus08psidextract.dta, clear
xtivreg lwage exp exp2 (wks = ms), fe 

// Panel bootstrap SEs
* xtivreg lwage exp exp2 (wks = ms), fe vce(boot, reps(400) seed(10101))

******* 9.3: HAUSMAN-TAYLOR ESTIMATOR

* Hausman-Taylor example of Baltagi and Khanti-Akom (1990)
use mus08psidextract.dta, clear
xthtaylor lwage occ south smsa ind exp exp2 wks ms union fem blk ed,  ///
  endog(exp exp2 wks ms union ed)

// Hausman-Taylor with panel bootstrap SEs or with jackknife
* xthtaylor lwage occ south smsa ind exp exp2 wks ms union fem blk ed, ///
*  endog(exp exp2 wks ms union ed) vce(boot, reps(400) nodots seed(10101))
* xthtaylor lwage occ south smsa ind exp exp2 wks ms union fem blk ed, ///
*  endog(exp exp2 wks ms union ed) vce(jackknife)

******* 9.4: ARELLANO-BOND ESTIMATOR

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

// Following not included in book
estat abond

* Use of xtdpd to exactly reproduce the previous xtdpdsys command
xtdpd L(0/2).lwage L(0/1).wks occ south smsa ind ms union, ///
  div(occ south smsa ind) dgmmiv(lwage, lagrange(2 4))          ///
  dgmmiv(ms union, lagrange(2 3)) dgmmiv(L.wks, lagrange(1 2))      ///
  lgmmiv(lwage wks ms union) twostep vce(robust) artests(3)

* Previous command if model error is MA(1)
xtdpd L(0/2).lwage L(0/1).wks occ south smsa ind ms union, ///
  div(occ south smsa ind) dgmmiv(lwage, lagrange(3 4))          ///
  dgmmiv(ms union, lagrange(2 3)) dgmmiv(L.wks, lagrange(1 2))      ///
  lgmmiv(L.lwage wks ms union) twostep vce(robust) artests(3)

******* 9.5: MIXED LINEAR MODELS

* Random intercept model estimated using xtmixed
use mus08psidextract.dta, clear
xtmixed lwage exp exp2 wks ed || id:, mle

* Cluster robust standard errors after xtmixed using bootstrap
xtset id
bootstrap, reps(400) seed(10101) cluster(id) nodots: ///
  xtmixed lwage exp exp2 wks ed || id:, mle

* Random-slopes model estimated using xtmixed
xtmixed lwage exp exp2 wks ed || id: exp wks, covar(unstructured) mle

* Random-coefficients model estimated using xtrc
quietly set matsize 600
xtrc lwage exp wks, i(id)

* List the estimated variance matrix
matrix list e(Sigma)

* Two-way random-effects model estimated using xtmixed
xtmixed lwage exp exp2 wks ed || _all: R.t || id: , mle

// Following not included in book
* gllamm
* gllamm lwage exp exp2 wks ed, i(id) nip(10) adapt
* xtmixed lwage exp exp2 wks ed || id:, mle

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
quietly xtreg pharvis, mle
display "Intra-class correlation for household: " e(rho)
quietly correlate pharvis L1.pharvis
display "Correlation for adjoining household:   " r(rho)

* OLS, RE, and FE estimation with clustering on household and on village
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

********** CLOSE OUTPUT
