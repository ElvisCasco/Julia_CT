* ct_panel.do based on mus08p1panlin.do mus09p1panlin2.do May 2013 for Stata version 12.0

********** OVERVIEW OF ct_panel.do **********

log using ct_panel.txt, text replace

* STATA Program 
* For A. Colin Cameron and Pravin Trivedi "Lectures in Microeconometrics"
* Binary models

* To run you need files
*   mus08psidextract.dta    
* in your directory
* Stata add-on
*   xtscc
* is used

********** SETUP **********

set more off
version 12.0
set mem 10m
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Source: The Panel Study of Income Dynamics
* Same as Stata website file psidextract.dta
*  Data due to  Baltagi and Khanti-Akom (1990) 
* This is corrected version of data in Cornwell and Rupert (1988).
* 595 individuals for years 1976-82

******* 8.3: PANEL DATA SUMMARY

* Read in data set
use mus08psidextract.dta, clear
drop tdum*

* Describe dataset
describe

* Summarize dataset 
summarize

* Organization of data set
list id t exp wks occ in 1/3, clean

* Declare individual identifier and time identifier
xtset id t

* Panel description of data set
xtdescribe

* Panel summary statistics: within and between variation
xtsum lwage exp ed t 

* Panel tabulation for a variable
xttab south

* Transition probabilities for a variable
xttrans south, freq

* Simple time series plot for each of 10 individuals
* xtline lwage if id<=10, overlay

* First-order autocorrelation in a variable
sort id t  
correlate lwage L.lwage L2.lwage L3.lwage L4.lwage L5.lwage L6.lwage

* Pooled OLS with incorrect default standard errors
regress lwage exp exp2 wks ed, noheader

* Pooled OLS with cluster-robust standard errors
regress lwage exp exp2 wks ed, noheader vce(cluster id)

* Pooled feasible GLS estimator with AR(2) error
xtgee lwage exp exp2 wks ed, corr(ar 2) vce(robust)

* Within or FE estimator with cluster-robust standard errors
xtreg lwage exp exp2 wks ed, fe vce(robust)

* Between estimator with default standard errors
xtreg lwage exp exp2 wks ed, be

* Random effects estimator with cluster-robust standard errors
xtreg lwage exp exp2 wks ed, re vce(robust) theta

* First difference estimator with cluster-robust standard errors
global xlist exp exp2 wks ed
regress D.(lwage $xlist), vce(cluster id)

* Compare various estimators all with cluster-robust se's except BE 
global xlist exp exp2 wks ed 
quietly regress lwage $xlist, vce(cluster id)
estimates store OLS
quietly xtgee lwage exp exp2 wks ed, corr(ar 2) vce(robust)
estimates store PFGLS
quietly xtreg lwage $xlist, be
estimates store BE
quietly xtreg lwage $xlist, re vce(robust)
estimates store RE
quietly xtreg lwage $xlist, fe vce(robust)
estimates store FE
estimates table OLS PFGLS BE RE FE,  b(%9.4f) se stats(N)

/*

* Panel IV example: FE with wks instrumented by external instrument ms
xtivreg lwage exp exp2 (wks = ms), fe 

* use mus08psidextract.dta, clear
* Hausman-taylor example
xthtaylor lwage occ south smsa ind exp exp2 wks ms union fem blk ed,  ///
  endog(exp exp2 wks ms union ed)

* Optimal or two-step GMM for a dynamic panel model
xtabond lwage occ south smsa ind, lags(2) maxldep(3)     ///
  pre(wks,lag(1,2)) endogenous(ms,lag(0,2))              ///
  endogenous(union,lag(0,2)) twostep vce(robust) artests(3)

* Test whether error is serially correlated
estat abond

* Test of overidentifying restrictions (first estimate with no vce(robust)
quietly xtabond lwage occ south smsa ind, lags(2) maxldep(3)     ///
  pre(wks,lag(1,2)) endogenous(ms,lag(0,2))              ///
  endogenous(union,lag(0,2)) twostep artests(3)
estat sargan

* Arellano/Bover or Blundell/Bond for a dynamic panel model
xtdpdsys lwage occ south smsa ind, lags(2) maxldep(3)     ///
  pre(wks,lag(1,2)) endogenous(ms,lag(0,2))              ///
  endogenous(union,lag(0,2)) twostep vce(robust) artests(3)
estat abond

* Use of xtdpd to exactly reproduce the previous xtdpdsys command
xtdpd L(0/2).lwage L(0/1).wks occ south smsa ind ms union, ///
  div(occ south smsa ind) dgmmiv(lwage, lag(2 4))          ///
  dgmmiv(ms union, lag (2 3)) dgmmiv(L.wks, lag(1 2))      ///
  lgmmiv(lwage wks ms union) twostep vce(robust) artests(3)

* Previous command if model error is MA(1)
xtdpd L(0/2).lwage L(0/1).wks occ south smsa ind ms union, ///
  div(occ south smsa ind) dgmmiv(lwage, lag(3 4))          ///
  dgmmiv(ms union, lag (2 3)) dgmmiv(L.wks, lag(1 2))      ///
  lgmmiv(L.lwage wks ms union) twostep vce(robust) artests(3)

* Random slopes model estimated using xtmixed
* xtmixed lwage exp exp2 wks ed || id: exp wks, covar(unstructured) mle

********** CLOSE OUTPUT
* log close
* clear
* exit

