* mus209panlinext.do  for Stata 17

cap log close

********** OVERVIEW OF mus209panlinext.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press  

* Chapter 9
*  9.2: PANEL IV ESTIMATORS
*  9.3: HAUSMAN-TAYLOR ESTIMATOR
*  9.4: ARELLANO-BOND ESTIMATOR 
*  9.5: LONG PANELS 

* To run you need files
*   mus208psid.dta
*   mus209cigar.dta
* in your directory

* community-contributed command
*   xtscc
* is used

********** SETUP **********

clear all
set linesize 80
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

*  mus210psid.dta
* PSID. Same as Stata website file psidextract.dta
* Data due to  Baltagi and Khanti-Akom (1990) 
* This is corrected version of data in Cornwell and Rupert (1988).
* 595 individuals for years 1976-82

* mus209cigar.dta
* Due to Baltagi et al. (2001)
* Panel on 46 U.S. states over 30 years 1963-92

******* 9.2: PANEL IV ESTIMATOR

* Panel IV example: FE with wks instrumented by external instrument ms
qui use mus208psid
xtivreg lwage exp exp2 (wks = ms), fe vce(robust)

* xtivreg lwage exp exp2 (wks = ms), fe vce(bootstrap, reps(400) seed(10101))

******* 9.3: HAUSMAN-TAYLOR ESTIMATOR

* Hausman-Taylor example of Baltagi and Khanti-Akom (1990)
xthtaylor lwage occ south smsa ind exp exp2 wks ms union fem blk ed,  ///
  endog(exp exp2 wks ms union ed) vce(robust)

******* 9.4: ARELLANO-BOND ESTIMATOR

* 2SLS or one-step GMM for a pure time-series AR(2) panel model
qui use mus208psid, clear 
xtabond lwage, lags(2) vce(robust)
 
* Optimal or two-step GMM for a pure time-series AR(2) panel model
xtabond lwage, lags(2) twostep vce(robust)

* Reduce the number of instruments for a pure time-series AR(2) panel model
xtabond lwage, lags(2) vce(robust) maxldep(1)

* Optimal or two-step GMM for a dynamic panel model
xtabond lwage occ south smsa ind, lags(2) maxldep(3)       ///
    pre(wks,lag(1,2)) endogenous(ms,lag(0,2))              ///
    endogenous(union,lag(0,2)) twostep vce(robust) artests(3)

* Test whether error is serially correlated
estat abond

* Test of overidentifying restrictions (first estimate with no vce(robust))
qui xtabond lwage occ south smsa ind, lags(2) maxldep(3) ///
    pre(wks,lag(1,2)) endogenous(ms,lag(0,2))            ///
    endogenous(union,lag(0,2)) twostep artests(3)
estat sargan

* xtreg lwage occ south smsa ind wks L1.wks ms union, re vce(robust)

* Arellano/Bover or Blundell/Bond for a dynamic panel model
xtdpdsys lwage occ south smsa ind, lags(2) maxldep(3)    ///
    pre(wks,lag(1,2)) endogenous(ms,lag(0,2))            ///
    endogenous(union,lag(0,2)) twostep vce(robust) artests(3)

estat abond

* Use of xtdpd to exactly reproduce the previous xtdpdsys command
xtdpd L(0/2).lwage L(0/1).wks occ south smsa ind ms union,       ///
    div(occ south smsa ind) dgmmiv(lwage, lagrange(2 4))         ///
    dgmmiv(ms union, lagrange(2 3)) dgmmiv(L.wks, lagrange(1 2)) ///
    lgmmiv(lwage wks ms union) twostep vce(robust) artests(3)

// Command given in text but output not provided 
* Previous command if model error is MA(1)
xtdpd L(0/2).lwage L(0/1).wks occ south smsa ind ms union,       ///
    div(occ south smsa ind) dgmmiv(lwage, lagrange(3 4))         ///
    dgmmiv(ms union, lagrange(2 3)) dgmmiv(L.wks, lagrange(1 2)) ///
    lgmmiv(L.lwage wks ms union) twostep vce(robust) artests(3)

* xtabond2: Two-step GMM for a pure time-series AR(2) panel model
xtabond2 lwage L.lwage L2.lwage, gmmstyle(L2.lwage) h(1) small twostep robust

* xtabond2: Two-step GMM for dynamic model with exogenous & endogenous vars.
xtabond2 lwage  L.lwage L2.lwage occ south smsa ind wks L.wks ms union ///
     tdum3-tdum6, gmmstyle(L.(L.lwage wks ms union))                   ///
     ivstyle(t, equation(level)) twostep small          
 
* xtabond2: Two-step GMM for dynamic model with limit on the lags generating IV
xtabond2 lwage  L.lwage L2.lwage occ south smsa ind wks L.wks ms union ///
   tdum3-tdum6, gmmstyle(L.(L.lwage wks ms union), laglimits(2 3)) twostep small             

******* CHAPTER 9.5: LONG PANEL

* Description of cigarette dataset
qui use mus209cigar, clear
describe

* Summary of cigarette dataset
summarize, separator(6)

* Pooled GLS with error correlated across states and state-specific AR(1) 
xtset state year 
xtgls lnc lnp lny lnpmin year, panels(correlated) corr(psar1)

* Comparison of various pooled OLS and GLS estimators
qui xtpcse lnc lnp lny lnpmin year, corr(ind) independent nmk
estimates store OLS_iid
qui xtpcse lnc lnp lny lnpmin year, corr(ind)
estimates store OLS_cor
qui xtscc lnc lnp lny lnpmin year, lag(4)
estimates store OLS_DK
qui xtpcse lnc lnp lny lnpmin year, corr(ar1)
estimates store AR1_cor
qui xtgls lnc lnp lny lnpmin year, corr(ar1) panels(iid)
estimates store FGLSAR1
qui xtgls lnc lnp lny lnpmin year, corr(ar1) panels(correlated)
estimates store FGLSCAR
estimates table OLS_iid OLS_cor OLS_DK AR1_cor FGLSAR1 FGLSCAR, b(%7.3f) se

* Comparison of various RE and FE estimators
qui use mus209cigar, clear
qui xtscc lnc lnp lny lnpmin, lag(4)
estimates store OLS_DK
qui xtreg lnc lnp lny lnpmin, fe
estimates store FE_REG
qui xtreg lnc lnp lny lnpmin, re
estimates store RE_REG
qui xtregar lnc lnp lny lnpmin, fe
estimates store FE_REGAR
qui xtregar lnc lnp lny lnpmin, re
estimates store RE_REGAR
qui xtscc lnc lnp lny lnpmin, fe lag(4)
estimates store FE_DK
estimates table OLS_DK FE_REG RE_REG FE_REGAR RE_REGAR FE_DK, b(%7.3f) se

* Run separate regressions for each state
statsby, by(state) clear: regress lnc lnp lny lnpmin year

* Report regression coefficients for each state
format _b* %9.2f
list, clean

********** END
