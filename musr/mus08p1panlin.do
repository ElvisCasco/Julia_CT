* mus08p1panlin.do Oct 2009 for Stata version 11

cap log close

********** OVERVIEW OF mus08p1panlin.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Chapter 8
* 8.3: PANEL-DATA SUMMARY
* 8.4: POOLED OR POPULATION-AVERAGED ESTIMATORS
* 8.5: WITHIN ESTIMATOR
* 8.6: BETWEEN ESTIMATOR
* 8.7: RANDOM EFFECTS ESTIMATOR
* 8.8: COMPARISON OF ESTIMATORS 
* 8.9: FIRST DIFFERENCE ESTIMATOR
* 8.10: LONG PANELS
* 8.11: PANEL-DATA MANAGEMENT

* To run you need files
*   mus08psidextract.dta
*   mus08cigar.dta
*   mus08cigarwide.dta
* in your directory
* Stata user-written command
*   xtscc
* is used

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

* mus08cigar.dta
* Due to Baltagi et al. (2001)
* Panel on 46 U.S. states over 30 years 1963-92

* mus08cigarwide.dta is a smaller wide from version of mus08cigar.dta

******* 8.3: PANEL-DATA SUMMARY

* Read in dataset and describe
use mus08psidextract.dta, clear
describe

* Summary of dataset
summarize

* Organization of dataset
list id t exp wks occ in 1/3, clean

* Declare individual identifier and time identifier
xtset id t

* Panel description of dataset
xtdescribe 

* Panel summary statistics: within and between variation
xtsum id t lwage ed exp exp2 wks south tdum1

* Panel tabulation for a variable
xttab south

* Transition probabilities for a variable
xttrans south, freq

// Following simpler command not included in book
* Simple time-series plot for each of 20 individuals
quietly xtline lwage if id<=20, overlay 

* Simple time-series plot for each of 20 individuals
quietly xtline lwage if id<=20, overlay legend(off) saving(lwage, replace)
quietly xtline wks if id<=20, overlay legend(off) saving(wks, replace)
graph combine lwage.gph wks.gph, iscale(1)
quietly graph export mus08timeseriesplot.eps, replace

* Scatterplot, quadratic fit and nonparametric regression (lowess)
graph twoway (scatter lwage exp, msize(small) msymbol(o))              ///
  (qfit lwage exp, clstyle(p3) lwidth(medthick))                       ///
  (lowess lwage exp, bwidth(0.4) clstyle(p1) lwidth(medthick)),        ///
  plotregion(style(none))                                              ///
  title("Overall variation: Log wage versus experience")               ///
  xtitle("Years of experience", size(medlarge)) xscale(titlegap(*5))   /// 
  ytitle("Log hourly wage", size(medlarge)) yscale(titlegap(*5))       ///
  legend(pos(4) ring(0) col(1)) legend(size(small))                    ///
  legend(label(1 "Actual Data") label(2 "Quadratic fit") label(3 "Lowess"))
graph export mus08scatterplot.eps, replace

* Scatterplot for within variation
preserve
xtdata, fe
graph twoway (scatter lwage exp) (qfit lwage exp) (lowess lwage exp),  ///
  plotregion(style(none)) title("Within variation: Log wage versus experience")
restore
graph export mus08withinplot.eps, replace

* Pooled OLS with cluster-robust standard errors
use mus08psidextract.dta, clear
regress lwage exp exp2 wks ed, vce(cluster id)

* Pooled OLS with incorrect default standard errors
regress lwage exp exp2 wks ed 

* First-order autocorrelation in a variable
sort id t  
correlate lwage L.lwage

* Autocorrelations of residual 
quietly regress lwage exp exp2 wks ed, vce(cluster id)
predict uhat, residuals
forvalues j = 1/6 {
     quietly corr uhat L`j'.uhat
     display "Autocorrelation at lag `j' = " %6.3f r(rho) 
     }

* First-order autocorrelation differs in different year pairs
forvalues s = 2/7 {
     quietly corr uhat L1.uhat if t == `s'
     display "Autocorrelation at lag 1 in year `s' = " %6.3f r(rho) 
     }

******* 8.4: POOLED OR POPULATION-AVERAGED ESTIMATORS

* Population-averaged or pooled FGLS estimator with AR(2) error
xtreg lwage exp exp2 wks ed, pa corr(ar 2) vce(robust) nolog

* Estimated error correlation matrix after xtreg, pa
matrix list e(R)

******* 8.5: WITHIN ESTIMATOR

* Within or FE estimator with cluster-robust standard errors
xtreg lwage exp exp2 wks ed, fe vce(cluster id)

* LSDV model fit using areg with cluster-robust standard errors
areg lwage exp exp2 wks ed, absorb(id) vce(cluster id)

* LSDV model fit using factor variables with cluster-robust standard errors
set matsize 800
quietly regress lwage exp exp2 wks ed i.id, vce(cluster id)
estimates table, keep(exp exp2 wks ed _cons) b se b(%12.7f)

******* 8.6: BETWEEN ESTIMATOR

* Between estimator with default standard errors
xtreg lwage exp exp2 wks ed, be

// Following gives heteroskedasrtic-robust se's for between estimator
* xtreg lwage exp exp2 wks ed, be vce(boot, reps(400) seed(10101) nodots)

******* 8.7: RANDOM EFFECTS ESTIMATORS

* Random-effects estimator with cluster-robust standard errors
xtreg lwage exp exp2 wks ed, re vce(cluster id) theta

* Calculate theta
* display "theta = "  1 - sqrt(e(sigma_e)^2 / (7*e(sigma_u)^2+e(sigma_e)^2))

******* 8.8: COMPARISON OF ESTIMATORS

use mus08psidextract.dta, clear

* Compare OLS, BE, FE, RE estimators, and methods to compute standard errors
global xlist exp exp2 wks ed 
quietly regress lwage $xlist, vce(cluster id)
estimates store OLS_rob
quietly xtreg lwage $xlist, be
estimates store BE
quietly xtreg lwage $xlist, fe 
estimates store FE
quietly xtreg lwage $xlist, fe vce(robust)
estimates store FE_rob
quietly xtreg lwage $xlist, re
estimates store RE
quietly xtreg lwage $xlist, re vce(robust)
estimates store RE_rob
estimates table OLS_rob BE FE FE_rob RE RE_rob,  ///
  b se stats(N r2 r2_o r2_b r2_w sigma_u sigma_e rho) b(%7.4f)

* Hausman test assuming RE estimator is fully efficient under null hypothesis
hausman FE RE, sigmamore

* Robust Hausman test using method of Wooldridge (2002)
quietly xtreg lwage $xlist, re
scalar theta = e(theta)
global yandxforhausman lwage exp exp2 wks ed
sort id
foreach x of varlist $yandxforhausman {
  by id: egen mean`x' = mean(`x')
  generate md`x' = `x' - mean`x'
  generate red`x' = `x' - theta*mean`x'
  }
quietly regress redlwage redexp redexp2 redwks reded mdexp mdexp2 mdwks, vce(cluster id)
test mdexp mdexp2 mdwks

* Prediction after OLS and RE estimation
quietly regress lwage exp exp2 wks ed, vce(cluster id)
predict xbols, xb
quietly xtreg lwage exp exp2 wks ed, re  
predict xbre, xb
predict xbure, xbu
summarize lwage xbols xbre xbure
correlate lwage xbols xbre xbure

******* 8.9: FIRST DIFFERENCE ESTIMATOR

sort id t
* First-differences estimator with cluster-robust standard errors
regress D.(lwage exp exp2 wks ed), vce(cluster id) noconstant

******* 8.10: LONG PANEL

* Description of cigarette dataset
use mus08cigar.dta, clear
describe

* Summary of cigarette dataset
summarize, separator(6)

* Pooled GLS with error correlated across states and state-specific AR(1) 
xtset state year 
xtgls lnc lnp lny lnpmin year, panels(correlated) corr(psar1)

* Comparison of various pooled OLS and GLS estimators
quietly xtpcse lnc lnp lny lnpmin year, corr(ind) independent nmk
estimates store OLS_iid
quietly xtpcse lnc lnp lny lnpmin year, corr(ind)
estimates store OLS_cor
quietly xtscc lnc lnp lny lnpmin year, lag(4)
estimates store OLS_DK
quietly xtpcse lnc lnp lny lnpmin year, corr(ar1)
estimates store AR1_cor
quietly xtgls lnc lnp lny lnpmin year, corr(ar1) panels(iid)
estimates store FGLSAR1
quietly xtgls lnc lnp lny lnpmin year, corr(ar1) panels(correlated)
estimates store FGLSCAR
estimates table OLS_iid OLS_cor OLS_DK AR1_cor FGLSAR1 FGLSCAR, b(%7.3f) se

* Run separate regressions for each state
statsby, by(state) clear: regress lnc lnp lny lnpmin year

* Report regression coefficients for each state
format _b* %9.2f
list, clean

* Comparison of various RE and FE estimators
use mus08cigar.dta, clear
quietly xtscc lnc lnp lny lnpmin, lag(4)
estimates store OLS_DK
quietly xtreg lnc lnp lny lnpmin, fe
estimates store FE_REG
quietly xtreg lnc lnp lny lnpmin, re
estimates store RE_REG
quietly xtregar lnc lnp lny lnpmin, fe
estimates store FE_REGAR
quietly xtregar lnc lnp lny lnpmin, re
estimates store RE_REGAR
quietly xtscc lnc lnp lny lnpmin, fe lag(4)
estimates store FE_DK
estimates table OLS_DK FE_REG RE_REG FE_REGAR RE_REGAR FE_DK, b(%7.3f) se

******* 8.11: PANEL-DATA MANAGEMENT

* Wide form data (observation is a state) 
use mus08cigarwide.dta, clear
list, clean

* Convert from wide form to long form (observation is a state-year pair) 
reshape long lnp lnc, i(state) j(year)

* Long-form data (observation is a state) 
list in 1/6, sepby(state)

* Reconvert from long form to wide form (observation is a state) 
reshape wide lnp lnc, i(state) j(year)

list, clean

* Create alternative wide-form data (observation is a year)
quietly reshape long lnp lnc, i(state) j(year) 
reshape wide lnp lnc, i(year) j(state)
list year lnp1 lnp2 lnc1 lnc2, clean

* Convert from wide form (observation is year) to long form (year-state)
reshape long lnp lnc, i(year) j(state)
list in 1/6, clean

********** CLOSE OUTPUT
