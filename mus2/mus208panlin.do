* mus208panlin.do  for Stata 17

cap log close

********** OVERVIEW OF mus208panlin.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* CHAPTER 8
* 8.2: PANEL-DATA METHODS OVERVIEW
* 8.3: SUMMARY OF PANEL DATA
* 8.4: POOLED OR PA ESTIMATORS
* 8.5: FE OR WITHIN ESTIMATOR
* 8.6: BETWEEN ESTIMATOR
* 8.7: RE ESTIMATOR
* 8.8: COMPARISON OF ESTIMATORS 
* 8.9: FIRST-DIFFERENCE ESTIMATOR
* 8.10: PANEL-DATA MANAGEMENT 
 
* To run you need files
*   mus208psid.dta
*   mus208cigarwide.dta
* in your directory

* No community-contributed command are used

********** SETUP **********

clear all
set linesize 80
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus208psid.dta
* PSID. Same as Stata website file psidextract.dta
* Data due to  Baltagi and Khanti-Akom (1990) 
* This is corrected version of data in Cornwell and Rupert (1988).
* 595 individuals for years 1976-82

* mus208cigarwide.dta is a smaller wide form version of mus211cigar.dta
* Due to Baltagi and James M. Griffin (2001), "The Econometrics of Rational Addiction: 
* The Case of Cigarettes", J. of Bus. and  Econ.  Statistics, 449-454. 

******* CHAPTER 8.3: SUMMARY OF PANEL DATA

* Read in dataset and describe
qui use mus208psid
describe

* Summary of dataset
summarize

* Organization of dataset
list id t exp wks occ in 1/3, clean

* Declare individual identifier and time identifier
xtset id t

* Panel description of dataset
xtdescribe 

* Panel summary statistics: Within and between variation
xtsum id t lwage ed exp exp2 wks south tdum1

* id t xtsum lwage occ south smsa ind exp wks ms union fem blk ed tdum1

* Panel tabulation for a variable
xttab south

* Transition probabilities for a variable
xttrans south, freq

* Simple time-series plot for each of 20 individuals
qui xtline lwage if id<=20, overlay 

* Simple time-series plot for each of 20 individuals
qui xtline lwage if id<=20, overlay legend(off) saving(graph1.gph, replace)
qui xtline wks if id<=20, overlay legend(off) saving(graph2.gph, replace)
graph combine graph1.gph graph2.gph, iscale(1.4) ysize(2.5) xsize(6.0)

* Scatterplot, quadratic fit, and nonparametric regression (lowess)
graph twoway (scatter lwage exp, msize(small) msymbol(o))              ///
    (qfit lwage exp, clstyle(p3) lwidth(medthick))                     ///
    (lowess lwage exp, bwidth(0.4) clstyle(p1) lwidth(medthick)),      ///
    plotregion(style(none)) scale(1.2)                                 ///
    title("Overall variation: log wage versus experience")             ///
    xtitle("Years of experience", size(medlarge)) xscale(titlegap(*5)) /// 
    ytitle("log hourly wage", size(medlarge)) yscale(titlegap(*5))     ///
    legend(pos(4) ring(0) col(1)) legend(size(small))                  ///
    legend(label(1 "Actual data") label(2 "Quadratic fit") label(3 "Lowess"))

* Scatterplot for within variation
preserve
xtdata, fe
graph twoway (scatter lwage exp, msize(small) msymbol(o))                ///
    (qfit lwage exp, clstyle(p3) lwidth(medthick))                       ///
    (lowess lwage exp, bwidth(0.4) clstyle(p1) lwidth(medthick)),        ///
    plotregion(style(none)) scale(1.2)                                   ///
    title("Within variation: log wage versus experience")                ///
    xtitle("Years of experience", size(medlarge)) xscale(titlegap(*5))   /// 
    ytitle("log hourly wage", size(medlarge)) yscale(titlegap(*5))       ///
    legend(pos(11) ring(0) col(1)) legend(size(small))                   ///
    legend(label(1 "Actual data") label(2 "Quadratic fit") label(3 "Lowess"))
restore

* graph twoway (scatter lwage exp) (qfit lwage exp) (lowess lwage exp),  ///
* plotregion(style(none)) title("Within variation: log wage versus experience")

* Pooled OLS with cluster–robust standard errors
use mus208psid, clear
regress lwage exp exp2 wks ed, vce(cluster id)

* Pooled OLS with incorrect default standard errors
regress lwage exp exp2 wks ed 

* First-order autocorrelation in a variable
sort id t  
correlate lwage L.lwage

* Autocorrelations of residual 
qui regress lwage exp exp2 wks ed, vce(cluster id)
predict uhat, residuals
forvalues j = 1/6 {
    qui corr uhat L`j'.uhat
    display "Autocorrelation at lag `j' = " %6.3f r(rho) 
    }

* First-order autocorrelation differs in different year pairs
forvalues s = 2/7 {
    qui corr uhat L1.uhat if t == `s'
    display "Autocorrelation at lag 1 in year `s' = " %6.3f r(rho) 
    }

******* CHAPTER 8.4: POOLED ESTIMATORS

* Population-averaged or pooled FGLS estimator with AR(2) error
xtreg lwage exp exp2 wks ed, pa corr(ar 2) vce(robust) nolog

* Estimated error correlation matrix after xtreg, pa
matrix list e(R)

******* CHAPTER 8.5: WITHIN ESTIMATOR

* Within or FE estimator with cluster–robust standard errors
xtreg lwage exp exp2 wks ed, fe vce(cluster id)

* LSDV model fit using areg with cluster–robust standard errors
areg lwage exp exp2 wks ed, absorb(id) vce(cluster id)

* LSDV model fit using factor variables with cluster–robust standard errors
qui regress lwage exp exp2 wks ed i.id, vce(cluster id)
estimates table, keep(exp exp2 wks ed _cons) b se b(%12.7f)

// Not included but mentioned in text
* Two-way fixed effects
xtreg lwage i.t exp exp2 wks ed, fe vce(cluster id)
reg2hdfe lwage exp exp2 wks ed, id1(id) id2(t) cluster(id)
felsdvreg lwage exp exp2 wks ed, ivar(id) jvar(t) cluster(id) peff(peff) ///
    feff(feff) xb(xn) res(res) mover(mover) group(group) mnum(mnum) pobs(pobs)
a2reg lwage exp exp2 wks ed, individual(id) unit(t) 

******* CHAPTER 8.6: BETWEEN ESTIMATOR

* Between estimator with default standard errors
xtreg lwage exp exp2 wks ed, be

* xtreg lwage exp exp2 wks ed, be vce(bootstrap, reps(400) seed(10101) nodots)

******* CHAPTER 8.7: RANDOM EFFECTS ESTIMATORS

* RE estimator with cluster–robust standard errors
xtreg lwage exp exp2 wks ed, re vce(cluster id) theta

* Calculate theta
* display "theta = "  1 - sqrt(e(sigma_e)^2 / (7*e(sigma_u)^2+e(sigma_e)^2))

* Mundlak correction: RE with individual-specific means added as regressors
sort id
foreach x of varlist exp exp2 wks ed {
    by id: egen mean`x' = mean(`x')
    }
xtreg lwage exp exp2 wks ed meanexp meanexp2 meanwks meaned, re vce(cluster id)

******* CHAPTER 8.8: COMPARISON OF ESTIMATORS


* Compare OLS, BE, FE, RE estimators, and methods to compute standard errors
global xlist exp exp2 wks ed 
qui regress lwage $xlist, vce(cluster id)
estimates store OLS_rob
qui xtreg lwage $xlist, be
estimates store BE
qui xtreg lwage $xlist, fe 
estimates store FE
qui xtreg lwage $xlist, fe vce(robust)
estimates store FE_rob
qui xtreg lwage $xlist, re
estimates store RE
qui xtreg lwage $xlist, re vce(robust)
estimates store RE_rob
estimates table OLS_rob BE FE FE_rob RE RE_rob,  ///
    b se stats(N r2 r2_o r2_b r2_w sigma_u sigma_e rho) b(%7.4f)

* Hausman test assuming RE estimator is fully efficient under null hypothesis
hausman FE RE, sigmamore

* Cluster-robust Hausman test using method of Wooldridge (2010)
qui regress lwage $xlist meanexp meanexp2 meanwks, vce(cluster id)
test meanexp meanexp2 meanwks

* Cluster-robust Hausman test using xtovierid command
qui xtreg lwage $xlist, re vce(robust)
xtoverid

* Prediction after OLS and RE estimation
qui regress lwage exp exp2 wks ed, vce(cluster id)
predict xbols, xb
qui xtreg lwage exp exp2 wks ed, re  
predict xbre, xb
predict xbure, xbu
summarize lwage xbols xbre xbure
correlate lwage xbols xbre xbure

*********** CHAPTER 8.9: FIRST DIFFERENCE ESTIMATOR

sort id t
* FD estimator with cluster–robust standard errors
regress D.(lwage exp exp2 wks ed), vce(cluster id)

// Not included - uses community-contributed xtivreg2
// xtivreg2 lwage exp exp2 wks ed, fd cluster(id) small

******* CHAPTER 8.10: WIDE FORM AND LONG FORM DATA

* Wide-form data (observation is a state) 
qui use mus208cigarwide, clear
list, clean

* Convert from wide form to long form (observation is a state-year pair) 
reshape long lnp lnc, i(state) j(year)

* Long-form data (observation is a state) 
list in 1/6, sepby(state)

* Reconvert from long form to wide form (observation is a state) 
reshape wide lnp lnc, i(state) j(year)

list, clean

* Create alternative wide-form data (observation is a year)
qui reshape long lnp lnc, i(state) j(year) 
reshape wide lnp lnc, i(year) j(state)
list year lnp1 lnp2 lnc1 lnc2, clean

* Convert from wide form (observation is year) to long form (year-state)
reshape long lnp lnc, i(year) j(state)
list in 1/6, clean

*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph

*********** END
