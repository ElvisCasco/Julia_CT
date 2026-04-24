* MSL_2021.do based on mus218multinomial.do  May 2021 for Stata version 17

capture log close

log using MSL_2021.txt, text replace

* MSL_2021.do

********** OVERVIEW OF MSL_2021.do **********

* Stata program 
* copyright C 2020 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2021)
* Stata Press

* To run you need files
*   mus218hklong.dta 
* in your directory

* No Stata addons are used

********** DATA DESCRIPTION **********

* mus218hk.dta is from J. A. Herriges and C. L. Kling, 
* "Nonlinear Income Effects in Random Utility Models", 
* Review of Economics and Statistics, 81(1999): 62-72
* Also analyzed in Cameron and Trivedi (2005) chapter 15

********** SETUP **********

set more off
* version 17
clear all
set linesize 84
set scheme s1mono  /* Graphics scheme */

********** MULTINOMIAL EXAMPLE: CHOICE OF FISHING MODEL

/*
use mus218hk.dta, clear
summarize, separator(0) 
tabulate mode
* Data are in wide form
list mode price pbeach ppier pprivate pcharter in 1, clean
* Convert data from wide form to long form
generate id = _n
reshape long d p q, i(id) j(fishmode beach pier private charter) string
label variable d "= 1 if this mode is chosen"
label variable p "price of this mode" 
label variable q "catch rate of this mode"
drop price crate
save mus218hklong.dta, replace
*/

* Data used
use mus218hklong.dta, clear
drop if fishmode=="charter" | mode == 4
describe
summarize
list id fishmode d p q income in 1/3, clean noobs

* Conditional logit with alternative-specific and case-specific regressors
* cmset before use of cm commands for alternative specific regressors 
cmset id fishmode
cmclogit d p q, casevars(income) basealternative(beach) nolog vce(robust) 
* Predicted probabilities of choice of each mode and compare to actual freqs
predict pasclogit, pr
* Average marginal effect of change in price
margins, dydx(p)

* Alternative-specific mixed logit or random parameters logit estimation
cmset id fishmode
cmmixlogit d q, casevars(income) random(p) basealternative(pier) ///
    vce(robust) nolog
* Average marginal effects with respect to price
margins, dydx(p)

********** END
