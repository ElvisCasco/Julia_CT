* ct_ldv.do  based on MUS files May 2013 for Stata version 12.0

* Extracts of mus14p1bin.do       binary
*             mus15p1censored.do  multinomial
*             mus16censored.do    censored

log using ct_ldv.txt, text replace

********** OVERVIEW OF ct_ldv.do **********

* STATA Program 
* For A. Colin Cameron  "Lectures in Microeconometrics"
* Brief summary: Binary models, multinomial models, censored models

* BINARY CHOICE: GENERATED DATA EXAMPLE
* BINARY CHOICE: EXAMPLE: PRIVATE HEALTH INSURANCE
* MULTINOMIAL CHOICE: EXAMPLE: FISHING MODE CHOICE
* CENSORED DATA: SIMULATED DATA TOBIT EXAMPLE
* CENSORED DATA: EXAMPLE: AMBULATORY EXPENDITURE

* To run you need files 
*   mus14data.dta 
*   mus15data.dta 
*   mus16data.dta 
* in your directory

********** SETUP **********

set more off
version 12
set mem 10m
set scheme s1mono  /* Graphics scheme */
 
********** BINARY CHOICE: GENERATED DATA EXAMPLE

* Generated data example with one regressor
set seed 10101
quietly set obs 200
quietly generate x = rnormal(0,1)
quietly generate y = 1 + 1*x + rnormal(0,1) > 0
summarize y x
quietly probit y x
quietly predict probit
quietly logit y x
quietly predict plogit
quietly regress y x
quietly predict pols
quietly sort x
graph twoway (scatter y x, msize(small) jitter(3))                  ///
  (line plogit x, clstyle(p1)) (line probit x, clstyle(p2))         ///
  (line pols x, clstyle(p3)), scale (1.2) plotregion(style(none))   ///
  xti("Regressor x", size(medlarge)) xsca(titlegap(*5))             /// 
  yti("Predicted Pr[y=1|x]", size(medlarge)) yscale(titlegap(*5))   ///
  legend(pos(4) ring(0) col(1)) legend(size(small))                 ///
  legend( label(1 "Data (jittered)") label(2 "Logit")               ///
          label(3 "Probit") label(4 "OLS"))
graph export ct_binarygraph.eps, replace

********** BINARY CHOICE: EXAMPLE: PRIVATE HEALTH INSURANCE

* Data Set comes from the 2000 Health and Retirement Survey
* Medicare benificiaries (mostly elderly)

* Describe and summarize dependent variable and regressors
use mus14data.dta, clear
label variable ins "1 if have private health insurance"
label variable retire "1 if retired"
label variable age "age in years"
label variable hstatusg "1 if health status good of better"
label variable hhincome "household annual income in $000's"
label variable educyear "years of education"
label variable married "1 if married"
label variable hisp "1 if hispanic"

describe ins retire age hstatusg hhincome educyear married hisp
summarize ins retire age hstatusg hhincome educyear married hisp

bysort ins: summarize retire age hstatusg hhincome educyear married hisp

* Logit regression
logit ins retire age hstatusg hhincome educyear married hisp

* Average marginal effect
margins, dydx(*)

* Logit, probit and OLS
logit ins retire age hstatusg hhincome educyear married hisp
probit ins retire age hstatusg hhincome educyear married hisp
regress ins retire age hstatusg hhincome educyear married hisp, vce(robust)

* Comparison of coefficient estimates from logit, probit and OLS 
* Estimation of several models
quietly logit ins retire age hstatusg hhincome educyear married hisp
estimates store blogit
quietly probit ins retire age hstatusg hhincome educyear married hisp
estimates store bprobit
quietly regress ins retire age hstatusg hhincome educyear married hisp
estimates store bols
quietly logit ins retire age hstatusg hhincome educyear married hisp, vce(robust)
estimates store blogitr
quietly probit ins retire age hstatusg hhincome educyear married hisp, vce(robust)
estimates store bprobitr
quietly regress ins retire age hstatusg hhincome educyear married hisp, vce(robust)
estimates store bolsr
* Compare coefficient estimates across models with default and robust standard errors
estimates table blogit bprobit bols blogitr bprobitr bolsr, /// 
  stats(N ll) b(%7.3f) t(%7.2f) stfmt(%8.2f)

* Comparison of predicted probabilities from logit, probit and OLS
quietly logit ins retire age hstatusg hhincome educyear married hisp
predict plogit, p
quietly probit ins retire age hstatusg hhincome educyear married hisp
predict pprobit, p
quietly regress ins retire age hstatusg hhincome educyear married hisp
quietly predict pOLS
summarize ins plogit pprobit pOLS

* Marginal effects for logit: AME differs from MEM
quietly logit ins retire age hstatusg hhincome educyear married hisp
margins, dydx(*) 
margins, dydx(*) atmean

* Old Stata comamnds superceded by margins
* margeff
* mfx

********** MULTINOMIAL CHOICE: EXAMPLE: FISHING MODE CHOICE

* Data Set comes from :
* J. A. Herriges and C. L. Kling, 
* "Nonlinear Income Effects in Random Utility Models", 
* Review of Economics and Statistics, 81(1999): 62-72

* Read in data and summarize
use mus15data.dta, clear
describe
list mode d* p* income in 1/2, clean
summarize d* p* q* income, separator(4)
preserve
bysort mode: summarize
restore

* Multinomial logit with base outcome alternative 1
mlogit mode income, baseoutcome(1)

* Compare average predicted probabilities to sample average frequencies
predict pmlogit1 pmlogit2 pmlogit3 pmlogit4, pr
summarize pmlogit* dbeach dpier dprivate dcharter, separator(4)

* AME of income change for outcome 3
margins, dydx(*) predict(outcome(3))
* MEM of income change for outcome 3
margins, dydx(*) predict(outcome(3)) atmean


********** CENSORED DATA: SIMULATED DATA TOBIT EXAMPLE

* Generate x, y*, ycensored, ytruncated and binary dy
* y* = -2500 + 1000*x + e  where e ~ N(0,1000^2) x ~ N(2.75, 0.6^2)
clear
set seed 10101
set obs 200 
generate e = rnormal(0,1000)
generate x = rnormal(2.75,0.6)
generate w = exp(x)   // not necessary
generate ystar = -2500 + 1000*x + e
generate ytruncated = ystar
replace ytruncated = . if (ystar < 0)
generate ycensored = ystar
replace ycensored = 0 if (ystar < 0)
generate dy = ycensored
replace dy = 1 if (ycensored>0)
summarize

* Compare various OLS estimates
quietly regress ystar x      // In practice not possible as y* not observed
estimates store COMPLETE  
quietly regress ycensored x  // Inconsistent  
estimates store CENSORED
quietly regress ytruncated x if dy==1  // Inconsistent 
estimates store TRUNCATED
estimates table COMPLETE CENSORED TRUNCATED, keep(_cons x) stats(N) b(%10.0f) se 

* Compute the theoretical (using d.g.p. betas) censored and truncated means
generate xb = -2500 + 1000*x
generate sigma = 1000
generate PHIxb = normal(xb/sigma)
generate phixb = normalden(xb/sigma)
generate lamda = phixb/PHIxb
generate Eytruncated = xb + sigma*lamda
generate Eycensored = PHIxb*Eytruncated
summarize

* Plot with scatterplot of data plus conditional means 
sort x
graph twoway (scatter ystar x, msize(small)) /// 
  (scatter Eytruncated x, c(l) msize(vtiny) clstyle(p3) clwidth(medthick)) ///
  (scatter Eycensored x, c(l) msize(vtiny) clstyle(p2) clwidth(medthick))  ///
  (scatter xb x, c(l) msize(vtiny) clstyle(p1) clwidth(medthick)),         ///
  scale (1.2) plotregion(style(none))                                      ///
  title("Tobit: Censored and Truncated Means") ///
  xtitle("x (natural logarithm of wage)", size(medlarge)) xscale(titlegap(*5)) /// 
  ytitle("Different Conditional Means", size(medlarge)) yscale(titlegap(*5)) ///
  legend(pos(5) ring(0) col(1)) legend(size(small)) ///
  legend( label(1 "Actual Latent Variable") label(2 "Truncated Mean") ///
          label(3 "Censored Mean") label(4 "Uncensored Mean"))
graph export ct_censoredcondmeans.wmf, replace

********** CENSORED DATA: EXAMPLE: AMBULATORY EXPENDITURE

/*  Subset of data in P. Deb, M. Munkin and P.K. Trivedi (2006)
    "Bayesian Analysis of Two-Part Model with Endogeneity", Journal of Applied Econometrics, 21, 1081-1100
    Only the data for year 2001 are used
    ambexp   Ambulatory medical expenditures (excluding dental and outpatient mental) 
    lambexp  Ln(ambexp) given ambexp > 0 ; missing otherwise
    dambexp  1 if ambexp > 0 and 0 otherwise
    lnambexp  ln(ambexp) if ambexp>0 and 0 if ambexp=0  
    age     age in years/10 
    female   1 for females, zero otherwise
    educ     years of schooling of decision maker
    blhisp   either black or hispanic 
    totchr     number of chronic diseases
    ins either PPO or HMO type insurance  */

* Raw data summary
use mus16data.dta, clear
summarize ambexp dambexp age female educ blhisp totchr ins

* Tobit on censored data
tobit ambexp age female educ blhisp totchr ins, ll(0)

* Set censoring point for data in logs (see MUS p.532 for explanation)
use mus16data.dta, clear
generate y = ambexp
generate dy = ambexp > 0
quietly generate lny = ln(y)                // Zero values will become missing
quietly summarize lny            
scalar gamma = r(min)               // This could be negative
quietly replace lny = gamma - 0.0000001 if lny == .

* Now do tobit on lny
tobit lny age female educ blhisp totchr ins, ll

* Heckman 2-step without exclusion restrictions
global xlist age female educ blhisp totchr ins
heckman lny $xlist, select(dy = $xlist) twostep

********** CLOSE OUTPUT
* log close
* clear
* exit

