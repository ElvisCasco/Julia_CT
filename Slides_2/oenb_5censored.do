* oenb_5censored.do  based on mus16p1tobit.do  Sept 2010 for Stata version 10.0

log using oenb_5censored.txt, text replace

********** OVERVIEW OF oenb_5censored.do **********

* STATA Program 
* For A. Colin Cameron and Pravin Trivedi "Lectures in Microeconometrics"
* Tobit models

* To run you need files 
*   mus16data.dta 
*   mus03data.dta
* in your directory
* and you need Stata user-written commands ??

********** SETUP **********

set more off
version 10
set mem 10m
set scheme s1mono  /* Graphics scheme */

*** TOBIT EXAMPLE WITH SIMULATED DATA 

* Generate x, y*, ycensored, ytruncated and binary dy
* y* = -2500 + 1000*x + e  where e ~ N(0,1000^2) x ~ N(2.75, 0.6^2)
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
* graph export ct_censoredcondmeans.wmf, replace

*** TOBIT MLE: REAL DATA EXAMPLE 

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

* (Inconsistent) OLS estimates on censored data 
regress ambexp age female educ blhisp totchr ins

* Tobit on censored data
tobit ambexp age female educ blhisp totchr ins, ll(0)

* Marginal effects for censored conditional mean evaluated at x = xbar
mfx compute, predict(ystar(0,.))

* ASIDE: Same AME using Stata 11 command margins
* margins, dydx(*) predict(ystar(0,.)) atmean

*** TOBIT MLE IN LOGS (NOT IN OVERHEADS)

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

* Prediction of censored conditional mean of y based on tobit for lny
predict xb, xb
matrix btobit = e(b)
scalar sigma = btobit[1,e(df_m)+2]   // sigma is estimate of sigma   
generate yfromlny = exp(xb+(sigma^2)/2) * (1-normal((gamma-xb-sigma^2)/sigma))
predict lnycensored, ystar(0,.)
summarize y yfromlny lny lnycensored

*** TWO-PART MODEL

* Two-part model
* First part is probit
probit dy age female educ blhisp totchr ins, nolog

* Second part is lognormal regression for positives
regress lny age female educ blhisp totchr ins if dy==1

* Two-part model prediction
quietly probit dy age female educ blhisp totchr ins
predict dyhat, pr
quietly regress lny age female educ blhisp totchr ins if dy==1
predict xbpos, xb
generate yhatpos = exp(xbpos+0.5*e(rmse)^2)
generate yhat2step = dyhat*yhatpos
summarize yhat2step y

* Heckman MLE without exclusion restrictions
global xlist age female educ blhisp totchr ins
heckman lny $xlist, select(dy = $xlist) nolog

* Heckman 2-step without exclusion restrictions
heckman lny $xlist, select(dy = $xlist) twostep

*** OLS IN MAtA using made-up data

clear
input x y
 1 1
 2 2
 3 2
 4 2
 5 3
end
summarize x y
regress y x, noheader
regress y x, noheader vce(robust)
 
mata
 X = (1,1 \ 1,2 \ 1,3 \ 1,4 \ 1,5)
 y = (1 \ 2 \ 2 \ 2 \ 3)
 X
 y
 b = luinv(X'*X)*X'y
 b
 u = y - X*b
 n = rows(X)
 k = cols(X)
 s2 = u'u/(n-k)
 Vdefault = s2*luinv(X'*X)
 seslope = sqrt(Vdefault[2,2])
 seslope
 XuhatsqX = (u:*X)'(u:*X)*n/(n-k)
 Vwhite = luinv(X'*X)*XuhatsqX*luinv(X'*X)
 seslopewhite = sqrt(Vwhite[2,2])
 seslopewhite
end 

regress y x, noheader vce(robust)

*** OLS IN MATA using actual data

use mus03data.dta, clear
keep if totexp > 0   // Analysis for positive medical expenditures only 
generate cons = 1
local y ltotexp
local xlist suppins phylim actlim totchr age female income cons

mata
  // Create y vector and X matrix from Stata dataset
  st_view(y=., ., "`y'")             // y is nx1
  st_view(X=., ., tokens("`xlist'")) // X is nxk
  XXinv = cholinv(cross(X,X))         // XXinv is inverse of X'X 
  b = XXinv*cross(X,y)               // b = [(X'X)^-1]*X'y
  e = y - X*b
  n = rows(X)
  k = cols(X)
  s2 = (e'e)/(n-k)
  vdef = s2*XXinv               // default VCE not used here
  vwhite = XXinv*((e:*X)'(e:*X)*n/(n-k))*XXinv  // robust VCE
  st_matrix("b",b')             // pass results from Mata to Stata
  st_matrix("V",vwhite)         // pass results from Mata to Stata
end

* Use Stata ereturn display to present nicely formatted results
matrix colnames b = `xlist'
matrix colnames V = `xlist'
matrix rownames V = `xlist'
ereturn post b V
ereturn display

* Compare to Stata regress
regress `y' `xlist', vce(robust)

********** CLOSE OUTPUT
* log close
* clear
* exit

