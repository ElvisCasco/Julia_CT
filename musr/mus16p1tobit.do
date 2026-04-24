* mus16p1tobit.do  Oct 2009 for Stata version 11 

cap log close

********** OVERVIEW OF mus16p1tobit.DO **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press 

* Chapter 16
* 16.3 TOBIT MODEL EXAMPLE
* 16.4 TOBIT FOR LOGNORMAL DATA
* 16.5 TWO-PART MODEL IN LOGS
* 16.6 SELECTION MODEL
* 16.7 PREDICTION FROM MODEL WITH OUTCOME IN LOGS

* Tobit and selection models example
* It provides
*   (1) Tobit estimation 
*   (2) Tobit analysis after log transformation  
*   (3) Two-part model estimation (similar to Table 16.1 in MMA)
*   (4) Selection model estimation
*     (4A) ML estimates  (Table 16.1)
*     (4B) Heckman 2-step estimates  (Table 16.1)
*     (4C) Check for possible collinearity problems in Heckman 2-Step 
*   (5) Prediction analysis

* To use this program you need file
*   mus16data.dta 

********** SETUP **********

set more off
version 11
clear all
set scheme s1mono   /* Used for graphs */

********** DATA DESCRIPTION **********

* Subset of data in P. Deb, M. Munkin and P.K. Trivedi (2006)
* "Bayesian Analysis of Two-Part Model with Endogeneity", Journal of Applied Econometrics, 21, 1081-1100
* Only the data for year 2001 are used

* Dependent variable is 
*      ambexp   Ambulatory medical expenditures (excluding dental and outpatient mental) 
*      lambexp  Ln(ambexp) given ambexp > 0 ; missing otherwise
*      dambexp  1 if ambexp > 0 and 0 otherwise
*      lnambexp  ln(ambexp) if ambexp>0 and 0 if ambexp=0  

* Regressors are 
*  - Health insurance measures
*       ins either PPO or HMO type insurance
*  - Health status measures
*       totchr     number of chronic diseases
*  - Socioeconomic characteristics
*       age     age in years/10 
*       female   1 for females, zero otherwise
*       educ     years of schooling of decision maker
*       blhisp   either black or hispanic 
*        
 
************ 16.3 TOBIT MODEL EXAMPLE

* Raw data summary
use mus16data.dta, clear
summarize ambexp age female educ blhisp totchr ins 

* Detailed summary to show skewness and kurtosis
summarize ambexp, detail

* Summary for positives only
summarize ambexp if ambexp >0, detail 

* Tobit analysis for ambexp using all expenditures
global xlist age female educ blhisp totchr ins   //Define regressor list $xlist
tobit ambexp $xlist, ll(0) 

* Tobit prediction and summary
predict yhatlin 
summarize yhatlin
 
** Tobit postestimation

* (1) ME on censored expected value E(y|x,y>0)
quietly tobit ambexp age i.female educ i.blhisp totchr i.ins, ll(0)
margins, dydx(*) predict(e(0,.)) atmean noatlegend

* (2) ME without censoring on E(y|x)
margins, dydx(*) predict(ystar(0,.)) atmean noatlegend

* (3) ME when E(y|0<y<535)
margins, dydx(*) predict(e(0,535)) atmean noatlegend

* Direct computation of marginal effects for E(y|x)
predict xb1, xb // xb1 is estimate of x'b
matrix btobit = e(b)
scalar sigma = btobit[1,11] // sigma is estimate of sigma
matrix bcoeff = btobit[1,1..9] // bcoeff is betas excl. constant
quietly summarize xb1 
scalar meanxb = r(mean) // mean of x'b equals (mean of x)'b 
scalar PHI = normal(meanxb/sigma) 
matrix deriv = PHI*bcoeff 
matrix list deriv 
* The following gives nicer looking results
ereturn post deriv 
ereturn display

* Compute margins on Pr(5000<ambexp<10000)  
quietly tobit ambexp age i.female educ i.blhisp totchr i.ins, ll(0) vce(robust)
margins, dydx(*) predict(pr(5000,10000)) atmean noatlegend

************ 16.4 TOBIT FOR LOGNORMAL DATA

* ASIDE: IN BELOW ....
* We observe  y = 0       if dy=0   and  y = y  if dy=1
* Let gamma be the smallest of the positive y's
* We define lny = gamma-0.000001 if dy=0   and lny = ln(y) if dy=1 
* This ensures that tobit, ll will correctly sort out
* the censored and uncensored observations
* Note that can't set lny = 0 if dy=1 as ln(y) < 0 if gamma
* drop xb

* Summary of log(expenditures) for positives only
summarize lambexp, detail

* "Tricking" Stata to handle log transformation
generate y = ambexp
generate dy = ambexp > 0
generate lny = ln(y)                // Zero values will become missing
quietly summarize lny            
scalar gamma = r(min)               // This could be negative
display "gamma = " gamma
replace lny = gamma - 0.0000001 if lny == .
tabulate y if y < 0.02              // .02 is arbitrary small value
tabulate lny if lny < gamma + 0.02 
* Label the variables
label variable y "ambexp"
label variable lny "lnambexp"
label variable dy "dambexp"

// Not included in book  
summarize y lny dy $xlist
summarize y lny dy, detail

* Now do tobit on lny and calculate threshold and lambda
tobit lny $xlist, ll

* OLS, not tobit
regress lny $xlist, noheader

* Now do two-limit tobit 
scalar upper = log(10000)
display upper
tobit lny $xlist, ll ul(9.2103404)

* Mills' ratio
quietly tobit lny $xlist, ll
predict xb, xb                            // xb is estimate of x'b
matrix btobit = e(b)
scalar sigma = btobit[1,e(df_m)+2]        // sigma is estimate of sigma
generate threshold = (gamma-xb)/sigma     // gamma: lower censoring point 
generate lambda = normalden(threshold)/normal(threshold)

* Generalized residuals
* gres1 and gres2 should have mean zero by the first-order conditions
* gres3 and gres4 have mean zero if model correctly specified
* Residual (scaled by sigma) for positive values
quietly generate uifdyeq1 = (lny - xb)/sigma if dy == 1
* First-sample moment 
quietly generate double gres1 = uifdyeq1 
quietly replace gres1 = -lambda if dy == 0
summarize gres1

* Second- to fourth-sample moments
quietly generate double gres2 = uifdyeq1^2 - 1 
quietly replace gres2 = -threshold*lambda if dy == 0
quietly generate double gres3 = uifdyeq1^3 
replace gres3 = -(2 + threshold^2)*lambda if dy == 0
quietly generate double gres4 = uifdyeq1^4 - 3 
quietly replace gres4 = -(3*threshold + threshold^3)*lambda if dy == 0

* Generate the scores to use in the LM test
* These are actual scores up to a multiple constant over i
* So get same LM test as would if added the multiple
* The scores wrt b are gres1 times the relevant component of x

* Generate the scores to use in the LM test
foreach var in $xlist {
   generate score`var' = gres1*`var'
}
global scores score* gres1 gres2

// Following not included in book
* The score wrt sigma can be shown to be gres2
* So the scores are the score`var', gres1 (for the intercept) and gres2
* Test that calculations done correct
* All the score components should have a mean of zero
* gres3 and gres4 may not if model misspecified and these are the basis of the test
summarize $scores gres3 gres4 threshold lambda

* Test of normality in tobit regression
* NR^2 from the uncentered regression has chi-squared distribution
generate one = 1
quietly regress one gres3 gres4 $scores, noconstant
display "N R^2 = " e(N)*e(r2) " with p-value = " chi2tail(2,e(N)*e(r2))

* Test of homoskedasticity in tobit regression
foreach var in $xlist {
   generate score2`var' = gres2*`var'
}
global scores2 score* score2* gres1 gres2
* summarize $scores2
quietly regress one gres3 gres4 $scores2, noconstant
display "N R^2 = " e(N)*e(r2) " with p-value = " chi2tail(2,e(N)*e(r2))

// Following not included in book
* Add-on command tobcm is for left-censored at 0 with 0 a censored observation
* For our example lny is left-censored at 0- but 0 is a noncensored observation
* This leads to slightly different result due to this single observation
tobit lny $xlist, ll(0)
tobcm

************ 16.5 TWO-PART MODEL IN LOGS

* Part 1 of the two-part model
probit dy $xlist, nolog
scalar llprobit = e(ll)

* Part 2 of the two-part model
regress lny $xlist if dy==1
scalar lllognormal = e(ll)
predict rlambexp, residuals

* Create two-part model log likelihood 
scalar lltwopart = llprobit + lllognormal  //two-part model log likelihood
display "lltwopart = " lltwopart

* hettest and sktest commands
quietly regress lny $xlist if dy==1
hettest
sktest rlambexp

************ 16.6 SELECTION MODEL

* Heckman MLE without exclusion restrictions
heckman lny $xlist, select(dy = $xlist) nolog

* Heckman 2-step without exclusion restrictions
heckman lny $xlist, select(dy = $xlist) twostep
 
* Heckman MLE with exclusion restriction
heckman lny $xlist, select(dy = $xlist income) nolog

// Following not included in book
* Heckman 2-step with exclusion restrictions
heckman lny $xlist, select(dy = $xlist income) twostep

************ 16.7 PREDICTION FROM MODELS WITH OUTCOME IN LOGS

// The following yields results that differ from book - book is in error
* Prediction from tobit on lny
generate yhat = exp(xb+0.5*sigma^2)*(1-normal((gamma-xb-sigma^2)/sigma))
generate ytrunchat = yhat / (1 - normal(threshold)) if dy==1
summarize y yhat 
summarize y yhat ytrunchat if dy==1

* Two-part model predictions
quietly probit dy $xlist
predict dyhat, pr
quietly regress lny $xlist if dy==1
predict xbpos, xb
generate  yhatpos = exp(xbpos+0.5*e(rmse)^2)

* Unconditional prediction from two-part model
generate yhat2step = dyhat*yhatpos
summarize yhat2step y
summarize yhatpos y if dy==1

* Heckman model predictions
quietly heckman lny $xlist, select(dy = $xlist) 
predict probpos, psel 
predict x1b1, xbsel 
predict x2b2, xb 
scalar sig2sq = e(sigma)^2 
scalar sig12sq = e(rho)*e(sigma)^2 
display "sigma1sq = 1" " sigma12sq = " sig12sq " sigma2sq = " sig2sq 
generate yhatheck = exp(x2b2 + 0.5*(sig2sq))*(1 - normal(-x1b1-sig12sq)) 
generate yhatposheck = yhatheck/probpos 
summarize yhatheck y probpos dy 
summarize yhatposheck probpos dy y if dy==1

// Following not included in book 
* THE FOLLOWING SHOWS THAT THE POOR PREDICTIONS ARE DUE TO HIGH SIGMA
* THE ABOVE MODEL HAD SIGMA = 2.78
* BY CONTRAST IF WE JUST REGRESS LNY on $XLIST WITH POSITIVE OBS THEN SIGMA = 1.27 
scalar sigma = 1 
generate yhatsigma1 = exp(xb+0.5*sigma^2)*(1-normal((gamma-xb-sigma^2)/sigma))
scalar sigma = 2
generate yhatsigma2 = exp(xb+0.5*sigma^2)*(1-normal((gamma-xb-sigma^2)/sigma))
summarize yhat*
summarize yhat* if dy==1

********** CLOSE OUTPUT
