* mus17p1cnt.do Oct 2009 for Stata version 11

cap log close

********** OVERVIEW OF mus17p1cnt.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Chapter 17
* This program analyzes the data used in chapter 17
* The original data is in file ??

* 17.2: FEATURES OF COUNT DATA
* 17.3: EMPIRICAL EXAMPLE 1
* 17.4: EMPIRICAL EXAMPLE 2
* 17.5: MODELS WITH ENDOGENOUS REGRESSORS

* To run you need file mus17data.dta    
* in your directory
* Stata user-written commands
*   spost9_ado     // For prvalue, prcount, listcoef in section 17.3, 17.4
*   hplogit        // For hurdle model in section 17.3
*   hnblogit       // For hurdle model in section 17.3
*   fmm            // For finite mixture models in section 17.3
* are used

********** SETUP **********

set more off
version 11
clear all
set scheme s1mono  /* Graphics scheme */

********** 17.2: FEATURES OF COUNT DATA

* Poisson (mu=1) generated data
quietly set obs 10000
set seed 10101              // set the seed !
generate xpois= rpoisson(1) // draw from Poisson(mu=1)    
summarize xpois
tabulate xpois

histogram xpois, discrete xtitle("Poisson")
graph save mus17xp, replace

* Negative binomial (mu=1 var=2) generated data
set seed 10101                   // set the seed !
generate xg = rgamma(1,1)
generate xnegbin = rpoisson(xg)  // NB generated as a Poisson-gamma mixture 
summarize xnegbin
tabulate xnegbin

histogram xnegbin, discrete xtitle("Negative binomial")
graph save mus17negbin, replace


********** 17.3: EMPIRICAL EXAMPLE 1

drop _all 

* Summary statistics for doctor visits data
use mus17data.dta              
global xlist private medicaid age age2 educyr actlim totchr 
summarize docvis $xlist

* Tabulate docvis after recoding values > 10 to ranges 11-40 or 41-143
generate dvrange = docvis
recode dvrange (11/40 = 40) (41/143 = 143)
tabulate dvrange

* Poisson with default ML standard errors
poisson docvis $xlist, nolog         

* Squared correlation between y and yhat
drop yphat
predict yphat, n
quietly correlate docvis yphat
display "Squared correlation between y and yhat = " r(rho)^2

* Poisson with robust standard errors
poisson docvis $xlist, vce(robust) nolog  // Poisson robust SEs

* Overdispersion test against V(y|x) = E(y|x) + a*{E(y|x)^2}
quietly poisson docvis $xlist, vce(robust)
predict muhat, n
quietly generate ystar = ((docvis-muhat)^2 - docvis)/muhat
regress ystar muhat, noconstant noheader

* Average marginal effect for Poisson
quietly poisson docvis i.private i.medicaid age age2 educyr i.actlim totchr, vce(robust)
margins, dydx(*)

// Following not included in book
* Poisson as glm
glm docvis $xlist, family(poisson) link(log) robust  // glm

* Standard negative binomial (NB2) with default SEs
nbreg docvis $xlist, nolog   

* Squared correlation between y and yhat
predict ynbhat, n
quietly correlate docvis ynbhat
display "Squared correlation between y and yhat = " r(rho)^2

* Poisson: Sample vs avg predicted probabilities of y = 0, 1, ..., 5 
countfit docvis $xlist, maxcount(5) prm nograph noestimates nofit

* NB2: Sample vs average predicted probabilities of y = 0, 1, ..., 5 
countfit docvis $xlist, maxcount(5) nbreg nograph noestimates nofit

* NB2: Predicted NB2 probabilities at x = x* of y = 0, 1, ..., 5 
quietly nbreg docvis $xlist
prvalue, x(private=1 medicaid=1) max(5) brief

* Generalized negative binomial with alpha parameterized
gnbreg docvis $xlist, lnalpha(female bh) nolog

* Nonlinear least squares
nl (docvis = exp({xb: $xlist one})), vce(robust) nolog

* Hurdle: Pr(y=0)=pi and Pr(y=k)=(1-pi) x Poisson(2) truncated at 0
quietly set obs 10000
set seed 10101             // set the seed !
scalar pi=1-(1-exp(-2))/2  // Probability y=0
generate xhurdle = 0
scalar minx = 0
while minx == 0 {
 generate xph = rpoisson(2) 
 quietly replace xhurdle = xph if xhurdle==0 
 drop xph
 quietly summarize xhurdle
 scalar minx = r(min)
}
replace xhurdle = 0 if runiform() < pi
summarize xhurdle

// Following not included in book
histogram xhurdle, discrete xtitle("Hurdle Poisson")
graph save mus17hurdle, replace 

* Hurdle logit-nb model manually
logit docvis $xlist, nolog

* Second step uses positives only
summarize docvis if docvis > 0

* Zero-truncated negative binomial 
ztnb docvis $xlist if docvis>0, nolog 

// Following not included in book
* Same hurdle logit-nb model using the user-written hnblogit command
hnblogit docvis $xlist, nolog

* margins for marginal effects of first part
quietly logit docvis i.private i.medicaid age age2 educyr i.actlim totchr, vce(robust)
margins, dydx(*) atmean noatlegend

* margins for marginal effects of second part
quietly ztnb docvis i.private i.medicaid age age2 educyr i.actlim totchr if docvis>0, vce(robust)
margins, dydx(*) atmean noatlegend

// Following repeats earlier code as xpois and xnegbin have since been dropped
* Poisson (mu=1) and NB(1,1) generated data
quietly set obs 10000
set seed 10101                   // set the seed !
generate xpois= rpoisson(1) 
generate xg = rgamma(1,1)
generate xnegbin = rpoisson(xg) // NB draws  

* Mixture: Poisson(.5) with prob .9 and Poisson(5.5) with prob .1
set seed 10101        // set the seed !
generate xp1= rpoisson(.5)
generate xp2= rpoisson(5.5)
summarize xp1 xp2
rename xp1 xpmix
quietly replace xpmix = xp2 if runiform() > 0.9
summarize xpmix

* Tabulate the mixture data
tabulate xpmix

// Following not included in book
histogram xpmix, discrete xtitle("Finite mixture Poisson")
graph save mus17pmix, replace

* Compare the four distributions, all with mean 1
graph combine mus17xp.gph mus17negbin.gph mus17pmix.gph ///
mus17hurdle.gph, title("Four different distributions with mean = 1") ///
ycommon xcommon

graph export mus17countdist.eps, replace

* Compare the four distributions, all with mean 1
summarize xpois xnegbin xpmix xhurdle
 
* Finite-mixture model using fmm command with constant probabilities
use mus17data.dta, clear
fmm docvis $xlist, vce(robust) components(2) mixtureof(poisson)

* Predict y for two components
quietly fmm docvis $xlist, vce(robust) components(2) mixtureof(poisson)
predict yfit1, equation(component1)
predict yfit2, equation(component2)
summarize yfit1 yfit2

* Marginal effects for component 1
margins, dydx(*) predict(eq(component1)) atmean noatlegend
* Marginal effects for component 2 
margins, dydx(*) predict(eq(component2)) atmean noatlegend

set scheme s1mono  /* Graphics scheme */
* Create histograms of fitted values  
quietly histogram yfit1, name(_comp_1, replace)
quietly histogram yfit2, name(_comp_2, replace)
quietly graph combine _comp_1 _comp_2

quietly graph export fitted_comp.eps, replace

* 2-component mixture of NB1 
fmm docvis $xlist, vce(robust) components(2) mixtureof(negbin1)

* Fitted values for 2-component NB1 mixture
drop yfit1 yfit2 
predict yfit1, equation(component1)
predict yfit2, equation(component2)
summarize yfit1 yfit2

********* 17.4 EMPIRICAL EXAMPLE 2
 
********* 
drop _all

set more off
use mus17data_z.dta 
* Summary stats for ER use model
use mus17data_z.dta      
global xlist1 age actlim totchr 
summarize er $xlist1
tabulate er 

* NB2 for er
nbreg er $xlist1, nolog

* Sample average fitted probabilities of y = 0 to max() 
prcounts erpr, max(3)
summarize erpr*

* Zero-inflated negative binomial for er
zinb er $xlist1, inflate($xlist1) vuong nolog

* Comparison of NB and ZINB using countfit
countfit er $xlist1, nbreg zinb nograph noestimates

********* 17.5: MODELS WITH ENDOGENOUS REGRESSORS
clear
drop _all 
set more off

* First-stage linear regression
use mus17data.dta
global xlist2 medicaid age age2 educyr actlim totchr
regress private $xlist2 income ssiratio, vce(robust)
predict lpuhat, residual

* Second-stage Poisson with robust SEs
poisson docvis private $xlist2 lpuhat, vce(robust) nolog

* Program and bootstrap for Poisson two-step estimator
program endogtwostep, eclass
  version 11
  tempname b
  capture drop lpuhat2
  regress private $xlist2 income ssiratio
  predict lpuhat2, residual
  poisson docvis private $xlist2 lpuhat2
  matrix `b' = e(b)
  ereturn post `b'
end
bootstrap _b, reps(400) seed(10101) nodots nowarn: endogtwostep

* Command gmm for one-step GMM (nonlinear IV) of overidentified Poisson model
gmm (docvis - exp({xb:private medicaid age age2 educyr actlim totchr}+{b0})), ///
 instruments(income ssiratio medicaid age age2 educyr actlim totchr) onestep nolog

* Command gmm for two-step GMM (nonlinear IV) of overidentified Poisson model
gmm (docvis - exp({xb:private medicaid age age2 educyr actlim totchr}+{b0})), ///
 instruments(income ssiratio medicaid age age2 educyr actlim totchr) twostep nolog

* Test of overidentifying restriction following two-step GMM
estat overid

********** CLOSE OUTPUT
