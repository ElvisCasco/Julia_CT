*mus14p1bin.do  Oct 2009 for Stata version 11

drop _all 
cap log close

********** OVERVIEW OF mus14bin.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
  
* 14.4 EXAMPLE
* 14.5 HYPOTHESIS AND SPECIFICATION TESTS
* 14.6 GOODNESS OF FIT AND PREDICTION
* 14.7 MARGINAL EFFECTS 
* 14.8 ENDOGENOUS REGRESSORS
* 14.9 GROUPED DATA

* To run this program you need data file 
*   mus14data.dta 
* in your directory
* No Stata user-written commands are used

* Dataset comes from HRS 2000

********** SETUP

set more off
version 11
clear all
set scheme s1mono  /* Graphics scheme */
  
********** 14.3 DATA DESCRIPTION

* Dataset comes from HRS 2000
 
*********** 14.4 EXAMPLE

clear
* Load data
use mus14data.dta
* Interaction variables
drop age2 agefem agechr agewhi
* Summary statistics of variables
global xlist age hstatusg hhincome educyear married hisp
generate linc = ln(hhinc)
global extralist linc female white chronic adl sretire
summarize ins retire $xlist $extralist

* 14.4.2 Logit regression

* Logit regression
logit ins retire $xlist
 
* Estimation of several models
quietly logit ins retire $xlist
estimates store blogit
quietly probit ins retire $xlist 
estimates store bprobit
quietly regress ins retire $xlist 
estimates store bols
quietly logit ins retire $xlist, vce(robust)
estimates store blogitr
quietly probit ins retire $xlist, vce(robust)
estimates store bprobitr
quietly regress ins retire $xlist, vce(robust)
estimates store bolsr

* Table for comparing models 
estimates table blogit blogitr bprobit bprobitr bols bolsr, /*
   */ t stats(N ll) b(%7.3f) stfmt(%8.2f)

********** 14.5 HYPOTHESIS AND SPECIFICATION TESTS

* Wald test for zero interactions
generate age2 = age*age
generate agefem = age*female
generate agechr = age*chronic
generate agewhi = age*white
global intlist age2 agefem agechr agewhi
quietly logit ins retire $xlist $intlist
test $intlist 

* Likelihood-ratio test
quietly logit ins retire $xlist $intlist
estimates store B 
quietly logit ins retire $xlist
lrtest B 

*** 14.5.3  HETEROSKEDASTIC PROBIT
* Heteroskedastic probit model
hetprob ins retire $xlist, het(chronic) nolog // Heteroskedastic Probit 

* Stukel score or LM test for asymmetric h-family logit
quietly logit ins retire $xlist
predict xbhat, xb
generate xbhatsq = xbhat^2
quietly logit ins retire $xlist xbhatsq
test xbhatsq


********** 14.6 GOODNESS OF FIT AND PREDICTION

* Hosmer-Lemeshow gof test with 4 groups
quietly logit ins retire $xlist
estat gof, group(4)  // Hosmer-Lemeshow gof test

quietly logit ins retire $xlist
* Hosmer-Lemeshow gof test with 10 groups
estat gof, group(10)  // Hosmer-Lemeshow gof test


**FIGURE 1: PLOT PREDICTED PROBABILITY AGAINST hhincome FOR MODELS

* Calculate and summarize fitted probabilities
quietly logit ins hhincome
predict plogit, pr
quietly probit ins hhincome  
predict pprobit, pr
quietly regress ins hhincome
predict pols, xb
summarize ins plogit pprobit pols

* Following gives Figure mus14fig1.eps
sort hhincome
graph twoway (scatter ins hhincome, msize(vsmall) jitter(3)) /*
  */ (line plogit hhincome, clstyle(p1)) /*
  */ (line pprobit hhincome, clstyle(p2)) /*
  */ (line pols hhincome, clstyle(p3)), /*
  */ scale (1.2) plotregion(style(none)) /*
  */ title("Predicted Probabilities Across Models") /*
  */ xtitle("HHINCOME (hhincome)", size(medlarge)) xscale(titlegap(*5)) /* 
  */ ytitle("Predicted probability", size(medlarge)) yscale(titlegap(*5)) /*
  */ legend(pos(1) ring(0) col(1)) legend(size(small)) /*
  */ legend(label(1 "Actual Data (jittered)") label(2 "Logit") /*
  */         label(3 "Probit") label(4 "OLS"))
graph export mus14fig1.eps, replace

* Fitted probabilities for selected baseline
quietly logit ins retire $xlist
prvalue, x(age=65 retire=0 hstatusg=1 hhincome=50 educyear=17 married=1 hisp=0)

prchange hhincome

* Comparing fitted probability and dichotomous outcome
quietly logit ins retire $xlist
estat classification

********** 14.7 MARGINAL EFFECTS 

* Marginal effects (MER) after logit
quietly logit ins i.retire age i.hstatusg hhincome educyear i.married i.hisp
margins, dydx(*) at (retire=1 age=75 hstatusg=1 hhincome=35 educyear=12   ///
 married=1 hisp=1) noatlegend   // (MER)

* Marginal effects (MEM) after logit
quietly logit ins i.retire age i.hstatusg hhincome educyear i.married i.hisp
margins, dydx(*) atmean noatlegend  // (MEM)

* Marginal effects (AME) after logit
quietly logit ins i.retire age i.hstatusg hhincome educyear i.married i.hisp
margins, dydx(*) noatlegend        // (AME)

* Computing change in probability after logit
quietly logit ins retire $xlist 
prchange hhincome

quietly logit ins retire $xlist
estat classification

********** 14.8  ENDOGENOUS REGRESSORS

use mus14data.dta, clear

* Endogenous probit using inconsistent probit MLE
generate linc = log(hhincome)
global xlist2 female age age2 educyear married hisp white chronic adl hstatusg
probit ins linc $xlist2, vce(robust) nolog

* Endogenous probit using ivprobit ML estimator
global ivlist2 retire sretire
ivprobit ins $xlist2 (linc = $ivlist2), vce(robust) nolog

* Endogenous probit using ivprobit 2-step estimator
ivprobit ins $xlist2 (linc = $ivlist2), twostep first

* Endogenous probit using ivregress to get 2SLS estimator
ivregress 2sls ins $xlist2 (linc = $ivlist2), vce(robust) noheader
estat overid


* 2SLS estimation
quietly ivregress 2sls ins $xlist (retire= female sretire seprhi adl), robust 
estimates store TSLS_IV4

********** 14.9 GROUPED DATA

* Using mus14data.dta to generate grouped data 
sort age
collapse av_ret=retire av_hhinc=hhincome av_educyear=educyear	      ///
	av_mar=married av_adl=adl av_hisp=hisp av_hstatusg=hstatusg   ///
	av_ins=ins, by(age)
generate logins = log(av_ins/(1-av_ins))
save mus14gdata.dta, replace

*Summarize grouped data
summarize logins  $xlistg

* Regressions with grouped data
regress logins av_ret av_hstatusg av_hhinc av_educyear av_mar av_hisp, vce(robust)

********** CLOSE OUTPUT **********
