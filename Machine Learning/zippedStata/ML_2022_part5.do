* ML_2022_part5.do  April 2022 

capture log close

log using ML_2022_part5.txt, text replace

********** OVERVIEW OF ML_2022_part5.do **********

* MACHINE LEARNING AND CAUSAL ANALYSIS


********** SETUP **********

set more off
* version 15
clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** PROGRAM AND DATA DESCRIPTION **********


* Data for inference on suppins example: 5 continuous and 13 binary variables
use mus203mepsmedexp.dta, clear
keep if ltotexp != .
describe ltotexp suppins
summarize ltotexp suppins

* Continuous variables
global xlist2 income educyr age famsze totchr
describe $xlist2
summarize $xlist2

* Discrete binary variables
global dlist2 female white hisp marry northe mwest south ///
    msa phylim actlim injury priolist hvgg
describe $dlist2	
	
* OLS on small model and full model
global rlist2 c.($xlist2)##c.($xlist2) i.($dlist2) c.($xlist2)#i.($dlist2)
qui regress ltotexp suppins $xlist2 $dlist2, vce(robust)
estimates store OLSSMALL
qui regress ltotexp suppins $rlist2, vce(robust)
estimates store OLSFULL
estimates table OLSSMALL OLSFULL, keep(suppins) b(%9.4f) se stats(N df_m r2)

* Partialing-out partial linear model using default plugin lambda
poregress ltotexp suppins, controls($rlist2)
estimates store POREG
estimates table POREG

* Standard heterogeneous effects estimates of ATE
teffects ra (ltotexp $rlist2) (suppins), vce(robust)
estimates store RA
teffects ipw (ltotexp) (suppins $rlist2), vce(robust)
estimates store IPW
teffects aipw (ltotexp $rlist2) (suppins $rlist2), vce(robust)
estimates store AIPW
estimates table RA IPW AIPW, keep(ATE: POmean:)

/* Not needed
teffects aipw (ltotexp $rlist2) (suppins $rlist2), vce(robust)
scalar bAIPW = r(table)[1,1]
scalar seAIPW = r(table)[2,1]
di "Estimates     RA     IPW    AIPW " _n ///
   "          " %7.3f  bRA  %7.3f bIPW %7.3f  bAIPW  _n ///
   "          " %7.3f  seRA  %7.3f seIPW %7.3f  seAIPW 
   */
   
* AIPW using lasso estimate of ATE
telasso (ltotexp $rlist2) (suppins $rlist2), selection(plugin) vce(robust)
estimates store PLUG
telasso (ltotexp $rlist2) (suppins $rlist2), selection(bic) vce(robust)
estimates store BIC
telasso (ltotexp $rlist2) (suppins $rlist2), selection(cv) xfolds(5) ///
   rseed(10101) vce(robust)
estimates store CV
estimates table PLUG BIC CV

* Did not work attempt to just report ATE
teffects ra (ltotexp $rlist2) (suppins), vce(robust)
estimates store RA
estimates table RA, b keep(ATE:1vs0.suppins)

********** CLOSE OUTPUT **********

* log close
* clear 
* exit



