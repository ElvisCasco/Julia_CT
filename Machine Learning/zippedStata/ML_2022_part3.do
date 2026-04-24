* ML_2022_part3.do  April 2021 for Stata version 16

capture log close

log using ML_2022_part3.txt, text replace

********** OVERVIEW OF ML_2022_part3.do **********

* To run you need files
*  mus203mepsmedexp.dta
*  mus228ajr.dta
* in your directory

* And Stata user-written commands 
*   none
* are used

* 2. PARTIALLING OUT ESTIMATOR
* 4. CROSSFIT PARTIALLING OUT ESTIMATOR
* 5. DOUBLE SELECTION ESTIMATOR APPLICATION
* 6. OTHER MODELS

********** SETUP **********

set more off
version 16
clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* File mus203mepsmedexp.dta is aothurs' extract from MEPS 
* (Medical Expenditure Panel Survey)
* for individuals 65 years and older in U.S. Medicare in 2003

* mus228ajr.dta is from Acemoglu, Johnson and Robinson (2001), "The colonial 
* origins of comparative development: an empirical investigation," AER, 1369-1401.

********* 2. PARTIALLING OUT ESTIMATOR 

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

* Lasso information
lassoinfo

* Variables selected
lassoknots, for(ltotexp)
lassoknots, for(suppins)
 
* Various postestimation options
lassocoef (., for(ltotexp))
lassocoef (., for(suppins))
// lassogof does not apply after po, xpo, ds 

* Partialing out done manually
qui lasso linear suppins $rlist2, selection(plugin) 
qui predict suppins_lasso, postselection
qui generate u_suppins = suppins - suppins_lasso
qui lasso linear ltotexp $rlist2, selection(plugin) 
qui predict ltotexp_lasso, postselection
qui generate u_ltotexp = ltotexp - ltotexp_lasso
regress u_ltotexp u_suppins, vce(robust) noconstant noheader

* Cross validation instead
poregress ltotexp suppins, controls($rlist2) selection(cv) rseed(10101)
lassoinfo

********* 4. CROSSFIT PARTIALLING OUT ESTIMATOR 
 
* Crossfit partialing out (double/debiased) using default plugin
xporegress ltotexp suppins, controls($rlist2) rseed(10101) nolog

* Summarize the number of selected variables across the ten folds
lassoinfo

/* 
This takes a long time so is commented out
* Multiple sample splits for xporegress
. xporegress ltotexp suppins, controls($rlist2) rseed(10101) nolog resample(10)
*/

********* 5. DOUBLE SELECTION ESTIMATOR 

* Double selection partial linear model using default plugin lambda
dsregress ltotexp suppins, controls($rlist2)

// Not included - do dsregress manually
qui lasso linear suppins $rlist2, selection(plugin) 
lassocoef
qui lasso linear ltotexp $rlist2, selection(plugin) 
lassocoef
regress ltotexp suppins income age c.income#c.totchr 1.marry#c.income   ///
    0.northe#c.income 1.hvgg#  c.income 1.white#c.educy 0.hisp#c.educyr ///
    0.marry#c.famsze ///
    totchr c.educyr#c.totchr c.age#c.totchr 0.actlim 1.phylim#c.educyr  /// 
    1.priolist#c.educyr 0.phylim#c.famsze 0.actlim#c.famsze             ///
    0.female#c.totchr 1.white#c.totchr 0.hisp#c.totchr 0.hvgg#c.totchr
dsregress ltotexp suppins, controls($rlist2)

********* 6. OTHER MODELS 

*** LOGIT PARTIALLING-OUT ESTIMATOR

* Logit variant of partial linear model and partialing-out estimator
generate dy = totexp > 4000
qui logit dy suppins $rlist2, or vce(robust)
estimates store FULL
qui pologit dy suppins, controls($rlist2) selection(plugin) coef
estimates store PARTIALOUT
qui dslogit dy suppins, controls($rlist2) coef
estimates store DOUBSEL
estimates table FULL PARTIALOUT DOUBSEL, keep(suppins) b(%9.4f) se ///
    stats(N df_m k_controls_sel)

tabulate dy
	
*** INSTRUMENTAL VARIAVBLES APPLICATION

* Read in Acemoglu-Johnson-Robinson data and define globals  
qui use mus228ajr.dta, clear
global xlist lat_abst edes1975 avelf temp* humid* steplow deslow ///
    stepmid desmid drystep  drywint goldm iron silv zinc oilres landlock
describe logpgp95 avexpr logem4
summarize logpgp95 avexpr logem4, sep(0)
correlate logpgp95 avexpr logem4

* Partialling-out IV using plugin for lambda
poivregress logpgp95 (avexpr=logem4), controls($xlist) selection(plugin, hom)

// Not included 
// The following shows selected instrument is logem4
// and selected controls are edes1975 avelf temp2 iron zinc
lassocoef (., for(avexpr))       // chooses logem4 and edes1975 zinc
lassocoef (., for(logpgp95))     // chooses edes1975 avelf
lassocoef (., for(pred(avexpr))) // chooses edes1975 avelf temp2 iron zinc  

********** CLOSE OUTPUT **********

* log close
* clear 
* exit
