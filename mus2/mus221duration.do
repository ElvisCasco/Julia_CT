* mus221duration.do  for Stata 17

capture log close

********** OVERVIEW OF mus221duration.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 21  SURVIVAL ANALYSIS FOR DURATION DATA
* 21.2: DATA AND DATA SUMMARY
* 21.3: SURVIVOR AND HAZARD FUNCTIONS
* 21.4: SEMIPARAMETRIC REGRESSION MODEL
* 21.5: FULLY PARAMETRIC REGRESSION MODELS
* 21.6: MULTIPLE-RECORDS DATA
* 21.7: DISCRETE-TIME HAZARDS LOGIT MODEL
* 21.8: TIME-VARYING REGRESSORS

* To run you need files
*   mus221mccall.dta
* in your directory

* No community-contributed commands are used

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus221mccall.dta is from B.P. McCall (1996), 
* "Unemployment Insurance Rules, Joblessness, and Part-time Work," 
* Econometrica, 64, 647-682.

********** 21.2: DATA AND DATA SUMMARY

* Unemployment spell length: Describe and summarize key variables
qui use mus221mccall, clear
describe spell censor1 ui logwage 
summarize spell censor1 ui logwage 

* List the first six observations of key variables
list spell censor1 ui logwage in 1/6, clean

* Command stset defines the dependent and censoring variables
stset spell, fail(censor1=1)

* Survival description of dataset
stdescribe

* Variation in regressors over time - relevant for multiple-record data
stvary ui logwage

* Summary of survival data by insurance status
stsum, by(ui)

********** 21.3: SURVIVOR AND HAZARD FUNCTIONS

* Graph histogram and density of survival data by ui status
qui graph twoway (hist spell if censor1==1)                   /// 
    (kdensity spell if censor1==1, lwidth(thick) lstyle(p1)), ///
    legend(pos(1) col(1) ring(0)) title("Completed spells")   ///
    saving(graph1.gph, replace)   
qui graph twoway (hist spell if censor1==0)                   /// 
    (kdensity spell if censor1==0, lwidth(thick) lstyle(p1)), ///
    legend(pos(1) col(1) ring(0)) title("Incomplete spells")  ///
    saving(graph2.gph, replace)   

graph combine graph1.gph graph2.gph, iscale(1.2) ycommon xcommon ///
    ysize(2.5) xsize(6) rows(1)

* Compute survivor function 
sts list

* Graph survivor function over all and by ui
qui sts graph, survival ci legend(pos(8) col(1) ring(0)) ///
    saving(graph1.gph, replace)
qui sts graph, by(ui) survival ci legend(pos(8) col(1) ring(0)) ///
    saving(graph2.gph, replace)

graph combine graph1.gph graph2.gph, iscale(1.2) ysize(2.5) xsize(6) rows(1)


* Graph cumulative hazard function and smoothed hazard function
local endash = ustrunescape("\u2013")
qui sts graph, cumhaz ci legend(pos(11) col(1) ring(0)) ///
    saving(graph1.gph, replace) title("Nelson`endash'Aalen cumulative hazard")
qui sts graph, hazard ci legend(pos(8) col(1) ring(0)) ///
    saving(graph2.gph, replace)

* Graph raw hazard rate 
qui sts list, cumhaz saving(cumhazard, replace)
preserve
use cumhazard, clear
qui generate haz = cumhaz - cumhaz[_n-1]
qui graph twoway (scatter haz time) (lfit haz time),           ///
    saving(graph3.gph, replace) legend(pos(11) col(1) ring(0)) ///
    title("Raw hazard estimate")
restore

graph combine graph1.gph graph2.gph graph3.gph, iscale(1.2) ///
    ysize(2) xsize(6) rows(1)

********** 21.4 SEMIPARAMETRIC REGRESSION MODEL

* Cox PH regression with raw coefficients
stcox ui logwage, nohr vce(robust) nolog

* Cox PH regression with exponentiated coefficients
stcox ui logwage, vce(robust) nolog

* PH model curves: Survivor, cumulative hazard, and hazard functions 
qui stcurve, survival at(ui=1) title("Survivor function") ///
    saving(graph1.gph, replace)
qui stcurve, cumhaz at(ui=1) title("Cumulative hazard function")   ///
    saving(graph2.gph, replace)
qui stcurve, hazard at(ui=1) title("Hazard function")     ///
    saving(graph3.gph, replace)

graph combine graph1.gph graph2.gph graph3.gph, iscale(1.2)   ///
    ysize(2) xsize(6) rows(1)

* PH model diagnostics: Check for parallel log-log survival curves
qui stphplot, by(ui) adjust(logwage) legend(pos(1) col(1) ring(0)) ///
    saving(graph1.gph, replace) title("Log`endash'log survival by UI")
qui summarize logwage, d
qui generate highwage = logwage > r(p50)
qui stphplot, by(highwage) adjust(ui) legend(pos(1) col(1) ring(0)) ///
    saving(graph2.gph, replace) title("Log`endash'log survival by wage")

graph combine graph1.gph graph2.gph, iscale(1.2) ysize(2.5) xsize(6) rows(1)

* PH model diagnostics: Check for good prediction of survival curve
qui stcoxkm, by(ui) legend(pos(1) col(1) ring(0)) ///
    saving(graph1.gph, replace) title("Survival predicted by UI")
qui stcoxkm, by(highwage) legend(pos(1) col(1) ring(0)) ///
    saving(graph2.gph, replace) title("Survival predicted by wage")

graph combine graph1.gph graph2.gph, iscale(1.2) ysize(2.5) xsize(6) rows(1)

* PH model diagnostics: Test of PH assumption
qui stcox ui logwage, vce(robust) nolog
estat phtest, detail

* PH model diagnostics: Graph the Schoenfeld residual against time
qui stcox ui logwage, vce(robust) nolog
qui predict double rui rlogwage, schoenfeld 
graph twoway (scatter rui _t) (lfit rui _t),  legend(off) ///
    ytitle("Schoenfeld residual for ui") legend(off) saving(graph1.gph, replace) 
graph twoway (scatter rlogwage _t) (lfit rlogwage _t),  legend(off) ///
    ytitle("Schoenfeld residual for logwage") saving(graph2.gph, replace) 

graph combine graph1.gph graph2.gph, iscale(1.2) ysize(2.5) xsize(6) rows(1)

* Semiparametric competing risks example - full-time job versus not full-time
qui generate event = 1 if censor1 == 1
qui replace event = 2 if (censor2 == 1 | censor3 == 1)
qui stset spell, fail(event=1) 
stcrreg ui logwage, compete(event == 2) nolog vce(robust)

// not included
stcurve, cif at1(ui=0) at2(ui=1)
gen one = 1
stcrreg one, compete(event == 2) nolog vce(robust)
stcurve, cif

********** 21.5 FULLY PARAMETRIC REGRESSION MODELS

* Parametric Weibull regression 
stset spell, fail(censor1=1) 
streg ui logwage, nohr vce(robust) dist(weibull) nolog

// Not included - Weibull with default standard errors
streg ui logwage, nohr dist(weibull) nolog

* Parametric Weibull model: Prediction of the conditional mean duration time
qui streg i.ui logwage, nohr vce(robust) dist(weibull) nolog
predict muspell, mean time
qui generate completedspell = spell if censor1==1
qui generate mucompletedspell = muspell if censor1==1
summarize muspell completedspell mucompletedspell

* Parametric Weibull model: AMEs for the conditional mean duration time
margins, dydx(*) predict(mean time)

* Parametric models (plus Cox) with PH parameterization 
qui streg ui logwage, dist(exponential) nohr vce(robust) 
qui stcurve, hazard title("Exponential") saving(graph1.gph, replace)
estimates store Exponential
qui streg ui logwage, dist(weibull) nohr vce(robust)
qui stcurve, hazard title("Weibull") saving(graph2.gph, replace)
estimates store Weibull
qui streg ui logwage, dist(gompertz) nohr vce(robust)
qui stcurve, hazard title("Gompertz") saving(graph3.gph, replace)
estimates store Gompertz
qui stcox ui logwage, nohr vce(robust)
qui stcurve, hazard title("Cox") saving(graph4.gph, replace)
estimates store Cox

* Estimates of parametric models with PH parameterization 
estimates table Exponential Weibull Gompertz Cox, eq(1) b(%11.3f) se ///
    stats(ll aic bic) 

* Parametric models with AFT parameterization 
qui streg ui logwage, dist(loglogistic) vce(robust)
qui stcurve, hazard title("Loglogistic") saving(graph5.gph, replace)
estimates store Loglogistic
qui streg ui logwage, dist(lognormal) vce(robust)
qui stcurve, hazard title("Lognormal") saving(graph6.gph, replace)
estimates store Lognormal
qui streg ui logwage, dist(exponential) nohr vce(robust) time
estimates store ExponAFT
qui streg ui logwage, dist(weibull) nohr vce(robust) time
estimates store WeibullAFT

* Estimates of parametric models with AFT parameterization 
estimates table Loglogistic Lognormal ExponAFT WeibullAFT, eq(1) b(%11.3f) ///
   se stats(ll aic bic) 

* Fitted hazard function for six parametric models
graph combine graph1.gph graph2.gph graph3.gph graph4.gph  ///
    graph5.gph graph6.gph, iscale(0.9) ysize(4) xsize(6)   ///
    rows(2) ycommon xcommon

* Predicted means and AMEs for four parametric models
foreach dist in exponential weibull loglogistic lognormal {
    qui streg i.ui logwage, dist(`dist') vce(robust) 
    qui margins, predict(mean time)
    display "Model `dist':" _col(21) "Ave mean time =" %6.2f r(b)[1,1] 
    qui margins, dydx(*) predict(mean time) 
    dis _col(21) "AME ui =" %6.2f r(b)[1,1] %6.2f " AME logwage= " r(b)[1,2]
    }

* Weibull model with gamma frailty
streg ui logwage, dist(weibull) frailty(invgau) nolog nohr vce(robust)

* Hazard curves conditional and unconditional on alpha
qui stcurve, hazard alpha1 ///
    title("Conditional on {&alpha}=1") saving(graph1.gph, replace)
qui stcurve, hazard unconditional ///
    title("Unconditional on {&alpha}") saving(graph2.gph, replace)

* AME at median for Weibull without and with frailty
qui streg 1.ui logwage, dist(weibull) vce(robust)
margins, dydx(*)
qui streg 1.ui logwage, dist(weibull) frailty(invgau) vce(robust)
margins, dydx(*)

graph combine graph1.gph graph2.gph, iscale(1.2) ysize(2.5) xsize(6) ///
    rows(1) xcommon

// Not included - gamma frailty  
streg ui logwage, dist(lognormal) frailty(gamma) nolog vce(robust)

********** 21.6 MULTIPLE-RECORDS DATA

* Expand data to have # observations per person = spell length
qui use mus221mccall, clear
generate id = _n                     // Need an individual identifier
stset spell, fail(censor1=1) id(id)  // stset must include id() option
stsplit t, every(1)             

* Data created by stsplit
summarize id spell censor1 ui logwage _* t
list id spell censor1 ui logwage _* t if id==9, clean
list id spell censor1 ui logwage _* t if id==10, clean 

* Cox PH with multiple-records data gives same results as single-record data
stcox ui logwage, nolog vce(robust) nohr

stcox ui logwage, nolog vce(robust)
stcox ui logwage, nolog vce(cluster id)
stcox ui logwage, nolog 

********** 21.7: DISCRETE-TIME HAZARDS MODEL

* Discrete-time hazards: Create indicator variable for getting employed
generate demploy = 0
replace demploy = 1 if censor1 == 1

* Discrete-time hazards: logit with quadratic time trend 
logit demploy ui logwage c.t##c.t, vce(cluster id) nolog

* Discrete-time hazards: cloglog with time dummies
cloglog demploy ui logwage i.t, vce(cluster id) nolog

* Discrete-time hazards logit with time dummies
logit demploy ui logwage i.t, vce(cluster id) nolog

* Compare AMEs for discrete-time hazards cloglog and logit
qui cloglog demploy 1.ui logwage i.t, vce(cluster id) nolog
margins, dydx(1.ui logwage)
qui logit demploy 1.ui logwage i.t, vce(cluster id) nolog
margins, dydx(1.ui logwage)

********** 21.8 TIME-VARYING REGRESSORS 

* Cox PH regression with time-varying regressors
generate tvlogwage = logwage + 0.5*rnormal(0,1)   // Create time-varying x
stcox ui tvlogwage, vce(robust) nohr nolog

* Not included. Test quadratic versus time dummies
generate tsq = t^2
logit demploy ui logwage t tsq i.t, vce(cluster id) nolog
testparm i(1/27).t    // Requires knowing that max(t) = 27

streg ui logwage, dist(weibull) nohr vce(robust) 

streg ui logwage, dist(weibull) nohr vce(robust) time


*** Erase files created by this program and not used elsewhere in the book

erase cumhazard.dta
erase graph1.gph
erase graph2.gph
erase graph3.gph
erase graph4.gph
erase graph5.gph
erase graph6.gph

********** END
