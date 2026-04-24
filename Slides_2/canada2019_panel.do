* canada2019_panel.do based on trcountimf2010.do May 2018 for Stata version 15
* and after that based on count2015sweden.do 
* and also racd29.do

capture log close
log using canada2019_panel.txt, text replace

********** OVERVIEW OF canada2019_panel.do **********

* STATA Program to demonstrate Cound Data Regression Models
* Based on mus17p1cnt.do
* A. Colin Cameron and Pravin K. Trivedi (2009)
* "Microeconometrics using Stata", Stata Press

* 1A. COUNT REGRESSION: BASIC CROSS-SECTION
* 1B. COUNT REGRESSION: INFERENCE (STANDARD ERRORS AND BOOTSTRAP) 
* 2: COUNT REGRESSION: EXTRAS (TRUNCATED/CENSORED, HURDLE, ZERO-INFLATED, 
*                              FINITE MIXTURE, ENDOGENOUS)
* 3: COUNT REGRESSION: BASIC PANEL

* To run you need files
*      mus18data.dta
*      racd09data.dta
* in your directory

********** SETUP **********

set more off
version 15
clear all
set scheme s1mono  /* Graphics scheme */

*********** 3: PANEL COUNT REGRESSION

* Describe dependent variables and regressors
use mus18data.dta, clear
describe mdu lcoins ndisease female age lfam child id year

* Summarize dependent variables and regressors
summarize mdu lcoins ndisease female age lfam child id year

*** PANEL DATA SUMMARY 

* Panel description of dataset 
xtset id year
xtdescribe

* Panel summary of dependent variable
xtsum mdu

* Year-to-year transitions in doctor visits
generate mdushort = mdu
replace mdushort = 4 if mdu >= 4
xttrans mdushort
corr mdu L.mdu

* Simple time-series plot for first 20 individuals (= first 85 obs)
quietly xtline mdu if _n<=85, overlay legend(off)
graph export countfig31panelplot.wmf, replace

* Autocorrelations - not the best way 
sort id year
corr mdu L.mdu L2.mdu L3.mdu

* Autocorrelations using all available data
pwcorr mdu L.mdu L2.mdu L3.mdu L4.mdu

*** 3. PANEL POISSON: POOLED AND POPULATION AVERAGED

* Pooled Poisson estimator with cluster-robust standard errors
poisson mdu lcoins ndisease female age lfam child, vce(cluster id)

* Without cluster-robust
poisson mdu lcoins ndisease female age lfam child

* Poisson PA estimator with unstructured error correlation and robust VCE
xtpoisson mdu lcoins ndisease female age lfam child, pa corr(unstr) vce(robust)

* Print out the correlation matrix
matrix list e(R)

* Allow for AR(1) correlation in y
xtpoisson mdu lcoins ndisease female age lfam child, pa corr(ar1) vce(robust)

*** 3. PANEL POISSON: RANDOM EFFECTS

* Poisson random-effects estimator (gamma) with default standard errors
xtpoisson mdu lcoins ndisease female age lfam child, re

* Poisson random-effects estimator (gamma) with cluster-robust standard errors
xtpoisson mdu lcoins ndisease female age lfam child, re vce(cluster id)
* Preceding same as xtpoisson , re vce(robust)

* Poisson random-effects estimator with normal intercept and default standard errors
xtpoisson mdu lcoins ndisease female age lfam child, re normal

* Poisson random-effects estimator with normal intercept and normal slope for one parameter

*** 3. PANEL POISSON: MIXED MODELS

* These take a very long time and are commented out

/* 
* Intercept only varies
meqrpoisson mdu lcoins ndisease female age lfam child || id: 
* Previous is the same as 
mepoisson mdu lcoins ndisease female age lfam child || id: 
* And is the same as 
xtpoisson mdu lcoins ndisease female age lfam child, re normal

* Now also allow a slope parameter to vary
mepoisson mdu lcoins ndisease female age lfam child || id: ndisease

* Now also allow slope parameter to covary with intercept parameter
mepoisson mdu lcoins ndisease female age lfam child || id: ndisease, ///
    cov(unstructured) intpoints(9)

*/

*** 3. PANEL POISSON: FIXED EFFECTS

* Poisson fixed-effects estimator with WRONG default standard errors
xtpoisson mdu lcoins ndisease female age lfam child, fe

* Poisson fixed-effects estimator with panel-robust standard errors
xtpoisson mdu lcoins ndisease female age lfam child, fe vce(robust)
* This supplants earlier add-on xtpqml
* xtpqml mdu lcoins ndisease female age lfam child, fe i(id) cluster(id)

*** 3. STATIC PANEL POISSON: ESTIMATOR COMPARISON

* Comparison of Poisson panel estimators using panel-robust standard errors
global xlist lcoins ndisease female age lfam child
quietly poisson mdu $xlist, vce(cluster id)
estimates store POOLED
quietly xtpoisson mdu $xlist, pa corr(unstr) vce(robust)
estimates store POPAVE
quietly xtpoisson mdu $xlist, re vce(cluster id)
estimates store RE_GAMMA
quietly xtpoisson mdu $xlist, re normal vce(cluster id)
estimates store RE_NORMAL
quietly xtpoisson mdu $xlist, fe vce(robust)
estimates store FIXED
estimates table POOLED POPAVE RE_GAMMA RE_NORMAL FIXED, equations(1) ///
   b(%8.4f) se stats(N ll) stfmt(%8.0f)

/* In early versions of Stata xtpoisson re  does not do panel-robust se's and instead 
   bootstrap as follows 
xtpoisson mdu lcoins ndisease female age lfam child, re vce(boot, reps(100) seed(10101))
*/

*** 3. DYNAMIC PANEL POISSON: ESTIMATOR COMPARISON

* Variable descriptions and summary statistics
use racd09data.dta, clear
global XLISTD PAT1 LOGR LOGR1 LOGR2 LOGK SCISECT dyear2 dyear3 dyear4 dyear5
describe PAT $XLISTD
summarize PAT1

* Summarize all
summarize $XLISTD

*** 3. DYNAMIC MODELS USING EXPONENTIAL FEEDBACK MODEL

*** 3. NON-FIXED EFFECTS DYNAMIC MODELS 

* Pooled Poisson
poisson PAT $XLISTD, vce(cluster id) nolog
estimates store DPCS
display "Table 9.5: first column Sum ln R"
lincom LOGR + LOGR1 + LOGR2

* Population Averaged Poisson with exchangeable errors
* Following standard errors are preferred
xtpoisson PAT $XLISTD, pa vce(robust)
display "Table 9.5: second column Sum ln R"
lincom LOGR + LOGR1 + LOGR2

* Poisson Random Effects - gamma
estimates store DPPA
xtpoisson PAT $XLISTD, re vce(robust)
estimates store DPRE
display "Table 9.5: third column Sum ln R"
lincom LOGR + LOGR1 + LOGR2

generate PAT_YEAR0 = PAT1 if YEAR==1    // = PAT1 in YEAR 1 and missing in YEARS 1-5
bysort id: egen PAT_INITIAL = mean(PAT_YEAR0)  // Replaces missings with PAT1 in YEAR1

* Conditionally Correlated Random Effects
sort id
by id: egen LOGRMEAN = mean(LOGR)
by id: egen LOGR1MEAN = mean(LOGR1)
by id: egen LOGR2MEAN = mean(LOGR2)
by id: egen LOGR3MEAN = mean(LOGR3)
by id: egen LOGR4MEAN = mean(LOGR4)
by id: egen LOGR5MEAN = mean(LOGR5)
global MEANS LOGRMEAN LOGR1MEAN LOGR2MEAN LOGR3MEAN LOGR4MEAN LOGR5MEAN

global XLISTD2 PAT1 LOGR LOGR1 LOGR2 LOGK SCISECT PAT_INITIAL LOGRMEAN LOGR1MEAN ///
   LOGR2MEAN dyear2 dyear3 dyear4 dyear5 

* Correlated random effects versions of the same
* Cross-section Poisson
poisson PAT $XLISTD2, vce(cluster id)
estimates store DPCS2
* Population averaged Poisson with exchangeable errrors
xtpoisson PAT $XLISTD2, pa vce(robust)
estimates store DPPA2

* Poisson Random Effects - gamma
xtpoisson PAT $XLISTD2, re vce(robust)
estimates store DPCCRE
display "Table 9.6: first column Sum ln R"
lincom LOGR + LOGR1 + LOGR2

estimates table DPCS DPCS2 DPPA DPPA2 DPRE DPCCRE, b(%7.4f) se(%7.3f) ///
   stats(N ll) stfmt(%9.1f) modelwidth(9) equations(1)

*** 3. FIXED EFFECTS DYNAMIC MODELS

* Fixed effects GMM using Chamberlain transformation
* This program is the same as gmm_poipre in Stata manual [r]gmm
program gmm_poipre
   version 11
   syntax varlist if, at(name) myrhs(varlist) ///
   mylhs(varlist) myidvar(varlist)
   quietly {
   tempvar mu mubar ybar
   gen double `mu' = 0 `if'
   local j = 1
   foreach var of varlist `myrhs' {
      replace `mu' = `mu' + `var'*`at'[1,`j'] `if'
      local j = `j' + 1
      }
   replace `mu' = exp(`mu')
   replace `varlist' = L.`mylhs' - L.`mu'*`mylhs'/`mu' `if'
   }
end

* Only include time-varying regressors
* Also here year 1 is dropped, so drop the year 2 dummy
* Regressors 
global XLISTTV PAT1 LOGR LOGR1 LOGR2 dyear3 dyear4 dyear5
* Instruments in just-identified case
global IVLISTTV PAT2 LOGR1 LOGR2 LOGR3 dyear3 dyear4 dyear5
* Instruments in over-identified case
global IVLISTTV2 PAT2 PAT3 PAT4 LOGR1 LOGR2 LOGR3 dyear3 dyear4 dyear5

* Just-identified
gmm gmm_poipre, mylhs(PAT) myrhs($XLISTTV) myidvar(id) nequations(1) ///
   parameters($XLISTTV) instruments($IVLISTTV, noconstant) onestep vce(cluster id)
estimates store DPGMM

* Overidentified
gmm gmm_poipre, mylhs(PAT) myrhs($XLISTTV) myidvar(id) nequations(1) ///
  parameters($XLISTTV) instruments($IVLISTTV2, noconstant) twostep vce(cluster id)
estimates store DPGMMOID
display "Table 9.6: second column Sum ln R"
lincom _b[LOGR:_cons]+_b[LOGR1:_cons]+_b[LOGR2:_cons]
estat overid
* predict residGMM
* correlate residGMM L.residGMM L2.residGMM

*** CHECK: THIS DOES POISSON FE USING GMM COMMAND

* This program is the same as gmm_poi in Stata manual [r]gmm
program gmm_poi2
   version 11
   syntax varlist if, at(name) myrhs(varlist) ///
   mylhs(varlist) myidvar(varlist)
   quietly {
   tempvar mu mubar ybar
   gen double `mu' = 0 `if'
   local j = 1
   foreach var of varlist `myrhs' {
      replace `mu' = `mu' + `var'*`at'[1,`j'] `if'
      local j = `j' + 1
      }
   replace `mu' = exp(`mu')
   egen double `mubar' = mean(`mu') `if', by(`myidvar')
   egen double `ybar' = mean(`mylhs') `if', by(`myidvar')
   replace `varlist' = `mylhs' - `mu'*`ybar'/`mubar' `if'
   }
end
global XLISTTIMEVARYING LOGR LOGR1 LOGR2 LOGR3 LOGR4 LOGR5 dyear2 dyear3 dyear4 dyear5
gmm gmm_poi2, mylhs(PAT) myrhs($XLISTTIMEVARYING)          ///
 myidvar(id) nequations(1) parameters($XLISTTIMEVARYING)   ///
 instruments($XLISTTIMEVARYING, noconstant) onestep vce(cluster id)
estimates store PFEGMM


********** CLOSE OUTPUT

* log close
* clear
* exit

