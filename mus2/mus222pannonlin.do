 * mus222pannonlin.do  for Stata 17

cap log close

********** OVERVIEW OF mus222pannonlin.do **********
* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press  
 
* Chapter 22  NONLINEAR PANEL MODELS
* 22.3 NONLINEAR PANEL-DATA EXAMPLE
* 22.4 BINARY OUTCOME AND ORDERED OUTCOME MODELS
* 22.5 TOBIT AND INTERVAL-DATA MODELS
* 22.6 COUNT-DATA MODELS
* 22.7 PANEL QR

* To run you need files
*   mus218rhie.dta    
*   mus208psid.dta
* in your directory

* community-contributed commands
*   logitfe
*   qreg2
* are used

********** SETUP **********

clear all
set linesize 80
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus218rhie.dta    
* Rand Health Insurance Experiment data 
* Essentially same data as in P. Deb and P.K. Trivedi (2002)
* "The Structure of Demand for Medical Care: Latent Class versus
* Two-Part Models", Journal of Health Economics, 21, 601-625
* except that article used different outcome (counts rather than $)
* Each observation is for an individual over a year.
* Individuals may appear in up to five years.
* All available sample is used except only fee for service plans included.
* If panel data used then clustering is on id (person id)

* mus210psid.dta
* PSID. Same as Stata website file psidextract.dta
* Data due to  Baltagi and Khanti-Akom (1990) 
* This is corrected version of data in Cornwell and Rupert (1988).
* 595 individuals for years 1976-82

********** CHAPTER 22.3 NONLINEAR PANEL DATA EXAMPLE

* Describe dependent variables and regressors
qui use mus218rhie
describe dmdu med mdu lcoins ndisease female age lfam child id year

* Summarize dependent variables and regressors
summarize dmdu med mdu lcoins ndisease female age lfam child id year

* Panel description of dataset 
xtset id year
xtdescribe 

* Panel summary of time-varying regressors
xtset id year 
xtsum age lfam child

********** CHAPTER 22.4 BINARY OUTCOME AND ORDERED OUTCOME MODELS

* Logit: Panel summary of dependent variable
xtsum dmdu

* Year-to-year transitions in whether visit doctor
xttrans dmdu

* Correlations over time in the dependent variable
pwcorr dmdu l.dmdu l2.dmdu l3.dmdu l4.dmdu

* Logit cross-section with panel-robust standard errors
logit dmdu lcoins ndisease female age lfam child, vce(cluster id) nolog

* Pooled logit cross-section with exchangeable errors and panel-robust VCE
xtlogit dmdu lcoins ndisease female age lfam child, pa corr(exch) ///
    vce(robust) nolog

* Within correlation of the exchangeable errors
matrix list e(R)

* Logit RE estimator
xtlogit dmdu lcoins ndisease female age lfam child, re nolog vce(robust)

* Logit FE estimator
xtlogit dmdu lcoins ndisease female age lfam child, fe nolog 

* Logit mixed-effects estimator (same as xtlogit, re)
* xtmelogit dmdu lcoins ndisease female age lfam child || id:

// Following code appears in text with output omitted
* Logit-correlated RE estimator
bysort id: egen aveage = mean(age)
by id: egen avelfam = mean(lfam)
by id: egen avechild = mean(child)
xtlogit dmdu lcoins ndisease female age lfam child ///
    aveage avelfam avechild, re vce(robust) nolog    
estimates store CRE

// logitfe takes a long time

* Panel logit estimator comparison
global xlist lcoins ndisease female age lfam child
qui logit dmdu $xlist, vce(cluster id)
estimates store POOLED
qui xtlogit dmdu $xlist, pa corr(exch) vce(robust)
estimates store PA
qui xtlogit dmdu $xlist, re vce(robust)    
estimates store RE
qui xtlogit dmdu $xlist, fe         // vce(robust) not available
estimates store FE
qui logitfe dmdu $xlist, teffects(no) nocorrection
estimates store FEDV
qui logitfe dmdu $xlist, teffects(no) analytical
estimates store FEDVBC
qui xtlogit dmdu $xlist aveage avelfam avechild, re vce(robust) 
estimates store CRE

* Panel logit estimator comparison results
estimates table PA RE FE FEDV FEDVBC CRE, equations(1) se b(%7.4f) ///
    stats(N ll) stfmt(%7.0f) varwidth(7)

* AMEs for PA, RE, and FE estimates
foreach model in pa re fe {
    qui xtlogit dmdu $xlist, `model'  
    qui margins, dydx(*)
    display "Model `model':" _col(11) "lcoins" _col(21) "ndisease"    ///
	    _col(31) "female" _col(41) "age" _col(51) "lfam" _col(61) "child" 
    display _col(11) %7.4f r(b)[1,1] _col(21) %7.4f r(b)[1,2]         /// 
            _col(31) %7.4f r(b)[1,3] _col(41) %7.4f r(b)[1,4]         ///
            _col(51) %7.4f r(b)[1,5] _col(61) %7.4f r(b)[1,6] 
    }

* AMEs for bias-corrected dummy variable FE estimator 
logitfe dmdu age lfam child, teffects(no) analytical

* AMEs for CREs estimator 
qui xtlogit dmdu $xlist aveage avelfam avechild, re vce(robust)
margins, dydx(*) 

// Not included
qui xtlogit dmdu $xlist aveage avelfam avechild, re vce(robust)
margins, dydx(*) predict(pu0)

// Not included
* Marginal effects (AME) in different time periods for RE estimates
qui xtlogit dmdu $xlist, pa vce(robust)  
forvalues i = 1/5 {
    qui margins if year==`i' , dydx(*)
    display "year `year':" _col(11) "lcoins" _col(21) "ndisease"    ///
	    _col(31) "female" _col(41) "age" _col(51) "lfam" _col(61) "child" 
    display _col(11) %7.4f r(b)[1,1] _col(21) %7.4f r(b)[1,2]       /// 
            _col(31) %7.4f r(b)[1,3] _col(41) %7.4f r(b)[1,4]       ///
            _col(51) %7.4f r(b)[1,5] _col(61) %7.4f r(b)[1,6] 
    }

* Tabulate cmdu generated by recoding count into seven ordered categories
recode mdu (0 = 1) (1/2 = 2) (3/4 = 3) (5/6 = 4) (7/8 = 5) ///
    (8/12 = 6) (13/999 = 7), gen(cmdu)
tabulate cmdu

* RE ordered logit: Estimates
xtologit cmdu lcoins ndisease female age lfam child, vce(cluster id) nolog

* RE ordered logit: AMEs of coinsurance rate on probability of category outcomes
margins, dydx(lcoins)

* RE ordered logit: Estimate fitted average category probabilities
predict pr*, pr 
summarize pr*

********** CHAPTER 22.5 TOBIT AND INTERVAL-DATA MODELS

* Tobit: Panel summary of dependent variable
xtsum med

* Tobit RE estimator
xttobit med lcoins ndisease female age lfam child, ll(0) nolog

// The following takes time
* Heckman sample-selection panel estimator
gen doutpdol = outpdol > 0
xtheckman outpdol lcoins ndisease age if id > 629388, intpoints(10) nolog ///
    select(doutpdol = lcoins ndisease age female lfam child hlthf hlthp linc) 

* Predictive margins for outpdol given outpdol > 0 at different ages 
margins, at(age=(20(40)60)) predict(ycond)

********** CHAPTER 22.6 COUNT-DATA MODELS

* Poisson: Panel summary of dependent variable
xtsum mdu

* Year-to-year transitions in doctor visits
generate mdushort = mdu
replace mdushort = 4 if mdu >= 4
xttrans mdushort

* Correlations over time in the dependent variable 
pwcorr mdu L1.mdu L2.mdu L3.mdu L4.mdu 

* Pooled Poisson estimator with cluster–robust standard errors
poisson mdu lcoins ndisease female age lfam child, vce(cluster id)

* Poisson PA estimator with unstructured error correlation and robust VCE
xtpoisson mdu lcoins ndisease female age lfam child, pa corr(unstr) vce(robust)

* Correlations over time 
matrix list e(R)

* Poisson RE estimator with default standard errors
xtpoisson mdu lcoins ndisease female age lfam child, re

* Poisson RE estimator with cluster–robust standard errors
xtpoisson mdu lcoins ndisease female age lfam child, re vce(robust) nolog

// Not included
* Poisson random-effects estimator with normal intercept and default standard errors
* xtpoisson mdu lcoins ndisease female age lfam child, re normal

// Not included
* Poisson RE estimator with normal intercept and normal slope for one parameter
* xtmepoisson mdu lcoins ndisease female age lfam child || id: NDISEASE

* Poisson FE estimator with default standard errors
xtpoisson mdu lcoins ndisease female age lfam child, fe i(id)  

* Poisson FE estimator with cluster–robust standard errors
xtpoisson mdu lcoins ndisease female age lfam child, fe vce(robust)

* Compute various Poisson panel estimators
global xlist age lfam child lcoins ndisease female
qui xtpoisson mdu $xlist, pa corr(unstr) vce(robust)
estimates store PPA_ROB
qui xtpoisson mdu $xlist, re vce(robust)
estimates store PRE
qui xtpoisson mdu $xlist, re normal vce(robust)
estimates store PRE_NORM
qui xtpoisson mdu $xlist, fe vce(robust)
estimates store PFE
qui xtpoisson mdu $xlist aveage avelfam avechild, re normal vce(robust)
estimates store CRE

* Report various Poisson panel estimators - results
estimates table PPA_ROB PRE PRE_NORM PFE CRE, equations(1) b(%8.4f) ///
   se stats(N ll) stfmt(%8.0f)

* Compute AMEs for PA, RE, FE, and CREs models
qui xtpoisson mdu $xlist, pa corr(unstr) vce(robust)
qui margins, dydx(*) predict(mu)
matrix PA = r(b)
qui xtpoisson mdu $xlist, re vce(robust) 
qui margins, dydx(*) expression(exp(predict(xb)))
matrix RE = r(b)
qui xtpoisson mdu $xlist, re normal vce(robust)
qui margins, dydx(*) predict(n)
matrix REnorm = r(b)
qui xtpoisson mdu $xlist, fe vce(robust)
qui margins, dydx(*) predict(nu0)
matrix FE_0 = r(b)
qui xtpoisson mdu $xlist aveage avelfam avechild, re normal vce(robust)
qui margins, dydx(*) predict(n)
matrix CREall = r(b)
matrix CRE = CREall[1..1, 1..6]
qui margins, dydx(*) predict(nu0)
matrix CRE_0all = r(b)
matrix CRE_0 = CRE_0all[1..1, 1..6]

* Report AMEs for PA, RE, FE, and CREs models
matrix rowjoinbyname ME = PA RE REnorm FE_0 CRE CRE_0
matrix rownames ME = PA RE REnorm FE_0 CRE CRE_0
matrix list ME, format(%8.4f)

* Comparison of negative binomial panel estimators
qui xtpoisson mdu lcoins ndisease female age lfam child, pa ///
    corr(exch) vce(robust)
estimates store PPA_ROB
qui xtnbreg mdu lcoins ndisease female age lfam child, pa   ///
    corr(exch) vce(robust)
estimates store NBPA_ROB
qui xtnbreg mdu lcoins ndisease female age lfam child, re 
estimates store NBRE
estimates table PPA_ROB NBPA_ROB NBRE, eq(1) b(%8.4f) se ///
    stats(N ll) stfmt(%8.0f)

********** CHAPTER 22.7 PANEL QUANTILE REGRESSION

* Pooled data quantile regression ignoring individual effects
qui use mus208psid, clear 
qui qreg2 lwage exp exp2 wks ed, quant(.25) cluster(id)
estimates store Pool25
qui qreg2 lwage exp exp2 wks ed, quant(.75) cluster(id)
estimates store Pool75

// Not included - compare standard errors with and without cluster–robust
qreg2 lwage exp exp2 wks ed, quant(.25) cluster(id) 
qreg2 lwage exp exp2 wks ed, quant(.25)
qreg2 lwage exp exp2 wks ed, quant(.75) cluster(id) 
qreg2 lwage exp exp2 wks ed, quant(.75)

// Not included - can do a cluster bootstrap
* xtset id
* bootstrap, reps(400) seed(10101) cluster(id):   ///
*     qreg lwage exp exp2 wks ed, quant(.25)

* FE QR: (1) Filter out the individual FE from the y variable
qui regress lwage exp exp2 wks ed
predict residuals_i, resid
sort id
by id: egen alphahat_i=mean(residuals_i)
summarize alphahat_i
gen lwagehat=lwage-alphahat_i

* FE QR: (2) QR of new y with bootstrap covariance estimates
set seed 10101 
qui bsqreg lwagehat exp exp2 wks ed, quant(.25) reps(400) 
estimates store FE25boot
qui bsqreg lwagehat exp exp2 wks ed, quant(.75) reps(400) 
estimates store FE75boot
xtset id
qui bootstrap, reps(400) seed(10101) cluster(id):  ///
    qreg lwagehat exp exp2 wks ed, quant(.25)
estimates store FE25clu 
qui bootstrap, reps(400) seed(10101) cluster(id):  ///
    qreg lwagehat exp exp2 wks ed, quant(.75)
estimates store FE75clu 

* Compare OLS, FE, Pooled_QR (0.22., 0.75), and Panel_QR(0.25, 0.75)
estimates table Pool25 FE25boot FE25clu Pool75 FE75boot FE75clu, b(%8.4f) se  

********** END
