* mus213nonlin.do   for Stata 17

capture log close

********** OVERVIEW OF mus213nonlin.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 13
* 13.2: NONLINEAR EXAMPLE: DOCTOR VISITS
* 13.3: NONLINEAR REGRESSION METHODS
* 13.4: DIFFERENT ESTIMATES OF THE VCE
* 13.5: PREDICTION
* 13.6: PREDICTIVE MARGINS
* 13.7: MES
* 13.8: MODEL DIAGNOSTICS
* 13.9: CLUSTERED DATA

* To run you need files
*   mus210mepsdocvisyoung.dta 
*   mus206vlss.dta   
* in your directory

* No community-contributed proagrams are used 

* To speed up program reduce reps in vce(bootstrap, reps(400))

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus210mepsdocvisyoung.dta authors' extract from 2002 MEPS
* Medical Expenditure Panel Survey
* U.S. data on office-based physician visits by persons aged 25-64 years
* Same as Deb, Munkin and Trivedi, "Bayesian analysis of the two-part model
* with endogeneity: application to health care expenditure", JBES, 1081-1099  
* Excludes those receiving public insurance (Medicare and Medicaid)
* Restricted to those working in the private sector but not self-employed.

* mus206vlss.dta
* Authors' extract from World Bank Vietnam Living Standards Survey 1997-98. 
* Data from Cameron and Trivedi (2005, p.848)

******* 13.2: NONLINEAR EXAMPLE: DOCTOR VISITS

* Read in dataset, select one year of data, and describe key variables
qui use mus210mepsdocvisyoung
keep if year02==1
describe docvis private chronic female income

* Summary of key variables
summarize docvis private chronic female income

******* 13.3: NONLINEAR REGRESSION METHODS

* Poisson regression (command poisson)
poisson docvis private chronic female income, vce(robust)

* ML estimation (command mlexp) for Poisson model
mlexp (-exp({xb: private chronic female income _cons}) + docvis*{xb:} - ///
       lnfactorial(docvis)), vce(robust)

* Nonlinear least-squares regression (command nl) with default standard errors
nl (docvis = exp({private}*private + {chronic}*chronic                  ///
      + {female}*female + {income}*income + {intercept})) 

* Nonlinear least-squares regression - alternative form of commannd nl
nl (docvis = exp({xb: private chronic female income} + {intercept})),  ///
    vce(robust) nolog

* Generalized linear models regression for poisson (command glm)
glm docvis private chronic female income,  ///
    family(poisson) link(log) vce(robust) nolog

* Command gmm for GMM estimation (nonlinear IV) for Poisson model
gmm (docvis - exp({xb:private chronic female income _cons})), ///
    instruments(firmsize chronic female income) onestep nolog 

// Output not included
* Computation of AMEs following for this gmm example
margins, dydx(*) expression(exp(predict(xb)))

* Data and globals for gmm estimation of chapter 7 control function example 
qui use mus207mepspresdrugs, clear
global y1 ldrugexp             // Dependent variable
global y2 hi_empunion          // Endogenous regressor
global x2list totchr age female blhisp linc  // Exogenous regressors
global xlist ssiratio $x2list  // Structural equation regressors
global zlist $y2 $x2list       // First-stage regressors 

* Command gmm for estimation of control function linear model
gmm (eq1: ($y2 - {zpi:ssiratio $x2list _cons} ) )                      ///
    (eq2: ($y1 - {xb:$zlist _cons} - {gamma}*($y2 - {zpi:}) ) )        /// 
    (eq3: (($y2 -{zpi:}) * ($y1 - {xb:} -{gamma}*($y2 - {zpi:})) ) ) , ///
    instruments(eq1: $xlist) instruments(eq2: $zlist)                  ///  
    instruments(eq3: )	onestep winitial(unadjusted, independent) nolog

******* 13.4: DIFFERENT ESTIMATES OF THE VCE

* Different VCE estimates after Poisson regression
qui use mus210mepsdocvisyoung, clear
keep if year02==1
qui poisson docvis private chronic female income
estimates store VCE_oim
qui poisson docvis private chronic female income, vce(opg)
estimates store VCE_opg
qui poisson docvis private chronic female income, vce(robust)
estimates store VCE_rob
qui poisson docvis private chronic female income, vce(cluster age)
estimates store VCE_clu
set seed 10101
qui poisson docvis private chronic female income, vce(bootstrap, reps(400))
estimates store VCE_boot
estimates table VCE_oim VCE_opg VCE_rob VCE_clu VCE_boot, b(%8.4f) se

******* 13.5: PREDICTION

* Predicted mean number of doctor visits using predict and predictnl 
qui poisson docvis private chronic female income, vce(robust)
predict muhat if e(sample), n
predictnl muhat2 = exp(_b[private]*private + _b[chronic]*chronic  ///
    + _b[female]*female + _b[income]*income + _b[_cons]), se(semuhat2)
summarize docvis muhat muhat2 semuhat2

predict resid, score
generate resid2 = docvis - muhat
summarize resid resid2 

* Out-of-sample prediction for year01 data using year02 estimates
qui use mus210mepsdocvisyoung, clear
qui poisson docvis private chronic female income if year02==1, vce(robust)
keep if year01 == 1
predict muhatyear01, n    
summarize docvis muhatyear01

* Prediction at a particular value of one of the regressors
qui use mus210mepsdocvisyoung, clear
keep if year02 == 1
qui poisson docvis private chronic female income, vce(robust)
preserve
replace private = 1
predict muhatpeq1, n
summarize muhatpeq1
restore

* Predict at a specified value of all the regressors using nlcom
nlcom exp(_b[_cons]+_b[private]*1+_b[chronic]*0+_b[female]*1+_b[income]*10)

* Predict at a specified value of all the regressors using lincom
lincom _cons + private*1 + chronic*0 + female*1 + income*10, eform

******* 13.7: PREDICTIVE MARGINS

* margins: Sample average of predicted number of events
qui poisson docvis private chronic female income, vce(robust)
margins

// Not included - use expression
margins, expression(exp(predict(xb)))

* margins: Sample average prediction at a particular value of a regressor
margins, at(private=1)

* margins: Prediction at a specified value of all regressors
margins, at(private=1 chronic=0 female=1 income=10)

* margins: Sample avg. prediction at different values of an indicator variable
qui poisson docvis i.private chronic female income, vce(robust)
margins private

******* 13.7: MARGINAL EFFECTS

* AMEs using margins command and finite differences
qui poisson docvis i.private i.chronic i.female income, vce(robust)
margins, dydx(*)

* AMEs using margins command and only calculus method
qui poisson docvis private chronic female income, vce(robust)
margins, dydx(*) 

* MEMs using margins command and finite differences
qui poisson docvis i.private i.chronic i.female income, vce(robust)
margins, dydx(*) atmeans

* MERs using margins command
qui poisson docvis i.private i.chronic i.female income, vce(robust)
margins, dydx(*) at(private=1 chronic=0 female=1 income=10) noatlegend

// Not included - same AME with vce(unconditional)
qui poisson docvis i.private i.chronic i.female income, vce(robust)
margins, dydx(*) vce(unconditional)

// Not included - example using expression()
margins, dydx(*) expression(exp(predict(xb)))

* AME computed manually for a single regressor
qui use mus210mepsdocvisyoung, clear
keep if year02 == 1
qui poisson docvis private chronic female income, vce(robust)
preserve
predict mu0, n
qui replace income = income + 0.01
predict mu1, n
generate memanual = (mu1-mu0)/0.01
summarize memanual
restore

* AME computed manually for all regressors
global xlist private chronic female income
preserve
predict mu0, n
foreach var of varlist $xlist {
    qui summarize `var'
    generate delta = r(sd)/1000
    qui generate orig = `var'
    qui replace `var' = `var' + delta
    predict mu1, n
    qui generate me_`var' = (mu1 - mu0)/delta
    qui replace `var' = orig
    drop mu1 delta orig
 }
summarize me_*
restore  

* AME for a polynomial regressor: Manual computation
generate inc2 = income^2
generate inc3 = income^3
qui poisson docvis private chronic female income inc2 inc3, vce(robust)
predict muhat, n
generate me_income = muhat*(_b[income]+2*_b[inc2]*income+3*_b[inc3]*inc2)
summarize me_income

* AME for a polynomial regressor: Computation using factor variables
qui poisson docvis private chronic female c.income c.income#c.income ///
   c.income#c.income#c.income, vce(robust)
margins, dydx(income)

* Specify model with interacted regressors using factor variables
poisson docvis private chronic i.female c.income i.female#c.income, ///
    vce(robust) nolog

* AME with interacted regressors given model specified using factor variables
margins, dydx(female income)

* AME computed manually for a complex model
preserve
predict mu0, n
qui summarize income
generate delta = r(sd)/100
qui replace income = income + delta
qui replace inc2 = income^2
qui replace inc3 = income^3
predict mu1, n 
generate me_inc = (mu1 - mu0)/delta
summarize me_inc
restore

* Usual ME evaluated at means of regressors
qui poisson docvis private chronic female income, vce(robust)
margins, dydx(income) atmeans noatlegend

* Elasticity evaluated at means of regressors
margins, eyex(income) atmeans noatlegend

* Semielasticity evaluated at means of regressors
margins, eydx(income) atmeans noatlegend

* Other semielasticity evaluated at means of regressors
margins, dyex(income) atmeans noatlegend

******* 13.8: MODEL DIAGNOSTICS

* Compute pseudo-R-squared after Poisson regression
qui poisson docvis private chronic female income, vce(robust)
display "Pseudo-R^2 = " 1 - e(ll)/e(ll_0)

* Report information criteria
estat ic

* Various residuals after command glm
qui glm docvis private chronic female income, family(poisson) vce(robust)
predict mu, mu
generate uraw = docvis - mu
predict upearson, pearson
predict udeviance, deviance
predict uanscombe, anscombe
summarize uraw upearson udeviance uanscombe

correlate uraw upearson udeviance uanscombe

******* 13.9: CLUSTERED DATA

* Read in Vietnam clustered data and delete one household in two communes
qui use mus206vlss, clear
drop if lnhhexp > 2.579681 & lnhhexp < 2.579683
drop if missing(lnhhexp) 
summarize pharvis lnhhexp illness commune

// Not included 
* Poisson estimation with default standard errors
poisson pharvis lnhhexp illness

// Not included - Poisson with heteroskedastic-robust standard errors
poisson pharvis lnhhexp illness, nolog vce(robust)

* Poisson estimation with cluster–robust standard errors clustered on commune
poisson pharvis lnhhexp illness, nolog vce(cluster commune)

nbreg pharvis lnhhexp illness, nolog
nbreg pharvis lnhhexp illness, nolog vce(robust)
nbreg pharvis lnhhexp illness, nolog vce(cluster commune)

* Generalized estimating equation estimation with cluster–robust st. errors clustered on commune
xtset commune
xtgee pharvis lnhhexp illness, nolog family(poisson) ///
  corr(exchangeable) vce(robust)

* RE estimation with cluster–robust st. errors clustered on commune
xtpoisson pharvis lnhhexp illness, re normal nolog vce(robust)

* Mixed estimation with cluster–robust st. errors clustered on commune
mepoisson pharvis lnhhexp illness || commune: , nolog vce(robust)

* FE estimation with cluster–robust st. errors clustering on commune
xtpoisson pharvis lnhhexp illness, fe nolog vce(robust)

* Dummy variables estimation with robust st. errors clustered on commune
qui poisson pharvis lnhhexp illness i.commune, vce(cluster commune)
estimates store POISSONDV
estimates table POISSONDV, keep(lnhhexp illness) b(%10.6f) se stats(N)

* Correlated RE estimation with cluster-specific effects
bysort commune: egen avelnhhexp = mean(lnhhexp)
by commune: egen aveillness = mean(illness)
xtpoisson pharvis lnhhexp illness avelnhhexp aveillness, re nolog vce(robust)

********** END
