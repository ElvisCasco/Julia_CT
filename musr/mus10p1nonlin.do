* mus10p1nonlin.do  Oct 2009 for Stata version 11

cap log close

********** OVERVIEW OF mus10p1nonlin.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Chapter 10
* 10.2: NONLINEAR EXAMPLE: DOCTOR VISITS
* 10.3: NONLINEAR REGRESSION METHODS
* 10.4: DIFFERENT ESIMATES OF THE VCE
* 10.5: PREDICTION
* 10.6: MARGINAL EFFECTS
* 10.7: MODEL DIAGNOSTICS

* To run you need files
*   mus10data.dta    
* in your directory

* To speed up program reduce reps in vce(boot,reps(400))

********** SETUP **********

set more off
version 11
clear all
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* 2002 Medical Expenditure Panel Survey (MEPS)
* U.S. individuals aged 25-64 years working in private sector 
* but not self-employed
* and not receiving public insurance (Medicare and Medicaid) 
* Data due to Deb, Munkin and Trivedi (2006)

******* 10.2: NONLINEAR EXAMPLE: DOCTOR VISITS

* Read in dataset, select one year of data, and describe key variables
use mus10data.dta, clear
keep if year02==1
describe docvis private chronic female income

* Summary of key variables
summarize docvis private chronic female income

******* 10.3: NONLINEAR REGRESSION METHODS

* Poisson regression (command poisson)
poisson docvis private chronic female income, vce(robust)

nl (docvis = exp({private}*private + {chronic}*chronic     ///
            + {female}*female + {income}*income + {intercept}))

* Nonlinear least-squares regression (command nl)
generate one = 1
nl (docvis = exp({xb: private chronic female income one})), vce(robust) nolog

* Generalized linear models regression for poisson (command glm)
glm docvis private chronic female income,  ///
  family(poisson) link(log) vce(robust) nolog

* Command gmm for GMM estimation (nonlinear IV) for Poisson model
gmm (docvis - exp({xb:private chronic female income}+{b0})),   ///
 instruments(firmsize chronic female income) onestep nolog

******* 10.4: DIFFERENT ESTIMATES OF THE VCE

* Different VCE estimates after Poisson regression
quietly poisson docvis private chronic female income
estimates store VCE_oim
quietly poisson docvis private chronic female income, vce(opg)
estimates store VCE_opg
quietly poisson docvis private chronic female income, vce(robust)
estimates store VCE_rob
quietly poisson docvis private chronic female income, vce(cluster age)
estimates store VCE_clu
set seed 10101
quietly poisson docvis private chronic female income, vce(boot,reps(400))
estimates store VCE_boot
estimates table VCE_oim VCE_opg VCE_rob VCE_clu VCE_boot, b(%8.4f) se

******* 10.5: PREDICTION

* Predicted mean number of doctor visits using predict and predictnl 
quietly poisson docvis private chronic female income, vce(robust)
predict muhat if e(sample), n
predictnl muhat2 = exp(_b[private]*private + _b[chronic]*chronic  ///
   + _b[female]*female + _b[income]*income + _b[_cons]), se(semuhat2)
summarize docvis muhat muhat2 semuhat2

predict resid, score
generate resid2 = docvis - muhat
summarize resid resid2 

* Out-of-sample prediction for year01 data using year02 estimates
use mus10data.dta, clear
quietly poisson docvis private chronic female income if year02==1, vce(robust)
keep if year01 == 1
predict muhatyear01, n    
summarize docvis muhatyear01

* Prediction at a particular value of one of the regressors
use mus10data.dta, clear
keep if year02 == 1
quietly poisson docvis private chronic female income, vce(robust)
preserve
replace private = 1
predict muhatpeq1, n
summarize muhatpeq1
restore

* Prediction at a specified value of all the regressors
nlcom exp(_b[_cons]+_b[private]*1+_b[chronic]*0+_b[female]*1+_b[income]*10)

* Predict at a specified value of all the regressors
lincom _cons + private*1 + chronic*0 + female*1 + income*10, eform

* Sample average of predicted number of events
quietly poisson docvis private chronic female income, vce(robust)
margins

* Sample average prediction at a particular value of one of the regressors
margins, at(private=1)

* Prediction at a specified value of all regressors
margins, at(private=1 chronic=0 female=1 income=10)

* Sample average prediction at different values of indicator variable private
quietly poisson docvis i.private chronic female income, vce(robust)
margins private

******* 10.6: MARGINAL EFFECTS

* Marginal effects at mean (MEM) using margins command and finite differences
quietly poisson docvis i.private i.chronic i.female income, vce(robust)
margins, dydx(*) atmean

* Marginal effects at mean (MEM) using margins command and only calculus method
quietly poisson docvis private chronic female income, vce(robust)
margins, dydx(*) atmean

* Marginal effects at representative value (MER) using margins command
quietly poisson docvis i.private i.chronic i.female income, vce(robust)
margins, dydx(*) at(private=1 chronic=0 female=1 income=10) noatlegend
* Average marginal effects (AVE) using margins command and finite differences
quietly poisson docvis i.private i.chronic i.female income, vce(robust)
margins, dydx(*)

* Usual ME evaluated at mean of regressors
quietly poisson docvis private chronic female income, vce(robust)
margins, dydx(income) atmean noatlegend

* Elasticity evaluated at mean of regressors
margins, eyex(income) atmean noatlegend

* Semielasticity evaluated at mean of regressors
margins, eydx(income) atmean noatlegend

* Other semielasticity evaluated at mean of regressors
margins, dyex(income) atmean noatlegend

* AME computed manually for a single regressor
use mus10data.dta, clear
keep if year02 == 1
quietly poisson docvis private chronic female income, vce(robust)
preserve
predict mu0, n
quietly replace income = income + 0.01
predict mu1, n
generate memanual = (mu1-mu0)/0.01
summarize memanual
restore

* AME computed manually for all regressors
global xlist private chronic female income
preserve
predict mu0, n
foreach var of varlist $xlist {
  quietly summarize `var'
  generate delta = r(sd)/1000
  quietly generate orig = `var'
  quietly replace `var' = `var' + delta
  predict mu1, n
  quietly generate me_`var' = (mu1 - mu0)/delta
  quietly replace `var' = orig
  drop mu1 delta orig
 }
summarize me_*
restore  

* AME for a polynomial regressor: manual computation
generate inc2 = income^2
generate inc3 = income^3
quietly poisson docvis private chronic female income inc2 inc3, vce(robust)
predict muhat, n
generate me_income = muhat*(_b[income]+2*_b[inc2]*income+3*_b[inc3]*inc2)
summarize me_income

* AME for a polynomial regressor: computation using factor variables
quietly poisson docvis private chronic female c.income c.income#c.income c.income#c.income#c.income, vce(robust)
margins, dydx(income)

* Specify model with interacted regressors using factor variables
poisson docvis private chronic i.female c.income i.female#c.income, vce(robust) nolog

* AME with interacted regressors given model specified using factor variables
margins, dydx(female income)

* AME computed manually for a complex model
preserve
predict mu0, n
quietly summarize income
generate delta = r(sd)/100
quietly replace income = income + delta
quietly replace inc2 = income^2
quietly replace inc3 = income^3
predict mu1, n 
generate me_inc = (mu1 - mu0)/delta
summarize me_inc
restore

******* 10.7: MODEL DIAGNOSTICS

* Compute pseudo-R-squared after Poisson regression
quietly poisson docvis private chronic female income, vce(robust)
display "Pseudo-R^2 = " 1 - e(ll)/e(ll_0)

* Report information criteria
estat ic

* Various residuals after command glm
quietly glm docvis private chronic female income, family(poisson)
predict mu, mu
generate uraw = docvis - mu
predict upearson, pearson
predict udeviance, deviance
predict uanscombe, anscombe
summarize uraw upearson udeviance uanscombe

// Following not presented in book
correlate uraw upearson udeviance uanscombe

********** CLOSE OUTPUT
