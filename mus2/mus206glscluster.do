* mus206glscluster.do for Stata 17

capture log close

********** OVERVIEW OF mus206glscluster.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 6
* 6.3: MODELING HETEROSKEDASTIC ERRORS
* 6.4: OLS FOR CLUSTEREED DATA
* 6.5: FGLS ESTIMATORS FOR CLUSTERED DATA
* 6.6: FE ESTIMATOR FOR CLUSTERED DATA
* 6.7: LINEAR MIXED MODELS FOR CLUSTERED DATA
* 6.8: SYSTEMS OF LINEAR REGRESSIONS
* 6.9: SURVEY DATA: WEIGHTING, CLUSTERING, AND STRATIFICATION

* To run you need files
*   mus206vlss.dta
*   mus206mepssur.dta
*   mus206nhanes2.dta (same as http://www.stata-press.com/data/r10/nhanes2.dta)
* in your directory

* Community-contributed commands
*    ivreg2
*    vcemway
* are used

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus206vlss.dta
* Authors' extract from World Bank Vietnam Living Standards survey 1997-98. 
* Data from Cameron and Trivedi (2005, p.848)

* mus206mepssur.dta is authors' extract from MEPS 
* (Medical Expenditure Panel Survey)
* for individuals 65 years and older in U.S. Medicare

* mus206nhanes2.dta  Same as Stata website file psidextract.dta
* Data from second National health and Nutrition Examination Survey (NHANES II)
* U.S. survey 1976-1980.
  
********** 6.3 MODELING HETEROSKEDASTIC ERORS

* This uses generated data
* Model is  y = 1 + 1*x2 + 1*x3 + u
* where     u = sqrt(exp(-1+0.2*x2))*e
*           x1 ~ N(0, 5^2)
*           x2 ~ N(0, 5^2)
*           e ~ N(0, 5^2)
* Errors are conditionally heteroskedastic with V[u|x]=exp(-1+1*x2)

* Generated data for heteroskedasticity example
set seed 10101
qui set obs 500
generate double x2 = 5*rnormal(0)
generate double x3 = 5*rnormal(0)
generate double e  = 5*rnormal(0)
generate double u  = sqrt(exp(-1+0.2*x2))*e
generate double y  = 1 + 1*x2 + 1*x3 + u
summarize

* OLS regression with default standard errors
regress y x2 x3

* OLS regression with heteroskedasticity-robust standard errors
regress y x2 x3, vce(robust) 

* Heteroskedasticity diagnostic scatterplot
qui regress y x2 x3
predict double uhat, resid
generate double absu = abs(uhat)
qui twoway (scatter absu x2) (lowess absu x2, bw(0.4) lw(thick)),       ///
    legend(off) ytitle("Absolute value of residual") yscale(titleg(*5)) ///
    xscale(titleg(*5)) plotr(style(none)) name(gls1, replace)
qui twoway (scatter absu x3) (lowess absu x3, bw(0.4) lw(thick)),       ///
    legend(off) ytitle("Absolute value of residual") yscale(titleg(*5)) ///
    xscale(titleg(*5)) plotr(style(none)) name(gls2, replace)
graph combine gls1 gls2, iscale(1.2) ysize(2.5) xsize(6.0)

* Test heteroskedasticity depending on x2, x3, and x2 and x3
estat hettest x2 x3, mtest

* Separate tests of heteroskedasticity using iid version of hettest
estat hettest x2, iid
estat hettest x3, iid
estat hettest x2 x3, iid

* FGLS: First step get estimate of skedasticity function
drop uhat
qui regress y x2 x3                       // Get bols
predict double uhat, resid
generate double uhatsq = uhat^2           // Get squared residual 
nl (uhatsq = exp({xb: x2}+{b0})), nolog   // NLS of uhatsq on exp(z'a)
predict double varu, yhat                 // Get sigmahat^2

* FGLS: Second step get estimate of skedasticity function
regress y x2 x3 [aweight=1/varu]  

* WLS estimator is FGLS with robust estimate of VCE
regress y x2 x3 [aweight=1/varu], vce(robust)

******* 6.4: OLS FOR CLUSTERED DATA

* Read in Vietnam clustered data, and delete one household in two communes
qui use mus206vlss, clear
drop if lnhhexp > 2.579681 & lnhhexp < 2.579683
drop if missing(lnhhexp)

* Create a unique household identifier and summarize data
egen hh = group(lnhhexp)
summarize pharvis lnhhexp illness hh commune

* OLS estimation with cluster–robust standard errors clustering on household
regress pharvis lnhhexp illness, vce(cluster hh)

* OLS estimation with various cluster–robust standard errors
qui regress pharvis lnhhexp illness
estimates store OLS_iid
qui regress pharvis lnhhexp illness, vce(robust)
estimates store OLS_het
qui regress pharvis lnhhexp illness, vce(cluster hh)
estimates store OLS_hh
qui regress pharvis lnhhexp illness, vce(bootstrap, cluster(hh) seed(10101) reps(400))
estimates store OLS_boot
qui regress pharvis lnhhexp illness, vce(cluster commune) 
estimates store OLS_comm
estimates table OLS_iid OLS_het OLS_hh OLS_boot OLS_comm,  ///
    b(%10.4f) se stats(r2 N)

* Within commune intracluster correlation for lnhhexp
loneway lnhhexp commune

* Standard-error inflation factor for b_lnhhexp when clustering on commune
qui loneway lnhhexp commune
scalar rho_x = r(rho)
qui regress pharvis lnhhexp illness
scalar numobs = e(N)
predict uhat, resid 
qui loneway uhat commune
scalar rho_uhat = r(rho)
qui sum 
scalar M = numobs/194
display "rho_x = " rho_x "  rho_u = " rho_uhat "  M = " M  ///
    _n "Standard error inflation = " sqrt(1+rho_x*rho_uhat*(M-1))

// Not included - alternative standard error inflation factor 
* Alternative inflation factor based on rho(x*u) and not on rho(x) times rho(u)
generate lnhhexpbyuhat = lnhhexp*uhat
loneway lnhhexpbyuhat commune   
scalar rho_xbyuhat = r(rho)
display "rho_xbyuhat = " rho_xbyuhat "  M = " m  ///
    _n "Standard error inflation = " sqrt(1+rho_xbyuhat*(m-1))
   

* Two-way cluster–robust standard errors using vcemway command
vcemway regress pharvis lnhhexp illness, cluster(commune illdays)


* Effect of two-way versus separate one-way clustering
qui regress pharvis lnhhexp illness, vce(robust)
estimates store het
qui regress pharvis lnhhexp illness, cluster(commune)
estimates store commune
qui regress pharvis lnhhexp illness, cluster(illdays)
estimates store illdays
qui vcemway regress pharvis lnhhexp illness, cluster(commune illdays)
estimates store twoway
qui vcemway regress pharvis lnhhexp illness, cluster(commune illdays) ///
    vmcfactor(minimum)
estimates store twowayalt
estimates table het commune illdays twoway twowayalt, b(%10.4f) se stats(r2 N)



******* 6.5: FGLS ESTIMATORS FOR CLUSTERED DATA

* Define the cluster identifier variable
xtset hh

* FGLS with equicorrelated errors using using xtreg, pa
xtreg pharvis lnhhexp illness, pa corr(exchangeable) vce(robust)

* Random hh effects FGLS using xtreg, re
xtreg pharvis lnhhexp illness, re vce(robust)

* Random hh effects FGLS with various standard errors
qui xtreg pharvis lnhhexp illness, re mle 
estimates store RE_MLE
qui xtreg pharvis lnhhexp illness, re
estimates store RE_def
qui xtreg pharvis lnhhexp illness, re vce(robust)
estimates store RE_rob
qui xtreg pharvis lnhhexp illness, re vce(cluster hh)
estimates store RE_hh
qui xtreg pharvis lnhhexp illness, re vce(cluster commune)
estimates store RE_comm
estimates table RE_MLE RE_def RE_rob RE_hh RE_comm RE_hh, ///
    keep(lnhhexp illness _cons) b(%10.4f) se stats(N) eq(1)

******* 6.6: FIXED EFFECTS ESTIMATOR FOR CLUSTERED DATA

* FE estimation with hh FE using xtreg, fe
xtreg pharvis lnhhexp illness, fe vce(robust)

* FE estimation with hh RE and various standard errors
qui xtreg pharvis lnhhexp illness, fe
estimate store FE_def
qui xtreg pharvis lnhhexp illness, fe vce(robust)
estimate store FE_rob
qui xtreg pharvis lnhhexp illness, fe vce(cluster hh)
estimate store FE_hh
qui xtreg pharvis lnhhexp illness, fe vce(cluster commune)
estimate store FE_comm
qui regress pharvis lnhhexp illness, vce(cluster hh)
estimates store OLS_hh
estimates table FE_def FE_rob FE_hh FE_comm OLS_hh, b(%10.4f) se stats(N)

* Generate integer-valued person identifiers
sort hh
by hh: generate numpersonsinhh = _N

* Various standard errors for hh FE estimation using xtreg, reg, and areg
preserve
keep if numpersonsinhh == 2
qui regress pharvis lnhhexp illness i.hh
estimates store LSDV_def
qui areg pharvis lnhhexp illness, absorb(hh)
estimates store AREG_def
qui xtreg pharvis lnhhexp illness, fe
estimates store XTREG_def
qui regress pharvis lnhhexp illness i.hh, vce(cluster hh)
estimates store LSDV_clu
qui areg pharvis lnhhexp illness, absorb(hh) vce(cluster hh)
display "AREG se times square-root(2) = " _se[illness]/sqrt(2)
estimates store AREG_clu
qui xtreg pharvis lnhhexp illness, fe vce(cluster hh)
estimates store XTREG_clu
estimates table LSDV_def AREG_def XTREG_def LSDV_clu AREG_clu XTREG_clu, ///
    keep(illness lnhhexp  _cons) b(%8.4f) se stats(N r2 df_m) varwidth(8)
restore

* RE estimation with hh fixed
sort hh
by hh: egen avelnhhexp = mean(lnhhexp)
by hh: egen aveillness = mean(illness)
xtreg pharvis avelnhhexp aveillness lnhhexp illness, re vce(cluster hh) 

// Not included
* Generate integer-valued person identifiers
sort hh
by hh: generate person = _n
sum person
* Define the cluster identifier variable and the within cluster identifier
xtset hh person
* Illustrate use
xtdescribe

******* 6.7: LINEAR MIXED MODELS FOR CLUSTERED DATA

* Random hh intercept model estimated using mixed
mixed pharvis lnhhexp illness || hh: , mle vce(robust)

// Not included as reml does not allow robust standard errors
mixed pharvis lnhhexp illness || hh: , reml 

* Robust standard errors after mixed (and equivalence with xtreg, mle)
qui xtreg pharvis lnhhexp illness, mle
estimates store remle_def
qui mixed pharvis lnhhexp illness || hh: , mle 
estimates store mixed_def
qui mixed pharvis lnhhexp illness || hh: , mle vce(robust)
estimates store mixed_rob
qui mixed pharvis lnhhexp illness || hh: , mle vce(cluster hh)
estimates store mixed_hh
qui mixed pharvis lnhhexp illness || hh: , mle vce(cluster commune)
estimates store mixed_comm
estimates table remle_def mixed_def mixed_rob mixed_hh mixed_comm, ///
    b(%10.4f) se stats(N)

* Random-slopes model estimated using mixed
mixed pharvis lnhhexp illness || hh: illness, mle covar(unstructured) ///
    difficult noheader nolog 

// Not included - default standard errors 
mixed pharvis lnhhexp illness || hh: illness, mle covar(unstructured) ///
    difficult noheader nolog 
  
* Hierarchical linear model with household and commune variance components
mixed pharvis lnhhexp illness || commune: || hh:, mle difficult ///
    nolog vce(cluster commune)

* Two-way nonnested errors - illness and commune
mixed pharvis lnhhexp illness || _all: R.illness || commune: , ///
    mle difficult nolog

// Not included - same results but slower
* mixed pharvis lnhhexp illness || _all: R.commune || illness: , mle difficult 

********** * 6.8: SYSTEMS OF LINEAR EQUATIONS

* Summary statistics for seemingly unrelated regressions example
clear all
qui use mus206mepssur
summarize ldrugexp ltotothr age age2 educyr actlim totchr medicaid private

summarize ldrugexp if ldrugexp!=. & ltotothr!=.

* SUR estimation of a seemingly unrelated regressions model
sureg (ldrugexp age age2 actlim totchr medicaid private) ///
    (ltotothr age age2 educyr actlim totchr private), corr

* Bootstrap to get heteroskedasticity-robust standard errors for SUR estimator 
bootstrap, reps(400) seed(10101) nodots: sureg          ///
    (ldrugexp age age2 actlim totchr medicaid private)  ///
    (ltotothr age age2 educyr actlim totchr private) 

* Test of variables in both equations
qui sureg (ldrugexp age age2 actlim totchr medicaid private) ///
    (ltotothr age age2 educyr actlim totchr private)
test age age2

* Test of variables in just the first equation
test [ldrugexp]age [ldrugexp]age2

* Test of a restriction across the two equations
test [ldrugexp]private = [ltotothr]private

* Specify a restriction across the two equations
constraint 1 [ldrugexp]private = [ltotothr]private

//  Out short version
// sureg (ldrugexp actlim totchr medicaid private)  ///
//   (ltotothr educyr actlim totchr private), constraints(1) noheader

* Estimate subject to the cross-equation constraint
sureg (ldrugexp age age2 actlim totchr medicaid private)        ///
    (ltotothr age age2 educyr actlim totchr private), constraints(1) 

* suest used for cross-equations test of separately estimated models
qui regress ldrugexp age age2 actlim totchr medicaid private
estimates store DRUG
qui regress ltotothr age age2 educyr actlim totchr private
estimates store OTHER
suest DRUG OTHER
test _b[DRUG_mean:private]  = _b[OTHER_mean:private]

// Not included - example that gives cluster–robust standard errors 
suest, coeflegend
suest DRUG OTHER, vce(cluster age)

******** 6.9 SURVEY DATA: WEIGHTING, CLUSTERING AND STRATIFICATION

* Survey data example: NHANES II data
clear all
qui use mus206nhanes
qui keep if age >= 21 & age <= 65
describe sampl finalwgt strata psu
summarize sampl finalwgt strata psu

* Declare survey design
svyset psu [pweight=finalwgt], strata(strata)

* Describe the survey design
svydescribe

* Estimate the population mean using svy:
svy: mean hgb

* Estimate the population mean using no weights and no cluster
mean hgb

* Regression using svy:
svy: regress hgb age female

* Regression using weights and cluster on PSU
generate uniqpsu = 2*strata + psu  // Make unique identifier for each PSU
regress hgb age female [pweight=finalwgt], vce(cluster uniqpsu)

* Regression using no weights and no cluster
regress hgb age female

********** END
