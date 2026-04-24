* mus212bootstrap.do  for Stata 17

capture log close

********** OVERVIEW OF mus212bootstrap.do **********

* Stata program
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi
* used for "Microeconometrics Using Stata, Second Edition"
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 12
* 12.3: BOOTSTRAP PAIRS USING THE VCE(BOOTSTRAP) OPTION
* 12.4: BOOTSTRAP PAIRS USING THE BOOTSTRAP COMMAND
* 12.5: BOOTSTRAPS WITH ASYMPTOTIC REFINEMENT
* 12.6: WILD BOOTSTRAP WITH ASYMPTOTIC REFINEMENT
* 12.7: BOOTSTRAP PAIRS USING BSAMPLE AND SIMULATE
* 12.8: ALTERNATIVE RESAMPLING SCHEMES
* 12.9: THE JACKKNIFE

* To run you need files
*   mus212bootdata.dta
*   mus207mepspresdrugs.dta
*   mus212occfewcluster.dta
*   mus219mepsambexp.dta
* in your directory

* community-contributed command
*     boottest
* is used

* To speed up program reduce reps in
*     vce(bootstrap, reps(400)) and vce(bootstrap, reps(999))
* or in bootstrap reps(999)

********** SETUP **********

clear all
set linesize 81
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus212bootdata.dta an extract of 50 observations from mus210mepsdocvisyoung.dta

* mus212occfewclusters is extract from Cameron, Gelbach and Miller (2011),
* "Robust Inference with Multi-way Clustering," JBES, 238-249.
* Based in turn on Hersch (1998), "Compensating Wage Differentials for "
* Gender-Specific Job Injury Rates", AER, 598-607.

* Other datsets used are
*   mus207mepspresdrugs.dta
*   mus219mepsambexp.dta


****** 12.3: BOOTSTRAP PAIRS USING THE CE(BOOT) OPTION

/* How mus212bootdata.dta was created
use mus210mepsdocvisyoung, clear
qui keep if year02 == 1
qui drop if _n > 50
qui keep docvis chronic age
qui save mus212bootdata.dta, replace
*/

* Sample is 50 observations and a few variables from chapter 10 data
qui use mus212bootdata
summarize

* Option vce(bootstrap) to compute bootstrap standard errors
poisson docvis chronic, vce(bootstrap, reps(400) seed(10101) nodots)

* Bootstrap standard errors for different reps and seeds
qui poisson docvis chronic, vce(bootstrap, reps(50) seed(10101))
estimates store boot50
qui poisson docvis chronic, vce(bootstrap, reps(50) seed(20202))
estimates store boot50diff
qui poisson docvis chronic, vce(bootstrap, reps(2000) seed(10101))
estimates store b2000
qui poisson docvis chronic, vce(bootstrap, reps(2000) seed(20202))
estimates store b2000diff
qui poisson docvis chronic, vce(robust)
estimates store robust
estimates table boot50 boot50diff b2000 b2000diff robust, b(%8.5f) se(%8.5f)

* Option vce(bootstrap, cluster) to compute cluster-bootstrap standard errors
poisson docvis chronic, vce(bootstrap, cluster(age) reps(400) seed(10101) nodots)

poisson docvis chronic, vce(cluster age)

* Bootstrap confidence intervals: Normal based, percentile, BC, and BCa
qui poisson docvis chronic, vce(bootstrap, reps(999) seed(10101) bca)
estat bootstrap, all

matrix list e(b_bs)

****** 12.4: BOOTSTRAP PAIRS USING THE BOOTSTRAP COMMAND

* bootstrap command applied to Stata estimation command
bootstrap, reps(400) seed(10101) nodots noheader: poisson docvis chronic

* Bootstrap estimate of the standard error of a coefficient estimate
bootstrap _b _se, reps(400) seed(10101) nodots: poisson docvis chronic

* Program to return b and robust estimate V of the VCE
program poissrobust, eclass
    version 17
    tempname b V
    poisson docvis chronic, vce(robust)
    matrix `b' = e(b)
    matrix `V' = e(V)
    ereturn post `b' `V'
end

* Check preceding program by running once
poissrobust
ereturn display

* Bootstrap standard-error estimate of robust standard errors
bootstrap _b _se, reps(400) seed(10101) nodots nowarn: poissrobust

* Set up the selection model two-step estimator data of the tobit chapter
qui use mus219mepsambexp, clear
global xlist age female educ blhisp totchr ins

* Program to return b for Heckman two-step estimator of selection model
program hecktwostep, eclass
    version 17
    tempname b V
    tempvar xb
    capture drop invmills
    probit dy $xlist
    predict `xb', xb
    generate invmills = normalden(`xb')/normprob(`xb')
    regress lny $xlist invmills if dy==1
    matrix `b' = e(b)
    ereturn post `b'
end

// Not included - Check preceding program by running once
hecktwostep
ereturn display

* Bootstrap for Heckman two-step estimator using tobit chapter example
bootstrap _b, reps(400) seed(10101) nodots nowarn: hecktwostep

// Not included - Check results - same coefficients
heckman lny $xlist, select(dy = $xlist) twostep

* Program to return (b1-b2) for Hausman test of endogeneity
program hausmantest, eclass
    version 17
    tempname b bols biv
    regress ldrugexp hi_empunion totchr age female blhisp linc
    matrix `bols' = e(b)
    ivregress 2sls ldrugexp (hi_empunion = ssiratio) totchr age female blhisp linc, vce(robust)
    matrix `biv' = e(b)
    matrix `b' = `bols' - `biv'
    ereturn post `b'
end

* Check preceding program by running once
qui use mus207mepspresdrugs, clear
hausmantest
ereturn display

* Bootstrap estimates for Hausman test using IV chapter example
qui use mus207mepspresdrugs, clear
bootstrap _b, reps(400) seed(10101) nodots nowarn: hausmantest

* Perform Hausman test on coefficient of the potentially endogenous regressor
test hi_empunion

* Perform Hausman test on the coefficients of all regressors
test hi_empunion totchr age female blhisp linc _cons

* Bootstrap estimate of the standard error of the coefficient of variation
qui use mus212bootdata, clear
bootstrap coeffvar=(r(sd)/r(mean)), reps(400) seed(10101) nodots   ///
    nowarn: summarize docvis

****** 12.5: BOOTSTRAPS WITH ASYMPTOTIC REFINEMENT

* Percentile-t for a single coefficient: Bootstrap the t statistic
qui use mus212bootdata, clear
qui poisson docvis chronic, vce(robust)
local theta = _b[chronic]
local setheta = _se[chronic]
bootstrap tstar=((_b[chronic]-`theta')/_se[chronic]), seed(10101) nodots ///
    reps(999) saving(percentilet, replace): poisson docvis chronic,  ///
    vce(robust)

* percentile-t: Plot the density of tstar
use percentilet, clear
tabstat tstar, stats(count mean sd skew kurt)	
kdensity tstar, bw(0.2) normal legend(off) xtitle("tstar") ///
    note(" ") scale(1.2) plotregion(style(none)) title(" ")

* Percentile-t p-value for symmetric two-sided Wald test of H0: theta = 0
qui count if abs(`theta'/`setheta') < abs(tstar)
display "p-value = " r(N)/_N

* Percentile-t confidence interval
_pctile tstar, p(2.5,97.5)
scalar lb = `theta' + r(r1)*`setheta'
scalar ub = `theta' + r(r2)*`setheta'
display "2.5 and 97.5 percentiles of t* distn: " r(r1) ", " r(r2) _n ///
    "95 percent percentile-t confidence interval is  (" lb ","  ub ")"

* Not included - compare to negative binomial default standard errors
qui use mus212bootdata, clear
nbreg docvis chronic

****** 12.6: WILD BOOTSTRAP WITH ASYMPTOTIC REFINEMENT

* Few clusters OLS: Heteroskedastic-robust standard errors
qui use mus212occfewcluster, clear
global xlist potexp potexpsq educ union nonwhite northe midw west
regress lnw occrate $xlist, vce(robust)

* Few clusters OLS: Cluster-robust standard errors
regress lnw occrate $xlist, vce(cluster occ_id)

* Few clusters OLS: CSS conservative effective number of clusters
clusteff lnw occrate $covars, cluster(occ_id) test(occrate)

* Few clusters OLS: Wild cluster bootstrap with Webb weights
boottest occrate, seed(10101) reps(9999) weight(webb)

// Not included - default Rademacher weights
boottest occrate, seed(10101)

* Few clusters logit: Cluster-robust standard errors
drop if occ_id==63 |  occ_id==113 | occ_id==133
generate dlnw = lnw > ln(10)
logit dlnw occrate $xlist, nolog vce(cluster occ_id)

* Few clusters logit: Score wild cluster bootstrap with Webb weights
boottest occrate, seed(10101) weight(webb) reps(999)

* Weak IV with independent observations: Wild bootstrap of Wald test
qui use mus207mepspresdrugs, clear
keep if _n <= 1500 & linc != .
global x2list totchr age female blhisp linc
ivregress 2sls ldrugexp $x2list (hi_empunion = ssiratio firmsz), ///
     noheader vce(robust)
boottest hi_empunion, seed(101010) reps(999)

* Weak IV with independent observations: Wild bootstrap of AR test
boottest, ar reps(999)
qui regress ldrugexp ssiratio firmsz $x2list, vce(robust)
test ssiratio firmsz    // The standard AR test for this example

******* 12.7: BOOTSTRAP USING BSAMPLE AND SIMULATE

* Program to do one bootstrap replication
program onebootrep, rclass
    version 17
    drop _all
    use mus212bootdata
    bsample
    poisson docvis chronic, vce(robust)
    return scalar tstar = (_b[chronic]-$theta)/_se[chronic]
end

* Now do 999 bootstrap replications
qui use mus212bootdata, clear
qui poisson docvis chronic, vce(robust)
global theta = _b[chronic]
global setheta = _se[chronic]
simulate tstar=r(tstar), seed(10101) reps(999) nodots  ///
  saving(percentilet2, replace): onebootrep

* Analyze the results to get the p-value
use percentilet2, clear
qui count if abs($theta/$setheta) < abs(tstar)
display "p-value = " r(N)/_N

* Program to do one bootstrap of B replications
program mybootstrap, rclass
    version 17
    use mus212bootdata, clear
    qui poisson docvis chronic, vce(robust)
    global theta = _b[chronic]
    global setheta = _se[chronic]
    simulate tstar=r(tstar), reps(999) nodots  ///
        saving(percentilet2, replace): onebootrep
    use percentilet2, clear
    qui count if abs($theta/$setheta) < abs(tstar)
    return scalar pvalue =  r(N)/_N
end

******* 12.8 ALTERNATIVE RESAMPLING SCHEMES (TO GET SE'S)

* Program to resample using bootstrap pairs
program bootpairs
    version 17
    drop _all
    use mus212bootdata
    bsample
    poisson docvis chronic
end

* Check the program by running once
bootpairs

* Bootstrap pairs for the parameters
simulate _b, seed(10101) reps(400) nodots: bootpairs
summarize

* Fit the model with original actual data and save estimates
qui use mus212bootdata, clear
qui nbreg docvis chronic
predict muhat
global alpha = e(alpha)

* Program for parametric bootstrap generating from negative binomial
program bootparametric, eclass
    version 17
    capture drop nu dvhat
    generate nu = rgamma(1/$alpha,$alpha)
    generate dvhat = rpoisson(muhat*nu)
    nbreg dvhat chronic
end

* Check the program by running once
set seed 10101
bootparametric

* Parametric bootstrap for the parameters
simulate _b, seed(10101) reps(400) nodots: bootparametric
summarize

* Program for residual bootstrap for OLS with iid errors
qui use mus212bootdata, clear
qui regress docvis chronic
predict uhat, resid
keep uhat
save residuals, replace
program bootresidual
    version 17
    drop _all
    use residuals
    bsample
    merge 1:1 _n using mus212bootdata
    regress docvis chronic
    predict xb
    generate ystar =  xb + uhat
    regress ystar chronic
end

* Check the program by running once
bootresidual

* Residual bootstrap for the parameters
simulate _b, seed(10101) reps(400) nodots: bootresidual
summarize

* Wild bootstrap for OLS with iid errors
qui use mus212bootdata, clear
program bootwild
    version 17
    drop _all
    use mus212bootdata
    regress docvis chronic
    predict xb
    predict u, resid
    generate ustar = -0.618034*u
    replace ustar = 1.618034*u if runiform() > 0.723607
    generate ystar =  xb + ustar
    regress ystar chronic
end

* Check the program by running once
bootwild

* Wild bootstrap for the parameters
simulate _b, seed(10101) reps(400) nodots: bootwild
summarize

******* 12.9 THE JACKKNIFE

* Jackknife estimate of standard errors
qui use mus212bootdata, replace
poisson docvis chronic, vce(jackknife, mse nodots)

* Not included - gives the same results
jacknife, mse nodots: poisson docvis chronic

*** Erase files created by this program and not used elsewhere in the book

erase residuals.dta
erase percentilet.dta
erase percentilet2.dta

********** END

	
	
	
