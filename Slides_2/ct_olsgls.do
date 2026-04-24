* ct_osgls.do  April 2013 for Stata version 12.0

********** OVERVIEW OF ct_olsgls.do **********

log using ct_olsgls.txt, text replace

* STATA Program 
* For A. Colin Cameron and Pravin Trivedi "Lectures in Microeconometrics"
* OLS and GLS
*    OLS with data
*    OLS simulation
*    FGLS with data

* To run you need file
*   mus10data.dta    
* in your directory

********** SETUP **********

set more off
version 12.0
set mem 10m
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* 2002 Medical Expenditure Panel Survey (MEPS)
* U.S. individuals aged 25-64 years working in private sector 
* but not self-employed
* and not receiving public insurance (Medicare and Medicaid) 
* Data due to Deb, Munkin and Trivedi (2006)

******* OLS WITH DOCTOR VISITS DATA

* Read in data, select 2002 data, describe and summarize key variables
use mus10data.dta, clear
quietly keep if year02==1
describe docvis private chronic female income
summarize docvis private chronic female income

* OLS regression with default standard errors
regress docvis private chronic female income

* OLS regression with robust standard errors
regress docvis private chronic female income, vce(robust)

* Comparison of standard errors
quietly regress docvis private chronic female income
estimates store DEFAULT
quietly regress docvis private chronic female income, vce(robust)
estimates store ROBUST
estimates table DEFAULT ROBUST, b(%9.4f) se(%9.3f) stats(N r2 F)

* Wald test of restrictions
quietly regress docvis private chronic female income, vce(robust) noheader
test (private = 0) (chronic = 0)

******* CONSISTENCY AND ASYMPTOTIC NORMAILITY FOR OLS

* Small sample: parameters differ from dgp values
clear all
quietly set obs 30
set seed 10101
quietly generate double x = rchi2(1)   
quietly generate y = 1 + 2*x + rchi2(1)-1     // demeaned chi^2 error 
regress y x, noheader

* Consistency: Large sample: parameters are very close to dgp values
clear all
quietly set obs 100000
set seed 10101
quietly generate double x = rchi2(1)   
quietly generate y = 1 + 2*x + rchi2(1)-1     // demeaned chi^2 error 
regress y x, noheader

* Central limit theorem
* Write program to obtain betas for one sample of size numobs (= 150)
program chi2data, rclass
    version 10.1  
    drop _all
    set obs $numobs
    generate double x = rchi2(1)   
    generate y = 1 + 2*x + rchi2(1)-1         // demeaned chi^2 error 
    regress y x
    return scalar b2 =_b[x]
    return scalar se2 = _se[x]
    return scalar t2 = (_b[x]-2)/_se[x]
    return scalar r2 = abs(return(t2))>invttail($numobs-2,.025)
    return scalar p2 = 2*ttail($numobs-2,abs(return(t2)))
end
* Run this program 1,000 times to get 1,000 betas etcetera
* Results differ from MUS (2008) as MUS did not reset the seed to 10101
* First define global macro numobs for sample size 
global numobs 150
set seed 10101 
quietly simulate b2f=r(b2) se2f=r(se2) t2f=r(t2) reject2f=r(r2) p2f=r(p2),  ///
   reps(1000) saving(chi2datares, replace) nolegend nodots: chi2data
* Summarize the 1,000 sample means 
summarize b2f se2f t2 reject2f p2f
mean b2f se2f t2 reject2f p2f

* Draw histogram of the t-test statistic of H0: b2 = 2
quietly histogram t2f, normal   ///
   xtitle("t-statistic for slope coeff from many samples")
graph export ct_olsglsclt.wmf, replace

******* FGLS WITH DOCTOR VISITS DATA

* Read in data, select 2002 data, describe and summarize key variables
use mus10data.dta, clear
quietly keep if year02==1

* FGLS: sigmahat^2 = exp(x'ghat) where ghat from NLS of uhat^2 on exp(x'g)
quietly regress docvis private chronic female income
quietly predict uhat, resid
quietly generate uhatsq = uhat^2               // compute squared residual
quietly generate one = 1                       
quietly nl (uhatsq = exp({xb: private chronic female income one})), nolog 
quietly predict varu, yhat                     // compute sigmahat^2
regress docvis private chronic female income [aweight=1/varu], noheader 

* WLS estimator is FGLS with robust estimate of VCE
regress docvis private chronic female income [aweight=1/varu], vce(robust) noheader

* Compare to OLS
regress docvis private chronic female income, vce(robust) noheader

********** CLOSE OUTPUT
* log close
* clear
* exit


