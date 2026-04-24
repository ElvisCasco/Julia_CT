* trpanel5.do for Stata 12.0  November 2013 

capture log close
log using trpanel5.txt, text replace

********** OVERVIEW OF trpanel5.do **********

* Based on "Microeconometrics using Stata" 
* by A. Colin Cameron and Pravin K. Trivedi (2008)

* Chapter 18
* NONLINEAR PANEL-DATA EXAMPLE
* PANEL BINARY OUTCOME MODELS
* PANEL TOBIT MODELS
* PANEL COUNT MODELS

* To run you need files
*   mus18data.dta    
* in your directory
* Stata user-written command
*   xtpqml
* is used

********** SETUP **********

set more off
version 12
clear all
set linesize 90
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus18data.dta    
* Rand Health Insurance Experiment data 
* Essentially same data as in P. Deb and P.K. Trivedi (2002)
* "The Structure of Demand for Medical Care: Latent Class versus
* Two-Part Models", Journal of Health Economics, 21, 601-625
* except that paper used different outcome (counts rather than $)
* Each observation is for an individual over a year.
* Individuals may appear in up to five years.
* All available sample is used except only fee for service plans included.
* If panel data used then clustering is on id (person id)

********** BINARY LOGIT EXAMPLE **********

* Begin with just analyse year 1
use mus18data.dta, clear
keep if year == 1
describe dmdu ndisease
sum dmdu ndisease

logit dmdu ndisease, nolog
margins, dydx(*)

predict plogit
sort ndisease
graph twoway (scatter dmdu ndisease, msize(small) jitter(8)) ///
   (line plogit ndisease, clstyle(p1) lwidth(thick))  ///
   (lfit dmdu ndisease, clstyle(p2) lwidth(thick)), scale (1.2)

********** NONLINEAR PANEL-DATA EXAMPLE

* Describe dependent variables and regressors
use mus18data.dta, clear
describe dmdu med mdu lcoins ndisease female age lfam child id year

* Summarize dependent variables and regressors
summarize dmdu med mdu lcoins ndisease female age lfam child id year

* Panel description of dataset 
xtset id year
xtdescribe

* Panel summary of time-varying regressors
xtsum age lfam child

* Panel summary of dependent variable
xtsum dmdu

* Year-to-year transitions in whether visit doctor
xttrans dmdu

* Correlations in the dependent variable
corr dmdu l.dmdu l2.dmdu

* Logit cross-section with panel-robust standard errors
logit dmdu lcoins ndisease female age lfam child, vce(cluster id) nolog

* Pooled logit cross-section with exchangeable errors and panel-robust VCE
xtlogit dmdu lcoins ndisease female age lfam child, pa corr(exch) vce(robust) nolog

* Logit random-effects estimator
xtlogit dmdu lcoins ndisease female age lfam child, re nolog

* Logit fixed-effects estimator
xtlogit dmdu lcoins ndisease female age lfam child, fe nolog 

* Logit mixed-effects estimator (same as xtlogit, re)
* xtmelogit dmdu lcoins ndisease female age lfam child || id:

* Panel logit estimator comparison
global xlist lcoins ndisease female age lfam child
quietly logit dmdu $xlist, vce(cluster id)
estimates store POOLED
quietly xtlogit dmdu $xlist, pa corr(exch) vce(robust)
estimates store PA
quietly xtlogit dmdu $xlist, re     // SEs are not cluster-robust
estimates store RE
quietly xtlogit dmdu $xlist, fe     // SEs are not cluster-robust
estimates store FE
estimates table POOLED PA RE FE, equations(1) se b(%8.4f) stats(N ll) stfmt(%8.0f)

********** PANEL TOBIT MODELS

* Panel summary of dependent variable
xtsum med

* Tobit random-effects estimator
xttobit med lcoins ndisease female age lfam child, ll(0) nolog

********** PANEL COUNT MODELS

* Panel summary of dependent variable
xtsum mdu

* Year-to-year transitions in doctor visits
generate mdushort = mdu
replace mdushort = 4 if mdu >= 4
xttrans mdushort

// Not included in book
corr mdu L.mdu

* Pooled Poisson estimator with cluster-robust standard errors
poisson mdu lcoins ndisease female age lfam child, vce(cluster id)

* Poisson PA estimator with unstructured error correlation and robust VCE
xtpoisson mdu lcoins ndisease female age lfam child, pa corr(unstr) vce(robust)

// Not included in book
* Poisson random-effects estimator with default standard errors
xtpoisson mdu lcoins ndisease female age lfam child, re

* Poisson random-effects estimator with cluster-robust standard errors
* xtpoisson mdu lcoins ndisease female age lfam child, re vce(boot, reps(400) seed(10101) nodots)
/* Previous command commented out as takes a long time
It gives the following output
Random-effects Poisson regression               Number of obs      =     20186
Group variable: id                              Number of groups   =      5908

Random effects u_i ~ Gamma                      Obs per group: min =         1
                                                               avg =       3.4
                                                               max =         5

                                                Wald chi2(6)       =    534.34
Log likelihood  = -43240.556                    Prob > chi2        =    0.0000

                                   (Replications based on 5908 clusters in id)
------------------------------------------------------------------------------
             |   Observed   Bootstrap                         Normal-based
         mdu |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      lcoins |  -.0878258   .0081916   -10.72   0.000     -.103881   -.0717706
    ndisease |   .0387629   .0024574    15.77   0.000     .0339466    .0435793
      female |   .1667192   .0376166     4.43   0.000     .0929921    .2404463
         age |   .0019159   .0016831     1.14   0.255     -.001383    .0052148
        lfam |  -.1351786   .0338651    -3.99   0.000     -.201553   -.0688042
       child |   .1082678   .0537636     2.01   0.044     .0028931    .2136426
       _cons |   .7574177   .0827935     9.15   0.000     .5951454      .91969
-------------+----------------------------------------------------------------
    /lnalpha |   .0251256   .0257423                     -.0253283    .0755796
-------------+----------------------------------------------------------------
       alpha |   1.025444   .0263973                      .9749897    1.078509
------------------------------------------------------------------------------
Likelihood-ratio test of alpha=0: chibar2(01) =  3.9e+04 Prob>=chibar2 = 0.000
*/

// Not included in book
* Poisson random-effects estimator with normal intercept and default standard errors
xtpoisson mdu lcoins ndisease female age lfam child, re normal

// Not included in book
* Poisson random-effects estimator with normal intercept and normal slope for one parameter
* xtmepoisson mdu lcoins ndisease female age lfam child || id: NDISEASE

* Poisson fixed-effects estimator with cluster-robust standard errors
* xtpoisson mdu lcoins ndisease female age lfam child, fe vce(boot, reps(400) seed(10101) nodots)
/* Previous command commented out as takes a long time
It gives the following output
Conditional fixed-effects Poisson regression    Number of obs      =     17791
Group variable: id                              Number of groups   =      4977

                                                Obs per group: min =         2
                                                               avg =       3.6
                                                               max =         5

                                                Wald chi2(3)       =      4.39
Log likelihood  = -24173.211                    Prob > chi2        =    0.2221

                                   (Replications based on 4977 clusters in id)
------------------------------------------------------------------------------
             |   Observed   Bootstrap                         Normal-based
         mdu |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
         age |  -.0112009   .0094595    -1.18   0.236    -.0297411    .0073394
        lfam |   .0877134   .1152712     0.76   0.447     -.138214    .3136407
       child |   .1059867   .0758987     1.40   0.163    -.0427721    .2547454
------------------------------------------------------------------------------
*/

// Not included in book
xtpqml mdu lcoins ndisease female age lfam child, fe

xtpoisson mdu lcoins ndisease female age lfam child, fe vce(robust)

// Not included in book
* Poisson fixed-effects estimator with default standard errors
xtpoisson mdu lcoins ndisease female age lfam child, fe i(id)  

* Comparison of Poisson panel estimators
quietly xtpoisson mdu lcoins ndisease female age lfam child, pa corr(unstr) vce(robust)
estimates store PPA_ROB
quietly xtpoisson mdu lcoins ndisease female age lfam child, re
estimates store PRE
quietly xtpoisson mdu lcoins ndisease female age lfam child, re normal
estimates store PRE_NORM
quietly xtpoisson mdu lcoins ndisease female age lfam child, fe vce(robust)
estimates store PFE_ROB
quietly xtpoisson mdu lcoins ndisease female age lfam child, fe
estimates store PFE
estimates table PPA_ROB PRE PRE_NORM PFE_ROB PFE, equations(1) b(%8.4f) se stats(N ll) stfmt(%8.0f)

// Not included in book
* Negative binomial pooled estimator with default standard errors
nbreg mdu lcoins ndisease female age lfam child
* Negative binomial pooled estimator with het robust standard errors
nbreg mdu lcoins ndisease female age lfam child, vce(robust)
* Negative binomial pooled estimator with cluster-robust standard errors
nbreg mdu lcoins ndisease female age lfam child, cluster(id)
* Negative binomial population-averaged estimator with equicorrelated errors
xtnbreg lcoins ndisease female age lfam child, pa corr(exch) vce(robust)
* Negative binomial random-effects estimator with default standard errors
xtnbreg mdu ndisease female age lfam child, re i(id)  
* Negative binomial fixed-effects estimator with default standard errors
xtnbreg mdu ndisease female age lfam child, fe i(id)  

* Comparison of negative binomial panel estimators
quietly xtpoisson mdu lcoins ndisease female age lfam child, pa corr(exch) vce(robust)
estimates store PPA_ROB
quietly xtnbreg mdu lcoins ndisease female age lfam child, pa corr(exch) vce(robust)
estimates store NBPA_ROB
quietly xtnbreg mdu lcoins ndisease female age lfam child, re
estimates store NBRE
quietly xtnbreg mdu lcoins ndisease female age lfam child, fe
estimates store NBFE
estimates table PPA_ROB NBPA_ROB NBRE NBFE, equations(1) b(%8.4f) se stats(N ll) stfmt(%8.0f)

********** CLOSE OUTPUT
* log close
* clear
* exit




