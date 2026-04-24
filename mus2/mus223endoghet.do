* mus223endoghet.do  for Stata 17
cap log close

********** OVERVIEW OF mus223endoghet.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press  
 
* Chapter 23.  PARAMETRIC MODELS FOR ENDOGENEITY AND HETEROGENEITY 
*   23.3: EMPIRICAL EXAMPLES OF FMMS
*   23.4: NONLINEAR MIXED-EFFECTS MODELS
*   23.5: STRUCTURAL EQUATION MODELS FOR LINEAR STRUCTURAL EQUATION MODEL
*   23.6: GSEM
*   23.7: ERM COMMANDS FOR ENDOGENEITY AND SELECTION 

* To run you need files
*    mus203mepsmedexp.dta
*    mus206vlss.dta
*    mus217hrs.dta
*    mus218hk.dta
*    mus219mepsambexp.dta
*    mus220mepsdocvis.dta
*    mus220mepsemergroom.dta
*    mus221mccall.dta 
*    mus223gsemfmmexample.dta
* in your directory

* No community-contributed command is used

********** SETUP **********

clear all
set linesize 82
set scheme s1mono   /* Used for graphs */

********** DATA DESCRIPTION **********

* mus203mepsmedexp.dta
* mus219mepsambexp.dta
* mus220mepsdocvis.dta
* mus220mepsemergroom.dta
* mus223gsemfmmexample.dta are authors' extracts 
* from MEPS (Medical Expenditure Panel Survey)
* for individuals 65 years and older in U.S. Medicare in 2003

* mus217hrs.dta comes from the Health and Retirement Study (HRS).
* wave 5 (2002) and is restricted to Medicare beneficiaries.

* mus218rhie.dta is from Rand Health Insurance Experiment data 
* Essentially same data as in P. Deb and P.K. Trivedi (2002)
* "The Structure of Demand for Medical Care: Latent Class versus
* Two-Part Models", Journal of Health Economics, 21, 601-625
* except that article used different outcome (counts rather than $)
* Each observation is for an individual over a year.
* Individuals may appear in up to five years.
* All available sample is used except only fee for service plans included.
* If panel data used then clustering is on id (person id)

* mus224mccall.dta is from B.P. McCall (1996), 
* "Unemployment Insurance Rules, Joblessness, and Part-time Work," 
* Econometrica, 64, 647-682.

************ CHAPTER 23.3: FMM: EMPIRICAL EXAMPLES OF FINITE MIXTURES MODELS

* EXAMPLE 1 : GAMMA REGRESSION MIXTURE FOR MEDICAL EXPENDITURE

* Gamma: Read chapter 3 example data and truncate expenditure at zero
qui use mus203mepsmedexp
keep if totexp > 0  
describe totexp totchr age female educyr private hvgg

* Gamma: Standard gamma regression (one component) as benchmark
glm totexp totchr age female educyr private hvgg, family(gamma) link(log) ///
    vce(robust) nolog

// not included
estat ic

// Not included - compare to default standard errors
glm totexp totchr age female educyr private hvgg, family(gamma) link(log) /// 
    nolog

// Not included - compare gamma to lognormal using same regressros as ch. 3
glm totexp suppins phylim actlim totchr age female income,  /// 
    family(gamma) link(log) vce(robust)
estimates store GAMMA
glm totexp suppins phylim actlim totchr age female income,  /// 
    family(gamma) link(log)
estimates store GAMDEF 
glm ltotexp suppins phylim actlim totchr age female income, vce(robust)
estimates store LNORMAL
glm ltotexp suppins phylim actlim totchr age female income
estimates store LNORMDEF
estimates table GAMMA GAMDEF LNORMAL LNORMDEF, b se eq(1)

* Gamma: Two-component finite mixture gamma regression 
fmm 2, nolog vce(robust): glm totexp totchr age female educyr private hvgg, ///
    family(gamma) link(log)

// not included
estat ic


* Gamma: Latent class marginal probabilities and latent class distribution means
estat lcprob
estat lcmean

* Gamma: ME of change in totchr in one- and two-component models
margins, dydx(totchr) noatlegend    // AME in two-component model
qui glm totexp totchr age female educyr private hvgg, family(gamma) ///
    link(log) vce(robust)
margins, dydx(totchr) noatlegend    // AME in one-component model

* Gamma: Generate and graph latent class predicted expenditures
qui fmm 2, vce(robust): glm totexp totchr age female educyr private hvgg, ///
    family(gamma) link(log)
estimates store FMM2
predict mu*
summarize mu1 mu2
qui histogram mu1, saving(graph1.gph, replace)
qui histogram mu2, saving(graph2.gph, replace) 

graph combine graph1.gph graph2.gph, iscale(1.1) rows(1) ysize(2.5) xsize(6) /// 
     ycommon xcommon

* Gamma: Separate AME of totchr for each component
margins, dydx(totchr) predict(mu class(1)) predict(mu class(2)) noatlegend


* Gamma: Generate and graph latent class posterior probabilities 
predict postprob*, classpost
summarize postprob*
qui histogram postprob1, saving(graph1.gph, replace)
qui histogram postprob2, saving(graph2.gph, replace) 

graph combine graph1.gph graph2.gph, iscale(1.2) rows(1) ysize(2.5) ///
    xsize(6) ycommon xcommon

// Not included
estat lcprob, classposteriorpr

* Gamma: Three-component finite mixture gamma regression 
qui fmm 3, vce(robust): glm totexp totchr age female educyr private hvgg, ///
    family(gamma)
estimates store FMM3 
estat lcprob
estat lcmean

* Gamma: Compare two-and three-component finite mixture gamma models
estimates stats FMM2 FMM3  

* Gamma: Two-component model with class probabilities varying with regressors 
fmm 2, nolog lcprob(female age) vce(robust): glm totexp totchr age female ///
    educyr private hvgg, family(gamma) link(log) 
estat lcmean 

// Not included
estat ic

* EXAMPLE 2: LOGIT REGRESSION MIXTURE FOR HEALTH INSURANCE CHOICE

clear all
* Logit: Read binary outcome chapter example data 
qui use mus217hrs, clear
describe ins hstatusg hhincome educyear married

// Not included 
logit ins hstatusg hhincome educyear married
estat ic

* Logit: Two-component finite mixture logit regression and predictions
fmm 2, nolog vce(robust): logit ins hstatusg hhincome educyear married
predict classpr*
twoway (histogram classpr1, width(.05)) (histogram classpr2,            ///
    fcolor(white) width(.05)), saving(graph1.gph, replace) legend(off)  ///
    xtitle("Predicted marginal probabilities for the two classes") scale(1.2)


// Not included 
estat ic

* Logit: Latent class marginal probabilities and latent class distribution means
estat lcmean
estat lcprob, classposteriorpr

* Logit: ME of change in totchr in one- and two-component models
qui fmm 2, nolog vce(robust): logit ins hstatusg hhincome educyear married
margins, dydx(hhincome)
qui logit ins hstatusg hhincome educyear married, nolog
margins, dydx(hhincome)

* Logit: Obtain coefficient legend necessary for tests
qui fmm 2, vce(robust): logit ins hstatusg hhincome educyear married
fmm, coeflegend

* Logit: Tests of coefficient equality and of regressor relevance
test (_b[ins:1.Class#c.hstatusg] = _b[ins:2.Class#c.hstatusg]) /// 
     (_b[ins:1.Class#c.hhincome] = _b[ins:2.Class#c.hhincome]) /// 
     (_b[ins:1.Class#c.educyear] = _b[ins:2.Class#c.educyear]) ///  
     (_b[ins:1.Class#c.married] =  _b[ins:2.Class#c.married]), mtest

// Not included - the Sidal correction for multiple testing
test (_b[ins:1.Class#c.hstatusg] = _b[ins:2.Class#c.hstatusg]) /// 
     (_b[ins:1.Class#c.hhincome] = _b[ins:2.Class#c.hhincome]) /// 
     (_b[ins:1.Class#c.educyear] = _b[ins:2.Class#c.educyear]) ///  
     (_b[ins:1.Class#c.married] =  _b[ins:2.Class#c.married]), mtest

* EXAMPLE 3: MULTINOMIAL LOGIT REGRESSION MIXTURE FOR FISHING MODE CHOICE
	 
* MNL: Read multinomial outcome chapter example data 
qui use mus218hk, clear
describe mode income
summarize mode income 

// Not included - MNL to get AIC, BIC
mlogit mode income, vce(robust)
estat ic 

* MNL: Two-component finite mixture logit regression 
fmm 2, nolog vce(robust): mlogit mode income

// Not incluided
estat ic

* MNL: Latent class marginal probabilities and distribution means
estat lcmean
estat lcprob, classposteriorpr

* MNL: ME of regressor changes in two-component model
margins, dydx(*)

* EXAMPLE 4: TOBIT REGRESSION MIXTURE FOR LOG HEALTH EXPENDITURES

* Tobit: Two-component finite mixture tobit
qui use mus219mepsambexp, clear
global xlist age female educ blhisp totchr ins
fmm 2, vce(robust) nolog: tobit ambexp $xlist, ll(0)
estimates store FMM2

* Tobit: AME for two-component finite mixture tobit for E[y*|x]
margins, dydx(*)

* Tobit: Compare FMM2 to regular tobit using AIC and BIC
qui tobit ambexp $xlist, ll(0) vce(robust)
estimates store FMM1
estimates stats FMM1 FMM2

* EXAMPLE 5: POISSON REGRESSION MIXTURE FOR DOCTOR VISITS

* Poisson: Two-component finite mixture Poisson using count chapter data
use  mus220mepsdocvis, clear
global xlist private medicaid age age2 educyr actlim totchr 
qui fmm 2, vce(robust): poisson docvis $xlist
estimates store fmm2pois
estat lcmean
estat lcprob

* Poisson: ME of regressor changes in two-component model
margins, dydx(totchr)
qui poisson docvis $xlist, vce(robust)
estimates store poisson 
margins, dydx(totchr)

// Output not given but results stated in text
* ME of totchr evaluated at sample means for FMM and Poisson models 
sum $xlist
qui fmm 2, vce(robust): poisson docvis $xlist
margins, dydx(totchr) at(medicaid=.166712 age=90 age2=8100 ///
    educyr=11.18031 actlim=.333152 totchr=1.843351)
qui poisson docvis $xlist, vce(robust)
margins, dydx(totchr) at(medicaid=.166712 age=90 age2=8100 ///
    educyr=11.18031 actlim=.333152 totchr=1.843351)

* Poisson: Three-component finite mixture Poisson regression
qui fmm 3, nolog vce(robust): poisson docvis $xlist
estat lcmean
estat lcprob
estimates store fmm3pois

* Poisson: Goodness-of-fit for two- and three-component models
estimates stats poisson fmm2pois fmm3pois 

// Not included
* Negative binomial: Finite mixture models and comparison to Poisson
qui nbreg docvis $xlist, vce(robust)
estimates store nb
qui fmm 2, nolog vce(robust): nbreg docvis $xlist
estimates store fmm2nb
qui fmm 3, nolog vce(robust): nbreg docvis $xlist
estimates store fmm3nb
estimates stats nb fmm2nb fmm3nb poisson fmm2pois fmm3pois 

// Not included - robust standard errors are different for lcmean and lcprob
fmm 2, vce(robust) nolog: poisson docvis totchr
estat lcmean
estat lcprob
fmm 2, nolog: poisson docvis totchr
* estimates store fmm2pdef 
estat lcmean
estat lcprob
* estimates table fmm2pdef fmm2pois, b se

* EXAMPLE 6: POINT-MASS COUNT REGRESSION WITH ZEROS

// Output not included
* fmm point-mass two-component command same as zero-inflated negative binomial
qui use mus220mepsemergroom, clear
fmm: (pointmass er, lcprob(age actlim totchr)) (nbreg er age actlim totchr) 
zinb er age actlim totchr, inflate(age actlim totchr) ltolerance(1e-9)

* EXAMPLE 7: WEIBULL REGRESSION MIXTURE FOR UNEMPLOYMENT DURATION 

* Weibull: Read in and describe unemployment duration data from survival chapter
qui use mus221mccall, clear   
describe spell ui logwage 
summarize spell ui logwage

* Weibull: Standard Weibull regression (one component) as benchmark
stset spell, fail(censor1=1)
qui streg i.ui logwage, nohr vce(robust) dist(weibull) nolog
margins, dydx(ui)

* Weibull: Two-component finite mixture Weibull regression
fmm 2, nolog vce(robust): streg i.ui logwage, distribution(weibull)

* Weibull: Latent class probabilities and means for two-component model
estat lcprob
estat lcmean

* Weibull: MEs for two-component model
margins, dydx(ui)

// Not included
qui streg i.ui logwage, nohr vce(robust) dist(weibull) nolog
estimates store WEIBULL
qui fmm 2, nolog vce(robust): streg i.ui logwage, distribution(weibull)
estimates store WEIBFMM2
qui fmm 2, nolog vce(robust): streg i.ui logwage, distribution(weibull)
estimates store WEIBFMM3
estimates stats WEIBULL WEIBFMM2 WEIBFMM2

*********** CHAPTER 23.4: ME: NONLINEAR MIXED EFFECTS MODELS

* Read in data and drop a few observations
qui use mus206vlss, clear
drop if missing(lnhhexp) | (lnhhexp > 2.579681 & lnhhexp < 2.579683)
summarize pharvis lnhhexp illness commune

* Mixed model: Poisson, 2nd level commune, coeff of intercept and illness vary
mepoisson pharvis lnhhexp illness || commune: illness, vce(robust) 

* Wald test of the random coefficients
test ([/]var(illness[commune])=0) ([/]var(_cons[commune])=0)

* Mixed model: Compute various models and standard errors
qui poisson pharvis lnhhexp illness, vce(cluster commune)
estimates store p_clu
qui mepoisson pharvis lnhhexp illness || commune: illness
estimates store mep_def
qui mepoisson pharvis lnhhexp illness || commune: illness, vce(cluster commune)
estimates store mep_rob
qui menbreg pharvis lnhhexp illness || commune: illness
estimates store menb_def
qui menbreg pharvis lnhhexp illness || commune: illness, vce(robust)
estimates store menb_rob

* Mixed model: Compare with standard Poisson and get robust standard errors
estimates table p_clu mep_def mep_rob menb_def menb_rob,     ///
   eq(1) b(%8.4f) se stfmt(%8.1f) stats(ll N aic chi2_c df_c) 

* Various predictions of random effects and of mean
estimates restore mep_rob
predict re_ebmeans*, reffects ebmeans 
predict mu_fixed, mu conditional(fixedonly)
predict mu_ebmeans, mu conditional(ebmeans)
predict mu_ebmode, mu conditional(ebmode)
predict mu_marg, mu marginal
summarize re_ebmeans* pharvis mu_fixed mu_ebmeans mu_ebmode mu_marg  

* Marginal effects after integrating out random effects
margins, dydx(*) predict(mu conditional(fixedonly))

*********** CHAPTER 23.5: SEM: LIMEAR STRUCTURAL EQUATION MODELS 

* SEM method ml: Linear regression example
qui use mus220mepsdocvis, clear
sem (docvis <- private medicaid age age2 educyr actlim totchr)

* SEM methods compared with OLS: Linear regression example
qui sem (docvis <- private medicaid age age2 educyr actlim totchr)
estimates store SEMMLE
qui sem (docvis <- private medicaid age age2 educyr actlim totchr), vce(robust)
estimates store SEMMLErob
qui sem (docvis <- private medicaid age age2 educyr actlim totchr), method(adf)
estimates store SEMADF
qui regress docvis private medicaid age age2 educyr actlim totchr
estimates store OLS
qui regress docvis private medicaid age age2 educyr actlim totchr, vce(robust)
estimates store OLSrob

* SEM methods compared with OLS: Table of estimates
estimates table SEMMLE SEMMLErob SEMADF OLS OLSrob , ///
   eq(1) b(%10.4f) se stfmt(%10.2f) stats(rmse r2 ll critvalue N)

// Not included - get precise rmse
qui regress docvis private medicaid age age2 educyr actlim totchr
di e(rmse) "   " e(rmse)^2

* SEM measurement error example: Generate the data
clear
qui set obs 1000
set seed 10101
gen x = rnormal(0,1)
gen z = x + rnormal(0,1)
gen y = 1 + 1*x + 1*z + rnormal(0,1)
gen x1 = x + rnormal(0,1)
gen x2 = x + rnormal(0,1)
gen x3 = x + rnormal(0,1)
gen x4 = x + rnormal(0,1)

* SEM measurement error: OLS (using SEM) with true regressor x consistent
sem y <- z x, nolog vce(robust)

* SEM measurement error: OLS with mismeasured regressor x1 inconsistent
sem y <- z x1, nolog  vce(robust)

* SEM measurement error: IV with x2, x3, x4 instruments with x1 consistent
sem (y <- z x1) (x1 <- z x2 x3 x4), covar(e.y*e.x1) nolog noheader vce(robust)

// Not included - liml and 2sls
ivregress liml y z (x1 = x2), vce(robust)
ivregress 2sls y z (x1 = x2), vce(robust)

* SEM measurement error: Measurement model for X using 3 measures x1, x2, x3
sem (x1 x2 x3 <- X), nolog

/* Variations not included
* SEM measurement error: Normalize X at 2 rescales some but not all parameters
sem (x1 <- X@2) (x2 <- X) (x3 <- X), nolog
* SEM measurement error: Similarly if set the variance
sem (x1 x2 x3 <- X), nolog var(X@0.25)
* SEM measurement error: SEM model with two measures converged but had problems
sem (x1 x2 <- X)
* SEM measurement error: SEM with four measures is over-identified
sem (x1 x2 x3 x4 <- X), nolog
*/

* SEM meas error: regress y on z and x controlling for measurement error in x
sem (x1 x2 x3 <- X) (y <- z X), nolog vce(robust)

/* Variations not included
* SEM measurement error: Rescale the latent variable
sem (x1 x2 x3 <- X) (y <- z X), nolog var(X@4)
* SEM measurement error: adf gives similar
* sem (x1 x2 x3 <- X) (y <- z X), adf
* SEM measurement error: Reduce to 2 measures still works as 4x5/2=10 > 9 parameters
sem (x1 x2 <- X) (y <- z X), nolog 
* SEM measurement error: With robust standard errors cannot do OIR test
sem (x1 x2 x3 <- X) (y <- z X), nolog vce(robust)
* SEM measurement error: sbentler gives the overidentifying restrictions test
sem (x1 x2 x3 <- X) (y <- z X), nolog vce(sbentler)
* SEM measurement error: Now just one measure
sem (x1 <- X) (y <- z X), nolog
* SEM measurement error: Now just two measures
sem (x1 x2 <- X) (y <- z X), nolog
*/

// Output not included
* SEM measurement error: Control for meas. error using known reliability ratio
sem (y <- z X) (x1 <- X), reliability(x1 0.5) nolog

* DGP for just-identified two-equation simultaneous equations example
clear all 
set seed 10101
drawnorm u1 u2, n(1000) corr(1, .7, 1) cstorage(lower)
drawnorm x1 x2, n(1000) corr(1, .3, 1) cstorage(lower) 
generate y1 =(1/(1-.25))*(x1 + 1*(x2+u2) + u1)   // Reduced form for y1
generate y2 = 0.25*y1 + 1*x2 + u2                // Generating y2 given y1

* SEM commands for 2SLS and 3SLS in just-identified model
qui ivregress 2sls y1 (y2=x2) x1, noheader
estimates store iv_2sls
qui sem (y1 <- y2 x1) (y2 <- x1 x2), covar(e.y1*e.y2)
estimates store sem_rf
qui reg3 (y1 = y2 x1) (y2 = y1 x2)
estimates store iv_3sls
qui sem (y1 <- y2 x1) (y2 <- y1 x2), covar(e.y1*e.y2)
estimates store sem_sf

* SEM estimators compared with 2SLS and 3SLS in just-identified model
estimates table iv_2sls sem_rf iv_3sls sem_sf, eq(1) b(%9.5f) se

*********** CHAPTER 23.6: GENERALIZED STRUCTURAL EQUATION MODELS
 
* GSEM: Poisson normal mixture model
qui use mus220mepsdocvis, clear
global xlist2 medicaid age age2 educyr actlim totchr
gsem (docvis <- private $xlist2 L, poisson), nolog vce(robust)

* GSEM: Compute Poisson, Poisson-normal, and negative binomial
qui poisson docvis private $xlist2
estimates store pois_def
qui poisson docvis private $xlist2, vce(robust) 
estimates store pois_rob
qui gsem (docvis <- private $xlist2 L, poisson)
estimates store normal 
qui gsem (docvis <- private $xlist2 L, poisson), vce(robust) 
estimates store norm_rob
qui nbreg docvis private $xlist2 
estimates store nb2_def
qui nbreg docvis private $xlist2, vce(robust) 
estimates store nb2_rob

* GSEM: Compare Poisson-normal with Poisson and negative binomial
estimates table pois_def pois_rob normal norm_rob nb2_def nb2_rob, ///
  eq(1) b(%8.4f) se stfmt(%8.1f) stats(ll N)

* GSEM: Endogenous regressor in a Poisson model
gsem (docvis <- private $xlist2 L, poisson) ///
   (private <- $xlist2 income ssiratio L), nolog vce(robust) 

* GSEM: Predict y controlling for latent variable L
predict mu_ebmeans, mu conditional(ebmeans)
predict mu_fixed, mu conditional(fixedonly)
predict mu_marg, mu marginal

* GSEM: Compare predictions of y
summarize docvis mu_ebmeans mu_fixed mu_marg
correlate docvis mu_ebmeans mu_fixed mu_marg

* GSEM: Read in, and describe data for bivariate dependent variables example 
qui use mus223gsemfmmexample, clear
global xlist privins medicaid numchron age male
summarize emr hosp $xlist

// Not included
correlate emr hosp
 
* GSEM: Fit two-component FMM bivariate negative binomial 
qui gsem (emr hosp <- $xlist), nbreg lclass(C 2)        ///
    startvalues(randomid, draws(5) seed(15)) vce(robust)
estat lcprob

* GSEM: Latent class means for two-component FMM bivariate negative binomial
estat lcmean

* GSEM: MEs of numchron
margins, dydx(numchron)

// Not included - the full output
gsem (emr hosp <- $xlist), nbreg lclass(C 2)           ///
    startvalues(randomid, draws(5) seed(15)) vce(robust)

// Not included - equivalence of single-equation fmm and gsem
fmm 2, nolog vce(robust): poisson emr $xlist
gsem (emr <- $xlist), poisson lclass(C 2) nolog vce(robust)


************ CHAPTER 23.7: ERM COMMANDS FOR ENDOGENEITY AND SELECTION 

* Read in data, define globals, and summarize key variables
qui use mus217hrs, clear
generate linc = log(hhincome)
global xlist female age age2 educyear married hisp white chronic adl hstatusg
global ivlist retire sretire

* Endogenous probit using eprobit ML estimator (an erm command)
eprobit ins $xlist, endogenous(linc = $xlist $ivlist) vce(robust) nolog

// not included
probit ins $xlist linc, vce(robust)

*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph

********** END
