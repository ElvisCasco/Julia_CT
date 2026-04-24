* mus224treat.do  for Stata 17

clear all 
cap log close

********** OVERVIEW OF mus224treat.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press  

********** SETUP **********

clear all
set scheme s1mono  /* Graphics scheme */
set linesize 82
 
* Chapter 24  TREATMENT EVALUATION FOR RCTS AND EXOGENOUS  

* 24.3: RCTS
* 24.4: REGRESSION IN AN RCT
* 24.8: OHIE EXAMPLE
* 24.9: TE ESTIMATES USING THE OHIE DATA
* 24.10: MULTILEVEL TES
* 24.11: CONDITIONAL QUANTILE TR

* To run you need files
*   mus224ohiesmallrecode.dta 
*   mus224mcbs.dta   
* in your directory

* community-contributed commands
*   poparms
*   qplot
* are used
* Additionally, pscore and atts appear but are commented out 
* due to long computational time

********** DATA DESCRIPTION **********

* mus224ohiesmallrecode.dta are from the Oregon Health Insurance Experiment
* downloaded from the MIT archive. 

* mus224mcbs.dta  are a subset of the data used in
* Qian Li and Pravin Trivedi (2015), "Adverse and advantageous selection in the 
* Medicare supplemental market: a Bayesian analysis of prescription drug
* expenditure, Health Economics, 192-211.

**** CHAPTER 24.3: RANDOMIZED CONTROLLED TRIALS

* Power calculation examples
* (1) Required sample size when m1=21; m2=23,24,25,26; sd=6
power twomeans 21 (23(1)26), sd(6)  
 
*power twomeans 21 (23(1)26), sd(6)   graph

* (2) Required sample size when m1=21; m2=24; sd=6; power=0.9
power twomeans 21 24, sd(6) power(.9)  

* (3) Required sample size when m1=21; m2=24; sd=6; power=0.9; one-sided
power twomeans 21 24, sd(6) onesided  

* (4) Required minimum TE size when m1=21; sd=6; N=200 
power twomeans 21, sd(6) power(0.8) n(200)

* (5) Power when m1=21; m2=25; sd=6; N=100
power twomeans 21 25, sd(6) n(100)

* (6) Power when m1=21; m2=25; N=100; sd=9; N=100 
power twomeans 21 25, sd(9) n(100)  

* (7) Required sample size when m1=21; m2=24; sd=6; cluster size 5; rho=0.2
power twomeans 21 24, sd(6) m1(5) m2(5) rho(0.2) 

******** CHAPTER 24.4: REGRESSION IN AN RCT

* RCT experiment based on generated data
clear all
* Simulated RCT: Exogenous covariate x; similar treat & control sample size
qui set obs 100
set seed 10101
generate x = rnormal(20,5)      // Exogenous regressor
set seed 10102
generate D = rbinomial(1,0.5)   // Treatment assignment
set seed 10103
generate u = rnormal(0,1)       // Model error
generate y = 1 + x + 2*D + u    // Outcome varies with treatment
summarize x D y u
pwcorr x D y u, star(0.05)

* ATE: Use ttest command (no controls)
ttest y, by(D) unequal

* ATE: OLS regress y on D (no controls)
regress y D, vce(robust)

* ATE: OLS regress y on D and x added as control
regress y x D, vce(robust)

* OLS regress y on D and x and interaction of D and x 
regress y i.D##c.x, vce(robust) noheader  

* Average predicted outcomes for D=0 and D=1 from interactive regression
margins D         // POMs

* ATE: Using margin, dydx command after interactive regression 
margins, dydx(D)  // ATE

* ATE: Equivalent results using teffects ra command
teffects ra (y x) (D), ate nolog

******** CHAPTER 24.8: OREGON HEALTH INSURANCE (OHIE)

* Variables: (1) outcomes (2) z: treatment related (3) x: outcome related
qui use mus224ohiesmallrecode, clear
global outcomes oop dowe ervisits  // Outcome variables
global y oop                       // Outcome variable for this chapter
global hh dhhsize2 dhhsize3        // Household size dummies
global wave dlotdraw* dsurvdraw*   // Lottery and survey draws
global zlist $hh $wave             // z variables for the treatment
global xlist dsmoke hhinc deduc2-deduc4 demploy2-demploy4  // x for outcome

* Variable descriptions 
describe $outcomes lottery medicaid household_id $zlist $xlist

* Summary statistics: Outcomes, treatments, and household size
summarize $outcomes lottery medicaid $hh

* Lottery dummy varies with household size and lottery and survey waves
regress lottery $zlist, vce(cluster household_id)

// Not included - lottery success by each variable in $zlist
foreach var in $zlist {
    bysort `var': summarize lottery
	}
	
* Summary statistics: xlist variables that may improve outcome model fit
summarize $xlist

* Quantile plot of outcome variable by treatment status
qui qplot $y if $y < 2000, over(lottery) clpattern(l _) recast(line)  ///
    ytitle("Out-of-pocket spending") legend(pos(11) ring(0) col(1))   ///
    title("Quantile plots for lottery") lwidth(medthick thick)        ///
    xtitle("Fraction  of the data")                                   ///
	lpattern(solid dash) saving(graph1.gph, replace) 
qui qplot $y if $y < 2000, over(medicaid) clpattern(l _) recast(line) ///
    ytitle("Out-of-pocket spending") legend(pos(11) ring(0) col(1))   ///
    title("Quantile plots for Medicaid") lwidth(medthick thick)       ///
    xtitle("Fraction  of the data")                                   ///
	lpattern(solid dash) saving(graph2.gph, replace) 

graph combine graph1.gph graph2.gph, ysize(2.5) xsize(6.0) iscale(1.2)

// Not included - simple difference in means test
ttest $y, by(lottery) unequal

* Regression of outcome on treatment with various controls
qui regress $y lottery, vce(robust)
estimates store Diff_rob
qui regress $y lottery, vce(cluster household_id)
estimates store Diff_clu
qui regress $y lottery $zlist, vce(cluster household_id)
estimates store zlist
qui regress $y lottery $xlist, vce(cluster household_id)
estimates store xlist
qui regress $y lottery $zlist $xlist, vce(cluster household_id)
estimates store Both
estimates table Diff_rob Diff_clu zlist xlist Both, keep(lottery) ///
   b(%8.4f) se stat(N r2) 

// Not included - test statistical significance of control regressors
regress $y lottery $zlist $xlist, vce(cluster household_id)
testparm $zlist
testparm $xlist

********* CHAPTER 24.9: TREATMENT EFFECTS ESTIMATES USING THE OHIE DATA

***** Regression adjustment 
* Regression-adjusted ATE using $zlist and $xlist
teffects ra ($y $xlist $zlist ) (lottery), nolog vce(cluster household_id)

* Regression-adjusted POMs using $zlist and $xlist
teffects ra ($y $xlist $zlist) (lottery), pomeans nolog vce(clu household_id)

* Regression-adjusted ATET using $zlist and $xlist
teffects ra ($y $xlist $zlist) (lottery), atet nolog vce(cluster household_id)

// Not included - ate
teffects ra ($y $xlist $zlist ) (lottery), nolog vce(cluster household_id) ///
    aequations

***** Inverse probability weighting, balance and overlap

* IPW ATE using $zlist
teffects ipw ($y) (lottery $zlist), nolog vce(cluster household_id)

tebalance density dhhsize2

/* tebalance overid takes a long time 91 iterations
So comment out but report the results in the text
* Formal test of balancing 
. tebalance overid 
Iteration 0:   criterion =  .00024449  
Iteration 192: criterion =  .00278518  (backed up)
Overidentification test for covariate balance
         H0: Covariates are balanced:
         chi2(16)     =  358.044
         Prob > chi2  =   0.0000
*/

* Overlap assumption - graph using teoverlap
teoverlap, ptlevel(1) kernel(triangle) bw(0.04)           ///
    title("Kernel density overlap") xtitle("Propensity score")   ///
    saving(graph1, replace)

* Overlap assumption - manual graph using logit prediction and histograms 
qui logit lottery $zlist, vce(cluster household_id)
qui predict lothat
summarize lothat
twoway (hist lothat if lottery==0, fcol(white)) (hist lothat if lottery==1), ///
    title("Histogram overlap") xtitle("Propensity score")               ///
    legend(label(1 "lottery=Not selected") label(2 "lottery=Selected")) ///
    saving(graph2, replace)

graph combine graph1.gph graph2.gph, ysize(2.5) xsize(6.0) iscale(1.1)

// Not included - count low and high predicted probabilities
count if lothat < 0.3 | lothat > 0.75

* Covariate balance summary
qui teffects ipw ($y) (lottery $zlist), nolog vce(cluster household_id)
tebalance summarize

***** Augmented IPW and IPW regression adjustment

* Augmented IPW ATE
teffects aipw ($y $xlist) (lottery $zlist), aequations nolog ///
    vce(cluster household_id)

* IPW with RA ATE
teffects ipwra ($y $xlist $zlist) (lottery $zlist), nolog ///
    vce(cluster household_id) 

***** Propensity score matching  

* PSM ATE
teffects psmatch ($y) (lottery $zlist)

***** Blocking and stratification

/* Comment out as takes a long time - selected output is included in the book
* Blocking and subsequent balancing tests using community-contributed command pscore
pscore lottery $zlist, logit blockid(myblock) pscore(pshat) ///
    numblo(8) level(0.001)
* Propensity score matching ATE using blocking and community-contributed command atts 
atts oop lottery, blockid(myblock) pscore(pshat) 
*/

***** Nearest neighbor matching

* NNM ATE
generate dhhbig = dhhsize2 + dhhsize3
teffects nnmatch ($y $wave) (lottery), ematch(dhhbig) metric(mahalanobis) 

teffects nnmatch ($y $zlist) (lottery), ///
    ematch(female_list) metric(euclidean) 

********* CHAPTER 24.10: MULTILEVEL TREATMENT EFFECTS

* Outcome - drug expenditure 
qui use mus224mcbs, clear
generate drugexp = aamttot
summarize drugexp

// count number of zero observations in dependent variable
count if drugexp == 0 

* Create treatment variable clevel - insurance ordered by increasing generosity
qui generate clevel= 0 if coverage == "Medicare" // Medicare 
qui replace clevel = 1 if coverage == "MMC"      // MMC (Medicare managed plan)
qui replace clevel = 2 if coverage == "Medigap"  // Medigap (Medicare suppl ins) 
qui replace clevel = 3 if coverage == "ESI"      // ESI (employer-sponsored)
tabulate clevel

* Variation in mean drug expenditure by type of insurance
generate dmedicare = clevel == 0
generate dmmc = clevel == 1
generate dmedigap = clevel == 2
generate desi = clevel == 3
regress drugexp dmmc dmedigap desi, noheader vce(robust)

****** Regression adjustment 

* Variation in mean drug expenditure by type of insurance
global xlist h_age h_male h_white income_c genhelth 
summarize $xlist

* RA ATE 
teffects ra (drugexp $xlist, poisson) (clevel), nolog

* RA ATET 
teffects ra (drugexp $xlist, poisson) (clevel), atet nolog

* RA POMs
teffects ra (drugexp $xlist, poisson) (clevel), pomeans nolog

* Contrasts of TEs after RA compared with base category
contrast r.clevel, nowald

* Contrasts of TEs after RA compared with adjacent category
contrast ar.clevel, nowald

****** Inverse probability weights

* Insurance regressors used for the MNL propensity scores 
global zlist prem1 prem2 prem3 prem4 penet
describe $zlist
summarize $zlist

* Predicted probabilities from the MNL model
qui mlogit clevel $zlist, base(0)
predict pshat1 pshat2 pshat3 pshat4
summarize dmedicare dmmc dmedigap desi pshat1 pshat2 pshat3 pshat4, sep(4)   

* Augmented IPW ATE 
teffects aipw (drugexp $xlist, poisson) (clevel $zlist)

* Check balance following use of the inverse-probability weights
tebalance summarize

******* CHAPTER 24.11 CONDITIONAL QUANTILE TREATMENT EFFECTS

// Code included in book without output
* IPW using community-contributed command poparms gives same results as teffects ipw
poparms (clevel $zlist) (drugexp), ipw
contrast r.clevel, nowald
teffects ipw (drugexp) (clevel $zlist)

// The following takes time
* IPW estimates at selected quantiles using community-contributed command poparms
set seed 10101
poparms (clevel $zlist) (drugexp $xlist), ipw quantiles(.50) ///
    vce(bootstrap, reps(400))   

* TEs for the 25th quantile
margins i.clevel, pwcompare predict(equation(#2)) 

// Not included
* Efficient influence function estimates at mean using the community-contributed command poparms
poparms (clevel $zlist) (drugexp $xlist) 
contrast r.clevel, nowald

*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph

********** END
