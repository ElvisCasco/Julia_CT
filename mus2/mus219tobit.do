* mus219tobit.do for Stata 17

cap log close

********** OVERVIEW OF mus219tobit.do **********
* Stata program
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi
* used for "Microeconometrics Using Stata, Second Edition"
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

*********************************************

* CHAPTER 19 TOBIT AND SELCTION MODELS
*   19.3 TOBIT MODEL EXAMPLE
*   19.4 TOBIT FOR LOGNORMAL DATA
*   19.5 TWO-PART MODEL IN LOGS
*   19.6 SELECTION MODEL
*   19.7 NONNORMAL MODELS OF SELECTION
*   19.8 PREDICTION FROM MODELS WITH OUTCOME IN LOGS
*   19.11 PANEL ATTRITION

* To use this program you need files
*   mus219mepsambexp.dta
*   mus219mabelunbalsmall.dta

* And you need community-contributed commands
*   heckmancopula
*   xtbalance

********** SETUP **********

clear all
set linesize 82
set scheme s1mono   /* Used for graphs */

********** DATA DESCRIPTION **********

* mus219mepsambexp.dta is a subset of data in
* P. Deb, M. Munkin and P.K. Trivedi (2006)
* "Bayesian Analysis of Two-Part Model with Endogeneity",
* Journal of Applied Econometrics, 21, 1081-1100
* Only the data for year 2001 are used

* mus219mabelbal.dta and mus219mabelunbalsmall.dta
* are extracts from MABEL EPS
* Medicine in Australia: Balancing Employment and Life (MABEL)
* An Australian longitudinal survey of physicians.

********** READ DATA **********

* Raw data summary
qui use mus219mepsambexp, clear
summarize ambexp age female educ blhisp totchr ins

******** CHAPTER 19.3 TOBIT MODEL EXAMPLE

* Detailed summary to show skewness and kurtosis
summarize ambexp, detail

* Summary for positives only
summarize ambexp if ambexp >0, detail

* Tobit analysis for ambexp using all expenditures
global xlist age female educ blhisp totchr ins  // Define regressor list $xlist
tobit ambexp $xlist, ll(0) vce(robust)

* Tobit prediction and summary
predict yhatlin
summarize yhatlin

* (1) MEMs for E(y|x,y>0)
qui tobit ambexp age i.female educ i.blhisp totchr i.ins, ll(0) vce(robust)
margins, dydx(*) predict(e(0,.)) atmeans noatlegend

* (2) MEMs for E(y|x)
margins, dydx(*) predict(ystar(0,.)) atmeans noatlegend

* Direct computation of MEMs for E(y|x)
predict xb1, xb                // xb1 is estimate of x'b
matrix btobit = e(b)
scalar sigmasq = btobit[1,11]  // sigmasq is estimate of sigma^2
matrix bcoeff = btobit[1,1..9] // bcoeff is betas excl. constant
qui summarize xb1
scalar meanxb = r(mean)        // Mean of x'b equals (mean of x)'b
scalar PHI = normal(meanxb / sqrt(sigmasq))
matrix deriv = PHI*bcoeff
matrix list deriv
ereturn post deriv             // Print nicer looking results
ereturn display

* Compute MEMs for Pr(5000<ambexp<10000)
qui tobit ambexp age i.female educ i.blhisp totchr i.ins, ll(0) vce(robust)
margins, dydx(*) predict(pr(5000,10000)) atmeans noatlegend

* Truncated tobit with truncation at zero
truncreg ambexp age female educ blhisp totchr ins, ll(0) nolog vce(robust)

* MEMs for E(y|x) following truncated tobit
margins, dydx(*) predict(e(0,.)) atmeans noatlegend

// Included as command without output
* gsem alternative for estimating tobit regression
gsem ambexp <- $xlist, family(gaussian, lcensored(0)) vce(robust)
tobit ambexp $xlist, ll(0) nolog vce(robust)

********** CHAPTER 19.4 TOBIT FOR LOGNORMAL DATA

* ASIDE: IN BELOW ....
* We observe  y = 0 if dy=0   and  y = y  if dy=1
* Let gamma be the smallest of the positive y's
* We define lny = gamma-0.000001 if dy=0   and lny = ln(y) if dy=1
* This ensures that tobit, ll will correctly sort out
* the censored and uncensored observations
* Note that can't set lny = 0 if dy=1 as ln(y) < 0 if gamma
* drop xb

* Summary of log(expenditures) for positives only
summarize lambexp, detail

/* Following code included in text explains how y, dy and lny created
* "Tricking" Stata to handle log transformation and labelling variables
generate y = ambexp
generate dy = ambexp > 0
generate lny = ln(y)                // Zero values will become missing
qui summarize lny
scalar gamma = r(min) - 0.0000001   // This could be negative
display "gamma = " gamma            // gamma was 0 for these data
replace lny = gamma  if lny == .
tabulate y if y < 0.02              // 0.02 is arbitrary small value
tabulate lny if lny < gamma + 0.02
*/

* Code below needs gamma which is the lower limit of lny (here 0) - 0.0000001
scalar gamma = -0.0000001

// Not included
summarize y lny dy $xlist
summarize y lny dy, detail

* Now do tobit on lny
tobit lny $xlist, ll(gamma) vce(robust)

* OLS, not tobit
regress lny $xlist, noheader vce(robust)

* Now do two-limit tobit
scalar upper = log(10000)
display upper
tobit lny $xlist, ll(gamma) ul(9.2103404) vce(robust)

* Compute Mills's ratio
qui tobit lny $xlist, ll(gamma) vce(robust)
predict xb, xb                            // xb is estimate of x'b
matrix btobit = e(b)
scalar sigma = sqrt(btobit[1,e(df_m)+2])  // sigma is estimate of sigma
generate threshold = (gamma-xb)/sigma     // gamma: lower censoring point
generate lambda = normalden(threshold)/normal(threshold)

* Generalized residuals: First create gres1, which has mean zero, by the f.o.c.
qui generate uifdyeq1 = (lny - xb)/sigma if dy == 1
qui generate double gres1 = uifdyeq1
qui replace gres1 = -lambda if dy == 0
summarize gres1

* Generalized residuals: Next create gres2, gres3, and gres4
qui generate double gres2 = uifdyeq1^2 - 1
qui replace gres2 = -threshold*lambda if dy == 0
qui generate double gres3 = uifdyeq1^3
qui replace gres3 = -(2 + threshold^2)*lambda if dy == 0
qui generate double gres4 = uifdyeq1^4 - 3
qui replace gres4 = -(3*threshold + threshold^3)*lambda if dy == 0

* Generate the scores to use in the Lagrange multiplier test
* These are the actual scores up to a multiple constant over i
* So get same Lagrange multiplier test as would if added the multiple
* The scores wrt b are gres1 times the relevant component of x

* Generate the scores to use in the Lagrange multiplier test of normality
foreach var in $xlist {
   generate score`var' = gres1*`var'
   }
global scores score* gres1 gres2

* The score wrt sigma can be shown to be gres2
* So the scores are the score`var', gres1 (for the intercept) and gres2
* Test that calculations done correct
* All the score components should have a mean of zero
* gres3 and gres4 may not if model misspecified and these are the basis of the test
summarize $scores gres3 gres4 threshold lambda

* Test of normality in tobit regression uses NR^2 from uncentered regression
generate one = 1
qui regress one gres3 gres4 $scores, noconstant
display "N R^2 = " e(N)*e(r2) " with p-value = " chi2tail(2,e(N)*e(r2))

* Test of homoskedasticity in tobit regression with w=x (aside from intercept)
foreach var in $xlist {
   generate gres2by`var' = gres2*`var'
}
qui regress one gres2by* $scores, noconstant
display "N R^2 = " e(N)*e(r2) " with p-value = " chi2tail(6,e(N)*e(r2))

********** CHAPTER 19.5 TWO-PART MODEL IN LOGS

* Part 1 of the two-part model
probit dy $xlist, nolog vce(robust)
scalar llprobit = e(ll)

* Part 2 of the two-part model
regress lny $xlist if dy==1, vce(robust)
scalar lllognormal = e(ll)
predict rlambexp, residuals

* Create two-part model log likelihood
scalar lltwopart = llprobit + lllognormal  //two-part model log likelihood
display "lltwopart = " lltwopart

* hettest and sktest commands require default standard errors
qui regress lny $xlist if dy==1
hettest
sktest rlambexp

********** CHAPTER 19.6 SELECTION MODELS

* Heckman MLE without exclusion restrictions
heckman lny $xlist, select(dy = $xlist) nolog vce(robust)

* Heckman two-step without exclusion restrictions
heckman lny $xlist, select(dy = $xlist) twostep

* Heckman MLE with exclusion restriction
heckman lny $xlist, select(dy = $xlist income) nolog vce(robust)

* Heckman MLE with exclusion restriction
* eregress lny $xlist, select(dy = $xlist income) nolog vce(robust)

// Not included
* Heckman two-step with exclusion restrictions
heckman lny $xlist, select(dy = $xlist income) twostep

********** CHAPTER 19.7  NON-NORMAL MODELS OF SELECTION

* Copula-based selection models with Gaussian copula (same as heckman)
qui use mus219mepsambexp, clear
global xlist age female educ blhisp totchr ins  // Regressor list $xlist
heckmancopula lny $xlist, select(dy = $xlist income) copula(gaussian) ///
   margsel(probit) margin1(normal) vce(robust)
estimates store Gaussian_N

* Copula-based selection models with several different copulas
qui heckmancopula lny $xlist, select(dy = $xlist income) copula(frank) ///
    margsel(probit) margin1(normal) vce(robust)
estimates store Frank_N
qui heckmancopula lny $xlist, select(dy = $xlist income) ///
    copula(fgm) margsel(probit) margin1(normal) vce(robust)
estimates store FGM_N
qui heckmancopula lny $xlist, select(dy = $xlist income) ///
    copula(plackett) margsel(probit) margin1(t) df(10) vce(robust)
estimates store Plackett_t10
estimates table Gaussian_N Frank_N FGM_N Plackett_t10, ///
    eq(1) se b(%13.4f) stats(N ll)

// Not included
heckmancopula lny $xlist, select(dy = $xlist income) copula(joe) ///
    margsel(probit) margin1(t) vce(robust)

********** CHAPTER 19.8 PREDICTION FROM MODELS WITH OUTCOME IN LOGS

* Prediction from tobit on lny
* The result is simple but derivation messy. Done in Carson and Sun.
*    E[y] = exp(xb+0.5*shat^2)(1 - PHI(((gamma-xb-shat^2)/shat)])
*    E[y | dy=1] = E[y] / Pr[ dy=1] = E[y] / (1-PHI(((xb-gamma)/shat)))
* drop yhat ytrunchat

* Prediction from tobit regression with dependent variable lny
qui use mus219mepsambexp, clear
scalar gamma = -0.0000001
qui tobit lny $xlist, ll(gamma) vce(robust)
predict xb, xb                            // xb is estimate of x'b
matrix btobit = e(b)
scalar sigma = sqrt(btobit[1,e(df_m)+2])  // sigma is estimate of sigma
generate threshold = (gamma-xb)/sigma     // gamma: lower censoring point
generate yhat = exp(xb+0.5*sigma^2)*(1-normal((gamma-xb-sigma^2)/sigma))
generate ytrunchat = yhat / (1 - normal(threshold)) if dy==1
summarize y yhat
summarize y yhat ytrunchat if dy==1

* Two-part model predictions
qui probit dy $xlist
predict dyhat, pr
qui regress lny $xlist if dy==1
predict xbpos, xb
generate yhatpos = exp(xbpos+0.5*e(rmse)^2)

* Unconditional prediction from two-part model
generate yhat2step = dyhat*yhatpos
summarize yhat2step y
summarize yhatpos y if dy==1

* Heckman model predictions
qui heckman lny $xlist, select(dy = $xlist) vce(robust)
predict probpos, psel
predict x1b1, xbsel
predict x2b2, xb
scalar sig2sq = e(sigma)^2
scalar sig12sq = e(rho)*e(sigma)^2
display "sigma1sq = 1" " sigma12sq = " sig12sq " sigma2sq = " sig2sq
generate yhatheck = exp(x2b2 + 0.5*(sig2sq))*(1 - normal(-x1b1-sig12sq))
generate yhatposheck = yhatheck/probpos
summarize yhatheck y probpos dy
summarize yhatposheck probpos dy y if dy==1

// Not included
* THE FOLLOWING SHOWS THAT THE POOR PREDICTIONS ARE DUE TO HIGH SIGMA
* THE ABOVE MODEL HAD SIGMA = 2.378
* BY CONTRAST IF WE JUST REGRESS LNY on $XLIST WITH POSITIVE OBS THEN SIGMA = 1.264
scalar sigma = 1
generate yhatsigma1 = exp(xb+0.5*sigma^2)*(1-normal((gamma-xb-sigma^2)/sigma))
scalar sigma = 2
generate yhatsigma2 = exp(xb+0.5*sigma^2)*(1-normal((gamma-xb-sigma^2)/sigma))
summarize yhat*
summarize yhat* if dy==1

********** CHAPTER 19.11 : PANEL ATTRITION

* Pattern of missing observations in panel dataset
qui use mus219mabelunbalsmall, clear
drop if lyearn==.|lyhrs==.|expr==.|exprsq==.|hospwork==.|selfempl==.
xtset
xtdescribe

* Summarize the data
summarize, sep(0)

* Selection indicators for each(id,wave) indicate which waves are available
generate s1 = 0
qui replace s1 = 1 if wave==1
qui replace s1 = 1 if (wave==2 & !missing(L.lyearn))
qui replace s1 = 1 if (wave==3 & !missing(L2.lyearn))
generate s2 = 0
qui replace s2 = 1 if wave==2
qui replace s2 = 1 if (wave==1 & !missing(F.lyearn))
qui replace s2 = 1 if (wave==3 & !missing(L.lyearn))
generate s3 = 0
qui replace s3 = 1 if wave==3
qui replace s3 = 1 if (wave==1 & !missing(F2.lyearn))
qui replace s3 = 1 if (wave==2 & !missing(F.lyearn))
generate balanced = (s1==1 & s2==1 & s3==1)
save mus219mabelfinal, replace

* Count of how many wave 2 also available in wave 1 and wave 3
tabulate s2 wave

* Compare wave 1 earnings by wave 3 attrition status
capture drop attrit*
twoway (kdensity lyearn if (wave==1 & s3==1), clstyle(p1))   ///
    (kdensity lyearn if (wave==1 & s3==0), clstyle(p2)),     ///
    title("Wave 1 earnings by wave 3 attrition status")      ///
    legend(label(1 "No attrition") label(2 "Attrition"))     ///
    legend(pos(10) ring(0) col(1)) saving(graph1.gph,replace)

// Not included
tabulate s3 if wave==1

* Compare wave 3 earnings by whether missed wave 2
generate wave2skip = .
qui replace wave2skip = 1 if (s1==1 & s2==0)
qui replace wave2skip = 0 if (s1==1 & s2==1)
twoway (kdensity lyearn if (wave==3 & wave2skip==0), clstyle(p1))  ///
    (kdensity lyearn if (wave==1 & wave2skip==1), clstyle(p2)),    ///
    title("Wave 3 earnings by missing wave 2")                     ///
    legend(label(1 "No attrition") label(2 "Attrition"))           ///
    legend(pos(10) ring(0) col(1)) saving(graph2.gph,replace)

* Combine the graphs
graph combine graph1.gph graph2.gph, scale(1.4) ysize(2.5) xsize(6)

// Not included
tabulate wave2skip if wave==3

* Difference in means test of wave 3 earnings by whether missed wave 2
ttest lyearn if wave==3, by(wave2skip) unequal

// Not included
regress lyearn wave2skip if wave==3, vce(robust)

* Manually create balanced dataset
keep if balanced == 1
xtdescribe

* xtbalance command creates a balenced dataset
qui use mus219mabelunbalsmall, clear
global xlist lyhrs expr exprsq hospwork
xtbalance, range(1 3) miss(lyearn $xlist)
xtdescribe

**** Analysis assuming missing at random

* OLS, random-effects, and fixed-effects regressions for unbalanced sample
* and for balanced sample
qui use mus219mabelfinal, clear   // Unbalanced sample
* drop if lyearn==.|lyhrs==.|expr==.|exprsq==.|hospwork==.
qui regress lyearn $xlist, vce(cluster id)
estimates store OLS_Unbal
qui xtreg lyearn $xlist, re vce(robust)
estimates store RE_Unbal
qui xtreg lyearn $xlist, fe vce(robust)
estimates store FE_Unbal
qui regress lyearn $xlist if balanced==1, vce(cluster id) // Balanced sample
estimates store OLS_Bal
qui xtreg lyearn $xlist if balanced==1, re vce(robust)
estimates store RE_Bal
qui xtreg lyearn $xlist if balanced==1, fe vce(robust)
estimates store FE_Bal

* Table comparing OLS, random-effects, and fixed-effects results for balanced
* and unbalanced data
estimates table OLS_Unbal OLS_Bal RE_Unbal RE_Bal FE_Unbal FE_Bal, ///
   b(%7.3f) se stats(N r2)

* Compare first-difference results for balanced and unbalanced data
global FDxlist D.lyhrs D.expr D.exprsq D.hospwork
sort id wave
qui regress D.lyearn $FDxlist, cluster(id)
estimates store FD_Unbal
qui regress D.lyearn $FDxlist if balanced==1, cluster(id) // Balanced sample
estimates store FD_Bal
estimates table FD_Unbal FD_Bal, b(%7.4f) se stats(N r2)

***** IPW

* Expand dataset to include missing waves 2 and 3 given present in wave 1
qui use mus219mabelfinal, clear
fillin id wave
sort id wave
list id wave s1 s2 s3 lyearn _fillin in 13/24, clean
generate swave = s1
replace swave = . if s1==0
drop s1 s2 s3
generate s1 = 0
qui replace s1 = 1 if (wave==1 & !missing(lyearn))
qui replace s1 = 1 if (wave==2 & !missing(L.lyearn))
qui replace s1 = 1 if (wave==3 & !missing(L2.lyearn))
generate s2 = 0
qui replace s2 = 1 if (wave==2 & !missing(lyearn))
qui replace s2 = 1 if (wave==1 & !missing(F.lyearn))
qui replace s2 = 1 if (wave==3 & !missing(L.lyearn))
generate s3 = 0
qui replace s3 = 1 if (wave==3 & !missing(lyearn))
qui replace s3 = 1 if (wave==1 & !missing(F2.lyearn))
qui replace s3 = 1 if (wave==2 & !missing(F.lyearn))
sum id wave s1 s2 s3 lyearn swave, sep(0)

* IPW: (1) Probit selection equation for wave 2 given in wave 1
global IPWxlist L.lyhrs L.expr L.exprsq L.hospwork
probit s2 $IPWxlist if wave==2
qui predict s2prob if e(sample)==1
tabulate s2 if e(sample)==1
tabstat s2prob if e(sample)==1, stats(count min p10 p50 p90 max)

// Not included
twoway (kdensity s2prob if e(sample)==1 & s2==1)  ///
    (kdensity s2prob if e(sample)==1 & s2==0),    ///
     legend(label(1 "s2=1") label(2 "s2=0"))

* IPW: (2) Generate weights and do weighted and unweighted first-difference OLS
qui generate s2weight = s2/s2prob if wave==2
tabstat s2weight if wave==2 & s2==1, stats(count min p10 p50 p90 max)
regress D.lyearn $FDxlist [pweight=s2weight] if wave==2, vce(robust) // Weighted
estimates store IPW2

// Not included - Manual check
preserve
generate manweight = s2/s2prob
generate wdlyearn = sqrt(manweight)*D.lyearn
generate wdlyhrs = sqrt(manweight)*D.lyhrs
generate wdlexpr = sqrt(manweight)*D.expr
generate wdlexprsq = sqrt(manweight)*D.exprsq
generate wdhospwork = sqrt(manweight)*D.hospwork
generate wdone = sqrt(manweight)
regress wdlyearn wdlyhrs wdlexpr wdlexprsq wdhospwork wdone if wave==2, ///
    noconstant vce(robust)
regress wdlyearn wdlyhrs wdlexpr wdlexprsq wdhospwork wdone if wave==2, ///
    noconstant vce(robust)
restore

* IPW: Now do wave 3 with selection on wave 1 variables
global IPWxlist3 L2.lyhrs L2.expr L2.exprsq L2.hospwork
qui probit s3 $IPWxlist3 if wave==3 & s2==1
qui predict s3prob if e(sample)
tabulate s3 if e(sample)==1
qui generate s3weight = s3/s3prob if wave==3
tabstat s3prob s3weight if e(sample)==1, stats(count min p10 p50 p90 max) col(s)
qui regress D.lyearn $FDxlist [pweight=s3weight] if wave==3, vce(robust)
estimates store IPW3

* IPW: Estimate waves 2 and waves 3 together
qui generate sweight = .
qui replace sweight = s2weight if wave==2
qui replace sweight = s3weight if wave==3
tabstat sweight if sweight!=0, stats(count min p10 p50 p90 max) col(s)
qui regress D.lyearn $FDxlist [pweight=sweight], vce(robust)
estimates store IPWboth

// Not included - time varying z for the selection probability
* IPW - wave 3 with selection on wave 2 variables
probit s3 $IPWxlist if wave==3 & s2==1
predict s3probnew if e(sample)
tabulate s3 if e(sample)==1
tabstat s3probnew if e(sample)==1, stats(count min p10 p50 p90 max) col(s)
generate s3weightnew = s3/(L.s2prob*s3probnew) if wave==3
regress D.lyearn $FDxlist [pweight=s3weightnew] if wave==3, vce(robust)
tabstat s3weightnew if e(sample)==1, stats(count min p10 p50 p90 max) col(s)

* IPW: Alternative weights based on auxiliary variables wave 2 versus wave 1
qui probit s2 $IPWxlist L.selfempl if wave==2
qui predict s2probaugment if e(sample)
qui generate s2relweight=s2prob/s2probaugment
tabstat s2relweight if e(sample)==1, stats(count min p10 p50 p90 max) col(s)
sum s2prob s2probaugment s2relweight
qui regress D.lyearn $FDxlist [pweight=s2relweight] if wave==2, vce(robust)
qui estimates store IPW2Aug

* Unweighted, IPW, and sample-selection estimates: First-difference in wave 2
qui regress D.lyearn $FDxlist if wave==2, vce(robust) noheader   // Unweighted
estimates store FD_UNW
estimates table FD_UNW IPW2 IPW3 IPWboth IPW2Aug, b(%7.4f) eq(1) se stats(N r2)

***** Selection model correction

* Sample-selection model: Pooled first-difference OLS using wave 2 and 3 data
generate dsel = .
replace dsel = s2 if wave==2
replace dsel = s3 if wave==3
tabulate dsel wave, missing
heckman D.lyearn $FDxlist, select(dsel = $IPWxlist L.selfempl) ///
    nolog vce(cluster id)
estimates store Heckpool

* Sample-selection model: First-difference OLS for wave 2 given wave 1
qui heckman D.lyearn $FDxlist if wave==2, select(s2 = $IPWxlist L.selfempl) ///
    mills(mills2) vce(robust)
display "N = " e(N) "  Independence test chi2(1) = " e(chi2_c) "  p = " e(p_c)
estimates store Heck2

* Sample-selection model: First-difference OLS for wave 3 given wave 1 and wave 2
qui heckman D.lyearn $FDxlist if wave==3, select(s3 = $IPWxlist L.selfempl) ///
    mills(mills3) vce(robust)
display "N = " e(N) "  Independence test chi2(1) = " e(chi2_c) "  p = " e(p_c)
estimates store Heck3

* Sample-selection model: Pooled first-difference OLS with different lambda for waves 2 and 3
generate millsp = mills2
replace millsp = mills3 if wave==3
regress D.lyearn $FDxlist i.wave#c.millsp, vce(cluster id)
test 2b.wave#c.millsp 3.wave#c.millsp
estimates store Heckboth

* Unweighted, IPW and sample-selection estiamtes: First-difference in wave 2
qui regress D.lyearn $FDxlist, vce(robust) noheader
estimates store FD_UNW
estimates table FD_UNW Heckpool Heck2 Heck3 Heckboth, b(%7.4f) eq(1) ///
   se stats(N r2) keep($FDxlist)


*** Erase files created by this program and not used elsewhere in the book

erase mus219mabelfinal.dta

*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph

********** END
