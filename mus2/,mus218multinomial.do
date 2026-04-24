* mus218multinomial.do  for Stata 17

capture log close

********** OVERVIEW OF mus218multinomial.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 18
* 18.3: MULTINOMIAL EXAMPLE: CHOICE OF FISHING MODEL
* 18.4: MNL MODEL
* 18.5: ALTERNATIVE-SPECIFIC CL MODEL
* 18.6: NL MODEL
* 18.7: MNP MODEL
* 18.8: ALTERNATIVE-SPECIFIC RANDOM-PARAMETERS LOGIT
* 18.9: ORDERED OUTCOME MODELS
* 18.11: MULTIVARIATE OUTCOMES

* To run you need files
*   mus218hk.dta 
*   mus218rhie.dta 
* in your directory

* No Stata addons are used

********** DATA DESCRIPTION **********

* mus218hk.dta is from J. A. Herriges and C. L. Kling, 
* "Nonlinear Income Effects in Random Utility Models", 
* Review of Economics and Statistics, 81(1999): 62-72
* Also analyzed in Cameron and Trivedi (2005) chapter 15

* mus218rhie.dta is from Rand Health Insurance Experiment data 
* Essentially same data as in P. Deb and P.K. Trivedi (2002)
* "The Structure of Demand for Medical Care: Latent Class versus
* Two-Part Models", Journal of Health Economics, 21, 601-625
* except that article used different outcome (counts rather than $)
* Each observation is for an individual over a year.
* Individuals may appear in up to five years.
* All available sample is used except only fee for service plans included.
* If panel data used then clustering is on id (person id)

********** SETUP **********

clear all
set linesize 80
set scheme s1mono  /* Graphics scheme */

********** 18.3: MULTINOMIAL EXAMPLE: CHOICE OF FISHING MODEL

* Read in dataset, and describe dependent variable and regressors
qui use mus218hk
describe

* Summarize dependent variable and regressors
summarize, separator(0) 

* Tabulate the dependent variable
tabulate mode

* Table of income by fishing mode
table mode, stat(count income) stat(mean income) stat(sd income)

* Table of fishing price by fishing mode
table (result) mode, stat(mean pbeach ppier pprivate pcharter) nformat(%6.0f)

* Table of fishing catch rate by fishing mode
table (result) mode, stat(mean qbeach qpier qprivate qcharter) nformat(%6.2f)


********** 18.4 MULTINOMIAL LOGIT MODEL

* Multinomial logit with base outcome alternative 1
mlogit mode income, baseoutcome(1) nolog vce(robust)

* Wald test of the joint significance of income
test income

* Relative-risk option reports exp(b) rather than b
mlogit mode income, rrr baseoutcome(1) nolog vce(robust)
estimates store MNL


* Predict probabilities of choice of each mode, and compare to actual freqs
predict pmlogit1 pmlogit2 pmlogit3 pmlogit4, pr
summarize pmlogit* dbeach dpier dprivate dcharter, separator(4)

* Sample average predicted probability of the third outcome
margins, predict(outcome(3)) noatlegend

* AME of income change for outcome 3
margins, dydx(*) predict(outcome(3)) noatlegend

* MEM of income change for outcome 3
margins, dydx(*) predict(outcome(3)) atmeans noatlegend

********** 18.5 ALTERNATIVE-SPECIFIC CONDITIONAL LOGIT MODEL

* Data are in wide form
list mode price pbeach ppier pprivate pcharter in 1, clean

* Convert data from wide form to long form
generate id = _n
reshape long d p q, i(id) j(fishmode beach pier private charter) string
save mus218hklong, replace

* List data for the first case after reshape
list in 1/4, clean noobs

* cmset before use of cm commands for alternative-specific regressors 
cmset id fishmode

* Summarize the data by chosen alternative specific
cmsummarize p q income, choice(d) statistics(N mean)

* Conditional logit with alternative-specific and case-specific regressors
cmclogit d p q, casevars(income) basealternative(beach) nolog vce(robust) 
estimates store CL

* Predicted probabilities of choice of each mode and compare to actual freqs
predict pasclogit, pr
encode fishmode, generate(fishmode2)
table fishmode2, stat(mean d pasclogit) stat(sd pasclogit) nototal nformat(%6.4f)


* AME of change in price
margins, dydx(p)

// Not included - AME on private of change in price
margins, dydx(p) alternative(private)

// Not included - MNL is CL with no alternative-specific regressors
cmclogit d, casevars(income) basealternative(beach) nolog vce(robust) 

********** 18.6 NESTED LOGIT MODEl

* Define the Tree for Nested logit
*       with nesting structure 
*             /     \
*           /  \   /  \

* Convert string variable fishmode to integer values and attach labels
* encode fishmode, gen(intmode)
* label list intmode     // Check ordering of intmode

* nlogitgen type = intmode(shore: 1 | 3, boat: 2 | 4)
* nlogitgen type = fishmode(shore: pier | private, boat: beach | charter)
* Define the tree for nested logit
nlogitgen type = fishmode(shore: pier | beach, boat: private | charter)

* Check the tree
nlogittree fishmode type, choice(d)

* Nested logit model estimate
nlogit d p q || type:, base(shore) || fishmode: income, case(id) notree ///
    vce(robust) nolog  
estimates store NL

* Predict level 1 and level 2 probabilities from NL model
predict plevel1 plevel2, pr
tabulate fishmode, summarize(plevel2)

list d plevel1 plevel2 in 1/4, clean
estat alternatives

* AME of beach price change computed manually
preserve
qui summarize p
generate delta = r(sd)/1000
qui replace p = p + delta if fishmode == "beach"
predict pnew1 pnew2, pr
generate dpdbeach = (pnew2 - plevel2)/delta
tabulate fishmode, summarize(dpdbeach)
restore

* Summary statistics for the logit models
estimates table MNL CL NL, keep(p q) stats(N ll aic bic) equation(1) ///
    b(%7.3f) stfmt(%7.0f)

*********** 18.7 MULTINOMIAL PROBIT MODEl 

// Output from following command was not included in the book
qui use mus218hk, clear
mprobit mode income, baseoutcome(1) vce(robust)

// Following with four choices did not converge or
// did not give SEs for variance parameters
// qui use mus218hklong, clear
// cmset id fishmode
// cmmprobit d p q, corr(exchangeable) basealt(charter) shownrtolerance 

* Multinomial probit with unstructured errors when charter is dropped  
qui use mus218hklong, clear
cmset id fishmode
drop if fishmode == "charter" | mode == 4
cmmprobit d p q, casevars(income) correlation(unstructured) structural ///
    vce(robust) nolog

* Show correlations and covariance
estat correlation
estat covariance 


// Not included
* Same estimator without option structural
cmmprobit d p q, casevars(income) correlation(unstructured) vce(robust) nolog

// Not included
* Compare to conditional logit with three choice
cmclogit d p q, casevars(income) basealternative(pier) nolog vce(robust) 

*********** 18.8 ALTERNATIVE-SPECIFIC RANDOM-PARAMETERS LOGIT

* Alternative-specific mixed logit or random parameters logit estimation
qui use mus218hklong, clear
drop if fishmode == "charter" | mode == 4
cmset id fishmode    // caseidvar casevars
cmmixlogit d q, casevars(income) random(p) basealternative(pier) ///
    vce(robust) nolog

* AMEs with respect to price
margins, dydx(p)

// Not included - default standard errors
cmmixlogit d q, casevars(income) random(p) basealternative(pier) nolog
margins, dydx(p)

// Not included - same example estimated as simpler conditional logit model
asclogit d q p, case(id)  alternatives(fishmode) casevars(income) ///
     basealternative(pier) vce(robust) nolog

/* The following requires the community-contributed program mixlogit
   This was used in the first and revised editions. Now use cmmixlogit 
* Data set up to include case-invariant regressors
qui use mus218hklong, clear
generate dbeach = fishmode=="beach"
generate dprivate = fishmode=="private"
generate dcharter = fishmode=="charter"
generate ybeach = dbeach*income
generate yprivate = dprivate*income
generate ycharter = dcharter*income
drop if fishmode=="charter" | mode == 4
mixlogit d q dbeach dprivate ybeach yprivate, group(id) rand(p)
*/

********** 18.9: ORDERED OUTCOME MODELS

* Create multinomial ordered outcome variables that take values y = 1, 2, 3
qui use mus218rhie, clear
qui keep if year == 2
generate hlthpf = hlthp + hlthf
generate hlthe = (1 - hlthpf - hlthg)
qui generate hlthstat = 1 if hlthpf == 1
qui replace hlthstat = 2 if hlthg == 1
qui replace hlthstat = 3 if hlthe == 1
label variable hlthstat "Health status"
label define hsvalue 1 "Poor or fair" 2 "Good" 3 "Excellent"
label values hlthstat hsvalue
tabulate hlthstat

* Summarize dependent and explanatory variables
summarize hlthstat age linc ndisease

* Ordered logit estimates
ologit hlthstat age linc ndisease, nolog vce(robust)

* Calculate predicted probability that y=1, 2, or 3 for each person
predict p1ologit p2ologit p3ologit, pr
summarize hlthpf hlthg hlthe p1ologit p2ologit p3ologit, separator(0)

* MEM for third outcome (health status excellent)
margins, dydx(*) predict(outcome(3)) atmeans noatlegend

// Not -included - ordered probit and heteroskedastic probit
oprobit hlthstat age linc ndisease, nolog vce(robust)
hetoprobit hlthstat age linc ndisease, het(age linc ndisease) vce(robust)

********** 18.11: MULTIVARIATE OUTCOMES

* Two binary dependent variables: hlthe and dmdu
tabulate hlthe dmdu
correlate hlthe dmdu

* Bivariate probit estimates
biprobit hlthe dmdu age linc ndisease, nolog vce(robust)

* Predicted probabilities
predict biprob1, pmarg1
predict biprob2, pmarg2
predict biprob11, p11
predict biprob10, p10
predict biprob01, p01
predict biprob00, p00
summarize hlthe dmdu biprob1 biprob2 biprob11 biprob10 biprob01 biprob00

* Separate probits
probit hlthe age linc ndisease, vce(robust)
probit dmdu age linc ndisease, vce(robust)

* Bivariate probit with different sets of regressors
biprobit (hlthe = age linc) (dmdu = age linc ndisease), nolog vce(robust)

* Nonlinear seemingly unrelated regressions estimator
nlsur (hlthe = normal({a1}*age+{a2}*linc+{a3}*ndisease+{a4})) ///
    (dmdu = normal({b1}*age+{b2}*linc+{b3}*ndisease+{b4})), vce(robust) nolog

*** Erase files created by this program and not used elsewhere in the book

erase mus218hklong.dta

********** END
