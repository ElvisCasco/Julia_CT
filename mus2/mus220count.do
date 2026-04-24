* mus220count.do for Stata 17

cap log close

********** OVERVIEW OF mus220count.do **********
* Stata program
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi
* used for "Microeconometrics Using Stata, Second Edition"
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 20
*  20.2: MODELING STRATEGIES FOR COUNT DATA
*  20.3: POISSON AND NB MODELS
*  20.4: HURDLE MODELS
*  20.5: FINITE-MIXTURE MODELS
*  20.6: ZERO-INFLATED MODELS
*  20.7: ENDOGENOUS REGRESSORS
*  20.9: QUANTILE REGRESSION FOR COUNT DATA

* To run you need file
*   mus220mepsdocvis.dta
*   mus220mepsemergroom.dta
*   mus220meps2003qr.dta
* in your directory

* community-contributed commands
*   spost13_ado   
*   prvalue and prcounts  (in spost13_ado package)
*   countfit  (in spost9_legacy package) 
*   chi2gof
*   hnblogit
*   qcount
*   qplot 
* are used

********** SETUP ****************

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Files mus220mepsdocvis.dta, mus220mepsemergroom.dta, mus220meps2003qr.dta
* are extracts from MEPS (Medical Expenditure Panel Survey)
* for individuals 65 years and older in U.S. Medicare in 2003


********** CHAPTER 20.2: MODELING STRATEGIES FOR COUNT DATA

* Poisson (mu=1) generated data
qui set obs 10000
set seed 10101                // Set the seed
generate xpois = rpoisson(1)  // Draw from Poisson(mu=1)
qui histogram xpois, discrete xtitle("Poisson") saving(histpois.gph, replace)
summarize xpois
tabulate xpois

* Negative binomial (mu=1 var=2) generated data
set seed 10101                   // Set the seed
generate xg = rgamma(1,1)
generate xnegbin = rpoisson(xg)  // NB generated as a Poisson-gamma mixture
qui histogram xnegbin, discrete xtitle("Negative binomial") ///
    saving(histnb2.gph, replace)
summarize xnegbin
tabulate xnegbin

********* CHAPTER 20.3 POISSON AND NEGATIVE BINOMIAL MODELS

* Summary statistics for doctor-visits data
qui use mus220mepsdocvis, clear
global xlist private medicaid age age2 educyr actlim totchr
summarize docvis $xlist

* Tabulate docvis after recoding values > 10 to ranges 11-40 or 41-143
generate dvrange = docvis
recode dvrange (11/40 = 40) (41/143 = 143)
tabulate dvrange

* Poisson with default ML standard errors
poisson docvis $xlist, nolog

* Poisson: Squared correlation between y and yhat
predict yphat, n
qui correlate docvis yphat
display "Squared correlation between y and yhat = " r(rho)^2

* Poisson with robust standard errors
poisson docvis $xlist, vce(robust) nolog  // Poisson robust SEs

// Not included
* Poisson as glm
glm docvis $xlist, family(poisson) link(log) robust  // glm

* Overdispersion test against V(y|x) = E(y|x) + a*{E(y|x)^2}
qui poisson docvis $xlist, vce(robust)
predict muhat, n
qui generate ystar = ((docvis-muhat)^2 - docvis)/muhat
regress ystar muhat, noconstant noheader

* Pearson and deviance measures after Poisson
qui poisson docvis $xlist, nolog
estat gof

* Poisson AMEs
qui poisson docvis i.private i.medicaid c.age##c.age educyr i.actlim ///
    totchr, vce(robust)
margins, dydx(*)

* NB2: Standard negative binomial with robust standard errors
nbreg docvis $xlist, vce(robust) nolog

* NB2: Squared correlation between y and yhat
predict ynbhat, n
qui correlate docvis ynbhat
display "Squared correlation between y and yhat = " r(rho)^2

// Not included
* Compute Pr(y_i = 2) for each observation and on average for the sample.
poisson docvis $xlist, nolog
predict p2hat, pr(2)
sum p2hat
drop p2hat
margins, predict(pr(2))

* Poisson: Sample versus avg predicted probabilities of y = 0, 1, ..., 5
countfit docvis $xlist, maxcount(5) prm nograph noestimates nofit

* NB2: Sample versus average predicted probabilities of y = 0, 1, ..., 5
countfit docvis $xlist, maxcount(5) nbreg nograph noestimates nofit

* NB2: Formal chisquare goodness of fit test
qui nbreg docvis $xlist
chi2gof, cells(0, 1, 2, 3, 4, 5) table

// Not included
qui nbreg docvis $xlist, vce(robust)
margins, predict(pr(2)) at(private=1 medicaid=1)

* NB2: Cumulative fitted probabilities and ave fitted probs of y = 0 to max()
qui nbreg docvis $xlist, vce(robust)
prcounts erpr, max(3)
summarize erp*

* NB2: Predicted NB2 probabilities at x = x* of y = 0, 1, ..., 5
qui nbreg docvis $xlist, vce(robust)
prvalue, x(private=1 medicaid=1) max(5) brief

* Generalized negative binomial with alpha parameterized
gnbreg docvis $xlist, lnalpha(female bh) vce(robust) nolog

* NLS with exponential conditional mean
nl (docvis = exp({xb: $xlist} + {constant})), vce(robust) nolog

********* CHAPTER 20.4 HURDLE MODELS

* Hurdle: Pr(y=0)=pi and Pr(y=k)=(1-pi) x Poisson(2) truncated at 0
clear
qui set obs 10000
set seed 10101             // Set the seed
scalar pi=1-(1-exp(-2))/2  // Probability y=0
generate xhurdle = 0
scalar minx = 0
while minx == 0 {
    generate xph = rpoisson(2)
    qui replace xhurdle = xph if xhurdle==0
    drop xph
    qui summarize xhurdle
    scalar minx = r(min)
}
replace xhurdle = 0 if runiform() < pi
qui histogram xhurdle, discrete xtitle("Hurdle Poisson") ///
   saving(histphurdle, replace)
summarize xhurdle

* Hurdle logit-nb2 model manually: (1) Logit for zeros
qui use mus220mepsdocvis, clear
logit docvis $xlist, nolog vce(robust)

* Hurdle logit-nb2 model manually: (2a) restrict to positives only
summarize docvis if docvis > 0

* ztnb (old) replaced by tnbreg (current)
* Hurdle logit-nb2 model manually: (2b) ZTNB for positives
tnbreg docvis $xlist if docvis>0, nolog  vce(robust)

// Not included
* Same hurdle logit-nb2 model using the community-contributed hnblogit command
hnblogit docvis $xlist, nolog vce(robust)

* Hurdle logit-nb2: AMEs for first part
global xlistfv i.private i.medicaid c.age##c.age educyr i.actlim totchr
qui logit docvis $xlistfv, vce(robust)
margins, dydx(*)
	
* Hurdle logit-nb2: AMEs for second part
qui tnbreg docvis $xlistfv if docvis>0, vce(robust)
margins, dydx(*)

* Hurdle logit-poisson: AMEs for combined model
generate visit = (docvis > 0)
generate novisit = 1 - visit
qui gsem (docvis <- $xlistfv, family(poisson, ltruncated(novisit)))  ///
    (visit <- $xlistfv, logit), vce(robust)
margins, dydx(*) expression( (exp(predict(eta))/(1-exp(-exp(predict(eta))))) ///
    * predict(equation(visit)) )

********* CHAPTER 20.5 FINITE MIXTURE MODELS

* Finite Mixture Models using official fmm command in Stata 16

* Poisson (mu=1) and NB(1,1) generated data
qui set obs 10000
set seed 10101                   // Set the seed
generate xpois= rpoisson(1)
generate xg = rgamma(1,1)
generate xnegbin = rpoisson(xg)  // NB draws

* Finite Mixture: Poisson(.5) with prob .9 and Poisson(5.5) with prob .1
set seed 10101        // Set the seed
generate xp1= rpoisson(.5)
generate xp2= rpoisson(5.5)
summarize xp1 xp2
rename xp1 xpmix
qui replace xpmix = xp2 if runiform() > 0.9
qui histogram xpmix, discrete xtitle("Finite-mixture Poisson") ///
    saving(histp2mix, replace)
summarize xpmix

* Finite Mixture: Relative frequencies
tabulate xpmix

* Histograms of four different distributions, all with mean 1
graph combine histpois.gph histnb2.gph histphurdle.gph histp2mix.gph, ///
    iscale(1.1) ycommon xcommon          ///
    title("Four different distributions with mean = 1")


// Need to re-create xhurdle for next table
set seed 10101             // Set the seed
scalar pi=1-(1-exp(-2))/2  // Probability y=0
generate xhurdle = 0
scalar minx = 0
while minx == 0 {
    generate xph = rpoisson(2)
    qui replace xhurdle = xph if xhurdle==0
    drop xph
    qui summarize xhurdle
    scalar minx = r(min)
}
replace xhurdle = 0 if runiform() < pi

* Means and standard deviations of four different distributions, all with mean 1
summarize xpois xnegbin xhurdle xpmix

* FMM two-component Poisson with constant probabilities
qui use mus220mepsdocvis, clear
global xlist i.private i.medicaid c.age##c.age c.educyr i.actlim c.totchr
fmm 2, nolog vce(robust): poisson docvis $xlist

* FMM two-component Poisson: Predicted y's and histogram for each component
qui fmm 2, vce(robust): poisson docvis $xlist
predict yhatp*
summarize yhatp1 yhatp2
qui histogram yhatp1, name(class_1, replace)
qui histogram yhatp2, name(class_2, replace)
qui graph combine class_1 class_2, iscale(1.2) rows(1) ycommon xcommon ///
    ysize(2.5) xsize(6)


* FMM two-component Poisson: Compute the component marginal predicted means
estat lcmean

* FMM two-component Poisson: Compute the component mean predicted probabilities
estat lcprob

* FMM two-component Poisson: MEs
margins, dydx(*)

* FMM two-component negative binomial with constant probabilities
fmm 2, nolog vce(robust): nbreg docvis $xlist

* FMM two-component NB2: Component marginal predicted probabilities and means
estat lcprob
estat lcmean

* Posterior probabilities of each latent class for Poisson and NB2
qui fmm 2, vce(robust): poisson docvis $xlist
predict postpois*, classposteriorpr
qui fmm 2, vce(robust): nbreg docvis $xlist
predict postnb*, classposteriorpr
sum postpois* postnb*

* Kernel densities for posterior probabilities of each class for Poisson and NB2
kdensity postpois1, lwidth(medthick) title(" ") name(postpois1, replace)
kdensity postpois2, lwidth(medthick) title(" ") name(postpois2, replace)
kdensity postnb1, lwidth(medthick) title(" ") name(postnb1, replace)
kdensity postnb2, lwidth(medthick) title(" ") name(postnb2, replace)
graph combine postpois1 postpois2 postnb1 postnb2, rows(2) ycommon ///
    xcommon ysize(5) xsize(6) iscale(0.8)


* Use coeflegend to find complete names of model coefficients
qui fmm 2, vce(robust): poisson docvis $xlist
fmm, coeflegend

* Testing equality of totchr coefficient across the two components
test (_b[docvis:1.Class#c.totchr] = _b[docvis:2.Class#c.totchr])
contrast c.totchr#a.Class, equation(docvis)

* FMM two-component Poisson: lcprob() option to allow varying mixture proportions
fmm 2, lcprob(actlim totchr) nolog vce(robust): poisson docvis private ///
   medicaid age age2 educyr

* FMM two-component Poisson: Postestimation summary statistics
estat lcmean
estat lcprob
estimates stats

* FMM two-component Poisson: Different regressors in the two components
fmm, lcprob(actlim totchr) nolog vce(robust): ///
    (poisson docvis private age age2 educyr) (poisson docvis private age age2)

** CHAPTER  20.6: Zero-inflated data and models(Empirical Example 2)

********* CHAPTER 20.6  ZERO-INFLATED MODELS

* Summary statistics for emergency room visits
qui use mus220mepsemergroom, clear
global xlist1 age i.actlim totchr
summarize er $xlist1
tabulate er
* NB2 for emergency room visits
nbreg er $xlist1, nolog vce(robust)

* Zero-inflated NB2 for emergency room visits
zinb er $xlist1, inflate($xlist1) nolog vce(robust)

* ZINB: AMEs
qui zinb er $xlist1, inflate($xlist1) vce(robust)
margins, dydx(*)

* NB2: AMEs
qui nbreg er $xlist1, vce(robust)
margins, dydx(*)

* Zero-inflated NB2 estimated using fmm command with point-mass zero
fmm, nolog vce(robust): (pointmass er, lcprob($xlist1)) (nbreg er $xlist1)

* Zero-inflated NB2 estimated using fmm command: Summary statistics
predict mu*
summarize mu*
estat lcmean
estat lcprob

* NB2 model: Average fitted frequencies and information criteria
countfit er age actlim totchr, nbreg nograph noestimates nofit
qui nbreg er $xlist1
estat ic

* ZINB model: Average fitted frequencies and information criteria
qui zinb er $xlist1, inflate($xlist1) vce(robust)
prcounts erpr, max(9)
summarize erprpr*
estat ic

********* CHAPTER 20.7: ENDOGENOUS REGRESSORS

* Poisson control function estimator: First-stage linear regression
qui use mus220mepsdocvis, clear
global xlist2 medicaid age age2 educyr actlim totchr
regress private $xlist2 income ssiratio, vce(robust)
predict lpuhat, residual

* Poisson control function estimator: Second-stage poisson regression
poisson docvis private $xlist2 lpuhat, vce(robust) nolog

* Poisson control function estimator: Obtain bootstrap standard errors for estimator
program endogtwostep, eclass
    version 17
    tempname b
    capture drop lpuhat2
    regress private $xlist2 income ssiratio
    predict lpuhat2, residual
    poisson docvis private $xlist2 lpuhat2
    matrix `b' = e(b)
    ereturn post `b'
end
bootstrap _b, reps(400) seed(10101) nodots nowarn: endogtwostep

* Poisson ML estimator using GSEM command
qui gsem (docvis <- private $xlist2 L, poisson) ///
    (private <- $xlist2 income ssiratio L), nolog vce(robust)
estimates table, keep(private) b(%10.7f) se stats(N)

* ivpoisson gmm: (1)-(2) additive one and two step; (3)-(4) multiplicative one and two
qui ivpoisson gmm docvis $xlist2 (private=income ssiratio), onestep
estimates store GMMadd1S
qui ivpoisson gmm docvis $xlist2 (private=income ssiratio)
estimates store GMMadd2S
qui ivpoisson gmm docvis $xlist2 (private=income ssiratio), onestep multiplic
estimates store GMMmult1S
qui ivpoisson gmm docvis $xlist2 (private=income ssiratio), multiplicative
estimates store GMMmult2S

* Table for GMM additive one step and two step and multiplicative one step and two step
estimates table GMMadd1S GMMadd2S GMMmult1S GMMmult2S,     ///
    b(%7.4f) se(%7.3f) stats(N) stfmt(%9.1f) modelwidth(9)

* Test of overidentifying restriction following two-step GMM
qui ivpoisson gmm docvis $XLISTEXOG (private=income ssiratio), multiplicative

// Output not included
* gmm command for additive and multiplicative one step
gmm (docvis - exp({xb:private $xlist2 _cons})),        ///
    instruments(income ssiratio $xlist2) onestep vce(robust)
gmm ( (docvis / exp({xb:private $xlist2 _cons})) - 1), ///
    instruments(income ssiratio $xlist2) onestep vce(robust)

* ivpoisson cfunction estimation for endogenous Poisson
ivpoisson cfunction docvis  medicaid age age2 educyr actlim totchr   ///
        (private=income ssiratio), vce(robust)

********** 20.9: QUANTILE REGRESSION FOR COUNT DATA

* Summary statistics for doctor-visits data
qui use mus220meps2003qr, clear
summarize docvis private totchr age female white, separator(0)

* QR: Generate jittered values and compare quantile plots
set seed 10101
generate docvisu = docvis + runiform()
qui qplot docvis if docvis < 40, recast(line) lwidth(medthick)     ///
    ytitle("Quantiles of doctor visits") saving(graph1, replace)
qui qplot docvisu if docvis < 40, recast(line) lwidth(medthick)    ///
    ytitle("Quantiles of jittered doctor vists") saving(graph2, replace)
graph combine graph1.gph graph2.gph, iscale(1.3) rows(1)    ///
    ycommon xcommon ysize(2.5) xsize(6)

* QCR: MEs from conventional negative binomial model
qui nbreg docvis private totchr age female white, vce(robust)
margins, dydx(*) atmeans noatlegend

* QCR for the median
set seed 10101
qcount docvis private totchr age female white, q(0.50) rep(500)

* QCR: MEs after QCR for the median
qcount_mfx

* QCR: MEs after QCR for q = 0.75
set seed 10101
qui qcount docvis private totchr age female white, q(0.75) rep(500)
qcount_mfx

*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph
erase histpois.gph
erase histnb2.gph
erase histphurdle.gph
erase histp2mix.gph

********** END
