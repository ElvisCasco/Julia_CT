* canada2019_crosssection.do for Stata 15
* based on trcountimf2010.do June 2015 for Stata version 13
* and after that based on count2015sweden.do 

capture log close
log using canada2019_crosssection.txt, text replace

********** OVERVIEW OF canada2019_crosssection.do **********

* STATA Program to demonstrate Count Data Regression Models
* Based on mus17p1cnt.do in Cameron and Trivedi(2009)
* "Microeconometrics using Stata", Stata Press

* 1A. COUNT REGRESSION: BASIC CROSS-SECTION
* 1B. COUNT REGRESSION: INFERENCE (STANDARD ERRORS AND BOOTSTRAP) 
* 2: COUNT REGRESSION: EXTRAS (TRUNCATED/CENSORED, HURDLE, ZERO-INFLATED, 
*                              FINITE MIXTURE, ENDOGENOUS)

* To run you need files
*      mus17data.dta
*      mus18data.dta    
*      musbootdata.dta
* in your directory

* Stata user-written commands
*   chi2gof        // For chis-squared goodness of fit test
*   boottest       // for wild bootstraps 
*   spost9_ado     // For countfit
*   hnblogit       // hurdle logit - negative binomial 
*   qcount         // quantile count regression  
* are used

********** SETUP **********

set more off
version 15
clear all
set scheme s1mono  /* Graphics scheme */

*********** 1A: COUNT REGRESSION: BASICS (POISSON, GLM, NLS, DIAGNOSTICS, NB)

*** 1A. BASICS: ANALYSIS WITHOUT REGRESSORS

* Summarize doctor visits data (MEPS 2003 data 65-90 years old: see MUS p.557) 
use mus17data.dta
summarize docvis
summarize docvis, detail
histogram docvis if docvis < 40, discrete ///
  xtitle("# Doctor visits (for < 40 visits)")
graph export countfig01histraw.wmf, replace

* Tabulate docvis after recoding values > 10 to ranges 11-40 and 41-60
generate dvrange = docvis
recode dvrange (11/40 = 40) (41/60 = 60)
label define dvcounts 40 "11-40" 60 "41-60"
label values dvrange dvcounts
tabulate dvrange

*** 1A. BASICS: POISSON REGRESSION

drop _all 

* Summary statistics for doctor visits data
use mus17data.dta              
global xlist private medicaid age age2 educyr actlim totchr 
describe docvis $xlist
summarize docvis $xlist, sep(10)

* Poisson with preffered robust standard errors
poisson docvis $xlist, vce(robust) nolog  // Poisson robust SEs

* Poisson with default ML standard errors
poisson docvis $xlist                     // Poisson default ML standard errors 

* Poisson coefficients as incidence ratios
poisson docvis $xlist, irr vce(robust) nolog  // Incidence ratios

* AME and MEM for Poisson 
quietly poisson docvis $xlist, vce(robust) 
margins, dydx(*)        // AME: Average marginal effect for Poisson
margins, dydx(*) atmean // MEM: ME for Poisson evaluated at average of x

* Older code uses add-ons mfx and margeff now superceded by margins     
* margeff  // AME: Average marginal effect for Poisson
* mfx      // MEM: Marginal effect for Poisson evaluated at average of x

* Also factor variables and noncaculus methods
* The i. are discrete and will calculate ME of one unit change (not derivative)
* The c.age##age means age and age-sauared appear and ME is w.r.t. age
poisson docvis i.private i.medicaid c.age##c.age educyr i.actlim totchr, vce(robust)
margins, dydx(*) // MEM: ME for Poisson evaluated at average of x

*** 1A. BASICS: GLM REGRESSION

* GLM for Poisson with robust standard errors
glm docvis $xlist, family(poisson) link(log) vce(robust) nolog

* GLM for Poisson with GLM corrected standard errors
* These are default ML standard errors multiplied by square root Pearson statistic 
glm docvis $xlist, family(poisson) link(log) scale(x2) nolog 

*** 1A. BASICS: NONLINEAR LEAST SQUARES  

* Nonlinear least squares with exponential conditional mean
nl (docvis = exp({xb: $xlist}+{b0})), vce(robust) nolog

*** 1A. BASICS: NEGATIVE BINOMIAL REGRESSION MOTIVATION

*** IMPORTANT: The following uses the old Stata 13 random number generator (kiss32)

set rng kiss32

* Generate Poisson data with sample mean of docvis
quietly summarize docvis
scalar dvmean = r(mean)
set seed 10101                      // set the seed !
generate ypoisson= rpoisson(dvmean) // Draw from Poisson(mu=dvmean) 
summarize docvis ypoisson
tabulate ypoisson
histogram docvis if docvis < 40, discrete addplot(hist ypoisson, fcolor(white) discrete) ///
  xtitle("# Doctor visits versus Poisson") legend(off)
graph export countfig02histpoisson.wmf, replace

* Generate negative binomial (mu=1 var=2) generated data
quietly nbreg docvis
scalar alpha = e(alpha)
display alpha
set seed 10101                      // set the seed !
generate mugamma = rgamma(1/alpha,alpha)
generate ynegbin = rpoisson(dvmean*mugamma) // NB generated as a Poisson-gamma mixture 
summarize docvis ynegbin mugamma 
tabulate ynegbin
histogram docvis if docvis < 40, discrete addplot(hist ynegbin if ynegbin < 40, fcolor(white) discrete) ///
  xtitle("# Doctor visits versus negative binomial") legend(off)
graph export countfig03histnegbin.wmf, replace

* Standard negative binomial (NB2) with default SEs
nbreg docvis $xlist, nolog 

* Newer command - use rather than countfit below  
chi2gof, cells(11) table

* NB2: Sample vs average predicted probabilities of y = 0, 1, ..., 10 
countfit docvis $xlist, maxcount(10) nbreg nograph noestimates nofit

* Generalized negative binomial with alpha parameterized
gnbreg docvis $xlist, lnalpha($xlist) nolog

* Comparison of Poisson and NB and various standard error estimates
quietly poisson docvis $xlist
estimates store PDEFAULT
quietly poisson docvis $xlist, vce(robust)
estimates store PROBUST
quietly glm docvis $xlist, family(poisson) link(log) scale(x2) 
estimates store PPEARSON
quietly nbreg docvis $xlist
estimates store NBDEFAULT
quietly nbreg docvis $xlist, vce(robust)
estimates store NBROBUST
esttab PDEFAULT PROBUST PPEARSON NBDEFAULT NBROBUST, b(%10.4f) se mtitles pr2

*********** 1B. COUNT REGRESSION: INFERENCE (STANDARD ERRORS AND BOOTSTRAP) 

*** 1B. INFERENCE: ROBUST STANDARD ERRORS

* Inference - Robust standard errors
clear
use mus17data.dta              
global xlist private medicaid age age2 educyr actlim totchr 

* Poisson with preferred robust standard errors
poisson docvis $xlist, vce(robust) nolog  // Poisson robust SEs

* Poisson with cluster robust standard errors - illustration
* Here cluster on age for illustration
* In practice the grouping variable would be, for example, village
poisson docvis $xlist, vce(cluster age) nolog  // Poisson robust SEs

*** 1B. BOOTSTRAP PAIRS STANDARD ERRORS

* Small data set for bootstraps so doesn't take long
clear
use musbootdata.dta 
describe 
summarize

* Default and robust standard errors
poisson docvis chronic, nolog
poisson docvis chronic, nolog vce(robust)

*** IMPORTANT: Stata 14 uses different random number generator (mt64) 
***            than earlier versions of Stata. 
***            My slides use the earlier random number generator (kiss32)
***            So here I use the command set rng kiss32

* Use the old Stata 13 and earlier random number generator
set rng kiss32

* Compute bootstrap standard errors using option vce(bootstrap) to 
poisson docvis chronic, vce(boot, reps(400) seed(10101) nodots)

* Same using bootstrap prefix command
bootstrap, seed(10101) reps(400) nodots: poisson docvis chronic

/* SKIP THE FOLLOWING AS TAKES TIME

* Bootstrap standard errors for different reps and seeds
quietly poisson docvis chronic, vce(boot, reps(50) seed(10101))
estimates store boot50
quietly poisson docvis chronic, vce(boot, reps(50) seed(20202))
estimates store boot50diff
quietly poisson docvis chronic, vce(boot, reps(2000) seed(10101))
estimates store boot2000
quietly poisson docvis chronic, vce(robust)
estimates store robust
estimates table boot50 boot50diff boot2000 robust, b(%8.5f) se(%8.5f)

*/

* The jackknife
poisson docvis chronic, vce(jackknife)

* Cluster bootstrap - cluster on age for illustration
poisson docvis chronic, vce(boot, cluster(age) reps(400) seed(10101) nodots)

* Bootstrap confidence intervals: normal-based, percentile, BC, and BCa
quietly poisson docvis chronic, vce(boot, reps(999) seed(10101) bca)
estat bootstrap, all

****** PERCENTILE-T BOOTSTRAP WITH ASYMPTOTIC REFINEMENT

* Percentile-t for a single coefficient: Bootstrap the t statistic
use musbootdata.dta, clear 
quietly poisson docvis chronic, vce(robust)
local theta = _b[chronic] 
local setheta = _se[chronic]
bootstrap tstar=((_b[chronic]-`theta')/_se[chronic]), seed(10101)        ///
  reps(999) nodots saving(percentilet, replace): poisson docvis chronic, ///
  vce(robust)

* Simple plot of the distribution of the tstar
use percentilet, clear
summarize
centile tstar, c(2.5,50,97.5)
kdensity tstar, normal

* Fancier plot of the distribution of the tstar
use percentilet, clear
summarize
centile tstar, c(2.5,50,97.5)
kdensity tstar, generate(evalpoint densityest) xtitle("tstar from the bootstrap replications")
generate phistnorm = normalden(evalpoint)
label variable phistnorm "Standard normal"
label variable densityest "Bootstrap density"
label variable evalpoint "tstar"
graph twoway (scatter densityest evalpoint, connect(l) msize(tiny)) ///
  (scatter phistnorm evalpoint, connect(l) msize(tiny))
* graph export ct_bootstrap1.wmf, replace

* Percentile-t p-value for symmetric two-sided Wald test of H0: theta = 0
use percentilet, clear
quietly count if abs(`theta'/`setheta') < abs(tstar)
display "p-value = " r(N)/(_N+1)

* Percentile-t critical values and confidence interval
_pctile tstar, p(2.5,97.5) 
scalar lb = `theta' + r(r1)*`setheta'
scalar ub = `theta' + r(r2)*`setheta'
display "2.5 and 97.5 percentiles of t* distn: " r(r1) ", " r(r2) _n ///
    "95 percent percentile-t confidence interval is:  (" lb   ","  ub ")"

/* Following did not work
* Command boottest
* Wild score bootstrap, Rademacher weights, null imposed, 999 replications
poisson docvis chronic
boottest chronic 
*/

*********** 2: COUNT REGRESSION: EXTRAS
*              (TRUNCATED/CENSORED, HURDLE, ZERO-INFLATED, FINITE MIXTURE, DIAGNOSTICS, ENDOGENOUS)

*** 2. EXTRAS: TRUNCATED AND CENSORED

* Summary statistics for doctor visits data
use mus17data.dta, clear    
global xlist private medicaid age age2 educyr actlim totchr 
describe docvis $xlist
summarize docvis $xlist, sep(10)

* Zero-truncated negative binomial 
tnbreg docvis $xlist if docvis>0, ll(0) nolog
estimates store ZTNB

* Replaces old command ztnb docvis $xlist if docvis>0, nolog 
* And only does lower truncationhelp ztpoisson

* Censored negative binomial at y <= 10
generate docviscens10 = docvis
replace docviscens10 = 10 if docvis > 10
 
* Censored Negbin ML program lfnbcens10 with censoring from right at 10
program lfnbcens10
  version 10.1
  args lnf theta1 a               // theta1=x'b, a=alpha, lnf=lnf(y)
  tempvar mu sumgammaratios probyle10 cens
  local y $ML_y1                  // Define y so program more readable
  generate double `mu'  = exp(`theta1')
  generate double `sumgammaratios'  = exp(lngamma(0+(1/`a')))/1 ///
    + exp(lngamma(1+(1/`a')))/1 + exp(lngamma(2+(1/`a')))/2     ///
    + exp(lngamma(3+(1/`a')))/6 + exp(lngamma(4+(1/`a')))/24    ///
    + exp(lngamma(5+(1/`a')))/120 + exp(lngamma(6+(1/`a')))/720 ///
    + exp(lngamma(7+(1/`a')))/5040 + exp(lngamma(8+(1/`a')))/40320 ///
    + exp(lngamma(9+(1/`a')))/362880 + exp(lngamma(10+(1/`a')))/3628800 
  generate double `probyle10'  = (`sumgammaratios'/exp(lngamma((1/`a')))) ///
     *((1/`a')^(1/`a'))*(`mu'^`y')/((1/`a'+`mu')^(1/`a'+`y'))
  generate double `cens' = `y' >= 10
  quietly replace `lnf' = (1-`cens') * (lngamma(`y'+(1/`a')) - lngamma(1/`a') ///
               -  lnfactorial(`y') - (`y'+(1/`a'))*ln(1+`a'*`mu')  ///
               + `y'*ln(`a') + `y'*ln(`mu') ) + `cens'*ln(`probyle10')
end
* Command lfnbcens10 implemented for negative binomial MLE
ml model lf lfnbcens10 (docviscens10 = $xlist) ()
ml maximize, nolog
estimates store NBCENS10 

/* Following did not work
  generate double `sumgammaratios'  = 0
  forvalues i = 0/10 {
    quietly replace `sumgammaratios' = `sumgammaratios' + exp(lngamma(`i'+(1/`a')))/exp(lngamma(`i'+1))
    }
*/

* Comparison with regular negative binomial
quietly nbreg docvis $xlist, nolog 
estimates store NBREG
estimates table NBREG ZTNB NBCENS10, equation(1) b(%10.4f) se  stats(N ll)

* Ordered logit
* Do for docviscens10 which has 10 categories
ologit docviscens10 $xlist, nolog 

*** 2. EXTRAS: HURDLE MODEL

* Hurdle logit-nb model manually
logit docvis $xlist, nolog
predict dvplogit, p       // Logit Pr[y=1] 
* Second step uses positives only
ztnb docvis $xlist if docvis>0, nolog 
predict dvztnbcm, cm     // ztnb E[y|y>0] = E[y](1-/Pr[|y=0]
/* Following computes this manually
scalar alpha = e(alpha)
predict dvztnb, n
generate pryeq0 = (1+alpha*dvztnb)^(-1/alpha)
generate dvztnbcmman = dvztnb/(1-pryeq0)
*/
generate dvhurdle = dvplogit*dvztnbcm

* Same hurdle logit-nb model using the user-written hnblogit command
hnblogit docvis $xlist, nolog
estimates store HURDLENB

*** 2. EXTRAS: ZERO-INFLATED

* Zero-inflated negative binomial
zinb docvis $xlist, inflate($xlist) nolog
estimates store ZINB
predict dvzinb, n

quietly nbreg docvis $xlist
estimates store NBREG
predict dvnbreg, n

*** 2. EXTRAS: CONTINUOUS MIXTURE AND HIERARCHICAL MODEL

* Poisson - normal mixture
gsem (docvis <- $xlist M1[educyr]), poisson

* Hierarchical model with vary by cluster on educyr
gsem (docvis <- $xlist M1[educyr]), poisson

*** 2. EXTRAS: MODEL COMPARISON

* Comparison of NB and ZINB using countfit
countfit docvis $xlist, nbreg prm nograph noestimates

* Compare LL, AIC, BIC across models
estimates table NBREG HURDLENB ZINB, equation(1) b(%12.4f) se  ///
   stats(N ll aic bic) stfmt(%12.1f)  newpanel

* Compare predicted means across models
summarize docvis dvnbreg dvhurdle dvzinb
correlate docvis dvnbreg dvhurdle dvzinb

*** 2. EXTRAS: QUANTILE COUNT REGRESSION - 2 COMPONENT NEGATIVE BINOMIAL

* To speed up number of reps is 50 but should be e.g. 1000
qcount docvis $xlist, q(0.5) rep(50)
* My data set had a variable one in it and this causes problems for qcount_mfx
capture drop one
* Marginal effects
qcount_mfx // MEM

*** 2. EXTRAS: FINITE MIXTURE MODEL - 2 COMPONENT NEGATIVE BINOMIAL

* Finite-mixture model using fmm command with constant probabilities
use mus17data.dta, clear
fmm 2, nolog: nbreg docvis $xlist

* Marginal predicted means for each component
estat lcmean

* Marginal predicted probabilities for each component
estat lcprob

* Marginal effects averaged across the two components
margins, dydx(*)

* Predict y for two components
estimates store FMM2
predict dvfit*
summarize dvfit1 dvfit2

* Create histograms of fitted values  
quietly histogram dvfit1 if dvfit1< 40, saving(graph1.gph, replace)
quietly histogram dvfit2 if dvfit1< 40, saving(graph2.gph, replace)
quietly graph combine graph1.gph graph2.gph, ysize(3) xsize(6) ///
    ycommon xcommon iscale(1.15)
quietly graph export countfig21finitemixture.wmf, replace

* Finite-mixture model poisson
use mus17data.dta, clear
fmm 2, nolog: poisson docvis $xlist

*** 2. EXTRAS: DIAGNOSTICS

* Residuals
quietly glm docvis $xlist, family(poisson) link(log)
predict rraw, response 
predict rpearson, pearson
predict rdeviance, deviance
predict ranscombe, anscombe
predict hat, hat
predict cooksd, cooksd
summarize rraw rpearson rdeviance ranscombe hat cooksd, sep(10)
correlate rraw rpearson rdeviance ranscombe

* An R-squared measure
capture drop yphat
quietly poisson docvis $xlist, vce(robust)
predict yphat, n
quietly correlate docvis yphat
display "Squared correlation between y and yhat = " r(rho)^2

* Overdispersion test against V[y|x] = E[y|x] + a*(E[y|x]^2)
quietly poisson docvis $xlist, vce(robust)
predict muhat, n
quietly generate ystar = ((docvis-muhat)^2 - docvis)/muhat
regress ystar muhat, noconstant noheader

* Newer command - use rather than countfit below  
quietly poisson docvis $xlist
chi2gof, cells(11) table

* Poisson: Sample vs avg predicted probabilities of y = 0, 1, ..., 10 
countfit docvis $xlist, maxcount(10) prm nograph noestimates nofit

*** 2. EXTRAS: ENDOGENOUS REGRESSORS

use mus17data.dta, clear
global xlist private medicaid age age2 educyr actlim totchr 

***** Approach 1: Conditional Moments

* GMM (Nonlinear IV) for Poisson: using command gmm
global zlist income ssiratio medicaid age age2 educyr actlim totchr
gmm (docvis - exp({xb:$xlist}+{b0})), instruments($zlist) onestep nolog

* GMM (Nonlinear IV) for Poisson: computation using command optimize
generate cons = 1
local y docvis
local xlist private medicaid age age2 educyr actlim totchr cons
local zlist income ssiratio medicaid age age2 educyr actlim totchr cons
mata
  void pgmm(todo, b, y, X, Z, Qb, g, H)
  {
    Xb = X*b'
    mu = exp(Xb)
    h = Z'(y-mu)
    W = cholinv(cross(Z,Z))
    Qb = h'W*h
    if (todo == 0) return
    G = -(mu:*Z)'X
    g = (G'W*h)'
    if (todo == 1) return
    H = G'W*G
    _makesymmetric(H)
   }
  st_view(y=., ., "`y'")
  st_view(X=., ., tokens("`xlist'"))
  st_view(Z=., ., tokens("`zlist'"))
  S = optimize_init()
  optimize_init_which(S,"min")
  optimize_init_evaluator(S, &pgmm())
  optimize_init_evaluatortype(S, "d2")
  optimize_init_argument(S, 1, y)
  optimize_init_argument(S, 2, X)
  optimize_init_argument(S, 3, Z)
  optimize_init_params(S, J(1,cols(X),0))
  optimize_init_technique(S,"nr")
  b = optimize(S)
  // Compute robust estimate of VCE
  Xb = X*b'
  mu = exp(Xb)
  h = Z'(y-mu)
  W = cholinv(cross(Z,Z))
  G = -(mu:*Z)'X
  Shat = ((y-mu):*Z)'((y-mu):*Z)*rows(X)/(rows(X)-cols(X))
  Vb = luinv(G'W*G)*G'W*Shat*W*G*luinv(G'W*G)
  st_matrix("b",b)
  st_matrix("Vb",Vb)
end

* Nonlinear IV estimator for Poisson: formatted results
matrix colnames b = `xlist'
matrix colnames Vb = `xlist'
matrix rownames Vb = `xlist'
ereturn post b Vb
ereturn display

***** Approach 2: Control Function but not MLE

* Control function approach to endogeneity
* First stage is reduced form to get predicted residuals
global xlist2 medicaid age age2 educyr actlim totchr
regress private $xlist2 income ssiratio, vce(robust)
predict lpuhat, residual

* Second-stage Poisson with robust SEs
poisson docvis private $xlist2 lpuhat, vce(robust) nolog

* Program and bootstrap for Poisson two-step estimator
program endogtwostep, eclass
  version 10.1
  tempname b
  capture drop lpuhat2
  regress private $xlist2 income ssiratio
  predict lpuhat2, residual
  poisson docvis private $xlist2 lpuhat2
  matrix `b' = e(b)
  ereturn post `b'
end
bootstrap _b, reps(400) seed(10101) nodots nowarn: endogtwostep

***** Approach 3: Fully structural MLE 

* Fully structural model of endogeneity
gsem (docvis <- private $xlist2 L, poisson) ///
   (private <- $xlist2 income ssiratio L), nolog

********** CLOSE OUTPUT

* log close
* clear
* exit

