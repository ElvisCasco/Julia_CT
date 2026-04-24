* count2017.do based on trcountimf2010.do August 2017 for Stata version 14

capture log close
log using count2017.txt, text replace

********** OVERVIEW OF count2017.do **********

* STATA Program to demonstrate Cound Data Regression Models
* Based on mus17p1cnt.do
* A. Colin Cameron and Pravin K. Trivedi (2009)
* "Microeconometrics using Stata", Stata Press

* 1. COUNT REGRESSION: BASIC CROSS-SECTION
* 2: COUNT REGRESSION: ADDITIONAL (TRUNCATED/CENSORED, HURDLE, ZERO-INFLATED, FINITE MIXTURE, ENDOGENOUS)
* 3: COUNT REGRESSION: BASIC PANEL

* To run you need files
*    mus17data.dta
*    mus18data.dta    
* in your directory

* Stata user-written commands
*   hnblogit       // for hurdle logit
*   chi2gof        // For chis-squared goodness of fit test
*   qcount         // quantile count regression   
* are used

********** SETUP **********

clear all
set more off
version 14
set linesize 82
set scheme s1mono  /* Graphics scheme */

*********** 1: COUNT REGRESSION: BASICS (POISSON, GLM, NLS, DIAGNOSTICS, NB)

*** 1. BASICS: ANALYSIS WITHOUT REGRESSORS

* Summarize doctor visits data (MEPS 2003 data 65-90 years old: see MUS p.557) 
use mus17data.dta
histogram docvis if docvis < 40, discrete ///
  xtitle("# Doctor visits (for < 40 visits)")
graph export countfig01histraw.wmf, replace

* Tabulate docvis after recoding values > 10 to ranges 11-40 and 41-60
generate dvrange = docvis
recode dvrange (11/40 = 40) (41/60 = 60)
label define dvcounts 40 "11-40" 60 "41-60"
label values dvrange dvcounts
tabulate dvrange

*** 1. BASICS: POISSON REGRESSION

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
quietly poisson docvis i.private i.medicaid c.age##c.age educyr i.actlim totchr, vce(robust)
margins, dydx(*) // MEM: ME for Poisson evaluated at average of x

*** 1. BASICS: EXAMPLES OF ML CODING

* Poisson providing the log-density and using mlexp for simple ML model 
mlexp (-exp({xb:$xlist _cons}) + docvis*({xb:}) - lnfactorial(docvis)), ///
   vce(robust)

* Poisson ML program lfpois to be called by command ml method lf
program lfpois
  version 15
  args lnf theta1                  // theta1=x'b, lnf=ln(y)
  tempvar lnyfact mu
  local y "$ML_y1"                 // Define y so program more readable
  generate double `lnyfact' = lnfactorial(`y')
  generate double `mu'      = exp(`theta1')
  quietly replace `lnf'     = -`mu' + `y'*`theta1' - `lnyfact'
end

* Command ml model including defining y and x, plus ml check
ml model lf lfpois (docvis = $xlist), vce(robust)
* ml check
* ml search    // search for better starting values
ml maximize, nolog  

global y docvis
generate cons = 1
global xlist2 private medicaid age age2 educyr actlim totchr cons
sum $y 
sum $xlist2

*** 1. BASICS: POISSON MLE AND NEWTON RAPHSON ITERATIONS IN MATA

* Complete Mata code for Poisson MLE NR iterations
mata:
  st_view(y=., ., "$y")            // read in stata data to y and X
  st_view(X=., ., tokens("$xlist")) 
  b = J(cols(X),1,0)                // compute starting values
  n = rows(X)   
  iter = 1                          // initialize number of iterations
  cha = 1                           // initialize stopping criterion
  do {
     mu = exp(X*b)                
     grad = X'(y-mu)                // k x 1 gradient vector
     hes = cross(X, mu, X)          // negative of the k x k Hessian matrix
     bold = b
     b = bold + cholinv(hes)*grad
     cha = (bold-b)'(bold-b)/(bold'bold)  
     iter = iter + 1
  } while (cha > 1e-16)             // end of iteration loops 
  mu = exp(X*b)
  hes = cross(X, mu, X)         
  vgrad = cross(X, (y-mu):^2, X)
  vb = cholinv(hes)*vgrad*cholinv(hes)*n/(n-cols(X))
  iter                              // number of iterations
  cha                               // stopping criterion
  st_matrix("b",b')                 // pass results from Mata to Stata
  st_matrix("V",vb)                 // pass results from Mata to Stata
end

* Present results, nicely formatted using Stata command ereturn
matrix colnames b = $xlist
matrix colnames V = $xlist
matrix rownames V = $xlist
ereturn post b V
ereturn display   

*** 1. BASICS: GLM REGRESSION

* GLM for Poisson with robust standard errors
glm docvis $xlist, family(poisson) link(log) vce(robust) nolog

* GLM for Poisson with GLM corrected standard errors
* These are default ML standard errors multiplied by square root Pearson statistic 
glm docvis $xlist, family(poisson) link(log) scale(x2) nolog 

*** 1. BASICS: NONLINEAR LEAST SQUARES  

* Nonlinear least squares with exponential conditional mean
nl (docvis = exp({xb: $xlist}+{b0})), vce(robust) nolog

*** 1. BASICS: DIAGNOSTICS

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

*** 1. BASICS: NEGATIVE BINOMIAL REGRESSION MOTIVATION

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

* Generate negative binomial (mu=ymean alpha=fromnbreg) generated data
quietly nbreg docvis           // Negative binomial on sclar givesalpha
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

*** 1. ADDITONAL: TRUNCATED AND CENSORED

use mus17data.dta, clear
global xlist private medicaid age age2 educyr actlim totchr 

* Zero-truncated negative binomial 
ztnb docvis $xlist if docvis>0, nolog 
estimates store ZTNB

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

*** 1. ADDITONAL: HURDLE MODEL

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

*** 1. ADDITONAL: ZERO-INFLATED

* Zero-inflated negative binomial
zinb docvis $xlist, inflate($xlist) vuong nolog
estimates store ZINB
predict dvzinb, n

quietly nbreg docvis $xlist
estimates store NBREG
predict dvnbreg, n

*** 1. ADDITIONAL: MODEL COMPARISON

* Comparison of NB and ZINB using countfit
countfit docvis $xlist, nbreg zinb nograph noestimates

* Compare LL, AIC, BIC across models
estimates table NBREG HURDLENB ZINB, equation(1) b(%12.4f) se  ///
   stats(N ll aic bic) stfmt(%12.1f)  newpanel

* Compare predicted means across models
summarize docvis dvnbreg dvhurdle dvzinb
correlate docvis dvnbreg dvhurdle dvzinb

*** 2. ADDITIONAL: FINITE MIXTURE MODEL

/* Stata 15 has a new command fmm
   The code below is commented out for those who only have Stata 14.
   The Stata command fmm supplants the earlier user-written fmm command
*/   

version 15

* Finite-mixture model using fmm command with constant probabilities
use mus17data.dta, clear
fmm 2, vce(robust): poisson docvis $xlist

* Obtain latent class probabilities 
estat lcprob 

* Obtain predicted mean for each individual by component 
predict yfit*
summarize yfit1 yfit2

* Obtain marginal effects in each component
margins, dydx(*) 
margins, dydx(*)  atmean noatlegend

* Create histograms of fitted values  
quietly histogram yfit1, name(class_1, replace)
quietly histogram yfit2, name(class_2, replace)
quietly graph combine class_1 class_2, iscale(1.1) ycommon xcommon ///
  ysize(2.5) xsize(5) 
graph export countfig21finitemixture.wmf, replace

* 2-component mixture of negative binomial 
fmm 2, nolog: nbreg docvis $xlist
estat lcprob 

*** 2. ADDITONAL: ENDOGENOUS REGRESSORS

use mus17data.dta, clear
global xlist private medicaid age age2 educyr actlim totchr 

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
estimates store twostep

/* Commented out as Stata 15 
* Use gsem which places greater parametric assumptions
global xlist2 medicaid age age2 educyr actlim totchr
gsem (docvis <- private $xlist2 L, poisson) ///
  (private <- $xlist2 income ssiratio L), var(L@1)
estimates store gsem1
gsem (docvis <- private $xlist2 L, poisson) ///
  (private <- $xlist2 income ssiratio L@1, probit), var(L@1)
estimates store gsem2
estimates table twostep gsem1 gsem2, b se 
* The next uses a logit model at the first stage
gsem (docvis <- private $xlist2 L, poisson) ///
  (private <- $xlist2 income ssiratio L, logit), var(L@1)
*/

*** 2. ADDITIONAL: PANEL COUNT REGRESSION

* Describe dependent variables and regressors
use mus18data.dta, clear
describe mdu lcoins ndisease female age lfam child id year

* Summarize dependent variables and regressors
summarize mdu lcoins ndisease female age lfam child id year

* Panel description of dataset 
xtset id year
xtdescribe

* Panel summary of dependent variable
xtsum mdu

* Year-to-year transitions in doctor visits
generate mdushort = mdu
replace mdushort = 4 if mdu >= 4
xttrans mdushort
corr mdu L.mdu

* Simple time-series plot for first 20 individuals (= first 85 obs)
quietly xtline mdu if _n<=85, overlay legend(off)
graph export countfig31panelplot.wmf, replace

*** PANEL POISSON: POOLED AND POPULATION AVERAGED

* Pooled Poisson estimator with cluster-robust standard errors
poisson mdu lcoins ndisease female age lfam child, vce(cluster id)

* Without cluster-robust
poisson mdu lcoins ndisease female age lfam child

* Poisson PA estimator with unstructured error correlation and robust VCE
xtpoisson mdu lcoins ndisease female age lfam child, pa corr(unstr) vce(robust)

* Print out the correlation matrix
matrix list e(R)

*** PANEL POISSON: RANDOM EFFECTS

* Poisson random-effects estimator with default standard errors
xtpoisson mdu lcoins ndisease female age lfam child, re

* Poisson random-effects estimator (gamma) with cluster-robust standard errors
xtpoisson mdu lcoins ndisease female age lfam child, re vce(cluster id)
* Preceding same as xtpoisson , re vce(robust)

* Poisson random-effects estimator with normal intercept and default standard errors
xtpoisson mdu lcoins ndisease female age lfam child, re normal

* Poisson random-effects estimator with normal intercept and normal slope for one parameter
* xtmepoisson mdu lcoins ndisease female age lfam child || id: NDISEASE
* Previous command commented out as takes a long time

*** PANEL POISSON: FIXED EFFECTS

* Poisson fixed-effects estimator with WRONG default standard errors
xtpoisson mdu lcoins ndisease female age lfam child, fe

* Poisson fixed-effects estimator with panel-robust standard errors
xtpoisson mdu lcoins ndisease female age lfam child, fe vce(robust)
* This supplants earlier add-on xtpqml
* xtpqml mdu lcoins ndisease female age lfam child, fe i(id) cluster(id)

*** PANEL POISSON: ESTIMATOR COMPARISON

* Comparison of Poisson panel estimators using panel-robust standard errors
global xlist lcoins ndisease female age lfam child
quietly poisson mdu $xlist, vce(cluster id)
estimates store POOLED
quietly xtpoisson mdu $xlist, pa corr(unstr) vce(robust)
estimates store POPAVE
quietly xtpoisson mdu $xlist, re vce(cluster id)
estimates store RE_GAMMA
quietly xtpoisson mdu $xlist, re normal vce(cluster id)
estimates store RE_NORMAL
quietly xtpoisson mdu $xlist, fe vce(robust)
estimates store FIXED
estimates table POOLED POPAVE RE_GAMMA RE_NORMAL FIXED, equations(1) b(%8.4f) se stats(N ll) stfmt(%8.0f)

********** CLOSE OUTPUT

* log close
* clear
* exit

