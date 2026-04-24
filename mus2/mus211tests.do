* mus211tests.do   for Stata 17

capture log close

********** OVERVIEW OF mus211tests.do **********

* Stata program
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi
* used for "Microeconometrics Using Stata, Second Edition"
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 11
* 11.2 CRITICAL VALUES AND P-VALUES
* 11.3 WALD TESTS AND CONDFIDENCE INTERVALS
* 11.4 LR TESTS
* 11.5 LM TEST (OR SCORE TEST)
* 11.6 MULTIPLE TESTING
* 11.7 TEST SIZE AND POWER
* 11.8 THE POWER ONEMEAN COMMAND FOR MULTIPLE REGRESSION
* 11.10 PERMUTATION TESTS AND RANDOMIZATION TESTS

* To run you need files
*   mus210mepsdocvisyoung.dta
* in your directory

* community-contributed command
*   multproc
* is used

********** SETUP **********

clear all
set linesize 81
set scheme s1mono   /* Used for graphs */

********** DATA DESCRIPTION **********

* mus210mepsdocvisyoung.dta from 2002 Medical Expenditure Panel Survey (MEPS)
* U.S. data on office-based physician visits by persons aged 25-64 years
* Same as Deb, Munkin and Trivedi, "Bayesian analysis of the two-part model
* with endogeneity: application to health care expenditure", JBES, 1081-1099
* Excludes those receiving public insurance (Medicare and Medicaid)
* Restricted to those working in the private sector but not self-employed.

********** 11.2 CRITICAL VALUES AND P-VALUES

* Create many draws from chi(5) and 5*F(5,30) distributions
set seed 10101
qui set obs 10000
generate chi5 = rchi2(5)           // Result xc ~ chisquared(10)
generate xfn = rchi2(5)/5          // Result numerator of F(5,30)
generate xfd = rchi2(30)/30        // Result denominator of F(5,30)
generate f5_30 = xfn/xfd           // Result xf ~ F(5,30)
generate five_x_f5_30 = 5*f5_30
summarize chi5 five_x_f5_30

* Plot the densities for these two distributions using kdensity
label var chi5 "chi(5)"
label var five_x_f5_30 "5*F(5,30)"
kdensity chi5, bw(1.0) generate(kx1 kd1) n(500)
kdensity five_x_f5_30, bw(1.0) generate(kx2 kd2) n(500)
qui drop if (chi5 > 25  |  five_x_f5_30 > 25)
graph twoway (line kd1 kx1) (line kd2 kx2, clstyle(p3)) if kx1 < 25, ///
    scale(1.2) plotregion(style(none))                               ///
    title("{&chi}{sup:2}(5) and 5*F(5,30) Densities")                ///
    xtitle("y", size(medlarge)) xscale(titlegap(*5))                 ///
    ytitle("Density f(y)", size(medlarge)) yscale(titlegap(*5))      ///
    legend(pos(1) ring(0) col(1)) legend(size(small))                ///
    legend(label(1 "{&chi}{sup:2}(5)") label(2 "5*F(5,30)"))

* p-values for t(30), F(1,30), Z, and chi(1) at y = 2
scalar y = 2
scalar p_t30 = 2*ttail(30,y)
scalar p_f1and30 = Ftail(1,30,y^2)
scalar p_z = 2*(1-normal(y))
scalar p_chi1 = chi2tail(1,y^2)
display "p-values" "  t(30) =" %7.4f p_t30 "  F(1,30)=" %7.4f ///
    p_f1and30 "  z =" %7.4f p_z "  chi(1)=" %7.4f p_chi1

* Critical values for t(30), F(1,30), Z, and chi(1) at level 0.05
scalar alpha = 0.05
scalar c_t30 = invttail(30,alpha/2)
scalar c_f1and30 = invFtail(1,30,alpha)
scalar c_z = -invnormal(alpha/2)
scalar c_chi1 = invchi2(1,1-alpha)
display "critical values" "  t(30) =" %7.3f c_t30 "  F(1,30)=" %7.3f ///
    c_f1and30 "  z =" %7.3f c_z "  chi(1)=" %7.3f c_chi1

********** 11.3 WALD TESTS AND CONFIDENCE INTERVALS

* Fit Poisson model used throughout this chapter
qui use mus210mepsdocvisyoung, clear
qui keep if year02==1
poisson docvis private chronic female income, vce(robust) nolog

* Test a single coefficient equal 0
test female

// Not included
test female
test female = 0
test [docvis]female = 0
test _b[female] = 0

* Test two hypotheses jointly using test
test (female) (private + chronic = 1)

* Test each hypothesis in isolation as well as jointly
test (female) (private + chronic = 1), mtest

test ([docvis]female = 0) ([docvis]private + [docvis]chronic = 1)
test (_b[female] = 0) (_b[private] + _b[chronic] = 1)

* Wald test of overall significance
test private chronic female income

* Manually compute overall test of significance using the formula for W
qui poisson docvis private chronic female income, vce(robust)
matrix b = e(b)'
matrix V = e(V)
matrix R = (1,0,0,0,0 \ 0,1,0,0,0 \ 0,0,1,0,0 \ 0,0,0,1,0 )
matrix r = (0 \ 0 \ 0 \ 0)
matrix W = (R*b-r)'*invsym(R*V*R')*(R*b-r)
scalar Wald = W[1,1]
scalar h = rowsof(R)
display "Wald test statistic: " Wald "  with p-value: " chi2tail(h,Wald)

* Test a nonlinear hypothesis using testnl
testnl _b[female]/_b[private] = 1

* Wald test is not invariant
test female = private

* Stepwise forward selection using statistical significance at 5%
stepwise, pe(.05): poisson docvis private chronic female income ///
    firmsize msa injury, vce(robust)

* Confidence interval for linear combinations using lincom
qui use mus210mepsdocvisyoung, clear
qui keep if year02==1
qui poisson docvis private chronic female income if year02==1, vce(robust)
lincom private + chronic - 1

* Confidence interval for nonlinear function of parameters using nlcom
nlcom _b[female] / _b[private] - 1

* Confidence interval for exp(b) using lincom option eform
lincom private, eform

* Confidence interval for exp(b) using lincom followed by exponentiate
lincom private

* Confidence interval for exp(b) using nlcom
nlcom exp(_b[private])

* Confidence interval from inverting test statistic
qui poisson docvis private chronic female income, vce(robust) nolog
postfile cifromtest b2 pvalue using pvalues, replace
forvalues i = 1/1000 {
    scalar b2 = (`i' - 500)/250
    qui test _b[private] = b2
    scalar p = r(p)
    post cifromtest (b2) (p)
}
postclose cifromtest
use pvalues, clear
sum if pvalue > 0.05

// Not included - verify visually that no gaps in the confidence interval
// scatter b2 p if pvalue > 0.05

********** 11.4 LIKELIHOOD RATIO TESTS

* LR tests output if estimate by ML with default estimate of VCE
qui use mus210mepsdocvisyoung, clear
qui keep if year02==1
nbreg docvis private chronic female income,  nolog

*testnl exp([lnalpha]_cons)=0

* LR test using command lrtest
qui nbreg docvis private chronic female income
estimates store unrestrict
qui nbreg docvis female income
estimates store restrict
lrtest unrestrict restrict

lrtest restrict unrestrict

* Wald test of the same hypothesis and using default standard errors (like LR)
qui nbreg docvis private chronic female income
test chronic private

* LR test using option force
qui nbreg docvis private chronic female income
estimates store nb
qui poisson docvis private chronic female income
estimates store poiss
lrtest nb poiss, force
display "Corrected p-value for LR-test = " r(p)/2

********** 11.5 LAGRANGE MULTIPLIER TEST (OR SCORE TEST)

* Perform LM test that b_private=0, b_chronic=0 using auxiliary regression
qui use mus210mepsdocvisyoung, clear
qui keep if year02==1
generate one = 1
constraint define 1 private = 0
constraint define 2 chronic = 0
qui nbreg docvis female income private chronic, constraints(1 2)
predict eqscore ascore, scores
generate s1restb = eqscore*one
generate s2restb = eqscore*female
generate s3restb = eqscore*income
generate s4restb = eqscore*private
generate s5restb = eqscore*chronic
generate salpha = ascore*one
qui regress one s1restb s2restb s3restb s4restb s5restb salpha, noconstant
scalar lm = e(N)*e(r2)
display "LM = N x uncentered Rsq = " lm " and p = " chi2tail(2,lm)

********** CHAPTER 11.6: MULTIPLE TESTING

* Multiple tests - form the age subgroups 25-34, 35-44, 45-54, 55-64
qui use mus210mepsdocvisyoung, clear
keep if _n <= 500
generate agegroup = round(age)-2   // Age is in tens of years
tabulate agegroup

* Using a scalar did not work below. Need to use a global.
* scalar nummodels = 4

* Multiple tests for multiple subgroups using p-value corrections
global nummodels 4
forvalues i = 1/$nummodels {
    qui poisson docvis private chronic female income ///
        if agegroup==`i', vce(robust)
    scalar beta = _b[private]
    scalar p = 2*(1-normal(abs(_b[private]/_se[private])))
    scalar pbonferroni = min(1,p*$nummodels)
    scalar psidak = 1-(1-p)^$nummodels
    di "Group " `i' " b =" %7.4f beta "  p-values:" " Usual =" %6.3f p ///
        " Bonferroni =" %6.3f pbonferroni " Sidak =" %6.3f psidak
}

clear
* Multiple tests using community-contributed command multproc and Holm's method
input pvalues
    .032
    .000
    .919
    .024
end
multproc, puncor(0.05) pvalue(pvalues) method(holm) rank(prank) ///
   gpuncor(alpha) gpcor(indivalpha) nhcred(credible) reject(reject)
list pvalues indivalpha reject credible prank alpha in 1/4, clean

* The following did not work. Why ?
* scalar nummodels = 4

* Multiple tests for multiple subgroups using mtests option of test
qui use mus210mepsdocvisyoung, clear
keep if _n <= 500
generate agegroup = round(age)-2   // Age is in tens of years
generate one = 1
qui poisson docvis i.agegroup#c.(one private chronic female income), ///
    noconstant vce(robust)
qui matrix list e(b)  // Yields the complicated names for the coefficients
test ([docvis]:1b.agegroup#c.private) ([docvis]:2.agegroup#c.private=0) ///
    ([docvis]:3.agegroup#c.private=0) ([docvis]:4.agegroup#c.private=0), mtest(s)

* Not included
* Alternative - have only private insurance interacted with age group
qui poisson docvis chronic female income agegroup#c.private, vce(robust)
qui matrix list e(b)  // Yields the complicated names for the coefficients
test 1b.agegroup#c.private 2.agegroup#private 3.agegroup#private ///
    4.agegroup#private, mtest(s)

* Multiple tests for multiple outcomes using p-value corrections
global nummodels 3
global ylist docvis injury educ
foreach var of varlist $ylist {
    qui regress `var' private chronic female income, vce(robust)
    scalar p = 2*ttail(e(df_r),abs(_b[income]/_se[income]))
    scalar pbonferroni = min(1,p*$nummodels)
    scalar psidak = 1-(1-p)^$nummodels
    di "Outcome = " "`var'" _col(19) " p-values:" " Usual =" %7.4f p ///
         " Bonferroni =" %7.4f pbonferroni " Sidak =" %7.4f psidak
}

// Not included - check the examples given in the text
input pvalues
 .002
 .004
 .008
 .012
 .015
 .016
 .034
 .058
 end
multproc, puncor(0.05) pvalue(pvalues) method(holm) rank(prank) ///
   gpuncor(alpha) gpcor(indivalpha) nhcred(credible) reject(reject)
list pvalues indivalpha reject credible prank alpha in 1/8, clean

********** 11.7 TEST SIZE AND POWER

* Do 1,000 simulations where each gets p-value of test of b2=2
set seed 10101
postfile sim pvalues using pvalues, replace
forvalues i = 1/1000 {
    drop _all
    qui set obs 150
    qui generate double x = rchi2(1)
    qui generate y = 1 + 2*x + rchi2(1)-1
    qui regress y x
    qui test x = 2
    scalar p = r(p)         // p-value for test this simulation
    post sim (p)
}
postclose sim

* Summarize the p-value from each of the 1,000 tests
use pvalues, clear
summarize pvalues

* Determine size of test at level 0.05
count if pvalues < .05
display "Test size from 1000 simulations = " r(N)/1000

* 95% simulation interval using exact binomial at level 0.05 with S=1000
cii proportions 1000 50

* Mentioned in text
cii proportions 10000 500

* Program to compute power of test given specified H0 and Ha values of b2
program mypower, rclass
    version 17
    args numsims numobs b2H0 b2Ha nominalsize
                                        // Setup before simulation loops
    drop _all
    set seed 10101
    postfile sim pvalues using power, replace
                                       // Simulation loop
    forvalues i = 1/`numsims' {
        drop _all
        qui set obs `numobs'
        qui generate double x = rchi2(1)
        qui generate y = 1 + `b2Ha'*x + rchi2(1)-1
        qui regress y x
        qui test x = `b2H0'
        scalar p = r(p)
        post sim (p)
    }
    postclose sim
    use power, clear
                                      // Determine the size or power
    qui count if pvalues < `nominalsize'
    return scalar power=r(N)/`numsims'
end

* Size = power of test of b2H0=2 when b2Ha=2, S=1000, N=150, alpha=0.05
mypower 1000 150 2.00 2.00 0.05
display r(power) " is the test size"

* Power of test of b2H0=2 when b2Ha=2.2, S=1000, N=150, alpha=0.05
mypower 1000 150 2.00 2.20 0.05
display r(power) " is the test power"

* Power of test of H0:b2=2 against Ha:b2=1.6,1.625, ..., 2.4
postfile simofsims b2Ha power using simresults, replace
forvalues i = 0/33 {
    drop _all
    scalar b2Ha = 1.6 + 0.025*`i'
    mypower 1000 150 2.00 b2Ha 0.05
    post simofsims (b2Ha) (r(power))
}
postclose simofsims
use simresults, clear
summarize

* Plot the power curve
twoway (connected power b2Ha), scale(1.2) plotregion(style(none))

* Change one example
* Power of chi(1) test when noncentrality parameter lambda = 5.634
display 1-nchi2(1,5.634,3.841)

********** CHAPTER 11.8 THE POWER COMMAND FOR MULTIPLE REGRESSION

* Power of Wald normal test: 2 versus 2.2, size 0.05, s_thetahat=0.0843
power onemean 2 2.2, knownsd alpha(0.05) n(1) sd(0.0843)

* Min effect size of normal test: 2 versus ?, size 0.05, power 0.80, s_t=.0843
power onemean 2, knownsd alpha(0.05) power(0.80) n(1) sd(0.0843) ///
    table(alpha power m0 ma sd N delta, formats(power "%7.3f"))

* Aside: 5% test value
di 2.0+1.960*0.0843

* Min sample size for normal test: 2 versus 2.2, size 0.05, power 0.80, independent
di sqrt(150)*0.0843
power onemean 2 2.2, knownsd alpha(0.05) power(0.80) sd(1.0325) ///
    table(alpha power m0 ma sd N delta, formats(power "%7.3f"))

* General power curve for N(0,1) test at size = 0.05
clear
power onemean 0 (-4(0.1)4), knownsd sd(1) n(1) alpha(0.05) ///
    graph(legend(off) scale(1.2) title("") subtitle("") ///
        note("") recast(line)  ///
    ylabel(, grid angle(0)) yline(0.05) xlabel(, grid)       ///
    xtitle("Ha: Number of standard errors from H0 value")    ///
    ytitle("Test power") saving(graph1.gph, replace))

* Used to create the second panel in the figure
clear
power onemean 0 (-4(0.1)4), knownsd sd(1) n(1) alpha(0.05 0.1) ///
    graph(legend(off) scale(1.2) title("") subtitle("") ///
        note("") recast(line)  ///
    ylabel(, grid angle(0)) yline(0.05 0.1) xlabel(, grid)   ///
    xtitle("Ha: Number of standard errors from H0 value")    ///
    ytitle("Test power") saving(graph2.gph, replace))


* Then combine
graph combine graph1.gph graph2.gph, iscale(1.2) rows(1) ycommon xcommon ///
    ysize(2.5) xsize(6)

// Not included
* Min effect size for t(148) test: 2 versus ?, size 0.05, power 0.80, s_t=0.0843
di sqrt(149)*0.0843
power onemean 2, alpha(0.05) power(0.80) n(12) sd(1.029) ///
    table(alpha power m0 ma sd N delta, formats(power "%7.3f"))
di "Size: " r(alpha) "  Power: " r(power) "  Min effect size: " r(delta) ///
    "  Target theta: " r(ma)

* Power of Wald t(11) test: 2 versus 2.2, size 0.05, s_thetahat=0.0843
di sqrt(12)*0.0843
power onemean 2 2.2, alpha(0.05) n(12) sd(0.29202) ///
  table(alpha power m0 ma sd N delta, formats(power "%7.3f"))

* Min effect size of t(11) test: 2 versus ?, size 0.05, power 0.80, s_t=.29202
power onemean 2, alpha(0.05) power(0.80) n(12) sd(0.29202) ///
    table(alpha power m0 ma sd N delta, formats(power "%7.3f"))

* Power for t(11) with specified test size and effect size
scalar df = 11
scalar se = 0.0843
scalar effectsize = 0.2
scalar ncp = effectsize/se
scalar size = 0.05
scalar tcrit = invttail(df, size/2)
display "power = " 1 - nt(df, ncp, tcrit) + - nt(df, ncp, -tcrit)

* Power with df increase from 11 to 23 (for example, double the number of clusters)
scalar df = 23
scalar se = 0.0843/sqrt(2)
scalar effectsize = 0.2
scalar ncp = effectsize/se
scalar size = 0.05
scalar tcrit = invttail(df, size/2)
display "power = " 1 - nt(df, ncp, tcrit) + - nt(df, ncp, -tcrit)

di sqrt(150)*0.0843

* Confidence interval width given estimator precision and sample size
ciwidth onemean, knownsd level(95) sd(1.0325) n(300) ///
    table(level N sd width, formats(level "%7.3f"))

* Minimum sample size for confidence interval of given width
ciwidth onemean, knownsd level(95) sd(1.0325) width(0.2) ///
    table(level N sd width, formats(level "%7.3f"))

**** 11.10 PERMUTATION TESTS AND RANDOMIZATION TESTS

* Permutation test example - DGP and usual test
clear
set obs 100
set seed 10101
gen d = runiform() > 0.5
gen y = rnormal(1,1) + 0.2*d
sum y d
regress y d, noheader

* Permutation test
capture program drop pols
program pols, rclass
    version 17
    args y d
    regress `y' `d'
    return scalar beta = _b[d]
    return scalar t = _b[d]/_se[d]
end
permute y t=r(t) beta=r(beta), nowarn nodots seed(10101) reps(1000): pols y d

*** Erase files created by this program and not used elsewhere in the book

erase power.dta
erase simresults.dta
erase pvalues.dta
erase graph1.gph
erase graph2.gph

********** END
