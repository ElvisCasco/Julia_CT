* mus12p1tests.do  Oct 2009 for Stata version 11

cap log close

********** OVERVIEW of mus12p1tests.do   **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Chapter 12
* 12.2 CRITICAL VALUES AND P-VALUES
* 12.3 WALD TESTS AND CONDFIDENCE INTERVALS
* 12.4 LIKELIHOOD RATIO TESTS
* 12.5 LAGRANGE MULTIPLIER TEST (OR SCORE TEST)
* 12.6 TEST SIZE AND POWER
* 12.7 SPECIFICATION TESTS

* To run you need files
*   mus10data.dta
* in your directory
* No Stata user-written commands are used

********** SETUP **********

version 11
clear all
set memory 10m
set linesize 81
set more off
set scheme s1mono   /* Used for graphs */

********** DATA DESCRIPTION **********

* 2002 Medical Expenditure Panel Survey (MEPS)
* U.S. individuals aged 25-64 years working in private sector 
* but not self-employed
* and not receiving public insurance (Medicare and Medicaid) 
* Data due to Deb, Munkin and Trivedi (2006)

********** 12.2 CRITICAL VALUES AND P-VALUES

* Create many draws from chi(5) and 5*F(5,30) distributions
set seed 10101
quietly set obs 10000 
generate chi5 = rchi2(5)           // result xc ~ chisquared(10) 
generate xfn = rchi2(5)/5          // result numerator of F(5,30)
generate xfd = rchi2(30)/30        // result denominator of F(5,30)
generate f5_30 = xfn/xfd           // result xf ~ F(5,30)
generate five_x_f5_30 = 5*f5_30
summarize chi5 five_x_f5_30

* Plot the densities for these two distributions using kdensity
label var chi5 "chi(5)"
label var five_x_f5_30 "5*F(5,30)"
kdensity chi5, bw(1.0) generate (kx1 kd1) n(500)
kdensity five_x_f5_30, bw(1.0) generate (kx2 kd2) n(500)
quietly drop if (chi5 > 25  |  five_x_f5_30 > 25)
graph twoway (line kd1 kx1) if kx1 < 25, name(chi)
graph twoway (line kd2 kx2) if kx2 < 25, name(F)
graph twoway (line kd1 kx1) (line kd2 kx2, clstyle(p3)) if kx1 < 25,  ///
  scale(1.2) plotregion(style(none))     ///
  title("Chisquared(5) and 5*F(5,30) Densities")                ///
  xtitle("y", size(medlarge)) xscale(titlegap(*5))              ///
  ytitle("Density f(y)", size(medlarge)) yscale(titlegap(*5))   ///
  legend(pos(1) ring(0) col(1)) legend(size(small))             ///
  legend(label(1 "Chi(5)") label(2 "5*F(5,30)"))
quietly graph export mus12fandchi5.eps, replace


* p-values for t(30), F(1,30), Z, and chi(1) at y = 2
scalar y = 2
scalar p_t30 = 2*ttail(30,y)
scalar p_f1and30 = Ftail(1,30,y^2)
scalar p_z = 2*(1-normal(y))
scalar p_chi1 = chi2tail(1,y^2)
display "p-values" "  t(30) =" %7.4f p_t30 "  F(1,30)=" %7.4f  ///
       p_f1and30 "  z =" %7.4f p_z "  chi(1)=" %7.4f p_chi1

* Critical values for t(30), F(1,30), Z, and chi(1) at level 0.05
scalar alpha = 0.05
scalar c_t30 = invttail(30,alpha/2)
scalar c_f1and30 = invFtail(1,30,alpha) 
scalar c_z = -invnormal(alpha/2)
scalar c_chi1 = invchi2(1,1-alpha)
display "critical values" "  t(30) =" %7.3f c_t30 "  F(1,30)=" %7.3f  ///
       c_f1and30 "  z =" %7.3f c_z "  chi(1)=" %7.3f c_chi1

********** 12.3 WALD TESTS AND CONFIDENCE INTERVALS

* Estimate Poisson model of chapter 10
use mus10data.dta, clear
quietly keep if year02==1
poisson docvis private chronic female income, vce(robust) nolog

* Test a single coefficient equal 0
test female

// Variations not included in book
test female
test female = 0
test [docvis]female = 0
test _b[female] = 0

* Test two hypotheses jointly using test
test (female) (private + chronic = 1)

* Test each hypothesis in isolation as well as jointly
test (female) (private + chronic = 1), mtest

// Following not included in book
test ([docvis]female = 0) ([docvis]private + [docvis]chronic = 1)
test (_b[female] = 0) (_b[private] + _b[chronic] = 1)

* Wald test of overall significance
test private chronic female income

* Manually compute overall test of significance using the formula for W
quietly poisson docvis private chronic female income, vce(robust)
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

* Confidence interval for linear combinations using lincom
use mus10data.dta, clear
quietly keep if year02==1
quietly poisson docvis private chronic female income if year02==1, vce(robust)
lincom private + chronic - 1

* Confidence interval for nonlinear function of parameters using nlcom
nlcom _b[female] / _b[private] - 1

* CI for exp(b) using lincom option eform
lincom private, eform

* CI for exp(b) using lincom followed by exponentiate
lincom private

* CI for exp(b) using nlcom
nlcom exp(_b[private])

********** 12.4 LIKELIHOOD RATIO TESTS

* LR tests output if estimate by ML with default estimate of VCE
nbreg docvis private chronic female income,  nolog

// Following not included in book
testnl exp([lnalpha]_cons)=0

* LR test using command lrtest
quietly nbreg docvis private chronic female income
estimates store unrestrict
quietly nbreg docvis female income
estimates store restrict
lrtest unrestrict restrict

// Following not included in book
lrtest restrict unrestrict

* Wald test of the same hypothesis
quietly nbreg docvis private chronic female income
test chronic private

// Following not included in book as can instead use option force
* LR test done manually
quietly nbreg docvis private chronic female income
scalar llunrestrict = e(ll)
quietly poisson docvis private chronic female income
scalar llrestrict = e(ll)
scalar lr = -2*(llrestrict - llunrestrict)
display "LR = " lr " with p = " chi2tail(1,lr)

* LR test using option force
quietly nbreg docvis private chronic female income
estimates store nb
quietly poisson docvis private chronic female income
estimates store poiss
lrtest nb poiss, force
display "Corrected p-value for LR-test = " r(p)/2

********** 12.5 LAGRANGE MULTIPLIER TEST (OR SCORE TEST)

* Perform LM test that b_private=0, b_chronic=0 using auxiliary regression
use mus10data.dta, clear
quietly keep if year02==1
generate one = 1
constraint define 1 private = 0
constraint define 2 chronic = 0
quietly nbreg docvis female income private chronic, constraints(1 2)
predict eqscore ascore, scores
generate s1restb = eqscore*one
generate s2restb = eqscore*female
generate s3restb = eqscore*income
generate s4restb = eqscore*private
generate s5restb = eqscore*chronic
generate salpha = ascore*one
quietly regress one s1restb s2restb s3restb s4restb s5restb salpha, noconstant
scalar lm = e(N)*e(r2)
display "LM = N x uncentered Rsq = " lm " and p = " chi2tail(2,lm)

********** 12.6 SIZE AND TEST POWER

// Following not included in book
* Generate one sample from DGP with b2 = 0 and estimate by OLS
clear
set seed 10101
quietly set obs 25
gen x = 12 + rnormal(0,2)
quietly generate y = 10 + 0.06*x + rnormal(0,0.3)
summarize x y
regress y x, noheader

* Do 1000 simulations where each gets p-value of test of b2=2
set seed 10101
postfile sim pvalues using pvalues, replace
forvalues i = 1/1000 {
  drop _all
  quietly set obs 150
  quietly generate double x = rchi2(1)
  quietly generate y = 1 + 2*x + rchi2(1)-1
  quietly regress y x
  quietly test x = 2
  scalar p = r(p)         // p-value for test this simulation
  post sim (p)
}
postclose sim

* Summarize the p-value from each of the 1000 tests
use pvalues, clear 
summarize pvalues

* Determine size of test at level 0.05
count if pvalues < .05
display "Test size from 1000 simulations = " r(N)/1000

* 95% simulation interval using exact binomial at level 0.05 with S=1000
cii 1000 50

// Following not included in book
cii 10000 500

* Program to compute power of test given specified H0 and Ha values of b2
program power, rclass
   version 10.1
   args numsims numobs b2H0 b2Ha nominalsize
                                        // Setup before iteration loops
   drop _all
   set seed 10101
   postfile sim pvalues using power, replace
                                        // Simulation loop
   forvalues i = 1/`numsims' {
     drop _all
     quietly set obs `numobs'
     quietly generate double x = rchi2(1)
     quietly generate y = 1 + `b2Ha'*x + rchi2(1)-1
     quietly regress y x
     quietly test x = `b2H0'
     scalar p = r(p)
     post sim (p)
   }
   postclose sim
   use power, clear
                                         // Determine the size or power
   quietly count if pvalues < `nominalsize'
   return scalar power=r(N)/`numsims'
end

* Size = power of test of b2H0=2 when b2Ha=2, S=1000, N=150, alpha=0.05
power 1000 150 2.00 2.00 0.05
display r(power) " is the test power"

* Power of test of b2H0=2 when b2Ha=2.2, S=1000, N=150, alpha=0.05
power 1000 150 2.00 2.20 0.05
display r(power) " is the test power"

* Power of test of H0:b2=2 against Ha:b2=1.6,1.625, ..., 2.4
postfile simofsims b2Ha power using simresults, replace
forvalues i = 0/33 {
  drop _all
  scalar b2Ha = 1.6 + 0.025*`i'
  power 1000 150 2.00 b2Ha 0.05
  post simofsims (b2Ha) (r(power))
}
postclose simofsims
use simresults, clear
summarize

* Plot the power curve
twoway (connected power b2Ha), scale(1.2) plotregion(style(none))
quietly graph export mus12power.eps, replace

* Power of chi(1) test when noncentrality parameter lambda = 5.67
display 1-nchi2(1,5.67,3.84)

********** CLOSE OUTPUT **********
