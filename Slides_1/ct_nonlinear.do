* ct_nonlinear.do  April 2013 for Stata version 12.0

********** OVERVIEW OF ct_nonlinear.do **********

log using ct_nonlinear.txt, text replace

* STATA Program 
* For A. Colin Cameron and Pravin Trivedi "Lectures in Microeconometrics"
* Nonlinear regression (M-estimation)

* To run you need file
*   mus10data.dta    
* in your directory

********** SETUP **********

set more off
version 12.0
set mem 10m
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* 2002 Medical Expenditure Panel Survey (MEPS)
* U.S. individuals aged 25-64 years working in private sector 
* but not self-employed
* and not receiving public insurance (Medicare and Medicaid) 
* Data due to Deb, Munkin and Trivedi (2006)

******* NONLINEAR REGRESSION - POISSON EXAMPLE

* Read in data, select 2002 data, describe and summarize key variables
use mus10data.dta, clear
quietly keep if year02==1
describe docvis private chronic female income
summarize docvis private chronic female income
histogram docvis if docvis < 40, width(1)

* Poisson regression with default standard errors
poisson docvis private chronic female income

* Poisson regression with robust standard errors
poisson docvis private chronic female income, vce(robust)

* Comparison of standard errors
quietly poisson docvis private chronic female income
estimates store DEFAULT
quietly poisson docvis private chronic female income, vce(robust)
estimates store ROBUST
estimates table DEFAULT ROBUST, b(%9.4f) se(%9.3f) stats(N r2 F)

* Marginal effects AME
margins, dydx(*)
* Marginal effects MEM
margins, dydx(*) atmean

* Old commands superceded by margins
* mfx
* margeff

* Marginal effects using finite difference for binary regressors
quietly poisson docvis i.private i.chronic i.female income, vce(robust)
margins, dydx(*)

* Wald test of nonlinear hypothesis
quietly poisson docvis private chronic female income, vce(robust)
testnl _b[female]/_b[private]=1

* Delta method confidence interval 
nlcom _b[female]/_b[private] - 1

* Nonlinear least-squares regression (command nl)
generate one = 1
nl (docvis = exp({xb: private chronic female income one})), vce(robust) nolog

* Newton-Raphson in Mata for Poisson MLE
* Set up data and local macros for dependent variable and regressors
generate cons = 1
local y docvis
local xlist private chronic female income cons
* Mata commands for Poisson MLE NR iterations
mata:
  st_view(y=., ., "`y'")            // read in stata data to y and X
  st_view(X=., ., tokens("`xlist'")) 
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
  iter                              // num iterations
  cha                               // stopping criterion
  st_matrix("b",b')                 // pass results from Mata to Stata
  st_matrix("V",vb)                 // pass results from Mata to Stata
end
* Present results, nicely formatted using Stata command ereturn
matrix colnames b = `xlist'
matrix colnames V = `xlist'
matrix rownames V = `xlist'
ereturn post b V
ereturn display   

********** CLOSE OUTPUT
* log close
* clear
* exit


