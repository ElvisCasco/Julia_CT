* mus11p1opt.do  Oct 2009 for Stata version 11

cap log close

********** OVERVIEW OF mus11p1opt.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Chapter 11
* 11.2: NEWTON-RAPHSON METHOD
* 11.3: GRADIENT METHODS
* 11.4: ML COMMAND: METHOD LF
* 11.5: CHECKING THE PROGRAM
* 11.6: ML COMMAND: METHODS D0, D1, D2, LF0, LF1, and LF2 METHODS
* 11.7: MATA OPTIMIZE() FUNCTION
* 11.8: GENERALIZED METHOD OF MOMENTS

* To run you need files
*   mus10data.dta 
* in your directory
* No Stata user-written commands are used

********** SETUP **********

set more off
version 11
clear all
set linesize 84
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* MEPS data described in chapter 10 program

************ 11.2: NEWTON-RAPHSON METHOD IN MATA

* Estimate model demonstrating some of the maximize options
* poisson docvis $xlist, robust

* Core Mata code for Poisson MLE NR iterations 
mata
  cha = 1                           // initialize stopping criterion
  do {
     mu = exp(X*b)                
     grad = X'(y-mu)                // kx1 gradient vector
     hes = cross(X, mu, X)          // negative of the kxk hessian matrix
     bold = b
     b = bold + cholinv(hes)*grad  
     cha = (bold-b)'(bold-b)/(bold'bold)  
     iter = iter + 1
  } while (cha > 1e-16)             // end of iteration loops 
end
*/

* Set up data and local macros for dependent variable and regressors
use mus10data.dta
keep if year02 == 1
generate cons = 1
local y docvis
local xlist private chronic female income cons

* Complete Mata code for Poisson MLE NR iterations
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
  iter                              // number of iterations
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

// Following not included in book is for comparison with above
poisson docvis private chronic female income, log trace gradient iterate(100) tol(1e-4) ltol(0) robust

************ 11.3: GRADIENT METHODS

* Objective function with multiple optima
graph twoway function                                        ///
  y=100-0.0000001*(x-10)*(x-30)*(x-50)*(x-50)*(x-70)*(x-80), ///
  range (5 90) plotregion(style(none))                       ///
  title("Objective function Q(theta) as theta varies")       ///
  xtitle("Theta", size(medlarge)) xscale(titlegap(*5))       ///
  ytitle("Q(theta)", size(medlarge)) yscale(titlegap(*5))
graph export mus11fig1.eps, replace


************ 11.4: ML COMMAND: METHOD LF

* Poisson ML program lfpois to be called by command ml method lf
program lfpois
  version 11
  args lnf theta1                  // theta1=x'b, lnf=ln(y)
  tempvar lnyfact mu
  local y "$ML_y1"                 // Define y so program more readable
  generate double `lnyfact' = lnfactorial(`y')
  generate double `mu'      = exp(`theta1')
  quietly replace `lnf'     = -`mu' + `y'*`theta1' - `lnyfact'
end

* Command ml model including defining y and x, plus ml check
ml model lf lfpois (docvis = private chronic female income), vce(robust)
ml check

* Search for better starting values
ml search

* Compute the estimator
ml maximize  

* Negbin ML program lfnb to be called by command ml method lf
program lfnb
  version 11
  args lnf theta1 a               // theta1=x'b, a=alpha, lnf=ln(y)
  tempvar mu
  local y $ML_y1                  // Define y so program more readable
  generate double `mu'  = exp(`theta1')
  quietly replace `lnf' = lngamma(`y'+(1/`a')) - lngamma((1/`a'))  ///
               -  lnfactorial(`y') - (`y'+(1/`a'))*ln(1+`a'*`mu')  ///
               + `y'*ln(`a') + `y'*ln(`mu') 
end

* Command lf implemented for negative binomial MLE
ml model lf lfnb (docvis = private chronic female income) ()
ml maximize, nolog

* NLS program lfnls to be called by command ml method lf
program lfnls
  version 11
  args lnf theta1                 // theta1=x'b, lnf=squared residual
  local y "$ML_y1"                // Define y so program more readable
  quietly replace `lnf' = -(`y'-exp(`theta1'))^2
end

* Command lf implemented for NLS estimator 
ml model lf lfnls (docvis = private chronic female income), vce(robust)
ml maximize

************ 11.5: CHECKING THE PROGRAM

* Example with high collinearity interpreted as perfect collinearity
generate extra = income + 0.001*runiform()
ml model lf lfpois (docvis = private chronic female income extra), vce(robust)
ml maximize

* Example with high collinearity not interpreted as perfect collinearity
generate extra2 = income + 0.01*runiform()
ml model lf lfpois (docvis = private chronic female income extra2), vce(robust)
ml maximize, nolog

* Detect multicollinearity using _rmcoll
_rmcoll income extra 
_rmcoll income extra2

* Generate dataset from Poisson DGP for large N
clear
set obs 10000               
set seed 10101              
generate x = rnormal(0,0.5)
generate mu = exp(2 + x)   
generate y = rpoisson(mu)
summarize mu x y            

* Consistency check: Run program lfpois and compare beta to DGP value
ml model lf lfpois (y = x)  
ml maximize, nolog

* Program to generate dataset, obtain estimate, and return beta and SEs
program secheck, rclass
  version 11 
  drop _all 
  set obs 500 
  generate x = rnormal(0,0.5)
  generate mu = exp(2 + x)   
  generate y = rpoisson(mu)
  ml model lf lfpois (y = x)
  ml maximize
  return scalar b1 =_b[_cons]
  return scalar se1 = _se[_cons]
  return scalar b2 =_b[x]
  return scalar se2 = _se[x]
end 

* Standard errors check: run program secheck
set seed 10101    
simulate "secheck" bcons=r(b1) se_bcons=r(se1) bx=r(b2) se_bx=r(se2), reps(2000) 
summarize

************ 11.6: ML COMMAND: METHODS D0, D1, D2

* Set up data and local macros for dependent variable and regressors
use mus10data.dta, clear
quietly keep if year02 == 1
generate cons = 1
local y docvis
local xlist private chronic female income cons

* Method d0: Program d0pois to be called by command ml method d0
program d0pois
  version 11
  args todo b lnf                  // todo is not used, b=b, lnf=lnL 
  tempvar theta1                   // theta1=x'b given in eq(1)
  mleval `theta1' = `b', eq(1)   
  local y $ML_y1                   // Define y so program more readable
  mlsum `lnf' = -exp(`theta1') + `y'*`theta1' - lnfactorial(`y')
end

* Method d0: implement Poisson MLE 
ml model d0 d0pois (docvis = private chronic female income)
ml maximize

* Method d1: Program d1pois to be called by command ml method d1
program d1pois
  version 11
  args todo b lnf g                // gradient g added to the arguments list
  tempvar theta1                   // theta1 = x'b given in eq(1)
  mleval `theta1' = `b', eq(1)   
  local y $ML_y1                   // Define y so program more readable
  mlsum `lnf' = -exp(`theta1') + `y'*`theta1' - lnfactorial(`y')
  if (`todo'==0 | `lnf'>=.) exit   // Extra code from here on
  tempname d1
  mlvecsum `lnf' `d1' = `y' - exp(`theta1')
  matrix `g' = (`d1')
end

// * ml check
* Method d1: implement Poisson MLE
ml model d1 d1pois (docvis = private chronic female income)  // Define y and x
ml maximize

* Method lf1: Program lf1poisrob is variation of program d1pois for robust se's
program lf1poisrob
  version 11
  args todo b lnfj g1 
  tempvar theta1                   // theta1 = x'b where x given in eq(1)
  mleval `theta1' = `b', eq(1)   
  local y $ML_y1                   // define y so program more readable
  quietly replace `lnfj' = -exp(`theta1') + `y'*`theta1' - lnfactorial(`y')
  if (`todo'==0) exit  
  quietly replace `g1' = `y' - exp(`theta1')  // extra code for robust 
end

* Method lf1: implement Poisson MLE with robust standard errors
ml model lf1 lf1poisrob (docvis = private chronic female income), vce(robust)
ml maximize, nolog

* Method d2: Program d2pois to be called by command ml method d2
program d2pois
  version 11
  args todo b lnf g H              // Add g and H to the arguments list
  tempvar theta1                   // theta1 = x'b where x given in eq(1)
  mleval `theta1' = `b', eq(1)     
  local y $ML_y1                   // Define y so program more readable
  mlsum `lnf' = -exp(`theta1') + `y'*`theta1' - lnfactorial(`y')
  if (`todo'==0 | `lnf'>=.) exit   // d1 extra code from here
  tempname d1
  mlvecsum `lnf' `d1' = `y' - exp(`theta1')
  matrix `g' = (`d1')
  if (`todo'==1 | `lnf'>=.) exit   // d2 extra code from here
  tempname d11
  mlmatsum `lnf' `d11' = exp(`theta1')
  matrix `H' = -`d11'
end

* Method d2: Poisson MLE with first and second derivatives provided
ml model d2 d2pois (docvis = private chronic female income)
ml maximize

************ 11.7: MATA OPTIMIZE() FUNCTION

local y docvis
local xlist private chronic female income cons

mata: mata set matastrict off
* Evaluator function for Poisson MLE using optimize gf2 evaluator
mata 
  void poissonmle(todo, b, y, X, lndensity, g, H)
  {
    Xb = X*b'
    mu = exp(Xb)
    lndensity = -mu + y:*Xb - lnfactorial(y)
    if (todo == 0) return
    g = (y-mu):*X   
    if (todo == 1) return
    H = - cross(X, mu, X)
  }
end

mata: mata drop poissonmle()
* Mata code to obtain Poisson MLE using command optimize
mata 
  void poissonmle(todo, b, y, X, lndensity, g, H)
  {
    Xb = X*b'
    mu = exp(Xb)
    lndensity = -mu + y:*Xb - lnfactorial(y)
    if (todo == 0) return
    g = (y-mu):*X   
    if (todo == 1) return
    H = - cross(X, mu, X)
  }
  st_view(y=., ., "`y'")  
  st_view(X=., ., tokens("`xlist'"))
  S = optimize_init()
  optimize_init_evaluator(S, &poissonmle())
  optimize_init_evaluatortype(S, "gf2")
  optimize_init_argument(S, 1, y)
  optimize_init_argument(S, 2, X)
  optimize_init_params(S, J(1,cols(X),0))
  b = optimize(S)  
  Vbrob = optimize_result_V_robust(S)
  serob = (sqrt(diagonal(Vbrob)))'
  b \ serob 
end

************ 11.8 GENERALIZED METHOD OF MOMENTS

mata: mata set matastrict off

local y docvis
local xlist private chronic female income cons
local zlist firmsize chronic female income cons

/*
mata
    Xb = X*b'                    // b for optimize is 1xk row vector
    mu = exp(Xb)
    h = (Z'(y-mu)                // h is rx1 column row vector
    W = cholinv(Z'Z)             // W is rxr wmatrix
    G = -(mu:*Z)'X               // G is rxk matrix
    S = ((y-mu):*Z)'((y-mu):*Z)  // S is rxr matrix
    Qb = h'W*h                   // Q(b) is scalar
    g = G'W*h                    // gradient for optimize is 1xk row vector
    H = G'W*G                    // hessian for optimize is kxk matrix
    V = luinv(G'W*G)*G'W*S*W*G*luinv(G'W*G)  
end
*/

* Mata code to obtain GMM estimator for Poisson using command optimize
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
  // Compute robust estimate of VCE and SEs
  Xb = X*b'
  mu = exp(Xb)
  h = Z'(y-mu)
  W = cholinv(cross(Z,Z))
  G = -(mu:*Z)'X
  Shat = ((y-mu):*Z)'((y-mu):*Z)*rows(X)/(rows(X)-cols(X))
  Vb = luinv(G'W*G)*G'W*Shat*W*G*luinv(G'W*G)
  seb = (sqrt(diagonal(Vb)))'
  b \ seb 
end

********** CLOSE OUTPUT
