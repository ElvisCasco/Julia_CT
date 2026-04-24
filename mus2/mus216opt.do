* mus216opt.do  for Stata 17

capture log close

********** OVERVIEW OF mus216opt.do **********

* Stata program
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi
* used for "Microeconometrics Using Stata, Second Edition"
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 16
* 16.2: NEWTON-RAPHSON METHOD
* 16.3: GRADIENT METHODS
* 16.5: THe ML COMMAND: LF METHOD
* 16.6: CHECKING THE PROGRAM
* 16.7: THE ML COMMAND: LF0-LF2, D0-D2, AND GF0 METHODS
* 16.8: NONLINEAR INSTRUMENTAL VARIABLES (GMM) EXAMPLE

* To run you need files
*   mus210mepsdocvisyoung.dta
* in your directory

* No community-contributed commands are used

* The simulations at the end of section 16.6 take a long time
* You may want to comment out

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus210mepsdocvisyoung.dta authors' extracvt from from 2002 MEPS (
* (Medical Expenditure Panel Survey)
* U.S. data on office-based physician visits by persons aged 25-64 years
* Same as Deb, Munkin and Trivedi, "Bayesian analysis of the two-part model
* with endogeneity: application to health care expenditure", JBES, 1081-1099
* Excludes those receiving public insurance (Medicare and Medicaid)
* Restricted to those working in the private sector but not self-employed.

************ 16.2: NEWTON-RAPHSON METHOD IN MATA

* Estimate model demonstrating some of the maximize options
* poisson docvis $xlist, robust

/* Following included as text in book
* Core Mata code for Poisson MLE NR iterations
mata
    p      = cols(X)                  // Number of regressors
    b      = J(p, 1, .02)             // Starting values are close to zero
    gHinvg = 1                        // Initialize scaled gradient 	
    iter   = 1                        // Initialize number of iterations
    do {
        mu   = exp(X*b)
        grad = X'*(y-mu)              // k x 1 gradient vector
        hes  = cross(X, mu, X)        // Negative of the kxk Hessian matrix
        diff = cholinv(hes)*grad      // Update amount
        bold = b
        b    = bold + diff
        iter = iter + 1
        gHinvg = grad'*diff           // Equals grad'*cholinv(hes)*grad
		printf("iter = %s, gHginv = %6.5g\n", strofreal(iter), gHinvg)
    } while (gHinvg > 1e-8)           // End of iteration loops
end
*/

* Set up data and local macros for dependent variable and regressors
qui use mus210mepsdocvisyoung
keep if year02 == 1
generate cons = 1
global y docvis
global xlist private chronic female income cons

* Complete Mata code for Poisson MLE NR iterations
mata:
    st_view(y=., ., "$y")             // Read in stata data to y and X
    st_view(X=., ., tokens("$xlist"))
    p      = cols(X)                  // Number of regressors
    n      = rows(X)
    b      = J(p, 1, .02)             // Start values are close to zero
    gHinvg = 1                        // Initialize scaled gradient 	
    iter   = 0                        // Initialize number of iterations
    do {
        mu   = exp(X*b)
        grad = X'*(y-mu)              // k x 1 gradient vector
        hes  = cross(X, mu, X)        // Negative of the kxk Hessian matrix
        diff = cholinv(hes)*grad      // Update amount
        bold = b
        b    = bold + diff
        iter = iter + 1
        gHinvg = grad'*diff           // Equals grad'*cholinv(hes)*grad
        printf("iter = %s, gHginv = %6.5g\n", strofreal(iter), gHinvg)
    } while (gHinvg > 1e-8)           // End of iteration loops
    mu = exp(X*b)
    hes = cross(X, mu, X)
    vgrad = cross(X, (y-mu):^2, X)
    vb = cholinv(hes)*vgrad*cholinv(hes)*n/(n-cols(X))
    st_matrix("b",b')                 // Pass results from Mata to Stata
    st_matrix("V",vb)                 // Pass results from Mata to Stata
end

* Present results, nicely formatted using Stata command ereturn
matrix colnames b = $xlist
matrix colnames V = $xlist
matrix rownames V = $xlist
ereturn post b V
ereturn display

// Not included - Example of various maximize options
poisson docvis private chronic female income, log trace gradient ///
  iterate(100) tol(1e-4) ltol(0) robust

************ 16.3: GRADIENT METHODS

* Objective function with multiple optima
graph twoway function                                          ///
    y=100-0.0000001*(x-10)*(x-30)*(x-50)*(x-50)*(x-70)*(x-80), ///
    range (5 90) plotregion(style(none)) scale(1.2)            ///
    title("Objective function Q({&theta}) as {&theta} varies") ///
    xtitle("{&theta}", size(medlarge)) xscale(titlegap(*5))    ///
    ytitle("Q({&theta})", size(medlarge)) yscale(titlegap(*5))

************ 16.5: THE ML COMMAND: METHOD LF

* Poisson ML program lfpois to be called by command ml method lf
program lfpois
    version 17
    args lnf theta1                  // theta1=x'b, lnf=ln(y)
    tempvar lnyfact mu
    local y "$ML_y1"                 // Define y so program more readable
    generate double `lnyfact' = lnfactorial(`y')
    generate double `mu'      = exp(`theta1')
    qui replace `lnf'     = -`mu' + `y'*`theta1' - `lnyfact'
end

* Following command only part of outpur included in book
* Command ml model defines y and x with heteroskedastic robust standard errors plus ml check
ml model lf lfpois (docvis = private chronic female income), vce(robust)
ml check

* Search for better starting values
ml search

* Command lfpois implemented for Poisson MLE
ml maximize

* Following in text but with output omitted
* Same model but with cluster–robust standard errors
ml model lf lfpois (docvis = private chronic female income), vce(cluster age)
ml maximize

* Negbin ML two-index program lfnb to be called by command ml method lf
program lfnb
    version 17
    args lnf theta1 a               // theta1=x'b, a=alpha, lnf=ln(y)
    tempvar mu
    local y $ML_y1                  // Define y so program more readable
    generate double `mu'  = exp(`theta1')
    qui replace `lnf' = lngamma(`y'+(1/`a')) - lngamma((1/`a'))    ///
               -  lnfactorial(`y') - (`y'+(1/`a'))*ln(1+`a'*`mu')  ///
               + `y'*ln(`a') + `y'*ln(`mu')
end

* Command lf implemented for negative binomial MLE
ml model lf lfnb (docvis = private chronic female income) (), vce(robust)
ml maximize, nolog

* NLS program lfnls to be called by command ml method lf
program lfnls
    version 17
    args lnf theta1                 // theta1=x'b, lnf=squared residual
    local y "$ML_y1"                // Define y so program more readable
    qui replace `lnf' = -(`y'-exp(`theta1'))^2
end

* Following in text but with output omitted
* Command lf implemented for NLS estimator
ml model lf lfnls (docvis = private chronic female income), vce(robust)
ml maximize

************ 16.6: CHECKING THE PROGRAM

* Example with high collinearity interpreted as perfect collinearity
generate extra = income + 0.001*runiform()
ml model lf lfpois (docvis = private chronic female income extra), vce(robust)
ml maximize

* Example with high collinearity not interpreted as perfect collinearity
generate extra2 = income + 0.01*runiform()
ml model lf lfpois (docvis = private chronic female income extra2), vce(robust)
ml maximize, nolog

* Detect multicollinearity by using _rmcoll
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

// The following takes a long time - may want to comment out

* Program to generate dataset, obtain estimate, and return beta and SEs
program secheck, rclass
    version 17
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

* Standard errors check: Run program secheck
set seed 10101
simulate "secheck" bcons=r(b1) se_bcons=r(se1) bx=r(b2) se_bx=r(se2), reps(2000)
summarize

/* The following takes a long time - so comment out

* Alternative using postfile
set seed 10101
postfile simsecheck bcons se_bcons bx se_bx using simresults.dta, replace
forvalues i = 1/2000 {
    qui {
        drop _all
        set obs 500
        generate x = rnormal(0,0.5)
        generate mu = exp(2 + x)
        generate y = rpoisson(mu)
        ml model lf lfpois (y = x)
        ml maximize
        post simsecheck (_b[_cons]) (_se[_cons]) (_b[x]) (_se[x])
    }
}
postclose simsecheck
use simresults.dta, clear
summarize

*/

************ 16.7: THE ML COMMAND: METHODS LF0, LF1, LF2


* ml method lf0: Program lf0pois gives lnf(y_i) in terms of index x_i'b
program lf0pois
    version 17
    args todo b lnfi
    tempvar theta1                   // theta1 = x'b where x given in eq(1)
    mleval `theta1' = `b', eq(1)
    local y $ML_y1                   // Define y so program more readable
    qui replace `lnfi' = -exp(`theta1') + `y'*`theta1' - lnfactorial(`y')
end

qui use mus210mepsdocvisyoung, clear
qui keep if year02 == 1
generate cons = 1
* ml method lf0: Obtain Poisson MLE with heteroskedastic-robust standard errors
ml model lf0 lf0pois (docvis = private chronic female income), vce(robust)
ml maximize

* Check
poisson docvis private chronic female income, nolog vce(robust)

* Not included Check cluster–robust
ml model lf0 lf0pois (docvis = private chronic female income), vce(cluster age)
ml maximize, nolog
poisson docvis private chronic female income, nolog vce(cluster age)

* ml method lf1: Program lf1pois adds analytical first derivatives
program lf1pois
    version 17
    args todo b lnfi gi
    tempvar theta1                   // theta1 = x'b where x given in eq(1)
    mleval `theta1' = `b', eq(1)
    local y $ML_y1                   // Define y so program more readable
    qui replace `lnfi' = -exp(`theta1') + `y'*`theta1' - lnfactorial(`y')
    if (`todo'==0) exit
    qui replace `gi' = `y' - exp(`theta1')  // Extra code for robust
end

* ml method lf1: Implement Poisson MLE with cluster–robust standard errors
ml model lf1 lf1pois (docvis = private chronic female income), vce(cluster age)
ml maximize, nolog

* Check
poisson docvis private chronic female income, nolog vce(cluster age)

* ml method lf2: Program lf1pois adds analytical first derivatives
program lf2pois
    version 17
    args todo b lnfi gi H
    tempvar theta1                   // theta1 = x'b where x given in eq(1)
    mleval `theta1' = `b', eq(1)
    local y $ML_y1                   // Define y so program more readable
    qui replace `lnfi' = -exp(`theta1') + `y'*`theta1' - lnfactorial(`y')
    if (`todo'==0) exit
    qui replace `gi' = `y' - exp(`theta1')  // Extra code for robust
    tempname d11
    mlmatsum `lnfi' `d11' = exp(`theta1')
    matrix `H' = -`d11'
end

* ml method lf2: Implement Poisson MLE with cluster–robust standard errors
ml model lf2debug lf2pois (docvis = private chronic female income), vce(cluster age)
ml check
ml maximize, nolog

* Check
poisson docvis private chronic female income, nolog vce(cluster age)

**** COMMAND ML: METHODS D0, D1, D2

* ml method d0: Program lf0pois gives lnf(y) in terms of index x_i'b
program d0pois
    version 17
    args todo b lnf                  // todo is not used, b=b, lnf=lnL
    tempvar theta1                   // theta1=x'b given in eq(1)
    mleval `theta1' = `b', eq(1)
    local y $ML_y1                   // Define y so program more readable
    mlsum `lnf' = -exp(`theta1') + `y'*`theta1' - lnfactorial(`y')
end

* ml method d0: Poisson MLE with default standard errors only possible
ml model d0 d0pois (docvis = private chronic female income)
ml maximize, nolog

* ml method d2: Program d2pois adds analytical first and second derivatives
program d2pois
    version 17
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

* ml method d2: Implement Poisson MLE with default standard errors only possible
ml model d2 d2pois (docvis = private chronic female income)
ml maximize, nolog

**** COMMAND ML: METHOD GF0

* ml method gf0: Program gf0pois gives lnf(y) in terms of index x_i'b
program gf0pois
    version 17
    args todo b lnf                  // todo is not used, b=b, lnf=lnL
    tempvar theta1                   // theta1=x'b given in eq(1)
    mleval `theta1' = `b', eq(1)
    local y $ML_y1                   // Define y so program more readable
    qui replace `lnf' = -exp(`theta1') + `y'*`theta1' - lnfactorial(`y')
end

* ml method gf0: Implement Poisson MLE with heteroskedastic-robust standard errors
ml model gf0 gf0pois (docvis = private chronic female income), vce(robust)
ml maximize, nolog

************ 16.8: NONLINEAR IV EXAMPLE

/* Following included as text in book
* Core Mata code for GMM example
mata
    Xb = X*b'                     // b for optimize is 1 x k row vector
    mu = exp(Xb)
    h  = (Z'(y-mu)                // h is r x 1 column row vector
    W  = cholinv(Z'Z)             // W is r x r wmatrix
    G  = -(mu:*Z)'X               // G is r x k matrix
    S  = ((y-mu):*Z)'((y-mu):*Z)  // S is r x r matrix
    Qb = h'W*h                    // Q(b) is scalar
    g  = G'W*h                    // Gradient for optimize is 1 x k row vector
    H  = G'W*G                    // Hessian for optimize is k x k matrix
    V  = luinv(G'W*G)*G'W*S*W*G*luinv(G'W*G)
end
*/

mata: mata set matastrict off

global y docvis
global xlist private chronic female income cons
global zlist firmsize chronic female income cons

* optimize() method d2: Evaluator and implement GMM estimator for Poisson
mata
    mata clear
    void pgmmd2(todo, b, y, X, Z, Qb, g, H)
    {
        Xb = X*b'
        mu = exp(Xb)
        h  = Z'(y-mu)
        W  = cholinv(cross(Z,Z))
        Qb = h'W*h
        if (todo == 0) return
        G  = -(mu:*Z)'X
        g  = (G'W*h)'
        if (todo == 1) return
        H = G'W*G
        _makesymmetric(H)
     }
    st_view(y=., ., "$y")
    st_view(X=., ., tokens("$xlist"))
    st_view(Z=., ., tokens("$zlist"))
    S = optimize_init()
    optimize_init_which(S,"min")
    optimize_init_evaluator(S, &pgmmd2())
    optimize_init_evaluatortype(S, "d2")
    optimize_init_argument(S, 1, y)
    optimize_init_argument(S, 2, X)
    optimize_init_argument(S, 3, Z)
    optimize_init_params(S, J(1,cols(X),0))
    optimize_init_technique(S,"nr")
    b = optimize(S)
    // Compute robust estimate of VCE and SEs
    Xb   = X*b'
    mu   = exp(Xb)
    h    = Z'(y-mu)
    W    = cholinv(cross(Z,Z))
    G    = -(mu:*Z)'X
    n  = rows(X)
    k  = cols(X)
    Shat = ((y-mu):*Z)'((y-mu):*Z)*rows(n)/(n-k)
    Vb   = luinv(G'W*G)*G'W*Shat*W*G*luinv(G'W*G)
    seb  = (sqrt(diagonal(Vb)))'
    b \ seb
end

*** Erase files created by this program and not used elsewhere in the book

********** END
