* mus2appendix.do  for Stata 17

cap log close

********** OVERVIEW OF mus2appendix.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Appendix A programs
*   A.1 STATA MATRIX COMMANDS
*   A.2 PROGRAMS
*   A.3 PROGRAM DEBUGGING
*   B.1 HOW TO RUN MATA
*   B.2 MATA MATRIX COMMANDS
*   B.3 PROGRAMMING IN MATA
*   C.1: MATA MOPTIMIZE() FUNCTION
*   C.2: MATA OPTIMIZE() FUNCTION

* To run you need file
*    mus203mepsmedexp.dta
*    mus210mepsdocvisyoung.dta
* in your directory

* No community-contributed commands are needed

********** SETUP **********

clear all
mata: mata set matastrict off
set linesize 82
set scheme s1mono  /* Graphics scheme */

**************** A.1 STATA MATRIX COMMANDS

* Define a matrix explicitly and list the matrix
matrix define A = (1,2,3 \ 4,5,6)
matrix list A

* Matrix row and column names
matrix rownames A = one two
matrix list A

* Read in data, summarize, and run regression
use mus203mepsmedexp
keep if _n <= 100
drop if ltotexp == . | totchr == .
summarize ltotexp totchr  
regress ltotexp totchr, noheader 

* Create a matrix from estimation results
matrix vbols = e(V)
matrix list vbols

* Create a matrix from variables
mkmat ltotexp, matrix(y)
generate intercept = 1
mkmat totchr intercept, matrix(X)

* Convert a matrix into variables
svmat X
summarize X* 

* Change value of an entry in matrix
matrix A[1,1] = A[1,2]
matrix list A

* Select part of matrix
matrix B = A[1...,2..3]
matrix list B

* Add columns to an existing matrix
matrix C = B, B
matrix list C

* Matrix operators
matrix D = C + 3*C
matrix list D

* Matrix functions
matrix r = rowsof(D)
matrix list r

* Can use scalar if 1 x 1 matrix
scalar ralt = rowsof(D)
display ralt

* Inverse of nonsymmetric square matrix
matrix Binv = inv(B)
matrix list Binv

* OLS estimator using X and y and matrix operators
matrix bols = (invsym(X'*X))*(X'*y)
matrix list bols

* OLS estimator using matrix accumulation operators
matrix accum XTX = totchr              // Form X'X including constant
matrix vecaccum yTX = ltotexp totchr   // Form y'X including constant
matrix cols = invsym(XTX)*(yTX)'
matrix list cols

* Illustrate Stata matrix commands: OLS with output
matrix accum XTX = totchr            // Form X'X including constant
matrix vecaccum yTX = ltotexp totchr // Form y'X including constant
matrix b = invsym(XTX)*(yTX)'
matrix accum yTy = ltotexp, noconstant
scalar k = rowsof(XTX)
scalar n = _N
matrix s2 = (yTy - b'*XTX'*b)/(n-k)
matrix V = s2*invsym(XTX)
matrix list b
matrix list V

* Stata matrix commands to compute standard errors and t statistics given b and V
matrix se = (vecdiag(cholesky(diag(vecdiag(V)))))'
matrix seinv = (vecdiag(cholesky(invsym(diag(vecdiag(V))))))'
matrix t = hadamard(b,seinv)
matrix results = b, se, t
matrix colnames results = coeff sterror tratio
matrix list results, format(%7.0g)

* Easier is to use ereturn post and display given b and V
matrix brow = b'
ereturn post brow V
ereturn display 

**************** A.2 STATA PROGRAMS

* Program with no arguments
program define time
  display c(current_time) c(current_date)
end

* Run the program
time

* Drop program if it already exists, write program, and run
capture program drop time
program time
  display c(current_time) c(current_date)
end
time
 
* Program with two positional arguments
program meddiff
    tempvar diff
    generate `diff' = `1' - `2'
    _pctile `diff', p(50)
    display "Median difference = " r(r1)
end

* Run the program with two arguments
meddiff ltotexp totchr

* Program with two named positional arguments
capture program drop meddiff
program meddiff
    args y x
    tempvar diff
    generate `diff' = `y' - `x'
    _pctile `diff', p(50)
    display "Median difference = " r(r1)
end
meddiff ltotexp totchr

* Program with results stored in r()
capture program drop meddiff
program meddiff, rclass
    args y x
    tempvar diff
    generate `diff' = `y' - `x'
    _pctile `diff', p(50)
    return scalar medylx = r(r1)
end

* Running the program does not immediately display the result
meddiff ltotexp totchr
return list
display r(medylx)

* Program that uses Stata commands syntax and gettoken to parse arguments
program myols
    syntax varlist [if] [in] [,vce(string)]
    gettoken y xvars : varlist
    display "varlist contains: "  "`varlist'" 
    display "and  if contains: "  "`if'" 
    display "and  in contains: "  "`in'"
    display "and vce contains: "  "`vce'"
    display "and   y contains: "  "`y'"
    display "& xvars contains: "  "`xvars'"
    regress `y' `xvars' `if' `in', `vce' noheader
end

* Execute program myols for an example 
myols ltotexp totchr if !missing(ltotexp) in 1/100, vce(robust)

capture program drop meddiff
*! version 2.1.0  15aug2021
program meddiff, rclass
    version 17
    args y x 
    tempvar diff
    qui {
        generate double `diff' = `y' - `x'
        _pctile `diff', p(50)
         return scalar medylx = r(r1)
    }
    display "Median of first variable - second variable = " r(r1)
end   

* Execute program meddiff for an example  
meddiff ltotexp totchr

**************** A.3 PROGRAM DEBUGGING 

* Display intermediate output to aid debugging
matrix accum XTX = totchr            // Recall constant is added
matrix list XTX                      // Should be 2 x 2 
matrix vecaccum yTX = ltotexp totchr
matrix list yTX                      // Should be 1 x 2
matrix bOLS = invsym(XTX)*(yTX)'
matrix list bOLS                     // Should be 2 x 1

* Debug an initial nonprogram version of a program
tempvar y x diff 
generate `y' = ltotexp
generate `x' = totchr
generate double `diff' = `y' - `x'
_pctile `diff', p(50)
scalar medylx = r(r1)
display "Median of first variable - second variable = " medylx

**************** B.1 HOW TO RUN MATA

* Read in data
clear all
use mus203mepsmedexp
keep if _n <= 100
generate cons = 1

* Sample Mata session
mata
I = I(2)  
I
end

* Mata commands issued from Stata
mata: I = I(2)
mata: I

// Stata commands issued from Mata
mata
stata("summarize ltotexp")
end

// Mata help
mata
// help mata det
// help m4 matrix
// help mata
end

****************   B.2 MATA MATRIX COMMANDS

mata  

// Create a matrix
A = (1,2,3 \ 4,5,6)

// List a matrix
A

// Create a 2 x 2 identity matrix
I = I(2)

// Create a 1 x 5 unit row vector with 1 in second entry and zeros elsewhere
e = e(2,5)
e

// Create a 2 x 5 matrix with entry 3 
J = J(2,5,3)
J

// Create a row vector with entries 8 to 15
a = 8..15
a

// Associate a Mata matrix with variables stored in Stata
st_view(y=., ., "ltotexp")
st_view(X=., ., ("totchr", "cons"))

// Create a Mata matrix from variables stored in Stata
Xloaded = st_data(., ("totchr", "cons"))

// Read Stata matrix (created in first line below) into Mata
stata("matrix define B = I(2)")
C = st_matrix("B")
C

// Transfer Mata matrix to Stata
st_matrix("D",C)
stata("matrix list D")

// Element by element multiplication of matrix by column vector
b = 2::3
J = J(2,5,3)
b:*J

// Matrix function that returns scalar
r = rows(A)
r

// Matrix function that returns matrix by element-by-element transformation
D = sqrt(A)
D

// Calculate eigenvalues and eigenvectors
E = (1, 2 \ 4, 3)
lambda = .
eigvecs = .
eigensystem(E,eigvecs,lambda)
lambda
eigvecs

// Matrix inversion: Use of makesymmetric() before cholinv()
F = 0.5*I(2)
G = makesymmetric(cholinv(F'F))
G

// Matrix cross product
beta = (cholinv(cross(X,X)))*(cross(X,y))
beta

// Matrix subscripts
A[1,2] = A[1,1]
A

// Combining matrices: add columns
M = A, A
M

// Combining matrices: add rows
N = A \ A
N

// Form submatrix using list subscripts
O = M[(1\2), (5::6)]
O

// Form submatrix using range subscripts
P = M[|1,5 \ 2,6|]
P

// Output mata matrix to Stata
st_matrix("Q", P)
stata("matrix list Q")

// Output mata matrix to Stata
yhat = X*beta
st_addvar("float", "ltotexphat")
st_store(.,"ltotexphat", yhat)
stata("summarize ltotexp ltotexphat")

end 

**************** B.3 MATA PROGRAMS

mata:
    void poissonmle(real scalar todo, 
    real rowvector b,
    real colvector y, 
    real matrix X,
    real colvector lndensity,
    real matrix g,
    real matrix H)
    {
        real colvector Xb
        real colvector mu
        Xb = X*b'
        mu = exp(Xb)
        lndensity = -mu + y:*Xb - lnfactorial(y)
        if (todo == 0) return
        g = (y-mu):*X   
        if (todo == 1) return
        H = - cross(X, mu, X)
    }
end

mata:
    void calcsum(varname, resultissum)  
    {
        st_view(x=., ., varname)
        resultissum = colsum(x)
    } 
    sum = .
    calcsum("ltotexp", sum)
    sum
end 

mata:
    void function calcsum2(varname)  
    {
        st_view(x=., ., varname)
        st_numscalar("r(sum)",colsum(x))
    } 
    calcsum2("ltotexp")
    stata("display r(sum)")
end 

program varsum
    version 17
    syntax varname
    mata: calcsum2("`varlist'")
    display r(sum)
end
varsum ltotexp

************ C.1: MATA COMMAND MOPTIMIZE

**** MOPTIMIZE: METHODS LF0, LF1, LF2

* Read in data
use mus210mepsdocvisyoung, clear
qui keep if year02 == 1

generate cons = 1

* moptimize() method lf: Program poissonlf gives lnf(y_i)  
mata:
    mata clear 
    function poissonlf(transmorphic M, real rowvector b, fv)
    {
        y1 = moptimize_util_depvar(M, 1)   // N x 1 vector y 
        Xb = moptimize_util_xb(M, b, 1)    // N x 1 vector Xb
        fv = -exp(Xb):+ (y1:*Xb) :- lnfactorial(y1)
    }
end

* Set up data and local macros for dependent variable and regressors
use mus210mepsdocvisyoung, clear 
keep if year02 == 1
generate cons = 1

* moptimize() method lf: Implement with default standard errors 
mata:
    M = moptimize_init()
    moptimize_init_evaluator(M, &poissonlf())
    moptimize_init_evaluatortype(M, "lf")
    moptimize_init_depvar(M, 1, "docvis")
    moptimize_init_eq_indepvars(M, 1, "private chronic female income")
    moptimize(M)
    moptimize_result_display(M)             // Default standard errors
end

* Check
poisson docvis private chronic female income, nolog

* moptimize() method lf: Implement with heteroskedastic-robust standard errors
mata:
    moptimize_result_display(M, "robust")   // Heteroskedastic-robust standard errors
end

* Check
poisson docvis private chronic female income, nolog  vce(robust)

* moptimize() method lf: Implement with cluster–robust standard errors
mata:
    M = moptimize_init()
    moptimize_init_evaluator(M, &poissonlf())
    moptimize_init_evaluatortype(M, "lf")
    moptimize_init_depvar(M, 1, "docvis")
    moptimize_init_eq_indepvars(M, 1, "private chronic female income")
    moptimize_init_cluster(M, "age")     // Define the cluster variable
    moptimize(M)
    moptimize_result_display(M)          // Now gives cluster–robust standard errors   
end

* Check
poisson docvis private chronic female income, nolog vce(cluster age)

* moptimize() method lf2: Add first and second derivatives to poissonlf  
mata:
    mata clear 
    function poissonlf2(transmorphic M, real scalar todo, ///
    real rowvector b, fv, S, H)
    {
        y1 = moptimize_util_depvar(M, 1)   // N x 1 vector y 
        Xb = moptimize_util_xb(M, b, 1)    // N x 1 vector Xb
        fv = -exp(Xb):+ (y1:*Xb) :- lnfactorial(y1)
        if (todo>=1) {
            s1 = (y1:-exp(Xb)) 
            S = s1
            if (todo==2) {
                h11 = -exp(Xb)
                H11 = moptimize_util_matsum(M, 1,1, h11, 0)
                H = H11
            }
        }
    }
end

* moptimize() method lf2: Implement with heteroskedastic-robust standard errors
mata:
M = moptimize_init()
moptimize_init_evaluator(M, &poissonlf2())
moptimize_init_evaluatortype(M, "lf2")
moptimize_init_depvar(M, 1, "docvis")
moptimize_init_eq_indepvars(M, 1, "private chronic female income")
moptimize(M)
moptimize_result_display(M, "robust")   // Robust standard errors
end

* Check
poisson docvis private chronic female income, vce(robust)

**** MOPTIMIZE: METHODS DF0, D1, D2

* moptimize() method d2: Program poissond2 
mata:
    mata clear 
    function poissond2(transmorphic M, real scalar todo, ///
        real rowvector b, fv, g, H)
    {
        y1 = moptimize_util_depvar(M, 1)    // N x 1 vector y1
        Xb = moptimize_util_xb(M, b, 1)     // N x 1 vector Xb
        fv = moptimize_util_sum(M, -exp(Xb):+ (y1:*Xb) :- lnfactorial(y1))
        if (todo>=1) {
            s1 = (y1:-exp(Xb)) 
            g1 = moptimize_util_vecsum(M, 1, s1, fv)
            g = g1
            if (todo==2) {
                h11 = -exp(Xb)
                H11 = moptimize_util_matsum(M, 1,1, h11, fv)
                H = H11
            }
        }
    }
end

* moptimize() method d2: Can implement only with default standard errors
mata:
    M = moptimize_init()
    moptimize_init_evaluator(M, &poissond2())
    moptimize_init_evaluatortype(M, "d2")
    moptimize_init_depvar(M, 1, "docvis")
    moptimize_init_eq_indepvars(M, 1, "private chronic female income")
    moptimize(M)
    moptimize_result_display(M)
end

* moptimize() method d2: Cannot obtain heteroskedastic-robust standard errors
mata:
    moptimize_result_display(M, "robust")
end

* moptimize() method d2: Compute heteroskedastic-robust standard errors manually
capture generate cons == 1         // We need to add a constant here
mata:
    b = moptimize_result_coefs(M)  // Estimates from previous moptimize()
    st_view(X=., ., tokens("private chronic female income cons"))
    st_view(y=., ., tokens("docvis"))
    N = rows(X)
    Xb = X*b'
    mu = exp(Xb)
    XmuXinv = cholinv(cross(X,mu,X))
    residsq = (y-mu):*(y-mu)
    Vrobust = (N/(N-1))*XmuXinv'(cross(X,residsq,X))*XmuXinv
    st_matrix("b",b)         // Pass results from Mata to Stata
    st_matrix("V",Vrobust)   // Pass results from Mata to Stata
end
matrix colnames b = private chronic female income cons
matrix colnames V = private chronic female income cons
matrix rownames V = private chronic female income cons
ereturn post b V
ereturn display

* Check
poisson docvis private chronic female income, nolog vce(robust)

**** MOPTIMIZE: METHODS GF*

// No example given as current example is cross section data

********** MOPTIMIZE: METHOD Q0

* Problems here

* moptimize() poisson q0: Poisson example
capture generate cons == 1         // We need to add a constant here
mata:
mata clear
function poissonq0(transmorphic M, real scalar todo, ///
    real rowvector b, r, S)
    {
       y = moptimize_util_depvar(M, 1)     // N x 1 vector y
           Xb = moptimize_util_xb(M, b, 1) // N x 1 vector Xb
       st_view(X=., ., tokens("private chronic female income cons"))
       mu = exp(Xb)                        // N x 1 vector mu
       r = X'(y-mu)                        // k x 1 vector r
    }
end

* moptimize() method gf1: Implement with heteroskedastic-robust standard errors
mata:
    M = moptimize_init()
    moptimize_init_evaluatortype(M, "q0")
    moptimize_init_evaluator(M, &poissonq0())
    moptimize_init_vcetype(M, "robust")
    moptimize_init_depvar(M, 1, "docvis")
    moptimize_init_eq_indepvars(M, 1, "private chronic female income")
    moptimize_init_technique(M,"gn")
    moptimize_init_conv_maxiter(M, 10)
        moptimize(M)
        b = moptimize_result_coefs(M)
    b
end


/*
W = I(5)
W
moptimize_init_gnweightmatrix(M, W)
moptimize_init_eq_coefs(M, 1, b0)
*/

* Not included
mata:
    moptimize_result_display(M)
end

poisson docvis private chronic female income, nolog

************   C.2: MATA COMAND OPTIMIZE

* mata: mata set matastrict off

mata: mata clear

* optimize() method gf2: Evaluator function 
mata 
    void pmlegf2(todo, b, y, X, lndensity, g, H)
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

// Need this to drop pmlegf2 which we redefine in next code
mata: mata clear

* optimize() method gf2: Implement with cluster–robust standard errors
mata 
    mata clear
    void pmlegf2(todo, b, y, X, lndensity, g, H)
    {
        Xb = X*b'
        mu = exp(Xb)
        lndensity = -mu + y:*Xb - lnfactorial(y)
        if (todo == 0) return
        g = (y-mu):*X   
        if (todo == 1) return
        H = - cross(X, mu, X)
    }
    st_view(y=., ., "docvis")  
    st_view(X=., ., tokens("private chronic female income cons"))
    S = optimize_init()
    optimize_init_evaluator(S, &pmlegf2())
    optimize_init_evaluatortype(S, "gf2")
    optimize_init_argument(S, 1, y)
    optimize_init_argument(S, 2, X)
    optimize_init_cluster(S, "age")  // Define cluster for cluster–robust standard errors 
    k = cols(X)
    optimize_init_params(S, J(1,k,0))
    b = optimize(S)  
    Vbrob = optimize_result_V_robust(S)
    serob = (sqrt(diagonal(Vbrob)))'
    b \ serob 
end

* Check
poisson docvis private chronic female income, vce(cluster age) nolog
  
********** END
