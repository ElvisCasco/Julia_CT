* musappxap1.do  Oct 2009 for Stata version 11

cap log close

********** OVERVIEW OF musappxap1.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Appendix A programs
*   A.1 STATA MATRIX COMMANDS
*   A.2 STATA PROGRAMS
*   A.3 PROGRAM DEBUGGING
*   B.1 HOW TO RUN MATA
*   B.2 MATA MATRIX COMMANDS
*   B.3 PROGRAMMING IN MATA
* To run you need file
*    mus03data.dta
* in your directory
* No Stata user-written command is needed

********** SETUP **********

set more off
version 11
clear all
mata: mata set matastrict off
set scheme s1mono  /* Graphics scheme */

**************** A.1 STATA MATRIX COMMANDS

* Define a matrix explicitly and list the matrix
matrix define A = (1,2,3 \ 4,5,6)
matrix list A

* Matrix row and column names
matrix rownames A = one two
matrix list A

* Read in data, summarize and run regression
use mus03data.dta
keep if _n <= 100
drop if ltotexp == . | totchr == .
summarize ltotexp totchr  
regress ltotexp totchr, noheader 

* Create a matrix from estimation results
matrix vbols = e(V)
matrix list vbols

// Following not in book
* Create a matrix from variables
mkmat ltotexp, matrix(y)
generate intercept = 1
mkmat totchr intercept, matrix(X)

// Following not in book
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

* Can use scalar if 1x1 matrix
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

* Stata matrix commands to compute SEs and t statistics given b and V
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
program time
  display c(current_time) c(current_date)
end

* Run the program
time

* Drop program if it already exists, write program and run
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
myols ltotexp totchr if ltotexp < . in 1/100, vce(robust)

* Example of an ado-file
capture program drop meddiff
*! version 1.1.0  05oct2009
program meddiff, rclass
  version 11
  args y x 
  tempvar diff
  quietly {
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
use mus03data.dta
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

* Stata commands issued from Mata
mata
stata("summarize ltotexp")
end


* Mata help with output not included in book
mata
help mata det
help m4 matrix
help mata
end

****************   B.2 MATA MATRIX COMMANDS

mata  

// Create a matrix
A = (1,2,3 \ 4,5,6)

// List a matrix
A

// Create a 2x2 identity matrix
I = I(2)

// Create a 1x5 unit row vector with 1 in second entry and zeros elsewhere
e = e(2,5)
e

// Create a 2x5 matrix with entry 3 
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

**************** B.3 PROGRAMMING IN MATA

mata
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
  version 11
  syntax varname
  mata: calcsum2("`varlist'")
  display r(sum)
end
varsum ltotexp
  
********** CLOSE OUTPUT
