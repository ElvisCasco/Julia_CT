* mus201stata.do  for Stata 17

capture log close

********** OVERVIEW OF mus201stata.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 1
*** 1.2 DOCUMENTATION
*** 1.3 COMMAND SYNTAX AND OPERATORS
*** 1.4 DO-FILE AND LOG FILE
*** 1.5 SCALARS AND MATRICES
*** 1.6 USING RESULTS FROM STATA COMMANDS
*** 1.7 GLOBAL AND LOCAL MACROS
*** 1.8 LOOPING COMMANDS
*** 1.11 TEMPLATE DO-FILE

* To run you need file access to file auto.dta 
* which is installed with Stata

* No community-contributed commands are used

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* File auto.dta is a dataset installed with Stata

********** 1.2 DOCUMENTATION

help contents

********** 1.3 COMMAND SYNTAX AND OPERATORS

* Use a dataset that was installed with Stata
sysuse auto
summarize

* Summary statistics of selected variables
summarize mpg price weight, separator(1)

* OLS regression
regress mpg price weight

* by foreign: regress mpg price weight if weight < 4000, vce(robust)

* Suppress output from command
quietly regress mpg price weight

* Factor variables for rep78 
summarize i.rep78

* Factor variables for rep78 in regression
regress mpg i.rep78

* Factor variables for rep78 - no category is omitted
regress mpg ibn.rep78, hascons

* Factor variables for interaction between two categorical variables
regress mpg i.rep78#i.foreign, allbaselevels

* Factor variables for interaction between categorical and continuous variables
regress mpg i.rep78#c.weight

* Factor variables for interaction between two continuous variables
regress mpg price c.weight c.weight#c.weight, noheader

* Wildcard (asterisk) example
summarize t*

* Arithmetic example and display of result
display -2*(9/(8+2-7))^2

* Missing value created by impossible arithmetic operation
display 2/0

* Example of logical operators - output omitted in book
* Commented out as autodata.dta needs to be re-read in first
* regress mpg price weight if weight <= 4000 & price <- 10000

* Error message example - following command will fail so commented out
* regress mpg notthere

********** 1.4 DO-FILE AND LOG FILE

/* Following commands and output are instead directly in chapter 1
* Change working directory 
cd "c:\Program Files\Stata16"
* List the contents of a file
type example.do
* Run the do-file
do example.do
*/

* Log file example saving output as monospaced (fixed-width) text 
*log using logs/example.txt, text replace

/* Following commands are instead directly in chapter 1
* Demonstrate use of comments
* This program reads in system file auto.dta and gets summary statistics
clear           // Remove data from memory
* The next two lines show how to split a single command across two lines
sysuse    ///
  auto.dta
summarize
*/

/* Following commands are instead directly in chapter 1
* Change delimiter from cr to semicolon and back to cr
#delimit ;
* More than one command per line and command spans more than one line;
clear; sysuse 
    auto.dta; summarize;
#delimit cr
*/


********** 1.5 SCALARS AND MATRICES

* Scalars: Example
scalar a = 2*3
scalar b = "2 times 3 = "
display b a

* Matrix commands: Example
matrix define A = (1,2,3 \ 4,5,6)
matrix list A
scalar c = A[2,3]
display c

********** 1.6 USING RESULTS

* Illustrate use of return list for r-class command summarize
summarize mpg
return list

* Illustrate use of r()
qui summarize mpg
scalar range = r(max) - r(min)
display "Sample range = " range

* Save a result in r() as a scalar 
scalar mpgmean = r(mean)

regress mpg price weight

* ereturn list after e-class command regress
ereturn list

* Use of e() where scalar
scalar r2 = e(mss)/(e(mss)+e(rss))
display "r-squared = " r2

* Use of e() where matrix
matrix best = e(b)
scalar bprice = best[1,1]
matrix Vest = e(V)
scalar Vprice = Vest[1,1]
scalar tprice = bprice/sqrt(Vprice)
display "t statistic for H0: b_price = 0 is " tprice 

* Direct use of matrix subscripts for e() matrices
display e(b)[1,1]/sqrt(e(V)[1,1])

********** 1.7 MACROS: GLOBAL AND LOCAL

* Global macro definition and use
global xlist price weight
regress mpg $xlist, noheader     // $ prefix for global macro is necessary

* Local macro definition and use
local xlist "price weight"
regress mpg `xlist', noheader     // Single quotes are necessary

* Local macro definition without double quotes
local y mpg           
regress `y' `xlist', noheader

* Local macro definition through function evaluation
local z = 2+2
display `z'

* Scalars disappear after clear all but macro does not
global b 3
local c 4
scalar d = 5
clear
display $b _skip(3) `c'   // Display macros
display d                 // Display the scalar
clear all  
display $b _skip(3) `c'   // Display macros
* Next line commented out as it will lead to error message
* display d               // Display the scalar

********** 1.8 LOOPS

* Make artificial dataset of 100 observations on 4 uniform variables
clear
set obs 100
set seed 10101
generate x1var = runiform()
generate x2var = runiform()
generate x3var = runiform()
generate x4var = runiform()

* Manually obtain the sum of four variables
generate sum = x1var + x2var + x3var + x4var
summarize sum

* foreach loop with a variable list
qui replace sum = 0
foreach var of varlist x1var x2var x3var x4var {
    qui replace sum = sum + `var'
}
summarize sum

* forvalues loop to create a sum of variables
qui replace sum = 0
forvalues i = 1/4  {
    qui replace sum = sum + x`i'var
}
summarize sum

* While loop and local macros to create a sum of variables
qui replace sum = 0
local i 1
while `i' <= 4 {
    qui replace sum = sum + x`i'var
    local i = `i' + 1
}
summarize sum

********** 1.10 TEMPLATE DO-FILE

/*
* 1. Program name
* mus01p2template.do written 12/01/2019 is a template do-file
* 2. Write output to a log file
log using logs/mus01p2template.txt, text replace
* 3. Stata version
version 17           // So will still run in a later version of Stata
* 4. Program explanation
* This illustrative program creates 100 uniform variates
* 5. Change Stata default settings - one example is given
set linesize 82     // Set the maximum width of Stata output
* 6. Set program parameters using logs/local and global macros
global numobs 100
local seed 10101
local xlist xvar
* 7. Generate data and summarize
set obs $numobs
set seed `seed'
generate xvar = runiform()
generate yvar = xvar^2
summarize
* 8. Demonstrate use of results stored in r()
summarize xvar
display "Sample range = " r(max)-r(min)
regress yvar `xlist', vce(robust)
scalar r2 = e(mss)/(e(mss)+e(rss))
display "r-squared = " r2
* 9. Close output file and exit Stata
log close
* exit, clear
*/

*** 1.11 STATA ADD-ON COMMANDS

search instrumental variables

********** END
