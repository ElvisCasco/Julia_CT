* mus203reg.do  for Stata 17

cap log close

********** OVERVIEW OF mus203reg.do **********

* Stata program
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi
* used for "Microeconometrics Using Stata, Second Edition"
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 3
* 3.2: DATA AND DATA SUMMARY
* 3.5: BASIC REGRESSION ANALYSIS
* 3.6: SPECIFICATION ANALYSIS
* 3.7: SPECIFICATION TESTS
* 3.8: SAMPLING WEIGHTS
* 3.9: OLS USING MATA

* To run you need files
*   mus203mepsmedexp.dta
* in your directory

* community-contributed commands
*   esttab 
* are used

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* File mus203mepsmedexp.dta is aothurs' extract from MEPS
* (Medical Expenditure Panel Survey)
* for individuals 65 years and older in U.S. Medicare in 2003

************ CHAPTER 3.2: DATA SUMMARY STATISTICS

* Variable description for medical expenditure dataset
use mus203mepsmedexp
describe totexp ltotexp posexp suppins phylim actlim totchr age female income

* Summary statistics for medical expenditure dataset
summarize totexp ltotexp posexp suppins phylim actlim totchr age female income

* Tabulate variable
tabulate income if income <= 0

* Detailed summary statistics of a single variable
summarize totexp, detail

* Two-way table of frequencies
table female totchr

* Two-way table with row and column percentages and Pearson chi-squared
tabulate female suppins, row col chi2

* Three-way table of frequencies
table female suppins totchr, nototals

* One-way table of summary statistics
table (result) female, stat(count totchr) stat(mean totchr) stat(sd totchr)  ///
    stat(p50 totchr)

* Two-way table of summary statistics
table female suppins, stat(count totchr) stat(mean totchr) nototals

* Summary statistics obtained using command tabstat
tabstat totexp ltotexp, statistics(count mean p50 sd skew kurt) columns(statistics)

* Kernel density plots with adjustment for highly skewed data
kdensity totexp if posexp==1, generate(kx1 kd1) n(500)
graph twoway (line kd1 kx1) if kx1 < 40000, name(levels, replace)
label variable ltotexp "Natural logarithm of expenditure"
kdensity ltotexp if posexp==1, generate(kx2 kd2) n(500)
graph twoway (line kd2 kx2) if kx2 < ln(40000), name(logs, replace)
graph combine levels logs, iscale(1.2) ysize(2.5) xsize(6.0)

*********** CHAPTER 3.5: BASIC REGRESSION ANALYSIS

* Pairwise correlations for dependent variable and regressor variables
correlate ltotexp suppins phylim actlim totchr age female income

* OLS regression with heteroskedasticity-robust standard errors
regress ltotexp suppins phylim actlim totchr age female income, vce(robust)

* Display stored results and list available postestimation commands
ereturn list
help regress postestimation

* Regression with constraints on the parameters
constraint 1 phylim = actlim
cnsreg ltotexp suppins phylim actlim totchr age female income, ///
    constraints(1) vce(robust)

* Wald test of equality of coefficients
qui regress ltotexp suppins phylim actlim totchr age female income, vce(robust)
test phylim = actlim

* Joint test of statistical significance of several variables
test phylim actlim totchr

* Store and then tabulate results from multiple regressions
qui regress ltotexp suppins phylim actlim totchr age female income, vce(robust)
estimates store REG1
qui regress ltotexp suppins phylim actlim totchr age female educyr, vce(robust)
estimates store REG2
estimates table REG1 REG2, b(%9.4f) se stats(N r2 F ll) stfmt(%9.2f) ///
    keep(suppins income educyr)

etable, estimates(REG1 REG2) keep(suppins income educyr) ///
        cstat(_r_b, nformat(%8.4f)) cstat(_r_se, nformat(%8.4f)) ///
       mstat(N) mstat(r2) mstat(F) mstat(ll) column(estimates) ///
        stars(0.1 "*" 0.05 "**" 0.01 "***") showstars showstarsnote

* Tabulate results using community-contributed command esttab to produce cleaner output
esttab REG1 REG2, b(%10.4f) se scalars(N r2 F ll) mtitles  ///
    keep(suppins income educyr) title("Model comparison of REG1-REG2")


* Tabulate results with some formatting using the table command
global xlist1 suppins phylim actlim totchr age female
table (colname[suppins income educyr] result) (command),             ///
    command(_r_b _r_se: regress ltotexp $xlist1 income, vce(robust)) ///
    command(_r_b _r_se: regress ltotexp $xlist1 educyr, vce(robust)) ///
    nformat(%8.4f) sformat("(%s)" _r_se) style(table-reg3)           ///
    stars(_r_p 0.01 "***" 0.05 "**" 0.1 "*", attach(_r_b))

* Tabulate results with even better formatting using the collect command
capture collect clear
collect title "Regression with income or education"
collect notes "Heteroskedastic-robust standard errors"
qui: collect _r_b _r_se, tag(model[(1): Income]):                 ///
    regress ltotexp suppins $xlist1 income, vce(robust)
qui: collect _r_b _r_se, tag(model[(2): Education]):              ///
    regress ltotexp suppins $xlist1 educyr, vce(robust)
collect label values colname suppins "Supplementary insurance", modify
collect stars _r_p 0.01 "***" 0.05 "**" 0.1 "*", attach(_r_b)
collect style cell, nformat(%9.4f) border(right, pattern(nil))
collect style cell result[_r_se], sformat("(%s)")
collect style header result[N r2 F ll], level(label)
collect label levels result r2 "R-squared", modify
collect style header result, level(hide)
collect style column, extraspace(1)
collect style row stack, spacer delimiter(" x ")
qui: collect layout (colname[suppins income educyr]#result result[N r2 F ll]) (model)
collect preview

* Write tabulated results to a file in Microsoft Word format
collect export mytable.docx, replace

etable, estimates(REG1 REG2) keep(suppins income educyr) ///
       cstat(_r_b, nformat(%8.4f)) cstat(_r_se, nformat(%8.4f)) ///
       mstat(N) mstat(r2) mstat(F) mstat(ll) column(estimates) ///
       stars(0.1 "*" 0.05 "**" 0.01 "***") showstars showstarsnote


* Factor variables for sets of indicator variables and interactions
regress ltotexp suppins phylim actlim totchr age female c.income ///
    i.famsze c.income#i.famsze, vce(robust) noheader allbaselevels

* Test joint significance of sets of indicator variables and interactions
testparm i.famsze c.income#i.famsze

* Factor variables for model with interactions and quadratic
regress ltotexp i.suppins i.phylim i.actlim c.totchr c.age##c.age ///
    i.female##c.income, vce(robust) noheader nofvlabel

* Average MEs in model with interactions and quadratic
margins, dydx(*) nofvlabel

* Cluster-robust standard errors
regress ltotexp suppins phylim actlim totchr age female income, vce(cluster age)

* Bootstrap to give heteroskedastic-robust standard errors
regress ltotexp suppins phylim actlim totchr age female income, ///
    vce(bootstrap, reps(999) seed(10101))

* Bootstrap to give cluster–robust standard errors
regress ltotexp suppins phylim actlim totchr age female income, ///
    vce(bootstrap, cluster(age) nodots reps(999) seed(10101))

* Bootstrap to compute confidence interval for a ratio of two coefficients
bootstrap _b[suppins]/_b[phylim], reps(999) seed(10101) nodots: ///
    regress ltotexp suppins phylim actlim totchr age female income

* Test of mean using ttest
ttest totexp = 7500

* Test of mean using regress
regress totexp, noheader
test _cons=7500

* Test of difference in means using ttest
ttest totexp, by(suppins) unequal

* Test of difference in means using regress
regress totexp suppins, vce(robust)

********** CHAPTER 3.6: SPECIFICATION ANALYSIS

* Plot of residuals against fitted values
qui regress ltotexp suppins phylim actlim totchr age female income, vce(robust)
rvfplot, msize(tiny) scale(1.2)

* Details on the outlier residuals
predict uhat, residual
predict yhat, xb
list totexp ltotexp yhat uhat if uhat < -5, clean

* Quantile-quantile plot of fitted against actual values
local endash = ustrunescape("\u2013")
qnorm uhat, msize(small) title("Q`endash'Q plot of residuals versus normal") ///
    name(graph1, replace)
kdensity uhat, normal legend(off) name(graph2, replace)
graph combine graph1 graph2, iscale(1.2) ysize(2.5) xsize(6.0)

// Not included in text - residuals from levels regression are not normal.
regress totexp suppins phylim actlim totchr age female income, vce(robust)
predict uhat2, residual
qnorm uhat2
drop uhat2

* Compute dfits that combines outliers and leverage
qui regress ltotexp suppins phylim actlim totchr age female income
predict dfits, dfits
scalar threshold = 2*sqrt((e(df_m)+1)/e(N))
display "dfits threshold = "  %6.3f threshold
tabstat dfits, statistics(min p1 p5 p95 p99 max) format(%9.3f) columns(statistics)
list dfits totexp ltotexp yhat uhat if abs(dfits) > 2*threshold & e(sample), ///
    clean

* Robust regression as a check on fat tails
qui use mus203mepsmedexp, replace
rreg ltotexp suppins phylim actlim totchr age female income, genwt(w)
estimates store ROB_DEF

* Weights used for robust regression
sum w, detail
drop w

* Median regression with default standard errors
qreg ltotexp suppins phylim actlim totchr age female income, nolog vce(robust)
estimates store MED_ROB

* Compare OLS, robust regression, and median regression
qui regress ltotexp suppins phylim actlim totchr age female income
estimates store OLS_DEF
qui regress ltotexp suppins phylim actlim totchr age female income, vce(robust)
estimates store OLS_ROB
qui qreg ltotexp suppins phylim actlim totchr age female income, nolog
estimates store MED_DEF
estimates table OLS_DEF MED_DEF ROB_DEF OLS_ROB MED_ROB, b(%10.4f) se stats(N)

* OLS, robust, and median regression for positive level of total expenditures
replace totexp = . if totexp <= 0
qui rreg totexp suppins phylim actlim totchr age female income, genwt(w)
estimates store ROB_LEVEL
qui qreg totexp suppins phylim actlim totchr age female income, vce(robust)
estimates store MED_LEVEL
qui regress totexp suppins phylim actlim totchr age female income, vce(robust)
estimates store OLS_LEVEL
estimates table OLS_LEVEL MED_LEVEL ROB_LEVEL, b(%10.4f) se stats(N)

********** CHAPTER 3.7: SPECIFICATION TESTS

* Variable augmentation test of conditional mean using estat ovtest
qui regress ltotexp suppins phylim actlim totchr age female income, vce(robust)
estat ovtest

* Link test of functional form of conditional mean
qui regress ltotexp suppins phylim actlim totchr age female income, vce(robust)
linktest

* Boxcox model with lhs variable transformed
boxcox totexp suppins phylim actlim totchr age female income if totexp>0, nolog

* Heteroskedasticity tests using estat hettest and option iid
qui regress ltotexp suppins phylim actlim totchr age female income
estat hettest, iid
estat hettest suppins phylim actlim totchr age female income, iid

* Information matrix test
qui regress ltotexp suppins phylim actlim totchr age female income
estat imtest

* Simulation to show tests have power in more than one direction
clear all
set obs 50
set seed 10101
generate x = runiform()                  // x ~ uniform(0,1)
generate u = rnormal()                   // u ~ N(0,1)
generate y = exp(1 + 0.25*x + 4*x^2) + u
generate xsq = x^2
regress y x xsq

* Test for heteroskedasticity
estat hettest

* Test for misspecified conditional mean
estat ovtest

********** CHAPTER 3.8 SAMPLING WEIGHTS

* Create artificial sampling weights
qui use mus203mepsmedexp, clear
generate swght = totchr^2 + 0.5
summarize swght

* Calculate the weighted mean
mean totexp [pweight=swght]

* Perform weighted regression
regress totexp suppins phylim actlim totchr age female income [pweight=swght]

* Weighted prediction
qui predict yhatwols
mean yhatwols [pweight=swght], noheader
mean yhatwols, noheader      // Unweighted prediction

******** CHAPTER 3.9 OLS USING MATA

* OLS with White robust standard errors using Mata
qui use mus203mepsmedexp, clear
keep if totexp > 0   // Analysis for positive medical expenditures only
generate cons = 1
local y ltotexp
local xlist suppins phylim actlim totchr age female income cons

mata
    // Create y vector and X matrix from Stata dataset
    st_view(y=., ., "`y'")             // y is nx1
    st_view(X=., ., tokens("`xlist'")) // X is nxk
    XXinv = cholinv(cross(X,X))        // XXinv is inverse of X'X
    b = XXinv*cross(X,y)               // b = [(X'X)^-1]*X'y
    e = y - X*b
    n = rows(X)
    k = cols(X)
    s2 = (e'e)/(n-k)
    vdef = s2*XXinv               // Default VCE not used here
    vwhite = XXinv*((e:*X)'(e:*X)*n/(n-k))*XXinv  // Robust VCE
    st_matrix("b",b')             // Pass results from Mata to Stata
    st_matrix("V",vwhite)         // Pass results from Mata to Stata
end

* Use Stata ereturn display to present nicely formatted results
matrix colnames b = `xlist'
matrix colnames V = `xlist'
matrix rownames V = `xlist'
ereturn post b V
ereturn display

rm mytable.docx
********** END
