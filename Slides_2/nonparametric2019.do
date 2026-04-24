* nonparametric2019.do  March 2019 for Stata version 15

log using nonparametric2019.txt, text replace

********** OVERVIEW OF nonparametric2019.do **********

* STATA Program by A. Colin Cameron
* Kernel density
* Nonparametric 

* To run you need files
*   nonparametric.dta  
* in your directory

* And you need Stata addons 
*    semipar
*    sls
*    gam   WHICH REQUIRES MS WINDOWS

* Also it uses npregress which is Stata 15
* This is commented out

********** SETUP **********

clear all
set more off
version 15
set scheme s1mono  /* Graphics scheme */
set linesize 82

********** DATA DESCRIPTION **********

* The original data are from the PSID Individual Level Final Release 1993 data
* From www.isr.umich.edu/src/psid  then choose Data Center 
* 4856 observations on 9 variables for Females 30 to 50 years 
* See mma09p1np.do for further description

******* OLS WITH DOCTOR VISITS DATA

* Read in data, select, describe and summarize key variables
use nonparametric.dta, clear
describe 

* Work with age 36 and nonmissing education data
keep if age == 36
drop if educatn == .
summarize

******* KERNEL DENSITY ESTIMATE

* Histogram
histogram lnhwage
histogram lnhwage, bin(30) scale(1.1)
* graph export nonparametricfig1.wmf, replace

* Kernel density 
kdensity lnhwage
kdensity lnhwage, bw(0.21)
graph twoway (kdensity lnhwage, bw(0.21))  ///
  (kdensity lnhwage, bw(0.07) clstyle(p2)) ///
  (kdensity lnhwage, bw(0.63) clstyle(p3)), legend( label(1 "Default") ///
  label(2 "Half default") label(3 "Twice default") ) scale(1.1)
* graph export nonparametricfig2.wmf, replace

* Histogram and kernel density
histogram lnhwage, kdensity  

* Kernel density and normal density with data mean and standard deviation
kdensity lnhwage, normal

******* NONPARAMETRIC REGRESSION 

* OLS 
regress lnhwage educatn

* Kernel (local constant) regression
lpoly lnhwage educatn, ci msize(small) scale(1.1)
* graph export nonparametricfig3.wmf, replace

* Local linear regression
lpoly lnhwage educatn, degree(1) ci 

* Lowess regression
lowess lnhwage educatn 

* Kernel for different bandwidths - default, halfdefault, twicedefault
graph twoway (lpoly lnhwage educatn, bw(1.5))   ///
  (lpoly lnhwage educatn, bw(0.75) clstyle(p2)) ///
  (lpoly lnhwage educatn, bw(3.0) clstyle(p3)), scale(1.1) ///
  legend(label(1 "Default") label(2 "Half default") label(3 "Twice default")) ///
  legend(pos(11) ring(0) col(1)) 

* Compare kernel, local linear, lowess with default bandwidths
graph twoway (lpoly lnhwage educ)                ///
  (lpoly lnhwage educ, degree(1) clstyle(p2))    ///
  (lowess lnhwage educ, clstyle(p3)), scale(1.1) ///
  legend( label(1 "Kernel") label(2 "Local linear") label(3 "lowess") ) ///
  legend(pos(11) ring(0) col(1))   
* graph export nonparametricfig5.wmf, replace

* OLS 
regress lnhwage educatn
regress lnhwage educatn, vce(robust)

******** NPREGRESS COMMAND

* npregress command - local linear
npregress kernel lnhwage educatn

* npregress with bootstrap standard errors
npregress kernel lnhwage educatn, vce(bootstrap, seed(10101) reps(50))

* 50 reps chosen to speed up program - should increase from 50 reps.

* Compute and plot predictions at various values with bootstrap st. errors
margins, at(educatn = (10(1)16)) vce(bootstrap, seed(10101) reps(50))
marginsplot, legend(off) scale(1.1)   ///
  addplot(scatter lnhwage educatn if lnhwage<50000, msize(tiny))
* graph export nonparametricfig11.wmf, replace

* Partial effects of changing education
margins, at(educatn = (10(1)16)) contrast(atcontrast(ar)) ///
    vce(bootstrap, seed(10101) reps(50))
marginsplot, legend(off)
* graph export nonparametricfig13.wmf, replace

******* SEMIPARAMETRIC REGRESSION  

* OLS 
regress lnhwage educatn hours, vce(robust) noheader

* Partial linear model - Robinson differencing estimator
semipar lnhwage educatn, nonpar(hours) robust ci title("Partial linear")
* graph export nonparametricfig16.wmf, replace

* Single index model - Ichimura semiparametric least squares
sls lnhwage hours educatn, trim(1,99)

* Plot predictions against the index x´b
predict yhat, ey
predict Index, xb
twoway (scatter y Index) (line yhat Index, sort lwidth(thick)), ///
   title("Single-index: yhat against x´b") scale(1.1) ///
   xtitle("Index") ytitle("y and yhat") legend(off)
* graph export nonparametricfig18.wmf, replace

* Generalized additive model - GAM REQUIRES MS WINDOWS
gam lnhwage educatn hours, df(3)

* Graphs
gamplot educatn, saving(graph1, replace)
gamplot hours, saving(graph2, replace)
graph combine graph1.gph graph2.gph, iscale(1.2) rows(1) ysize(2.5) xsize(5)
* graph export nonparametricfig21.wmf, replace

********** CLOSE OUTPUT
* log close
* clear
* exit


