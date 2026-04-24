* mus07p1qr.do  Oct 2009 for Stata version 11

cap log close

* To speed up program reduce number in reps() for bsqreg and sqreg
* the program usually uses 400
* and reduce number in rep() for qcount

********** OVERVIEW OF mus07p1qr.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Chapter 7
* 7.3: QUANTILE REGRESSION FOR MEDICAL EXPENDITURES DATA
* 7.4: QUANTILE REGRESSION FOR GENERATED HETEROSKEDASTIC DATA
* 7.5: QUANTILE REGRESSION FOR COUNT DATA

* To run you need files
*   mus03data.dta 
*   mus07qrcnt.dta   
* in your directory
* Stata user-written commands
*    qplot
*    grqreg  
*    qcount
* are used and need to be user-installed in Stata

********** SETUP **********

set more off
version 11
clear all
set scheme s1mono  /* Graphics scheme */

********** 7.3: QUANTILE REGRESSION FOR MEDICAL EXPENDITURES DATA

* Read in log of medical expenditures data and summarize
use mus03data.dta, clear
drop if ltotexp == . 
summarize ltotexp suppins totchr age female white, separator(0)

* Quantile plot for ltotexp using user-written command qplot
qplot ltotexp, recast(line) scale(1.5)
quietly graph export mus07fig1_qltotexp.eps, replace 


* Basic quantile regression for q = 0.5
qreg ltotexp suppins totchr age female white

* Obtain multiplier to convert QR coeffs in logs to AME in levels. 
quietly predict xb
generate expxb = exp(xb)
quietly summarize expxb
display "Multiplier of QR in logs coeffs to get AME in levels = " r(mean)

* Compare (1) OLS; (2-4) coeffs across quantiles; (5) bootstrap SEs
quietly regress ltotexp suppins totchr age female white   
estimates store OLS
quietly qreg ltotexp suppins totchr age female white, quantile(.25)  
estimates store QR_25
quietly qreg ltotexp suppins totchr age female white, quantile(.50) 
estimates store QR_50
quietly qreg ltotexp suppins totchr age female white, quantile(.75) 
estimates store QR_75
set seed 10101 
quietly bsqreg ltotexp suppins totchr age female white, quant(.50) reps(400) 
estimates store BSQR_50 
estimates table OLS QR_25 QR_50 QR_75 BSQR_50, b(%7.3f) se  

* Test for heteroskedasticity in linear model using estat hettest
quietly regress ltotexp suppins totchr age female white
estat hettest suppins totchr age female white, iid

* Simultaneous QR regression with several values of q
set seed 10101
sqreg ltotexp suppins totchr age female white, q(.25 .50 .75) reps(400)

* Test of coefficient equality across QR with different q
test [q25=q50=q75]: suppins

* Plots of each regressor's coefficients as quantile q varies
quietly bsqreg ltotexp suppins totchr age female white, quantile(.50) reps(400)
label variable suppins "=1 if supp ins"
label variable totchr "# of chronic condns"
grqreg, cons ci ols olsci scale(1.1)
quietly graph export mus07fig3_qrcoeff.eps, replace


********** 7.4: QUANTILE REGRESSION FOR GENERATED HETEROSKEDASTIC DATA

* Generated dataset with heteroskedastic errors
set seed 10101
set obs 10000
generate x2 = rchi2(1)
generate x3 = 5*rnormal(0)
generate e = 5*rnormal(0)
generate u = (.1+0.5*x2)*e
generate y = 1 + 1*x2 + 1*x3 + u
summarize e x2 x3 u y 

* Quantile regression for q = .25, .50 and .75
sqreg y x2 x3, quantile(.25 .50 .75)

* Predicted coefficient of x2 for q = .25, .50 and .75
quietly summarize e, detail
display "Predicted coefficient of x2 for q = .25, .50, and .75" ///
 _newline 1+.5*r(p25) _newline 1+.5*r(p50) _newline 1+.5*r(p75)

* Generate scatterplots and qplot
quietly kdensity u, scale(1.25) lwidth(medthick) saving(density, replace)
quietly qplot y, recast(line) scale(1.4) lwidth(medthick) saving(quanty, replace)
quietly scatter y x2, scale(1.25) saving(yversusx2, replace)
quietly scatter y x3, scale(1.25) saving(yversusx3, replace)
graph combine density.gph quanty.gph yversusx2.gph yversusx3.gph
quietly graph export mus07hetexample.eps, replace

* Quantile regression for q = .25, .50 and .75
sqreg y x2 x3, quantile(.25 .50 .75)

* Predicted coefficient of x2 for q = .25, .50 and .75
quietly summarize e, detail
display "Predicted coefficient of x2 for q = .25, .50, and .75" ///
 _newline 1+.5*r(p25) _newline 1+.5*r(p50) _newline 1+.5*r(p75)



* OLS and quantile regression for q = .25, .5, .75
quietly regress y x2 x3
estimates store OLS
quietly regress y x2 x3, vce(robust)
estimates store OLS_Rob
quietly bsqreg y x2 x3, quantile(.25) reps(400)
estimates store QR_25
quietly bsqreg y x2 x3, quantile(.50) reps(400)
estimates store QR_50
quietly bsqreg y x2 x3, quantile(.75) reps(400)
estimates store QR_75
estimates table OLS OLS_Rob QR_25 QR_50 QR_75, b(%7.3f) se  

* Test equality of coeff of x2 for q=.25 and q=.75
set seed 10101
quietly sqreg y x2 x3, q(.25 .75) reps(400)
test [q25]x2 = [q75]x2
test [q25]x3 = [q75]x3

********** 7.5: QUANTILE REGRESSION FOR COUNT DATA

* Read in doctor visits count data and summarize
use mus07qrcnt.dta, clear
summarize docvis private totchr age female white, separator(0)

* Generate jittered values and compare quantile plots
set seed 10101
generate docvisu = docvis + runiform()
quietly qplot docvis if docvis < 40, recast(line) scale(1.25) ///
  lwidth(medthick) saving(docvisqplot, replace)
quietly qplot docvisu if docvis < 40, recast(line) scale(1.25) ///
  lwidth(medthick) saving(docvisuqplot, replace)
graph combine docvisqplot.gph docvisuqplot.gph
graph export mus07dvqplot2.eps, replace

* Marginal effects from conventional negative binomial model
quietly nbreg docvis i.private totchr age i.female i.white, vce(robust)
margins, dydx(*) atmean noatlegend

* Quantile count regression
set seed 10101 
qcount docvis private totchr age female white, q(0.50) rep(500)

* Marginal effects after quantile regression for median
set linesize 81
qcount_mfx

* Quantile count regression for q = 0.75
set seed 10101
quietly qcount docvis private totchr age female white, q(0.75) rep(500)
qcount_mfx

* Save coefficients and generate graph for a range of quantiles
use mus03data.dta, clear
drop if ltotexp == . 
capture program drop mus07plot
program mus07plot
 postfile myfile percentile b1 upper lower using bsqrcoef1.dta, replace 
 forvalues tau1=0.10(0.1)0.9 {
 set seed 10101
 quietly bsqreg ltotexp suppins age female white totchr, quant(`tau1') reps(400)
 matrix b = e(b)
 scalar b1=b[1,1]
 matrix V = e(V)
 scalar v1=V[1,1]
 scalar df=e(df_r) 
 scalar upper = b1 + invttail(df,.025)*sqrt(v1)
 scalar lower = b1 - invttail(df,.025)*sqrt(v1)
 post myfile (`tau1') (b1) (upper) (lower) 
 matrix drop V b
 scalar drop b1 v1 upper lower df
 }
 postclose myfile
end
mus07plot
program drop mus07plot
use bsqrcoef1.dta, clear
twoway connected b1 percentile || line upper percentile||line lower percentile, ///
  title("Slope Estimates") subtitle("Coefficient of suppins") ///
  xtitle("Quantile", size(medlarge)) ///
  ytitle("Slope and confidence bands", size(medlarge)) ///
  legend(label(1 "Quantile slope coefficient") ///
  label(2 "Upper 95% bs confidence band")  ///
  label(3 "Lower 95% bs confidence band"))  
graph save bsqrcoef1.gph, replace


********** CLOSE OUTPUT **********
