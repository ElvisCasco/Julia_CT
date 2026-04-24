* mus215quantile.do for Stata 17

cap log close

********** OVERVIEW OF mus215quantile.do **********

* Stata program
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi
* used for "Microeconometrics Using Stata, Second Edition"
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* CHAPTER 15
* 15.2: CQR
* 15.3: CQR FOR MEDICAL EXPENDITURES DATA
* 15.4: CQR FOR GENERATED HETEROSKEDASTIC DATA
* 15.5: QUANTILE TREATMENT EFFECTS FOR A BINARY TREATMENT

* To run you need files
*   mus203mepsmedexp.dta
* in your directory

* community-contributed commands
*    qplot
*    grqreg
*    qreg2
* are used and need to be installed in Stata
* To speed up program reduce number in reps() for bsqreg and sqreg
* the program usually uses 400

********** SETUP **********

set linesize 82
clear all
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* File mus203mepsmedexp is authors' extract from MEPS
* (Medical Expenditure Panel Survey)
* for individuals 65 years and older in U.S. Medicare in 2003

********** CHAPTER 15.3: QUANTILE REGRESSION FOR MEDICAL EXPENDITURES DATA

* Read in log of medical expenditures data and summarize
qui use mus203mepsmedexp
drop if ltotexp == .
summarize ltotexp suppins totchr age female white, separator(0)

* Intercept-only quantile regression gives the raw quantile
centile ltotexp, centile(25)
qreg ltotexp, quantile(.25) nolog

* Quantile plot for ltotexp using the community-contributed command qplot
qplot ltotexp, recast(line) scale(1.5) ytitle("Quantiles of ln(totexp)") xtitle("Fraction of the data")

* Basic quantile regression for q = 0.5
qreg ltotexp suppins totchr age female white

* Obtain multiplier to convert QR coeffs in logs to AME in levels.
qui predict xb
generate expxb = exp(xb)
qui summarize expxb
display "Multiplier of QR in logs coeffs to get AME in levels = " r(mean)

* Compare (1) OLS; (2-4) coeffs across quantiles .25, .50, and .75
qui regress ltotexp suppins totchr age female white, vce(robust)
estimates store OLS
qui qreg ltotexp suppins totchr age female white, quantile(.25) vce(robust)
estimates store QR_25
qui qreg ltotexp suppins totchr age female white, quantile(.50) vce(robust)
estimates store QR_50
qui qreg ltotexp suppins totchr age female white, quantile(.75) vce(robust)
estimates store QR_75
estimates table OLS QR_25 QR_50 QR_75, b(%7.3f) se

* qreg2 provides heteroskedastic-robust standard errors and also tests for heteroskedasticity
qreg2 ltotexp suppins totchr age female white, quantile(0.5)

* Median regression: i.i.d., heteroskedastic-robust, and cluster–robust standard errors
qui qreg ltotexp suppins totchr age female white, quantile(.50)
estimates store IID
qui qreg ltotexp suppins totchr age female white, quantile(.50) vce(robust)
estimates store HETWHITE
qui qreg2 ltotexp suppins totchr age female white, quantile(.50)
estimates store HETANAL
set seed 10101
qui bsqreg ltotexp suppins totchr age female white, quant(.50) reps(400)
estimates store HETBOOT
qui qreg2 ltotexp suppins totchr age female white, quant(.50) cluster(educyr)
estimates store CLUANAL
qui bootstrap, cluster(educyr) seed(10101) reps(400):       ///
   qreg ltotexp suppins totchr age female white, quant(.50)
estimates store CLUBOOT
estimates table IID HETWHITE HETANAL HETBOOT CLUANAL CLUBOOT, b(%7.3f) se

// Not included - intracluster correlation tests
qreg2 totexp suppins totchr age female white, quantile(0.50) cluster(educyr)

* Test for heteroskedasticity in linear model using estat hettest
qui regress ltotexp suppins totchr age female white
estat hettest suppins totchr age female white, iid

* Simultaneous QR regression with several values of q
set seed 10101
sqreg ltotexp suppins totchr age female white, q(.25 .50 .75) reps(400) nodots

* Test of coefficient equality across QR with different q
test [q25=q50=q75]: suppins

* Plots of each regressor's coefficients as quantile q varies
set seed 10101
qui bsqreg ltotexp suppins totchr age female white, reps(400)
grqreg, cons ci ols olsci scale(1.1) seed(10101)

********** CHAPTER 15.4: QUANTILE REGRESSION FOR GENERATED HETEROSKEDASTIC DATA

* Generated dataset with heteroskedastic errors
clear all
set seed 10101
qui set obs 10000
generate x2 = rchi2(1)
generate x3 = 5*rnormal(0)
generate e = 5*rnormal(0)
generate u = (.1+0.5*x2)*e
generate y = 1 + 1*x2 + 1*x3 + u
summarize e x2 x3 u y

* Generate scatterplots and qplot
qui kdensity u, scale(1.25) lwidth(medthick) saving(graph1.gph, replace)
qui qplot y, recast(line) scale(1.4) lwidth(medthick) saving(graph2.gph, replace) xtitle("Fraction of the data")
qui scatter y x2, msize(tiny) scale(1.25) saving(graph3.gph, replace)
qui scatter y x3, msize(tiny) scale(1.25) saving(graph4.gph, replace)
graph combine graph1.gph graph2.gph graph3.gph graph4.gph

* OLS and quantile regression for q = .25, .5, .75
qui regress y x2 x3
estimates store OLS
qui regress y x2 x3, vce(robust)
estimates store OLS_Rob
set seed 10101
qui bsqreg y x2 x3, quantile(.25) reps(400)
estimates store QR_25
qui bsqreg y x2 x3, quantile(.50) reps(400)
estimates store QR_50
qui bsqreg y x2 x3, quantile(.75) reps(400)
estimates store QR_75
estimates table OLS OLS_Rob QR_25 QR_50 QR_75, b(%7.3f) se

* Predicted coefficients of x2 from quantile regression
qui summarize e, detail
display "Predicted coefficient of x2 for q = .25, .50, and .75" _newline ///
    "     are       " 1+.5*5*invnormal(0.25) ", " 1+.5*5*invnormal(0.5)  ///
    ", and " 1+.5*5*invnormal(0.75)

* Test equality of coeff of x2 and equality of coeff of x3 for q=.25 and q=.75
set seed 10101
qui sqreg y x2 x3, q(.25 .75) reps(400)
test [q25]x2 = [q75]x2
test [q25]x3 = [q75]x3

********** CHAPTER 15.5: QUANTILE TREATMENT EFFECT FOR BINARY TREATMENT

* QTE at q=0.25 by separate estimation for suppins==0 and suppins==1
qui use mus203mepsmedexp, clear
qui drop if ltotexp == .
qreg ltotexp if suppins==0, quantile(0.25) vce(robust) nolog
scalar q25_0 = _b[_cons]
qreg ltotexp if suppins==1, quantile(0.25) vce(robust) nolog
scalar q25_1 = _b[_cons]
display "QTE of suppins at q=25 = " q25_1 - q25_0

* QTE at q=0.25 by direct CQR on intercept and suppins
qreg ltotexp suppins, quantile(0.25) vce(robust) nolog

* Quantile plot of ltotexp for those with and without supplementary insurance
qplot ltotexp, over(suppins) clp(1 _) recast(line) scale(1.1)          ///
    ytitle("Quantiles of lntotexp") xtitle("Fraction of the data")    ///
    legend(pos(11) ring(0) col(1))       ///
    legend(label(1 "suppins=0") label(2 "suppins==1")) saving(graph1, replace)

* Plot of the QTE of supplementary insurance at different quantiles
set seed 10101
qui bsqreg ltotexp suppins, reps(400)
grqreg, ci ols olsci scale(1.1) seed(10101)
graph save graph2, replace

graph combine graph1.gph graph2.gph, ysize(2.5) xsize(6.0) iscale(1.2)

*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph
erase graph3.gph
erase graph4.gph

********** END
