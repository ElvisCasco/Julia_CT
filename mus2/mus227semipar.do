* mus227semipar.do   for Stata 17
capture log close

********** OVERVIEW OF mus227semipar.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 27.  SEMIPARAMETRIC REGRESSION
* 27.4: NONPARAMETRIC SINGLE REGRESSOR EXAMPLE
* 27.5: NONPARAMETRIC MULTIPLE REGRESSOR EXAMPLE
* 27.6: PARTIAL LINEAR MODEL
* 27.7: SINGLE-INDEX MODEL
* 27.8: GAM

* To run you need files
*   mus202psid92m.dta
* in your directory

* And you need community-contributed commands
*   semipar
*   sls
*   gam   (which only works in MS Windows)

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus202psid92m.txt is authors' extract from the 1992 PSID for males 30-50 years

********** 27.4: NONPARAMETRIC SINGLE REGRESSOR EXAMPLE

* Read in data and choose a 5% sample
qui use mus202psid92m
drop if earnings==0 | earnings==. | earnings > 145000
set seed 10101
keep if runiform() < 0.05
sum earnings hours age edcat ibn.edcat

* Single regressor: OLS 
regress earnings hours, vce(robust)

* Single regressor: Local linear with epan2 kernel and no standard errors
npregress kernel earnings hours, kernel(epan2) 

// To speed up program reduce reps to e.g. 40 in code below

* Single regressor: Local linear with epan2 kernel and bootstrap standard errors
npregress kernel earnings hours, kernel(epan2) ///
    vce(bootstrap, seed(10101) reps(400) nodots)

* The nonidentified observation 
qui npregress kernel earnings hours, kernel(epan2) 
list earnings hours if _unident_sample == 1, clean
list hours if hours > 3590, clean

* OLS with the nonidentified obervation dropped
regress earnings hours if _unident_sample == 0, vce(robust)

* Three ways to obtain identification for all observations 
qui npregress kernel earnings hours, kernel(epan2) meanbwidth(243, copy)
matrix list e(b)   // Wider bandwidth for same kernel
qui npregress kernel earnings hours, kernel(gaussian)
matrix list e(b)   // Different unbounded kernel
qui npregress kernel earnings hours, kernel(epan2) estimator(constant) ///
    noderivatives  // Local constant rather than local linear
matrix list e(b)

* Trimmed mean dropping 5% of observations with lowest f_hat(x)
qui npregress kernel earnings hours, kernel(epan2) meanbwidth(243, copy) 
qui predict m_earnings
kdensity hours, kernel(epan2) bwidth(243) at(hours) generate(hours_x hours_d)
qui centile hours_d, centile(5)
list hours_d hours earnings m_earnings if hours_d < r(c_1), clean
sum m_earnings if hours_d >= r(c_1) 
sum m_earnings

// To speed up program reduce reps to e.g. 40 in code below

* Compute predicted mean at hours = 1000, 2000, 3000
qui npregress kernel earnings hours, kernel(epan2) ///
   vce(bootstrap, seed(10101) reps(400))
margins, at(hours=(1000(1000)3000)) vce(bootstrap, seed(10101) reps(400) nodots)

// Not included - check that bootstrap prefix gives same results
* bootstrap, seed(10101) reps(400): npregress kernel earnings hours, kernel(epan2)

// Not included - check that can do cluster bootstrap
npregress kernel earnings hours, kernel(epan2) ///
    vce(bootstrap, cluster(age) seed(10101) reps(50))
   
* Compute predicted effect of change in hours of 1,000 hours
margins, at(hours=(1000(1000)3000)) contrast(atcontrast(ar)) ///
    vce(bootstrap, seed(10101) reps(400) nodots)

// Manual calculation
generate zkernel = (hours - 3000) / 237.4
generate earn3000 = (3/4)*(1-zkernel^2)*earnings if abs(zkernel) < 1
sum earn3000

// This takes a long time with reps = 400

* Plot predictions at many levels of hours along with 95% confidence intervals
qui npregress kernel earnings hours, kernel(epan2) 
qui margins, at(hours=(400(100)3200)) ///
    vce(bootstrap, seed(10101) reps(400) nodots) 
marginsplot, legend(off) saving(graph1.gph, replace) /// 
    addplot(scatter earnings hours if earnings < 80000, msize(vsmall))  graphregion(margin(r+2))

// This takes a long time with reps = 400

* Plot effect of change in hours of 100 hours along with 95% confidence intervals
qui margins, at(hours=(400(100)3200)) contrast(atcontrast(ar)) ///
   vce(bootstrap, seed(10101) reps(400) nodots)
marginsplot, yline(0) saving(graph2.gph, replace)

* Create the figure
graph combine graph1.gph graph2.gph, ///
    iscale(1.0) xcommon rows(1) ysize(2.5) xsize(6)

* Third-order natural spline with number of knots determined by cross-validation
npregress series earnings hours, spline criterion(cv)

* Third-order natural spline with number of knots determined by cross-validation
npregress series earnings hours, spline criterion(cv)

// Not included - compute R2
predict yhatnpseries
correlate earnings yhatnpseries
display r(rho)^2

* Compute predicted mean at hours = 1000, 2000, 3000
margins, at(hours=(1000(1000)3000))

// Not included - marginsplot after npregress series
marginsplot, yline(0) 

********** 27.5: NONPARAMETRIC MULTIPLE REGRESSOR EXAMPLE

* Multiple regressor: Local linear with epan2 kernel and bootstrap standard errors
npregress kernel earnings hours age i.edcat, kernel(epan2) ///
    vce(bootstrap, seed(10101) reps(400) nodots)
count if _unident_sample == 1

* Test joint statistical significance of education categories 
matrix list e(b)
test r2vs1.edcat r3vs1.edcat r4vs1.edcat

* Multiple regressor: OLS on same sample as npregress kernel
regress earnings hours age i.edcat if _unident_sample == 0, vce(robust)

// Not included
testparm i(2/4).edcat

* Multiple regressor: OLS on same sample as npregress kernel
npregress series earnings hours age, spline criterion(cv) asis(i.edcat)

// Not included - compute R2
predict yhatmultnpseries
correlate earnings yhatmultnpseries
display r(rho)^2

********** 27.6: PARTIAL LINEAR MODEL

* Generated data: y = 1 + 1*x * f(z) + u where f(z) = z + z^2
clear
set obs 200
set seed 10101
generate x1 = rnormal()
generate x2 = rnormal() + 0.5*x1
generate z = rnormal() + 0.5*x1
generate zsq = z^2
generate y = 1 + x1 + x2 + z + zsq + 2*rnormal()

* Robinson estimator given unknown g(z) using community-contributed command semipar 
semipar y x1 x2, nonpar(z) robust ci title("Partial linear: f(z) against z") ///
    ytitle("y-b*x and f(z)") xtitle("z") 
graph save graph1, replace

* Not included: I actually get a different R-squared using corr^2(y,yhat)
qui semipar y x1 x2, nonpar(z) generate(yhatsemipar) nograph
qui correlate y yhatsemipar
display "R-squared = " r(rho)^2  

* OLS estimation given knowledge of the DGP with g(z) = z + z^2
regress y x1 x2 z zsq, vce(robust)

// Not included: compare fitted values 
predict yhatols
sum y yhatsemipar yhatols
correlate y yhatsemipar yhatols

/*
* Output not included but results reported in text: 
  Hardle and Mammen test of whether g(z) is quadratic
* Takes time. To speed up reduce nsim()
set seed 10101
semipar y x1 x2, nonpar(z) test(1) nsim(1000) nograph
set seed 10101
semipar y x1 x2, nonpar(z) test(2) nsim(1000) nograph
set seed 10101
semipar y x1 x2, nonpar(z) test(3) nsim(1000) nograph
*/

********** 27.7: SINGLE INDEX MODEL

* Ichimura semiparametric least squares for single-index model
sls y x1 x2 z zsq, trim(1,99)

* Create plot of yhat against the index x'b
predict yhat, ey
predict Index, xb
qui correlate y yhat
display "R-squared = " r(rho)^2
twoway (scatter y Index) (line yhat Index, sort lwidth(thick)), ///
    title("Single-index: yhat against x'b")                     ///
    xtitle("Index") ytitle("y and yhat") legend(off) saving(graph2, replace)

graph combine graph1.gph graph2.gph, ///
    iscale(1.2) ycommon rows(1) ysize(2.5) xsize(6)

********** 27.8: GENERALIZED LINEAR MODEL

* The following will be commented out as it only works on a Windows computer

/*
* Generalized additive model - requires Windows version of Stata
gam y x1 x2 z,  df(3)

/* 
200 records merged.
Generalized Additive Model with family gauss, link ident.
Model df     =    10.004                           No. of obs =       200
Deviance     =   820.723                           Dispersion =   4.31968
-------------------------------------------------------------------------
           y |   df    Lin. Coef.  Std. Err.      z        Gain    P>Gain
-------------+-----------------------------------------------------------
          x1 |  3.000   .8474175   .1804832     4.695     1.998    0.3683
          x2 |  3.001   .9830467   .1459322     6.736     2.134    0.3444
           z |  3.002   1.143934   .1438065     7.955   131.768    0.0000
       _cons |      1     2.1644    .146964    14.727         .         .
-------------------------------------------------------------------------
Total gain (nonlinearity chisquare) =   135.900 (6.003 df), P = 0.0000
*/

* Plot fitted smooth (and partial residual) against x for each regressor x
qui gam y x1 x2 z, df(3)
gamplot x1, saving(graph1.gph, replace)
gamplot x2, saving(graph2.gph, replace)
gamplot z, saving(graph3.gph, replace)
graph combine graph1.gph graph2.gph graph3.gph, ///
    iscale(1.2) rows(1) ysize(2) xsize(6)


* Partial linear model y = b0+b1*x1+b2*x2+g(z) estimated using GAM
gam y x1 x2 z, df(x1:1, x2:1, z:3)
display "R2 = " 1 - e(disp)*e(tdf) / 2584.6
*/
*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph

********** END
