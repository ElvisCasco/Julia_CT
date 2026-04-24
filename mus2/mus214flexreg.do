* mus214flexreg.do  for Stata 17

cap log close

********** OVERVIEW OF mus214flexreg.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press 

* CHAPTER 14
* 14.2 : MODELS BASED ON FINITE MIXTURES
* 14.3 : FMM EXAMPLE: EARNINGS OF DOCTORS 
* 14.4 : GLOBAL POLYNOMIALS
* 14.5 : REGRESSION SPLINES 
* 14.6 : NONPARAMETRIC REGRESSION 

* To run you need files
*   mus214mabelfmm.dta 
* in your directory

* community-contributed command
*    gam   (which will work only in Windows version of Stata)
* is used
 
************ SETUP ***********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

************ DATA DESCRIPTION ***********

* File mus214mabelfmm.dta is authors' extract from MABEL waves 1-3 
* Medicine in Australia: Balancing Employment and Life (MABEL)
* An Australian longitudinal survey of physicians.

************* 14.2: MODELS BASED ON FINITE MIXTURES 

* Compare mixture of normal densities with r.v. that is weighted sum of normals 
set obs 10000
generate y_finite_mixture = rnormal(4,1)
replace y_finite_mixture = rnormal(8,1) if runiform()>0.75
kdensity y_finite_mixture, saving(graph1.gph, replace) ///
    title("Finite mixture of normals") note(" ")
generate y_weighted_sum = 0.75*rnormal(4,1)+0.25*rnormal(8,1) 
kdensity y_weighted_sum, saving(graph2.gph, replace)   /// 
    title("Weighted sum of normals") note(" ")
sum y_finite_mixture y_weighted_sum 

graph combine graph1.gph graph2.gph, ycommon xcommon ysize(2.5) xsize(6) ///
    iscale(1.2)

* DGP for regression 3-component mixture of normals: Well-separated components
clear all 
set obs 10000
set seed 10101
generate x = runiform() 
generate y = 1 + 1*x + rnormal() 
generate class = runiform()                      
replace y = 4 + 4*x + .8*rnormal() if class > 0.5
replace y = 8 + 8*x + .4*rnormal() if class > 0.8  

* Regression mixture of normals: Summary statistics and density
summarize y x 
kdensity y, kernel(gaussian) lwidth(medthick) note(" ") scale(1.2)

********* 14.3: FMM EXAMPLE: EARNINGS OF DOCTORS 
 
* Log earnings of doctors example: Read in data and give summary statistics
qui use mus214mabelfmm, clear
global xlist logyhrs female childu5 expr exprsq pgradoth pracsize 
summarize logyearn $xlist, sep(0)
tabstat logyearn, stat(mean p50 sd skewness kurtosis)  

// Not included but mentioned in text
kdensity logyearn, normal lwidth(medthick)

* Standard OLS log-earnings regression (= 1 component FMM) with robust standard errors
regress logyearn $xlist, vce(robust) 
estimates store fmm1

// Not included but mentioned in text
predict uhat, resid
tabstat uhat, stat(mean p50 sd skewness kurtosis)  
kdensity uhat

* Two-component mixture of (log) normals: ML estimates
fmm 2, nolog vce(robust): regress logyearn $xlist 
estimates store fmm2

// Not included but mentioned in text - regression output with log file
fmm 2, vce(robust): regress logyearn $xlist 

* Two-component mixture: Predicted conditional means for each observation
predict lyearnhat*, mu
summarize lyearnhat*

* Distribution of conditional means for the two components
twoway (kdensity lyearnhat1, lwidth(medthick) clstyle(p1) note(" "))  ///
    (kdensity lyearnhat2, lwidth(medthick) clstyle(p2) note(" ")),    ///
    legend(pos(11) ring(0) col(1)) legend(size(small))                /// 
    xtitle("Conditional mean of log earnings")                        /// 
    title("Conditional means") saving(graph1.gph, replace)
	     
* Two-component mixture: estat lcmean gives average predicted conditional means
estat lcmean

* Random draws from the two components
generate randomly1 = lyearnhat1 + rnormal(0,sqrt(.4314492)) 
generate randomly2 = lyearnhat2 + rnormal(0,sqrt(.1555543)) 
summarize randomly1 randomly2
twoway (kdensity randomly1, lwidth(medthick) clstyle(p1) note(" ")) ///
    (kdensity randomly2, lwidth(medthick) clstyle(p2) note(" "))    ///
    (kdensity logyearn, lwidth(medthick) clstyle(p3) note(" ")),    /// 
    legend(pos(11) ring(0) col(1)) legend(size(small))              ///
    xtitle("Random draw of log earnings")                           /// 
    title("Random draws") saving(graph2.gph, replace)
	  
graph combine graph1.gph graph2.gph, scale(1.4) ysize(2.5) xsize(6)

* Two-component mixture: Average marginal effects
margins, dydx(*)

* Two-component mixture: Predicted class probabilities for each observation
predict classprob*, classpr
summarize classprob*

* Two-component mixture: estat lcprob gives average predicted class probs
estat lcprob

* Two-component mixture: Predicted posterior class probabilities for each obs
predict postprob*, classposteriorpr
summarize postprob*

* Estimate one- two- three- and four-component mixture of normals
qui fmm 1, vce(robust): regress logyearn $xlist
estimates store fmm1   
qui fmm 2, vce(robust): regress logyearn $xlist
estimates store fmm2
qui fmm 3, vce(robust): regress logyearn $xlist
estimates store fmm3
estat lcprob
qui fmm 4, vce(robust): regress logyearn $xlist
estimates store fmm4   
estat lcprob 

* Model comparison: one- two- three- and four-component mixture of normals
estimates stats fmm1 fmm2 fmm3 fmm4

* Coeflegend gives full names of regression coefficients 
qui fmm 2, vce(robust): regress logyearn logyhrs female childu5 ///
    expr exprsq pgradoth pracsize 
fmm, coeflegend

* FMM2: Test of a single restriction across latent classes
test (_b[logyearn:1.Class#c.logyhrs] = _b[logyearn:2.Class#c.logyhrs])

* FMM2: Same test using the contrast command
contrast c.logyhrs#a.Class, equation(logyearn) 

* FMM2: Test of a joint restriction across latent classes
test (_b[logyearn:1.Class#c.female] = _b[logyearn:2.Class#c.female]) ///
    (_b[logyearn:1.Class#c.childu5]=  _b[logyearn:2.Class#c.childu5])

* FMM2 model with varying mixture probabilities
fmm 2, lcprob(female) lcbase(2) nolog vce(robust):  ///
    regress logyearn logyhrs expr exprsq pgradoth pracsize  

********** 14.4 GLOBAL POLYNOMIALS 

* Generated data: y = 1 + 1*x1 + 1*x2 + f(z) + u where f(z) = z + z^2
clear
set obs 200
set seed 10101
generate x1 = rnormal()
generate x2 = rnormal() + 0.5*x1
generate z = rnormal() + 0.5*x1
generate zsq = z^2
generate y = 1 + x1 + x2 + z + zsq + 2*rnormal()
summarize

// Not included - estimate same model as DGP
reg y x1 x2 z zsq

* Quartic global polynomial model 
reg y c.z##c.z##c.z##c.z, vce(robust)

* Graph comparing quartic model predictions to quadratic model predictions
predict yquartic, xb
sort z
twoway (scatter y z, msize(small)) (qfit y z, lwidth(medthick)clstyle(p2)) ///
    (line yquartic z, lwidth(medthick)), scale(1.2)                        ///
    legend(pos(11) ring(0) col(1)) legend(size(small))                     ///
    legend(label(1 "Actual data") label(2 "Quadratic") label(3 "Quartic")) 


// Not included - use npregress series comand instead
npregress series y z, polynomial(4)
predict yquarticnpseries
correlate yquartic yquarticnpseries

* Quartic model estimated using fractional polynomial command
fp <z>, fp(1 2 3 4) scale replace: regress y <z>
predict yfpquartic
correlate yquartic yfpquartic

* Orthogonalize the quartic polynomials 
orthpoly z, generate(pz*) deg(4)
correlate pz*

// Not included 
regress y pz*, vce(robust)

********** 14.5 REGRESSION SPLINES

* Create the basis function manually with three segments and knots at -1 and 1
generate zseg1 = z
generate zseg2 = 0
replace zseg2 = z - (-1) if z > -1
generate zseg3 = 0
replace zseg3 = z - 1 if z > 1

* Piecewise linear regression with three sections
regress y zseg1 zseg2 zseg3, vce(robust)
predict yhat
twoway (scatter y z) (line yhat z, sort lwidth(thick)),                    ///
    title("Piecewise linear: y=a+f(z)+u") ytitle("y and f(z)") xtitle("z") ///
    legend(off) saving(graph1.gph, replace)
 
* Repeat piecewise linear using command mkspline to create the basis functions
mkspline zmk1 -1 zmk2 1 zmk3 = z, marginal
summarize zseg1 zmk1 zseg2 zmk2 zseg3 zmk3, sep(8) 
regress y zmk1 zmk2 zmk3, vce(robust) noheader

// Not included - use npregress series comand instead
npregress series y z, knots(4) spline
matrix define mknots = (-1, 1)
matrix list mknots
npregress series y z, knotsmat(mknots) spline(1)
predict yhatspline
correlate yhat*

* Piecewise regression with additional regressors x1 and x2
regress y x1 x2 zmk1 zmk2 zmk3, vce(robust) 
generate partresid = y - _b[_cons] - _b[x1]*x1 - _b[x2]*x2
generate fz = _b[zmk1]*zmk1 + _b[zmk2]*zmk2 + _b[zmk3]*zmk3 
twoway (scatter partresid z) (line fz z, sort lwidth(thick)), xtitle("z") ///
    title("Piecewise linear: y=a+b*x1+c*x2+f(z)+u") legend(off)       ///
    ytitle("Partial residual and f(z)")  saving(graph2.gph, replace)

graph combine graph1.gph graph2.gph, xcommon ///
    iscale(1.2) rows(1) ysize(2.5) xsize(6)
  
* Natural or restricted cubic spline regression of y on z
mkspline zspline = z, cubic nknots(5) displayknots
regress y zspline*, vce(robust)

* Plot the predicted values from natural cubic spline regression
predict yhatnatural
twoway (scatter y z) (line yhatnatural z, sort lwidth(thick)), ///
    title("Natural cubic spline: y=a+f(z)+u") xtitle("z")      ///
    ytitle("f(z)") legend(off) saving(graph1.gph, replace)

/*
// Following with gam will work only in Windows version of Stata

* Smoothing spline regression of y on z - requires Windows version of Stata
gam y z, df(4)
display "R2 = " 1 - e(disp)*e(tdf) / 2584.6
twoway (scatter y z) (line GAM_mu z, sort lwidth(thick)), ///
    title("Smoothing spline: y=a+f(z)+u") xtitle("z")     ///
    ytitle("f(z)") legend(off) saving(graph2.gph, replace)

graph combine graph1.gph graph2.gph, ///
    iscale(1.2) rows(1) ysize(2.5) xsize(6)

*/
************* 14.6 NONPARAMETRIC REGRESSION 

* Local linear using lpoly
lpoly y z, degree(1) at(z) generate(yhatlpoly) kernel(epan2)     ///
    title("Local linear using lpoly") saving(graph1.gph, replace) 

* Local linear using npregress kernel
npregress kernel y z, estimator(linear) kernel(epan2) nolog  ///
    vce(bootstrap, seed(10101) reps(400) nodots)  
npgraph, title("Local linear using npregress") saving(graph2.gph, replace)  

graph combine graph1.gph graph2.gph, ycommon xcommon iscale(1.2) ///
    ysize(2.5) xsize(6)

* Compare predicted conditional means from lpoly, npregress, and lowess
predict yhatnp, mean              // Prediction from npregress
lowess y z, generate(yhatlowess)
sum y yhatnp yhatlpoly yhatlowess
correlate y yhatnp yhatlpoly yhatlowess

// Not included - npregress series comand
* Examples of the npregress series command with specified degree and knot
npregress series y z, polynomial(4)
npregress series y z, spline knots(2)
* Examples of npregress series command with data-determined degree and knot
npregress series y z, polynomial
npregress series y z, spline

*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph

**** END  

                                  
