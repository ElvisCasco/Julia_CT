* ML_2022_part4.do  April 2022 for Stata version 16

capture log close

log using ML_2022_part4.txt, text replace

********** OVERVIEW OF ML_2022_part4.do **********

* To run you need files
*  mus203mepsmedexp.dta
* in your directory

* And Stata user-written commands 
*   rforest
* are used

* 2.1 DINMENSION REDUCTION: PRINCIPAL COMPONENTS
* 3. BASIS FUNCTIONS - GLOBAL POLYNOMIALS, SPLINES
* 4. NEURAL NETWORK EXAMPLE
* 6. PREDICTION EXAMPLE

********** SETUP **********

set more off
version 16
clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Data for Principla Components are ggenerated

* Data for Prediction example
* File mus203mepsmedexp.dta is aothurs' extract from MEPS 
* (Medical Expenditure Panel Survey)
* for individuals 65 years and older in U.S. Medicare in 2003

********** 2.1 DINMENSION REDUCTION: PRINCIPAL COMPONENTS

* Generate three correlated variables (rho = 0.5) and y linear only in x1
clear
quietly set obs 40
set seed 12345
matrix MU = (0,0,0)
scalar rho = 0.5
matrix SIGMA = (1,rho,rho \ rho,1,rho \ rho,rho,1)
drawnorm x1 x2 x3, means(MU) cov(SIGMA)
generate y = 2 + 1*x1 + rnormal(0,3)
saveold ML_2022_part4, version(11) replace

* Standardize regressors and demean y
foreach var of varlist x1 x2 x3 {
     qui egen double z`var' = std(`var')
}
qui summarize y
qui generate double ydemeaned = y - r(mean)
summarize ydemeaned z*

* Summarize data
summarize
correlate

* OLS regression of y on x1-x3
regress y x1 x2 x3, vce(robust)

capture drop pc* yhat
* Principal components using default option that first standardizes the data
pca x1 x2 x3

* Not included - get same results using standardized data and covariance option
pca zx1 zx2 zx3, covariance

* Compute the 3 principal components and their means, st.devs., correlations
predict pc1 pc2 pc3
summarize pc1 pc2 pc3
correlate pc1 pc2 pc3

* Manually compute the first principal component and compare to pc1
generate double pc1manual = 0.6306*zx1 +  0.5712*zx2 + 0.5254*zx3
summarize pc1 pc1manual

capture drop yhat

* Compare R from OLS on all three regressors, on pc1, on x1, on x2, on x3
qui regress y x1 x2 x3
predict yhat
correlate y yhat pc1 x1 x2 x3

// Not included
* Compare OLS on x1 with OLS on first principal component
regress y x1
regress y pc1

********** 3. BASIS FUNCTIONS - GLOBAL POLYNOMIALS, SPLINES

*** GLOBAL POLYNOMIALS

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
    legend(label(1 "Actual data") label(2 "Quadratic") label(3 "quartic")) 

// Not included - use npregress series command instead
npregress series y z, polynomial(4)
predict yquarticnpseries
correlate yquartic yquarticnpseries

*** REGRESSION SPLINES

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
summarize zseg1 zmk1 zseg2 zmk2 zseg3 zmk3, sep (8) 
regress y zmk1 zmk2 zmk3, vce(robust) noheader

// Not included - use npregress series commmand instead
npregress series y z, knots(4) spline
matrix define mknots = (-1, 1)
matrix list mknots
npregress series y z, knotsmat(mknots) spline(1)
predict yhatspline
correlate yhat*

* Natural or restricted cubic spline regression of y on z
mkspline zspline = z, cubic nknots(5) displayknots
regress y zspline*, vce(robust)

* Plot the predicted values from natural cubic spline regression
predict yhatnatural
twoway (scatter y z) (line yhatnatural z, sort lwidth(thick)), ///
    title("Natural cubic spline: y=a+f(z)+u") xtitle("z")      ///
    ytitle("f(z)") legend(off) saving(graph1.gph, replace)

********** 4. NEURAL NETWORK EXAMPLE

/* Following gave me error message - so comment out
. brain define, input(x) output(y) hidden(20)
non-natively compiled windows plugin detected, e.g. cygwin/mingw
unable to load brainwin.plugin from directory C:\Users\ccameron\ado\plus/b/
perhaps additional dlls are required in that directory, e.g.:
brainwin.plugin
libgomp-1.dll
libwinpthread-1.dll
libgcc_s_seh-1.dll
r(999);

* Example from help file for user-written brain command
clear 
set obs 200
gen x = 4*_pi/200*_n
gen y = sin(x)
brain define, input(x) output(y) hidden(20)
quietly brain train, iter(200) eta(2)
brain think ybrain
sort x
twoway (scatter y x) (lfit y x) (line y x)
*/

********** 6. PREDICTION EXAMPLE

* Data for prediction example: 5 continuous and 14 binary variables
qui use mus203mepsmedexp.dta, clear
keep if ltotexp != .
global xlist income educyr age famsze totchr
global dlist suppins female white hisp marry northe mwest south ///
    msa phylim actlim injury priolist hvgg
global rlist c.($xlist)##c.($xlist) i.($dlist) c.($xlist)#i.($dlist)

* Summary statistics for full sample
summarize ltotexp $xlist $dlist

* OLS for full sample
regress ltotexp $xlist $dlist, vce(robust) noheader

// Not included - find model degrees of freedom
ereturn list

* Split the sample with 80% in training sample
splitsample ltotexp, generate(train) split(1 4) values(0 1) rseed(10101)
tabulate train

* OLS with 19 regressors
regress ltotexp $xlist $dlist if train==1, noheader vce(robust)
qui predict y_small

* OLS with 188 potential regressors and 104 estimated 
qui regress ltotexp $rlist if train==1
qui predict y_full

// Not included - find model degrees of freedom
ereturn list

* LASSO with 188 potential regressors leads to 32 selected
qui lasso linear ltotexp $rlist if train==1, selection(adaptive) ///
    rseed(10101) nolog
lassoknots
qui predict y_laspen                  // use penalized coefficients
qui predict y_laspost, postselection  // use post selection OLS coeffs

* Principal components using the first 5 principal components of 19 variables
qui pca $xlist $dlist if train==1
qui predict pc* 
qui regress ltotexp pc1-pc5 if train==1
qui predict y_pca

/* Following does not work so drop
* Neural network with 19 variables and one hidden layers with 10 units
brain define, input($xlist $dlist) output(ltotexp) hidden(10)
qui brain train if train==1, iter(500) eta(2)
brain think y_neural
*/    

* Random forest with 19 variables
qui rforest ltotexp $xlist $dlist if train==1, ///
    type(reg) iter(200) depth(10) lsize(5)
qui predict y_ranfor	

/* Boost requires questionable add-on so drop 
capture program drop boost_plugin
* Boosting linear regression with 19 variables
program boost_plugin, plugin using("C:\ado\personal\boost64.dll")
qui boost ltotexp $xlist $dlist if train==1, ///
    distribution(normal) trainfraction(0.8) maxiter(100) predict(y_boost)  
*/

* Training MSE and test MSE for the various methods
qui regress ltotexp
qui predict y_noreg
foreach var of varlist y_noreg y_small y_full y_laspen y_laspost y_pca ///
                       y_ranfor {
    qui gen `var'errorsq = (`var' - ltotexp)^2
    qui sum `var'errorsq if train == 1
    scalar mse`var'train = r(mean)
    qui sum `var'errorsq if train == 0
    qui scalar mse`var'test = r(mean)
    display "Predictor: " "`var'" _col(21) ///
	    " Train MSE = " %5.3f mse`var'train "  Test MSE = " %5.3f mse`var'test 
    }

********** CLOSE OUTPUT **********

* log close
* clear 
* exit
