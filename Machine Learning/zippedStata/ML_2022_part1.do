* ML_2022_part1.do  April 2021 2020 for Stata version 16

capture log close

log using ML_2022_part1.txt, text replace

********** OVERVIEW OF ML_2022_part1.do **********

* To run you need files
*  none as data is generated
* in your directory

* And Stata user-written commands 
*   crossfold, loocv, vselect
* are used

* 2.1 GENERATED DATA SAMPLE
* 2.3 MSE, INFORMATION CRITERIA AND RELATED PENALTY MEASURES
* 2.4 SPLIT SAMPLE COMMAND
* 2.5 SINGLE SPLIT CROSS VALIDATION
* 2.6 K-FOLD CROSS VALIDATION
* 2.7 LEAVE-ONE-OUT CROSS VALIDATION
* 2.8 BEST SUBSETS AND STEPWISE SELECTION
* 2.9 SELECTION USING STATISTICAL SIGNIFICANCE

********** SETUP **********

set more off
version 16
clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Data are generated
* But are nonetheless saved in case want to use a different program than Stata 

********** 2.1 GENERATED DATA

* Generate three correlated variables (rho = 0.5) and y linear only in x1
clear
quietly set obs 40
set seed 12345
matrix MU = (0,0,0)
scalar rho = 0.5
matrix SIGMA = (1,rho,rho \ rho,1,rho \ rho,rho,1)
drawnorm x1 x2 x3, means(MU) cov(SIGMA)
generate y = 2 + 1*x1 + rnormal(0,3)
saveold ML_2021_part1, version(11) replace

* Summarize data
summarize
correlate

* OLS regression of y on x1-x3
regress y x1 x2 x3, vce(robust)

regress y x1 
di e(rss) "  " e(rss)/e(N)  "  " e(N)

regress y x1 x3 
di e(rss) "  " e(rss)/e(N)  "  " e(N)


*********** 2.3 MSE, INFORMATION CRITERIA AND RELATED PENALTY MEASURES2.4 MODEL SELECTION BASED ON PENALIZED GOODNESS-OF-FIT

* Regressor lists for all possible models
global xlist1
global xlist2 x1
global xlist3 x2
global xlist4 x3
global xlist5 x1 x2
global xlist6 x2 x3
global xlist7 x1 x3
global xlist8 x1 x2 x3

* Full sample estimates with AIC, BIC, Cp, R2adj penalties
quietly regress y $xlist8
scalar s2full = e(rmse)^2  // Needed for Mallows Cp
forvalues k = 1/8 {
    quietly regress y ${xlist`k'}
    scalar mse`k' = e(rss)/e(N)
    scalar r2adj`k' = e(r2_a)
    scalar aic`k' = -2*e(ll) + 2*e(rank)
    scalar bic`k' = -2*e(ll) + e(rank)*ln(e(N))
    scalar cp`k' =  e(rss)/s2full - e(N) + 2*e(rank)
    display "Model " "${xlist`k'}" _col(15) " MSE=" %8.5f mse`k'  ///
     " R2adj=" %6.3f r2adj`k' "  AIC=" %7.2f aic`k'  ///
     " BIC=" %7.2f bic`k' " Cp=" %6.3f cp`k'
}

********** 2.4 SPLITSAMPLE COMMAND

* Split sample into five equal size parts using splitsample command
splitsample, nsplit(5) generate(snum) rseed(10101)
tabulate snum

********** 2.5 SINGLE SPLIT CROSS VALIDATION

* Split into training group with 32 observations (dtrain = 1) 
* and validation sample with 8 observations (dtrain = 0)
splitsample, split(1 4) values(0 1) generate(dtrain) rseed(10101)
tabulate dtrain

* Single split validation - training and test MSE for the 8 possible models
forvalues k = 1/8 {
    qui reg y ${xlist`k'} if dtrain==1
    qui predict y`k'hat
    qui gen y`k'errorsq = (y`k'hat - y)^2
    qui sum y`k'errorsq if dtrain == 1
    scalar mse`k'train = r(mean)
    qui sum y`k'errorsq if dtrain == 0
    qui scalar mse`k'test = r(mean)
    display "Model " "${xlist`k'}" _col(16)  ///
        " Training MSE = " %7.3f mse`k'train " Test MSE = " %7.3f mse`k'test 
}
drop y*hat y*errorsq 

********** 2.6 K-FOLD CROSS VALIDATION

* Five-fold cross validation example for model with all regressors
splitsample, nsplit(5) generate(foldnum) rseed(10101)
matrix allmses = J(5,1,.)
capture drop y*hat y*errorsq
forvalues i = 1/5 {
    qui reg y x1 x2 x3 if foldnum != `i'
    qui predict y`i'hat
    qui gen y`i'errorsq = (y`i'hat - y)^2
    qui sum y`i'errorsq if foldnum ==`i'
    matrix allmses[`i',1] = r(mean)
}
matrix list allmses

* Compute the average MSE over the five folds and standard deviation 
svmat allmses, names(vallmses) 
qui sum vallmses
display "CV5 = " %5.3f r(mean) " with st. dev. = " %5.3f r(sd)

* Five-fold cross validation measure for one model with ll regressors 
set seed 10101
crossfold regress y x1 x2 x3, k(5)
drop _est*  // Drop variables created 

* Five-fold cross validation measure for all possible models
forvalues k = 1/8 {
    set seed 10101
    qui crossfold regress y ${xlist`k'}, k(5)  
    matrix RMSEs`k' = r(est)
    svmat RMSEs`k', names(rmse`k') 
    qui generate mse`k' = rmse`k'^2
    qui sum mse`k'
    scalar cv`k' = r(mean)
    scalar sdcv`k' = r(sd)
    display "Model " "${xlist`k'}" _col(16) "  CV5 = " %7.3f cv`k' ///
        " with st. dev. = " %7.3f sdcv`k'
}

********** 2.7 LEAVE-ONE-OUT CROSS VALIDATION

* Leave-one-out cross validation (loocv sets same seed each time)
loocv regress y x1 
display "LOOCV MSE = " r(rmse)^2

* Not included
loocv regress y x1 x2
loocv regress y x1 x2 x3

********** 2.8 BEST SUBSETS SELECTION AND STEPWISE REGRESSION

* Best subset selection with user-written add-on vselect
vselect y x1 x2 x3, best

* Best subset selection with user-written add-on vselect
vselect y x1 x2 x3, best
* Stepwise forwards using AIC
vselect y x1 x2 x3, forward aic

* Not included
* Stepwise backwards using AIC
vselect y x1 x2 x3, backward aic

* Not included
* Best subsets with x1 always included
vselect y x2 x3, fix(x1) best

* Not included
* Add-on command gvselect for OLS regression with x1 always included
* Problem here
* gvselect <xlist> x2 x3: regress y <xlist> x1

********** 2.9 SELECTION USING STATISTICAL SIGNIFICANCE

* This is old school use same p = 0.05 regardless of number of regressors.
* Recent work says reduce p-value with number of potential regressors

* Stepwise forward using statistical significance at five percent
stepwise, pe(.05): regress y x1 x2 x3

* Stepwise backward using statistical significance at five percent
stepwise, pr(.05): regress y x1 x2 x3

* Not included
* Stepwise forward in specified order (hierarchical) testing at five percent
stepwise, pe(.05): regress y x1 x2 x3

********** CLOSE OUTPUT **********

* log close
* clear 
* exit



