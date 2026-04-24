* ML_2022_part2.do  April 2021 2020 for Stata version 16

capture log close

log using ML_2022_part2.txt, text replace

********** OVERVIEW OF ML_2022_part2.do **********

* To run you need files
*  none as data is generated
* in your directory

* And Stata user-written commands 
*   none
* are used

* 3. SHRINKAGE ESTIMATION THEORY: RIDGE, LASSO, ELASTICNET 
* 4. PREDICTION USING LASSO, RIDGE AND ELSTICNET
* 4.1 THE LASSO COMMAND
* 4.2 LASSO LINEAR REGRESSION EXAMPLE
* 4.3 LASSO POSTESTIMATION COMMANDS EXAMPLE
* 4.4 ADAPTIVE LASSO
* 4.5 ELASTICNET AND RIDGE REGRESSION
* 4.6 COMPARISON OF SHRINKAGE ESTIMATORS
* 4.7 SHRINKAGE FOR LOGIT, PROBIT AND POISSON MODELS

********** SETUP **********

set more off
version 16
clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Data are generated
* But are nonetheless saved in case want to use a different program than Stata 

********** GENERATED DATA

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

********** 3: SHRINKAGE ESTIMATION

* Standardize regressors and demean y
foreach var of varlist x1 x2 x3 {
     qui egen double z`var' = std(`var')
}
qui summarize y
qui generate double ydemeaned = y - r(mean)
summarize ydemeaned z*

********** 4.2 LASSO LINEAR REGRESSION EXAMPLE

// Long output so only is included 
* Lasso linear using 5-fold cross validation
lasso linear y x1 x2 x3, selection(cv) folds(5) rseed(10101)

********** 4.3 LASSO POSTESTIMATION COMMANDS EXAMPLE

* List the values of lambda at which variables are added or removed
lassoknots

// Not included
// lasso linear y x1 x2 x3, selection(none) 
// lassoknots, display(bic)
// lassoknots, display(bic) alllambdas

* Plot the change in the penalized objective function as lambda changes
cvplot, saving(graph1, replace)

* Plot how estimated coefficients change with lambda
coefpath, xunits(rlnlambda) saving(graph2, replace)

graph combine graph1.gph graph2.gph, iscale(1.25) ysize(2.5) xsize(6.0)
* graph export mus228fig1cvplots.eps, replace

* Provide a summary of the lasso
lassoinfo

* Lasso coefficients for the standardized regressors
lassocoef, display(coef, standardized)

*  Lasso coefficients for the unstandardized regressors
lassocoef, display(coef, penalized) nolegend

* Post-selection estimated coefficients for the unstandardized regressors
lassocoef, display(coef, postselection) nolegend

* Goodness-of-fit with penalized coefficients and postselection coefficients
lassogof, penalized
lassogof, postselection

* Compare to OLS with the lasso selected regressors 
regress y x1 x2, noheader

********* * 4.4 ADAPTIVE LASSO

* Lasso linear using 5-fold adaptive cross validation
qui lasso linear y x1 x2 x3, selection(adaptive) folds(5) rseed(10101)
lassoknots

* Lasso linear with no method for selecting lambda
qui lasso linear y x1 x2 x3, selection(none) folds(5)
lassoknots

// Not included - Lasso linear with plugin lambda
lasso linear y x1 x2 x3, selection(plugin) folds(5)

********* 4.5 ELASTICNET AND RIDGE REGRESSION

// Not included - Ridge with complete grid datasets
elasticnet linear y x1 x2 x3, alpha(0) rseed(10101)

* Ridge estimation using the elasticnet command and selected results
qui elasticnet linear y x1 x2 x3, alpha(0) rseed(10101) folds(5)
lassoknots
lassocoef, display(coef, penalized) nolegend
lassogof, penalized

* Elastic net estimation and selected results
qui elasticnet linear y x1 x2 x3, alpha(0.9(0.05)1) rseed(10101) folds(5)
lassoknots
lassocoef, display(coef, penalized) nolegend
lassogof, penalized

// Not included - Elastic net with complete grid datasets
elasticnet linear y x1 x2 x3, alpha(0.1(0.3)1) rseed(10101) folds(5)

********* 4.6 COMPARISON OF SHRINKAGE ESTIMATORS

* Estimate various models and store results
qui regress y x1 x2 x3
estimates store OLS
qui lasso linear y x1 x2 x3, selection(cv) folds(5) rseed(10101)
estimates store LASCV
qui lasso linear y x1 x2 x3, selection(adaptive) folds(5) rseed(10101)
estimates store LASADAPT
qui lasso linear y x1 x2 x3, selection(plugin) folds(5)
estimates store LASPLUG
qui elasticnet linear y x1 x2 x3, alpha(0) selection(cv) folds(5) rseed(10101)
estimates store RIDGECV
qui elasticnet linear y x1 x2 x3, alpha(0.9(0.05)1) rseed(10101) folds(5)
estimates store ELASTIC

* Compare in-sample fit and selected coefficients of various models
lassogof OLS LASCV LASADAPT LASPLUG RIDGECV ELASTIC
lassocoef OLS LASCV LASADAPT LASPLUG RIDGECV ELASTIC, display(coef) nolegend  

********** 4.7 SHRINKAGE FOR LOGIT, PROBIT AND POISSON MODELS

* Lasso for logit example
qui generate dy = y > 3
qui lasso logit dy x1 x2 x3, rseed(10101) folds(5)
lassoknots

* Lasso for count data example
qui generate ycount = rpoisson(exp(-1 + x1)) 
qui lasso poisson ycount x1 x2 x3, rseed(10101) folds(5)
lassoknots

********** CLOSE OUTPUT **********

* log close
* clear 
* exit
