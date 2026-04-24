* mus228prediction.do  for Stata 17
capture log close

********** OVERVIEW OF mus228prediction.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 28  MACHINE LEARNING FOR PREDICTION AND INFERENCE
* 28.2: MEASURING THE PREDICTIVE ABILITY OF A MODEL
* 28.3: SHRINKAGE ESTIMATORS
* 28.5: DIMENSION REDUCTION
* 28.7: PREDICTION APPLICATION
* 28.8: MACHINE LEARNING FOR INFERENCE IN PARTIAL LINEAR MODEL

* To run you need files
*   mus203mepsmedexp.dta
*   mus228ajr.dta
* in your directory

* And community-contributed commands 
*   vselect
*   crossfold
*   loocv 
*   rforest
*   boost
*   program boost_plugin, plugin using("C:\ado\personal\boost64.dll")
* are used

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Most data are generated

* File mus203mepsmedexp.dta is aothurs' extract from MEPS 
* (Medical Expenditure Panel Survey)
* for individuals 65 years and older in U.S. Medicare in 2003

* mus228ajr.dta is from Acemoglu, Johnson and Robinson (2001), "The colonial 
* origins of comparative development: an empirical investigation," AER, 1369-1401.

********** 28.2: MEASURING PREDICTIVE ABILITY

* Generate three correlated variables (rho = 0.5) and y linear only in x1
qui set obs 40
set seed 12345
matrix MU = (0,0,0)
scalar rho = 0.5
matrix SIGMA = (1,rho,rho \ rho,1,rho \ rho,rho,1)
drawnorm x1 x2 x3, means(MU) cov(SIGMA)
generate y = 2 + 1*x1 + rnormal(0,3)

* Summarize data
summarize
correlate

* OLS regression of y on x1-x3
regress y x1 x2 x3, vce(robust)
   
* Regressor lists for all possible models
global xlist1
global xlist2 x1
global xlist3 x2
global xlist4 x3
global xlist5 x1 x2
global xlist6 x2 x3
global xlist7 x1 x3
global xlist8 x1 x2 x3

* Full-sample estimates with AIC, BIC, Cp, R2adj penalties
qui regress y $xlist8
scalar s2full = e(rmse)^2  // Needed for Mallows Cp
forvalues k = 1/8 {
    qui regress y ${xlist`k'}
    scalar mse`k' = e(rss)/e(N)
    scalar r2adj`k' = e(r2_a)
    scalar aic`k' = -2*e(ll) + 2*e(rank)
    scalar bic`k' = -2*e(ll) + e(rank)*ln(e(N))
    scalar cp`k' =  e(rss)/s2full - e(N) + 2*e(rank)
    display "Model " "${xlist`k'}" _col(15) " MSE=" %6.3f mse`k'  ///
        " R2adj=" %6.3f r2adj`k' "  AIC=" %7.2f aic`k'  ///
        " BIC=" %7.2f bic`k' " Cp=" %6.3f cp`k'
}	

* Split sample into five equal-size parts using splitsample command
splitsample, nsplit(5) generate(snum) rseed(10101)
tabulate snum

* Form indicator for training data (80% of sample) and test data (20%)
splitsample, split(1 4) values(0 1) generate(dtrain) rseed(10101)
tabulate dtrain

* Single-split validation - training and test MSE for the 8 possible models
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

* Five-fold CV example for model with all regressors
splitsample, nsplit(5) generate(foldnum) rseed(10101)
matrix allmses = J(5,1,.)
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
qui sum vallmses1
display "CV5 = " %5.3f r(mean) " with st. dev. = " %5.3f r(sd)

/* 

tempvar k
scalar `k' = 1
di `k'

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

*/		


* Five-fold CV measure for all possible models
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

* LOOCV
loocv regress y x1 
display "LOOCV MSE = " r(rmse)^2

* Not included
loocv regress y x1 x2
loocv regress y x1 x2 x3

* Best subset selection with community-contributed command vselect
vselect y x1 x2 x3, best

********** 28.4: SHRINKAGE ESTIMATORS

* Standardize regressors and demean y
foreach var of varlist x1 x2 x3 {
     qui egen double z`var' = std(`var')
}
qui summarize y
qui generate double ydemeaned = y - r(mean)
summarize ydemeaned z*

// Long output so only is included 
* Lasso linear using 5-fold CV
lasso linear y x1 x2 x3, selection(cv, folds(5)) rseed(10101)

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

* Provide a summary of the lasso
lassoinfo

* Lasso coefficients for the standardized regressors
lassocoef, display(coef, standardized)

*  Lasso coefficients for the unstandardized regressors
lassocoef, display(coef, penalized) nolegend

* Postselection estimated coefficients for the unstandardized regressors
lassocoef, display(coef, postselection) nolegend

* Goodness of fit with penalized coefficients and postselection coefficients
lassogof, penalized
lassogof, postselection

* Compare with OLS with the lasso-selected regressors 
regress y x1 x2, noheader

* Lasso linear using 5-fold adaptive CV
qui lasso linear y x1 x2 x3, selection(adaptive, folds(5)) rseed(10101)
lassoknots

* Lasso linear with no method for selecting lambda
qui lasso linear y x1 x2 x3, selection(none) 
lassoknots

// Not included - Lasso linear with plugin lambda
lasso linear y x1 x2 x3, selection(plugin) folds(5)

// Not included - Ridge with complete grid datasets
elasticnet linear y x1 x2 x3, alpha(0) rseed(10101)

* Ridge estimation using the elasticnet command and selected results
qui elasticnet linear y x1 x2 x3, alpha(0) rseed(10101) selection(cv, folds(5))
lassoknots
lassocoef, display(coef, penalized) nolegend
lassogof, penalized

* Elastic net estimation and selected results
qui elasticnet linear y x1 x2 x3, alpha(0.9(0.05)1) rseed(10101) selection(cv, folds(5))
lassoknots
lassocoef, display(coef, penalized) nolegend
lassogof, penalized

// Not included - Elastic net with complete grid datasets
elasticnet linear y x1 x2 x3, alpha(0.1(0.3)1) rseed(10101) folds(5)

* Fit various models and store results
qui regress y x1 x2 x3
estimates store OLS
qui lasso linear y x1 x2 x3, selection(cv, folds(5)) rseed(10101)
estimates store LASCV
qui lasso linear y x1 x2 x3, selection(adaptive, folds(5)) rseed(10101)
estimates store LASADAPT
qui lasso linear y x1 x2 x3, selection(plugin) 
estimates store LASPLUG
qui elasticnet linear y x1 x2 x3, alpha(0) selection(cv, folds(5)) rseed(10101)
estimates store RIDGECV
qui elasticnet linear y x1 x2 x3, alpha(0.9(0.05)1) rseed(10101) selection(cv, folds(5))
estimates store ELASTIC

* Compare in-sample fit and selected coefficients of various models
lassogof OLS LASCV LASADAPT LASPLUG RIDGECV ELASTIC
lassocoef OLS LASCV LASADAPT LASPLUG RIDGECV ELASTIC, display(coef) nolegend  

* Lasso for logit example
qui generate dy = y > 3
qui lasso logit dy x1 x2 x3, rseed(10101) selection(cv, folds(5))
lassoknots

* Lasso for count data example
qui generate ycount = rpoisson(exp(-1 + x1)) 
qui lasso poisson ycount x1 x2 x3, rseed(10101) selection(cv, folds(5))
lassoknots

********** 28.5: DIMENSION REDUCTION

capture drop pc* yhat
* Principal components using default option that first standardizes the data
pca x1 x2 x3

* Not included - get same results using standardized data and covariance option
pca zx1 zx2 zx3, covariance

* Compute the three principal components and their means, st.devs., correlations
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

********** 28.7: PREDICTION EXAMPLE USING MACHINE LEARNING

* Data for prediction example: 5 continuous and 14 binary variables
qui use mus203mepsmedexp, clear
keep if !missing(ltotexp)
global xlist income educyr age famsze totchr
global dlist suppins female white hisp marry northe mwest south ///
    msa phylim actlim injury priolist hvgg
global rlist c.($xlist)##c.($xlist) i.($dlist) c.($xlist)#i.($dlist)

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
qui predict y_laspen                  // Use penalized coefficients
qui predict y_laspost, postselection  // Use post selection OLS coeffs

* Principal components using the first 5 principal components of 19 variables
qui pca $xlist $dlist if train==1
qui predict pc* 
qui regress ltotexp pc1-pc5 if train==1
qui predict y_pca

* Neural network with 19 variables and 2 hidden layers each with 10 units
brain define, input($xlist $dlist) output(ltotexp) hidden(10)
qui brain train if train==1, iter(500) eta(2)
brain think y_neural    

* Random forest with 19 variables
qui rforest ltotexp $xlist $dlist if train==1, ///
    type(reg) iter(200) depth(10) lsize(5)
qui predict y_ranfor	

capture program drop boost_plugin

* Boosting linear regression with 19 variables
program boost_plugin, plugin using("C:\ado\personal\boost64.dll")
qui boost ltotexp $xlist $dlist if train==1, ///
    distribution(normal) trainfraction(0.8) maxiter(100) predict(y_boost)  

* Training MSE and test MSE for the various methods
qui regress ltotexp
qui predict y_noreg
foreach var of varlist y_noreg y_small y_full y_laspen y_laspost y_pca ///
                       y_neural y_ranfor y_boost {
    qui gen `var'errorsq = (`var' - ltotexp)^2
    qui sum `var'errorsq if train == 1
    scalar mse`var'train = r(mean)
    qui sum `var'errorsq if train == 0
    qui scalar mse`var'test = r(mean)
    display "Predictor: " "`var'" _col(21) ///
	    " Train MSE = " %5.3f mse`var'train "  Test MSE = " %5.3f mse`var'test 
    }

********** 28.8: INFERENCE USING MACHINE LEARNING

* Data for inference on suppins example: 5 continuous and 13 binary variables
qui use mus203mepsmedexp, clear
keep if ltotexp != .
global xlist2 income educyr age famsze totchr
global dlist2 female white hisp marry northe mwest south ///
    msa phylim actlim injury priolist hvgg
global rlist2 c.($xlist2)##c.($xlist2) i.($dlist2) c.($xlist2)#i.($dlist2)

* OLS on small model and full model
qui regress ltotexp suppins $xlist2 $dlist2, vce(robust)
estimates store OLSSMALL
qui regress ltotexp suppins $rlist2, vce(robust)
estimates store OLSFULL
estimates table OLSSMALL OLSFULL, keep(suppins) b(%9.4f) se stats(N df_m r2)

* Partialing-out partial linear model using default plugin lambda
poregress ltotexp suppins, controls($rlist2)

* Lasso information
lassoinfo

* Cluster-robust partialing-out partial linear model using default plugin lambda
poregress ltotexp suppins, controls($rlist2) vce(cluster age)

* Lasso information
lassoinfo

lassocoef (., for(ltotexp))
lassocoef (., for(suppins))

// Not included - various postestimation options
poregress ltotexp suppins, controls($rlist2) selection(cv) rseed(10101)
lassoinfo
lassocoef (., for(ltotexp))
lassocoef (., for(suppins))
lassoknots, for(ltotexp)
lassoknots, for(suppins)
cvplot, for(ltotexp)
cvplot, for(suppins)
coefpath, for(ltotexp)
// lassogof does not apply after po, xpo, ds 

* Partialing out done manually
qui lasso linear suppins $rlist2, selection(plugin) 
qui predict suppins_lasso, postselection
qui generate u_suppins = suppins - suppins_lasso
qui lasso linear ltotexp $rlist2, selection(plugin) 
qui predict ltotexp_lasso, postselection
qui generate u_ltotexp = ltotexp - ltotexp_lasso
regress u_ltotexp u_suppins, vce(robust) noconstant noheader

* Cross-fit partialing-out (double/debiased) using default plugin
xporegress ltotexp suppins, controls($rlist2) rseed(10101) nolog

* Summarize the number of selected variables across the ten folds
lassoinfo

* Double-selection partial linear model using default plugin
dsregress ltotexp suppins, controls($rlist2)

// Not included - do dsregress manually
qui lasso linear suppins $rlist2, selection(plugin) 
lassocoef
qui lasso linear ltotexp $rlist2, selection(plugin) 
lassocoef
regress ltotexp suppins income age c.income#c.totchr 1.marry#c.income   ///
    0.northe#c.income 1.hvgg#  c.income 1.white#c.educy 0.hisp#c.educyr ///
    0.marry#c.famsze ///
    totchr c.educyr#c.totchr c.age#c.totchr 0.actlim 1.phylim#c.educyr  /// 
    1.priolist#c.educyr 0.phylim#c.famsze 0.actlim#c.famsze             ///
    0.female#c.totchr 1.white#c.totchr 0.hisp#c.totchr 0.hvgg#c.totchr
dsregress ltotexp suppins, controls($rlist2)

* Exponential variant of partial linear model and partialing-out estimator
generate ycount = floor(sqrt(totexp/500))
summarize ycount
qui poisson ycount suppins $xlist2 $dlist2, vce(robust)
estimates store PSMALL
qui poisson ycount suppins $rlist2, vce(robust)
estimates store PFULL
qui popoisson ycount suppins, controls($rlist2) coef
estimates store PPOLASSO
estimates table PSMALL PFULL PPOLASSO, keep(suppins) b(%9.4f) se ///
    stats(N df_m k_controls_sel)

// Not included - popoisson estimation for the level of health expenditure
poisson totexp suppins $xlist2 $dlist2, vce(robust)
popoisson totexp suppins, controls($rlist2) coef

* Logit variant of partial linear model and partialing-out estimator
generate dy = totexp > 4000
tabulate dy
qui logit dy suppins $xlist2 $dlist2, or vce(robust)
estimates store LSMALL
qui logit dy suppins $rlist2, or vce(robust)
estimates store LFULL
qui pologit dy suppins, controls($rlist2) coef
estimates store LPOLLASSO
estimates table LSMALL LFULL LPOLLASSO, keep(suppins) b(%9.4f) se ///
    stats(N df_m k_controls_sel)

* Read in Acemoglu-Johnson-Robinson data and define globals  
qui use mus228ajr, clear
global xlist lat_abst edes1975 avelf temp* humid* steplow deslow ///
    stepmid desmid drystep  drywint goldm iron silv zinc oilres landlock
describe logpgp95 avexpr logem4
summarize logpgp95 avexpr logem4, sep(0)

* Partialing-out IV using plugin for lambda
poivregress logpgp95 (avexpr=logem4), controls($xlist) selection(plugin, hom)

* poivregress estimator in just-identified model obtained manually
gen y = logpgp95
gen d = avexpr
global zlist logem4
qui lasso linear y $xlist, selection(plugin, hom)    // Lasso of y on x
qui predict yhat, postselection
generate yresid = y - yhat                           // Generate y residual
qui lasso linear d $xlist $zlist, selection(plugin, hom)  // Lasso d on x,z
qui predict dhat, postselection                      // Generate dhat 
qui lasso linear dhat $xlist, selection(plugin, hom) // Lasso dhat on x 
predict dhat_hat, postselection
generate dhatresid = dhat - dhat_hat                 // Generate dhat residual
generate dresid = d - dhat_hat                       // Generate d "residual"
ivregress 2sls yresid (dresid = dhatresid), noconstant vce(robust)

// Not included 
// The following shows selected instrument is logem4
// and selected controls are edes1975 avelf temp2 iron zinc
qui poivregress logpgp95 (avexpr=logem4), controls($xlist) selection(plugin, hom)
lassocoef (., for(avexpr))       // Chooses logem4 and edes1975 zinc
lassocoef (., for(logpgp95))     // Chooses edes1975 avelf
lassocoef (., for(pred(avexpr))) // Chooses edes1975 avelf temp2 iron zinc  

// Not included IV with the LASSO selected variables 
// leads to same estimate in this special just-identified example
ivregress 2sls logpgp95 (avexpr = logem4) edes1975 avelf temp2 iron zinc, ///
    vce(robust) noheader
 
********** END
