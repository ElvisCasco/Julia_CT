* trpanel1.do  September 2013 for Stata 12.0
* Based on trmus0809panel.do based on mus08p1panlin.do
* Earlier version using Stata 10.0 gave wrong se's for xtreg, re vce(robust)

********** OVERVIEW OF trpanel1.do **********

log using trpanel1.txt, text replace

* Based on "Microeconometrics using Stata" 
* by A. Colin Cameron and Pravin K. Trivedi (2008)

* 8.3: PANEL DATA SUMMARY
* 8.4: POOLED OR POPULATION-AVERAGED ESTIMATORS
* 8.5: WITHIN ESTIMATOR
* 8.6: BETWEEN ESTIMATOR
* 8.7: RANDOM EFFECTS ESTIMATOR
* 8.8: COMPARISON OF ESTIMATORS 
* 8.9: FIRST DIFFERENCE ESTIMATOR

* To run you need files
*   mus08psidextract.dta   
* in your directory

* Stata add-on
*   xtoverid
* is used

* To run faster drop the bootstraps at the end or rduce the number of bootstraps

********** SETUP **********

set more off
version 12.0  
set mem 10m
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Source: The Panel Study of Income Dynamics
* Same as Stata website file psidextract.dta
* Data due to  Baltagi and Khanti-Akom (1990) 
* This is corrected version of data in Cornwell and Rupert (1988).
* 595 individuals for years 1976-82.

******* PANEL DATA SUMMARY

* Read in data set
use mus08psidextract.dta, clear
drop tdum*

* Describe dataset
describe

* Summarize dataset 
summarize

* Organization of data set
list id t exp wks occ in 1/3, clean

* Declare individual identifier and time identifier
xtset id t

* Panel description of data set
xtdescribe

* Panel summary statistics: within and between variation
xtsum lwage exp ed t 

* Panel tabulation for a variable
xttab south

* Transition probabilities for a variable
xttrans south, freq

* Simple time series plot for each of 10 individuals
* xtline lwage if id<=10, overlay

* First-order autocorrelation in a variable
sort id t  
correlate lwage L.lwage L2.lwage L3.lwage L4.lwage L5.lwage L6.lwage

* The previous uses few observations as needs L6.lwage
* The following uses as much data as is available
pwcorr lwage L.lwage L2.lwage L3.lwage L4.lwage L5.lwage L6.lwage

* Following does the same - this was used in the Stata book
forvalues j = 1/6 {
     quietly corr lwage L`j'.lwage
     display "Autocorrelation at lag `j' = " %6.3f r(rho) 
     }

* First-order autocorrelation differs in different year pairs
forvalues s = 2/7 {
     quietly corr lwage L1.lwage if t == `s'
     display "Autocorrelation at lag 1 in year `s' = " %6.3f r(rho) 
     }

* Pooled OLS with incorrect default standard errors
regress lwage exp exp2 wks ed, noheader

* Pooled OLS with cluster-robust standard errors
regress lwage exp exp2 wks ed, noheader vce(cluster id)

******* POOLED OR POPULATION-AVERAGED ESTIMATORS

* Pooled feasible GLS estimator with AR(2) error and cluster-robust se's
xtreg lwage exp exp2 wks ed, pa corr(ar 2) vce(robust) nolog

* xtgee = xtreg, pa since xtgee has linear model as default and cluster-robust se's
* xtgee lwage exp exp2 wks ed, corr(ar 2) vce(robust) molog

*******  WITHIN ESTIMATOR

* Within or FE estimator with cluster-robust standard errors
xtreg lwage exp exp2 wks ed, fe vce(robust)

* FE model fitted as LSDV using areg with cluster-robust standard errors
* This gives different standard errors to xtreg, fe
* due to different degrees of freedom correction
areg lwage exp exp2 wks ed, absorb(id) vce(cluster id)

* FE model fitted using LSDV using regress with cluster-robust standard errors
* This gives same standard errors as areg
set matsize 800
quietly xi: regress lwage exp exp2 wks ed i.id, vce(cluster id)
estimates table, keep(exp exp2 wks ed _cons) b se b(%12.7f)

* FE model fitted by add mean of x as a regressor with cluster-robust se's
* This has wrong degrees of freedom as it does not allow for hte N dummies
* It uses NT-k not NT-N-k
global xlist exp exp2 wks ed
sort id
foreach x of varlist $xlist {
  by id: egen mean`x' = mean(`x')
  }
regress lwage exp exp2 wks ed mean*, vce(robust)
drop mean*

*******  BETWEEN ESTIMATOR

* Between estimator with default (nonrobust) standard errors
xtreg lwage exp exp2 wks ed, be

*******  RANDOM EFFECTS ESTIMATOR

* Random effects estimator with cluster-robust standard errors
xtreg lwage exp exp2 wks ed, re vce(robust) theta 

*******  ESTIMATOR COMPARISON

* Compare various estimators all with cluster-robust se's except BE 
global xlist exp exp2 wks ed 
quietly regress lwage $xlist, vce(cluster id)
estimates store OLS
quietly xtgee lwage exp exp2 wks ed, corr(ar 2) vce(robust)
estimates store PFGLS
quietly xtreg lwage $xlist, be
estimates store BE
quietly xtreg lwage $xlist, re vce(robust)
estimates store RE
quietly xtreg lwage $xlist, fe vce(robust)
estimates store FE
estimates table OLS PFGLS BE RE FE,  b(%9.4f) se stats(N)

* Compare default versus cluster-robust se's
quietly xtreg lwage $xlist, re
estimates store RE_def
quietly xtreg lwage $xlist, fe
estimates store FE_def
estimates table RE_def RE FE_def FE,  b(%9.4f) se stats(N)

* Wrong Hausman test requires that RE estimator is fully efficient under null hypothesis
hausman FE_def RE_def, sigmamore

* Correct Robust Hausman test for RE estimator not fully efficient using method of Wooldridge (2002)

* From Cameron and Trivedi (2010) the hard way
quietly xtreg lwage $xlist, re vce(cluster id)
scalar theta = e(theta)
global yandxforhausman lwage exp exp2 wks
sort id
foreach x of varlist $yandxforhausman {
  by id: egen mean`x' = mean(`x')
  generate md`x' = `x' - mean`x'
  generate red`x' = `x' - theta*mean`x'
  }
quietly regress redlwage redexp redexp2 redwks mdexp mdexp2 mdwks, vce(cluster id)
test mdexp mdexp2 mdwks

* Simpler with same result is to add the means as regressors
quietly regress lwage exp exp2 wks meanexp meanexp2 meanwks, vce(cluster id)
test meanexp meanexp2 meanwks
drop mean*

* Or simpler use xtoverid, but this gives a different result for some reason
quietly xtreg lwage exp exp2 wks, re 
xtoverid, cluster(id)

* Prediction after OLS and RE estimation
quietly reg lwage exp exp2 wks ed, vce(cluster id)
predict xbols, xb
quietly xtreg lwage exp exp2 wks ed, re  
predict xbre, xb
predict xbure, xbu
summarize lwage xbols xbre xbure
correlate lwage xbols xbre xbure

******* FIRST DIFFERENCE ESTIMATOR

sort id t
* First-differences estimator with cluster-robust standard errors
* To compare with FE estimates directly use noconstant option
regress D.(lwage exp exp2 wks ed), vce(cluster id) noconstant

******* PANEL BOOTSTRAPS

* Following cluster(id) bootstrap does not work
/*
. bootstrap _b, reps(400) seed(10101) cluster(id): ///
    regress lwage exp exp2 wks ed
(running regress on estimation sample)
Bootstrap replications (400)
----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 
repeated time values within panel
the most likely cause for this error is misspecifying 
 the cluster(), idcluster(), or group() option
r(451);
*/

* Does work if drop the time variable in xtset
xtset id

* OLS panel bootstrap using regress and cluster id
quietly bootstrap _b, reps(400) seed(10101) cluster(id) nodots: ///
   regress lwage exp exp2 wks ed
estimates store OLSboot1

* OLS panel bootstrap using xtreg, pa and vce(boot)
xtreg lwage exp exp2 wks ed, pa corr(ind) vce(boot, reps(400) ///
  seed(10101) nodots)
estimates store OLSboot2

* RE panel bootstrap using bootstrap: and cluster(id)
bootstrap _b, reps(400) seed(10101) cluster(id) nodots: ///
   xtreg lwage exp exp2 wks ed, re
estimates store REboot1

* RE panel bootstrap using vce(boot)
xtreg lwage exp exp2 wks ed, re vce(boot, reps(400) seed(10101) nodots)
estimates store REboot2

* RE panel bootstrap using bootstrap: with cluster(id) and idcluster
* Need to add idcluster for bootstrap of FE
bootstrap _b, reps(400) seed(10101) cluster(id) idcluster(newid) ///
    nodots: xtreg lwage exp exp2 wks ed, fe
estimates store FEboot1

* FE panel bootstrap using vce(boot)
xtreg lwage exp exp2 wks ed, fe vce(boot, reps(400) seed(10101) nodots)
estimates store FEboot2

* Expect bootstrap using method 1 to equal bootstrap using method 2
estimates table OLSboot1 OLSboot2 REboot1 REboot2 FEboot1 FEboot2, b(%9.4f) se stats(N)


********** CLOSE OUTPUT
* log close
* clear
* exit

