* mus204regext.do  for Stata 17

cap log close

********** OVERVIEW OF mus204regext.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 4
* 4.2: IN-SAMPLE PREDICTION 
* 4.3: OUT-OF-SAMPLE PREDICTION 
* 4.4: PREDICTIVE MARGINS
* 4.5: MEs
* 4.6: REGRESSION DECOMPOSITION ANALYSIS
* 4.7: SHAPLEY DECOMPOSITION OF RELATIVE REGRESSOR IMPORTANCE
* 4.8: DIFFERENCE-IN-DIFFERENCE ESTIMATES

* To run you need files
*   mus203mepsmedexp.dta 
*   mus204mabeldecomp.dta
*   mus204nswpsid.dta
* in your directory

* community-contributed commands 
*   oaxaca 
*   rego (net from http://www.marco-sunder.de/stata/  net install rego)
* are used
 
********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* File mus203mepsmedexp is authors' extract from MEPS 
* (Medical Expenditure Panel Survey)
* for individuals 65 years and older in U.S. Medicare in 2003

* File mus204mabeldecomp.dta is authors' extract from MABEL EPS
* Medicine in Australia: Balancing Employment and Life (MABEL)
* An Australian longitudinal survey of physicians.

* File mus204nswpsid.dta is orginally from R.H. Dehejia and S. Wahba (1999) 
* "Causal Effects in Nonexperimental Studies: reevaluating the 
* Evaluation of Training Programs", JASA, 1053-1062
* or R.H. Dehejia and S. Wahba (2002) "Propensity-score Matching Methods 
* for Nonexperimental Causal Studies", ReStat, 151-161

*********** 4.2 IN-SAMPLE PREDICTION

clear all

* OLS in levels for positive medical expenditures
qui use mus203mepsmedexp 
keep if totexp > 0   
regress totexp suppins phylim actlim totchr age female income, vce(robust)

* In-sample prediction following OLS in levels
predict yhatlevels
summarize totexp yhatlevels

* Compare median prediction and sample median actual value after OLS
tabstat totexp yhatlevels, stat(count p50) col(stat)

* In-sample prediction in levels from a logarithmic model
qui regress ltotexp suppins phylim actlim totchr age female income, vce(robust)
qui predict lyhat
generate yhatwrong = exp(lyhat)
generate yhatnormal = exp(lyhat)*exp(0.5*e(rmse)^2)
qui predict uhat, residual
generate expuhat = exp(uhat)
qui summarize expuhat
generate yhatduan = r(mean)*exp(lyhat) 
summarize totexp yhatwrong yhatnormal yhatduan yhatlevels 

* Effect of suppins: (1) compare actual means with suppins==1 and suppins==0
bysort suppins: summarize totexp

* Effect of suppins: (2) compare predicted means by suppins value
bysort suppins: summarize yhatlevels yhatduan

* Effect of suppins: (3) predict with all set to suppins==0 versus all set to 1
qui regress totexp suppins phylim actlim totchr age female income, vce(robust)
preserve
qui replace suppins = 1
qui predict yhat1
qui replace suppins = 0
qui predict yhat0 
generate treateffect = yhat1 - yhat0
summarize yhat1 yhat0 treateffect 
restore

* Effect of suppins model: (3b) same as previous but easier and gives se 
qui regress totexp i.suppins phylim actlim totchr age female income, vce(robust)
margins, dydx(suppins) nofvlabel

* Effect of suppins model: (4) different fitted models for suppins =0 and =1
qui regress totexp phylim actlim totchr age female income if suppins==0
qui predict yhatmodel0
qui regress totexp phylim actlim totchr age female income if suppins==1
qui predict yhatmodel1
generate treateffectra = yhatmodel1 - yhatmodel0
summarize yhatmodel1 yhatmodel0 treateffectra 

// Not included as output in text but mentioned
// Obtain previous estimate using using teffects ra command
teffects ra (totexp phylim actlim totchr age female income) (suppins), nolog 

* Compute standard errors of prediction and forecast with default VCE
qui regress totexp suppins phylim actlim totchr age female income
display "Estimated standard deviation of the error: s = " e(rmse)
predict yhatstdp, stdp
predict yhatstdf, stdf
summarize yhatlevels yhatstdp yhatstdf

*********** 4.3 OUT-OF-SAMPLE PREDICTION

* Out-of-sample OLS predictions using correctly specified model
clear all
qui set obs 100
set seed 10101
generate u = 3*rnormal()
generate x = 10 + 4*rnormal() 
generate y = 1 + x + u
qui regress y x if _n < 91
qui predict yhatos if _n > 90   // Predict 10 out-of-sample observations
qui regress y x
qui predict yhatis if _n > 90   // Predict 10 in-sample observations 
qui generate ferroros = y - yhatos 
summarize y yhatis yhatos ferroros if _n > 90

* Out-of-sample prediction under misspecification (omitted regressor)    
generate z = 10 + 4*rnormal()
generate ynew =  1 + x + z + u
qui regress ynew x if _n <91     // Regression with variable z omitted
qui predict ynewhatos if _n >90  // Predict 10 out-of-sample observations
qui regress ynew x
qui predict ynewhatis if _n > 90 // Predict 10 in-sample observations 
qui generate ferrornewos = ynew - ynewhatos if _n > 90
list ferrornewos in 91/100, clean
summarize ynew ynewhatis ynewhatos ferrornewos if _n > 90

********** 4.4: PREDICTIVE MARGINS

* Factor variables used to define model with interactions and quadratic
qui use mus203mepsmedexp, clear
regress totexp i.suppins i.phylim i.actlim c.totchr c.age##c.age ///
    i.female##c.income, vce(robust) noheader nofvlabel

* Predictive margin using sample average 
margins

// Not included in text - for OLS equals sample mean for estimation sample
sum totexp if e(sample)==1

* Predictive margin using sample average with population-average inference
margins, vce(unconditional)

* Predictive margins by categorical regressor gender using sample averages
margins female, nofvlabel

// Not included - for OLS margins result for female=1 less that for female=0
// equals the marginal effect from margins, dydx()
margins, dydx(female)

// Output not included in text but discussed
// margins after OLS differs from raw difference in means due to regressors
qui sum totexp if female==1 & e(sample)==1
scalar meanfemale = r(mean)
qui sum totexp if female==0 & e(sample)==1
display "Mean totexp for female = " meanfemale " and for male = " r(mean)

* Predictive margins by categorical regressor gender using sample means 
margins female, atmeans nofvlabel

* Predictive margins by categorical regressor gender at specified values 
margins female, at(suppins=1 phylim=1 actlim=1 totchr=2 age=70 income=50) nofvlabel

// Output not included but mentioned in text - the following fails
//   margins age
// gives error: gfactor 'age' not found in list of covariates  r(322)

* Margins at 3 different ages (averaged over other regressors' sample values)
margins, at(age=(65(10)85))

* Marginsplot for predictive margins by age
qui margins, at(age=(65(5)90))
marginsplot, saving(graph1, replace)

* Marginsplot for predictive margins by age for each gender
margins, at(age=(65(5)90)) over(female) nofvlabel
marginsplot, saving(graph2, replace)

graph combine graph1.gph graph2.gph, iscale(1.2) ///
    ycommon xcommon rows(1) ysize(2.5) xsize(6)

* Marginsplot for predictive margins by age, separate graphs by gender
qui margins, at(age=(65(5)90)) over(female)
marginsplot, by(female) saving(graph1, replace)

graph combine graph1.gph, iscale(1.2) rows(1) ysize(2.5) xsize(6)

* Predictive margins for three income values
margins, at(income=(10(10)30)) 

* Pairwise comparison of predictive margins for three income values
margins, at(income=(10(10)30)) pwcompare(effects)

* Not included - compare by gender
margins female
margins female, pwcompare 
margins, at(female=(0 1)) pwcompare
margins, over(female) pwcompare

* Pairwise comparison of predictive margins for six income-gender combinations
margins female, at(income=(10(10)30)) pwcompare(effects) nofvlabel

* Contrasts of predictive margins for gender
margins female, contrast(ci) nofvlabel

* Not included - ci option gives same as pwcompare
margins female, contrast(ci)
margins female, pwcompare

* Contrasts of predictive margins for gender at each of three income levels
margins female, at(income=(10(10)30)) contrast(nowald) nofvlabel

********** 4.5: MARGINAL EFFECTS

* AMEs computed with factor variables used to define model
qui regress totexp i.suppins i.phylim i.actlim c.totchr c.age##c.age ///
    i.female##c.income, vce(robust) noheader
margins, dydx(*) nofvlabel

* AMEs computed manually
generate meage = _b[age] + 2*_b[c.age#c.age]*age
generate mefemale = _b[1.female] + _b[1.female#c.income]*income
generate meincome = _b[income] + _b[1.female#c.income]*female
sum meage mefemale meincome

* Wrong way to compute (average) MEs if any nonlinearity present
generate  agesq = age^2
generate fembyinc = female*income
qui regress totexp suppins phylim actlim totchr age agesq ///
    female income fembyinc, vce(robust) 
margins, dydx(*)

// Not included as in this example it gives the same results as earlier 
* MER computed with factor variables 
qui regress totexp i.suppins i.phylim i.actlim c.totchr c.age##c.age     ///
    i.female##c.income, vce(robust) 
margins, dydx(*) atmeans noatlegend

* MER computed with factor variables 
qui regress totexp i.suppins i.phylim i.actlim c.totchr c.age##c.age     ///
    i.female##c.income, vce(robust) noheader
margins, dydx(*) at(suppins=1 phylim=1 actlim=1 totchr=2 age=70 female=1 ///
    income=20) nofvlabel

* Sample average elasticity with respect to income
margins, eyex(income)

* Manual computation of sample average elasticity with respect to income
predict yhat
generate eincome = (_b[income]+_b[1.female#c.income]*female) * income / yhat 
sum eincome

* Elasticity with respect to income evaluated at representative values
margins, eyex(income) at(suppins=1 phylim=1 actlim=1 totchr=2 age=70 ///
    female=1 income=20)

* Not included - why separate means for e.g. 0.female and 1.female
margins, eyex(income) atmeans

* Semielasticities: Proportionate change in y w.r.t. unit change in x
margins, eydx(*) at(suppins=1 phylim=1 actlim=1 totchr=2 age=70 ///
    female=1 income=20) nofvlabel

* Semielasticities: Change in y w.r.t. proportionate change in x
margins, dyex(totchr age income) at(suppins=1 phylim=1 actlim=1 totchr=2 ///
    age=70 female=1 income=20)

*********** 4.6: REGRESSION DECOMPOSITION ANALYSIS 

* t tests on difference in mean log-earnings by gender
qui use mus204mabeldecomp, clear
ttest logyearn, by(female) unequal

// Not included as output in text but mentioned - repeat in levels
generate earn = exp(logyearn)
ttest earn, by(female) unequal

* Gender differences in variables for the regression estimation sample
global xlist1 yhrs expr exprsq fellow pgradoth pracsize childu5 visa
qui regress logyearn $xlist1
table () (female) if e(sample), stat(mean logyearn $xlist1) ///
    stat(count logyearn) nototals


* Separate earnings regressions for male and female doctors
qui regress logyearn $xlist1 if female==0, vce(robust)
estimates store MALE
qui regress logyearn $xlist1 if female==1, vce(robust)
estimates store FEMALE
estimates table MALE FEMALE, b(%11.5f) t(%11.2f) stats(N r2 F) 

* Tests of coefficient equality in separate regressions
qui regress logyearn $xlist1 if female==0
estimates store MALE
qui regress logyearn $xlist1 if female==1
estimates store FEMALE
suest MALE FEMALE
test _b[MALE_mean:yhrs]  = _b[FEMALE_mean:yhrs]

* Blinder-Oaxaca decomposition in logarithmic units
oaxaca logyearn $xlist1, by(female) vce(robust)

// Not included as output but mentioned in text - repeat in transformed units
oaxaca, eform 

*********** 4.7: SHAPLEY DECOMPOSITION OF RELATIVE REGRESSOR IMPORTANCE

* Shapley-Owen decomposition of regressors' contributions to R-squared
qui rego logyearn yhrs expr exprsq pracsize fellow pgradoth     ///
    if female==0, vce(robust) 
rego logyearn yhrs \ expr exprsq \ pracsize \ fellow pgradoth \   ///
    if female==0, vce(robust)   

// Not included - repeat without grouping variables
rego logyearn yhrs expr exprsq pracsize fellow pgradoth if female==0, vce(robust) 

*********** 4.8:  DIFFERENCE-IN-DIFFERENCES ESTIMATOR

qui use mus204nswpsid, clear
* DID: Outcome, treatment, and other variables  
describe re78 re75 treat age educ nodegree black hisp marr u74 u75 re74 ///
    agesq educsq re74sq re75sq u74black u74hisp

* Descriptive statistics for treatment and control samples
table () (treat), stat(mean re78 re75 treat age educ nodegree black hisp ///
    marr u74 u75 re74 agesq educsq re74sq re75sq u74black u74hisp)      ///
    stat(count re75) nototals

   
// Not included - Descriptive statistics of whole sample
summarize re78 re75 treat age educ nodegree black hisp marr ///
    u74 u75 re74 agesq educsq re74sq re75sq u74black u74hisp

* Before–after comparison for the treated 
generate badiff = re78 - re75
mean badiff if treat==1
 
* Treatment-control comparison (difference in means)
regress re78 treat, vce(robust)

// Not included - alternative way to do same test
ttest re78, by(treat) unequal

* DID - first stack into panel format
generate id = _n
generate earns1 = re75
generate earns2 = re78
qui reshape long earns, i(id) j(year)
generate dpost = 0
replace dpost = 1 if year==2
rename treat dtreat
qui generate dpostbydtreat = dpost*dtreat
summarize id year dtreat earns dpost dpostbydtreat age, sep(0)

* DID estimation - no controls
regress earns dpost dtreat dpostbydtreat, vce(cluster id)

* DID estimation - with time-invariant controls
regress earns dpost dtreat dpostbydtreat age educ nodegree black hisp marr ///
    u74 u75 re74 agesq educsq re74sq re75sq u74black u74hisp, vce(cluster id) 

* xtdidregress command for DID - no controls
xtset id year
xtdidregress (earns) (dpostbydtreat), group(id) time(year)

*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph

*********** END
