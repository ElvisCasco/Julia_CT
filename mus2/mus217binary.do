* mus217binary.do  for Stata 17

capture log close

********** OVERVIEW OF mus217binary.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press
 
* Chapter 17. Binary outcome models.
* 17.4: EXAMPLE
* 17.5: GOODNESS OF FIT AND PREDICTION
* 17.6: MES
* 17.7: CLUSTERED DATA
* 17.8: ADDITIONAL MODELS
* 17.9: ENDOGENOUS REGRESSORS
* 17.10: GROUPED AND AGGREGATE DATA

* To run you need file
*   mus217hrs.dta
* in your directory

* community-contributed commands 
*   prvalue
*   prchange
*   locreg 
*   moremata (needed for locreg)
*   kdens    (needed for locreg)
* are used

********** SETUP

clear all
set scheme s1mono  /* Graphics scheme */
  
********** DATA DESCRIPTION

* mus217hrs.dta authors' extract from from the Health and Retirement Study (HRS).
* wave 5 (2002) and is restricted to Medicare beneficiaries. 
 
*********** 17.4 EXAMPLE (LOGIT, PROBIT, OLS AND GLM MODELS)

* Read in data, define globals, and summarize key variables
qui use mus217hrs
global xlist age hstatusg hhincome educyear married hisp
global extralist female white chronic adl sretire
summarize ins retire $xlist $extralist

* Logit regression
logit ins retire $xlist, vce(robust)

* Comparison of estimates for logit, probit and LPM models
qui logit ins retire $xlist
estimates store blogit
qui probit ins retire $xlist 
estimates store bprobit
qui regress ins retire $xlist 
estimates store bols
qui logit ins retire $xlist, vce(robust)
estimates store blogitr
qui probit ins retire $xlist, vce(robust)
estimates store bprobitr 
qui regress ins retire $xlist, vce(robust)
estimates store bolsr

* Table for comparing models 
estimates table blogit blogitr bprobit bprobitr bols bolsr, ///
    t stats(N ll) b(%7.3f) stfmt(%8.2f) eq(1)

* Wald test for no interactions
global intlist c.age#c.age c.age#i.hstatusg c.age#c.hhincome ///
    c.age#c.educyear c.age#i.married c.age#i.hisp
qui logit ins retire $xlist $intlist, vce(robust) nolog
testparm $intlist

* LR test for no interactions
qui logit ins retire $xlist $intlist
estimates store B 
qui logit ins retire $xlist
lrtest B 

// Output not included
* Logit command results duplicated using the glm command
glm ins retire $xlist, link(logit) family(binomial) vce(robust)

* NLS estimation of logit model using the glm command
glm ins retire $xlist, link(logit) vce(robust) nolog

*********** 17.5: GOODNESS OF FIT AND PREDICTION

* Hosmer–Lemeshow goodness-of-fit test with 4 groups
qui logit ins retire $xlist
estat gof, group(4) table  // Hosmer–Lemeshow goodness-of-fit test

* Hosmer–Lemeshow goodness-of-fit test with 10 groups
estat gof, group(10)  // Hosmer–Lemeshow goodness-of-fit test

* Comparing fitted probability and dichotomous outcome
qui logit ins retire $xlist
estat classification

* Calculate and summarize fitted probabilities for models with a single regressor
qui logit ins hhincome
predict plogit, pr
qui probit ins hhincome  
predict pprobit, pr
qui regress ins hhincome
predict pols, xb
summarize ins plogit pprobit pols

* Plot of predicted probabilities for models with a single regressor
sort hhincome
graph twoway (scatter ins hhincome, msize(vsmall) jitter(3))  ///
    (line plogit hhincome, clstyle(p1))                       ///
    (line pprobit hhincome, clstyle(p2))                      ///
    (line pols hhincome, clstyle(p3)),                        ///
    plotregion(style(none))                                   ///
    title("Predicted probabilities across models")                       ///
    xtitle("Household income (hhincome)", size(medlarge)) xscale(titlegap(*5))   /// 
    ytitle("Predicted probability", size(medlarge)) yscale(titlegap(*5)) ///
    legend(pos(1) ring(0) col(1)) legend(size(small))                    ///
    legend(label(1 "Actual data (jittered)") label(2 "Logit")            ///
    label(3 "Probit") label(4 "OLS")) saving(graph1.gph, replace)        

* ROC curve following logit estimation
qui logit ins hhincome
lroc, lwidth(thick) title("Receiver operator characteristics curve") ///
    msize(tiny) saving(graph2.gph, replace) 
 
graph combine graph1.gph graph2.gph, iscale(1.2) ysize(2.5) xsize(6.0)

* Fitted probabilities for selected baseline using margins
qui logit ins retire $xlist, vce(robust)
margins, at(age=65 retire=0 hstatusg=1 hhincome=50 educyear=17 married=1 hisp=0) ///
    noatlegend

* Fitted probabilities for selected baseline using prvalue
qui logit ins retire $xlist, vce(robust)
prvalue, x(age=65 retire=0 hstatusg=1 hhincome=50 educyear=17 married=1 hisp=0)

*********** 17.6: MARGINAL EFFECTS

* AME after logit
qui logit ins i.retire age i.hstatusg hhincome educyear i.married i.hisp, ///
    vce(robust)
margins, dydx(*) noatlegend        // (AME)

* MEM after logit
qui logit ins i.retire age i.hstatusg hhincome educyear i.married i.hisp, ///
    vce(robust)
margins, dydx(*) atmeans noatlegend  // (MEM)

* MER after logit
qui logit ins i.retire age i.hstatusg hhincome educyear i.married i.hisp, ///
    vce(robust)
margins, dydx(*) at (retire=1 age=75 hstatusg=1 hhincome=35 educyear=12   ///
    married=1 hisp=1) noatlegend   // (MER)

* Computing change in probability after logit
qui logit ins retire $xlist, vce(robust)
prchange hhincome

*********** 17.7 CLUSTERED DATA

/*  Not included but examples of commands for clustered data
xtset educyear
xtgee ins retire, family(binomial) link(logit) corr(exchangeable) vce(robust)
xtlogit ins retire, pa corr(exchangeable)
xtlogit ins retire $xlist, re
*/ 

*********** 17.8 ADDITIONAL MODELS

* Heteroskedastic probit model
hetprobit ins retire $xlist, het(age chronic) nolog vce(robust)  

* Stukel score or Lagrange multiplier test for asymmetric h-family logit
qui logit ins retire $xlist
predict xbhat, xb
generate xbhatsq = xbhat^2
qui logit ins retire $xlist xbhatsq, vce(robust)
test xbhatsq

* Local logit regression
locreg ins, dummy(retire hstatusg married hisp)                ///
    continuous(age hhincome educyear) logit bandwidth(0.5 0.8) ///
    lambda(0.8 1.0) generate(ploclog, replace)

* Compare local logit and logit ML predicted probabilities
qui logit ins retire $xlist
qui predict plogitml   
qui scatter ploclog plogitml, xtitle("Logit ML predicted probability") ///
    ytitle("Local logit predicted probability") msize(tiny)            ///
    scale(1.2) saving(graph1.gph,replace)


// Not included - more detailes in predicted probabilities
summarize ins ploclog plogitml
correlate ins ploclog plogitml
summarize ploclog, d

*********** 17.9 ENDOGENOUS REGRESSORS  

* Endogenous probit using inconsistent probit MLE
generate linc = log(hhincome)
global xlist2 female age age2 educyear married hisp white chronic adl hstatusg
probit ins linc $xlist2, vce(robust) nolog

* Endogenous probit using ivprobit ML estimator
global ivlist2 retire sretire
ivprobit ins $xlist2 (linc = $ivlist2), vce(robust) nolog

// Output was not included 
* Endogenous probit using eprobit ML estimator (an erm command)
eprobit ins $xlist2, endogenous(linc =  $xlist2 $ivlist2) vce(robust) nolog

* Endogenous probit using ivprobit two-step estimator
ivprobit ins $xlist2 (linc = $ivlist2), twostep first

* Endogenous probit using ivregress to get 2SLS estimator
ivregress 2sls ins $xlist2 (linc = $ivlist2), vce(robust) noheader
estat overid

// Output was not included 
* Endogenous probit using nonlinear IV or GMM estimation
gmm (ins - normal({xb:linc $xlist2 _cons})), instruments($xlist2 $ivlist2) 
margins, dydx(*) expression(normal(predict(xb)))
estat overid

*********** 17.10 GROUPED AND AGGREGATE DATA 

* Changed July 18,2018

* Create grouped data 
qui use mus217hrs, clear
bysort age: egen num_in_cell = count(ins)
bysort age: gen num_ins_in_cell = num_in_cell*ins
sort age
collapse n=num_in_cell r=num_ins_in_cell av_ins=ins ///
    av_educyear=educyear av_hstatusg=hstatusg, by(age)
drop if n < 5
summarize

* Grouped data: Binomial logit model for number insured (= n times y)
glm r av_educyear av_hstatusg, family(binomial n) link(logit) noheader

* Grouped data: NLS with actual y as dependent variable
nl (av_ins = 1 / (1 + exp(-{xb: av_educyear av_hstatusg}+{b0}))), vce(robust) nolog   

* Grouped data: "logit" with actual y as dependent variable
fracreg logit av_ins av_educyear av_hstatusg, vce(robust) noheader nolog

* Grouped data: least squares with transformation of y as dependent variable
generate logins = log(av_ins/(1-av_ins))
regress av_ins av_educyear av_hstatusg, vce(robust) noheader

*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph

********** END
