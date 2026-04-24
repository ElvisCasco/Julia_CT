* mus225treatext.do  for Stata 17

cap log close

********** OVERVIEW OF mus2285reatext.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press  

********** SETUP ******

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

* Chapter 25:  ENDOGENOUS TREATMENT EFFECTS
* 25.3: ERM COMMANDS FOR ET
* 25.4: ET COMMANDS FOR BINARY TREATMENT
* 25.5: THE LATE ESTIMATOR FOR HETEROGENEOUS EFFECTS
* 25.6: DID AND SYNTHETIC CONTROL
* 25.7: RD DESIGN
* 25.8: CONDITIONAL QUANTILE REGRESSION WITH ENDOGENOUS REGRESSORS
* 25.9: UNCONDITIONAL QUANTILES

* To run this program you need files
*    mus224mcbs.dta
*    mus207mepspresdrugs.dta
*    mus210mepsdocvisyoung.dta
*    mus224ohiesmallrecode.dta 
*    mus225smoking.dta
*    mus225rdsenate.dta
*    mus203mepsmedexp.dta

* community-contributed commands
*    synth
*    synth_runner (includes effects_graphs and single_treatment_graphs)
*    rdrobust (includes rdplot and rdbwselect)
*    ivqte    (needs moremata and kdens)
*    moremata (needed for ivqte)
*    kdens    (needed for ivqte)
*    rifreg   (at https://sites.google.com/view/nicole-m-fortin/data-and-programs)
*    cdeco    (net install counterfactual, from("https://sites.google.com/site/mellyblaise/")
* are used

********** DATA DESCRIPTION **********

* mus225smoking.dta
* From Galiani and Quistorff (2017), "The synthrunnerpackage: Utilities to 
* automate synthetic control stimation using synth, " Stata Journal, 834-849.
* based in turn on Abadie, Diamond and Hainmuller (2010), 
* "Synthetic control methods for comparative case studies: Estimating the 
* effect of California's tobacco control program, JASA, 493-505.

* mus225rdsenate.dta from Galiani and Quistorff (2017), "The synthrunnerpackage:
* Utilities to automate synthetic control estimation using synth," 
* Stata Journal, 834-849.

************ CHAPTER 25.3  ERM COMMANDS FOR ENDOGENOUS TREATMENT

* Read in data, drop zero drug spending, create ordered multivalued treatment 
qui use mus224mcbs
generate drugexp = aamttot
drop if drugexp == 0
qui generate ldrugexp = ln(drugexp) 
qui replace income_c = income_c/1000
qui tabulate coverage, generate(inslevel)
qui summarize inslevel*
qui generate clevel= 0
qui replace clevel=0 if inslevel3==1  // Base Medicare
qui replace clevel=1 if inslevel2==1  // MMC (Medicare managed plan)
qui replace clevel=2 if inslevel4==1  // Medigap (Medicare suppl ins)
qui replace clevel=3 if inslevel1==1  // ESI (Employer sponsored)

// Not included 
tabulate clevel

* eregress example: log drug expenditure by insurance type
regress ldrugexp ibn.clevel, noheader noconstant

* eregress: Summary stats for outcome, exogenous variables, and instruments
global xlist h_age h_male income_c genhelth // Exogenous in both
global zlist prem1 prem2 prem3 prem4 penet  // Instruments
summarize ldrugexp $xlist $zlist

* eregress: Endogenous multivalued treatment and regression-adjustment model
eregress ldrugexp $xlist, entreat(clevel = $xlist $zlist) vce(robust) nolog

* eregress: ATE for endogenous multivalued treatment with unconditional standard errors
estat teffects, ate

// Not included - ATET estimate
* eregress: ATET for endogenous multivalued treatment
estat teffects, atet 

* eregress: ATE following regression with default standard errors had smaller standard errors
qui eregress ldrugexp $xlist, entreat(clevel = $xlist $zlist) nolog
estat teffects, ate

* eregress: Endogenous multivalued treatment & simpler no interaction model
qui eregress ldrugexp $xlist, entreat(clevel = $xlist $zlist, nointeract) ///
    vce(robust)
estat teffects, ate 

eregress ldrugexp $xlist, entreat(clevel = $xlist $zlist, nointeract) vce(robust)
estat teffects

********* CHAPTER 25.4: ET COMMANDS FOR ENDOGENOUS BINARY TREATMENT

* etregress: ML for endogenous binary treatment (insurance)  
qui use mus207mepspresdrugs, clear
global x2list totchr age female blhisp linc
etregress ldrugexp $x2list, treat(hi_empunion = $x2list ssiratio) ///
   vce(robust) first nolog

* etregress: ATE of insurance status following ML estimation
margins, predict(cte)

// Not included - gives same result
margins r.hi_empunion, contrast(nowald)

* etregress: Control function estimator for endog binary treatment (insurance)
qui etregress ldrugexp $x2list, treat(hi_empunion = $x2list ssiratio) ///
    cfunction vce(robust)
margins, predict(cte)

* etregress: More flexible potential-outcomes model and ATET
global x3list c.totchr c.age i.female i.blhisp c.linc  
qui etregress ldrugexp i.hi_empunion#($x3list),     ///
    treat(hi_empunion = $x3list ssiratio) vce(robust) nolog poutcomes 
margins, predict(cte) subpop(if hi_empunion==1)

* etpoisson sample: Subsample of young MEPS sample
qui use mus210mepsdocvisyoung, clear
set seed 10101
qui sample 2000, count
keep if docvis < 50
regress docvis private, vce(robust) noheader

* poisson: AME = ATE of exogenous binary treatment using margins and finite diffs
qui poisson docvis i.private i.chronic i.female age, vce(robust)
margins, dydx(i.private)

* etpoisson: ML for endogenous binary treatment
etpoisson docvis i.chronic i.female age, ///
    treat(private= i.chronic i.female age income) vce(robust) nolog

* etpoisson: ATE following ML for endogenous binary treatment
margins, predict(cte) vce(unconditional)

* etpoisson: ATET following ML for endogenous binary treatment
margins, predict(cte) subpop(if private==1) vce(unconditional)

* etpoisson: Manual computation of ATE using predict command
qui etpoisson docvis i.chronic i.female age, ///
    treat(private= i.chronic i.female age income) vce(robust) nolog
preserve
qui replace private=1 if private==0
qui predict docvis1 if private==1, pomean
qui replace private=0 if private==1
qui predict docvis0 if private==0, pomean
generate ate = docvis1 - docvis0
summarize ate docvis1 docvis0
restore

// Not included - gives same estimate
margins r.private, contrast(nowald)

********* CHAPTER 25.5: THE LATE ESTIMATOR

* LATE: Treatment D is Medicaid and instrument z is lottery win/lose
qui use mus224ohiesmallrecode, clear
label variable medicaid "Enrolled in Medicaid"
tabulate medicaid lottery, cell nokey

* LATE: Treatment is Medicaid and instrument is lottery win/lose
global y oop                       // Outcome variable for this chapter
global xlist dhhsize2 dhhsize3 dlotdraw* dsurvdraw*  // x variables for treat
qui regress $y lottery $xlist, vce(cluster household_id)
estimates store intent
qui regress $y medicaid $xlist, vce(cluster household_id)
estimates store ols
qui regress medicaid lottery $xlist, vce(cluster household_id)
estimates store first
qui ivregress 2sls $y (medicaid = lottery) $xlist, vce(cluster household_id)
estimates store iv
estimates table intent ols first iv, keep(medicaid lottery) ///
    b(%10.3f) se stat(N r2)

// Not included ut results mentioned in text
* LATE: Without controls
global xlist
qui regress $y lottery $xlist, vce(cluster household_id)
estimates store intent
qui regress $y medicaid $xlist, vce(cluster household_id)
estimates store ols
qui regress medicaid lottery $xlist, vce(cluster household_id)
estimates store first
qui ivregress 2sls $y (medicaid = lottery) $xlist, vce(cluster household_id)
estimates store iv
estimates table intent ols first iv, keep(medicaid lottery) ///
    b(%10.3f) se stat(N r2)

******** 25.6 DIFFERENCE IN DIFFERENCES AND SYNTHETIC CONTROL

* Synthetic control: Smoking dataset and summary
qui use mus225smoking, clear
summarize, sep(7)

* Synthetic control: synth command
tsset state year
qui synth cigsale beer(1984(1)1988) lnincome(1972(1)1988) retprice ///
    age15to24 cigsale(1988) cigsale(1980) cigsale(1975),           ///
	trunit(3) trperiod(1989) figure
	
* Synthetic control: synth_runner command			
synth_runner cigsale beer(1984(1)1988) lnincome(1972(1)1988) retprice ///
    age15to24 cigsale(1988) cigsale(1980) cigsale(1975),              ///
	trunit(3) trperiod(1989) gen_vars

* Synthetic control: synth_runner postestimation command effect_graphs	
effect_graphs, trlinediff(-1) tc_options(ti("Treated & control outcomes")) ///
    effect_options(title("Difference between treated & control"))
graph combine tc effect, xcommon ysize(2.5) xsize(6) scale(1.4)


* Synthetic control: synth_runner postestimation command single_treatment_graphs	
single_treatment_graphs, trlinediff(-1) effects_ylabels(-50(10)50)     ///
    do_color(gs11) effects_ymax(50) effects_ymin(-50)                  ///
    raw_options(title("Treated & donor outcomes"))                 ///
    effects_options(title("Difference between treated & donor outcomes")) 
graph combine raw effects, xcommon ysize(2.5) xsize(6) scale(1.3)


******** CHAPTER 25.7: REGRESSION DISCONTINUITY DESIGN

* SRD DGP: Quadratic in running variable x and TE 80 at x = 0
clear all
qui set obs 500
set seed 10101
generate x = 0.25*(_n - 250) + rnormal(0,5)  // The running variable
generate xsq = x^2
generate D = x > 0                           // The sharp cutoff  
generate y = -10 + 80*D + 2*x - 0.025*xsq + rnormal(0,60)  // The outcome
summarize y x D

* SRD design example: Scatterplot with separate global quadratic fits
twoway (scatter y x, xline(0) yline(-10, lpat(dash)) yline(70, lpat(dash))   ///
    msize(vsmall) xtitle("Running variable x") ytitle("Outcome") leg(off))   ///
    (qfitci y x if D==1, lcolor(black)) (qfitci y x if D==0, lcolor(black)), ///
    title("RD: Scatterplot and global quadratic fits") saving(graph1, replace)	

* SRD design example: ATE estimated by global quadratic
generate Dx = D*x
generate Dxsq = D*xsq
regress y D x xsq Dx Dxsq, vce(robust)

* SRD design example: Scatterplot with separate local linear fits
twoway (scatter y x, xline(0) yline(-10, lpat(dash)) yline(70, lpat(dash)) ///
    msize(vsmall) xtitle("Running variable x") ytitle("Outcome") leg(off)) ///
    (lpolyci y x if D==1, kernel(triangle) bw(20) deg(1) lcolor(black))    ///
    (lpolyci y x if D==0, kernel(triangle) bw(20) deg(1) lcolor(black)),   ///
    title("RD: Scatterplot and local linear fits") saving(graph2, replace)

graph combine graph1.gph graph2.gph, xcomm ycomm ysize(2.5) xsize(6) scale(1.4)

qui lpoly y x if D==0, kern(tri) bw(20) deg(1) gen(xminus yminus)
qui lpoly y x if D==1, kern(tri) bw(20) deg(1) gen(xplus yplus)
di "ATE = " yplus[1] " - " yminus[50] " = " yplus[1]-yminus[50]

* SRD design example: Binned data with separate local linear fits
sum x
scalar lowbw = (0 - r(min))/20
scalar highbw = (0 + r(max))/20
generate xbin = lowbw*floor(x/lowbw)+ lowbw/2
replace xbin = highbw*floor(x/highbw)+ highbw/2
bysort xbin: egen ybinmean = mean(y)
twoway (scatter ybinmean xbin, xline(0) xtitle("Bins of running variable x") ///
    msize(vsmall) ytitle("Mean of outcome in each bin of x") legend(off))    ///
    (lpoly y x if D==1, kern(tri) bw(20) deg(1) lcol(black) fcol(none))      ///
    (lpoly y x if D==0, kern(tri) bw(20) deg(1) lcol(black) fcol(none)),     ///
    title("RD: Scatterplot of binned data") saving(graph1, replace)

* SRD design example: Visual check no discontinuity in x at x = 0
bysort xbin: egen xcount = count(x)
twoway (scatter xcount xbin, xline(0) xtitle("Bins of running variable x") ///
    msize(vsmall) ytitle("Count of x in each bin of x") legend(off))       ///
    (lpoly xcount xbin if D==1, kern(triangle) bw(20) deg(1) lcol(black))  ///
    (lpoly xcount xbin if D==0, kern(triangle) bw(20) deg(1) lcol(black)), ///
    title("RD: Check discontinuity in running variable") saving(graph2, replace)
	
graph combine graph1.gph graph2.gph, xcommon ysize(2.5) xsize(6) scale(1.4)

*****

* SRD: Cattaneo et al. U.S. Senate elections data 
qui use mus225rdsenate, clear
summarize

* SRD: Scatterplot of the data 
twoway scatter vote margin, msize(tiny) xline(0)                      ///
	ytitle("Vote share at t+6") xtitle("Margin of victory at t")      ///
	title("Simple scatterplot of the data") saving(graph1.gph, replace)

* SRD: rdplot command for the same data using command defaults
rdplot vote margin, c(0)                                          ///
    graph_options(title("rdplot command using command defaults")  ///
	ytitle("Vote share at t+6") xtitle("Margin of victory at t")) 
graph save graph2.gph, replace

graph combine graph1.gph graph2.gph, xcommon ysize(2.5) xsize(6) scale(1.4)

* SRD: ATE estimate using parametric fourth-order polynomials
generate d = margin > 0
qui regress vote i.d##c.margin##c.margin##c.margin##c.margin, vce(robust)
di "Parametric ATE = " %6.3f _b[1.d] " het-robust st. error = " %6.3f _se[1.d]
qui regress vote i.d##c.margin##c.margin##c.margin##c.margin, vce(cluster state)
di "Parametric ATE = " %6.3f _b[1.d] " clu-robust st. error = " %6.3f _se[1.d]

* SRD: rdrobust command using command defaults aside from option all
rdrobust vote margin, c(0) all

* SRD: rdplot of the local linear estimate of ATE obtained by rdrobust
qui rdrobust vote margin, c(0)
qui rdplot vote margin if -e(h_l)<= margin & margin <= e(h_r),    ///
    binselect(esmv) kernel(triangular) h(`e(h_l)' `e(h_r)') p(1)  ///
    graph_options(title("RD plot of the default ATE estimate")    ///
	ytitle("Vote share at t+6") xtitle("Margin of victory at t")) 
graph save graph1.gph, replace

* SRD: rdplot of the local quadratic estimate with narrower bandwidth
qui rdplot vote margin if -10 <= margin & margin <= 10,           ///
    binselect(esmv) kernel(uniform) h(`-5' `5') p(2)              ///
    graph_options(title("RD plot with quadratic and narrower bw") ///
	ytitle("Vote share at t+6") xtitle("Margin of victory at t")) 
graph save graph2.gph, replace

graph combine graph1.gph graph2.gph, xcommon ysize(2.5) xsize(6) scale(1.3)

* SRD: rdrobust estimates with different estimation settings
qui rdrobust vote margin, c(0) all h(10 10) rho(0.5) p(1) kernel(triangular)
estimates store P1tri	
qui rdrobust vote margin, c(0) all h(10 10) rho(0.5) p(1) kernel(tri) vce(hc0)
estimates store sewhite
qui rdrobust vote margin, c(0) all h(10 10) rho(0.5) p(1) kernel(uniform)
estimates store P1unif
qui rdrobust vote margin, c(0) all h(10 10) rho(0.5) p(2) kernel(triangular)
estimates store P2tri
qui rdrobust vote margin, c(0) all h(15 15) rho(0.5) p(2) kernel(triangular)
estimates store wide
qui rdrobust vote margin, c(0) all h(5 5) rho(0.5) p(2) kernel(triangular)
estimates store narrow
estimates table P1tri sewhite P1unif P2tri wide narrow, b(%6.3f) se ///
   stfmt(%6.0f) stats(N_h_l N_h_r )

* SRD: rdbwselect command showing all available bandwidth selection methods
rdbwselect vote margin, c(0) all

* FRD DGP: 10% below cutoff are treated and 40% above cutoff are treated
gen dtreat = margin > 0
set seed 10101
qui replace dtreat = 1 if (runiform() < 0.1 & margin < 0)
qui replace dtreat = 0 if (runiform() < 0.4 & margin > 0)

* FRD: rdplot command for treatment indicator against margin
qui rdplot dtreat margin, c(0)                                          ///
    graph_options(title("Fuzzy rd: rdplot of treatment on running variable") ///
    ytitle("Treatment indicator at t") xtitle("Margin of victory at t")) 


* Fuzzy RD: rdrobust with fuzzy() option
rdrobust vote margin, c(0) fuzzy(dtreat) all

* SRD: Placebo test
qui rdplot population margin, c(0)                              ///
graph_options(title("Placebo rdplot using population")  ///
	ytitle("Population at t + 6") xtitle("Margin of victory at t")) 
graph save graph1.gph, replace
tsset state year
generate votelagged = 0
replace votelagged = l6.vote   // Six two-year periods ago = 12 years 
qui rdplot votelagged margin, c(0)                               /// 
graph_options(title("Placebo rdplot using lagged vote")  ///
	ytitle("Vote share at t - 6") xtitle("Margin of victory at t")) 
graph save graph2.gph, replace

graph combine graph1.gph graph2.gph, xcommon ysize(2.5) xsize(6) scale(1.3)

********** CHAPTER 25.8: CONDITIONAL QR WITH ENDOGENOUS TREATMENTS 

* Conditional QTE of exogenous suppins using ivqte (same as qreg) 
qui use mus203mepsmedexp, clear
drop if ltotexp == . 
ivqte ltotexp totchr age female white (suppins), quantile(0.5) variance

// Not included
correlate suppins marry

* Conditional QTE of endogenous suppins using ivqte with aai option 
ivqte ltotexp (suppins=marry), dummy(white female) continuous(totchr age) ///
   quantile(0.5) var trim(0.01) generate_p(predprob) aai  

// Not included
summarize predprob

// Also get ivqte aai coefficient at 0.25 and 0.75
ivqte ltotexp (suppins=marry), dummy(white female) continuous(totchr age) ///
   quantile(0.25) variance trim(0.01) aai  
ivqte ltotexp (suppins=marry), dummy(white female) continuous(totchr age) ///
   quantile(0.75) variance trim(0.01) aai  

// Not included is an ivqreg command for Chernozhukov and Hansen (2018)
// It takes a long time here and may be problematic
// ivqreg ltotexp totchr age female white (suppins=marry), quantile(0.5)

// Not included uses ivqreg2 for location-scale model and takes time
// First compare location-scale to qreg in exogenous case
// ivqreg2 ltotexp suppins totchr age female white, quantile(0.5)
//   qreg ltotexp suppins totchr age female white, quantile(0.5)
// Then do location-scale in endogenous case
// ivqreg2 ltotexp suppins totchr age female white,  ///
//    inst(marry totchr age female white) quantile(0.5)

********** CHAPTER 25.12: UNCONDITIONAL QUANTILES

* Unconditional QTE of suppins on ltotexp using community-contributed command ivqte
ivqte ltotexp (suppins), dummy(white female) continuous(totchr age) ///
   quantile(0.25(0.25)0.75) variance trim(0.1)

* Corresponding conditional QTE of suppins on ltotexp using qreg command
forvalues i = 1/3 {
    local j = `i'/4
    qui qreg ltotexp suppins totchr age white female, quant(`j') vce(robust)
    di "q = " `j' "   b = " _b[suppins]  "   se = " _se[suppins]
}

* Unconditional QTE of suppins at median using community-contributed command rifreg
rifreg ltotexp suppins totchr age female white, quantile(.25) ///
    generate(yval kdensval) 

// Not included - check the kernel density estimates
scatter kdensval yval

* Unconditional QTE of suppins at quartiles using community-contributed command rifreg
forvalues i = 1/3 {
    local j = `i'/4
    qui rifreg ltotexp suppins totchr age white female, quantile(`j')
    di "q = " `j' "   b = " _b[suppins] 
}

* Unconditional QTE of suppins on ltotexp using community-contributed command cdeco
set seed 10101
cdeco ltotexp totchr age white female, group(suppins) method(locsca) ///
    nreg(100) reps(400) quantile(0.25(0.25)0.75)

// Not included as takes a while - use method(qr) instead
// cdeco ltotexp totchr age white female, group(suppins) method(qr) ///
//    nreg(100) reps(400) quantile(0.25(0.25)0.75)

* Unconditional QTE of endogenous suppins using ivqte 
ivqte ltotexp (suppins=marry), dummy(white female) continuous(age totchr) ///
  variance trim(0.01) quantile(0.5)

// Not included - Also get unconditional ivqte estimates at 0.25 and 0.75
ivqte ltotexp (suppins=marry), dummy(white female) continuous(age totchr) ///
  variance trim(0.01) quantile(0.25) 
ivqte ltotexp (suppins=marry), dummy(white female) continuous(age totchr) ///
  variance trim(0.01) quantile(0.75)  
	
*** Erase files created by this program and not used elsewhere in the book

erase graph1.gph
erase graph2.gph

********** END
 
