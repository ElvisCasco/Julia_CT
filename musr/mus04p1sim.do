* MUS04P1SIM.DO  Oct 2009 for Stata version 11

cap log close

********** OVERVIEW OF mus04sim.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Chapter 4
* This program provides analysis for the second half of chapter 4
* Simulations for regression

* 4.2: PSEUDORANDOM NUMBER GENERATORS: INTRODUCTION
* 4.3; DISTRIBUTION OF THE SAMPLE MEAN
* 4.4: PSEUDORANDOM NUMBER GENERATORS: FURTHER DETAILS
* 4.5: COMPUTING INTEGRALS
* 4.6: SIMULATION FOR REGRESSION: INTRODUCTION

* To run, you need updated version of Stata 11 with the new random
* number generators runiform(), rnormal(), ...
* All data are generated.

********** SETUP **********

set more off
version 11
clear all
set memory 10m
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Generated data thoughout

********** 4.2: PSEUDORANDOM NUMBER GENERATORS: INTRODUCTION

* Single draw of a uniform number 
set seed 10101
scalar u = runiform()
display u

* 1000 draws of uniform numbers
quietly set obs 1000
set seed 10101
generate x = runiform()
list x in 1/5, clean
summarize x

// Output not included in book
histogram x, start(0) width(0.1)

* First three autocorrelations for the uniform draws
generate t = _n
tsset t
pwcorr x L.x L2.x L3.x, star(0.05)      

// Output not included in book
ac x      // Autocorrelations with 95% confidence band

* normal and uniform
clear
quietly set obs 1000
set seed 10101                           // set the seed
generate uniform = runiform()            // uniform(0,1)
generate stnormal = rnormal()            // N(0,1)
generate norm5and2 = rnormal(5,2)    
tabstat uniform stnormal norm5and2, stat(mean sd skew kurt min max) col(stat)

* t, chi-squared, and F with constant degrees of freedom
clear
quietly set obs 2000
set seed 10101
generate xt = rt(10)               // result xt ~ t(10)
generate xc = rchi2(10)            // result xc ~ chisquared(10) 
generate xfn = rchi2(10)/10        // result numerator of F(10,5)
generate xfd = rchi2(5)/5          // result denominator of F(10,5)
generate xf = xfn/xfd              // result xf ~ F(10,5) 
summarize xt xc xf

* Discrete rv's: binomial
set seed 10101
generate p1 = runiform()              // here p1~uniform(0,1) 
generate trials = ceil(10*runiform()) // here # trials varies btwn 1 & 10
generate xbin = rbinomial(trials,p1)  // draws from binomial(n,p1)  
summarize p1 trials xbin

* Discrete rv's: independent poisson and negbin draws
set seed 10101 
generate xb= 4 + 2*runiform()
generate xg = rgamma(1,1)           // draw from gamma;E(v)=1
generate xbh = xb*xg                // apply multiplicative heterogeneity
generate xp = rpoisson(5)           // result xp ~ Poisson(5)
generate xp1 = rpoisson(xb)         // result xp1 ~ Poisson(xb)
generate xp2 = rpoisson(xbh)        // result xp2 ~ NB(xb)
summarize xg xb xp xp1 xp2

* Example of histogram and kernel density plus graph combine
quietly twoway (histogram xc, width(1)) (kdensity xc, lwidth(thick)), ///
   title("Draws from chisquared(10)") 
quietly graph save mus04cdistr.gph, replace
quietly twoway (histogram xp, discrete) (kdensity xp, lwidth(thick) ///
   w(1)),  title("Draws from Poisson(mu) for 5<mu<6") 
quietly graph save mus04poissdistr.gph, replace
graph combine mus04cdistr.gph mus04poissdistr.gph, ///
   title("Random-number generation examples", margin(b=2) size(vlarge)) 
graph export mus04fig1.eps, replace


* Simulation using postfile
set seed 10101
postfile sim_mem xmean using simresults, replace
forvalues i = 1/10000 {
     drop _all
     quietly set obs 30
     tempvar x
     generate `x' = runiform()
     quietly  summarize `x'
     post sim_mem (r(mean))
}
postclose sim_mem

* See the results stored in simresults
use simresults, clear
summarize



********** 4.3 DISTRIBUTION OF THE SAMPLE MEAN

clear all
* Draw 1 sample of size 30 from uniform distribution 
quietly set obs 30
set seed 10101
generate x = runiform()

* Summarize x and produce a histogram
summarize x
quietly histogram x, width(0.1) xtitle("x from one sample")

quietly graph export mus04fig1clt.eps, replace

* Program to draw 1 sample of size 30 from uniform and return sample mean
program onesample, rclass
    drop _all
    quietly set obs 30
    generate x = runiform()
    summarize x
    return scalar meanforonesample = r(mean)
end

* Run program onesample once as a check
set seed 10101
onesample
return list

* Run program onesample 10,000 times to get 10,000 sample means 
simulate xbar = r(meanforonesample), seed(10101) reps(10000) nodots: onesample

* Summarize the 10,000 sample means and draw histogram
summarize xbar
quietly histogram xbar, normal xtitle("xbar from many samples")

quietly graph export mus04fig1clt2.eps, replace

********** 4.4: RANDOM NUMBER GENERATION: FURTHER DETAILS

clear all
* Inverse probability transformation example: standard normal
quietly set obs 2000
set seed 10101
generate xstn = invnormal(runiform())     

* Inverse probability transformation example: unit exponential
generate xue = -ln(1-runiform())

* Inverse probability transformation example: Bernoulli (p = 0.6)
generate xbernoulli = runiform() > 0.6   // Bernoulli(0.6)
summarize xstn xue xbernoulli
 
* Draws from truncated normal x ~ N(mu,sigma^2) in [a,b]
quietly set obs 2000
set seed 10101
scalar a = 0                             // lower truncation point
scalar b = 12                            // upper truncation point
scalar mu = 5                            // mean
scalar sigma = 4                         // standard deviation
generate u = runiform()
generate w=normal((a-mu)/sigma)+u*(normal((b-mu)/sigma)-normal((a-mu)/sigma))
generate xtrunc = mu + sigma*invnormal(w)
summarize xtrunc

// Note: The following may yield different output to the book

* Bivariate normal example: 
* means 10, 20; variances 4, 9; and correlation 0.5
clear
quietly set obs 1000
set seed 10101
matrix MU = (10,20)                   // MU is 2 x 1
scalar sig12 = 0.5*sqrt(4*9)
matrix SIGMA = (4, sig12 \ sig12, 9)  // SIGMA is 2 x 2 
drawnorm y1 y2, means(MU) cov(SIGMA)
summarize y1 y2
correlate y1 y2

* MCMC example: Gibbs for bivariate normal mu's=0 v's=1 corr=rho=0.9
set seed 10101    
clear all
set obs 1000
generate double y1 =.
generate double y2 =.
mata:
  s0 = 10000             // Burn-in for the Gibbs sampler (to be discarded)
  s1 = 1000              // Actual draws used from the Gibbs sampler
  y1 = J(s0+s1,1,0)      // Initialize y1 
  y2 = J(s0+s1,1,0)      // Initialize y2 
  rho = 0.90             // Correlation parameter
  for(i=2; i<=s0+s1; i++) {
      y1[i,1] = ((1-rho^2)^0.5)*(rnormal(1, 1, 0, 1)) + rho*y2[i-1,1]
      y2[i,1] = ((1-rho^2)^0.5)*(rnormal(1, 1, 0, 1)) + rho*y1[i,1]
  }
  y = y1,y2
  y = y[|(s0+1),1 \ (s0+s1),.|]  // Drop the burn-ins
  mean(y)                        // Means of y1, y2
  variance(y)                    // Variance matrix of y1, y2
  correlation(y)                 // Correlation matrix of y1, y2
end 

mata:
  y2 = y[|2,2 \ s1,2|]
  y2lag1 = y[|1,2 \ (s1-1),2|]
  y2andlag1 = y2,y2lag1
  correlation(y2andlag1,1)       // Correlation between y2 and y2 lag 1
end

************ 4.5: COMPUTiNG INTEGRALS

* Integral evaluation by Monte Carlo simulation with S=100 
clear all
quietly set obs 100
set seed 10101
generate double y = invnormal(runiform())
generate double gy = exp(-exp(y))
quietly summarize gy, meanonly
scalar Egy = r(mean)
display "After 100 draws the MC estimate of E[exp{-exp(x)}] is " Egy

* Program mcintegration to compute E{g(y)} numsims times
program mcintegration, rclass
   version 11 
   args numsims        // Call to program will include value for numsims
   drop _all
   quietly set obs `numsims'
   set seed 10101
   generate double y = rnormal(0) 
   generate double gy = exp(-exp(y)) 
   quietly summarize gy, meanonly
   scalar Egy = r(mean)
   display "#simulations: " %9.0g `numsims'  ///
     "  MC estimate of  E[exp{-exp(x)}] is " Egy
end

* Run program mcintegration S = 10, 100, ...., 100000 times
mcintegration 10 
mcintegration 100
mcintegration 1000
mcintegration 10000
mcintegration 100000


// * Following code was not included in the book
// * It compares use of simulate and postfile for a simple example   

* Program to generate one sample and return the sample median program 

* Simulate distribution of sample median
program simulateexample, rclass
version 11
drop _all
set obs 9
tempvar y
generate `y' = rnormal(0)
centile `y'
return scalar median = r(c_1)
end

* Run program once
set seed 10101
simulateexample
display "Median is " r(median)

* Run program 1000 times
set seed 10101
simulate mutilde = r(median), reps(1000) nodots: simulateexample 

* See the results of the program.  
summarize 

* Simulation using postfile  
set seed 10101
tempname simfile
postfile `simfile' mutilde using simresults, replace
forvalues i = 1/1000 {
     drop _all
     quietly set obs 9
     tempvar y
     generate `y' = rnormal(0)
     quietly  centile `y'
     post `simfile' (r(c_1))
     }
postclose `simfile'

* See the results stored in simresults
use simresults, clear
summarize

// End example not in book 

************ 4.6: SIMULATION FOR REGRESSION: INTRODUCTION

* Define global macros for sample size and number of simulations
global numobs 150             // sample size N
global numsims "1000"         // number of simulations 

* Program for finite-sample properties of OLS
program chi2data, rclass
    version 11  
    drop _all
    set obs $numobs
    generate double x = rchi2(1)   
    generate y = 1 + 2*x + rchi2(1)-1     // demeaned chi^2 error 
    regress y x
    return scalar b2 =_b[x]
    return scalar se2 = _se[x]
    return scalar t2 = (_b[x]-2)/_se[x]
    return scalar r2 = abs(return(t2))>invttail($numobs-2,.025)
    return scalar p2 = 2*ttail($numobs-2,abs(return(t2)))
end

* Show that test gives same result as doing test manually
set seed 10101
quietly chi2data
return list
quietly test x=2
return list

// Note: ideally the seed would have been reset to 10101 before the following

* Simulation for finite-sample properties of OLS
simulate b2f=r(b2) se2f=r(se2) t2f=r(t2) reject2f=r(r2) p2f=r(p2),  ///
reps($numsims) saving(chi2datares, replace) nolegend nodots: chi2data
summarize b2f se2f reject2f

* Summarize results
mean b2f se2f reject2f

// histogram t2f
// histogram p2

* t-statistic distribution
kdensity t2f,  n(1000) gen(t2_x t2_d) nograph
generate double t2_d2 = tden(148, t2_x)
graph twoway (line t2_d t2_x) (line t2_d2 t2_x)
graph export mus04ttest.eps, replace

// Following was mislabelled as fixed regressors rather than power

* Program for finite-sample properties of OLS: power
program chi2datab, rclass
    version 11  
    drop _all
    set obs $numobs
    generate double x = rchi2(1)   
    generate y = 1 + 2*x + rchi2(1)-1     // demeaned chi^2 error 
    regress y x
    return scalar b2  =_b[x]
    return scalar se2 =_se[x]
    test x=2.1
    return scalar r2 = (r(p)<.05)
end

* Power simulation for finite-sample properties of OLS
simulate b2f=r(b2) se2f=r(se2) reject2f=r(r2),  ///
reps($numsims) saving(chi2databres, replace) nolegend nodots: chi2datab
mean b2f se2f reject2f

* Inconsistency of OLS in errors-in-variables model (measurement error)
clear
quietly set obs 10000
set seed 10101
matrix mu = (0,0,0)
matrix sigmasq = (9,0,0\0,1,0\0,0,1)
drawnorm xstar u v, means(mu) cov(sigmasq)
generate y = 1*xstar + u   // DGP for y depends on xstar
generate x = xstar + v     // x is mismeasured xstar 
regress y x, noconstant

* Program for finite-sample properties of OLS: random regressors
clear
set seed 10101
program randomreg, rclass
    version 11 
    drop _all
    set obs $numobs
    generate x = rnormal(0)                  // random regressors 
    generate y = 10 + 5*x + 1*rnormal(0)     // normal error 
    regress y x
    return scalar b2 =_b[x]
    return scalar se2 = _se[x]
    return scalar t2 = (_b[x]-5)/_se[x]
    return scalar r2 = abs(return(t2))>invttail($numobs-2,.025)
    return scalar p2 = 2*ttail($numobs-2,abs(return(t2)))
end

* Simulation for finite-sample properties of OLS: random regressors
simulate b2r=r(b2) se2r=r(se2) t2r=r(t2) reject2r=r(r2) p2r=r(p2),  ///
     reps($numsims) nolegend nodots: randomreg
summarize
mean b2r se2r t2r reject2r p2r, noheader

* Postfile and post example: random regressors again
clear
set seed 10101
program simbypost
    version 11 
    tempname simfile
    postfile `simfile' b2 se2 t2 reject2 p2 using simresults, replace
    quietly {
       forvalues i = 1/$numsims { 
         drop _all
         set obs $numobs
         generate x = rnormal(0)                 // random regressors 
         generate y = 10 + 5*x + 1*rnormal(0)    // normal error 
         regress y x
         scalar b2 =_b[x]
         scalar se2 = _se[x]
         scalar t2 = (_b[x]-5)/_se[x]
         scalar reject2 = abs(t2) > invttail($numobs-2,.025)
         scalar p2 = 2*ttail($numobs-2,abs(t2))
         post `simfile' (b2) (se2) (t2) (reject2) (p2)
      }
    }
    postclose `simfile'
end
simbypost
use simresults, clear
summarize
mean b2 se2 t2 reject2 p2, noheader
 
* Endogenous regressor
clear
set seed 10101
program endogreg, rclass
    version 11 
    drop _all
    set obs $numobs
    generate u = rnormal(0) 
    generate x = 0.5*u + rnormal(0)         // endogenous regressors
    generate y = 10 + 2*x + u 
    regress y x
    return scalar b2 =_b[x]
    return scalar se2 = _se[x]
    return scalar t2 = (_b[x]-2)/_se[x]
    return scalar r2 = abs(return(t2))>invttail($numobs-2,.025)
    return scalar p2 = 2*ttail($numobs-2,abs(return(t2)))
end

simulate b2r=r(b2) se2r=r(se2) t2r=r(t2) reject2r=r(r2) p2r=r(p2),  ///
     reps($numsims) nolegend nodots: endogreg
mean b2r se2r reject2r 

********** CLOSE OUT
