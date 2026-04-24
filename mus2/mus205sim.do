* mus205sim.do  for Stata 17

cap log close

********** OVERVIEW OF mus205sim.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* 5.2: PSEUDORANDOM-NUMBER GENERATORS
* 5.3: DISTRIBUTION OF THE SAMPLE MEAN
* 5.4: PSEUDORANDOM-NUMBER GENERATION: FURTHER DETAILS
* 5.5: COMPUTING INTEGRALS 
* 5.6: SIMULATION FOR REGRESSION: INTRODUCTION

* No community-contributed commands are used

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Generated data thoughout

* Results will differ if random number generators from 
* versions of Stata prior to version 14 are used.

********** CHAPTER 5.2: PSEUDORANDOM NUMBER GENERATORS

* Single draw of a uniform number 
set seed 10101
scalar u = runiform()
display u

* 1,000 draws of uniform numbers
qui set obs 1000
set seed 10101
generate x = runiform()
list x in 1/5, clean
summarize x

// Not included as output in text but mentioned
histogram x, start(0) width(0.1)

* First three autocorrelations for the uniform draws
generate t = _n
tsset t
pwcorr x L.x L2.x L3.x, star(0.05)      

// Not included as output - autocorrelations with 95% confidence band
ac x      
* Normal and uniform and uniform draws
clear
qui set obs 1000
set seed 10101                           // Set the seed
generate uniform = runiform()            // uniform(0,1)
generate stnormal = rnormal()            // N(0,1)
generate norm5and2 = rnormal(5,2)    
tabstat uniform stnormal norm5and2, stat(mean sd skew kurt min max) col(stat)

* Student's t, chi-squared, and F draws with constant degrees of freedom
clear
qui set obs 2000
set seed 10101
generate xt = rt(10)               // Result xt ~ t(10)
generate xc = rchi2(10)            // Result xc ~ chisquared(10) 
generate xfn = rchi2(10)/10        // Result numerator of F(10,5)
generate xfd = rchi2(5)/5          // Result denominator of F(10,5)
generate xf = xfn/xfd              // Result xf ~ F(10,5) 
summarize xt xc xf

* Discrete random variables: Binomial draws with n and p varying over trials
set seed 10101
generate p1 = runiform()              // Here p1~uniform(0,1) 
generate trials = runiformint(1,10)   // Here # trials varies btwn 1 & 10
generate xbin = rbinomial(trials,p1)  // Draws from binomial(n,p1)  
summarize p1 trials xbin

* Discrete random variables: Poisson and negative binomial draws i.i.d. 
set seed 10101 
generate xb= 4 + 2*runiform()
generate xg = rgamma(1,1)           // Draw from gamma;E(v)=1
generate xbh = xb*xg                // Apply multiplicative heterogeneity
generate xp = rpoisson(5)           // Result xp ~ Poisson(5)
generate xp1 = rpoisson(xb)         // Result xp1 ~ Poisson(xb)
generate xp2 = rpoisson(xbh)        // Result xp2 ~ NB(xb)
summarize xg xb xp xp1 xp2

* Example of histogram and kernel density plus graph combine
qui twoway (histogram xc, width(1)) (kdensity xc, lwidth(thick)), ///
    legend(off) title("Draws from chisquared(10)") saving(graph1.gph, replace) 
qui twoway (histogram xp, discrete) (kdensity xp, lwidth(thick) w(1)),  ///
    legend(off) title("Draws from Poisson(mu=5)") saving(graph2.gph, replace)
graph combine graph1.gph graph2.gph, iscale(1.2) ysize(2.5) xsize(6.0) 

********** 5.3: DISTRIBUTION OF THE SAMPLE MEAN

* Draw 1 sample of size 30 from uniform distribution 
clear all
qui set obs 30
set seed 10101
generate x = runiform()

* Summarize x and produce a histogram
summarize x
histogram x, width(0.1) xtitle("x's from one sample") scale(1.2)

* Program to draw 1 sample of size 30 from uniform and return sample mean
program onesample, rclass
    drop _all
    qui set obs 30
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
histogram xbar, normal xtitle("Sample mean xbar from many samples") scale(1.2) 


* Simulation using postfile
set seed 10101
postfile sim_mem xmean using simresults, replace
forvalues i = 1/10000 {
    drop _all
    qui set obs 30
    tempvar x
    generate `x' = runiform()
    qui  summarize `x'
    post sim_mem (r(mean))
}
postclose sim_mem

* See the results stored in simresults
use simresults, clear
summarize

********** CHAPTER 5.4: RANDOM NUMBER GENERATION: FURTHER DETAILS

* Inverse-probability transformation example: Standard normal
clear all
qui set obs 2000
set seed 10101
generate xstn = invnormal(runiform())     

* Inverse-probability transformation example: Unit exponential
generate xue = -ln(1-runiform())

* Inverse-probability transformation example: Bernoulli (p = 0.6)
generate xbernoulli = runiform() > 0.6   // Bernoulli(0.6)
summarize xstn xue xbernoulli
 
* Draws from truncated normal x ~ N(mu,sigma^2) in [a,b]
qui set obs 2000
set seed 10101
scalar a = 0                             // Lower truncation point
scalar b = 12                            // Upper truncation point
scalar mu = 5                            // Mean
scalar sigma = 4                         // Standard deviation
generate u = runiform()
generate w=normal((a-mu)/sigma)+u*(normal((b-mu)/sigma)-normal((a-mu)/sigma))
generate xtrunc = mu + sigma*invnormal(w)
summarize xtrunc

* Bivariate normal example: Means 10, 20; variances 4, 9; and correlation 0.5
clear
qui set obs 1000
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

************ CHAPTER 5.5: COMPUTiNG INTEGRALS

* Integral evaluation by Monte Carlo simulation with S=100 
clear all
qui set obs 100
set seed 10101
generate double y = runiform()
generate double gy = exp(-exp(y))
qui summarize gy, meanonly
scalar Egy = r(mean)
display "After 100 draws the MC estimate of E[exp{-exp(x)}] is " Egy

* Program mcintegration to compute E{g(y)} numsims times
program mcintegration, rclass
    version 17 
    args numsims        // Call to program will include value for numsims
    drop _all
    qui set obs `numsims'
    set seed 10101
    generate double y = rnormal(0) 
    generate double gy = exp(-exp(y)) 
    qui summarize gy
    scalar Egy = r(mean)
    scalar seEgy = r(sd)/sqrt(`numsims')
    display "#sims:" %7.0f `numsims'  "  MC estimate is " Egy  ///
            "  Standard error is " seEgy
end

* Run program mcintegration S = 10, 100, ...., 100000 times
mcintegration 10 
mcintegration 100
mcintegration 1000
mcintegration 10000
mcintegration 100000

* Generate 2-dimensional Halton sequences of length 200
clear
qui set obs 200
mata
    h = halton(200,2)
    st_addvar("float", ("halton1", "halton2"))
    st_store(., ("halton1", "halton2"),h[.,.])
end
list halton1 halton2 in 1/8, clean

* Summary statistics and serial correlation for Halton sequences
summarize halton1 halton2
qui generate laghalt1 = halton1[_n-1]
correlate halton1 laghalt1

// Not included but mentioned: How to generate Hammersley with start at 1
//   h = halton(200,2,1,1)  // .,.,1,1 is start=1, hammersley

* Compare deterministic Halton draws to pseudo-random uniform draws 
set seed 10101
gen uniform1 = runiform()
gen uniform2 = runiform()
qui graph twoway (scatter halton1 halton2), title("Halton sequences") ///
    saving(graph1.gph, replace)
qui graph twoway (scatter uniform1 uniform2), title("Uniform draws") ///
    saving(graph2.gph, replace)
qui graph combine graph1.gph graph2.gph, ysize(2.5) xsize(6) iscale(1.5)

// Following code was not included in the book
// It compares use of simulate and postfile for a simple example  

* Program to generate one sample and return the sample median program 

* Simulate distribution of sample median
program simulateexample, rclass
    version 17
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
    qui set obs 9
    tempvar y
    generate `y' = rnormal(0)
    qui  centile `y'
    post `simfile' (r(c_1))
    }
postclose `simfile'

* See the results stored in simresults
use simresults, clear
summarize

// End example not in book 

********** CHAPTER 5.6: SIMULATION FOR REGRESSION: INTRODUCTION

* Define global macros for sample size and number of simulations
global numobs 150             // Sample size N
global numsims "1000"         // Number of simulations 

* Program for finite-sample properties of OLS
program chi2data, rclass
    version 17  
    drop _all
    set obs $numobs
    generate double x = rchi2(1)   
    generate y = 1 + 2*x + rchi2(1)-1     // Demeaned chi^2 error 
    regress y x
    return scalar b2 =_b[x]
    return scalar se2 = _se[x]
    return scalar t2 = (_b[x]-2)/_se[x]
    return scalar r2 = abs(return(t2))>invttail($numobs-2,.025)
    return scalar p2 = 2*ttail($numobs-2,abs(return(t2)))
end

* An F test gives same p-value as the manual t test in program chi2data
set seed 10101
qui chi2data
return list
qui test x=2
return list

* Simulation for finite-sample properties of OLS
simulate b2f=r(b2) se2f=r(se2) t2f=r(t2) reject2f=r(r2) p2f=r(p2),  ///
    seed(10101) reps($numsims) nolegend nodots: chi2data
summarize b2f se2f reject2f

* Report results for simulation averages
mean b2f se2f reject2f

* Plot the density of the t statistics and compare with theoretical t(148)
kdensity t2f, student(148) legend(off) scale(1.2) ///
    title("Density of the {it:t} statistics versus {it:t}(148)")

* Get 2.5 and 97.5 percentils for the t statistic
centile t2f, centile(2.5 97.5)
display invttail(148, 0.025)

* Histogram of p2f
histogram p2f

/* Comment out as takes time
* Increase to 10,000 simulations
simulate b2f=r(b2) se2f=r(se2) t2f=r(t2) reject2f=r(r2) p2f=r(p2),  ///
    seed(10101) reps(10000) nolegend nodots: chi2data
mean reject2f
*/

* Program for finite-sample properties of OLS: power
program chi2datab, rclass
    version 17
    drop _all
    set obs $numobs
    generate double x = rchi2(1)
    generate y = 1 + 2.1*x + rchi2(1)-1     // Demeaned chi^2 error
    regress y x
    return scalar b2  =_b[x]
    return scalar se2 =_se[x]
    test x=2
    return scalar r2 = (r(p)<.05)
end

* Power simulation for finite-sample properties of OLS
simulate b2f=r(b2) se2f=r(se2) reject2f=r(r2),  ///
    seed(10101) reps($numsims) nolegend nodots: chi2datab
mean b2f se2f reject2f

* Consistency of OLS in preceding simulation setup
clear
qui set obs 10000
set seed 10101
generate double x = rchi2(1)
generate y = 1 + 2*x + rchi2(1)-1     // Demeaned chi^2 error
regress y x, noheader

* Inconsistency of OLS in errors-in-variables model (measurement error)
clear
qui set obs 10000
set seed 10101
matrix mu = (0,0,0)
matrix sigmasq = (9,0,0\0,1,0\0,0,1)
drawnorm xstar u v, means(mu) cov(sigmasq)
generate y = 1*xstar + u   // DGP for y depends on xstar
generate x = xstar + v     // x is mismeasured xstar 
regress y x, noconstant noheader

* Program for OLS with endogenous regressor
clear
program endogreg, rclass
    version 17 
    drop _all
    set obs $numobs
    generate u = rnormal(0) 
    generate z = rnormal()
    generate v = rnormal(0)
    generate x = z + u + v   // Endogenous regressor
    generate y = 10 + 2*x + u 
    regress y x
    return scalar b2 =_b[x]
    return scalar se2 = _se[x]
    test x=0
    return scalar r2 = (r(p)<.05)
end

* Simulation for OLS with endogenous regressor
simulate b2r=r(b2) se2r=r(se2) reject2r=r(r2),              ///
     seed(10101) reps($numsims) nolegend nodots: endogreg
mean b2r se2r reject2r 

*** Erase files created by this program and not used elsewhere in the book

erase simresults.dta
erase graph1.gph
erase graph2.gph

********** END
