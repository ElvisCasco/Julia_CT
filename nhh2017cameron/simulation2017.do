* simulation2017.do  based on crete2016_simulation August 2017
* By A. Colin Cameron 

log using simulation2017.txt, text replace

********** OVERVIEW OF simulation2017.do **********

* RANDOM NUMBER DRAWS
* COMPUTING INTEGRALS
* MONTE CARLO EXPERIMENT OLS with Error chisquare(1) - 1

* To run you need files
*   no files    
* in your directory
* No Stata user-written commands are used

********** SETUP **********

clear all
set more off
version 14
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Only simulated data used

****  RANDOM NUMBER DRAWS 

*** PLot of fifty uniform draws - do they look independent
set obs 50
set seed 10101
generate xuniform = runiform()
generate drawnumber = _n
summarize 
scatter xuniform drawnumber, title("Uniform draw vs.draw number")
graph export simulation01fig.wmf, replace

***  Now do 1,000 draws and compare to uniform
clear
set obs 10000
set seed 10101
generate x = runiform()
list x in 1/5, clean
summarize x
display "Theoretical mean = 0.5 and standard deviaion = " 1/sqrt(12)
histogram x, start(0) width(0.1)

*** Autocorrelations for the uniform draws should be zero
generate t = _n
tsset t
* line x t if t <= 100
pwcorr x L.x L2.x L3.x, star(0.05)
ac x

*** Inverse Transformation Method 

* N(5, 2^2) draws using old inverse transformation method and new rnormal()
clear
set obs 1000
set seed 10101 
generate x1 = 5 + 2*invnorm(runiform())
generate x2 = rnormal(5,2)    
graph twoway (kdensity x1) (kdensity x2)
summarize x1 x2
tabstat x1 x2, stat(mean sd skew kurt min max) col(stat)

* Unit exponential  u = F(x) = 1 - exp(-x) so x = F^-1(u_ = -ln(1-u)
clear
set obs 10000
set seed 10101 
generate u = runiform()
generate x = -ln(1-u)
summarize u x
sort x
line u x, yline(0.64) xline(1.0216) ytitle("Cdf F(x)") plotreg(style(none)) ///
    title("Inverse Transformation Method") xtitle("Random variable x")
graph export simulation03fig.wmf, replace	

*** Gibbs sampler example

* MCMC example: Gibbs for bivariate normal mu's=0 v's=1 corr=rho=0.9
clear all
set seed 10101    
set obs 1000
generate double y1 =.
generate double y2 =.

clear all
set seed 10101
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
  // Skip view in mata:  mean(y), variance(y), correlation(y)
  stata("quietly set obs 1000")          // This requires s1 = 1000
  st_addvar("float", ("y1", "y2"))
  st_store(., ("y1", "y2"), y)
end 

summarize
correlate y1 y2
gen s = _n
tsset s
corrgram y1, lag(5)

****** COMPUTING INTEGRALS

* Find the mean of standard normal using Riemann sum midpoint method
* Evaluate 

clear
set obs 100
set seed 01101
generate x = rnormal(0,1)
generate f1x = x
generate f2x = exp(-exp(x)) 
summarize f1x f2x
mean f1x f2x

clear
set obs 10000
set seed 01101
generate x = rnormal(0,1)
* We want E[X] for N(0,1) = integral xf(x) where f is normal density
* We want E[exp(-exp(X))] for N(0,1) = integral exp(-exp(x)f(x) where f is normal density
generate f1x = x
generate f2x = exp(-exp(x)) 
summarize f1x f2x
mean f1x f2x

clear
* x is the evaluation point - ranges from -2.95 to 3.95
set obs 81
generate x = 3.95 - 0.1*_n  
* We want E[X] for N(0,1) = integral xf(x) where f is normal density
generate f1xtimesbase = 0.1*x*normalden(x)
summarize f1xtimesbase
display "E[X] computed by midpoint method = " r(sum)
* We want E[exp(-exp(X))] for N(0,1) = integral exp(-exp(x)f(x) where f is normal density
generate f2xtimesbase = 0.1*exp(-exp(x))*normalden(x) 
summarize f2xtimesbase
display "E[exp(-exp(X))] computed by midpoint method = " r(sum)

***** MONTE CARLO EXPERIMENT OLS with Error chisquare(1) - 1

* Small sample - estimates will differ from d.g.p. values
clear
set seed 101
quietly set obs 30
generate double x = rchi2(1)   
generate y = 1 + 2*x + rchi2(1)-1   
regress y x, noheader

* Consistency - sample size is set large at e.g. 100000
clear
set seed 10101
quietly set obs 100000
generate double x = rchi2(1)   
generate y = 1 + 2*x + rchi2(1)-1   
regress y x, noheader

* Set up for Monte Carlo experiment
clear
clear programs
global numobs 150 
global numsims "1000"

* Postfile example for simulation
clear
set seed 54321
postfile simalternative b se t r p using simresults, replace
forvalues i = 1/$numsims {
    drop _all
    quietly set obs $numobs
    quietly generate double x = rchi2(1)   
    quietly generate y = 1 + 2*x + rchi2(1)-1     // demeaned chi^2 error 
    quietly regress y x
    scalar b2 =_b[x]
    scalar se2 = _se[x]
    scalar t2 = (_b[x]-2)/_se[x]
    scalar r2 = abs(t2)>invttail($numobs-2,.025)
    scalar p2 = 2*ttail($numobs-2,abs(t2))
    post simalternative (b2) (se2) (t2) (r2) (p2)
}
postclose simalternative

* Analyze the simulation results
use simresults, clear
summarize
mean b se t r p 

/* For some reason problem here, so comment out 
* Simulation interval for test size when test at 0.05 and 1000 simulations
cii proportions 1000 0.05
*/

* t-statistic distribution
summarize t 
kdensity t,  n(1000) gen(t_x t_d) nograph
generate double t_d2 = tden(148, t_x)
graph twoway (line t_d t_x, clstyle(p1)) (line t_d2 t_x, clstyle(p2)), ///
   xtitle("t-statistic for slope coeff from many samples") ///
   legend( label(1 "t from simulation") label(2 "t with N-k dof"))
graph export simulation10fig.wmf, replace

* To check power
* We tested H0: beta2 = 2 when this is d.g.p. value
* To test power against H0: beta2 = 2.1 change d.g.p. value to 2.1
* Choose 2.1 as a bit more than one standard error away from 2
clear
set seed 54321
postfile simalternative b se t r p using simresults, replace
forvalues i = 1/$numsims {
    drop _all
    quietly set obs $numobs
    quietly generate double x = rchi2(1)   
    quietly generate y = 1 + 2.1*x + rchi2(1)-1     // demeaned chi^2 error 
    quietly regress y x
    scalar b2 =_b[x]
    scalar se2 = _se[x]
    scalar t2 = (_b[x]-2)/_se[x]
    scalar r2 = abs(t2)>invttail($numobs-2,.025)
    scalar p2 = 2*ttail($numobs-2,abs(t2))
    post simalternative (b2) (se2) (t2) (r2) (p2)
}
postclose simalternative
* Analyze the simulation results
use simresults, clear
summarize
mean b se t r p 

***** USE SIMULATE INSTEAD

* Program for finite-sample properties of OLS
program chi2data, rclass
    version 14 
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

* Run the program once
set seed 54321
chi2data

* Simulation for finite-sample properties of OLS
set seed 54321
simulate b2f=r(b2) se2f=r(se2) t2f=r(t2) reject2f=r(r2) p2f=r(p2),  ///
reps($numsims) saving(chi2datares, replace) nolegend nodots: chi2data
summarize b2f se2f t2 reject2f p2f
* Summarize results
mean b2f se2f t2 reject2f p2f

* Power against 2.1
clear
program chi2datapower, rclass
    version 14 
    drop _all
    set obs $numobs
    generate double x = rchi2(1)   
    generate y = 1 + 2.1*x + rchi2(1)-1     // Change to 2.1
    regress y x
    return scalar b2 =_b[x]
    return scalar se2 = _se[x]
    return scalar t2 = (_b[x]-2)/_se[x]
    return scalar r2 = abs(return(t2))>invttail($numobs-2,.025)
    return scalar p2 = 2*ttail($numobs-2,abs(return(t2)))
end
set seed 54321
simulate b2f=r(b2) se2f=r(se2) t2f=r(t2) reject2f=r(r2) p2f=r(p2),  ///
reps($numsims) saving(chi2datares, replace) nolegend nodots: chi2datapower
summarize b2f se2f t2 reject2f p2f
mean b2f se2f t2 reject2f p2f

********** CLOSE OUTPUT
* log close
* clear
* exit

