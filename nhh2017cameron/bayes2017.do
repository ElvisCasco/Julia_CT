* bayes2017.do  August 2017

log using bayes2017.txt, text replace

********** OVERVIEW OF bayes2017.do **********

* PROBIT MCMC
*  1. Generate Data and MLE
*  2. Sata 14 bayesmh command
*  3. Produce graph for normal-normal example 
*  4. Own code with METROPOLIS HASTINGS RANDOM WALK CHAIN & UNINFORMATIVE PRIOR
*  5. Own code with DATA AUGMENTATION AND GIBBS SAMPLER

* To run you need files
*   no files    
* in your directory
* No Stata user-written commands are used

********** SETUP **********

set more off
version 14
clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Only simulated data used
* Stata program based in part on Koop's MATLAB program chapter9b.m
* Probit example of Koop (2003) chapter 9.3 "Bayesian Econometrics"

********** 1. GENERATE THE DATA FOR PROBIT WITH ONE REGRESSOR **********

* Generate artificial data set for probit illustration 
*  - explanatory variable x ~ N[0,1]
*  - dependent variable   y = 1(0.5*x + e > 0) 
*                     for e ~ N[0,1]

* Generate data N = 100  Pr[y=1|x] = PHI(0 + 0.5*x)
clear 
set obs 100
set seed 1234567
gen x = rnormal(0,1)
gen ystar = 0.5 + 1*x + rnormal(0,1)
gen y = (ystar > 0)
gen cons = 1
summarize
save bayes2017sample.dta, replace

* Estimate model by MLE
probit y x

********** 2. ANALYZE USING BAYES or BAYESMH COMMAND

* Following the same as version 15 command bayes, rseed(10101): probit y x
bayesmh y x, likelihood(probit) prior({y: }, normal(0,10000)) rseed(10101)

* Now different prior for each parameter and a tight prior
bayesmh y x, likelihood(probit) prior({y:x}, normal(0.5,0.01)) ///
  prior({y:_cons}, normal(0.5,0.01)) rseed(10101)

* Diagostics
bayesmh y x, likelihood(probit) prior({y: }, normal(0,10000)) rseed(10101)
bayesgraph diagnostics {y:x}
graph export bayesdiagnostics.wmf, replace

* Summary statistics
bayesstats summary {y:x}

* Probability that slope is in range 0.4 to 0.6
bayestest interval {y:x}, lower(0.4) upper(0.6)

* Effective sample size
bayesstats ess

* Save retained draws
bayesmh y x, likelihood(probit) prior({y: }, normal(0,10000)) rseed(10101) ///
   saving(bayesmhout.dta, replace)
use bayesmhout.dta, clear
summarize

********** 3. BAYES NORMAL-NORMAL GRAPH

clear
set obs 1000
gen x = -2 + _n*20/1000   // from -2 to 18
generate prior = normalden(x,5,sqrt(3))
generate posterior = normalden(x,8,sqrt(1.2))
generate likelihood = normalden(x,10,sqrt(2))
graph twoway (scatter prior x, msize(small) lstyle(p4)) ///
  (scatter likelihood x, msize(small) lstyle(p2)) (scatter posterior x, msize(small) lstyle(p3)) 
graph export bayesnormalnormal.wmf, replace

********** 4. PROBIT BAYESIAN WITH RANDOM WALK CHAIN METROPOLIS HASTNGS **********

clear
use bayes2017sample.dta

* Start the Bayesian MCMC

* Globals for number of reps and a key tuning parameter
 
global s1 10000      // number of retained reps
global s0 10000      // number of burnin reps
global sdscale 0.25  // use random walk b + $sdscale * N(0, I)

* Mata to obtain the posterior draws of b

mata
  // Create y vector and X matrix from Stat data set using st_view command
  st_view(y=., ., "y")            // dependent
  st_view(X=., ., ("cons", "x"))  // regressors
  Xnames = ("cons", "x")          // used to label output

  // Calculate a few quantities outside the loop for later use
  n = rows(X)
  k = cols(X)
  ones = J(n,1,1)

  // Specify the number of replications 
  s0 = $s0     // number of burnin reps
  s1 = $s1     // number of retained reps
  s = s0+s1    // total reps

  // Store all draws and MH acceptance ratein the following matrices
  b_all = J(k,s1,0)
  accept_all = J(1,s1,0) 

  // Initialization
  bdraw = J(k,1,0)      // starting b value is vector of zeroes
  lpostdraw = -1*10^10  // starting value of ln(posterior) is small 
                        // so accept initial MH draw
   
  // Now do Metropolis-Hastings loop and make the posterior draws

  for (irep=1; irep<=s; irep++) {

     // Draw new candidate value of b from MH random walk chain
     bcandidate = bdraw + $sdscale*rnormal(k,1,0,1)

     // Note: For different data you may need to change the global sdscale
     // And best is bcandidate = bdraw + z  where z ~ N(0, posterior variance of b)
 
     // Compute the log posterior at the candidate value of b
     // The assumed prior for b is uninformative 
     // so the posterior is proportional to the usual probit liklihood
     probitprob = normal(X*bcandidate)
     lpostcandidate = ones'( y:*ln(probitprob) + (ones - y):*ln(ones - probitprob) )
 
     // Accept the candidate draw on basis of posterior probability ratio 
     // if  uniform > (posterior(bcandidate) / posterior(bdraw))
     // where bcandidate is current b  and  bdraw is previous b
     // Taking logs the rule is the same as
     // if  ln(uniform) > (lpostcandidate - lpostdraw)
     laccprobability = lpostcandidate - lpostdraw
     accept = 0
     if ( ln(runiform(1,1)) < laccprobability ) {
       lpostdraw = lpostcandidate
       bdraw = bcandidate    
       accept = 1 
       }

    // Store the draws after burn-in of b and whether accept draw 
    if (irep>s0) {
        // after discarding burnin, store all draws
        j = irep-s0
        b_all[.,j] = bdraw         // These are the posterior draws
        accept_all[.,j] = accept   // These are one if new draw accepted
    } 

  }                
  
  // End MH loop 

  // Pass results back to Stata
  // This bit works only for k = 2 (intercept plus one slope)
  beta = b_all'
  accept = accept_all'
  st_addvar("float", ("beta1", "beta2", "accept"))
  // The following line is needed for conformability 
  // It assumes that s1 (number of draws) exceeds the original sample size 
  stata("set obs $s1")
  st_store(., ("beta1", "beta2"), beta)  
  st_store(., ("accept"), accept)  
end

** Analyze the results

summarize
save bayesprobitmh.dta, replace
rename beta2 b

* The posterior mean and posterior standard deviation of the slope coefficient
quietly summarize b
display "posterior mean = " r(mean) "  and posterior standard deviation = " r(sd)
summarize b
centile b, centile(2.5, 97.5)
kdensity b, normal
graph export bayesprobitmhposterior.wmf, replace

* Plot the posterior draws of b2 and acf
generate s= _n
tsset s
line b s if s < 100
graph export bayeprobitmhdraws.wmf, replace

* Give the corrleations and the acceptance rate in the random walk chain MH 
corrgram b, lags(10)
quietly summarize accept
display "MH acceptance rate = " r(mean) "


********** 5. PROBIT BAYESIAN WITH DATA AUGMENTATION AND GIBBS SAMPLER **********

clear
use bayes2017sample.dta

* Estimate model by MLE
probit y x

* Start the Bayesian MCMC

* Globals for number of reps 

global s1 10000 // number of retained reps

mata
  // Create y vector and X matrix from Stat data set using st_view command
  st_view(y=., ., "y")            // dependent
  st_view(X=., ., ("cons", "x"))  // regressors
  Xnames = ("cons", "x")          // used to label output

  // Calculate a few quantities outside the loop for later use
  n = rows(X)
  k = cols(X)
  Xsquare = X'*X
  Xtxinv = invsym(Xsquare)
  Xtxinvchol = cholesky(Xtxinv)

  v1 = n

  // Specify the number of replications 
  s0 = 10000  // number of burnin reps
  s1 = $s1    // number of retained reps
  s = s0+s1   // total reps

  // store all draws in the following matrices
  y1 = J(s0+s1,1,0)
  b_all = J(k,s1,0)

  // Beta prior is Noninformative

  // choose a starting value for latent data
  ystar = y 

  // Now do Gibbs loop

  for (irep=1; irep<=s; irep++) {

     // Posterior-step:  draw from beta | y* ~ N[bols*, (X'X)^-1]
     // This is using noninformative prior for beta
     bols = Xtxinv*X'*ystar
     b1 = bols
     bdraw = b1 + Xtxinvchol*rnormal(k,1,0,1)  //invnormal(uniform(k,1))

     // Imputation step: make one draw of vector ystar 
     // where for ith observation ystar_i | y,b  is truncated normal
     // Right: If y = 1 we need draw from truncated N[0,1] with ystar > -mu
     // Left:  If y = 0 we need draw from truncated N[0,1] with ystar < -mu
     // The method used here is explained in 
     // Cameron and Trivedi (2009, section 4.4.4) Microeconometrics using Stata

     for (i=1; i<=n; i++) {
     mu = X[i,.]*bdraw
       if (y[i,1]==0) {
           uright = normal(-mu)*uniform(1,1)        
           ystar[i,1] = mu + invnormal(uright)
       }
       else {
           uleft = normal(-mu) + (1-normal(-mu))*uniform(1,1)
           ystar[i,1] = mu + invnormal(uleft)
       }
     }
    
    // Store the draws of b after burn-in plus a diagnostic used in Koop
    if (irep>s0) {
        // after discarding burnin, store all draws
        j = irep-s0
        b_all[.,j] = bdraw    // These are the posterior draws
    } 

  }                
  
  // End Gibbs loop 

  beta = b_all'
  st_addvar("float", ("beta1", "beta2"))
  // The following line is needed for conformability 
  // It assumes that s1 (number of draws) exceeds the original sample size 
  stata("set obs $s1")
  st_store(., ("beta1", "beta2"), beta)

end

** Analyze the results

summarize
save bayesprobitgibbsaugment.dta, replace
rename beta2 b

* Plot the draws and acf
generate s= _n
tsset s
line b s if s < 100
graph export bayesprobitgibbsaugmentdraws.wmf, replace
corrgram b, lags(10)

* Posterior mean etc.
kdensity b, normal
graph export bayesprobitgibbsaugmentposterior.wmf, replace
summarize b
centile b, centile(2.5, 97.5)

********** CLOSE OUTPUT **********

* log close
* clear 
* exit
