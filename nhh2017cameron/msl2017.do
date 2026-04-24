* msl2017.do  August 2017
* By A. Colin Cameron 

log using msl2017.txt, text replace

********** OVERVIEW OF msl2017.do **********

* 1. LOGIT MLE USING STATA ML COMMAND
* 2. RANDOM PARAMETERS LOGIT BY MAX SIMULATED LIKELIHOOD USING STATA ML COMMAND
* 3. RANDOM PARAMETERS LOGIT BY MSL USING STATA ADDON MIXEDLOGIT
* 4. RANDOM PARAMETERS LOGIT BY MSL USING STATA 15 COMMAND ASMIXLOGIT

* To run you need files
*   no files    
* in your directory

* Stata user-written commands 
*     mixlogit
* is used

********** SETUP **********

clear all
set more off
version 14
set linesize 81
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Only simulated data used

**** 1. LOGIT MLE USING STATA ML COMMAND

* Generate the data - Pr[y=1|x] = LAMDA(1 + 1*x)
clear all

* Use Stata 13 random kiss32 number generator 
set rng kiss32

set obs 1000
set seed 10101
gen u = runiform()
gen ulogistic = ln(u) - ln(1-u)
gen x = rnormal(0,2)
gen y = 1 + 1*x + ulogistic > 0
summarize

* ASIDE: Alternative way to draw y, perhaps simpler
* gen y2 = runiform() < exp(1+1*x) / (1 + exp(1+1*x)) 
* summarize y y2

***  Logit ml using logit command

logit y x, nolog

** Logit ml using user-written program and command ml

* ML program lflogit to be called by command ml method lf
program lflogit
  args lnf theta1              // theta1=x'b, lnf=lnf(y)
  tempvar p                    // Will define p to make program more readable
  local y "$ML_y1"             // Define y so program more readable
  generate double `p'   = exp(`theta1')/(1+exp(`theta1')) 
  quietly replace `lnf' = `y'*ln(`p') + (1-`y')*ln(1-`p')
end

* Command ml model including defining y and x
ml model lf lflogit (y = x)
ml maximize

** Logit ml using user-written program and command ml 

* The following program is a variation
* that will be extended to random parameters binary logit 
*    b1 and b2 are separate parameters
*    alternative way to define lnf
*    get robust se's 
program lflogitnew
  args lnf b1 b2        // 1 is intercept and b2 is slope 
  tempvar p f
  local y "$ML_y1"
  gen double `p' = exp(`b1' + `b2') / (1 + exp(`b1' + `b2'))
  quietly generate `f' = `p' if `y'==1
  quietly replace `f' = 1 - `p' if `y'==0
  quietly replace `lnf' = ln(`f')
end

ml model lf lflogitnew (b1: y = ) (b2: x, nocons), vce(robust)
ml init 1 1, copy
ml maximize


**** 2. RANDOM PARAMETERS LOGIT BY MAX SIMULATED LIKELIHOOD USING STATA ML COMMAND

* Generate the data - Pr[y=1] = LAMDA(1 + (1+e)*x)
clear all
set obs 1000
set seed 10101
gen u = runiform()
gen ulogistic = ln(u) - ln(1-u)
gen x = rnormal(0,2)
gen e = rnormal(0,1)
gen y = 1 + (1+e)*x + ulogistic > 0
save msl2017sample.dta, replace
summarize

* Create 100 draws (S=100) from the uniform for each observation (n=1000)
* These will be used to in turn get draws from the normal distribution
set seed 10101
forvalues i = 1/100 {
   gen draws`i' = runiform()
   }

* Program to calculate the log-density using Monte Carlo integration
program lflogitmsl
  args lnf b1 b2 ln_sd       // if use sd then problems if sd < 0
  tempvar p sim_f sim_avef
  local y "$ML_y1"
  local sd = exp(`ln_sd')    // convert back to sd 
  qui gen `sim_avef' = 0 
  set seed 10101
  forvalues d = 1/100 {
    gen double `p' = exp(`b1' + `b2' + `sd'*invnormal(draws`d')*x) ///
                 / (1 + exp(`b1' + `b2' + `sd'*invnormal(draws`d')*x))
    qui gen `sim_f' = `p' if `y'==1
    qui replace `sim_f' = 1 - `p' if `y'==0
    qui replace `sim_avef' = `sim_avef' + `sim_f'/100
    drop `p' `sim_f' 
    }
    qui replace `lnf' = ln(`sim_avef')
end

* The following takes time
* Now calculate the maximum simulated likelihood estimator
ml model lf lflogitmsl (b1: y = ) (b2: x, nocons) (ln_sd:), vce(robust)
ml init 1 1 0, copy
ml maximize, difficult

* And convert back to sd = exp(ln_sd)
nlcom exp(_b[ln_sd:_cons])

* NOTE: The dgp is b1=1 b2=1 sd=1
* Note sure why ln_sd and the se's vary from slides even though b1, b2, lnL the same

**** 3. RANDOM PARAMETERS LOGIT BY MSL USING USER ADDON COMMAND MIXLOGIT 

* Continue using the same dataset as above except drop draws*
drop draws*

* Do regular logit
logit y x, vce(robust)
estimates store logitold

* mixlogit requires having separate data for each of the two choices

* Data before conversion
list y x in 1/5, clean
* Now convert to dataset with data for each alternative
gen id = _n
gen x1 = 0
rename x x2
rename y y2
gen y1 = 1 - y2
reshape long y x, i(id) j(alt)

* See what expanded data set looks like
sum id alt y x
list id alt y x in 1/10, clean

* Confirm that asclogit gives same results as earlier using logit
asclogit y x, case(id) alternatives(alt) vce(robust)
estimates store logitnew

* Confirm that clogit gives same results as earlier using logit
* Need to manually provide a dummy variable for alternative 2
gen d2 = alt==2
clogit y d2 x, group(id)
estimates store logitasc
estimates table logitold logitnew logitasc, b(%10.4f) se

* Now do mixlogit which has similar command structure to clogit
mixlogit y d2, group(id) rand(x) nrep(50)

**** 4. RANDOM PARAMETERS LOGIT BY MSL USING STATA 15 COMMAND ASMIXLOGIT

/*
* The following requires Stata 15 and is commented out 
version 15 
* asmixlogit has similar command structure to asclogit
asmixlogit y, case(id) alternatives(alt) random(x) nolog
*/ 

********** CLOSE OUTPUT
* log close
* clear
* exit
