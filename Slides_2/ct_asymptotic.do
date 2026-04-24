* ct_asymptotic.do  based on mus04p1sim  January 2013 for Stata version 12.0

********** OVERVIEW OF ct_asymptotic.do **********

log using ct_asymptotic.txt, text replace

* STATA Program 
* For A. Colin Cameron and Pravin Trivedi "Lectures in Microeconometrics"
* Asymptotic Theory - LLN and CLT

* Requires no dataset

********** SETUP **********

set more off
version 12.0
set mem 10m
set scheme s1mono  /* Graphics scheme */

******* ASYMPTOTIC THEORY - LLN and CLT in IID CASE

* Draw from uniform with population mean 0.5
* Demonstrate LLN by finding average for a very large sample
* Demonstrate CLT by simulating to obtain many averages

* Small sample: sample mean differs from population mean
quietly set obs 30 
set seed 10101
quietly generate x = runiform()
mean x

* Consistency: Large sample: sample mean is very close to population mean
clear all 
quietly set obs 100000
set seed 10101
quietly generate x = runiform()
mean x

* Central limit theorem
* Write program to obtain sample mean for one sample of size numobs (= 30)
program onesample, rclass
    args numobs
    drop _all
    quietly set obs `numobs'
    generate x = runiform()
    summarize x
    return scalar meanforonesample = r(mean)
end
* Run this program 10,000 times to get 10,000 sample means
quietly simulate xbar = r(meanforonesample), seed(10101) reps(10000) nodots: ///
   onesample 30
* Summarize the 10,000 sample means
summarize xbar

* Draw histogram and kernel density estimate for z=(xbar-mu)/(sigma/sqrt(N))
generate z = (xbar-0.500) / (sqrt(1/12)/sqrt(30))
quietly histogram z, normal xtitle("z from many samples")
* graph export ct_asymptoticclt.wmf, replace

********** CLOSE OUTPUT
* log close
* clear
* exit


