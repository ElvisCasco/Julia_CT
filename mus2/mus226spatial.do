* mus226spatial.do  for Stata 17
* Based on mus2_spatial.do 

capture log close
* log using mus226spatial.txt, text replace

********** OVERVIEW OF mus226spatial.do **********

* Stata program 
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Second Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 26 - SPATIAL REGRESSION
* 26.3: GEOSPATIAL DATA
* 26.4: THE SPATIAL WEIGHTING MATRIX
* 26.5: ORDINARY LEAST-SQUARES REGRESSION AND TEST FOR SPATIAL CORRELATION
* 26.7: SAR MODELS
* 26.8: SPATIAL IVS

* To run you need files
*  mus226georgia.dta 
*  mus226georgia_shp.dta
* in your directory

* community-contributed command
*    x_ols  
*  (at http://economics.uwo.ca/people/conley_docs/code_to_download_gmm.html) 
* is used

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* Files mus226georgia.dta and mus226georgia_shp.dta
* are extracts for state of Georgia from 
* Stata example dataset files homicide1990.dta and homicide1990_shp.dta

/* Not included - The following created the Georgia data from Stata's data

* Read in homicide1990.dta and pick only Georgia 
use homicide1990, clear
summarize
list _ID if sname == "Georgia", clean noobs
keep if _ID==2342 | _ID== 2161 |  _ID==2694 |  _ID==2648 |  _ID==2542 | ///
 _ID==2555 |  _ID==2620 |  _ID==2311 |  _ID==2106 |  _ID==2499 |  _ID==2356 |  ///
 _ID==2309 |  _ID==2127 |  _ID==2586 |  _ID==2199 |  _ID==2680 |  _ID==2105 |  ///
 _ID==2387 |  _ID==2176 |  _ID==2693 |  _ID==2517 |  _ID==2190 |  _ID==2533 |  ///
 _ID==2462 |  _ID==2287 |  _ID==2210 |  _ID==2476 |  _ID==2168 |  _ID==2494 |  ///
 _ID==2692 |  _ID==2358 |  _ID==2167 |  _ID==2430 |  _ID==2109 |  _ID==2285 |  ///
 _ID==2483 |  _ID==2600 |  _ID==2624 |  _ID==2598 |  _ID==2353 |  _ID==2193 |  ///
 _ID==2512 |  _ID==2250 |  _ID==2587 |  _ID==2331 |  _ID==2588 |  _ID==2159 |  ///
 _ID==2394 |  _ID==2194 |  _ID==2200 |  _ID==2329 |  _ID==2532 |  _ID==2605 |  ///
 _ID==2326 |  _ID==2722 |  _ID==2275 |  _ID==2486 |  _ID==2521 |  _ID==2649 |  ///
 _ID==2469 |  _ID==2134 |  _ID==2591 |  _ID==2461 |  _ID==2110 |  _ID==2413 |  ///
 _ID==2104 |  _ID==2475 |  _ID==2611 |  _ID==2696 |  _ID==2389 |  _ID==2261 |  ///
 _ID==2479 |  _ID==2637 |  _ID==2582 |  _ID==2103 |  _ID==2691 |  _ID==2271 |  ///
 _ID==2652 |  _ID==2497 |  _ID==2386 |  _ID==2396 |  _ID==2398 |  _ID==2318 |  ///
 _ID==2682 |  _ID==2590 |  _ID==2659 |  _ID==2678 |  _ID==2401 |  _ID==2254 |  ///
 _ID==2547 |  _ID==2695 |  _ID==2251 |  _ID==2400 |  _ID==2631 |  _ID==2129 |  ///
 _ID==2525 |  _ID==2516 |  _ID==2511 |  _ID==2664 |  _ID==2531 |  _ID==2596 |  ///
 _ID==2537 |  _ID==2155 |  _ID==2269 |  _ID==2549 |  _ID==2700 |  _ID==2338 |  ///
 _ID==2395 |  _ID==2383 |  _ID==2556 |  _ID==2107 |  _ID==2603 |  _ID==2676 |  ///
 _ID==2669 |  _ID==2226 |  _ID==2558 |  _ID==2458 |  _ID==2452 |  _ID==2239 |  ///
 _ID==2572 |  _ID==2623 |  _ID==2262 |  _ID==2653 |  _ID==2177 |  _ID==2575 |  ///
 _ID==2609 |  _ID==2447 |  _ID==2641 |  _ID==2150 |  _ID==2100 |  _ID==2595 |  ///
 _ID==2247 |  _ID==2344 |  _ID==2548 |  _ID==2671 |  _ID==2363 |  _ID==2428 |  ///
 _ID==2457 |  _ID==2651 |  _ID==2490 |  _ID==2446 |  _ID==2524 |  _ID==2316 |  ///
 _ID==2108 |  _ID==2402 |  _ID==2225 |  _ID==2336 |  _ID==2304 |  _ID==2274 |  ///
 _ID==2307 |  _ID==2328 |  _ID==2314 |  _ID==2322 |  _ID==2399 |  _ID==2644 |  ///
 _ID==2235 |  _ID==2455 |  _ID==2471 |  _ID== 2222 
summarize
save mus226georgia, replace

* Read in and summarize dataset with all coordinate boundaries for each county
use homicide1990_shp, clear
summarize 
keep if _ID==2342 | _ID== 2161 |  _ID==2694 |  _ID==2648 |  _ID==2542 | ///
 _ID==2555 |  _ID==2620 |  _ID==2311 |  _ID==2106 |  _ID==2499 |  _ID==2356 |  ///
 _ID==2309 |  _ID==2127 |  _ID==2586 |  _ID==2199 |  _ID==2680 |  _ID==2105 |  ///
 _ID==2387 |  _ID==2176 |  _ID==2693 |  _ID==2517 |  _ID==2190 |  _ID==2533 |  ///
 _ID==2462 |  _ID==2287 |  _ID==2210 |  _ID==2476 |  _ID==2168 |  _ID==2494 |  ///
 _ID==2692 |  _ID==2358 |  _ID==2167 |  _ID==2430 |  _ID==2109 |  _ID==2285 |  ///
 _ID==2483 |  _ID==2600 |  _ID==2624 |  _ID==2598 |  _ID==2353 |  _ID==2193 |  ///
 _ID==2512 |  _ID==2250 |  _ID==2587 |  _ID==2331 |  _ID==2588 |  _ID==2159 |  ///
 _ID==2394 |  _ID==2194 |  _ID==2200 |  _ID==2329 |  _ID==2532 |  _ID==2605 |  ///
 _ID==2326 |  _ID==2722 |  _ID==2275 |  _ID==2486 |  _ID==2521 |  _ID==2649 |  ///
 _ID==2469 |  _ID==2134 |  _ID==2591 |  _ID==2461 |  _ID==2110 |  _ID==2413 |  ///
 _ID==2104 |  _ID==2475 |  _ID==2611 |  _ID==2696 |  _ID==2389 |  _ID==2261 |  ///
 _ID==2479 |  _ID==2637 |  _ID==2582 |  _ID==2103 |  _ID==2691 |  _ID==2271 |  ///
 _ID==2652 |  _ID==2497 |  _ID==2386 |  _ID==2396 |  _ID==2398 |  _ID==2318 |  ///
 _ID==2682 |  _ID==2590 |  _ID==2659 |  _ID==2678 |  _ID==2401 |  _ID==2254 |  ///
 _ID==2547 |  _ID==2695 |  _ID==2251 |  _ID==2400 |  _ID==2631 |  _ID==2129 |  ///
 _ID==2525 |  _ID==2516 |  _ID==2511 |  _ID==2664 |  _ID==2531 |  _ID==2596 |  ///
 _ID==2537 |  _ID==2155 |  _ID==2269 |  _ID==2549 |  _ID==2700 |  _ID==2338 |  ///
 _ID==2395 |  _ID==2383 |  _ID==2556 |  _ID==2107 |  _ID==2603 |  _ID==2676 |  ///
 _ID==2669 |  _ID==2226 |  _ID==2558 |  _ID==2458 |  _ID==2452 |  _ID==2239 |  ///
 _ID==2572 |  _ID==2623 |  _ID==2262 |  _ID==2653 |  _ID==2177 |  _ID==2575 |  ///
 _ID==2609 |  _ID==2447 |  _ID==2641 |  _ID==2150 |  _ID==2100 |  _ID==2595 |  ///
 _ID==2247 |  _ID==2344 |  _ID==2548 |  _ID==2671 |  _ID==2363 |  _ID==2428 |  ///
 _ID==2457 |  _ID==2651 |  _ID==2490 |  _ID==2446 |  _ID==2524 |  _ID==2316 |  ///
 _ID==2108 |  _ID==2402 |  _ID==2225 |  _ID==2336 |  _ID==2304 |  _ID==2274 |  ///
 _ID==2307 |  _ID==2328 |  _ID==2314 |  _ID==2322 |  _ID==2399 |  _ID==2644 |  ///
 _ID==2235 |  _ID==2455 |  _ID==2471 |  _ID==2222 
summarize
save mus226georgia_shp, replace

* Need to modify spset to attach to correct shapefile
use mus226georgia, clear
spset
spset, modify shpfile(mus226georgia_shp)
spset
save mus226georgia, replace

*/ 

********** 26.3: GEOSPATIAL DATA

* Read in Georgia homicide data and summarize
qui use mus226georgia, clear
describe _ID _CX _CY cname hrate ln_population poverty
summarize _ID _CX _CY cname hrate ln_population poverty

* Not included - county 2108 
list _ID _CX _CY cname hrate if _ID==2108, clean

* Attach to correct shapefile using spset, modify, and save
spset, modify shpfile(mus226georgia_shp)
save mus226georgia, replace

* Read in shapefile with coordinates for Georgia counties and summarize
qui use mus226georgia_shp, clear
summarize

* List the county boundary coordinates for county 2108
list _ID _X _Y shape_order if _ID == 2108, clean

* Provide a heat map for homicide rates in Georgia counties
qui use mus226georgia, clear
grmap, activate
grmap hrate, clmethod(custom) clbreaks(0 10 20 40) legend(pos(1)) fcolor(Greys)

********** 26.4: THE SPATIAL WEIGHTING MATRIX

* Compute Euclidean distance between observations using planar coordinates
spdistance 2100 2103
display "Distance = " sqrt((-83.40281-(-84.96217))^2+(34.87811-34.80377)^2)

capture spmatrix drop W 

* Create & summarize weighting matrix W - contiguity with spectral normalization
spmatrix create contiguity W
spmatrix summarize W

* Not included - properties of W
* Shows each nonzero entry is 0.1651122936 and id=1 has 2 and id=2 has 4
spmatrix matafromsp Wmata id = W
mata: Wmata[1,1]
mata: rowsum(Wmata[|1,1 \ 1,159|])
mata: rowsum(Wmata[|2,1 \ 2,159|])
mata: Wmata[|1,1 \ 1,159|]  // Row 1
mata: Wmata[|2,1 \ 2,159|]  // Row 2

* Create spatial lag Wy 
spgenerate Whrate = W*hrate
summarize hrate Whrate
correlate hrate Whrate

* Not included
regress hrate ln_population poverty, vce(robust) noheader
regress hrate ln_population poverty Whrate, vce(robust) noheader

********** 26.5: OLS REGRESSION AND TEST FOR SPATIAL CORRELATION

* OLS estimation using regress
regress hrate ln_population poverty, vce(robust)

* Moran test for spatial correlation following OLS
qui regress hrate ln_population poverty
estat moran, errorlag(W)

capture drop window epsilon dis1 dis2
capture drop const cutoff1 cutoff2

* OLS with Conley spatial HAC standard errors using user addon x_ols
generate const = 1
generate cutoff1 = 1.04 // One standard deviation of _CX
generate cutoff2 = 1.19 // One standard deviation of _CY
x_ols _CX _CY cutoff1 cutoff2 hrate ln_population poverty const, xreg(3) coord(2)

* Not included - different cutoffs
replace cutoff1 = 2
replace cutoff2 = 2
capture drop window epsilon dis1 dis2
x_ols _CX _CY cutoff1 cutoff2 hrate const ln_population poverty, xreg(3) coord(2)

********** 26.7: SPATIAL AUTOREGRESSIVE MODELS

* OLS estimation using spregress
spregress hrate ln_population poverty, gs2sls heteroskedastic

* SAR(1) in mean: Errors independent
spregress hrate ln_population poverty, gs2sls dvarlag(W) heteroskedastic

capture drop one Wone Wln_population Wpoverty

* Next uses spgenerate Whrate = W*hrate created earlier

* SAR(1) in mean: Errors independent estimated using ivregress
generate one = 1
spgenerate Wone = W*one
spgenerate Wln_population = W*ln_population
spgenerate Wpoverty = W*poverty
qui ivregress 2sls hrate ln_population poverty ///
    (Whrate = W*one Wln_population Wpoverty), vce(robust)
estimates store IVREG1 
qui spregress hrate ln_population poverty, gs2sls dvarlag(W) ///
    heteroskedastic impower(1)
estimates store SPREG1
estimates table IVREG1 SPREG1, b(%9.4f) se eq(1)

* Not included Out two lags
spgenerate W2ln_population = W*Wln_population
spgenerate W2poverty = W*Wpoverty
spgenerate W2one = W*Wone
ivregress 2sls hrate ln_population poverty (Whrate = W*one Wln_population ///
    Wpoverty W2*one W2ln_population W2poverty), vce(robust)
spregress hrate ln_population poverty, gs2sls dvarlag(W) heteroskedastic

* SAR(1) in mean: ml with errors iid normal
spregress hrate ln_population poverty, ml dvarlag(W) 

* Impact multipliers following SAR(1) in mean
qui spregress hrate ln_population poverty, gs2sls dvarlag(W) heteroskedastic
estat impact

* SAR(1) in error: errors independent
spregress hrate ln_population poverty, gs2sls errorlag(W) heteroskedastic

// Not included
* SAR(1) in error: ml with errors iid normal
spregress hrate ln_population poverty, ml errorlag(W)

// Not included
* OLS and SAR(1) model comparison
qui spregress hrate ln_population poverty, gs2sls heteroskedastic
estimates store OLS
qui spregress hrate ln_population poverty, gs2sls dvarlag(W) heteroskedastic
estimates store MEAN
qui spregress hrate ln_population poverty, ml dvarlag(W) 
estimates store MEAN_ml
qui spregress hrate ln_population poverty, gs2sls errorlag(W) heteroskedastic
estimates store ERROR
qui spregress hrate ln_population poverty, ml errorlag(W)
estimates store ERROR_ml
estimates table OLS MEAN MEAN_ml ERROR ERROR_ml, b(%9.3f) se eq(1)

* SARAR(1,1) model with additionally SAR in X: gs2sls
spregress hrate ln_population poverty, gs2sls dvarlag(W) errorlag(W) ///
    ivarlag(W:ln_population poverty) heteroskedastic

* Not included - add Wxgini but not gini
spregress hrate ln_population poverty, gs2sls dvarlag(W) errorlag(W) ///
   ivarlag(W:ln_population poverty gini) heteroskedastic

// Not included
* Impact multipliers following SARAR(1,1)
estat impact

********** 26.8: SPATIAL INSTRUMENTAL VARIABLES

* Endogenous regressors in SAR in mean model
spivregress hrate ln_population (poverty=gini), gs2sls dvarlag(W) heteroskedastic

********** END
