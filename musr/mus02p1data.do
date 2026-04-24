* mus02p1data.do  Oct 2009 for Stata version 11

cap log close

********** OVERVIEW OF mus02p1data.do **********

* Stata program 
* copyright C 2010 by A. Colin Cameron and Pravin K. Trivedi 
* used for "Microeconometrics Using Stata, Revised Edition" 
* by A. Colin Cameron and Pravin K. Trivedi (2010)
* Stata Press

* Chapter 2
* 2.3: INPUTTING DATA
* 2.4: DATA MANAGEMENT
* 2.5: MANIPULATING DATASETS
* 2.6: GRAPHICAL DISPLAY OF DATA

* To run you need files
*   mus02file1.csv
*   mus02file2.csv
*   mus02file3.txt
*   mus02file4.txt
*   mus02psid92m.txt
*   mus02psid92m.do
* in your directory
* No Stata user-written commands are used

********** SETUP **********

set more off
version 11
clear all
set scheme s1mono  /* Graphics scheme */

************ 2.3: INPUTTING DATA

* Read in dataset from an Internet web site
use http://www.stata-press.com/data/r11/census.dta, clear
clear

* Data input from keyboard
input str20 name age female income
  "Barry" 25 0 40.990
  "Carrie" 30 1 37.000
  "Gary" 31 0 48.000
end

list, clean

* Read data from a csv file that includes variable names using insheet
clear 
insheet using mus02file1.csv
list, clean

* Read data from a csv file without variable names and assign names
clear 
insheet name age female income using mus02file2.csv

* Read data from free-format text file using infile
clear 
infile str20 name age female income using mus02file2.csv
list, clean

* Read data from fixed-format text file using infix
clear
infix str20 name 1-10 age 11-12 female 13 income 14-20 using mus02file3.txt
list, clean

* Read data using infix where an observation spans more than one line
clear
infix str20 name 1-10 age 11-12 female 13 2: income 1-7 using mus02file4.txt

************* 2.4: DATA MANAGEMENT

* Commands to read in data from PSID extract
type mus02psid92m.do

clear
#delimit ;
*  PSID DATA CENTER *****************************************************
   JOBID            : 10654                             
   DATA_DOMAIN      : PSID                              
   USER_WHERE       : ER32000=1 and ER30736 ge 30 and ER
   FILE_TYPE        : All Individuals Data              
   OUTPUT_DATA_TYPE : ASCII Data File                   
   STATEMENTS       : STATA Statements                  
   CODEBOOK_TYPE    : PDF                               
   N_OF_VARIABLES   : 12                                
   N_OF_OBSERVATIONS: 4290                              
   MAX_REC_LENGTH   : 56                                
   DATE & TIME      : November 3, 2003 @ 0:28:35
*************************************************************************
;
insheet 
   ER30001 ER30002 ER32000 ER32022 ER32049 ER30733 ER30734 ER30735 ER30736
    ER30748 ER30750 ER30754
using mus02psid92m.txt, delim("^") clear 
;
destring, replace ;
label variable er30001  "1968 INTERVIEW NUMBER"  ;
label variable er30002  "PERSON NUMBER                 68"  ;
label variable er32000  "SEX OF INDIVIDUAL"  ;
label variable er32022  "# LIVE BIRTHS TO THIS INDIVIDUAL"  ;
label variable er32049  "LAST KNOWN MARITAL STATUS"  ;
label variable er30733  "1992 INTERVIEW NUMBER"  ;
label variable er30734  "SEQUENCE NUMBER               92"  ;
label variable er30735  "RELATION TO HEAD              92"  ;
label variable er30736  "AGE OF INDIVIDUAL             92"  ;
label variable er30748  "COMPLETED EDUCATION           92"  ;
label variable er30750  "TOT LABOR INCOME              92"  ;
label variable er30754  "ANN WORK HRS                  92"  ;

#delimit cr;    //  Change delimiter to default cr

* Data description
describe

* Data summary
summarize

* save mus01psid92m, replace

* Rename variables
rename er32000 sex
rename er30736 age
rename er30748 education
rename er30750 earnings
rename er30754 hours

* Relabel some of the variables
label variable age "AGE OF INDIVIDUAL"
label variable education "COMPLETED EDUCATION"
label variable earnings "TOT LABOR INCOME"
label variable hours "ANN WORK HRS"

* Define the label gender for the values taken by variable sex
label define gender 1 male 2 female
label values sex gender
list sex in 1/2, clean

* Data summary of key variables after renaming
summarize sex age education earnings hours

* List first 2 observations of two of the variables
list age hours in 1/2, clean

* Tabulate all values taken by a single variable
tabulate education

* Replace missing values with missing-data code
replace education = . if education == 0 | education == 99

* Listing of variable including missing value
list education in 46/48, clean  

* Example of data analysis with some missing values
summarize education age

* Replace missing values with missing-data code
replace earnings = . if education >= .
replace hours = . if education >= .

* Assign more than one missing code
replace education = .a if education == 0
replace education = .b if education == 99

* This command will include missing values
list education in 40/60 if education > 16, clean

* This command will not include missing values
list education in 40/60 if education > 16 & education < . , clean

* Summarize cleaned up data
summarize sex age education earnings

* Create identifier using generate command
generate id = _n

* Create new variable using generate command
generate lnearns = ln(earnings)

* Create new variable using egen command
egen aveearnings = mean(earnings) if earnings < .

* Replace existing data using the recode command
recode education (1/11=1) (12=2) (13/15=3) (16/17=4), generate(edcat) 

* Create new variable using bysort: prefix
bysort education: egen aveearnsbyed = mean(earnings)
sort id

* Create indicator variable using generate command with logical operators
generate d1 = earnings > 0 if earnings < .

* Create indicator variable using generate and replace commands
generate d2 = 0
replace d2 = 1 if earnings > 0
replace d2 = . if earnings >= .

* Create indicator variable using recode commands
recode earnings (0=0) (1/999999=1), generate(d3)

* Verify that three methods create the same indicator variable
// summarize d1 d2 d3
summarize d1 

* Create a set of indicator variables using tabulate with generate() option
quietly tabulate edcat, generate(eddummy)
summarize eddummy*

* Set of indicator variables using factor variables - no category is omitted
summarize ibn.edcat


* Create interactive variable using generate commands
generate d1education = d1*education

* Set of interactions using factor variables
summarize i.edcat#c.earnings


* Create demeaned variables
egen double aveage = mean(age)
generate double agedemean = age - aveage
generate double agesqdemean = agedemean^2
summarize agedemean agesqdemean

* Save as Stata data file
save mus02psid92m.dta, replace

* Save as comma-separated values spreadsheet
outsheet age education eddummy* earnings d1 hours using mus02psid92m.csv, comma replace

* Save as formatted text (ascii) file
outfile age education eddummy* earnings d1 hours using mus02psid92m.asc, replace

* Select the sample used in a single command using the if qualifier
summarize earnings lnearns if age >= 35 & age <= 44 

* Select the sample using command keep
keep if (lnearns != .) & (age >= 35 & age <= 44) 
summarize earnings lnearns

* Select the sample using keep and drop commands
use mus02psid92m.dta, clear
keep lnearns age
drop in 1/1000

************* 2.5: MANIPULATING DATASETS

* Commands preserve and restore illustrated
use mus02psid92m.dta, clear
list age in 1/1, noheader clean
preserve
replace age = age + 1000
list age in 1/1, noheader clean
restore
list age in 1/1, noheader clean

* Create first dataset with every third observation
use mus02psid92m.dta, clear
keep if mod(_n,3) == 0 
keep id education earnings
list in 1/4, clean
quietly save merge1.dta, replace

* Create second dataset with every second observation
use mus02psid92m.dta, clear
keep if mod(_n,2) == 0 
keep id education hours
list in 1/4, clean
quietly save merge2.dta, replace

* Merge two datasets with some observations and variables different
clear
use merge1.dta
sort id
merge 1:1 id using merge2.dta    
sort id
list in 1/4, clean

* Append two datasets with some observations and variables different
clear
use merge1.dta
append using merge2.dta
sort id
list in 1/5, clean

************* 2.6: GRAPHICAL DISPLAY OF DATA

use mus02psid92m.dta, clear
twoway scatter lnearns hours
graph export mus02two1.eps, replace

* More advanced graphics command with two plots and with several options
graph twoway (scatter lnearns hours, msize(small))     ///
   (lfit lnearns hours, lwidth(medthick)),             ///
   title("Scatterplot and OLS fitted line")       
graph export mus02two2.eps, replace

/*
* More advanced graph example 
graph twoway (scatter lnearns hours, msize(small))    ///
  (lfit lnearns hours, lwidth(medthick)),             ///
  title("Scatterplot and OLS fitted line") 
graph export mus02examplefigure.eps, replace
*/

use mus02psid92m.dta, clear
label define edtype 1 "< High School" 2 "High School" 3 "Some College" 4 "College Degree"
label values edcat edtype

* Box and whisker plot of single variable over several categories
graph box hours, over(edcat) scale(1.2) marker(1,msize(vsmall))   ///
  ytitle("Annual hours worked by education") yscale(titlegap(*5)) 
quietly graph export mus02boxfig.eps, replace

histogram lnearns

* Histogram with bin width and start value set
histogram lnearns, width(0.25) start(4.0)
graph export mus02hist.eps, replace

kdensity lnearns
* Kernel density plot with bandwidth set and fitted normal density overlaid
kdensity lnearns, bwidth(0.20) normal n(4000)
graph export mus02kd1.eps, replace

* Histogram and nonparametric kernel density estimate 
histogram lnearns if lnearns > 0, width(0.25) kdensity             ///
  kdenopts(bwidth(0.2) lwidth(medthick))                               ///
  plotregion(style(none)) scale(1.2)                                   ///
  title("Histogram and density for log earnings")                      ///
  xtitle("Log annual earnings", size(medlarge)) xscale(titlegap(*5))   ///
  ytitle("Histogram and density", size(medlarge)) yscale(titlegap(*5)) 
graph export mus02histdens.eps, replace

* Simple two-way scatterplot
scatter lnearns hours

* Two-way scatterplot and quadratic regression curve with 95% ci for y|x  
twoway (qfitci lnearns hours, stdf)  (scatter lnearns hours, msize(small))
graph export mus02scatter2.eps, replace


* Scatterplot with lowess and local linear nonparametric regression
graph twoway (scatter lnearns hours, msize(tiny))                 ///
  (lowess lnearns hours, clstyle(p1) lwidth(medthick))            /// 
  (lpoly lnearns hours, kernel(epan2) degree(1) lwidth(medthick)  ///
  bwidth(500)), plotregion(style(none))                           ///
  title("Scatterplot, lowess, and local linear regression")        ///
  xtitle("Annual hours", size(medlarge))                          /// 
  ytitle("Natural logarithm of annual earnings", size(medlarge)) ///
  legend(pos(4) ring(0) col(1)) legend(size(small))              ///
  legend(label(1 "Actual Data") label(2 "Lowess") label(3 "Local linear"))
quietly graph export mus02twoway.eps, replace

* Multiple scatterplots
label variable age "Age"
label variable lnearns "Log earnings"
label variable hours "Annual hours"
graph matrix lnearns hours age, by(edcat) msize(small) 

*********** CLOSE OUTPUT
