* mus202data.do  for Stata 17

capture log close

********** OVERVIEW OF mus202data.do **********

* Stata program
* copyright C 2022 by A. Colin Cameron and Pravin K. Trivedi
* used for "Microeconometrics Using Stata, Second Edition"
* by A. Colin Cameron and Pravin K. Trivedi (2022)
* Stata Press

* Chapter 2
* 2.3: INPUTTING DATA
* 2.4: DATA MANAGEMENT
* 2.5: MANIPULATING DATASETS
* 2.6: GRAPHICAL DISPLAY OF DATA

* To run you need files
*   mus202file1.csv
*   mus202file2.csv
*   mus202file3.txt
*   mus202file4.txt
*   mus202psid92m.txt
* in your directory

* No community-contributed commands are used

********** SETUP **********

clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

********** DATA DESCRIPTION **********

* mus202psid92m.txt is authors' extract from the 1992 PSID for males 30-50 years

************ 2.2: TYPES OF DATA

* Example of numerical error
display %25.20f 0.3+0.6+0.1          ///
    _n %25.20f 0.3+0.1+0.6           ///
    _n %25.20f (0.3+0.6+0.1)-(0.3+0.1+0.6)

* Examples of formats
display %15.2f 123.456 %15.2g 123.456 %15.2e 123.456
display %15.5f 123.456 %15.5g 123.456 %15.5e 123.456

************ 2.3: INPUTTING DATA

* Read in dataset from an Internet website
use http://www.stata-press.com/data/r17/census

* Data input from keyboard using input
clear
input str20 name age female income
    "Barry" 25 0 40.990
    "Carrie" 30 1 37.000
    "Gary" 31 0 48.000
end

list, clean

* Read data from a .csv file that includes variable names using import delimited
clear
import delimited using mus202file1.csv
list, clean

* Read data from a .csv file without variable names, and assign names
clear
import delimited name age female income using mus202file2.csv

* Read data from free-format text file using infile
clear
infile str20 name age female income using mus202file2.csv
list, clean

* Read data from fixed-format text file using infix
clear
infix str20 name 1-10 age 11-12 female 13 income 14-20 using mus202file3.txt
list, clean

* Read data using infix where an observation spans more than one line
clear
infix str20 name 1-10 age 11-12 female 13 2: income 1-7 using mus202file4.txt

************* 2.4: DATA MANAGEMENT

* Commands to read in data from PSID extract in a delimited file
type mus202psid92m.do

* Commands to read in data from PSID extract
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
import delimited
   er30001 er30002 er32000 er32022 er32049 er30733 er30734 er30735 er30736
    er30748 er30750 er30754
using mus202psid92m.txt, delim("^") clear 
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

* Rename variables
rename er32000 sex
rename er30736 age
rename er30748 education
rename er30750 earnings
rename er30754 hours

* Relabel some of the variables
label variable age "Age of individual"
label variable education "Completed education"
label variable earnings "Total labor income"
label variable hours "Annual work hours"

* Define the label gender for the values taken by variable sex
label define gender 1 male 2 female
label values sex gender
list sex in 1/2, clean

* Data summary of key variables after renaming
summarize sex age education earnings hours

* List first 2 observations of two of the variables
list age hours in 1/2, clean

* Change display format of variable hours
format %15.2f hours
list age hours in 1/2, clean

* Count number of observations satisfying a given condition
count if education < 12

* Tabulate all values taken by a single variable
tabulate education

* Replace missing values with missing-data code
replace education = . if education == 0 | education == 99

* Listing of variable including missing value
list education in 46/48, clean

* Example of data analysis with some missing values
summarize education age

* Replace missing values with missing-data code
replace earnings = . if missing(education)
replace hours = . if missing(education)

/*
The next piece appeared at text in the book. It is an an alternative
to earlier command replace education = . if education == 0 | education == 99
* Assign more than one missing code
replace education = .a if education == 0
replace education = .b if education == 99
*/

* This command will include missing values
list education in 40/60 if education > 16, clean

* This command will not include missing values
list education in 40/60 if education > 16 & !missing(education), clean

* Summarize cleaned-up data
summarize sex age education earnings

* Create identifier using generate command
generate id = _n

* Create new variable using generate command
generate lnearns = ln(earnings)

* Create new variable using egen command
egen aveearnings = mean(earnings) if !missing(earnings)

* Replace existing data using the recode command
recode education (1/11=1) (12=2) (13/15=3) (16/17=4), generate(edcat)

* Convert numeric variable to string using decode and preexisting label
list sex in 1, nolabel clean    // List the numeric value
list sex in 1, clean            // List the associated label
decode sex, generate(str_sex)   // Make a new string variable
list str_sex in 1, clean        // List the new string variable

* Convert string variable to numeric with label using encode
encode str_sex, generate(num_sex)  // Make a new numeric variable
list num_sex in 1, nolabel clean   // List the numeric value
list num_sex in 1, clean           // List the associated label

* Create new variable using bysort prefix
bysort education: egen aveearnsbyed = mean(earnings)
sort id

* Create indicator variable using generate command with logical operators
generate d1 = earnings > 0 if !missing(earnings)

summarize d1

* Following was not included - alternative eaysto create indicator variable
* Create indicator variable using generate and replace commands
generate d2 = 0
replace d2 = 1 if earnings > 0
replace d2 = . if earnings >= .
* Create indicator variable using recode commands
recode earnings (0=0) (1/999999=1), generate(d3)
* Verify that three methods create the same indicator variable
summarize d1 d2 d3
drop d2 d3

* Create a set of indicator variables using tabulate with generate() option
qui tabulate edcat, generate(eddummy)
summarize eddummy*

* Set of indicator variables using factor variables - no category is omitted
summarize i.edcat


* Create interactive variable using generate commands
generate d1education = d1*education

* Set of interactions using factor-variable operators
summarize i.edcat#c.earnings

* Create demeaned variables
egen double aveage = mean(age)
generate double agedemean = age - aveage
generate double agesqdemean = agedemean^2
summarize agedemean agesqdemean

* Save as Stata data file (also used in semiparametric regression chapter)
save cleaneddata, replace

* Save as Stata data file readable by versions 12 and above
saveold cleaneddata, version(12) replace

* Save as comma-separated values spreadsheet
export delimited age education eddummy* earnings d1 hours ///
    using cleaneddata.csv, replace

* Save as formatted text (ascii) file
outfile age education eddummy* earnings d1 hours using cleaneddata.asc, replace

* Select the sample used in a single command using the if qualifier
summarize earnings lnearns if age >= 35 & age <= 44

* Select the sample using command keep
keep if (!missing(lnearns)) & (age >= 35 & age <= 44)
summarize earnings lnearns

* Select the sample using keep and drop commands
use cleaneddata, clear
keep lnearns age
drop in 1/1000

************* 2.5: MANIPULATING DATASETS

* Commands preserve and restore illustrated
use mus202psid92m, clear
list age in 1/1, noheader clean
preserve
replace age = age + 1000
list age in 1/1, noheader clean
restore
list age in 1/1, noheader clean

* Create a new data frame
frame
frame rename default first
frame copy first second

* Change to the new data frame and manipulate
frame change second
replace age = age + 1000
list age in 1/1, noheader clean
frame first: list age in 1/1, noheader clean

* Revert back to the original data frame
frame change first
list age in 1/1, noheader clean

* Collapse to dataset of medians of earnings and age for each value of edcat
preserve
keep if !missing(earnings) & !missing(age)
collapse (median) earnings age, by(edcat)
list, clean
restore

* Expand dataset with number of duplicate observations = edcat if edcat < 4
preserve
qui keep if !missing(earnings) & !missing(age)
qui collapse (median) earnings age, by(edcat)
expand edcat if edcat < 4
list, clean
restore

* Create first dataset with every third observation
use mus202psid92m, clear
keep if mod(_n,3) == 0
keep id education earnings
list in 1/4, clean
qui save merge1, replace

* Create second dataset with every second observation
use mus202psid92m, clear
keep if mod(_n,2) == 0
keep id education hours
list in 1/4, clean
qui save merge2, replace

* Merge 1:1 two datasets with some observations and variables different
clear
use merge1
sort id
merge 1:1 id using merge2
sort id
list in 1/4, clean

* Not included - drop sort id
clear
use merge1
merge 1:1 id using merge2
list in 1/4, clean

* Not included - reverse the roles of master and using data sets
clear
use merge2
sort id
merge 1:1 id using merge1
sort id
list in 1/4, clean

* Append two datasets with some observations and variables different
clear
use merge1
append using merge2
sort id
list in 1/5, clean

************* 2.6: GRAPHICAL DISPLAY OF DATA

*!! sls do not include graphregion(margin(r+2)) in logs -- used bc we have a
*problem with graph export

use mus202psid92m, clear
twoway scatter lnearns hours, graphregion(margin(r+2))

* More advanced graphics command with two plots and with several options
graph twoway (scatter lnearns hours, msize(small))     ///
    (lfit lnearns hours, lwidth(medthick)),            ///
    graphregion(margin(r+2))                           ///
    title("Scatterplot and OLS fitted line")

use mus202psid92m, clear
label define edtype 1 "< high school" 2 "High school" 3 "Some college" 4 "College degree"
label values edcat edtype

* Box-and-whisker plot of single variable over several categories
graph box hours, over(edcat) scale(1.2) marker(1,msize(vsmall))   ///
    ytitle("Annual hours worked by education") yscale(titlegap(*5))

histogram lnearns

* Histogram with bin width and start value set
histogram lnearns, width(0.25) start(4.0) scale(1.2)

// Output not included in text but this provides the default bandwidth
kdensity lnearns

* Kernel density plot with bandwidth set and fitted normal density overlaid
kdensity lnearns, bwidth(0.12) normal n(4000) scale(1.2)

* Histogram and nonparametric kernel density estimate
histogram lnearns if lnearns > 0, width(0.25) kdensity       ///
    kdenopts(bwidth(0.2) lwidth(medthick))                   ///
    plotregion(style(none)) scale(1.2)                       ///
    title("Histogram and density for log earnings")          ///
    xtitle("Log annual earnings", size(medlarge)) xscale(titlegap(*5))  ///
    ytitle("Histogram and density", size(medlarge)) yscale(titlegap(*5))

* Simple two-way scatterplot
scatter lnearns hours

* Two-way scatterplot and quadratic regression curve with 95% ci for y|x
twoway (qfitci lnearns hours, stdf) (scatter lnearns hours, msize(small))

* Local constant with epan2 kernel and 95% confidence bands
use mus202psid92m, clear
lpoly lnearns hours, kernel(epan2) ci msize(tiny) lwidth(medthick) ///
    plotregion(style(none)) xtitle("Annual hours", size(medlarge)) ///
    title("Local constant smooth") scale(1.1)                      ///
    ytitle("Natural log of annual earnings", size(medlarge))       ///
    legend(pos(4) ring(0) col(1)) legend(size(small))              ///
    legend(label(1 "CI") label(2 "Actual data") label(3 "Local constant"))

* Scatterplot with lowess and local linear nonparametric regression
graph twoway (scatter lnearns hours, msize(tiny))                    ///
    (lpoly lnearns hours, kernel(epan2) degree(1) clstyle(p1) lwidth(thick) ///
    bwidth(500)) (lowess lnearns hours, clstyle(p2) lwidth(thick)),  ///
    plotregion(style(none)) title("Local linear and lowess fits")    ///
    xtitle("Annual hours", size(medlarge)) scale(1.1)                ///
    ytitle("Natural log of annual earnings", size(medlarge))         ///
    legend(pos(4) ring(0) col(1)) legend(size(small))                ///
    legend(label(1 "Actual data") label(2 "Local linear") label(3 "Lowess")) ///
    graphregion(margin(r+2))

* Multiple scatterplots
label variable age "Age"
label variable lnearns "Log earnings"
label variable hours "Annual hours"
graph matrix lnearns hours age, by(edcat) msize(tiny)

*** Erase files created by this program and not used elsewhere in the book

erase merge1.dta
erase merge2.dta
erase cleaneddata.dta
erase cleaneddata.csv
erase cleaneddata.asc

*********** END
