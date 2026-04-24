* ML_2022_part6.do   April 2022
* based on machinelearn2022_part2.do  

capture log close

log using ML_2022_part6.txt, text replace

********** OVERVIEW OF ML_2022_part6.do **********

* To run you need files
*    mus203mepsmedexp.dta
* in your directory

* Stata user-written command
*    svmachines (for support vector machines) 
* are used

*  1. CLASSIFICATION
*       LOGIT
*       K NEAREST NEIGHBORS
*       LINEAR DISCRIMINANT ANALYSIS
*       QUADRATIC DISCRIMINANT ANALYSIS
*       SUPPORT VECTOR MACHINES
*  2. UNSUPERVISED LEARNING
*       CLUSTER ANALYSIS
 
************ SETUP ***********

set more off
* version 17
clear all
set linesize 82
set scheme s1mono  /* Graphics scheme */

************ DATA DESCRIPTION ***********

**************** CATEGORICAL DATA

* Data for 65-90 year olds on supplementary insurance indicator and regressors
use mus203mepsmedexp.dta, clear
global xlist income educyr age female white hisp marry ///
   totchr phylim actlim hvgg 
describe suppins $xlist

* Summary statistics   
summarize suppins $xlist  

* logit model
logit suppins $xlist, nolog

* Classification table
estat classification

* Classification table manually
predict ph_logit
generate yh_logit = ph_logit >= 0.5
generate err_logit = (suppins==0 & yh_logit==1) | (suppins==1 & yh_logit==0)
summarize suppins ph_logit yh_logit err_logit
tabulate suppins yh_logit 

* K-nearest neighbors 
discrim knn $xlist, group(suppins) k(11) notable
predict yh_knn  
estat classtable, nototals nopercents looclass

* K-nn classification table with leave-one out cross validation not as good
estat classtable, nototals nopercents  // without LOOCV

* Linear discriminant analysis
discrim lda $xlist, group(suppins) notable
predict yh_lda
estat classtable, nototals nopercents

* Quadratic discriminant analysis
discrim qda $xlist, group(suppins) notable
predict yh_qda
estat classtable, nototals nopercents

* Support vector machines - need y to be byte not float and matsize > n
* set matsize 3200     // Newer versions of Stata do this automatically
global xlistshort income educyr age female marry totchr
generate byte ins = suppins
svmachines ins income
svmachines ins $xlist
predict yh_svm
tabulate ins yh_svm

* Compare various in-sample predictions
correlate suppins yh_logit yh_knn yh_lda yh_qda yh_svm

**************** CLUSTER ANALYSIS

* Cluster analysis
* Generated data: y = 1 + 1*x1 + 1*x2 + f(z) + u where f(z) = z + z^2
clear
set obs 200
set seed 10101
generate x1 = rnormal()
generate x2 = rnormal() + 0.5*x1
generate z = rnormal() + 0.5*x1
generate zsq = z^2
generate y = 1 + x1 + x2 + z + zsq + 2*rnormal()
summarize

* k-means clustering with defaults and three clusters

* use machlearn_part2_spline.dta, replace

* Generated data: y = 1 + 1*x1 + 1*x2 + f(z) + u where f(z) = z + z^2
clear
set obs 200
set seed 10101
generate x1 = rnormal()
generate x2 = rnormal() + 0.5*x1
generate z = rnormal() + 0.5*x1
generate zsq = z^2
generate y = 1 + x1 + x2 + z + zsq + 2*rnormal()
summarize

// Not included - estimate same model as DGP
reg y x1 x2 z zsq

graph matrix x1 x2 z     // matrix plot of the three variables
cluster kmeans x1 x2 z, k(3) name(myclusters)
tabstat x1 x2 z, by(myclusters) stat(mean)
* graph export machlearn_us214fig7_locallinear.eps, replace

**** CLOSE OUTPUT  

                                  
