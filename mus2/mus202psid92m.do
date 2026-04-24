* mus202psid92m.do
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
