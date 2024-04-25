
/*******************************************************************
Project:  Metformin use in the first trimester of pregnancy and risk of non-live birth and congenital malformations:
emulating a target trial using real-world data."

Authors: Yu-Han Chiu

Date created: April, 2024
codes run by Jennifer Yland, Helen Mogun

********************************************************************/

/*******************************************************************

Table 2: Baseline characteristics of women with type 2 diabetes who were included for emulating a target trial of insulin monotherapy vs. insulin plus metformin at the end of the grace period, 
                    
Figure 1: Flow chart of eligible individuals of the target trial restricted to metformin users before pregnancy

********************************************************************/


libname MAX '/PHShome/hm020/app4a/DM_Target_Trial_Emul'; 
libname outfile '/PHShome/jjy13/yuhan/';


/*flow chart: the target trial*/

data flowchart; set MAX.master;
where   cohort2=1;
if exposure0_num_new in (1,3,4,5,6,7,8) then delete;
ID=pregnancy_id;
age_lmpsq=age_lmp*age_lmp;
exposure0_num =SUBSTR(exposure0,1,1)*1;
exposure1_num =SUBSTR(exposure1,1,1)*1;
exposure2_num =SUBSTR(exposure2,1,1)*1;
exposure3_num =SUBSTR(exposure3,1,1)*1;
gestage_w=round(Gestage/7,1);
if cohort in (0,1) then live=1; else live=0;
insulin_w=round(GA_first_Insulin/7,1); 
metformin_w=round(GA_first_metformin/7,1); 
othermeds_w=round(GA_first_othermeds/7,1); 
if exposure1_num=4 then main_exposed=1;
if exposure1_num=3 then main_exposed=0;
run;



proc freq data=flowchart; where exposure1_num=1; tables Gestage; run;
proc freq data=flowchart; where exposure1_num=2; tables Gestage; run;
proc freq data=flowchart; where exposure1_num=3; tables Gestage; run;
proc freq data=flowchart; where exposure1_num=4; tables Gestage; run;
proc freq data=flowchart; where exposure1_num in (5,6); tables Gestage; run;


proc freq data=flowchart; where exposure1_num=1; tables live*Gestage; run;
proc freq data=flowchart; where exposure1_num=2; tables live*Gestage; run;
proc freq data=flowchart; where exposure1_num=3; tables live*Gestage; run;
proc freq data=flowchart; where exposure1_num=4; tables live*Gestage; run;
proc freq data=flowchart; where exposure1_num in (5,6); tables live*Gestage; run;


data  flowchart_grace0; set flowchart;
censor=0;
if exposure1_num in (2,4,5,6) then censor=1;
if exposure1_num in (1) AND GestAge>90 then censor=1;
remain=0;
if exposure1_num in (3) AND GestAge>90 then remain=1; 
/*those remained pregnancy+ use insulin in 1st trimester)*/
run;

data  flowchart_grace1; set flowchart;
censor=0;
if exposure1_num in (5,6) then censor=1;
if exposure1_num in (1,2,3) AND GestAge>90  then censor=1;

remain=0;
if exposure1_num in (4) AND GestAge>90 then remain=1; 
/*those remained pregnancy+ use insulin+ metformin in 1st trimester remain=1)*/
run;


/*********************Table 2 (target trial with strict eligibility)**********************************************/
proc freq data=flowchart_grace0; title "main: grace0 (uncensored or no loess by 90 days); n=670"; 
where remain=1; tables 
dm_complications hbp hyperlipidemia hyperglycemia hypoglycemia  smoking pcos 
Race_WHITE  Race_BLACK_OR_AFRICAN_AMERICAN Race_HISPANIC_OR_LATINO
Race_ASIAN       Race_UNKNOWN_OTHER   ; run;

proc means mean std data=flowchart_grace0;  where remain=1; var age_lmp; run;



proc freq data=flowchart_grace1; title "main: grace1 (uncensored or no loess by 90 days); n=1286"; 
where remain=1; tables 
dm_complications hbp hyperlipidemia hyperglycemia hypoglycemia  smoking pcos 
Race_WHITE  Race_BLACK_OR_AFRICAN_AMERICAN Race_HISPANIC_OR_LATINO
Race_ASIAN       Race_UNKNOWN_OTHER   ; run;

proc means mean std data=flowchart_grace1;  where remain=1; var age_lmp; run;


/*do a chis square test for baseline characteristics by treatment groups*/
data flowchart2_grace0; set flowchart2_grace0;
where remain=1;
tx=0;
run;
data flowchart2_grace1; set flowchart2_grace1;
where remain=1;
tx=1;
run;

data compare; set 
flowchart2_grace1 flowchart_grace0;
run;

proc freq data=compare; title "main: grace1 (uncensored or no loess by 90 days); n=1660"; 
tables 
(dm_complications hbp hyperlipidemia hyperglycemia hypoglycemia  smoking pcos 
Race_WHITE  Race_BLACK_OR_AFRICAN_AMERICAN Race_HISPANIC_OR_LATINO
Race_ASIAN       Race_UNKNOWN_OTHER   obesity age_group exposure0_num_new live  MAIN_Malform_Overall cohort3)*tx/chisq ; run;






/********************************************************
Supplemental Figure 1
flow chart2: the target trial with broad eligibility*/


data flowchart2; set MAX.master;
where   cohort2=1;
ID=pregnancy_id;
age_lmpsq=age_lmp*age_lmp;
exposure0_num =SUBSTR(exposure0,1,1)*1;
exposure1_num =SUBSTR(exposure1,1,1)*1;
exposure2_num =SUBSTR(exposure2,1,1)*1;
exposure3_num =SUBSTR(exposure3,1,1)*1;
gestage_w=round(Gestage/7,1);
if cohort in (0,1) then live=1; else live=0;
insulin_w=round(GA_first_Insulin/7,1); 
metformin_w=round(GA_first_metformin/7,1); 
othermeds_w=round(GA_first_othermeds/7,1); 
if exposure1_num=4 then main_exposed=1;
if exposure1_num=3 then main_exposed=0;
/*add pre-pregnancy exposure*/
if exposure0_num=1 then pre_tx1=1; else pre_tx1=0;
if exposure0_num=2 then pre_tx2=1; else pre_tx2=0;
if exposure0_num=3 then pre_tx3=1; else pre_tx3=0;
if exposure0_num=4 then pre_tx4=1; else pre_tx4=0;
if exposure0_num in (5,6) then pre_tx5=1; else pre_tx5=0;
/*REVISION: exclude women if they use other ADM before pregnancy */
if exposure0_num_new in (5,6,7,8) then delete;
run;




proc freq data=flowchart2; where exposure1_num=1; tables Gestage; run;
proc freq data=flowchart2; where exposure1_num=2; tables Gestage; run;
proc freq data=flowchart2; where exposure1_num=3; tables Gestage; run;
proc freq data=flowchart2; where exposure1_num=4; tables Gestage; run;
proc freq data=flowchart2; where exposure1_num in (5,6); tables Gestage; run;

proc freq data=flowchart2; where exposure1_num=1; tables live*Gestage; run;
proc freq data=flowchart2; where exposure1_num=2; tables live*Gestage; run;
proc freq data=flowchart2; where exposure1_num=3; tables live*Gestage; run;
proc freq data=flowchart2; where exposure1_num=4; tables live*Gestage; run;
proc freq data=flowchart2; where exposure1_num in (5,6); tables live*Gestage; run;


proc sort data=flowchart2; by id;

data  flowchart2_grace0; set flowchart2;
censor=0;
if exposure1_num in (2,4,5,6) then censor=1;
if exposure1_num in (1) AND GestAge>90 then censor=1;
remain=0;
if exposure1_num in (3) AND GestAge>90 then remain=1; 
/*those remained pregnancy+ use insulin in 1st trimester
remain=1: )*/
run;


data  flowchart2_grace1; set flowchart2;
censor=0;
if exposure1_num in (5,6) then censor=1;
if exposure1_num in (1,2,3) AND GestAge>90  then censor=1;

remain=0;
if exposure1_num in (4) AND GestAge>90 then remain=1; 
/*those remained pregnancy+ use insulin+ metformin in 1st trimester
remain=1:)*/
run;




/***************Table 2 (target trial with broad eligibility)*********************************************/
proc freq data=flowchart2_grace0; title "main: grace0 (uncensored or no loess by 90 days); n=5648"; 
where remain=1; tables 
dm_complications hbp hyperlipidemia hyperglycemia hypoglycemia  smoking pcos 
Race_WHITE  Race_BLACK_OR_AFRICAN_AMERICAN Race_HISPANIC_OR_LATINO
Race_ASIAN       Race_UNKNOWN_OTHER  obesity  exposure0_num_new live MAIN_Malform_Overall  ; run;

proc means mean std data=flowchart2_grace0;  where remain=1; var age_lmp; run;

/*metformin + insulin group*/
proc freq data=flowchart2_grace1; tables censor remain; run; 


proc means mean std data=flowchart2_grace1;  where censor=0; var age_lmp; run;*/

proc freq data=flowchart2_grace1; title "main: grace1 (uncensored or no loess by 90 days); n=2880"; 
where remain=1; tables 
dm_complications hbp hyperlipidemia hyperglycemia hypoglycemia  smoking pcos 
Race_WHITE  Race_BLACK_OR_AFRICAN_AMERICAN Race_HISPANIC_OR_LATINO
Race_ASIAN       Race_UNKNOWN_OTHER   obesity  exposure0_num_new live  MAIN_Malform_Overall ; run;

proc means mean std data=flowchart2_grace1;  where remain=1; var age_lmp; run;







