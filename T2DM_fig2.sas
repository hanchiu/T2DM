/*******************************************************************
Project:  Metformin use in the first trimester of pregnancy and risk of non-live birth and congenital malformations:
emulating a target trial using real-world data."

Authors: Yu-Han Chiu

Date created: March, 2024
codes run by Jennifer Yland, Helen Mogun

********************************************************************/

/*******************************************************************

                      Figure 2:
Comparison of estimated risk ratio of congenital malformations among live births using different analytic approaches

********************************************************************/



libname outfile '/PHShome/hm020/app4a/DM_Target_Trial_Emul/Temp';
/*libname outfile '/PHShome/jjy13/yuhan/';*/
libname MAX '/PHShome/hm020/app4a/DM_Target_Trial_Emul'; 


/************Cohort 1: overall cohort***************/
data cohort_1; set MAX.master; 
where infant_linkage=1 ;
age_lmpsq=age_lmp*age_lmp;
exposure0_num =SUBSTR(exposure0,1,1)*1;
exposure1_num =SUBSTR(exposure1,1,1)*1;
exposure2_num =SUBSTR(exposure2,1,1)*1;
exposure3_num =SUBSTR(exposure3,1,1)*1;
id=pregnancy_id;
if exposure3_num=1   then mi_exposed=1;
if exposure3_num=2   then mi_exposed=0;
if exposure1_num=4 then main_exposed=1;
if exposure1_num=3 then main_exposed=0;
/*add pre-pregnancy exposure*/
if exposure0_num_new=1 then pre_tx1=1; else pre_tx1=0;
if exposure0_num_new=2 then pre_tx2=1; else pre_tx2=0;
if exposure0_num_new=3 then pre_tx3=1; else pre_tx3=0;
if exposure0_num_new=4 then pre_tx4=1; else pre_tx4=0;
if exposure0_num_new=5 then pre_tx5=1; else pre_tx5=0;
if exposure0_num_new=6 then pre_tx6=1; else pre_tx6=0;
if exposure0_num_new=7 then pre_tx7=1; else pre_tx7=0;
if exposure0_num_new=8 then pre_tx8=1; else pre_tx8=0;
run;



/*restricting to metformin vs insulin*/
data cohort_11; set cohort_1;
where exposure3_num in (1,2);
run;

/*restricting to metformin+insulin vs insulin*/
data cohort_12; set cohort_1;
where exposure1_num in (3,4);
run;



/*Cohort 2: restrict to cohort2=1 (women with pregestational T2DM)*/
data cohort_2; set MAX.master; 
where infant_linkage=1 AND cohort2=1;
exposure0_num =SUBSTR(exposure0,1,1)*1;
exposure1_num =SUBSTR(exposure1,1,1)*1;
exposure2_num =SUBSTR(exposure2,1,1)*1;
exposure3_num =SUBSTR(exposure3,1,1)*1;
id=pregnancy_id;
age_lmpsq=age_lmp*age_lmp;
/*add pre-pregnancy exposure*/
if exposure0_num_new=1 then pre_tx1=1; else pre_tx1=0;
if exposure0_num_new=2 then pre_tx2=1; else pre_tx2=0;
if exposure0_num_new=3 then pre_tx3=1; else pre_tx3=0;
if exposure0_num_new=4 then pre_tx4=1; else pre_tx4=0;
if exposure0_num_new=5 then pre_tx5=1; else pre_tx5=0;
if exposure0_num_new=6 then pre_tx6=1; else pre_tx6=0;
if exposure0_num_new=7 then pre_tx7=1; else pre_tx7=0;
if exposure0_num_new=8 then pre_tx8=1; else pre_tx8=0;
if exposure3_num=1   then mi_exposed=1;
if exposure3_num=2   then mi_exposed=0;
if exposure1_num=4 then main_exposed=1;
if exposure1_num=3 then main_exposed=0;
if     age_lmp<30  then age_gp1=1; else age_gp1=0;
if 30=<age_lmp<35  then age_gp2=1; else age_gp2=0;
if 35=<age_lmp<40  then age_gp3=1; else age_gp3=0;
if     age_lmp>=40 then age_gp4=1; else age_gp4=0;
/*REVISION: exclude women if they use other ADM before pregnancy */
if exposure0_num_new in (5,6,7,8) then delete;
run;

data cohort_21; set cohort_2;
where exposure3_num in (1,2);
run;

data cohort_22; set cohort_2;
where exposure1_num in (3,4);
run;


data cohort_23; set cohort_2;
where exposure0_num in (3);
run;

 data cohort_24; set cohort_2;
where exposure0_num in (1);
run;


/*REVISION: exclude women if they use other ADM before pregnancy */
data cohort_221; set cohort_22;
if exposure0_num_new in (5,6,7,8) then delete;
run;

/*REVISION: metformin alone at baseline */
data cohort_226; set cohort_22;
if exposure0_num_new in (1,3,4,5,6,7,8) then delete;
run;



%macro ormodel(data=, outdata=, covar=);
proc genmod data= &data  descending ;
class metformin_exposed(ref="0") id;
model MAIN_Malform_Overall =  metformin_exposed  &covar /DIST=POISSON LINK=LOG;
lsmeans metformin_exposed /exp diff cl; 
title "&outdata";
repeated subject = id/ type = unstr;
run;
%mend;

/*#1 All pregnant women: Met-exposed vs unexposed: 1a, 1d*/
%ormodel(data=cohort_1, outdata=DM_01a_RR, covar=);
%ormodel(data=cohort_1, outdata=DM_01b_RR, covar=age_lmp age_lmp*age_lmp  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity);
%ormodel(data=cohort_1, outdata=DM_01c_RR, covar=PreGest_Diabetes_Type2  age_lmp age_lmp*age_lmp  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity);
%ormodel(data=cohort_1, outdata=DM_01d_RR, covar=pre_tx1 pre_tx2 pre_tx3 pre_tx4  PreGest_Diabetes_Type2  age_lmp age_lmp*age_lmp  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity);



/*#2 Type II DM: Met-exposed vs unexposed: 2a, 2c */
%ormodel(data=cohort_2, outdata=DM_02a_RR, covar=);
%ormodel(data=cohort_2, outdata=DM_02c_RR, covar=pre_tx1 pre_tx2 pre_tx3 pre_tx4  age_lmp age_lmp*age_lmp  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity);





%macro ormodel2(data=, outdata=, covar=);
proc genmod data= &data  descending ;
class mi_exposed(ref="0") id;
model MAIN_Malform_Overall =  mi_exposed  &covar /DIST=POISSON LINK=LOG;
lsmeans mi_exposed /exp diff cl; 
title "&outdata";
repeated subject = id/ type = unstr;
run;
%mend;



/*#3 Type II DM: Metformin vs insulin): 3a, 3c*/
%ormodel2(data=cohort_21, outdata=DM_03a_RR, covar=);
%ormodel2(data=cohort_21, outdata=DM_03c_RR, covar=pre_tx1 pre_tx2 pre_tx3 pre_tx4 age_lmp age_lmp*age_lmp  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity);




%macro ormodel3(data=, outdata=, covar=);
proc genmod data= &data  descending ;
class main_exposed(ref="0") id;
model MAIN_Malform_Overall =  main_exposed  &covar /DIST=POISSON LINK=LOG;
lsmeans main_exposed /exp diff cl; 
title "&outdata";
repeated subject = id/ type = unstr;
run;
%mend;



/*#4 Type II DM: Metformin+insulin vs insulin): 4a, 4c*/
%ormodel3(data=cohort_22, outdata=DM_04a_RR, covar=);
%ormodel3(data=cohort_22, outdata=DM_04c_RR, covar=pre_tx1 pre_tx2 pre_tx3 pre_tx4  age_lmp age_lmp*age_lmp  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity);




/*#5 Type II DM + excluding women who used other antidiabetics before LMP: Metformin+insulin vs insulin): 5a, 5b*/
%ormodel3(data=cohort_221, outdata=DM_041a_RR, covar=);
%ormodel3(data=cohort_221, outdata=DM_041_RR, covar=pre_tx1 pre_tx2 pre_tx3  age_lmp age_lmp*age_lmp  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity);


/*#6 Type II DM + pre-LMP metformin monotherapy: Metformin+insulin vs insulin): 6a, 6b*/
%ormodel3(data=cohort_226, outdata=DM_046a_RR, covar=);
%ormodel3(data=cohort_226, outdata=DM_046_RR, covar=          age_lmp age_lmp*age_lmp  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity);







