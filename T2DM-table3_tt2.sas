/*******************************************************************
Project:  Metformin use in the first trimester of pregnancy and risk of non-live birth and congenital malformations:
emulating a target trial using real-world data."

Authors: Yu-Han Chiu

Date created: March, 2024
codes run by Jennifer Yland, Helen Mogun

********************************************************************/

/*******************************************************************

       Table 3

Estimated risk of non-live births and non-chromosomal major congenital malformations 
under insulin monotherapy or insulin plus metformin therapy during the first trimester 
among women with type 2 diabetes (excluding women who used non-metformin and non-insulin antidiabetics before LMP)

Broad eligibility: Include T2DM women (excluding women who used any other antidiabetics within 180 days before LMP)            
********************************************************************/


libname MAX '/PHShome/hm020/app4a/DM_Target_Trial_Emul'; /*change the location of the path*/

%include '/PHShome/hm020/app4a/DM_Target_Trial_Emul/Temp/step2.1_bootstrap12n.sas';

/*%include '/PHShome/jjy13/yuhan/step2.1_bootstrap12n.sas';*/

data master1; set MAX.master;
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
if exposure0_num_new=1 then pre_tx1=1; else pre_tx1=0;
if exposure0_num_new=2 then pre_tx2=1; else pre_tx2=0;
if exposure0_num_new=3 then pre_tx3=1; else pre_tx3=0;
if exposure0_num_new=4 then pre_tx4=1; else pre_tx4=0;
if exposure0_num_new=5 then pre_tx5=1; else pre_tx5=0;
if exposure0_num_new=6 then pre_tx6=1; else pre_tx6=0;
if exposure0_num_new=7 then pre_tx7=1; else pre_tx7=0;
if exposure0_num_new=8 then pre_tx8=1; else pre_tx8=0;

/*REVISION: exclude women if they use other ADM before pregnancy */
if exposure0_num_new in (5,6,7,8) then delete;
/*
excluding women with use of non-metformin and non-insulin agent before pregnancy.

exposure0_num_new
1. no treatment
2. metformin only
3. insulin only
4. metformin and insulin
5. any other antidiabetic only
6. metformin in combination (other than insulin)
7. insulin in combination (other than metformin)
8. Any other combination.
*/
run;



proc sort data=master1; by id;



%macro computemean1;
data wide; set &data;
array KK time1-time11;
do i=1 to dim(KK);
KK{i}=i;
end;
insulin_week=round(GA_first_Insulin/7,1); 
metformin_week=round(GA_first_metformin/7,1); 
othermeds_week=round(GA_first_othermeds/7,1); 
if othermeds_week=. then othermeds_week=9999;
run;


data long; set wide;
time=time1; output;
time=time2; output;
time=time3; output;
time=time4; output;
time=time5; output;
time=time6; output;
time=time7; output;
time=time8; output;
time=time9; output;
time=time10; output;
time=time11; output;
run;
data long; set long;
drop time1-time11;
run;



data longdata; set long;
insulin=0; metformin=0;
if exposure1_num in (3,4) then do;
if time>=insulin_week then insulin=1;
end;
if exposure1_num in (2,4) then do;
if time>=metformin_week then metformin=1;
end;

if exposure1_num in (5,6) then do;
otherdrug=0;
if time>=othermeds_week then otherdrug=1;
end;
run;

/*exclusive insulin*/
proc sort data=longdata; by id time;
data long_txt0; set longdata;
C=0; L=0; 
if metformin=1     then C=1;
if otherdrug=1     then C=1;
if time>=gestage_w then L=1;
if live=0 AND gestage_w>10 AND time=11 then L=1;
if  exposure1_num=1 AND time in (11) then C=1;
lagL=lag(L);
if time=1 then lagL=0;
totL=lagL+L;
if totL>1 then delete;
if C=1 then delete;
if time=1 then time1=1; else time1=0;
if time=2 then time2=1; else time2=0;
if time=3 then time3=1; else time3=0;
if time=4 then time4=1; else time4=0;
if time=5 then time5=1; else time5=0;
if time=6 then time6=1; else time6=0;
if time=7 then time7=1; else time7=0;
if time=8 then time8=1; else time8=0;
if time=9 then time9=1; else time9=0;
if time=10 then time10=1; else time10=0;
if time=11 then time11=1; else time11=0;
run;




******************;
proc logistic data= long_txt0  descending outest=LBoutc0  noprint;
model L =time8 time9 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications 
pre_tx1 pre_tx2 pre_tx3   obesity
;
freq numberhits;
run;


proc logistic data= wide  descending outest=maloutc0  noprint;
where exposure1_num in (3,4) AND infant_linkage=1;
model MAIN_Malform_Overall = main_exposed  age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications 
pre_tx1 pre_tx2 pre_tx3   obesity
;
freq numberhits;
run;





data outc0;
set LBoutc0;
if _type_='PARMS';
_sample_ = 0;
array avar intercept  time8 time9 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications pre_tx1 pre_tx2 pre_tx3   obesity
;
array abeta boutc00-boutc14;
do i=1 to dim(avar);
abeta(i)=avar(i);
end;
keep _sample_ boutc00-boutc14;
run;


data mal_outc0;
set maloutc0;
if _type_='PARMS';
_sample_ = 0;
array avar intercept main_exposed  age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications pre_tx1 pre_tx2 pre_tx3  obesity
;
array abeta bmal00-bmal12 ;
do i=1 to dim(avar);
abeta(i)=avar(i);
end;
keep _sample_ bmal00-bmal12;
run;



data simul0 (keep=ID  _sample_  time time1 time2 time3 time4 time5 time6 time7 time8 time9 time10 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications pre_tx1 pre_tx2 pre_tx3 pre_tx4  obesity
 ); 
set long_txt0 ;
 _sample_ = 0;
 where time=1;
run;


data sim0; 
merge outc0 simul0 mal_outc0 end=_end_;
by _sample_;
uno=1;
main_exposed=0;
array aboutc boutc00-boutc14;
array asS    uno  time8 time9 time11  age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications pre_tx1 pre_tx2 pre_tx3  obesity
;
array abmal  bmal00-bmal12;
array asY    uno main_exposed age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications pre_tx1 pre_tx2 pre_tx3  obesity
;

 mY = 0 ;
do i=1 to dim(abmal);
  mY =sum(mY,abmal(i)*asY(i));
  end;
pY = 1/(1+exp(-mY));

array sS {1:11 } ;
do time=1 to 11;
if time=1 then time1=1; else time1=0;
if time=2 then time2=1; else time2=0;
if time=3 then time3=1; else time3=0;
if time=4 then time4=1; else time4=0;
if time=5 then time5=1; else time5=0;
if time=6 then time6=1; else time6=0;
if time=7 then time7=1; else time7=0;
if time=8 then time8=1; else time8=0;
if time=9 then time9=1; else time9=0;
if time=10 then time10=1; else time10=0;
if time=11 then time11=1; else time11=0;
  mS = 0 ;
  do i=1 to dim(aboutc);
  mS =sum(mS,aboutc(i)*asS(i));
  end;
pS = 1/(1+exp(-mS));
sS[time] = pS ;
end;
cumsurv= 1;
cuminc =0;
array cumincr{1:11} ;
array cumsurvr{1:11} ;
do t=1 to 11;
inc = cumsurv *sS[t];
cuminc = cuminc + inc;
cumincr[t] = cuminc ;
surv=(1.0 - sS[t]);
cumsurv=cumsurv*surv;
cumsurvr[t] = cumsurv ;
end;
drop t inc  ;
totY=cumsurv*pY;
run;

proc means data=sim0 noprint;
output out=stats_holder0
mean (pY  totY cumsurv  )
=     cY0  totY0  live0 
;
run;


data long_txt1; set longdata;
C=0; L=0; 
if otherdrug =1     then C=1;
if time>=gestage_w  then L=1;
if live=0 AND gestage_w>10 AND time=11 then L=1;
if  exposure1_num=2 AND time in (11) then C=1;
if  exposure1_num=3 AND time in (11) then C=1;
if C=1 then delete;
lagL=lag(L);
if time=1 then lagL=0;
totL=lagL+L;
if totL>1 then delete;
if time=1 then time1=1; else time1=0;
if time=2 then time2=1; else time2=0;
if time=3 then time3=1; else time3=0;
if time=4 then time4=1; else time4=0;
if time=5 then time5=1; else time5=0;
if time=6 then time6=1; else time6=0;
if time=7 then time7=1; else time7=0;
if time=8 then time8=1; else time8=0;
if time=9 then time9=1; else time9=0;
if time=10 then time10=1; else time10=0;
if time=11 then time11=1; else time11=0;
run;





proc logistic data= long_txt1  descending outest=LBoutc1 noprint;
model L =time8 time9 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications pre_tx1 pre_tx2 pre_tx3  obesity;
freq numberhits;
run;


data outc1;
set LBoutc1;
if _type_='PARMS';
_sample_ = 0;
array avar intercept  time8 time9 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications pre_tx1 pre_tx2 pre_tx3  obesity
;
array abeta boutc00-boutc14;
do i=1 to dim(avar);
abeta(i)=avar(i);
end;
keep _sample_ boutc00-boutc14;
run;

data simul1 (keep=ID  _sample_  time time1 time2 time3 time4 time5 time6 time7 time8 time9 time10 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications pre_tx1 pre_tx2 pre_tx3  obesity
 ); 
set long_txt1 ;
 _sample_ = 0;
 where time=1;
run;


data sim1; 
merge outc1 simul1 mal_outc0 end=_end_;
by _sample_;
uno=1;
main_exposed=1;
array aboutc boutc00-boutc14;
array asS    uno  time8 time9 time11  age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications pre_tx1 pre_tx2 pre_tx3  obesity
;
array abmal  bmal00-bmal12;
array asY    uno  main_exposed age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications pre_tx1 pre_tx2 pre_tx3  obesity
;

 mY = 0 ;
do i=1 to dim(abmal);
  mY =sum(mY,abmal(i)*asY(i));
  end;
pY = 1/(1+exp(-mY));

array sS {1:11 } ;
do time=1 to 11;
if time=1 then time1=1; else time1=0;
if time=2 then time2=1; else time2=0;
if time=3 then time3=1; else time3=0;
if time=4 then time4=1; else time4=0;
if time=5 then time5=1; else time5=0;
if time=6 then time6=1; else time6=0;
if time=7 then time7=1; else time7=0;
if time=8 then time8=1; else time8=0;
if time=9 then time9=1; else time9=0;
if time=10 then time10=1; else time10=0;
if time=11 then time11=1; else time11=0;
  mS = 0 ;
  do i=1 to dim(aboutc);
  mS =sum(mS,aboutc(i)*asS(i));
  end;
pS = 1/(1+exp(-mS));
sS[time] = pS ;
end;
cumsurv= 1;
cuminc =0;
array cumincr{1:11} ;
array cumsurvr{1:11} ;
do t=1 to 11;
inc = cumsurv *sS[t];
cuminc = cuminc + inc;
cumincr[t] = cuminc ;
surv=(1.0 - sS[t]);
cumsurv=cumsurv*surv;
cumsurvr[t] = cumsurv ;
end;
drop t inc  ;
totY=cumsurv*pY;
run;


proc means data=sim1 noprint;
output out=stats_holder1
mean (pY  totY cumsurv  )
=    cY1  totY1  live1 
;
run;


data final;
merge stats_holder1 stats_holder0;
by _type_;
RR_live=live1/live0;
RD_live=live1-live0;
RR_mal_t=totY1/totY0;
RD_mal_t=totY1-totY0;
RR_mal_c=cY1/cY0;
RD_mal_c=cY1-cY0;
OR_mal_c=(cY1/(1-cY1))/(cY0/(1-cY0));
drop _type_;
run;

%mend;


%bootstrap12n(data=master1,nsamples=500,compmacro=computemean1,resultdatatemp=final, 
dataout=tt2,
coeff1=live0,      coeff2=live1, /*live birth*/
coeff3=totY0,      coeff4=totY1, /*total effect of malformation*/
coeff5=cY0,        coeff6=cY1,    /*direct effect of malformations among live birth*/
coeff7=RD_live,    coeff8=RR_live,    /*RD. RR for live birth*/
coeff9=RD_mal_c,   coeff10=RR_mal_c, /*RD. RR for malformation among live birth*/
coeff11=RD_mal_t,  coeff12=RR_mal_t  /*RD, RR for total effects*/
);






