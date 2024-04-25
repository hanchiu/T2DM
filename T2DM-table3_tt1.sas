/*******************************************************************
Project:  Metformin use in the first trimester of pregnancy and risk of non-live birth and congenital malformations:
emulating a target trial using real-world data."

Authors: Yu-Han Chiu

Date created: March, 2023
codes run by Jennifer Yland, Helen Mogun

********************************************************************/

/*******************************************************************

       Table 3 (target trial with stricter eligibility)

Estimated risk of non-live births and non-chromosomal major congenital malformations 
under insulin monotherapy or insulin plus metformin therapy during the first trimester 
among women with type 2 diabetes and used metformin monotherapy before LMP.

 
Strict eligibility: Include T2DM women with metformin only within 180 days before LMP (requiring only 1 prescription)
                   
********************************************************************/


libname MAX '/PHShome/hm020/app4a/DM_Target_Trial_Emul'; 


%include '/PHShome/jjy13/yuhan/step2.1_bootstrap12n.sas';


data master1; set MAX.master;
where   cohort2=1; /*restrict to women with pre-gestational T2DM*/
if exposure0_num_new in (1,3,4,5,6,7,8) then delete; /*restrict to women with metformin only before pregnancy*/
/*****************************************
exposure0_num_new (exposure between (LMP-180) and LMP-1)
1. no treatment
2. metformin only
3. insulin only
4. metformin and insulin
5. any other antidiabetic only
6. metformin in combination (other than insulin)
7. insulin in combination (other than metformin)
8. Any other combination.
**********************************************/
ID=pregnancy_id;
age_lmpsq=age_lmp*age_lmp;
exposure0_num =SUBSTR(exposure0,1,1)*1;
exposure1_num =SUBSTR(exposure1,1,1)*1;
exposure2_num =SUBSTR(exposure2,1,1)*1;
exposure3_num =SUBSTR(exposure3,1,1)*1;
gestage_w=round(Gestage/7,1);
if cohort in (0,1) then live=1; else live=0; /*live birth*/
insulin_w=round(GA_first_Insulin/7,1);       /*week at the first insulin prescription*/
metformin_w=round(GA_first_metformin/7,1);   /*week at the first metformin prescription*/
othermeds_w=round(GA_first_othermeds/7,1);   /*week at the other antidiabetics prescription*/
if exposure1_num=4 then main_exposed=1; /*use of metfomin + insulin during the first trimester*/
if exposure1_num=3 then main_exposed=0; /*use of  insulin monotherapy during the first trimester*/
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

/insulin monotherapy group*/
proc sort data=longdata; by id time;
data long_txt0; set longdata;
C=0; L=0; /*C: indicator of censoring if they deviate the protocol; L: indicator of non-live birth*/
if metformin=1     then C=1; 
if otherdrug=1     then C=1;
if time>=gestage_w then L=1;
if live=0 AND gestage_w>10 AND time=11 then L=1;
if  exposure1_num=1 AND time in (11) then C=1;
lagL=lag(L);
if time=1 then lagL=0;
totL=lagL+L;
if totL>1 then delete;
if C=1 then delete;   /*remove women who is artificially censored*/
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




**********Estimating the probability of non-live birth ********;

******************;
proc logistic data= long_txt0  descending outest=LBoutc0  noprint;
model L =time8 time9 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications 
  obesity
;
freq numberhits;
run;


proc logistic data= wide  descending outest=maloutc0  noprint;
where exposure1_num in (3) AND infant_linkage=1;
model MAIN_Malform_Overall =  age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications 
 obesity
;
freq numberhits;
run;


proc logistic data= wide  descending outest=maloutc1  noprint;
where exposure1_num in (4) AND infant_linkage=1;
model MAIN_Malform_Overall = age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications 
  obesity
;
freq numberhits;
run;





data outc0;
set LBoutc0;
if _type_='PARMS';
_sample_ = 0;
array avar intercept  time8 time9 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications   obesity
;
array abeta boutc00-boutc11;
do i=1 to dim(avar);
abeta(i)=avar(i);
end;
keep _sample_ boutc00-boutc11;
run;


data mal_outc0;
set maloutc0;
if _type_='PARMS';
_sample_ = 0;
array avar intercept   age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications  obesity
;
array abeta bmal00-bmal08 ;
do i=1 to dim(avar);
abeta(i)=avar(i);
end;
keep _sample_ bmal00-bmal08;
run;


data mal_outc1;
set maloutc1;
if _type_='PARMS';
_sample_ = 0;
array avar intercept   age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity
;
array abeta bmal00-bmal08 ;
do i=1 to dim(avar);
abeta(i)=avar(i);
end;
keep _sample_ bmal00-bmal08;
run;


data simul0 (keep=ID  _sample_  time time1 time2 time3 time4 time5 time6 time7 time8 time9 time10 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity
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
array aboutc boutc00-boutc11;
array asS    uno  time8 time9 time11  age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity
;
array abmal  bmal00-bmal08;
array asY    uno age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity
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
model L =time8 time9 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity;
freq numberhits;
run;


data outc1;
set LBoutc1;
if _type_='PARMS';
_sample_ = 0;
array avar intercept  time8 time9 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity
;
array abeta boutc00-boutc11;
do i=1 to dim(avar);
abeta(i)=avar(i);
end;
keep _sample_ boutc00-boutc11;
run;

data simul1 (keep=ID  _sample_  time time1 time2 time3 time4 time5 time6 time7 time8 time9 time10 time11 age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity
 ); 
set long_txt1 ;
 _sample_ = 0;
 where time=1;
run;


data sim1; 
merge outc1 simul1 mal_outc1 end=_end_;
by _sample_;
uno=1;
main_exposed=1;
array aboutc boutc00-boutc11;
array asS    uno  time8 time9 time11  age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity
;
array abmal  bmal00-bmal08;
array asY    uno  age_lmp age_lmpsq  hbp hyperlipidemia hyperglycemia hypoglycemia dm_complications obesity
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
dataout=tt1,
coeff1=live0,      coeff2=live1, /*live birth*/
coeff3=totY0,      coeff4=totY1, /*total effect of malformation*/
coeff5=cY0,        coeff6=cY1,    /*direct effect of malformations among live birth*/
coeff7=RD_live,    coeff8=RR_live,    /*RD. RR for live birth*/
coeff9=RD_mal_c,   coeff10=RR_mal_c, /*RD. RR for malformation among live birth*/
coeff11=RD_mal_t,  coeff12=RR_mal_t  /*RD, RR for total effects*/
);



