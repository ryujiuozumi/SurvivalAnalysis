*======================================================================*;
/*

大橋靖雄・浜田知久馬・魚住龍史 (2022)
生存時間解析 第2版 --SASによる生物統計
東京大学出版会

第3章　コックス回帰とその応用 (PHREGプロシジャ)

Program : Chap3.sas
Made by : 魚住龍史 (e-mail : uozumi at kuhp.kyoto-u.ac.jp)
Date : 2022年8月

*/
*======================================================================*;

*==================================================*;
* 表3.1.1 ;
*==================================================*;

proc phreg data = Gehan; model Week*Remiss(0) = Drug / ties = exact;
proc phreg data = Gehan; model Week*Remiss(0) = Drug / ties = breslow;
proc phreg data = Gehan; model Week*Remiss(0) = Drug / ties = efron;
proc phreg data = Gehan; model Week*Remiss(0) = Drug / ties = discrete;
run;

*==================================================*;
* プログラム3.1.1 ;
*==================================================*;

data Gast;
do Center = 1 to 3; do Group = 0 to 1; do i = 1 to 3;
  input Time Censor w @@;
  if Center = 2 then C1 = 1; else C1 = 0;
  if Center = 3 then C2 = 1; else C2 = 0;
output; end; end; end;
cards;
1 1 15	3 1 17	3 0 7			1 1 15	3 1 17	3 0 2
1 1 12	3 1 13	3 0 15		1 1 17	3 1 17	3 0 10
1 1 3		3 1 17	3 0 18		1 1 7		3 1 17	3 0 16
;
proc phreg data = Gast;
  model Time*Censor(0) = Group C1 C2 / ties = exact; freq w;
proc phreg data = Gast;
  model Time*Censor(0) = Group C1 C2 / ties = breslow; freq w;
proc phreg data = Gast;
  model Time*Censor(0) = Group C1 C2 / ties = efron; freq w;
proc phreg data = Gast;
  model Time*Censor(0) = Group C1 C2 / ties = discrete; freq w;
run;

*==================================================*;
* プログラム3.3.1 ;
*==================================================*;

proc phreg data = Pcancer simple;
  model Time*Censor(0) = Age Sex Treat / rl;
run;

*==================================================*;
* プログラム3.4.1 ;
*==================================================*;

data Dummy; do Drug = 0 to 1; output; end; run;
proc phreg data = Gehan;
  model Week*Remiss(0) = Drug;
  baseline out = Out covariates = Dummy survival = S loglogs = LLS / nomean method = pl;
run;
proc print data = Out; by Drug; run;

*==================================================*;
* 図3.4.1 ;
*==================================================*;

data Out; set Out;
  LogWeek = log(Week);
run;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc sgplot data = Out noautolegend;
  step x = Week y = S / group = Drug name = "Survival";
  yaxis values = (0 to 1 by 0.1) label = "S"; 
  xaxis values = (0 to 30 by 10) label = "Week";
  keylegend "Survival" / noborder;
  format Drug Drugf.;
run;
proc sgplot data = Out noautolegend;
  step x = LogWeek y = LLS / group = Drug name = "Survival";
  yaxis values = (-5 to 2 by 1) label = "LLS"; 
  xaxis values = (0 to 3.2 by 0.2) label = "log(Week)";
  keylegend "Survival" / noborder;
  format Drug Drugf.;
run;
ods graphics off;
ods html close;

*==================================================*;
* プログラム3.4.2 ;
*==================================================*;

data Dummy; do Drug = 0 to 1; output; end; run;
proc phreg data = Gehan;
  model Week*Remiss(0) = Drug;
  baseline out = Outsurv covariates = Dummy survival = S upper = U lower = L / method = pl cltype = log nomean;
run;

*==================================================*;
* 図3.4.2 ;
*==================================================*;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc sgpanel data = Outsurv noautolegend;
  panelby Drug / novarname;
  step x = Week y = S / lineattrs = (pattern = 1);
  step x = Week y = L / lineattrs = (pattern = 2);
  step x = Week y = U / lineattrs = (pattern = 3);
  refline 0.5 / axis = y lineattrs = (pattern = 34 thickness = 1.5);
  rowaxis values = (0 to 1 by 0.1);
  colaxis values = (0 to 30 by 10);
  format Drug Drugf.;
run;
ods graphics off;
ods html close;

*==================================================*;
* プログラム3.5.1 ;
*==================================================*;

proc phreg data = Pcancer;
  model Time*Censor(0) = Age Sex Treat BUI CH P Stage PS 
              / selection = stepwise slentry = 0.20 slstay = 0.20 details;
run;

*==================================================*;
* プログラム3.5.2 ;
*==================================================*;

proc phreg data = Pcancer;
  model Time*Censor(0) = Age Sex Treat BUI CH P Stage PS 
              / selection = score best = 2;
run;

*==================================================*;
* 表3.5.2 ;
*==================================================*;

proc phreg data = Pcancer;
  model Time*Censor(0) = Treat;
proc phreg data = Pcancer;
  model Time*Censor(0) = Treat BUI;
proc phreg data = Pcancer;
  model Time*Censor(0) = Treat BUI Stage;
proc phreg data = Pcancer;
  model Time*Censor(0) = Age Treat BUI Stage;
proc phreg data = Pcancer;
  model Time*Censor(0) = Age Treat BUI P Stage;
proc phreg data = Pcancer;
  model Time*Censor(0) = Age Treat BUI P Stage PS;
proc phreg data = Pcancer;
  model Time*Censor(0) = Age Treat BUI CH P Stage PS;
proc phreg data = Pcancer;
  model Time*Censor(0) = Age Sex Treat BUI CH P Stage PS;
run;

*==================================================*;
* プログラム3.5.3 ;
*==================================================*;

data Dummy;
  do Stage = 3; do Treat = 0 to 1; do BUI = 0; output; end; end; end;
  do Stage = 4; do Treat = 0 to 1; do BUI = 1; output; end; end; end;
run;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc phreg data = Pcancer plots(overlay) = survival;
  model Time*Censor(0) = Treat BUI Stage / rl;
  baseline out = Out covariates = Dummy survival = S;
run;
ods graphics off;
ods html close;

*==================================================*;
* 図3.5.1 ;
*==================================================*;

data Out; set Out;
  if Stage = 3 then do;
    if Treat = 0 then Group = 1;
    else if Treat = 1 then Group = 2;
  end;
  else if Stage = 4 then do;
    if Treat = 0 then Group = 3;
    else if Treat = 1 then Group = 4;
  end;
run;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc sgplot data = Out noautolegend;
  step x = Time y = S / group = Group;
  yaxis values = (0 to 1 by 0.1);
  xaxis values = (0 to 25 by 5);
  keylegend / noborder;
  refline 10 / axis = x lineattrs = (pattern = 34 thickness = 1.5);
  refline 0.25 0.5 0.75 / axis = y lineattrs = (pattern = 34 thickness = 1.5);
run;
ods graphics off;
ods html close;

*==================================================*;
* プログラム3.6.1 ;
*==================================================*;

data Work;
input ID Group Time Censor @@;
cards;
1 0 4 1  2 0 5 1  3 0 7 0  4 0 9 1  5 1 6 1  6 1 10 1
7 1 11 1  8 1 12 0
;
proc phreg data = Work;
  model Time*Censor(0) = Group T_G;
  T_G = Time*Group;
  put ID =  Group =  Time =  Censor =  T_G = ;
run;

*==================================================*;
* プログラム3.6.2 ;
*==================================================*;

proc phreg data = Gehan;
  model Week*Remiss(0) = Drug X;
  if Drug = 0 then X = 0;
  if Drug = 1 then X = log(Week);
run;

*==================================================*;
* プログラム3.6.3 ;
*==================================================*;

proc phreg data = Heart;
   model Time*Status(0) = Z Acc_Age;
   if WaitTime = . or Time < WaitTime then Z = 0; else Z = 1;
run;

*==================================================*;
* プログラム3.6.4 ;
*==================================================*;

/* SELECT-WHERE構文 */
proc phreg data = Tumor;
  model Time*Censor(0) = Dose NPap;
  select (Time);
    when (27) NPap = P1; when (34) NPap = P2;
    when (37) NPap = P3; when (41) NPap = P4;
    when (43) NPap = P5; when (45) NPap = P6;
    when (46) NPap = P7; when (47) NPap = P8;
    when (49) NPap = P9; when (50) NPap = P10;
    when (51) NPap = P11; when (53) NPap = P12;
    when (65) NPap = P13; when (67) NPap = P14;
    when (71) NPap = P15; otherwise NPap = 0;
  end;
run;

/* ARRAYとIF-THEN構文 */
proc phreg data = Tumor;
  model Time*Censor(0) = Dose NPap;
  array pp{*} P1-P14; array tt{*} t1-t15;
  t1 = 27; t2 = 34; t3 = 37; t4 = 41; t5 = 43;
  t6 = 45; t7 = 46; t8 = 47; t9 = 49; t10 = 50;
  t11 = 51; t12 = 53; t13 = 65; t14 = 67; t15 = 71;
  if Time <  tt[1]  then NPap = 0;
  else if time >= tt[15] then NPap = P15;
  else do i = 1 to dim(pp);
    if tt[i] <= Time < tt[i+1] then NPap= pp[i];
  end;
run;

*==================================================*;
* 表3.6.1 ;
*==================================================*;

data Tumor2; set Tumor;
  array pp{*} P1-P14; array qq{*} P2-P15;
  array tt{1:15} _temporary_ (27 34 37 41 43 45 46 47 49 50 51 53 65 67 71);
  T1 = 0; T2 = 0; Censor2 = 0;
  if ( Time = tt[1] ) then do;
    T2 = tt[1]; NPap = p1; Censor2 = Censor; output;
  end;
  else do j = 1 to dim(pp);
    if ( tt[j] = Time ) then do;
      T2= Time; NPap = pp[j]; Censor2 = Censor; output;
    end;
    else if ( tt[j]  < Time ) then do;
      if ( pp[j]  ^= qq[j] ) then do;
        if qq[j] = . then T2 = Time; else T2= tt[j];
        NPap= pp[j]; Censor2= 0; output; T1 = T2;
      end;
    end;
  end;
  if ( Time >= tt[15] ) then do;
    T2 = Time; NPap = P15; Censor2 = Censor; output;
  end;
  keep ID Time Censor Dose T1 T2 NPap Censor2;
run;

*==================================================*;
* プログラム3.6.5 ;
*==================================================*;

proc phreg data = Tumor2;
  model (T1,T2)*Censor2(0) = Dose NPap;
  id ID;
run;

*==================================================*;
* 表3.6.2 ;
*==================================================*;

proc transpose data = Tumor out = TDC; by ID; var P1-P15;
data TDC; set TDC;
  where COL1 ^= . ; rename COL1 = NPap;
  t = input(substr(_NAME_, 2, 4), best.); 
data TDC; set TDC; by ID NPap; if first.NPap=1 then output;
data TDC; set TDC;
  if t = 1 then Ontime = 27;
  else if t = 2 then Ontime = 34;
  else if t = 3 then Ontime = 37;
  else if t = 4 then Ontime = 41;
  else if t = 5 then Ontime = 43;
  else if t = 6 then Ontime = 45;
  else if t = 7 then Ontime = 46;
  else if t = 8 then Ontime = 47;
  else if t = 9 then Ontime = 49;
  else if t = 10 then Ontime = 50;
  else if t = 11 then Ontime = 51;
  else if t = 12 then Ontime = 53;
  else if t = 13 then Ontime = 65;
  else if t = 14 then Ontime = 67;
  else if t = 15 then Ontime = 71;
  else Ontime = . ;
  keep ID Ontime Npap;
run;

*==================================================*;
* プログラム3.6.6 ;
*==================================================*;

/* 時間に依存しない情報を含んだデータセット */
data Base; set Tumor; keep ID Time Censor Dose;
run;

/* 個体ごとに情報を統合 */
proc transpose data = TDC out = Ontime(drop = _NAME_) prefix = T;
  var Ontime; by ID;
run;
proc transpose data = TDC out = NPap(drop = _NAME_) prefix = P;
  var NPap; by ID;
data Tumor3; merge Base Ontime NPap; by ID;
run;

*==================================================*;
* プログラム3.6.7 ;
*==================================================*;

proc phreg data = Tumor3;
  model Time*Censor(0) = Dose NPap;
  array P(5) P1-P5; array T(5) T1-T5;
  do i = 2 to 5; if T{ i } = .  or  T{ i } > Time then goto exit; end;
  exit: NPap = P{ i - 1 };
run;

*==================================================*;
* プログラム3.7.1 ;
*==================================================*;

data Data; X = 1; Censor = 1;
  do Strata = 0 to 1; do Treat = 0 to 1; do i = 1 to 10;
    input Time @@;output;
  end; end; end;
cards;
2.10	2.89	2.08	1.64	1.95	2.06	2.05	2.72	2.30	1.69
2.13	2.71	2.09	2.86	2.10	2.07	1.99	1.95	3.77	2.42
2.55	2.56	2.63	4.27	2.44	2.62	2.38	2.80	3.16	4.60
3.54	3.52	3.53	4.69	3.52	3.35	5.81	5.83	5.52	6.07
;
data Dummy; X = 1;
  do Strata = 0 to 1; do Treat = 0 to 1; output; end; end;
proc sort data = Data out = Data; by Strata Treat; run;
/*---  KM MODEL	 ------------------	(i)	-----------------*/
proc phreg data = Data;
  model Time*Censor(0) = X;
  by Strata Treat;
  baseline out = Out1 covariates = Dummy survival = S loglogs = LLS / nomean;
/*---  PH MODEL	BY STRATA  --- (ii) -----------------*/
proc phreg data = Data;
  model Time*Censor(0) = Treat;
  by Strata;
  baseline out = Out2 covariates = Dummy survival = S loglogs = LLS / nomean;
/*---  PH MODEL STRATA STRATA  -- (iii)	-----------------*/
proc phreg data = Data;
  model Time*Censor(0) = Treat;
  strata Strata;
  baseline out = Out3 covariates = Dummy survival = S loglogs = LLS / nomean;
/*---  PH MODEL  ------------------- (iv) -----------------*/
proc phreg data = Data;
  model Time*Censor(0) = Strata Treat;
  baseline out = Out4 covariates = Dummy survival = S loglogs = LLS / nomean;
run;

*==================================================*;
* 表3.7.1 ;
*==================================================*;

data Data; set Data;
  if Strata = 0 and Treat = 0 then Group = 1;
  else if Strata = 0 and Treat = 1 then Group = 2;
  else if Strata = 1 and Treat = 0 then Group = 3;
  else if Strata = 1 and Treat = 1 then Group = 4;
run;

proc means data = Data;
  var Time; by Group;
run;

*==================================================*;
* 図3.7.1-3.7.2 ;
*==================================================*;

proc sort data = Out3 out = Out3;
  by Treat Strata Time;
data Out3; set Out3;
  by Treat Strata Time; if first.Time = 1 then output;
run;

data Out1; set Out1; Model=1;
data Out2; set Out2; Model=2;
data Out3; set Out3; Model=3;
data Out4; set Out4; Model=4;
data Out; set Out1-Out4;
run;

data Out; set Out;
  if Strata = 0 and Treat = 0 then Group = 1;
  else if Strata = 0 and Treat = 1 then Group = 2;
  else if Strata = 1 and Treat = 0 then Group = 3;
  else if Strata = 1 and Treat = 1 then Group = 4;
run;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc sgpanel data = Out noautolegend;
  panelby Model / novarname;
  step x = Time y = S / group = Group name = "Survival";
  rowaxis values = (0 to 1 by 0.1) label = "S"; 
  colaxis values = (0 to 7 by 1) label = "Time";
  keylegend "Survival" / noborder;
  format Model Modelf. Group Groupf.;
run;
proc sgpanel data = Out noautolegend;
  panelby Model / novarname;
  step x = Time y = LLS / group = Group name = "Survival";
  rowaxis values = (-6 to 4 by 1) label = "LLS"; 
  colaxis values = (1 to 8 by 1) label = "Time";
  keylegend "Survival" / noborder;
  format Model Modelf. Group Groupf.;
run;
ods graphics off;
ods html close;

*==================================================*;
* 図3.7.3 ;
*==================================================*;

/*---  PH MODEL  ------------------- (v) -----------------*/
proc phreg data = Data;
  model Time*Censor(0) = Strata Treat Strata*Treat;
  baseline out = Out covariates = Dummy survival = S loglogs = LLS / nomean;
run;

data Out; set Out;
  if Strata = 0 and Treat = 0 then Group = 1;
  else if Strata = 0 and Treat = 1 then Group = 2;
  else if Strata = 1 and Treat = 0 then Group = 3;
  else if Strata = 1 and Treat = 1 then Group = 4;
run;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc sgplot data = Out noautolegend;
  step x = Time y = S / group = Group name = "Survival";
  yaxis values = (0 to 1 by 0.1) label = "S"; 
  xaxis values = (0 to 7 by 1) label = "Time";
  keylegend "Survival" / noborder;
  format Group Groupf.;
run;
proc sgplot data = Out noautolegend;
  step x = Time y = LLS / group = Group name = "Survival";
  yaxis values = (-6 to 4 by 1) label = "LLS"; 
  xaxis values = (1 to 7 by 1) label = "Time";
  keylegend "Survival" / noborder;
  format Group Groupf.;
run;
ods graphics off;
ods html close;

*==================================================*;
* プログラム3.7.2 ;
*==================================================*;

/*  交互作用を含むモデル            */ 
proc phreg data = Pcancer;
  class Stage Treat / param = effect ref = first;
  model Time*Censor(0) = Stage Treat Stage*Treat / rl;
run;
/*  交互作用を含まないモデル         */ 
proc phreg data = Pcancer;
  class Stage Treat / param = effect ref = first;
  model Time*Censor(0) = Stage Treat / rl;
run;
proc phreg data = Pcancer;
  class Stage Treat / param = ref ref = first;
  model Time*Censor(0) = Stage Treat / rl;
run;

*==================================================*;
* プログラム3.7.3 ;
*==================================================*;

proc phreg data = Pcancer;
  model Time*Censor(0) = T1 T2 T3 T4 S1 S2 S3 S4 / rl;
  T1 = 0; T2 = 0; T3 = 0; T4 = 0; S1 = 0; S2 = 0; S3 = 0; S4 = 0;
  select;
    when( 0 < Time <= 2 ) do; T1 = Treat; S1 = Stage; end; 
    when( 2 < Time <= 4 ) do; T2 = Treat; S2 = Stage; end; 
    when( 4 < Time <= 8 ) do; T3 = Treat; S3 = Stage; end; 
    when( 8< Time ) do; T4 = Treat; S4 = Stage; end; 
  otherwise; end;
run;

*==================================================*;
* 図3.7.4 ;
*==================================================*;

/*---	KM	*/
ods listing close;
proc lifetest data = Pcancer outs = Out1;
  time Time*Censor(0); strata Treat;
  format Treat Treatf.;
run;
ods listing;

data Out1; set Out1;
  LLS = log(-log(Survival)); LogTime = log(Time); Method=1;
run;

data Dummy;	do Treat = 0 to 1; output; end;
run;

/*---	PH	*/
proc phreg data = Pcancer;
  model Time*Censor(0) = Treat;
  baseline out = Out2 covariates = Dummy survival = Survival loglogs = LLS / nomean;
run;

data Out2; set Out2;
  LogTime = log(Time); Method=2;
run;

data Out; set Out1-Out2;
run;

/*---	グラフ作成    */
ods html image_dpi=645 style=journal;
ods graphics / reset noborder noscale
               width=500px height=300px
               imagefmt=png
;
proc sgpanel data = Out noautolegend;
  panelby Method / novarname;
  step x = Time y = Survival / group = Treat  name = "Survival";
  rowaxis values = (0 to 1 by 0.1) label = "S"; 
  colaxis values = (0 to 30 by 10) label = "Time";
  keylegend "Survival" / noborder;
  format Method Method_survf. Treat Treatf.;
  refline 2 4 8 / axis = x lineattrs = (pattern = 34 thickness = 1.5);
run;
proc sgpanel data = Out noautolegend;
  panelby Method / novarname;
  step x = LogTime y = LLS / group = Treat  name = "Survival";
  rowaxis values = (-5 to 3 by 1) label = "LLS"; 
  colaxis values = (-3 to 3 by 1) label = "log(Time)";
  keylegend "Survival" / noborder;
  format Method Method_llsf. Treat Treatf.;
  refline 0.7 1.4 2.1 / axis = x lineattrs = (pattern = 34 thickness = 1.5);
run;
ods graphics off;
ods html close;

*==================================================*;
* 図3.7.5 ;
*==================================================*;

/*---	KM	*/
ods listing close;
proc lifetest data = Pcancer outs = Out1;
  time Time*Censor(0); strata Stage;
  format Stage Stagef.;
run;
ods listing;

data Out1; set Out1;
  LLS = log(-log(Survival)); LogTime = log(Time); Method=1;
run;

data Dummy;	do Stage = 3 to 4; output; end;
run;

/*---	PH	*/
proc phreg data = Pcancer;
  model Time*Censor(0) = Stage;
  baseline out = Out2 covariates = Dummy survival = Survival loglogs = LLS / nomean;
run;

data Out2; set Out2;
  LogTime = log(Time); Method=2;
run;

data Out; set Out1-Out2;
run;

/*---	グラフ作成    */
ods html image_dpi=645 style=journal;
ods graphics / reset noborder noscale
               width=500px height=300px
               imagefmt=png
;
proc sgpanel data = Out noautolegend;
  panelby Method / novarname;
  step x = Time y = Survival / group = Stage  name = "Survival";
  rowaxis values = (0 to 1 by 0.1) label = "S"; 
  colaxis values = (0 to 30 by 10) label = "Time";
  keylegend "Survival" / noborder;
  format Method Method_survf. Stage Stagef.;
  refline 2 4 8 / axis = x lineattrs = (pattern = 34 thickness = 1.5);
run;
proc sgpanel data = Out noautolegend;
  panelby Method / novarname;
  step x = LogTime y = LLS / group = Stage  name = "Survival";
  rowaxis values = (-5 to 3 by 1) label = "LLS"; 
  colaxis values = (-3 to 3 by 1) label = "log(Time)";
  keylegend "Survival" / noborder;
  format Method Method_llsf. Stage Stagef.;
  refline 0.7 1.4 2.1 / axis = x lineattrs = (pattern = 34 thickness = 1.5);
run;
ods graphics off;
ods html close;

*==================================================*;
* プログラム3.7.4 ;
*==================================================*;

data Scancer2; set Scancer; LogDose = log(Dose);
  select (Dose);
    when (10)	do; X1 = 0; X2 = 0; end;
    when (30)	do; X1 = 1; X2 = 0; end;
    when (90)	do; X1 = 0; X2 = 1; end;
  end;
run;
/*-----------	(i)	------------------------------------------*/
proc phreg data = Scancer2;
  model Time*Censor(0) = X1 X2;
run;
/*-----------	(ii) ------------------------------------------*/
proc phreg data = Scancer2;
  model Time*Censor(0) = Dose;
run;
/*-----------	(iii) ------------------------------------------*/
proc phreg data = Scancer2;
  model Time*Censor(0) = LogDose;
run;

*==================================================*;
* 図3.7.6-3.7.7 ;
*==================================================*;

/*---	KM	*/
ods listing close;
proc lifetest data = Scancer2 outs = Out0;
  time Time*Censor(0); strata Dose;
  format Dose Dosecf.;
run;
ods listing;

data Out0; set Out0; LLS = log(-log(Survival)); LogTime = log(Time);
run;

/*---	PH	*/
/*-----------	(i) ------------------------------------------*/
data Dummy; do Dose = 10, 30, 90; LogDose = log(Dose); output; end;
proc phreg data = Scancer2;
  class Dose;
  model Time*Censor(0) = Dose;
  baseline out = Out1 covariates = Dummy survival = Survival loglogs = LLS / nomean;
run;
/*----------- (ii) ------------------------------------------*/
proc phreg data = Scancer2;
  model Time*Censor(0) = Dose;
  baseline out = Out2 covariates = Dummy survival = Survival loglogs = LLS / nomean;
run;
/*----------- (iii) ------------------------------------------*/
proc phreg data = Scancer2;
  model Time*Censor(0) = LogDose;
  baseline out = Out3 covariates = Dummy survival = Survival loglogs = LLS / nomean;
run;

/*---	グラフ作成    */
data Out0; set Out0; Method = 0;
data Out1; set Out1; Method = 1;
data Out2; set Out2; Method = 2;
data Out3; set Out3; Method = 3;
data Out; set Out0-Out3;
run;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc sgpanel data = Out noautolegend;
  panelby Method / novarname;
  step x = Time y = Survival / group = Dose name="Survival";
  rowaxis values = (0 to 1 by 0.1) label = "S"; 
  colaxis values = (0 to 70 by 10) label = "Time";
  keylegend "Survival" / noborder;
  format Method Model_Scancerf. Dose Dosecf.;
run;
proc sgpanel data = Out noautolegend;
  panelby Method / novarname;
  step x = Time y = LLS / group = Dose name="Survival";
  rowaxis values = (-6 to 2 by 1) label="LLS"; 
  colaxis values = (0 to 70 by 10) label = "Time";
  keylegend "Survival" / noborder;
  format Method Model_Scancerf. Dose Dosecf.;
run;
ods graphics off;
ods html close;

*==================================================*;
* プログラム3.7.5 ;
*==================================================*;

data Pcancer; set Pcancer;
  select;
    when (         Age < 50 )  do; Cat = 1; A1 = 0; A2 = 0; A3 = 0; end;
    when ( 50 <= Age <60 )  do; Cat = 2; A1 = 1; A2 = 0; A3 = 0; end;
    when ( 60 <= Age <70 )  do; Cat = 3; A1 = 0; A2 = 1; A3 = 0; end;
    when ( 70 <= Age       )  do; Cat = 4; A1 = 0; A2 = 0; A3 = 1; end;
  otherwise; end;
run;

/* ダミー変数A1, A2, A3を用いた解析 */
proc phreg data = Pcancer;
  model Time*Censor(0) = Sex Treat A1 A2 A3; 
  Age: test A1 = A2 = A3 = 0;
run;
/* ダミー変数A1, A2, A3を含めない解析 */
proc phreg data = Pcancer;
  model Time*Censor(0) = Sex Treat;
run;
/* CLASS文による変数Catのダミー変数化を用いた解析 */
proc phreg data = Pcancer;
  class Cat / param = ref ref = first;
  model Time*Censor(0) = Sex Treat Cat;
run;

*==================================================*;
* プログラム3.7.6 ;
*==================================================*;

proc phreg data = Pcancer;
  model Time*Censor(0) = Treat BUI Stage;
  output out = Residual dfbeta = DFB1-DFB3 ld = LD lmax = Lmax;
run;

*==================================================*;
* 図3.7.8 ;
*==================================================*;

data Residual; set Residual;
  if ID in(72) then ID2 = ID; else ID2 = . ;
  label DFB1="Treatに対するDFBETA統計量"
  DFB2="BUIに対するDFBETA統計量"
  DFB3="Stageに対するDFBETA統計量"
  LD="LD統計量"
  Lmax="LMAX統計量"
  ;
run;

data Residual; set Residual;
  if ID in(7, 17, 80, 84) then ID2 = ID; else ID2 = . ;
run;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc sgscatter data = Residual;
  matrix DFB1-DFB3 LD Lmax / datalabel = ID2;
run;
ods graphics off;
ods html close;

*==================================================*;
* 表3.7.8 ;
*==================================================*;

proc phreg data = Pcancer;
  model Time*Censor(0) = Treat BUI Stage;
  output out = Residual dfbeta = DFB1-DFB3 ld = LD lmax = Lmax;
run;

proc phreg data = Pcancer;
  where ID notin(7);
  model Time*Censor(0) = Treat BUI Stage;
  output out = Residual dfbeta = DFB1-DFB3 ld = LD lmax = Lmax;
run;

*==================================================*;
* プログラム3.8.1 ;
*==================================================*;

ods html image_dpi=1000 style=journal;
ods graphics / reset noborder noscale
               width=250px height=250px
               imagefmt=png
;
proc phreg data = Liver plots(overlay = individual) = (roc) rocoptions(method = km at = (5, 7));
  model Time*Status(0) = Bilirubin Age Edema;
run;
ods graphics off;
ods html close;

*==================================================*;
* 図3.8.2 ;
*==================================================*;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc phreg data = Liver plots = auc rocoptions(method = km);
  model Time*Status(0) = Bilirubin Age Edema;
run;
ods graphics off;
ods html close;

*==================================================*;
* プログラム3.8.2 ;
*==================================================*;

proc phreg data = Liver concordance = harrell;
   model Time*Status(0) = Bilirubin Age Edema;
run;

*==================================================*;
* END ;
*==================================================*;
