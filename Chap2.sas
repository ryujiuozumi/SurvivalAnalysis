*======================================================================*;
/*

大橋靖雄・浜田知久馬・魚住龍史 (2022)
生存時間解析 第2版 --SASによる生物統計
東京大学出版会

第2章　生存関数のノンパラメトリックな推定と検定 (LIFETEST プロシジャ)

Program : Chap2.sas
Made by : 魚住龍史 (e-mail : uozumi at kuhp.kyoto-u.ac.jp)
Date : 2022年8月

*/
*======================================================================*;

*==================================================*;
* 図2.1.4 ;
*==================================================*;

data Work;
input Group Time Censor @@;
cards;
1 4 1  1 5 1  1 9 1
2 6 1  2 10 1  2 11 1
;
run;

proc lifetest data = Work notable;
  time Time*Censor(0);
  strata Group;
run;

*==================================================*;
* 図2.1.5 ;
*==================================================*;

data Work;
input Group Time Censor @@;
cards;
1 4 1  1 5 1  1 9 1  1 7 0
2 6 1  2 10 1  2 11 1  2 12 0
;
run;

proc lifetest data = Work notable;
  time Time*Censor(0);
  strata Group;
run;

*==================================================*;
* プログラム2.3.1 ;
*==================================================*;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc lifetest data = Gehan plots = (s(atrisk outside)) atrisk outs = OutS;
  time Week*Remiss(0);
  strata Drug;
  format Drug Drugf.;
run;
ods graphics off;
ods html close;

proc sort data = OutS out = OutS; by Drug; run;
proc print data = OutS; by Drug; run;

*==================================================*;
* プログラム2.3.2 ;
*==================================================*;

data OutS; set OutS;
  retain SURVIVAL_t; by Drug;
  if SURVIVAL ^= . then SURVIVAL_t = SURVIVAL;
  else SURVIVAL = SURVIVAL_t;
run;

goptions reset = all;
goptions vsize = 10 in hsize = 15 in htitle = 3.5cells htext = 3.5cells;
options linesize = 98 pagesize = 200;
/*filename gplotout "C:\Users\Ryuji\Documents\SAS_GRAPH\Gplot.bmp";*/
goptions device = bmp gsfname = gplotout gsfmode = replace;
goptions ftext = "Times New Roman";
proc gplot data = OutS;
  plot SURVIVAL*Week = Drug /  noframe legend = legend1; 
  legend1 value = (tick = 2 justify = l); label SURVIVAL = 'S';
  symbol1 i = steplj l = 1 v = none c = green w = 3;
  symbol2 i = steplj l = 2 v = none c = red w = 3;
run; quit;

*==================================================*;
* プログラム2.3.3 ;
*==================================================*;

data OutS2; set OutS;
  Surv = SURVIVAL; output;
  if _CENSOR_ = 1 then do;
    Surv = SURVIVAL + 0.05; output;  /* 0.05: Length of whisker*/
    Surv = SURVIVAL; output;
  end;
run;

data OutS3; set OutS;
  LLS = log(-log(SURVIVAL)); LogWeek = log(Week);
run;

/* カプラン・マイヤープロット作成 (ひげ付き) */
proc gplot data = OutS2;
  plot Surv*Week = Drug / noframe legend = legend1; 
  legend1 value = (tick = 2 justify = l); label Surv = 'S';
  symbol1 i = steplj l = 1 v = none c = green w = 3;
  symbol2 i = steplj l = 2 v = none c = red w = 3;
run; quit;

/* 2重対数プロット作成 */
proc gplot data = Outs3;
  plot LLS*LogWeek = Drug / noframe legend = legend1;
  legend1 value = (tick = 2 justify = l);
  symbol1 i = steplj l = 1 v = none c = green w = 3;
  symbol2 i = steplj l = 2 v = none c = red w = 3;
run; quit;

*==================================================*;
* プログラム2.4.1 ;
*==================================================*;

proc lifetest data = Pcancer outt = OutT notable;
  time Time*Censor(0);
  test Age Sex Treat BUI CH P Stage PS;
run;
proc print data = OutT; run;

*==================================================*;
* プログラム2.5.1 ;
*==================================================*;

proc lifetest data = Scancer notable;
  time Time*Censor(0); strata Dose;
run;

*==================================================*;
* プログラム2.5.2 ;
*==================================================*;

proc lifetest data = Scancer notable;
  where Dose in(10, 30);
  time Time*Censor(0); strata Dose;
run;

*==================================================*;
* 表2.5.1 ;
*==================================================*;

proc lifetest data = Scancer notable;
  where Dose in(10, 30);
  time Time*Censor(0); strata Dose;
run;
proc lifetest data = Scancer notable;
  where Dose in(10, 90);
  time Time*Censor(0); strata Dose;
run;
proc lifetest data = Scancer notable;
  where Dose in(30, 90);
  time Time*Censor(0); strata Dose;
run;

*==================================================*;
* プログラム2.5.3 ;
*==================================================*;

proc iml; reset print nolog;
  u = {-13.8626, 4.814431, 9.048202};
  V = {10.2810 -5.7934 -4.4875,
        -5.7934 9.0072 -3.2138,
        -4.4875 -3.2138 7.7013};
  CT = {10, 30, 90};
  C1 = {-1, 0, 1};
  C2 = {-1, 2, -1};
  X2L = u`*ginv(V)*u;
  X2T = (CT`*u)##2/(CT`*V*CT);
  X21 = (C1`*u)##2/(C1`*V*C1);
  X22 = (C2`*u)##2/(C2`*V*C2); 
  PL = (1-probchi(X2L, 2));
  PT = (1-probchi(X2T, 1));
  P1 = (1-probchi(X21, 1));
  P2 = (1-probchi(X22, 1));
run; quit;

*==================================================*;
* 図2.5.2 ;
*==================================================*;

proc iml; reset print nolog;
  u = {-13.8626, 4.814431, 9.048202};
  V = {10.2810 -5.7934 -4.4875,
        -5.7934 9.0072 -3.2138,
        -4.4875 -3.2138 7.7013};
  IV = inv(V);
  GIV = ginv(V);
  DET = det(V);
  V={10.2809 -5.7934 -4.4875,
      -5.7934 9.0072 -3.2138,
      -4.4875 -3.2138 7.7013};
  IV = inv(V);
  GIV = ginv(V);
  DET = det(V);
  X2L = u`*ginv(V)*u;
run; quit;

*==================================================*;
* 図2.5.3 ;
*==================================================*;

proc iml; reset print nolog;
  u = {4.814431, 9.048202};
  V = { 9.0072 -3.2138,
        -3.2138 7.7013};
  IV = inv(V);
  GIV = ginv(V);
  X2L = u`*inv(V)*u;
run; quit;

*==================================================*;
* プログラム2.5.4 ;
*==================================================*;

data Scancer2; set Scancer; if Censor = 1 then Censor = 2; run;
proc multtest data = Scancer2 out = Out noprint;
  class Dose; test peto(Censor / time = Time);
  contrast 'Tarone'   10 30 90;
  contrast 'Linear'   -1 0 1;
  contrast 'Quadratic'   -1 2 -1;
run;
data Out; set Out; where _tstrat_ = .; X2 = cinv(1 - raw_p, 1); run;
proc print data = Out label;
  var _test_ _var_ _contrast_ _exp_ _se_ raw_p X2;
run;

*==================================================*;
* 図2.6.1 ;
*==================================================*;

proc iml; reset print nolog;
  u = {4.814431, 9.048202};
  V = {9.007 -3.2138,
        -3.2138 7.7013};
  IV = inv(V);
  R = IV*u;
  EXPR = exp(R);
run; quit;

*==================================================*;
* プログラム2.6.1 ;
*==================================================*;

data Scancer2; set Scancer;
  select(Dose);
    when(10) do; X1 = 0; X2 = 0; end;
    when(30) do; X1 = 1; X2 = 0; end;
    when(90) do; X1 = 0; X2 = 1; end;
  end;
run;

proc phreg data = Scancer2;
  model Time*Censor(0) = X1 X2 / itprint;
  Dose: test X1 = X2 = 0;
run;

/* CLASS文でダミー変数化 */
proc phreg data = Scancer;
  class Dose / param = ref ref = first;
  model Time*Censor(0) = Dose / itprint;
run;

*==================================================*;
* プログラム2.7.1 ;
*==================================================*;

proc sort data = Pcancer; by Stage; run;
/* 層別解析 */
proc lifetest data = Pcancer notable;
  time Time*Censor(0);
  strata Treat; by Stage;
run;
/* 平均的な効果の解析 */
proc lifetest data = Pcancer notable;
  time Time*Censor(0);
  strata Stage / group = Treat;
run;

*==================================================*;
* プログラム2.7.2 ;
*==================================================*;

/* 層別解析 */
proc phreg data = Pcancer;
  model Time*Censor(0) = Treat; 
run;
/* 平均的な効果の解析 */
proc phreg data = Pcancer;
  model Time*Censor(0) = Treat; strata Stage;
run;

*==================================================*;
* 図2.7.1-2.7.2 ;
*==================================================*;

ods listing close;
proc lifetest data = Pcancer outs = OutS;
  time Time*Censor(0);
  strata Treat;
  by Stage;
  format Treat Treatf.;
run;
ods listing;

data OutS; set OutS;
  retain SURVIVAL_t; by Stage;
  if SURVIVAL ^= . then SURVIVAL_t = SURVIVAL;
  else SURVIVAL = SURVIVAL_t;
run;

data OutS; set OutS;
 LLS = log(-log(SURVIVAL));
 LogTime = log(Time);
run;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc sgpanel data = OutS;
  panelby Stage / novarname;
  step x = Time y = SURVIVAL / group = Treat;
  keylegend / noborder;
  rowaxis values = (0 to 1 by 0.1) label = "　"; 
  colaxis values = (0 to 30 by 10); 
  format stage stagef.;
run;
proc sgpanel data = OutS;
  panelby Stage / novarname;
  step x = LogTime y = LLS / group = Treat;
  keylegend / noborder;
  rowaxis values = (-4 to 2 by 1) label = "　"; 
  colaxis values = (-1 to 4 by 1) label = "log(Time)"; 
  format stage stagef.;
run;
ods graphics off;
ods html close;

*==================================================*;
* プログラム2.8.1 ;
*==================================================*;

data Work; 
input Time Status @@; 
cards;
5 0  6 1  7 0  9 2  11 1
;
run;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc lifetest data = Work;
  time Time*Status(0) / eventcode = 1;
run;
ods graphics off;
ods html close;

*==================================================*;
* プログラムB.1 ;
*==================================================*;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc lifetest data = Gehan method = life intervals = 0 to 24 by 3;
  time Week*Remiss(0);
  strata Drug;
  format Drug Drugf.;
run;
ods graphics off;
ods html close;

*==================================================*;
* プログラムD.1 ;
*==================================================*;

proc lifetest data = Gehan rmst(tau=20 bc cl);
  time Week*Remiss(0); strata Drug;
run;

*==================================================*;
* END ;
*==================================================*;
