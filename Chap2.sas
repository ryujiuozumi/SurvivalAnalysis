*======================================================================*;
/*

�勴���Y�E�l�c�m�v�n�E���Z���j (2022)
�������ԉ�� ��2�� --SAS�ɂ�鐶�����v
������w�o�ŉ�

��2�́@�����֐��̃m���p�����g���b�N�Ȑ���ƌ��� (LIFETEST �v���V�W��)

Program : Chap2.sas
Made by : ���Z���j (e-mail : uozumi at kuhp.kyoto-u.ac.jp)
Date : 2022�N7��

*/
*======================================================================*;

*==================================================*;
* �}2.1.4 ;
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
* �}2.1.5 ;
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
* �t�H�[�}�b�g ;
*==================================================*;

proc format;
  value Drugf 0='CONTROL' 1='6-MP';
  value Treatf 0='�p���Ǝ˂Ȃ�' 1='�p���Ǝ˂���'; 
  value Stagef 3='�X�e�[�W3' 4='�X�e�[�W4';
run;

*==================================================*;
* Gehan�̔����a�̃f�[�^ ;
*==================================================*;

data Gehan;
do Drug=0,1;
   do i=1 to 21; input Week Remiss @@;output;end;end;
cards;
 1 1   1 1   2 1   2 1   3 1   4 1   4 1   5 1   5 1   8 1   8 1   8 1
 8 1  11 1  11 1  12 1  12 1  15 1  17 1  22 1  23 1
 6 0   6 1   6 1   6 1   7 1   9 0  10 0  10 1  11 0  13 1  16 1   17 0
 19 0  20 0  22 1  23 1  25 0  32 0  32 0  34 0  35 0
;
run;

*==================================================*;
* �v���O����2.3.1 ;
*==================================================*;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
;
proc lifetest data = Gehan plots = (s(atrisk outside)) atrisk outs=OutS;
  time Week*Remiss(0);
  strata Drug;
  format Drug Drugf.;
run;
ods graphics off;
ods html close;

proc sort data = OutS out = OutS; by Drug; run;
proc print data = OutS; by Drug; run;

*==================================================*;
* �v���O����2.3.2 ;
*==================================================*;

data OutS; set OutS;
  retain SURVIVAL_t; by Drug;
  if SURVIVAL ^= . then SURVIVAL_t = SURVIVAL;
  else SURVIVAL = SURVIVAL_t;
run;

goptions reset = all;
goptions vsize = 10 in hsize = 15 in htitle = 3.5cells htext = 3.5cells;
options linesize = 98 pagesize = 200;
filename gplotout "C:\Users\Ryuji\Documents\SAS_GRAPH\Gplot.bmp";
goptions device = bmp gsfname = gplotout gsfmode = replace;
goptions ftext = "Times New Roman";
proc gplot data = OutS;
 plot SURVIVAL*Week = Drug /  noframe; label SURVIVAL = 'S';
 symbol1 i = steplj l = 1 v = none c = green w = 3;
 symbol2 i = steplj l = 2 v = none c = red w = 3;
run; quit;

*==================================================*;
* �v���O����2.3.3 ;
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

/* �J�v�����E�}�C���[�v���b�g�쐬 (�Ђ��t��) */
proc gplot data = OutS2;
 plot Surv*Week = Drug / noframe; label Surv = 'S';
 symbol1 i = steplj l = 1 v = none c = green w = 3;
 symbol2 i = steplj l = 2 v = none c = red w = 3;
run; quit;

/* 2�d�ΐ��v���b�g�쐬 */
proc gplot data = outs3;
 plot LLS*LogWeek = Drug / noframe;
 symbol1 i = steplj l = 1 v = none c = green w = 3;
 symbol2 i = steplj l = 2 v = none c = red w = 3;
run; quit;

*==================================================*;
* �X�����̃f�[�^ ;
*==================================================*;

data Pcancer;
input ID Time Censor Age Sex Treat BUI CH P Stage PS; x=1;
cards;
 1    2.4   1     66   0   0   1   4   1   4   3
 2    1.7   1     69   0   0   1   4   1   4   3
 3    0.1   1     48   0   0   1   1   0   3   2
 4    1.0   1     73   0   0   1   4   0   3   4
 5    4.8   1     65   0   0   1   4   0   4   2
 6    6.4   1     38   0   0   0   4   0   3   3
 7   10.8   1     62   1   0   1   4   0   3   3
 9    5.1   1     59   1   0   1   1   1   4   3
10    1.1   1     53   0   0   1   1   1   4   3
11    0.5   1     70   0   0   1   1   1   4   3
12    0.8   1     71   0   0   0   3   0   4   3
14    4.0   1     61   1   0   1   4   0   4   3
15    4.0   1     69   0   0   0   4   0   3   3
16    4.0   1     41   1   0   1   4   0   4   4
17    8.5   1     49   0   0   1   4   0   3   2
18    3.6   1     56   0   0   0   4   0   3   2
19    6.9   1     59   1   0   0   4   1   3   4
20    6.2   1     53   0   0   0   3   1   4   4
21    1.0   1     72   1   0   0   3   0   4   2
22    6.2   1     57   1   0   0   3   0   3   2
23    4.3   1     49   0   0   0   3   0   4   2
24    3.1   1     74   0   0   1   4   0   3   3
25    8.3   1     43   1   1   1   4   0   4   2
26   12.7   1     60   1   1   1   4   1   3   3
28    4.9   1     55   1   1   1   4   0   3   3
30    2.7   1     70   0   1   1   3   0   3   3
33   10.6   1     63   0   1   1   1   0   3   2
35   18.2   1     69   1   1   0   4   0   3   2
36    1.4   1     66   1   1   1   1   1   4   2
37    5.8   1     58   1   1   1   1   0   4   2
38    3.0   1     67   1   1   1   1   1   4   3
39    1.5   1     74   0   1   1   1   1   4   2
40    2.4   1     77   1   1   1   1   1   3   2
41    2.0   1     70   1   1   1   4   0   3   4
42    1.1   1     75   0   1   1   4   0   4   2
43    2.5   1     65   1   1   1   2   0   3   2
44    5.4   1     71   1   1   1   4   0   3   2
45    4.4   1     50   0   1   1   1   1   4   2
46    4.8   1     56   0   1   1   3   0   3   2
47    3.1   1     68   1   1   1   1   0   3   1
48    5.6   1     65   1   1   1   3   0   3   3
49    3.1   1     65   0   1   0   3   0   3   2
50    1.3   1     43   1   1   1   3   0   3   2
51   11.5   1     83   0   1   0   3   0   3   2
53    3.8   1     65   0   1   1   2   0   3   2
54    2.9   1     63   0   1   1   1   0   4   3
55    2.2   1     47   1   1   1   2   0   4   2
56    1.7   1     75   1   1   1   2   1   4   2
58    3.5   1     63   1   1   1   1   1   4   2
62   11.3   1     54   0   1   1   2   0   3   2
63    9.0   1     56   1   1   0   4   0   4   3
64   12.5   1     50   0   1   0   1   0   3   1
65    6.8   1     62   0   1   1   1   0   3   1
66   10.8   1     53   1   1   1   1   0   3   3
67    3.0   1     63   0   1   1   1   1   4   4
68    1.8   1     59   0   1   1   4   0   4   4
69    5.0   1     66   1   1   1   4   0   4   2
70    8.0   1     62   0   1   0   3   0   3   2
71    6.8   1     72   0   1   1   1   0   3   2
72   11.1   0     54   0   1   0   3   0   3   1
73    9.4   1     68   1   1   0   3   0   3   2
74    3.9   1     63   1   1   0   4   0   3   2
75    2.1   1     68   1   1   1   2   0   3   2
76    4.3   1     48   1   1   1   4   0   3   3
77    9.3   1     68   0   1   1   1   0   4   2
78    8.8   1     75   1   1   1   1   0   4   2
79    2.4   1     49   0   1   1   1   1   4   3
80   21.6   1     62   0   1   1   1   0   4   3
81    5.6   1     56   1   1   1   1   0   4   2
83   11.4   1     56   0   1   1   1   0   3   2
84   18.3   1     59   1   1   1   1   1   4   3
85    9.2   1     59   0   1   1   1   1   3   2
86    4.5   1     48   0   1   1   2   1   3   3
87    8.2   1     64   0   1   1   2   0   3   2
88   15.0   1     60   0   1   1   4   0   3   4
89    6.9   1     75   0   1   1   1   0   4   2
90    3.5   1     46   0   1   1   1   1   4   2
91    2.1   1     53   1   1   1   1   1   4   2
92    3.1   1     62   0   1   1   3   0   3   3
93    3.2   1     47   0   1   1   3   0   4   3
94    1.9   1     62   0   1   1   1   0   4   3
95    2.1   1     55   0   1   1   1   1   3   2
96    7.0   1     80   0   1   1   1   0   4   3
;
run;

*==================================================*;
* �v���O����2.4.1 ;
*==================================================*;

proc lifetest data = Pcancer outt = OutT notable;
  time Time*Censor(0);
  test Age Sex Treat BUI CH P Stage PS;
run;
proc print data = OutT; run;

*==================================================*;
* �畆���̃f�[�^ ;
*==================================================*;

data Scancer;
do Dose=10,30,90; do i=1 to 30;
   input Time Censor @@; output; end; end;
cards;
40 1  76 0  76 0  76 0  64 0  66 1  76 0  76 0  76 0  76 0
32 1  40 1  60 1  72 0  76 0  44 1  62 1  60 0  76 0  76 0
40 1  42 1  60 1  76 0  76 0  48 1  76 0  76 0  76 0  76 0
26 0  46 1  32 1  49 0  44 1  44 0  43 0  40 1  44 1  45 1
22 1  43 0  48 1  44 1  44 1  36 1  44 1  42 1  45 1  49 0
33 0  38 1  48 0  48 0  47 0  41 1  46 1  46 1  38 1  35 0
36 1  40 1  44 1  44 1  49 0  29 1  28 1  34 0  48 1  49 0
40 1  42 1  40 1  38 1  38 1  32 1  38 1  32 1  49 0  22 1
32 1  38 1  48 0  23 0  32 1  49 0  44 0  45 1  49 0   1 0
;
run;

*==================================================*;
* �v���O����2.5.1 ;
*==================================================*;

proc lifetest data = Scancer notable;
  time Time*Censor(0); strata Dose;
run;

*==================================================*;
* �v���O����2.5.2 ;
*==================================================*;

proc lifetest data = Scancer notable;
  where Dose in(10, 30);
  time Time*Censor(0); strata Dose;
run;

*==================================================*;
* �\2.5.1 ;
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
* �v���O����2.5.3 ;
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
* �}2.5.2 ;
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
* �}2.5.3 ;
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
* �v���O����2.5.4 ;
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
* �}2.6.1 ;
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
* �v���O����2.6.1 ;
*==================================================*;

data Scancer2; set Scancer;
  select(Dose);
    when(10) do; X1 = 0; X2 = 0; end;
    when(30) do; X1 = 1; X2 = 0; end;
    when(90) do; X1 = 0; X2 = 1; end;
  end;
run;

proc phreg data=Scancer2;
  model Time*Censor(0) = X1 X2 / itprint;
  Dose: test X1 = X2 = 0;
run;

/* CLASS���Ń_�~�[�ϐ��� */
proc phreg data = Scancer;
  class Dose / param = ref ref = first;
  model Time*Censor(0) = Dose / itprint;
run;

*==================================================*;
* �v���O����2.7.1 ;
*==================================================*;

proc sort data = Pcancer; by Stage; run;
/* �w�ʉ�� */
proc lifetest data = Pcancer notable;
  time Time*Censor(0);
  strata Treat; by Stage;
run;
/* ���ϓI�Ȍ��ʂ̉�� */
proc lifetest data = Pcancer notable;
  time Time*Censor(0);
  strata Stage / group = Treat;
run;

*==================================================*;
* �v���O����2.7.2 ;
*==================================================*;

/* �w�ʉ�� */
proc phreg data = Pcancer;
  model Time*Censor(0) = Treat; 
run;
/* ���ϓI�Ȍ��ʂ̉�� */
proc phreg data = Pcancer;
  model Time*Censor(0) = Treat; strata Stage;
run;

*==================================================*;
* �}2.7.1-2.7.2 ;
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
  rowaxis values = (0 to 1 by 0.1) label = "�@"; 
  colaxis values = (0 to 30 by 10); 
  format stage stagef.;
run;
proc sgpanel data = OutS;
  panelby Stage / novarname;
  step x = LogTime y = LLS / group = Treat;
  keylegend / noborder;
  rowaxis values = (-4 to 2 by 1) label = "�@"; 
  colaxis values = (-1 to 4 by 1) label = "log(Time)"; 
  format stage stagef.;
run;
ods graphics off;
ods html close;

*==================================================*;
* �v���O����2.8.1 ;
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
* END ;
*==================================================*;
