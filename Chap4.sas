*======================================================================*;
/*

大橋靖雄・浜田知久馬・魚住龍史 (2022)
生存時間解析 第2版 --SASによる生物統計
東京大学出版会

第4章　生存時間解析に関連した他のプロシジャ

Program : Chap4.sas
Made by : 魚住龍史 (e-mail : uozumi at kuhp.kyoto-u.ac.jp)
Date : 2022年8月

*/
*======================================================================*;

*==================================================*;
* プログラム4.1.1 ;
*==================================================*;

proc lifereg data = Gehan;
  model Week*Remiss(0) = Drug / dist = weibull;
run;

*==================================================*;
* プログラム4.1.2 ;
*==================================================*;

proc nlmixed data = Gehan;
  parms lambda = 0.02 gamma = 1 beta = 0;
  bounds lambda > 0, gamma > 0;
  expxb = exp(-beta*Drug);
  G_t = exp(-lambda*(Week*expxb)**gamma);
  g = expxb*gamma*lambda*(Week*expxb)**(gamma-1)*G_t;
  LogL = (Remiss = 1)*log(g) + (Remiss = 0)*log(G_t);
  model Week ~ general(LogL);
run;

*==================================================*;
* 図4.1.3 ;
*==================================================*;

ods output ParameterEstimates = Est00;
proc lifereg data = Gehan;
  model Week*Remiss(0) = Drug / dist=weibull;
run;

data Est00; set Est00; keep Parameter Estimate;
proc transpose data = Est00 out = Est00;
  id Parameter; var Estimate;
data Est00; set Est00; Lambda=exp(-Intercept/Scale); Gamma=1/Scale;
run;

data Surv00; set Est00;
  do t = 0 to 35 by 0.01;
    logt = log(t);
    do g = 0;
      S = exp(-lambda*(t*exp(-Drug*g))**Gamma); LLS = log(-log(S)); output;
    end;
    do g = 1;
      S = exp(-lambda*(t*exp(-Drug*g))**Gamma); LLS = log(-log(S)); output;
    end;
end;
run;

data Surv00; set Surv00;
  if round(t, 1e-2) = 7.4 then do; if g = 0 then y = LLS; end;
run;

ods html image_dpi=400 style=journal;
ods graphics / reset noborder noscale
               width=600px height=400px
               imagefmt=png
			   antialiasmax=7100
;
proc sgplot data = Surv00 noautolegend;
  series x = t y = S / group = g;
  yaxis values = (0 to 1 by 0.1) label = "生存関数の推定値"; 
  xaxis values = (0 to 35 by 5) label = "Week"; 
run;
proc sgplot data = Surv00 noautolegend;
  series x = logt y = LLS / group = g;
  yaxis values = (-5 to 2 by 1) label = "2重対数プロット"; 
  xaxis values = (0 to 4 by 1) label = "log(Week)"; 
  vector x = logt y = y / xorigin = 2.0014800002 yorigin = -2.068039496 arrowdirection = both;
  inset "2群間の垂直方向の距離：β・γ = β/σ = 1.7309"
        / position = bottomright textattrs = GraphAnnoText
  ;
run;
ods graphics off;
ods html close;

proc lifereg data = TD50 outest = Est(drop = NAME TYPE MODEL SHAPE1);
  model (TimeS, TimeE) = LogDose / d = weibull scale=0.33333 noscale;
  output out = Out q = 0.5 p = Quantile cdf = CDF std = STD;
  weight W; by Tumor;
run;
data Est; set Est; 
  Beta = LogDose; Gamma = 1/_SCALE_; Lambda = exp(-Intercept*Gamma);
  TD50 = (log(2)/(Lambda*104**3))**(-1/Beta/3);
  keep Tumor Lambda Gamma Beta TD50;
run;
proc print data = Est; run;

*==================================================*;
* プログラム4.1.4 ;
*==================================================*;

data Gast2; input C1 C2 Group Time TimeS TimeE W @@;
cards;
0 0 0 0 . 2 15   0 0 0 0 2 . 24   0 0 0 1 . 2 17   0 0 0 1 2 . 7
0 0 1 0 . 2 15   0 0 1 0 2 . 19   0 0 1 1 . 2 17   0 0 1 1 2 . 2
1 0 0 0 . 2 12   1 0 0 0 2 . 28   1 0 0 1 . 2 13   1 0 0 1 2 . 15
1 0 1 0 . 2 17   1 0 1 0 2 . 27   1 0 1 1 . 2 17   1 0 1 1 2 . 10
0 1 0 0 . 2 3    0 1 0 0 2 . 35   0 1 0 1 . 2 17   0 1 0 1 2 . 18
0 1 1 0 . 2 7    0 1 1 0 2 . 33   0 1 1 1 . 2 17   0 1 1 1 2 . 16
;
proc lifereg data = Gast2;
  model (TimeS, TimeE) = Group Time C1 C2 / d = exponential;
  weight W;
run;

*==================================================*;
* プログラム4.2.1 ;
*==================================================*;

data Gast3;
input C1 C2 Group Time Y R @@;
R1 = R*10000; R2 = R*100000;
cards;
0 0 0 0 15 63   0 0 0 1 17 31   0 0 1 0 15 53   0 0 1 1 17 21
1 0 0 0 12 68   1 0 0 1 13 43   1 0 1 0 17 71   1 0 1 1 17 37
0 1 0 0  3 73   0 1 0 1 17 53   0 1 1 0  7 73   0 1 1 1 17 49
;
/* LOGISTICプロシジャ (Rを10,000倍) */
proc logistic data = Gast3;
  model Y / R1 = Group Time C1 C2;
run;
/* LOGISTICプロシジャ (Rを100,000倍) */
proc logistic data = Gast3;
  model Y / R2 = Group Time C1 C2 / expb stb;
run;
/* PROBITプロシジャ (Rを10,000倍) */
proc probit data = Gast3;
  model Y / R1 = Group Time C1 C2 / d = logistic;
run;

*==================================================*;
* プログラム4.2.2 ;
*==================================================*;

data Gast4; set Gast3; LogR = log(R);
proc genmod data = Gast4;
  model Y = Group Time C1 C2 / d = poisson offset = LogR;
run;

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

*==================================================*;
* プログラム4.3.1 ;
*==================================================*;

proc phreg data = Work;
  model Time*Censor(0) = Group;
run;

*==================================================*;
* プログラム4.3.2 ;
*==================================================*;

data Data;
  input Time Group $ Death @@;
cards;
1 A 1   1 A 0   1 A 0   1 B 0   1 B 0   1 B 0
2 A 1   2 A 0   2 B 0   2 B 0   2 B 0
3 A 0   3 B 1   3 B 0   3 B 0
4 A 1   4 B 0   4 B 0
5 B 1   5 B 0
6 B 1
;
proc freq data = Data;
  tables Time*Death*Group / chisq cmh norow nocol nopercent;
run;

*==================================================*;
* END ;
*==================================================*;
