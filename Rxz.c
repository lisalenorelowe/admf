#include "DeclareFunctions.h"

void Rxz(double ***Ricci24,Tensor gij,Tensor d2gxx,Tensor d2gyy,Tensor d2gzz,Tensor d2gxy,Tensor d2gxz,Tensor d2gyz,Vector dgxx,Vector dgyy, Vector dgzz, Vector dgxy, Vector dgxz, Vector dgyz,Params Par){

 int i,j,k;
 double t1,t2,t3,t4,t5,t6,t7,t8,t11,t12,t13,t14,t15,t16,t19,t20,t21,t23,t24,t25,t28,t32,t33,t34,t36,t40,t41,t43,t47,t48,t49,t53,t57,t58,t59,t62,t65,t66,t69,t71,t75,t76,t77,t78,t81,t82,t85,t86,t87,t90,t93,t96,t99,t100,t102,t103,t107,t112,t116,t120,t124,t128,t132,t137,t138,t139,t142,t143,t147,t150,t154,t155,t158,t159,t163,t167,t171,t172,t179,t183,t187,t188,t192,t193,t194,t199,t200,t203,t204,t207,t223,t224,t233,t236,t240,t241,t242,t244,t250,t252,t256,t263,t267,t273,t274,t280,t285,t290,t291,t295,t296,t299,t303,t309,t317,t323,t327,t332,t336,t337,t338,t341,t345,t355,t360,t363,t364,t369,t370,t371,t374,t375,t378,t381,t382,t385,t387,t400,t408,t411,t412,t416,t427,t431,t438,t439,t441,t443,t450,t457,t458,t466,t473,t480,t484,t485,t491,t495,t511,t522,t526,t529,t550,t551,t554,t574,t577,t580,t585,t588,t617,t625,t628,t633,t639,t646,t649,t657,t671,t673,t679,t683,t687,t691,t700,t705,t711,t714,t715,t728,t737,t741,t747,t756,t761,t762,t764,t766,t775,t779,t789,t804,t806,t815,t819,t823,t826,t844,t868,t879,t882,t886,t897,t911,t935,t951,t957,t983,t986,t1006,t1015,t1023,t1049,t1056,t1088,t1123,t1125,t1153,t1156,t1158,t1185,t1204,t1228,t1245,t1263,t1302,t1339,t1347,t1351;

double t9,t17,t18,t26,t27,t29,t35,t37,t39,t44,t54,t55,t56,t60,t67,t83,t84,t88,t91,t95,t97,t101,t109,t110,t113,t117;
double t118,t119,t125,t126,t127,t131,t136,t141,t144,t149,t156,t162,t166,t170,t173,t174,t177,t178,t181,t185,t189,t196;
double t205,t209,t212,t213,t217,t227,t228,t232,t235,t246,t257,t265,t269,t276,t279,t283,t286,t289,t300,t304,t308,t311;
double t312,t321,t325,t339,t342,t346,t357,t362,t366,t379,t380,t384,t389,t392,t395,t399,t402,t405,t409,t413,t419,t423;
double t426,t434,t442,t445,t449,t452,t455,t460,t463,t464,t468,t478,t492,t496,t504,t506,t508,t524,t538,t542,t545,t571;
double t586,t589,t626,t655,t660,t663,t668,t674,t677,t684,t695,t702,t709,t742,t758,t769,t777,t803,t808,t811,t813,t827,t830;
double t869,t876,t908,t914,t926,t931,t944,t947,t954,t969,t977,t984,t995,t1003,t1014,t1021,t1030,t1054,t1097,t1136;
double t1146,t1150,t1177,t1180,t1211,t1235,t1251,t1280,t1315,t1359,t1363;

double t22,t30,t31,t38,t42,t50,t51,t63,t72,t73,t74,t79,t80,t89,t92,t98,t104,t108,t115,t121,t129,t130,t133,t140;
double t146,t148,t151,t160,t164,t165,t168,t169,t180,t184,t186,t197,t202,t208,t216,t225,t229,t234,t239,t248,t249;
double t260,t275,t281,t288,t292,t302,t305,t306,t315,t319,t322,t349,t350,t354,t358,t359,t367,t368,t383,t388,t396;
double t403,t410,t414,t420,t421,t425,t429,t459,t482,t497,t498,t500,t523,t528,t532,t533,t536,t543,t544,t547,t548;
double t558,t562,t563,t575,t578,t581,t584,t587,t590,t591,t596,t600,t603,t606,t613,t632,t635,t638,t658,t662,t666;
double t680,t718,t725,t730,t740,t768,t772,t773,t784,t802,t814,t817,t825,t838,t859,t867,t895,t906,t912,t929,t932;
double t939,t949,t950,t959,t960,t974,t998,t1000,t1019,t1026,t1035,t1041,t1072,t1079,t1082,t1086,t1087,t1124,t1128;
double t1129,t1132,t1152,t1175,t1222,t1223,t1230,t1240,t1259,t1265,t1269,t1301,t1348,t1392,t1424,t1427,t1455,t1485;
double t1492,t1509,t1527,t1570,t1578,t1582;

for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
//> C(Rxz,optimized);
      t1 = gij.xx[i][j][k];
      t2 = t1*t1;
      t3 = gij.zz[i][j][k];
      t4 = t2*t3;
      t5 = gij.yy[i][j][k];
      t6 = dgzz.x[i][j][k];
      t7 = t5*t6;
      t8 = dgyz.y[i][j][k];
      t9 = t7*t8;
      t12 = t1*t3;
      t13 = gij.yz[i][j][k];
      t14 = t13*t13;
      t15 = dgxx.z[i][j][k];
      t16 = t14*t15;
      t17 = dgxy.y[i][j][k];
      t21 = gij.xz[i][j][k];
      t22 = t21*t21;
      t23 = t22*t5;
      t24 = t13*t15;
      t25 = dgyz.x[i][j][k];
      t29 = gij.xy[i][j][k];
      t30 = t29*t3;
      t31 = d2gxx.yz[i][j][k];
      t36 = t2*t13;
      t37 = d2gxy.zz[i][j][k];
      t38 = t37*t5;
      t42 = d2gxz.xy[i][j][k];
      t47 = t21*t5;
      t48 = t47*t29;
      t49 = dgxy.x[i][j][k];
      t50 = t3*t49;
      t51 = t50*t6;
      t54 = t13*t49;
      t55 = dgxz.z[i][j][k];
      t62 = dgxx.y[i][j][k];
      t63 = t13*t62;
      t67 = dgzz.z[i][j][k];
      t72 = t47*t1;
      t73 = dgxz.x[i][j][k];
      t74 = t13*t73;
      t75 = dgyz.z[i][j][k];
      t79 = t22*t29;
      t80 = dgyy.x[i][j][k];
      t85 = t29*t29;
      t86 = t21*t85;
      t87 = dgxz.y[i][j][k];
      t88 = t3*t87;
      t89 = t88*t25;
      t92 = dgxx.x[i][j][k];
      t93 = t13*t92;
      t97 = d2gxz.yz[i][j][k];
      t98 = t97*t5;
      t102 = -2.0*t4*t9+2.0*t12*t16*t17+2.0*t23*t24*t25+2.0*t30*t31*t22*t5-2.0*
t36*t38*t3-2.0*t30*t42*t22*t5-2.0*t48*t51-4.0*t23*t54*t55+2.0*t23*t54*t6+2.0*
t23*t63*t55-2.0*t23*t29*t67*t49+4.0*t72*t74*t75+2.0*t79*t13*t80*t55+2.0*t86*t89
-2.0*t23*t93*t75+2.0*t36*t98*t3;
      t103 = t14*t1;
      t104 = t21*t62;
      t108 = t1*t13;
      t109 = t108*t21;
      t110 = t29*t87;
      t115 = t5*t25*t55;
      t118 = t5*t5;
      t119 = t21*t118;
      t120 = t92*t3;
      t121 = t120*t55;
      t124 = t3*t15;
      t125 = dgxy.z[i][j][k];
      t126 = t124*t125;
      t129 = t21*t13;
      t130 = d2gyz.xx[i][j][k];
      t131 = t129*t130;
      t132 = t1*t5;
      t133 = t3*t132;
      t136 = dgyy.z[i][j][k];
      t140 = t21*t14;
      t141 = dgzz.y[i][j][k];
      t142 = t1*t141;
      t146 = t30*t21;
      t147 = t5*t73;
      t148 = t147*t87;
      t151 = t30*t1;
      t155 = t85*t3;
      t156 = t21*t73;
      t160 = t5*t141;
      t164 = t3*t62;
      t165 = t164*t75;
      t168 = t29*t13;
      t169 = t22*t49;
      t173 = t29*t14;
      t178 = d2gyy.xz[i][j][k];
      t180 = t29*t21;
      t181 = t180*t13;
      t184 = -2.0*t103*t104*t75-4.0*t109*t110*t75+2.0*t109*t115-2.0*t121*t119
-3.0*t48*t126+2.0*t131*t133+2.0*t109*t50*t136-2.0*t140*t142*t49+2.0*t146*t148+
2.0*t151*t74*t136-2.0*t155*t156*t136-2.0*t151*t160*t73+2.0*t72*t165-4.0*t168*
t169*t75+2.0*t173*t21*t92*t75-4.0*t12*t178*t181;
      t186 = d2gxy.yz[i][j][k];
      t192 = t3*t3;
      t193 = t29*t192;
      t194 = t5*t15;
      t197 = t85*t6;
      t202 = t85*t141;
      t208 = t29*t73;
      t212 = t29*t141;
      t216 = t130*t1;
      t224 = t22*t13;
      t225 = t87*t87;
      t229 = d2gxz.xz[i][j][k];
      t234 = d2gxx.zz[i][j][k];
      t239 = d2gzz.xx[i][j][k];
      t240 = t239*t1;
      t248 = 4.0*t12*t186*t181-t151*t24*t136+t193*t194*t62-t129*t197*t25+t129*
t197*t87+t129*t202*t15-2.0*t72*t24*t75-4.0*t23*t208*t75+2.0*t23*t212*t73+2.0*
t30*t216*t14+2.0*t30*t130*t22*t5+2.0*t224*t29*t225-4.0*t47*t229*t85*t3+2.0*t47*
t234*t85*t3+2.0*t47*t240*t14+2.0*t47*t239*t85*t3;
      t249 = t12*t21;
      t250 = t29*t25;
      t257 = t5*t92;
      t260 = t1*t80;
      t265 = t13*t25;
      t269 = t13*t87;
      t275 = dgyy.y[i][j][k];
      t279 = t180*t1;
      t280 = t3*t136;
      t281 = t280*t87;
      t283 = t22*t80;
      t288 = t1*t67;
      t292 = t129*t42;
      t295 = d2gxy.xz[i][j][k];
      t296 = t129*t295;

      t299 = -2.0*t249*t250*t8+2.0*t249*t110*t8+t193*t257*t25-t193*t260*t125+
t193*t260*t25-2.0*t249*t265*t17+2.0*t249*t269*t17+t23*t164*t25-t79*t3*t275*t15+
t279*t281+t30*t283*t125-t30*t283*t25+t119*t288*t15-t36*t281-2.0*t292*t133-2.0*
t296*t133;
      t302 = t129*t31;
      t305 = t13*t125;
      t306 = t305*t8;
      t309 = t265*t8;
      t312 = t269*t8;
      t315 = t29*t125;
      t319 = t14*t92;
      t322 = t1*t192;
      t323 = t5*t62;
      t337 = t160*t25;
      t339 = t21*t15;
      t342 = t6*t6;
      t349 = 2.0*t302*t133+2.0*t4*t306+2.0*t4*t309-2.0*t4*t312-2.0*t249*t315*t8
-t47*t319*t6+t322*t323*t125+t322*t323*t87-t322*t323*t25+t322*t29*t275*t15-t30*
t16*t62+t322*t194*t80+t4*t337-t173*t339*t125-t168*t132*t342-2.0*t30*t295*t22*t5
;
      t350 = t31*t1;
      t354 = t85*t13;
      t358 = t168*t1;
      t359 = t5*t67;
      t363 = t229*t1;
      t367 = d2gyz.xz[i][j][k];
      t368 = t367*t5;
      t382 = t22*t21;
      t383 = t382*t5;
      t384 = t25*t25;
      t387 = t85*t29;
      t388 = t387*t192;
      t395 = t3*t80;
      t396 = t395*t55;
      t399 = t21*t87;
      t403 = t125*t125;
      t410 = t295*t1;
      t414 = 2.0*t30*t350*t14-2.0*t354*t339*t75+2.0*t358*t359*t73-4.0*t47*t363*
t14+2.0*t36*t368*t3-4.0*t146*t54*t87+4.0*t146*t54*t25-2.0*t193*t1*t49*t136-2.0*
t383*t384+2.0*t388*t130-2.0*t103*t21*t80*t55+2.0*t358*t396+2.0*t103*t399*t25+
2.0*t12*t47*t403+2.0*t146*t24*t80-2.0*t30*t410*t14;
      t420 = t168*t5;
      t421 = t120*t6;
      t423 = t7*t87;
      t425 = t359*t25;
      t427 = t359*t125;
      t429 = t359*t87;
      t431 = t7*t141;
      t441 = t29*t6;
      t459 = -t193*t257*t125-t193*t260*t87-t420*t421+t151*t423-t36*t425+t36*
t427-t36*t429+t36*t431-2.0*t79*t3*t25*t17+2.0*t72*t13*t67*t49+2.0*t23*t441*t125
-t23*t164*t125-t358*t359*t15-2.0*t193*t216*t5+2.0*t224*t208*t136-2.0*t108*t367*
t22*t5;
      t482 = t367*t29*t21;
      t497 = t22*t1;
      t498 = t7*t136;
      t500 = t160*t125;
      t508 = -2.0*t108*t367*t85*t3-2.0*t108*t97*t22*t5-2.0*t108*t97*t85*t3+4.0*
t181*t147*t6+4.0*t181*t194*t55-t249*t13*t275*t15+6.0*t103*t482+2.0*t354*t126
-4.0*t129*t202*t73-t193*t257*t87+t279*t425+t140*t441*t62-t249*t194*t136-t497*
t498+t497*t500+2.0*t181*t194*t6-4.0*t168*t229*t133;
      t523 = t124*t25;
      t526 = t280*t125;
      t528 = t22*t136;
      t532 = t13*t6;
      t533 = t532*t17;
      t536 = t21*t141;
      t543 = t22*t14;
      t544 = t49*t25;
      t547 = t14*t13;
      t548 = t21*t547;
      t551 = t42*t1;
      t554 = t382*t13;
      t558 = t21*t387;
      t562 = 2.0*t168*t234*t133-6.0*t155*t131-2.0*t12*t283*t136-2.0*t86*t305*
t55+2.0*t354*t523-t279*t526-t108*t528*t125+t36*t526-2.0*t151*t533-2.0*t155*t536
*t49-2.0*t155*t74*t125-2.0*t543*t544+2.0*t548*t410+2.0*t548*t551-2.0*t554*t49*
t136-2.0*t558*t37*t3;
      t563 = t382*t29;
      t574 = t22*t85;
      t575 = t25*t75;
      t578 = t125*t75;
      t581 = t87*t75;
      t584 = t22*t178;
      t587 = t22*t186;
      t590 = d2gyz.xy[i][j][k];
      t591 = t22*t590;
      t596 = t387*t13;
      t600 = t29*t547;
      t603 = t234*t1;
      t606 = t85*t14;
      t613 = 2.0*t563*t368+2.0*t558*t367*t3+2.0*t563*t98+2.0*t558*t97*t3+2.0*
t574*t575-2.0*t574*t578+2.0*t574*t581-2.0*t584*t103+2.0*t587*t155+2.0*t591*t103
+2.0*t591*t155-2.0*t596*t67*t73+4.0*t600*t363-2.0*t600*t603-8.0*t606*t229*t21+
4.0*t606*t234*t21;
      t632 = t73*t55;
      t635 = t73*t6;
      t638 = t15*t55;
      t658 = t2*t192;
      t662 = 4.0*t606*t239*t21+4.0*t596*t229*t3-2.0*t596*t234*t3-2.0*t600*t240
-2.0*t596*t239*t3-2.0*t600*t92*t55+4.0*t606*t632-2.0*t606*t635-2.0*t606*t638+
4.0*t383*t49*t75-2.0*t383*t62*t75-2.0*t584*t155-2.0*t4*t186*t14-2.0*t383*t141*
t49+2.0*t383*t87*t25-2.0*t658*t178*t5;
      t666 = d2gxz.yy[i][j][k];
      t674 = t5*t125*t55;
      t680 = t29*t15;
      t684 = t7*t25;
      t715 = 2.0*t658*t590*t5+2.0*t4*t666*t14+2.0*t322*t666*t85+2.0*t79*t674
-3.0*t224*t441*t80-3.0*t224*t680*t136-3.0*t109*t684-3.0*t129*t197*t125-2.0*t146
*t93*t136+2.0*t193*t551*t5+6.0*t155*t292+2.0*t193*t410*t5+4.0*t224*t212*t49+2.0
*t140*t208*t125+2.0*t140*t208*t87-2.0*t140*t208*t25+2.0*t354*t51;
      t718 = t21*t49;
      t725 = t129*t5;
      t730 = t1*t73;
      t740 = d2gzz.xy[i][j][k];
      t768 = t108*t5;
      t769 = t164*t55;
      t772 = 4.0*t173*t718*t55-2.0*t173*t718*t6-t725*t120*t25+t109*t395*t125
-4.0*t173*t730*t75+2.0*t173*t142*t73-2.0*t103*t110*t55+2.0*t108*t740*t22*t5+2.0
*t12*t197*t8-4.0*t12*t587*t5+2.0*t168*t239*t133+2.0*t420*t121+12.0*t168*t229*
t22*t5-6.0*t168*t234*t22*t5-6.0*t168*t239*t22*t5-2.0*t768*t769;
      t773 = t740*t5;
      t784 = t5*t87*t55;
      t802 = t1*t15;
      t814 = t85*t67;
      t817 = -2.0*t36*t773*t3+2.0*t108*t283*t75+2.0*t103*t315*t55+2.0*t109*t784
-t108*t197*t141-t279*t427+t279*t429+2.0*t23*t680*t75+2.0*t108*t37*t22*t5+2.0*
t108*t37*t85*t3+2.0*t140*t802*t136-4.0*t173*t104*t55+2.0*t79*t88*t17-t768*t523+
t768*t126-t108*t814*t125;
      t825 = t160*t87;
      t838 = t24*t8;
      t859 = t108*t814*t25+t108*t814*t87+t108*t528*t87+t4*t825-t109*t395*t25-
t12*t197*t136+t140*t680*t87+2.0*t497*t9+2.0*t23*t124*t17+2.0*t79*t838-2.0*t79*
t3*t125*t17+6.0*t155*t296-2.0*t193*t350*t5-6.0*t155*t302-2.0*t109*t29*t80*t67
-4.0*t12*t666*t181;
      t867 = t120*t75;
      t895 = t147*t25;
      t906 = -4.0*t168*t22*t87*t25-3.0*t173*t339*t25-2.0*t354*t867-2.0*t155*t74
*t87+2.0*t155*t74*t25-t47*t814*t15+t768*t3*t6*t62+t768*t124*t87+2.0*t173*t802*
t75+4.0*t354*t156*t75-2.0*t497*t309+2.0*t497*t312+2.0*t79*t533-2.0*t146*t895+
2.0*t155*t536*t62-2.0*t497*t306+2.0*t322*t315*t17;
      t912 = t740*t29*t21;
      t929 = t147*t125;
      t932 = t666*t22;
      t939 = t22*t275;
      t949 = t22*t22;
      t950 = t949*t275;
      t959 = -6.0*t103*t912-2.0*t79*t115-2.0*t249*t160*t62+2.0*t322*t250*t17
-2.0*t322*t110*t17-2.0*t249*t305*t17+2.0*t146*t929+4.0*t12*t932*t5-4.0*t12*t591
*t5+2.0*t12*t939*t25-2.0*t12*t939*t87+2.0*t12*t47*t384+t950*t87-2.0*t383*t403
-2.0*t72*t88*t125-2.0*t72*t89;
      t960 = t50*t75;
      t974 = t7*t125;
      t995 = t180*t37;
      t998 = -4.0*t72*t960+2.0*t72*t3*t141*t49-2.0*t103*t250*t55+t725*t120*t87+
t109*t395*t87-t109*t974-t725*t124*t62-2.0*t119*t603*t3-2.0*t119*t240*t3-t109*
t423+2.0*t912*t133-2.0*t79*t784+2.0*t86*t265*t55-t279*t431-2.0*t109*t674+2.0*
t995*t133;
      t1000 = t180*t97;
      t1019 = t532*t275;
      t1021 = t160*t15;
      t1023 = t265*t136;
      t1026 = t395*t75;
      t1035 = -2.0*t1000*t133-2.0*t482*t133-t497*t825-t23*t124*t80-t23*t164*t87
-4.0*t224*t212*t62-t12*t202*t87+t12*t202*t125-t249*t7*t80+t4*t1019-t109*t1021+
t497*t1023-t497*t1019-2.0*t36*t1026-4.0*t109*t250*t75+4.0*t109*t315*t75;
      t1041 = t14*t62;
      t1072 = t22*t118;
      t1079 = 2.0*t109*t441*t136-t12*t16*t80-t12*t1041*t125-t12*t1041*t87+t12*
t1041*t25-t497*t337+2.0*t108*t740*t85*t3+2.0*t563*t125*t8-2.0*t563*t87*t8-2.0*
t383*t6*t17-2.0*t932*t103-2.0*t932*t155+2.0*t658*t186*t5+2.0*t383*t87*t125+4.0*
t1072*t632-2.0*t1072*t635-2.0*t1072*t638;
      t1082 = t1*t547;
      t1086 = t2*t14;
      t1087 = t141*t25;
      t1124 = t85*t192;
      t1125 = t62*t125;
      t1128 = 2.0*t1082*t62*t55-2.0*t1086*t1087-2.0*t1086*t578+2.0*t1086*t581+
2.0*t1086*t575-4.0*t574*t367*t13+4.0*t574*t37*t13+2.0*t224*t29*t403-4.0*t574*
t97*t13+4.0*t574*t740*t13-2.0*t563*t80*t75-2.0*t563*t773-2.0*t558*t740*t3-2.0*
t563*t38-2.0*t658*t666*t5-2.0*t1124*t1125;
      t1129 = t49*t87;
      t1132 = t49*t125;
      t1152 = t387*t3;
      t1175 = 2.0*t1124*t1129+2.0*t1124*t1132-2.0*t1124*t544-2.0*t155*t21*t225+
2.0*t554*t42*t5-2.0*t554*t31*t5-2.0*t140*t1*t225+2.0*t224*t29*t384+2.0*t1152*
t141*t73-2.0*t548*t216-2.0*t554*t130*t5+2.0*t554*t295*t5-2.0*t548*t350-4.0*t543
*t42*t29-4.0*t543*t295*t29+4.0*t543*t31*t29;
      t1222 = 4.0*t543*t130*t29+2.0*t543*t1132+2.0*t543*t1129+2.0*t563*t25*t8+
4.0*t168*t22*t62*t75+6.0*t103*t1000+2.0*t12*t939*t125+t30*t283*t87-4.0*t224*
t110*t125+2.0*t249*t7*t17+2.0*t155*t21*t6*t80+2.0*t155*t339*t136-2.0*t86*t165+
4.0*t12*t590*t181-2.0*t322*t194*t17+2.0*t47*t319*t55;
      t1223 = t50*t55;
      t1230 = t1*t6;
      t1240 = t92*t136;
      t1259 = 4.0*t48*t1223-6.0*t103*t995+4.0*t358*t960+2.0*t173*t1230*t125-
t1152*t141*t15-t1124*t15*t80+t558*t67*t125+t1124*t1240-t554*t80*t87-t548*t92*
t87+t554*t80*t25+t543*t1240+t548*t15*t62+t558*t6*t141+t543*t62*t87-t543*t62*t25
+t596*t67*t15;
      t1263 = t15*t6;
      t1265 = t92*t67;
      t1269 = t15*t15;
      t1301 = -t606*t1263+t606*t1265+t600*t92*t6+t119*t3*t1269-t47*t85*t342-
t658*t275*t125-t658*t275*t25+t658*t275*t87+t658*t80*t136+t548*t92*t25-t554*t80*
t125+2.0*t4*t178*t14+2.0*t322*t178*t85-2.0*t322*t186*t85-2.0*t4*t590*t14-2.0*
t322*t590*t85;
      t1348 = 2.0*t587*t103-2.0*t543*t15*t17-2.0*t574*t6*t8+4.0*t382*t666*t168+
4.0*t382*t178*t168-4.0*t382*t186*t168-4.0*t382*t590*t168-2.0*t383*t15*t8+2.0*
t383*t141*t62+2.0*t554*t125*t17+2.0*t554*t25*t17-2.0*t554*t87*t17-2.0*t119*t288
*t73-2.0*t86*t396+2.0*t23*t441*t25+2.0*t86*t269*t55;
      t1392 = 2.0*t30*t108*t225-2.0*t30*t108*t403+4.0*t146*t63*t125+t151*t684+
t151*t1021+t30*t319*t87+2.0*t173*t1230*t25-2.0*t173*t288*t49+2.0*t354*t21*t67*
t49-8.0*t181*t147*t55+2.0*t140*t1*t87*t125-2.0*t151*t838+2.0*t249*t194*t8-2.0*
t224*t929+4.0*t12*t584*t5+2.0*t279*t1026;
      t1424 = -t548*t92*t125-t574*t1087+t574*t141*t87-t574*t141*t125+t554*t275*
t15+t383*t15*t136-t563*t25*t136+t383*t6*t80+t563*t6*t275+t119*t1*t342-t1072*
t1263-t47*t14*t1269+t1072*t1265-t1082*t15*t125-t1082*t6*t62+t1082*t15*t25-t1082
*t15*t87;
      t1427 = t80*t67;
      t1455 = t1086*t1427-t1086*t6*t136-t543*t1125+t563*t136*t125-t558*t67*t25-
t558*t67*t87-t563*t136*t87+t574*t1427+t1152*t6*t125-t1152*t6*t25-t1152*t6*t87-
t4*t1023-t151*t974+t249*t250*t136-t146*t7*t62-t151*t532*t80;
      t1485 = t2*t547;
      t1492 = -2.0*t181*t257*t67+2.0*t354*t769-4.0*t354*t1223+2.0*t155*t399*
t125-2.0*t30*t108*t384+2.0*t47*t814*t73-4.0*t146*t54*t125-t249*t441*t275+t119*
t421-t48*t523+t725*t120*t125-t950*t25-2.0*t949*t666*t5+2.0*t1485*t37+2.0*t388*
t31-2.0*t388*t42;
      t1509 = t382*t118;
      t1527 = 2.0*t949*t590*t5+2.0*t949*t186*t5+t949*t80*t136-t950*t125-2.0*
t1485*t367-2.0*t388*t295+t600*t1269+t596*t342+2.0*t1509*t239-4.0*t1509*t229-2.0
*t949*t178*t5+2.0*t1485*t740-2.0*t1485*t97+2.0*t1509*t234-2.0*t224*t148+2.0*
t224*t895;
      t1570 = 2.0*t140*t1230*t80-2.0*t140*t730*t136-2.0*t30*t551*t14+2.0*t30*
t169*t136+4.0*t109*t212*t25+t30*t319*t125-t146*t194*t87-t30*t319*t25+2.0*t140*
t142*t62+2.0*t23*t24*t125+2.0*t48*t867+2.0*t47*t603*t14+4.0*t119*t363*t3-t168*
t5*t3*t1269+t4*t498-t12*t202*t25-t4*t500;
      t1578 = pow(t133-t103-t155+2.0*t181-t23,2.0);
      t1582 = -(t1492+t998+t1175+t1035+t1570+t414+t1301+t1527+t1259+t562+t1128+
t772+t508+t959+t349+t906+t1392+t1455+t817+t459+t1424+t662+t859+t299+t1222+t248+
t102+t1079+t184+t1348+t715+t613)/t1578/4.0;

 Ricci24[i][j][k] = -t1582;

}}}

return;

}
