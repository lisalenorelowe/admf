#include "DeclareFunctions.h"

void Rzz(double ***Ricci44,Tensor gij,Tensor d2gxx,Tensor d2gyy,Tensor d2gzz,Tensor d2gxy,Tensor d2gxz,Tensor d2gyz,Vector dgxx,Vector dgyy, Vector dgzz, Vector dgxy, Vector dgxz, Vector dgyz,Params Par){

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

for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){

//> C(Rzz,optimized);
      t1 = gij.yy[i][j][k];
      t2 = gij.zz[i][j][k];
      t3 = t1*t2;
      t4 = gij.xy[i][j][k];
      t5 = t3*t4;
      t6 = gij.yz[i][j][k];
      t7 = dgxx.z[i][j][k];
      t8 = t6*t7;
      t9 = dgxz.z[i][j][k];
      t13 = gij.xz[i][j][k];
      t14 = t13*t13;
      t15 = t4*t4;
      t16 = t14*t15;
      t17 = dgzz.y[i][j][k];
      t18 = t17*t17;
      t20 = dgzz.z[i][j][k];
      t21 = t15*t20;
      t24 = t6*t6;
      t25 = t13*t24;
      t26 = gij.xx[i][j][k];
      t27 = dgxy.z[i][j][k];
      t29 = dgyz.z[i][j][k];
      t33 = t26*t26;
      t34 = t2*t2;
      t35 = t33*t34;
      t36 = dgyy.z[i][j][k];
      t37 = t36*t36;
      t39 = t24*t14;
      t40 = t27*t27;
      t43 = t14*t14;
      t44 = d2gyz.yz[i][j][k];
      t48 = t24*t24;
      t49 = d2gxz.xz[i][j][k];
      t53 = t4*t2;
      t54 = t6*t13;
      t55 = dgyz.x[i][j][k];
      t56 = t55*t55;
      t60 = t14*t6;
      t62 = t1*t27*t9;
      t65 = t54*t26;
      t66 = dgyy.x[i][j][k];
      t67 = t2*t66;
      t71 = t15*t2;
      t76 = t1*t1;
      t77 = t76*t2;
      t78 = t13*t7;
      t82 = 2.0*t5*t8*t9+t16*t18-t3*t21*t7-4.0*t25*t26*t27*t29+t35*t37+4.0*t39*
t40+4.0*t43*t44*t1+4.0*t48*t49*t26+4.0*t53*t54*t56-4.0*t60*t62-2.0*t65*t67*t29
-4.0*t71*t13*t17*t55-2.0*t77*t78*t9;
      t83 = dgxz.y[i][j][k];
      t84 = t83*t83;
      t88 = t15*t34;
      t91 = dgxz.x[i][j][k];
      t95 = t3*t13;
      t96 = dgxy.x[i][j][k];
      t97 = t6*t96;
      t101 = t26*t83;
      t109 = t1*t20;
      t110 = t109*t83;
      t113 = t4*t83;
      t117 = t4*t20;
      t118 = dgxx.y[i][j][k];
      t119 = t117*t118;
      t125 = t24*t4;
      t126 = t2*t96;
      t127 = dgzz.x[i][j][k];
      t131 = d2gxy.zz[i][j][k];
      t132 = t131*t26;
      t136 = d2gyz.xz[i][j][k];
      t141 = 4.0*t53*t54*t84-2.0*t88*t56+2.0*t3*t21*t91-4.0*t95*t97*t9+4.0*t25*
t101*t29+4.0*t25*t26*t55*t29-2.0*t65*t110-4.0*t60*t113*t29-t95*t119+4.0*t65*t2*
t36*t27+2.0*t125*t126*t127-4.0*t53*t132*t24+4.0*t53*t136*t14*t1;
      t143 = dgyy.y[i][j][k];
      t144 = t43*t143;
      t149 = t26*t17;
      t155 = t24*t15;
      t156 = t127*t127;
      t158 = t1*t34;
      t159 = t4*t127;
      t162 = t4*t7;
      t163 = t162*t29;
      t166 = t4*t17;
      t167 = t166*t91;
      t170 = dgxx.x[i][j][k];
      t171 = t48*t170;
      t173 = t76*t34;
      t174 = t7*t7;
      t177 = t15*t4;
      t178 = t177*t34;
      t181 = d2gzz.yy[i][j][k];
      t185 = -2.0*t144*t29-2.0*t88*t84-4.0*t25*t149*t55+2.0*t88*t40+t155*t156-
t158*t159*t118+2.0*t95*t163+2.0*t95*t167+t171*t127+t173*t174+t144*t17-4.0*t178*
t131-2.0*t43*t181*t1;
      t189 = t109*t55;
      t192 = t26*t2;
      t193 = t24*t66;
      t196 = d2gxz.yz[i][j][k];
      t199 = t159*t83;
      t203 = t26*t1*t2;
      t204 = t6*t20;
      t205 = t204*t118;
      t209 = d2gzz.xy[i][j][k];
      t212 = t192*t13;
      t213 = t1*t17;
      t217 = d2gxx.zz[i][j][k];
      t223 = d2gzz.xx[i][j][k];
      t227 = t14*t13;
      t228 = t227*t4;
      t232 = -4.0*t125*t126*t9-2.0*t65*t189+t192*t193*t127+4.0*t178*t196+6.0*
t25*t199-t203*t205+4.0*t178*t136-4.0*t178*t209-2.0*t212*t213*t27-2.0*t48*t217*
t26-2.0*t171*t9-2.0*t48*t223*t26-2.0*t228*t36*t29;
      t235 = t227*t1;
      t236 = dgyz.y[i][j][k];
      t246 = t227*t6;
      t250 = t26*t34;
      t257 = d2gyy.zz[i][j][k];
      t265 = t14*t66;
      t269 = t1*t55*t9;
      t276 = t4*t91*t29;
      t279 = t1*t7;
      t283 = -4.0*t235*t236*t9+2.0*t235*t236*t127+2.0*t235*t36*t9+2.0*t246*t143
*t9-4.0*t250*t44*t15-2.0*t173*t170*t9-2.0*t43*t257*t1+4.0*t53*t14*t36*t27+t53*

t265*t17+4.0*t60*t269+2.0*t192*t21*t236+4.0*t25*t276+2.0*t25*t279*t9;
      t285 = t4*t13;
      t286 = t285*t6;
      t289 = t4*t36;
      t299 = t53*t13;
      t300 = t6*t17;
      t304 = t209*t26;
      t308 = t6*t170;
      t311 = dgxy.y[i][j][k];
      t312 = t4*t311;
      t321 = t44*t14;
      t325 = t4*t143;
      t332 = -4.0*t3*t223*t286+2.0*t212*t289*t29-2.0*t71*t13*t127*t36-t192*t21*
t36-4.0*t299*t300*t118-4.0*t53*t304*t24+t95*t308*t17+4.0*t250*t312*t29-2.0*t250
*t312*t17-4.0*t299*t269-8.0*t192*t321*t1+2.0*t250*t325*t9+2.0*t158*t149*t96;
      t337 = t33*t2;
      t338 = t6*t236;
      t339 = t338*t29;
      t342 = t26*t66;
      t346 = t4*t96;
      t357 = t338*t17;
      t362 = t192*t4;
      t363 = t300*t27;
      t366 = t1*t66;
      t369 = t166*t7;
      t374 = t4*t34;
      t375 = t136*t26;
      t379 = -4.0*t299*t8*t36+4.0*t337*t339+2.0*t25*t342*t20-2.0*t158*t346*t127
-4.0*t158*t26*t96*t29-4.0*t158*t101*t55-2.0*t337*t357+t212*t289*t17+2.0*t362*
t363-t250*t366*t127+t95*t369-4.0*t53*t54*t40-4.0*t374*t375*t1;
      t380 = t6*t127;
      t381 = t380*t83;
      t384 = t14*t1;
      t385 = t2*t311;
      t389 = t338*t9;
      t392 = t109*t236;
      t395 = t6*t66;
      t399 = t14*t26;
      t402 = t1*t127;
      t405 = t1*t91;
      t409 = t159*t55;
      t412 = t6*t36;
      t413 = t412*t29;
      t416 = t54*t209;
      t419 = t4*t236;
      t423 = t117*t96;
      t426 = -4.0*t71*t381+4.0*t384*t385*t9-4.0*t362*t389-2.0*t337*t392+8.0*
t299*t395*t9+2.0*t399*t392+t212*t402*t36-4.0*t25*t405*t9-2.0*t25*t409-2.0*t337*
t413-4.0*t416*t203-4.0*t212*t419*t29-2.0*t95*t423;
      t434 = t6*t91;
      t441 = t24*t6;
      t442 = t441*t26;
      t445 = t1*t18;
      t449 = t66*t127;
      t452 = t441*t13;
      t455 = t17*t118;
      t457 = t26*t24;
      t460 = t441*t4;
      t463 = -t212*t325*t20+2.0*t25*t423+2.0*t5*t434*t127-4.0*t5*t434*t9+t442*
t20*t118-t399*t445+t235*t66*t20+t39*t449+t337*t445-t452*t127*t118+t39*t455-t457
*t1*t156-t460*t7*t127;
      t464 = t223*t26;
      t468 = t2*t17;
      t478 = t1*t236;
      t485 = t4*t55;
      t492 = t6*t311;
      t496 = t54*t136;
      t504 = t109*t36;
      t506 = 4.0*t3*t464*t24-2.0*t457*t468*t96+4.0*t457*t126*t29+4.0*t158*t346*
t9+4.0*t212*t478*t9-2.0*t212*t213*t83-4.0*t60*t485*t29+2.0*t54*t21*t83-4.0*t212
*t492*t29+4.0*t496*t203+2.0*t203*t381+2.0*t212*t213*t55+t337*t504;
      t508 = t338*t127;
      t524 = t66*t9;
      t529 = t7*t36;
      t538 = t196*t26;
      t542 = t109*t27;
      t545 = 2.0*t362*t508-2.0*t54*t15*t127*t17+8.0*t39*t136*t4-8.0*t39*t131*t4
+8.0*t39*t196*t4-2.0*t39*t524+2.0*t88*t455+2.0*t88*t529+2.0*t95*t199+2.0*t95*
t409-2.0*t25*t163+4.0*t53*t538*t24+2.0*t65*t542;
      t551 = t14*t4;
      t571 = t300*t83;
      t574 = t6*t27;
      t586 = -2.0*t25*t167+4.0*t53*t375*t24-2.0*t551*t363+2.0*t212*t109*t311
-4.0*t374*t538*t1+t125*t2*t127*t118+2.0*t374*t342*t29+2.0*t362*t189-2.0*t60*
t213*t7-2.0*t551*t571-4.0*t71*t574*t9+2.0*t71*t13*t66*t20-2.0*t212*t1*t36*t9;
      t589 = t24*t311;
      t626 = 4.0*t192*t589*t9-12.0*t71*t496-2.0*t3*t24*t174-4.0*t77*t49*t14+2.0
*t77*t223*t14-2.0*t3*t14*t84+4.0*t246*t131*t1+4.0*t452*t132-4.0*t246*t136*t1
-4.0*t452*t375-4.0*t246*t196*t1-4.0*t452*t7*t27+2.0*t442*t7*t29;
      t655 = t14*t181;
      t660 = 2.0*t442*t17*t91-4.0*t442*t91*t29+t43*t37+4.0*t321*t457+2.0*t250*
t257*t15+2.0*t39*t529-2.0*t442*t127*t83+2.0*t442*t127*t55+2.0*t442*t127*t27+2.0
*t337*t181*t24+4.0*t321*t71-2.0*t655*t457-2.0*t655*t71;
      t663 = t4*t6;
      t668 = t217*t26;
      t674 = t54*t131;
      t677 = t6*t83;
      t684 = t6*t55;
      t695 = t6*t143;
      t702 = -8.0*t227*t44*t663-2.0*t173*t464-2.0*t173*t668-2.0*t192*t14*t37
-4.0*t674*t203-4.0*t362*t677*t29+2.0*t203*t204*t96-4.0*t362*t684*t29-2.0*t203*
t300*t91-t158*t149*t118+t362*t380*t36+t212*t695*t127+4.0*t53*t196*t14*t1;
      t709 = t177*t2;
      t742 = -2.0*t460*t7*t9+4.0*t460*t91*t9+2.0*t709*t127*t17+2.0*t235*t17*t27

+2.0*t235*t17*t83-2.0*t235*t17*t55-2.0*t235*t20*t311+2.0*t246*t66*t29-4.0*t246*
t36*t27-8.0*t39*t209*t4+4.0*t452*t304+4.0*t246*t209*t1-2.0*t35*t143*t29;
      t747 = t24*t217;
      t758 = t24*t170;
      t766 = t412*t9;
      t769 = t1*t311;
      t777 = t4*t170;
      t789 = -2.0*t747*t384-2.0*t460*t91*t127-4.0*t60*t402*t83-t250*t325*t127+
4.0*t3*t758*t9+2.0*t158*t4*t118*t9+2.0*t362*t766+2.0*t250*t769*t127+4.0*t457*t2
*t83*t55+2.0*t158*t777*t29+t65*t67*t17-2.0*t54*t21*t27+2.0*t250*t366*t9;
      t803 = t49*t26;
      t808 = t1*t83*t9;
      t811 = t695*t20;
      t813 = t118*t29;
      t827 = t24*t49;
      t830 = 2.0*t212*t419*t17-8.0*t299*t677*t55-2.0*t3*t758*t127-2.0*t25*t26*
t127*t36-8.0*t3*t803*t24-4.0*t299*t808-t399*t811-4.0*t88*t813+2.0*t88*t449+2.0*
t709*t20*t27+4.0*t246*t311*t29-2.0*t246*t311*t17+4.0*t827*t384;
      t869 = 2.0*t250*t181*t15-4.0*t39*t311*t9+2.0*t39*t311*t127+4.0*t227*t181*
t663-t399*t504-2.0*t203*t8*t29-4.0*t71*t13*t27*t29+4.0*t362*t574*t29+3.0*t60*
t159*t36+4.0*t71*t684*t9+4.0*t71*t677*t9+4.0*t374*t304*t1+2.0*t71*t205;
      t876 = t24*t223;
      t908 = -2.0*t452*t96*t127+4.0*t452*t96*t9-2.0*t876*t384+2.0*t452*t170*t29
+4.0*t441*t223*t285-8.0*t441*t49*t285-4.0*t452*t538+2.0*t77*t217*t14+2.0*t158*
t223*t15-4.0*t158*t49*t15+2.0*t158*t217*t15-2.0*t3*t14*t56-2.0*t3*t14*t40;
      t914 = t14*t17;
      t926 = t412*t17;
      t931 = t54*t196;
      t944 = t159*t27;
      t947 = 2.0*t25*t405*t127+t3*t914*t118-4.0*t158*t162*t27-2.0*t212*t695*t9+
4.0*t95*t8*t27+t399*t926-4.0*t192*t257*t286+4.0*t931*t203+4.0*t71*t13*t83*t29+
4.0*t551*t389-2.0*t551*t508+2.0*t362*t571+2.0*t95*t944;
      t951 = t2*t143;
      t954 = t26*t20;
      t969 = t204*t311;
      t977 = t402*t17;
      t984 = 2.0*t212*t492*t17+t551*t951*t127-2.0*t77*t954*t91+4.0*t3*t14*t96*
t29+t48*t174+2.0*t60*t109*t118+8.0*t192*t44*t286+2.0*t551*t969-t212*t366*t20-
t25*t1*t170*t20-2.0*t362*t977+12.0*t71*t416+t25*t279*t127;
      t995 = t2*t170;
      t1003 = t4*t27;
      t1014 = t257*t14;
      t1021 = -t457*t109*t7-4.0*t399*t339-4.0*t25*t113*t9+2.0*t65*t977+t125*
t995*t17-2.0*t125*t995*t29+t77*t954*t7+4.0*t60*t1003*t29-2.0*t203*t380*t55-t246
*t143*t127+2.0*t362*t110+4.0*t192*t1014*t1-2.0*t192*t589*t127;
      t1030 = t2*t118;
      t1054 = -2.0*t25*t944-2.0*t551*t766+12.0*t71*t674-12.0*t71*t931-2.0*t125*
t1030*t9+2.0*t362*t300*t55-2.0*t876*t71-4.0*t39*t96*t29-2.0*t39*t813+2.0*t39*
t17*t96+4.0*t827*t71-2.0*t747*t71-2.0*t442*t20*t96;
      t1097 = -4.0*t337*t44*t24+4.0*t441*t217*t285+4.0*t25*t1003*t9-4.0*t53*
t131*t14*t1+4.0*t3*t14*t83*t55-t77*t78*t127+4.0*t203*t434*t29-2.0*t203*t380*t27
+4.0*t125*t2*t7*t27+2.0*t457*t109*t91+t95*t380*t118+4.0*t60*t808-4.0*t3*t217*
t286;
      t1136 = -2.0*t551*t951*t9-2.0*t362*t542-4.0*t53*t209*t14*t1+4.0*t374*t132
*t1-2.0*t95*t308*t29-2.0*t212*t478*t127-2.0*t192*t193*t9+t77*t13*t170*t20-2.0*
t384*t385*t127-2.0*t384*t67*t9-4.0*t25*t485*t9+4.0*t299*t62-4.0*t192*t181*t286;
      t1146 = t6*t118;
      t1150 = t14*t143;
      t1177 = 8.0*t3*t49*t286+t384*t67*t127+2.0*t54*t21*t55-2.0*t95*t1146*t9
-2.0*t192*t1150*t17+2.0*t551*t385*t17-2.0*t457*t1030*t29-2.0*t3*t14*t118*t29
-2.0*t3*t914*t96+4.0*t3*t668*t24+t457*t468*t118-3.0*t25*t119+8.0*t299*t1146*t29
;
      t1180 = t300*t7;
      t1211 = t5*t8*t127+t203*t1180+4.0*t192*t1150*t29+t35*t143*t17-t228*t36*
t17+t228*t143*t20-t235*t127*t36-4.0*t374*t26*t36*t27+6.0*t60*t166*t55-t192*t15*
t18+4.0*t35*t44*t1-2.0*t35*t257*t1-2.0*t35*t181*t1;
      t1235 = t13*t91;
      t1251 = -2.0*t1014*t457+4.0*t227*t257*t663+4.0*t228*t236*t29-2.0*t228*
t236*t17-2.0*t16*t20*t236+2.0*t337*t257*t24-t374*t342*t17-2.0*t71*t1180-2.0*t77
*t1235*t127-2.0*t53*t265*t29-3.0*t60*t4*t66*t20-4.0*t95*t276-4.0*t551*t385*t29;
      t1280 = 2.0*t95*t97*t127-4.0*t299*t395*t127+3.0*t25*t369-t362*t395*t20+
2.0*t399*t357-t442*t17*t7-t452*t170*t17+t77*t26*t156-t3*t15*t156+t173*t170*t127
+t460*t170*t20+t16*t20*t36-t246*t66*t17;
      t1315 = t155*t20*t7+4.0*t71*t13*t55*t29-4.0*t250*t769*t9-t5*t308*t20+2.0*
t158*t26*t118*t29-t337*t926+t337*t811+4.0*t173*t803-2.0*t457*t2*t40-2.0*t457*t2
*t84-2.0*t457*t2*t56-2.0*t155*t20*t91-4.0*t88*t524;
      t1351 = -2.0*t709*t20*t83-2.0*t709*t20*t55+4.0*t88*t83*t55+2.0*t158*t26*
t40+2.0*t158*t26*t84+2.0*t158*t26*t56-2.0*t1014*t71+2.0*t452*t118*t9+2.0*t399*
t413-t158*t777*t17-2.0*t362*t969+4.0*t77*t1235*t9+4.0*t192*t655*t1;
      t1359 = pow(t203-t457-t71+2.0*t286-t384,2.0);
      t1363 = -(t141+t379+t506+t232+t1351+t1251+t1136+t545+t1211+t185+t660+t908
+t1097+t1054+t702+t984+t1177+t426+t82+t830+t1021+t947+t283+t626+t742+t1280+
t1315+t586+t332+t789+t869+t463)/t1359/4.0;

 Ricci44[i][j][k] = -t1363;

}}}

return;

}


