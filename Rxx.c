#include "DeclareFunctions.h"

void Rxx(double ***Ricci22,Tensor gij,Tensor d2gxx,Tensor d2gyy,Tensor d2gzz,Tensor d2gxy,Tensor d2gxz,Tensor d2gyz,Vector dgxx,Vector dgyy, Vector dgzz, Vector dgxy, Vector dgxz, Vector dgyz,Params Par){

 int i,j,k;
 double t1,t2,t3,t4,t5,t6,t7,t10,t13,t14,t15,t16,t18,t21,t22,t23,t24,t25,t26,t29,t30,t31,t34,t35,t36,t38,t42;
 double t45,t47,t50,t54,t57,t62,t63,t65,t66,t67,t68,t72,t76,t80,t81,t82,t86,t87,t88,t89,t91,t94,t95,t96,t97;
 double t98,t101,t103,t106,t107,t109,t110,t111,t117,t118,t122,t128,t129,t130,t133,t134,t135,t138,t139,t140;
 double t144,t148,t149,t150,t153,t161,t162,t164,t165,t170,t173,t174,t175,t178,t180,t185,t191,t193,t198,t205,t208;
 double t209,t210,t213,t216,t225,t228,t230,t231,t234,t243,t244,t245,t255,t267,t286,t297,t306,t310,t316,t317,t323;
 double t330,t331,t350,t355,t358,t366,t370,t385,t388,t398,t404,t409,t416,t448,t453,t465,t469,t478,t481,t482,t484;
 double t487,t489,t490,t491,t493,t496,t502,t503,t505,t509,t512,t515,t521,t528,t529,t533,t536,t543,t548,t550,t551;
 double t553,t554,t555,t557,t561,t562,t567,t573,t574,t580,t584,t588,t598,t601,t602,t605,t608,t611,t614,t617,t623;
 double t624,t632,t635,t638,t639,t642,t645,t650,t656,t659,t662,t665,t673,t680,t683,t684,t687,t690,t693,t696,t706;
 double t713,t717,t747,t756,t769,t788,t791,t794,t830,t856,t857,t860,t863,t866,t869,t870,t873,t876,t907,t916,t953;
 double t973,t993,t1014,t1029,t1033,t1037,t1038,t1041,t1044,t1049,t1068,t1080,t1090,t1112,t1134,t1137,t1140,t1146;
 double t1149,t1152,t1175,t1188,t1223,t1263,t1296,t1304,t1308;


for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){

//> C(Rxx,optimized);

      t1 = gij.xx[i][j][k];
      t2 = t1*t1;


      t3 = gij.yz[i][j][k];
      t4 = t3*t3;
      t5 = t4*t3;
      t6 = t2*t5;


      t7 = d2gxy.xz[i][j][k];
      t10 = d2gyz.xx[i][j][k];


      t13 = gij.xz[i][j][k];
      t14 = t13*t13;
      t15 = t14*t14;
      t16 = d2gxx.yy[i][j][k];

      t18 = gij.yy[i][j][k];

      t21 = gij.xy[i][j][k];
      t22 = t21*t21;
      t23 = t22*t22;


      t24 = dgzz.z[i][j][k];
      t25 = t23*t24;
      t26 = dgxz.x[i][j][k];
      t29 = t2*t4;
      t30 = dgyz.x[i][j][k];
      t31 = t30*t30;
      t34 = dgyy.y[i][j][k];
      t35 = t15*t34;
      t36 = dgxx.y[i][j][k];
      t38 = d2gyy.xx[i][j][k];
      t42 = t14*t22;
      t45 = d2gxx.zz[i][j][k];
      t47 = gij.zz[i][j][k];
      t50 = d2gxy.xy[i][j][k];
      t54 = d2gxz.xy[i][j][k];
      t57 = t1*t3;
      t62 = t14*t4;
      t63 = t36*t36;
      t65 = 4.0*t6*t7-4.0*t6*t10-2.0*t15*t16*t18-2.0*t25*t26+2.0*t29*t31+t35*

t36-2.0*t15*t38*t18+4.0*t42*t31-2.0*t23*t45*t47+4.0*t15*t50*t18+4.0*t6*t54+4.0*
t57*t54*t14*t18+t62*t63;
      t66 = t13*t21;
      t67 = dgxy.z[i][j][k];
      t68 = t67*t67;
      t72 = d2gxz.xz[i][j][k];
      t76 = d2gzz.xx[i][j][k];
      t80 = dgxx.x[i][j][k];
      t81 = t4*t80;
      t82 = dgxz.y[i][j][k];
      t86 = t18*t18;
      t87 = t2*t86;
      t88 = dgzz.x[i][j][k];
      t89 = t88*t88;
      t91 = t82*t82;
      t94 = t22*t1;
      t95 = dgzz.y[i][j][k];
      t96 = t47*t95;
      t97 = dgxy.x[i][j][k];
      t98 = t96*t97;
      t101 = dgxx.z[i][j][k];
      t103 = d2gxx.yz[i][j][k];
      t106 = t22*t4;
      t107 = t101*t101;
      t109 = t1*t47;
      t110 = t109*t13;
      t111 = t18*t36;
      t117 = t1*t4;
      t118 = t21*t80;
      t122 = 4.0*t57*t66*t68+4.0*t23*t72*t47-2.0*t23*t76*t47+2.0*t66*t81*t82+
t87*t89-2.0*t29*t91-2.0*t94*t98+t25*t101-4.0*t6*t103+t106*t107-2.0*t110*t111*
t82-2.0*t29*t68+2.0*t117*t118*t95;
      t128 = t2*t3;
      t129 = t18*t95;
      t130 = t129*t26;
      t133 = t14*t21;
      t134 = t3*t36;
      t135 = t134*t30;
      t138 = t1*t86;
      t139 = t13*t26;
      t140 = dgxz.z[i][j][k];
      t144 = t13*t97;
      t148 = t22*t13;
      t149 = t3*t101;
      t150 = t149*t30;
      t153 = t149*t67;
      t161 = dgyy.x[i][j][k];
      t162 = t161*t161;
      t164 = t47*t47;
      t165 = t1*t164;
      t170 = t2*t164;
      t173 = -4.0*t57*t103*t22*t47+2.0*t128*t130-2.0*t133*t135+4.0*t138*t139*
t140-4.0*t117*t144*t30-2.0*t148*t150-2.0*t148*t153-4.0*t57*t10*t22*t47+t23*t89+
t15*t162+t165*t18*t63-t109*t4*t63+t170*t34*t36;
      t174 = t22*t21;

      t175 = t174*t13;
      t178 = t47*t107;
      t180 = t1*t18;
      t185 = t95*t36;
      t191 = t22*t18;
      t193 = t174*t3;
      t198 = t174*t47;
      t205 = -t175*t24*t36+t138*t178-t180*t4*t107+t62*t80*t161+t42*t185-t175*
t95*t101+t87*t24*t101-t191*t178-t193*t101*t88+t193*t80*t24+t198*t80*t95+t106*
t80*t88-t198*t88*t36;
      t208 = t14*t13;
      t209 = t208*t18;
      t210 = dgyy.z[i][j][k];
      t213 = t208*t21;
      t216 = t208*t3;
      t225 = t14*t18;
      t228 = t101*t210;
      t230 = t2*t47;
      t231 = t18*t31;
      t234 = t18*t68;
      t243 = t209*t80*t210-t213*t36*t210+t216*t80*t34-t216*t36*t161-t209*t101*
t161-t213*t34*t101-t225*t47*t63+t42*t228+2.0*t230*t231+2.0*t230*t234-2.0*t109*
t22*t31-2.0*t35*t97+t170*t162;
      t244 = t180*t21;
      t245 = t47*t101;
      t255 = t18*t91;
      t267 = dgyz.z[i][j][k];
      t286 = -2.0*t244*t245*t67-2.0*t109*t22*t91+2.0*t230*t16*t4+2.0*t230*t255+
4.0*t193*t26*t140-2.0*t193*t26*t88-2.0*t193*t101*t140+4.0*t175*t26*t267-2.0*
t175*t101*t267+2.0*t175*t24*t97+2.0*t165*t16*t22-2.0*t170*t34*t97+4.0*t170*t50*
t18;
      t297 = t14*t1;
      t306 = dgyz.y[i][j][k];
      t310 = t14*t16;
      t316 = t14*t38;
      t317 = t22*t47;
      t323 = dgxy.y[i][j][k];
      t330 = -2.0*t170*t16*t18-2.0*t170*t38*t18-2.0*t109*t22*t68-2.0*t297*t234+
2.0*t213*t34*t26+2.0*t213*t97*t210+4.0*t213*t97*t306-2.0*t310*t117+2.0*t42*t101
*t306-2.0*t316*t317-4.0*t42*t26*t306+2.0*t209*t101*t323+4.0*t175*t103*t47;
      t331 = t2*t18;
      t350 = t10*t18;
      t355 = t26*t210;
      t358 = t14*t50;
      t366 = t1*t5;
      t370 = -4.0*t331*t72*t4-2.0*t180*t22*t89+2.0*t331*t45*t4+4.0*t87*t72*t47
-2.0*t87*t45*t47-2.0*t87*t76*t47+4.0*t213*t350-2.0*t316*t117-2.0*t42*t355+4.0*
t358*t317+4.0*t358*t117-4.0*t213*t161*t30-2.0*t366*t80*t82;
      t385 = t95*t97;
      t388 = t88*t161;
      t398 = t54*t18;
      t404 = t103*t18;
      t409 = -2.0*t366*t80*t67+2.0*t29*t185+2.0*t366*t80*t30+4.0*t29*t82*t67
-4.0*t29*t385+2.0*t29*t388-4.0*t29*t355+2.0*t29*t228+2.0*t366*t101*t36-4.0*t213
*t398+2.0*t175*t95*t26+4.0*t213*t404+2.0*t42*t388;
      t416 = t7*t18;
      t448 = 4.0*t175*t10*t47-4.0*t175*t88*t30-4.0*t213*t416+8.0*t42*t54*t3+8.0
*t42*t7*t3-8.0*t42*t103*t3-8.0*t42*t10*t3-4.0*t175*t54*t47-2.0*t42*t385-4.0*
t175*t7*t47-4.0*t138*t72*t14+2.0*t331*t76*t4+2.0*t138*t76*t14;
      t453 = t22*t72;
      t465 = t22*t45;
      t469 = t13*t3;
      t478 = t18*t101;
      t481 = t47*t80;
      t482 = t481*t88;
      t484 = t18*t80;
      t487 = t96*t36;
      t489 = -2.0*t87*t24*t26+4.0*t453*t117-4.0*t42*t97*t267+2.0*t42*t36*t267
-2.0*t198*t80*t267-2.0*t465*t117-8.0*t174*t72*t469+4.0*t174*t45*t469+4.0*t174*

t76*t469+t148*t478*t88-t191*t482-t148*t484*t24+t94*t487;
      t490 = t3*t24;
      t491 = t490*t36;
      t493 = t47*t36;
      t496 = t134*t210;
      t502 = t3*t34;
      t503 = t502*t101;
      t505 = t478*t210;
      t509 = t22*t95;
      t512 = t21*t36;
      t515 = t3*t80;
      t521 = t94*t491+t133*t493*t161+t297*t496-t225*t481*t161-t133*t481*t34-
t230*t503-t230*t505-t109*t81*t161+t57*t509*t101-t165*t512*t161-t110*t515*t34-
t230*t496+t165*t484*t161;
      t528 = t109*t21;
      t529 = t515*t210;
      t533 = t21*t34;
      t536 = t22*t101;
      t543 = t149*t161;
      t548 = t129*t101;
      t550 = t165*t118*t34+t110*t134*t161-t528*t529+t110*t478*t161+t110*t533*
t101+t109*t536*t210-t110*t484*t210+t110*t512*t210+t528*t543+t297*t503+t297*t505
-t331*t487-t128*t548;
      t551 = t66*t1;
      t553 = t57*t13;
      t554 = t18*t88;
      t555 = t554*t36;
      t557 = t484*t95;
      t561 = t47*t97;
      t562 = t561*t30;
      t567 = t14*t95;
      t573 = t180*t13;
      t574 = t21*t24;
      t580 = t13*t101;
      t584 = t551*t548+t553*t555-t553*t557+t244*t149*t88-4.0*t148*t562-t244*
t515*t24+t180*t567*t36-t244*t481*t95-t331*t491+t573*t574*t36+t244*t47*t88*t36-
t138*t580*t88+t138*t482;
      t588 = t13*t80;
      t598 = t481*t140;
      t601 = t3*t97;
      t602 = t601*t140;
      t605 = t601*t88;
      t608 = t134*t140;
      t611 = t561*t267;
      t614 = t493*t267;
      t617 = t490*t97;
      t623 = -t180*t81*t88+t138*t588*t24-3.0*t133*t529+2.0*t133*t557+4.0*t244*
t561*t140+2.0*t191*t598+4.0*t148*t602-2.0*t148*t605-2.0*t148*t608+4.0*t94*t611
-2.0*t94*t614-2.0*t94*t617+4.0*t109*t310*t18;
      t624 = t21*t26;
      t632 = t3*t161*t30;
      t635 = t601*t210;
      t638 = t3*t26;
      t639 = t638*t323;
      t642 = t638*t161;
      t645 = t149*t323;
      t650 = t601*t306;
      t656 = t134*t82;
      t659 = t134*t67;
      t662 = t515*t306;
      t665 = 4.0*t117*t624*t67+2.0*t133*t493*t323+4.0*t297*t632-2.0*t297*t635+
4.0*t133*t639-2.0*t133*t642-2.0*t133*t645+3.0*t133*t543-4.0*t297*t650+2.0*t225*
t481*t323-2.0*t133*t656+6.0*t133*t659+2.0*t133*t662;
      t673 = t502*t26;
      t680 = t18*t82*t67;
      t683 = t18*t26;
      t684 = t683*t306;
      t687 = t683*t210;
      t690 = t478*t306;
      t693 = t134*t306;
      t696 = t21*t97;
      t706 = t22*t26;
      t713 = -2.0*t110*t478*t323+2.0*t230*t673-2.0*t110*t533*t26-4.0*t230*t680
-4.0*t230*t684+2.0*t230*t687+2.0*t230*t690-2.0*t230*t693-2.0*t110*t696*t210-4.0
*t110*t696*t306-2.0*t109*t536*t306+4.0*t109*t706*t306-8.0*t109*t358*t18;
      t717 = t14*t34;
      t747 = t66*t3;
      t756 = 4.0*t109*t316*t18+4.0*t109*t717*t97-2.0*t109*t706*t210+4.0*t109*
t22*t82*t67-2.0*t109*t717*t36+4.0*t110*t21*t161*t30-4.0*t133*t561*t323+2.0*t133
*t561*t161-2.0*t110*t683*t161+4.0*t230*t650+8.0*t109*t50*t747-4.0*t109*t16*t747
-4.0*t109*t38*t747;
      t769 = t22*t24;
      t788 = t638*t267;
      t791 = t149*t267;
      t794 = -2.0*t165*t484*t323+2.0*t528*t656+2.0*t528*t135+2.0*t528*t659-2.0*
t528*t662+4.0*t180*t769*t26-2.0*t110*t111*t30+2.0*t110*t111*t67-4.0*t148*t683*
t140+2.0*t148*t683*t88+2.0*t148*t478*t140-4.0*t94*t788+2.0*t94*t791;
      t830 = 2.0*t109*t81*t323+2.0*t110*t484*t306+4.0*t165*t696*t323-2.0*t165*
t696*t161-2.0*t165*t512*t323-4.0*t110*t601*t323+2.0*t110*t601*t161+2.0*t110*
t134*t323-4.0*t230*t632+2.0*t230*t635+2.0*t110*t512*t306-4.0*t528*t639+2.0*t528
*t642;
      t856 = t66*t54;
      t857 = t180*t47;
      t860 = t66*t7;
      t863 = t66*t103;
      t866 = t66*t10;
      t869 = t66*t18;
      t870 = t481*t67;
      t873 = t481*t82;
      t876 = 2.0*t528*t645+4.0*t110*t683*t323+4.0*t57*t22*t88*t30-4.0*t57*t10*
t14*t18+4.0*t57*t7*t14*t18+2.0*t117*t588*t210-4.0*t128*t398*t47+4.0*t856*t857+
4.0*t860*t857-4.0*t863*t857-4.0*t866*t857-2.0*t869*t870-2.0*t869*t873;
      t907 = t22*t76;
      t916 = -2.0*t148*t245*t161-4.0*t133*t601*t67-4.0*t133*t601*t82+4.0*t133*
t601*t30-4.0*t148*t638*t67-2.0*t66*t81*t30+2.0*t198*t101*t67+2.0*t198*t101*t30
-2.0*t198*t101*t82+2.0*t138*t45*t14-2.0*t907*t117-2.0*t106*t80*t140+2.0*t198*
t36*t140;
      t953 = -4.0*t198*t97*t140-2.0*t907*t225+4.0*t453*t225-2.0*t465*t225+2.0*

t198*t97*t88+2.0*t165*t38*t22+2.0*t230*t38*t4-4.0*t165*t50*t22-4.0*t230*t50*t4
-2.0*t109*t14*t162-2.0*t310*t317-2.0*t62*t80*t323-2.0*t209*t80*t306;
      t973 = t21*t3;
      t993 = 4.0*t216*t97*t323-2.0*t216*t97*t161-2.0*t216*t36*t323-2.0*t213*t36
*t306-4.0*t209*t26*t323+2.0*t209*t26*t161-8.0*t208*t50*t973+4.0*t208*t16*t973+
4.0*t208*t38*t973+2.0*t209*t36*t82+2.0*t209*t36*t30-2.0*t209*t36*t67-2.0*t297*
t255;
      t1014 = t57*t18;
      t1029 = t21*t95;
      t1033 = -2.0*t297*t231-4.0*t128*t416*t47-12.0*t117*t860+4.0*t128*t404*t47
+12.0*t117*t863+4.0*t128*t350*t47+12.0*t117*t866+2.0*t1014*t870+2.0*t1014*t873+
4.0*t117*t144*t67+4.0*t117*t144*t82+4.0*t57*t54*t22*t47+8.0*t553*t1029*t97;
      t1037 = t57*t21;
      t1038 = t561*t82;
      t1041 = t478*t82;
      t1044 = t561*t67;
      t1049 = t21*t101;
      t1068 = 4.0*t117*t624*t82-4.0*t1037*t1038+2.0*t553*t1041-4.0*t1037*t1044+
4.0*t1037*t562-4.0*t117*t1049*t82+4.0*t57*t66*t91-2.0*t297*t673+4.0*t297*t680+
4.0*t297*t684-2.0*t297*t687-2.0*t297*t690+2.0*t297*t693;
      t1080 = t21*t88;
      t1090 = t554*t30;
      t1112 = -2.0*t57*t509*t26-4.0*t57*t103*t14*t18-4.0*t553*t1029*t36-2.0*
t117*t1080*t36-4.0*t553*t1080*t161-4.0*t553*t1049*t210-4.0*t128*t1090-4.0*t57*
t66*t31-4.0*t117*t13*t36*t67-8.0*t553*t21*t82*t67-12.0*t117*t856+3.0*t148*t3*
t88*t36+4.0*t551*t1090;
      t1134 = t683*t67;
      t1137 = t683*t82;
      t1140 = t683*t30;
      t1146 = t245*t36;
      t1149 = 2.0*t66*t81*t67-4.0*t148*t493*t67+4.0*t148*t1038-4.0*t133*t1041
-3.0*t148*t515*t95+4.0*t148*t1044+6.0*t148*t149*t82-4.0*t117*t624*t30-4.0*t553*
t1134-4.0*t553*t1137+4.0*t553*t1140+8.0*t553*t624*t210-2.0*t1014*t1146;
      t1152 = t481*t30;
      t1175 = t515*t267;
      t1188 = -2.0*t1014*t1152-2.0*t117*t580*t161+4.0*t57*t7*t22*t47+2.0*t573*
t150+2.0*t573*t153+2.0*t244*t481*t267-4.0*t573*t624*t267+2.0*t573*t1049*t267+
2.0*t148*t1175+2.0*t331*t614+2.0*t331*t617-2.0*t573*t574*t97-4.0*t244*t638*t140
;
      t1223 = 2.0*t244*t638*t88+2.0*t244*t149*t140-2.0*t138*t139*t88-2.0*t138*
t580*t140-2.0*t244*t561*t88-4.0*t573*t602+2.0*t573*t605+2.0*t573*t608-4.0*t331*
t611+2.0*t331*t98+8.0*t180*t72*t747-4.0*t180*t45*t747-4.0*t180*t76*t747;
      t1263 = -2.0*t138*t598-2.0*t244*t245*t30+2.0*t244*t245*t82-8.0*t180*t453*
t47+4.0*t180*t465*t47+4.0*t331*t788-2.0*t331*t791-2.0*t66*t4*t101*t36+4.0*t180*
t14*t97*t267-2.0*t573*t1175-2.0*t180*t14*t36*t267-2.0*t180*t769*t101-2.0*t180*
t567*t97;
      t1296 = -4.0*t148*t638*t82+4.0*t148*t638*t30-2.0*t133*t555+4.0*t133*t1134
+4.0*t133*t1137-4.0*t133*t1140-2.0*t551*t130+2.0*t869*t1146+2.0*t869*t1152+2.0*
t148*t481*t210+4.0*t180*t907*t47+2.0*t180*t81*t140-2.0*t244*t493*t140;
      t1304 = pow(t857-t117-t317+2.0*t747-t225,2.0);
      t1308 = -(t916+t993+t623+t876+t1149+t1296+t794+t1033+t953+t65+t122+t1188+
t1068+t173+t205+t243+t286+t1223+t330+t370+t665+t409+t448+t1263+t489+t521+t550+
t584+t713+t756+t1112+t830)/t1304/4.0;

   Ricci22[i][j][k] = -t1308;
}}}

return;
}


