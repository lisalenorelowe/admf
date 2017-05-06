#include "DeclareFunctions.h"

double MomCon(Tensor gij,Tensor gup,Tensor Kij, Tensor Kup,Vector mc, Params Par, double DxK[XMAX][YMAX][ZMAX], double DyK[XMAX][YMAX][ZMAX], double DzK[XMAX][YMAX][ZMAX],Connection C2,Tensor Ricci) {

//Call ginv and HamCon before calling MomCon
 
 int i,j,k,m;
 double ***dummy,Kc_rms,***trK,***Kr_yx,***Kr_zx,***Kr_zy;
 Scalar detg;
 Vector ddetg,dtrK,dKrxx,dKryy,dKrzz,dKrxy,dKrxz,dKryz,dKryx,dKrzx,dKrzy;
 Tensor Kr;
 Vector dKcx, dKcy, dKcz;
 double DaMa, vol;
 FILE *fp3;

dKcx.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKcy.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKcz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKcx.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKcy.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKcz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKcx.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKcy.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKcz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Kc.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
//Kc.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
//Kc.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

Kr.xx=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kr.yy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kr.zz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kr_yx=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kr_zx=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kr_zy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kr.xy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kr.xz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kr.yz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

dKrxx.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrxx.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrxx.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKryy.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKryy.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKryy.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrzz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrzz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrzz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrxy.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrxy.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrxy.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKryz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKryz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKryz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrxz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrxz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrxz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrzx.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrzx.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrzx.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrzy.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrzy.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKrzy.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKryx.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKryx.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKryx.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

trK=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dtrK.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dtrK.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dtrK.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
detg.s=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dummy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
ddetg.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
ddetg.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
ddetg.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

/*----------------------------------------------------------------*/

//Get K = TraceK = g^{ab}K_{ab} 
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
    trK[i][j][k] = gup.xx[i][j][k]*Kij.xx[i][j][k] + gup.yy[i][j][k]*Kij.yy[i][j][k] +gup.zz[i][j][k]*Kij.zz[i][j][k] +
       2.*gup.xy[i][j][k]*Kij.xy[i][j][k] +2.*gup.xz[i][j][k]*Kij.xz[i][j][k] +2.*gup.yz[i][j][k]*Kij.yz[i][j][k] ;
}}}

//Get K^b_a = g^{bc}K_{ac}, (==Kr.ab;  Kr.ab != Kr.ba)
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
   Kr.xx[i][j][k] = gup.xx[i][j][k]*Kij.xx[i][j][k] + gup.xy[i][j][k]*Kij.xy[i][j][k] + gup.xz[i][j][k]*Kij.xz[i][j][k]; 
   Kr.yy[i][j][k] = gup.xy[i][j][k]*Kij.xy[i][j][k] + gup.yy[i][j][k]*Kij.yy[i][j][k] + gup.yz[i][j][k]*Kij.yz[i][j][k];
   Kr.zz[i][j][k] = gup.xz[i][j][k]*Kij.xz[i][j][k] + gup.yz[i][j][k]*Kij.yz[i][j][k] + gup.zz[i][j][k]*Kij.zz[i][j][k];
   Kr.xy[i][j][k] = gup.xy[i][j][k]*Kij.xx[i][j][k] + gup.yy[i][j][k]*Kij.xy[i][j][k] + gup.yz[i][j][k]*Kij.xz[i][j][k];
   Kr.xz[i][j][k] = gup.xz[i][j][k]*Kij.xx[i][j][k] + gup.yz[i][j][k]*Kij.xy[i][j][k] + gup.zz[i][j][k]*Kij.xz[i][j][k];
   Kr.yz[i][j][k] = gup.xz[i][j][k]*Kij.xy[i][j][k] + gup.yz[i][j][k]*Kij.yy[i][j][k] + gup.zz[i][j][k]*Kij.yz[i][j][k];
   Kr_yx[i][j][k] = gup.xx[i][j][k]*Kij.xy[i][j][k] + gup.xy[i][j][k]*Kij.yy[i][j][k] + gup.xz[i][j][k]*Kij.yz[i][j][k];
   Kr_zx[i][j][k] = gup.xx[i][j][k]*Kij.xz[i][j][k] + gup.xy[i][j][k]*Kij.yz[i][j][k] + gup.xz[i][j][k]*Kij.zz[i][j][k];
   Kr_zy[i][j][k] = gup.xy[i][j][k]*Kij.xz[i][j][k] + gup.yy[i][j][k]*Kij.yz[i][j][k] + gup.yz[i][j][k]*Kij.zz[i][j][k];
}}}

//get determinant of metric 
GetDetg(gij,detg,Par);

//Find partial derivs of K^b_a (dKrab) and of TrK, detg
FunctionFirstDerivs(Kr.xx,&dKrxx,Par);
FunctionFirstDerivs(Kr.yy,&dKryy,Par);
FunctionFirstDerivs(Kr.zz,&dKrzz,Par);
FunctionFirstDerivs(Kr.xy,&dKrxy,Par);
FunctionFirstDerivs(Kr.xz,&dKrxz,Par);
FunctionFirstDerivs(Kr.yz,&dKryz,Par);
FunctionFirstDerivs(Kr_yx,&dKryx,Par);
FunctionFirstDerivs(Kr_zx,&dKrzx,Par);
FunctionFirstDerivs(Kr_zy,&dKrzy,Par);
FunctionFirstDerivs(trK,&dtrK,Par);
 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
   for(k=0;k<Par.nzb;k++){
    dummy[i][j][k] = sqrt(detg.s[i][j][k]);
    DxK[i][j][k] = dtrK.x[i][j][k];
    DyK[i][j][k] = dtrK.y[i][j][k];
    DzK[i][j][k] = dtrK.z[i][j][k];
 }}}
FunctionFirstDerivs(dummy,&ddetg,Par);

//$M_a = D_b K^b_a - D_a K$
//$D_b K^b_a = \partial_bK^b_a + \Gamma^b_{db} K^d_a - \Gamma^d_{ab}K^b_d $

for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
// add trK part, $D_a K$
   mc.x[i][j][k] = -dtrK.x[i][j][k];  
   mc.y[i][j][k] = -dtrK.y[i][j][k];  
   mc.z[i][j][k] = -dtrK.z[i][j][k];

//add $\partial_b K^b_a$ part
   mc.x[i][j][k] = mc.x[i][j][k] + dKrxx.x[i][j][k] + dKrxy.y[i][j][k] + dKrxz.z[i][j][k];
   mc.y[i][j][k] = mc.y[i][j][k] + dKryx.x[i][j][k] + dKryy.y[i][j][k] + dKryz.z[i][j][k];
   mc.z[i][j][k] = mc.z[i][j][k] + dKrzx.x[i][j][k] + dKrzy.y[i][j][k] + dKrzz.z[i][j][k];

//add $\Gamma^b_{db} K^d_a$ part, $\Gamma^b_{db} = \partial_d\sqrt(g)$
   mc.x[i][j][k] = mc.x[i][j][k] + (ddetg.x[i][j][k]*Kr.xx[i][j][k] + ddetg.y[i][j][k]*Kr.xy[i][j][k] 
	+ ddetg.z[i][j][k]*Kr.xz[i][j][k])/dummy[i][j][k];
   mc.y[i][j][k] = mc.y[i][j][k] + (ddetg.x[i][j][k]*Kr_yx[i][j][k] + ddetg.y[i][j][k]*Kr.yy[i][j][k] 
	+ ddetg.z[i][j][k]*Kr.yz[i][j][k])/dummy[i][j][k];
   mc.z[i][j][k] = mc.z[i][j][k] + (ddetg.x[i][j][k]*Kr_zx[i][j][k] + ddetg.y[i][j][k]*Kr_zy[i][j][k] 
	+ ddetg.z[i][j][k]*Kr.zz[i][j][k])/dummy[i][j][k];

//add $-\Gamma^d_{ab}K^b_d$ part (== C2dab * Kr.db)
   mc.x[i][j][k] = mc.x[i][j][k] -
    (C2.xxx[i][j][k]*Kr.xx[i][j][k] + C2.xxy[i][j][k]*Kr.xy[i][j][k] + C2.xxz[i][j][k]*Kr.xz[i][j][k] +
     C2.yxx[i][j][k]*Kr_yx[i][j][k] + C2.yxy[i][j][k]*Kr.yy[i][j][k] + C2.yxz[i][j][k]*Kr.yz[i][j][k] +
     C2.zxx[i][j][k]*Kr_zx[i][j][k] + C2.zxy[i][j][k]*Kr_zy[i][j][k] + C2.zxz[i][j][k]*Kr.zz[i][j][k]  );

   mc.y[i][j][k] = mc.y[i][j][k] - 
    (C2.xxy[i][j][k]*Kr.xx[i][j][k] + C2.xyy[i][j][k]*Kr.xy[i][j][k] + C2.xyz[i][j][k]*Kr.xz[i][j][k] +
     C2.yxy[i][j][k]*Kr_yx[i][j][k] + C2.yyy[i][j][k]*Kr.yy[i][j][k] + C2.yyz[i][j][k]*Kr.yz[i][j][k] +
     C2.zxy[i][j][k]*Kr_zx[i][j][k] + C2.zyy[i][j][k]*Kr_zy[i][j][k] + C2.zyz[i][j][k]*Kr.zz[i][j][k]  );

   mc.z[i][j][k] = mc.z[i][j][k] - 
    (C2.xxz[i][j][k]*Kr.xx[i][j][k] + C2.xyz[i][j][k]*Kr.xy[i][j][k] + C2.xzz[i][j][k]*Kr.xz[i][j][k] + 
     C2.yxz[i][j][k]*Kr_yx[i][j][k] + C2.yyz[i][j][k]*Kr.yy[i][j][k] + C2.yzz[i][j][k]*Kr.yz[i][j][k] +
     C2.zxz[i][j][k]*Kr_zx[i][j][k] + C2.zyz[i][j][k]*Kr_zy[i][j][k] + C2.zzz[i][j][k]*Kr.zz[i][j][k]  );

}}}

//printf("MC at origin = %2.10f %2.10f %2.10f \n", mcx[Par.nxb/2][Par.nyb/2][Par.nzb/2], mcy[Par.nxb/2][Par.nyb/2][Par.nzb/2], 
//	mcz[Par.nxb/2][Par.nyb/2][Par.nzb/2]);

Kc_rms = 0.0;
vol = 0.0;

for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
//        Kc_rms  += mcx[i][j][k]*mcx[i][j][k] + mcy[i][j][k]*mcy[i][j][k] +mcz[i][j][k]*mcz[i][j][k];
	Kc_rms += dummy[i][j][k]*(mc.x[i][j][k]*gup.xx[i][j][k]*mc.x[i][j][k] + mc.y[i][j][k]*gup.yy[i][j][k]*mc.y[i][j][k] 
		+ mc.z[i][j][k]*gup.zz[i][j][k]*mc.z[i][j][k] + 2.0*mc.x[i][j][k]*gup.xy[i][j][k]*mc.y[i][j][k]
		+ 2.0*mc.x[i][j][k]*gup.xz[i][j][k]*mc.z[i][j][k] + 2.0*mc.y[i][j][k]*gup.yz[i][j][k]*mc.z[i][j][k]) ;
	vol += dummy[i][j][k];
}}}

//Kc_rms = sqrt(Kc_rms/(Par.nxb*Par.nyb*Par.nzb));
Kc_rms = sqrt(Kc_rms/vol);

/*
FunctionFirstDerivs(mcx,&dKcx,Par);
FunctionFirstDerivs(mcy,&dKcy,Par);
FunctionFirstDerivs(mcz,&dKcz,Par);
DaMa = dKcx.x[Par.nxb/2][Par.nyb/2][Par.nzb/2]*gup.xx[Par.nxb/2][Par.nyb/2][Par.nzb/2] 
	+ dKcy.y[Par.nxb/2][Par.nyb/2][Par.nzb/2]*gup.yy[Par.nxb/2][Par.nyb/2][Par.nzb/2]
	+ dKcz.z[Par.nxb/2][Par.nyb/2][Par.nzb/2]*gup.zz[Par.nxb/2][Par.nyb/2][Par.nzb/2]
	+ dKcx.y[Par.nxb/2][Par.nyb/2][Par.nzb/2]*gup.xy[Par.nxb/2][Par.nyb/2][Par.nzb/2] 
	+ dKcy.x[Par.nxb/2][Par.nyb/2][Par.nzb/2]*gup.xy[Par.nxb/2][Par.nyb/2][Par.nzb/2] 
	+ dKcx.z[Par.nxb/2][Par.nyb/2][Par.nzb/2]*gup.xz[Par.nxb/2][Par.nyb/2][Par.nzb/2] 
	+ dKcz.x[Par.nxb/2][Par.nyb/2][Par.nzb/2]*gup.xz[Par.nxb/2][Par.nyb/2][Par.nzb/2] 
	+ dKcy.z[Par.nxb/2][Par.nyb/2][Par.nzb/2]*gup.yz[Par.nxb/2][Par.nyb/2][Par.nzb/2] 
	+ dKcz.y[Par.nxb/2][Par.nyb/2][Par.nzb/2]*gup.yz[Par.nxb/2][Par.nyb/2][Par.nzb/2];

printf("DaMa at origin = %2.10f\n",DaMa); 
fp3 = fopen("DaMafile","a");
fprintf(fp3,"%2.10f \n",DaMa);
fclose(fp3);
*/

free_dVector3D(dKcx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKcx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKcx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKcy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKcy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKcy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKcz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKcz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKcz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

free_dVector3D(Kr.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kr.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kr.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kr.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kr.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kr.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kr_yx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kr_zy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kr_zx,0,Par.nxb,0,Par.nyb,0,Par.nzb);

free_dVector3D(dKrxx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrxx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrxx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKryy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKryy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKryy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrzz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrzz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrzz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrxy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrxy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrxy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKryz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKryz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKryz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrxz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrxz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrxz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKryx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKryx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKryx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrzy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrzy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrzy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrzx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrzx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKrzx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(trK,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dtrK.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dtrK.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dtrK.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(detg.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dummy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(ddetg.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(ddetg.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(ddetg.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

return Kc_rms;

}
