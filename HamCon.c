#include "DeclareFunctions.h"

//Find Hamiltonian Constraint

// Call ginv before calling HamCon

double HamCon(Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,Scalar hc,Params Par,double KupKij[XMAX][YMAX][ZMAX], Connection C2,Tensor Ricci){
 int i,j,k,m;
 double hc_rms,***R,***K;
 Scalar detg;
 double vol;
 FILE *fp4;

detg.s = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
R  = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
K  = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
/*
Ricci.xx     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Ricci.yy     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Ricci.zz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Ricci.xy     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Ricci.xz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Ricci.yz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gup.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gup.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gup.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gup.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gup.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gup.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
*/
/*Kup.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kup.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kup.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kup.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kup.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Kup.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);*/
/*C2xxx= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2xxy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2xxz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2yxy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2xyy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2xyz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2yxz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2zyy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2xzz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2zxy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2yyy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2yyz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2yxx= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2zxx= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2yzz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2zxz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2zyz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C2zzz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
*/

//ginv(gij,&gup,Par);

GetRicci(gij,gup,&Ricci,Par,C2);

//*kc_rms = MomCon(gij,gup,Kij,Par,C2,mc,DxK,DyK,DzK,Ricci);
//printf("in HamCon, kc_rms = %2.14f \n",*kc_rms);

for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){

K[i][j][k] = gup.xx[i][j][k]*Kij.xx[i][j][k] + gup.yy[i][j][k]*Kij.yy[i][j][k] +gup.zz[i][j][k]*Kij.zz[i][j][k] +
       2.*gup.xy[i][j][k]*Kij.xy[i][j][k] +2.*gup.xz[i][j][k]*Kij.xz[i][j][k] +2.*gup.yz[i][j][k]*Kij.yz[i][j][k] ;
R[i][j][k] = gup.xx[i][j][k]*Ricci.xx[i][j][k] + gup.yy[i][j][k]*Ricci.yy[i][j][k] +gup.zz[i][j][k]*Ricci.zz[i][j][k] +
       2.*gup.xy[i][j][k]*Ricci.xy[i][j][k] +2.*gup.xz[i][j][k]*Ricci.xz[i][j][k] +2.*gup.yz[i][j][k]*Ricci.yz[i][j][k] ;


Kup.xx[i][j][k] = 
   gup.xx[i][j][k]*gup.xx[i][j][k]*Kij.xx[i][j][k] + 2.*gup.xx[i][j][k]*gup.xy[i][j][k]*Kij.xy[i][j][k] + 2.*gup.xx[i][j][k]*gup.xz[i][j][k]*Kij.xz[i][j][k]
 + gup.xy[i][j][k]*gup.xy[i][j][k]*Kij.yy[i][j][k] + 2.*gup.xy[i][j][k]*gup.xz[i][j][k]*Kij.yz[i][j][k] + gup.xz[i][j][k]*gup.xz[i][j][k]*Kij.zz[i][j][k];

Kup.yy[i][j][k] = 
   gup.xy[i][j][k]*gup.xy[i][j][k]*Kij.xx[i][j][k] + 2.*gup.xy[i][j][k]*gup.yy[i][j][k]*Kij.xy[i][j][k] + 2.*gup.xy[i][j][k]*gup.yz[i][j][k]*Kij.xz[i][j][k]
 + gup.yy[i][j][k]*gup.yy[i][j][k]*Kij.yy[i][j][k] + 2.*gup.yy[i][j][k]*gup.yz[i][j][k]*Kij.yz[i][j][k] + gup.yz[i][j][k]*gup.yz[i][j][k]*Kij.zz[i][j][k];

Kup.zz[i][j][k] = 
   gup.xz[i][j][k]*gup.xz[i][j][k]*Kij.xx[i][j][k] + 2.*gup.xz[i][j][k]*gup.yz[i][j][k]*Kij.xy[i][j][k] + 2.*gup.xz[i][j][k]*gup.zz[i][j][k]*Kij.xz[i][j][k]
 + gup.yz[i][j][k]*gup.yz[i][j][k]*Kij.yy[i][j][k] + 2.*gup.yz[i][j][k]*gup.zz[i][j][k]*Kij.yz[i][j][k] + gup.zz[i][j][k]*gup.zz[i][j][k]*Kij.zz[i][j][k];

Kup.xy[i][j][k] = 
    gup.xx[i][j][k]*gup.xy[i][j][k]*Kij.xx[i][j][k] + gup.xx[i][j][k]*gup.yy[i][j][k]*Kij.xy[i][j][k] + gup.xx[i][j][k]*gup.yz[i][j][k]*Kij.xz[i][j][k]
 +  gup.xy[i][j][k]*gup.xy[i][j][k]*Kij.xy[i][j][k] + gup.xy[i][j][k]*gup.yy[i][j][k]*Kij.yy[i][j][k] + gup.xy[i][j][k]*gup.yz[i][j][k]*Kij.yz[i][j][k]
 +  gup.xz[i][j][k]*gup.xy[i][j][k]*Kij.xz[i][j][k] + gup.xz[i][j][k]*gup.yy[i][j][k]*Kij.yz[i][j][k] + gup.xz[i][j][k]*gup.yz[i][j][k]*Kij.zz[i][j][k];

Kup.xz[i][j][k] = 
    gup.xx[i][j][k]*gup.xz[i][j][k]*Kij.xx[i][j][k] + gup.xx[i][j][k]*gup.yz[i][j][k]*Kij.xy[i][j][k] + gup.xx[i][j][k]*gup.zz[i][j][k]*Kij.xz[i][j][k]
  + gup.xy[i][j][k]*gup.xz[i][j][k]*Kij.xy[i][j][k] + gup.xy[i][j][k]*gup.yz[i][j][k]*Kij.yy[i][j][k] + gup.xy[i][j][k]*gup.zz[i][j][k]*Kij.yz[i][j][k]
  + gup.xz[i][j][k]*gup.xz[i][j][k]*Kij.xz[i][j][k] + gup.xz[i][j][k]*gup.yz[i][j][k]*Kij.yz[i][j][k] + gup.xz[i][j][k]*gup.zz[i][j][k]*Kij.zz[i][j][k];

Kup.yz[i][j][k] = 
    gup.xy[i][j][k]*gup.xz[i][j][k]*Kij.xx[i][j][k] + gup.xy[i][j][k]*gup.yz[i][j][k]*Kij.xy[i][j][k] + gup.xy[i][j][k]*gup.zz[i][j][k]*Kij.xz[i][j][k]
  + gup.yy[i][j][k]*gup.xz[i][j][k]*Kij.xy[i][j][k] + gup.yy[i][j][k]*gup.yz[i][j][k]*Kij.yy[i][j][k] + gup.yy[i][j][k]*gup.zz[i][j][k]*Kij.yz[i][j][k]
  + gup.yz[i][j][k]*gup.xz[i][j][k]*Kij.xz[i][j][k] + gup.yz[i][j][k]*gup.yz[i][j][k]*Kij.yz[i][j][k] + gup.yz[i][j][k]*gup.zz[i][j][k]*Kij.zz[i][j][k];

KupKij[i][j][k] = Kup.xx[i][j][k]*Kij.xx[i][j][k] + 2.*Kup.xy[i][j][k]*Kij.xy[i][j][k] + 2.*Kup.xz[i][j][k]*Kij.xz[i][j][k]
		+ Kup.yy[i][j][k]*Kij.yy[i][j][k] + 2.*Kup.yz[i][j][k]*Kij.yz[i][j][k] + Kup.zz[i][j][k]*Kij.zz[i][j][k];

hc.s[i][j][k] = R[i][j][k] - KupKij[i][j][k] + K[i][j][k]*K[i][j][k];

}}}

GetDetg(gij,detg,Par);
hc_rms = 0.0;
vol = 0.0;

//Get hc_rootmeansquare 
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
   hc_rms += sqrt(detg.s[i][j][k])*hc.s[i][j][k]*hc.s[i][j][k]; 
   vol += sqrt(detg.s[i][j][k]);
}}}

//hc_rms = sqrt(hc_rms/(Par.nxb*Par.nyb*Par.nzb));
hc_rms = sqrt(hc_rms/vol);

/*free_dVector3D(Ricci.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricci.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricci.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricci.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricci.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricci.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
*/
/*free_dVector3D(Kup.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kup.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kup.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kup.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kup.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kup.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);*/
/*free_dVector3D(C2xxx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2xxy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2xxz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2yxy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2xyy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2xyz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2yxz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2zyy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2xzz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2zxy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2yyy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2yyz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2yxx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2zxx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2yzz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2zxz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2zyz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2zzz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
*/
free_dVector3D(R,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(K,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(detg.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);

return hc_rms;

}
