#include "DeclareFunctions.h"

void Loperator(Vector Eks,Tensor LEks, Tensor gij, Tensor gup, Connection C2, Params Par){

 int i, j, k;
 Vector dEksx, dEksy, dEksz;
 double trace3[XMAX][YMAX][ZMAX];

dEksx.x     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dEksx.y     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dEksx.z     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dEksy.x     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dEksy.y     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dEksy.z     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dEksz.x     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dEksz.y     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dEksz.z     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

 FunctionFirstDerivs(Eks.x,&dEksx,Par);
 FunctionFirstDerivs(Eks.y,&dEksy,Par);
 FunctionFirstDerivs(Eks.z,&dEksz,Par);

 for(i=0;i<Par.nxb;i++){
   for(j=0;j<Par.nyb;j++){
      for(k=0;k<Par.nzb;k++){
         LEks.xx[i][j][k] = dEksx.x[i][j][k] + dEksx.x[i][j][k] - 2.0*(C2.xxx[i][j][k]*Eks.x[i][j][k]
                                                                      +C2.yxx[i][j][k]*Eks.y[i][j][k]
                                                                      +C2.zxx[i][j][k]*Eks.z[i][j][k]);
         LEks.yy[i][j][k] = dEksy.y[i][j][k] + dEksy.y[i][j][k] - 2.0*(C2.xyy[i][j][k]*Eks.x[i][j][k]
                                                                      +C2.yyy[i][j][k]*Eks.y[i][j][k]
                                                                      +C2.zyy[i][j][k]*Eks.z[i][j][k]);
         LEks.zz[i][j][k] = dEksz.z[i][j][k] + dEksz.z[i][j][k] - 2.0*(C2.xzz[i][j][k]*Eks.x[i][j][k]
                                                                      +C2.yzz[i][j][k]*Eks.y[i][j][k]
                                                                      +C2.zzz[i][j][k]*Eks.z[i][j][k]);
         LEks.xy[i][j][k] = dEksx.y[i][j][k] + dEksy.x[i][j][k] - 2.0*(C2.xxy[i][j][k]*Eks.x[i][j][k]
                                                                      +C2.yxy[i][j][k]*Eks.y[i][j][k]
                                                                      +C2.zxy[i][j][k]*Eks.z[i][j][k]);
         LEks.xz[i][j][k] = dEksx.z[i][j][k] + dEksz.x[i][j][k] - 2.0*(C2.xxz[i][j][k]*Eks.x[i][j][k]
                                                                      +C2.yxz[i][j][k]*Eks.y[i][j][k]
                                                                      +C2.zxz[i][j][k]*Eks.z[i][j][k]);
         LEks.yz[i][j][k] = dEksy.z[i][j][k] + dEksz.y[i][j][k] - 2.0*(C2.xyz[i][j][k]*Eks.x[i][j][k]
                                                                      +C2.yyz[i][j][k]*Eks.y[i][j][k]
                                                                      +C2.zyz[i][j][k]*Eks.z[i][j][k]);
 }}}

 for(i=0;i<Par.nxb;i++){
   for(j=0;j<Par.nyb;j++){
      for(k=0;k<Par.nzb;k++){
        trace3[i][j][k] = (gup.xx[i][j][k]*LEks.xx[i][j][k] + gup.yy[i][j][k]*LEks.yy[i][j][k] 
                         + gup.zz[i][j][k]*LEks.zz[i][j][k] + 2.0*(gup.xy[i][j][k]*LEks.xy[i][j][k]
                         + gup.xz[i][j][k]*LEks.xz[i][j][k] + gup.yz[i][j][k]*LEks.yz[i][j][k]))/3.0;
 }}}

 for(i=0;i<Par.nxb;i++){
   for(j=0;j<Par.nyb;j++){
      for(k=0;k<Par.nzb;k++){
         LEks.xx[i][j][k] = LEks.xx[i][j][k] - trace3[i][j][k]*gij.xx[i][j][k];
         LEks.yy[i][j][k] = LEks.yy[i][j][k] - trace3[i][j][k]*gij.yy[i][j][k];
         LEks.zz[i][j][k] = LEks.zz[i][j][k] - trace3[i][j][k]*gij.zz[i][j][k];
         LEks.xy[i][j][k] = LEks.xy[i][j][k] - trace3[i][j][k]*gij.xy[i][j][k];
         LEks.xz[i][j][k] = LEks.xz[i][j][k] - trace3[i][j][k]*gij.xz[i][j][k];
         LEks.yz[i][j][k] = LEks.yz[i][j][k] - trace3[i][j][k]*gij.yz[i][j][k];
 }}}

free_dVector3D(dEksx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dEksx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dEksx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dEksy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dEksy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dEksy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dEksz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dEksz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dEksz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

 return;
}

