#include "DeclareFunctions.h"

void Deltaop(Tensor LEks,Vector DelEks, Tensor gij, Tensor gup, Connection C2, Params Par){

 int i, j, k;
 Vector dLEksxx, dLEksxy, dLEksxz, dLEksyy, dLEksyz, dLEkszz;

 dLEksxx.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksxx.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksxx.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksxy.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksxy.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksxy.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksxz.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksxz.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksxz.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksyy.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksyy.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksyy.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksyz.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksyz.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEksyz.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEkszz.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEkszz.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 dLEkszz.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

 FunctionFirstDerivs(LEks.xx,&dLEksxx,Par);
 FunctionFirstDerivs(LEks.xy,&dLEksxy,Par);
 FunctionFirstDerivs(LEks.xz,&dLEksxz,Par);
 FunctionFirstDerivs(LEks.yy,&dLEksyy,Par);
 FunctionFirstDerivs(LEks.yz,&dLEksyz,Par);
 FunctionFirstDerivs(LEks.zz,&dLEkszz,Par);

 for(i=0;i<Par.nxb;i++){
   for(j=0;j<Par.nyb;j++){
      for(k=0;k<Par.nzb;k++){
        DelEks.x[i][j][k] = gup.xx[i][j][k]*(dLEksxx.x[i][j][k] - C2.xxx[i][j][k]*LEks.xx[i][j][k] 
                                         - C2.yxx[i][j][k]*LEks.xy[i][j][k] - C2.zxx[i][j][k]*LEks.xz[i][j][k]
                                         - C2.xxx[i][j][k]*LEks.xx[i][j][k] - C2.yxx[i][j][k]*LEks.xy[i][j][k]
                                         - C2.zxx[i][j][k]*LEks.xz[i][j][k])
                          + gup.yy[i][j][k]*(dLEksxy.y[i][j][k] - C2.xxy[i][j][k]*LEks.xy[i][j][k] 
                                         - C2.yxy[i][j][k]*LEks.yy[i][j][k] - C2.zxy[i][j][k]*LEks.yz[i][j][k]
                                         - C2.xyy[i][j][k]*LEks.xx[i][j][k] - C2.yyy[i][j][k]*LEks.xy[i][j][k]
                                         - C2.zyy[i][j][k]*LEks.xz[i][j][k])
                          + gup.zz[i][j][k]*(dLEksxz.z[i][j][k] - C2.xxz[i][j][k]*LEks.xz[i][j][k] 
                                         - C2.yxz[i][j][k]*LEks.yz[i][j][k] - C2.zxz[i][j][k]*LEks.zz[i][j][k]
                                         - C2.xzz[i][j][k]*LEks.xx[i][j][k] - C2.yzz[i][j][k]*LEks.xy[i][j][k]
                                         - C2.zzz[i][j][k]*LEks.xz[i][j][k])
                          + gup.xy[i][j][k]*(dLEksxy.x[i][j][k] + dLEksxx.y[i][j][k] 
                                         - C2.xxx[i][j][k]*LEks.xy[i][j][k] - C2.xxy[i][j][k]*LEks.xx[i][j][k]
                                         - C2.yxx[i][j][k]*LEks.yy[i][j][k] - C2.yxy[i][j][k]*LEks.xy[i][j][k]
                                         - C2.zxx[i][j][k]*LEks.yz[i][j][k] - C2.zxy[i][j][k]*LEks.xz[i][j][k]
                                         - 2.0*(C2.xxy[i][j][k]*LEks.xx[i][j][k] + C2.yxy[i][j][k]*LEks.xy[i][j][k]
                                               +C2.zxy[i][j][k]*LEks.xz[i][j][k])) 
                          + gup.xz[i][j][k]*(dLEksxz.x[i][j][k] + dLEksxx.z[i][j][k] 
                                         - C2.xxx[i][j][k]*LEks.xz[i][j][k] - C2.xxz[i][j][k]*LEks.xx[i][j][k]
                                         - C2.yxx[i][j][k]*LEks.yz[i][j][k] - C2.yxz[i][j][k]*LEks.xy[i][j][k]
                                         - C2.zxx[i][j][k]*LEks.zz[i][j][k] - C2.zxz[i][j][k]*LEks.xz[i][j][k]
                                         - 2.0*(C2.xxz[i][j][k]*LEks.xx[i][j][k] + C2.yxz[i][j][k]*LEks.xy[i][j][k]
                                               +C2.zxz[i][j][k]*LEks.xz[i][j][k])) 
                          + gup.yz[i][j][k]*(dLEksxz.y[i][j][k] + dLEksxy.z[i][j][k] 
                                         - C2.xxy[i][j][k]*LEks.xz[i][j][k] - C2.xxz[i][j][k]*LEks.xy[i][j][k]
                                         - C2.yxy[i][j][k]*LEks.yz[i][j][k] - C2.yxz[i][j][k]*LEks.yy[i][j][k]
                                         - C2.zxy[i][j][k]*LEks.zz[i][j][k] - C2.zxz[i][j][k]*LEks.yz[i][j][k]
                                         - 2.0*(C2.xyz[i][j][k]*LEks.xx[i][j][k] + C2.yyz[i][j][k]*LEks.xy[i][j][k]
                                               +C2.zyz[i][j][k]*LEks.xz[i][j][k])) ;

        DelEks.y[i][j][k] = gup.xx[i][j][k]*(dLEksxy.x[i][j][k] - C2.xxy[i][j][k]*LEks.xx[i][j][k] 
                                         - C2.yxy[i][j][k]*LEks.xy[i][j][k] - C2.zxy[i][j][k]*LEks.xz[i][j][k]
                                         - C2.xxx[i][j][k]*LEks.xy[i][j][k] - C2.yxx[i][j][k]*LEks.yy[i][j][k]
                                         - C2.zxx[i][j][k]*LEks.yz[i][j][k])
                          + gup.yy[i][j][k]*(dLEksyy.y[i][j][k] - C2.xyy[i][j][k]*LEks.xy[i][j][k] 
                                         - C2.yyy[i][j][k]*LEks.yy[i][j][k] - C2.zyy[i][j][k]*LEks.yz[i][j][k]
                                         - C2.xyy[i][j][k]*LEks.xy[i][j][k] - C2.yyy[i][j][k]*LEks.yy[i][j][k]
                                         - C2.zyy[i][j][k]*LEks.yz[i][j][k])
                          + gup.zz[i][j][k]*(dLEksyz.z[i][j][k] - C2.xyz[i][j][k]*LEks.xz[i][j][k] 
                                         - C2.yyz[i][j][k]*LEks.yz[i][j][k] - C2.zyz[i][j][k]*LEks.zz[i][j][k]
                                         - C2.xzz[i][j][k]*LEks.xy[i][j][k] - C2.yzz[i][j][k]*LEks.yy[i][j][k]
                                         - C2.zzz[i][j][k]*LEks.yz[i][j][k])
                          + gup.xy[i][j][k]*(dLEksyy.x[i][j][k] + dLEksxy.y[i][j][k] 
                                         - C2.xxy[i][j][k]*LEks.xy[i][j][k] - C2.xyy[i][j][k]*LEks.xx[i][j][k]
                                         - C2.yxy[i][j][k]*LEks.yy[i][j][k] - C2.yyy[i][j][k]*LEks.xy[i][j][k]
                                         - C2.zxy[i][j][k]*LEks.yz[i][j][k] - C2.zyy[i][j][k]*LEks.xz[i][j][k]
                                         - 2.0*(C2.xxy[i][j][k]*LEks.xy[i][j][k] + C2.yxy[i][j][k]*LEks.yy[i][j][k]
                                               +C2.zxy[i][j][k]*LEks.yz[i][j][k])) 
                          + gup.xz[i][j][k]*(dLEksyz.x[i][j][k] + dLEksxy.z[i][j][k] 
                                         - C2.xxy[i][j][k]*LEks.xz[i][j][k] - C2.xyz[i][j][k]*LEks.xx[i][j][k]
                                         - C2.yxy[i][j][k]*LEks.yz[i][j][k] - C2.yyz[i][j][k]*LEks.xy[i][j][k]
                                         - C2.zxy[i][j][k]*LEks.zz[i][j][k] - C2.zyz[i][j][k]*LEks.xz[i][j][k]
                                         - 2.0*(C2.xxz[i][j][k]*LEks.xy[i][j][k] + C2.yxz[i][j][k]*LEks.yy[i][j][k]
                                               +C2.zxz[i][j][k]*LEks.yz[i][j][k])) 
                          + gup.yz[i][j][k]*(dLEksyz.y[i][j][k] + dLEksyy.z[i][j][k] 
                                         - C2.xyy[i][j][k]*LEks.xz[i][j][k] - C2.xyz[i][j][k]*LEks.xy[i][j][k]
                                         - C2.yyy[i][j][k]*LEks.yz[i][j][k] - C2.yyz[i][j][k]*LEks.yy[i][j][k]
                                         - C2.zyy[i][j][k]*LEks.zz[i][j][k] - C2.zyz[i][j][k]*LEks.yz[i][j][k]
                                         - 2.0*(C2.xyz[i][j][k]*LEks.xy[i][j][k] + C2.yyz[i][j][k]*LEks.yy[i][j][k]
                                               +C2.zyz[i][j][k]*LEks.yz[i][j][k])) ;
                         
        DelEks.z[i][j][k] = gup.xx[i][j][k]*(dLEksxz.x[i][j][k] - C2.xxz[i][j][k]*LEks.xx[i][j][k] 
                                         - C2.yxz[i][j][k]*LEks.xy[i][j][k] - C2.zxz[i][j][k]*LEks.xz[i][j][k]
                                         - C2.xxx[i][j][k]*LEks.xz[i][j][k] - C2.yxx[i][j][k]*LEks.yz[i][j][k]
                                         - C2.zxx[i][j][k]*LEks.zz[i][j][k])
                          + gup.yy[i][j][k]*(dLEksyz.y[i][j][k] - C2.xyz[i][j][k]*LEks.xy[i][j][k] 
                                         - C2.yyz[i][j][k]*LEks.yy[i][j][k] - C2.zyz[i][j][k]*LEks.yz[i][j][k]
                                         - C2.xyy[i][j][k]*LEks.xz[i][j][k] - C2.yyy[i][j][k]*LEks.yz[i][j][k]
                                         - C2.zyy[i][j][k]*LEks.zz[i][j][k])
                          + gup.zz[i][j][k]*(dLEkszz.z[i][j][k] - C2.xzz[i][j][k]*LEks.xz[i][j][k] 
                                         - C2.yzz[i][j][k]*LEks.yz[i][j][k] - C2.zzz[i][j][k]*LEks.zz[i][j][k]
                                         - C2.xzz[i][j][k]*LEks.xz[i][j][k] - C2.yzz[i][j][k]*LEks.yz[i][j][k]
                                         - C2.zzz[i][j][k]*LEks.zz[i][j][k])
                          + gup.xy[i][j][k]*(dLEksyz.x[i][j][k] + dLEksxz.y[i][j][k] 
                                         - C2.xxz[i][j][k]*LEks.xy[i][j][k] - C2.xyz[i][j][k]*LEks.xx[i][j][k]
                                         - C2.yxz[i][j][k]*LEks.yy[i][j][k] - C2.yyz[i][j][k]*LEks.xy[i][j][k]
                                         - C2.zxz[i][j][k]*LEks.yz[i][j][k] - C2.zyz[i][j][k]*LEks.xz[i][j][k]
                                         - 2.0*(C2.xxy[i][j][k]*LEks.xz[i][j][k] + C2.yxy[i][j][k]*LEks.yz[i][j][k]
                                               +C2.zxy[i][j][k]*LEks.zz[i][j][k])) 
                          + gup.xz[i][j][k]*(dLEkszz.x[i][j][k] + dLEksxz.z[i][j][k] 
                                         - C2.xxz[i][j][k]*LEks.xz[i][j][k] - C2.xzz[i][j][k]*LEks.xx[i][j][k]
                                         - C2.yxz[i][j][k]*LEks.yz[i][j][k] - C2.yzz[i][j][k]*LEks.xy[i][j][k]
                                         - C2.zxz[i][j][k]*LEks.zz[i][j][k] - C2.zzz[i][j][k]*LEks.xz[i][j][k]
                                         - 2.0*(C2.xxz[i][j][k]*LEks.xz[i][j][k] + C2.yxz[i][j][k]*LEks.yz[i][j][k]
                                               +C2.zxz[i][j][k]*LEks.zz[i][j][k])) 
                          + gup.yz[i][j][k]*(dLEkszz.y[i][j][k] + dLEksyz.z[i][j][k] 
                                         - C2.xyz[i][j][k]*LEks.xz[i][j][k] - C2.xzz[i][j][k]*LEks.xy[i][j][k]
                                         - C2.yyz[i][j][k]*LEks.yz[i][j][k] - C2.yzz[i][j][k]*LEks.yy[i][j][k]
                                         - C2.zyz[i][j][k]*LEks.zz[i][j][k] - C2.zzz[i][j][k]*LEks.yz[i][j][k]
                                         - 2.0*(C2.xyz[i][j][k]*LEks.xz[i][j][k] + C2.yyz[i][j][k]*LEks.yz[i][j][k]
                                               +C2.zyz[i][j][k]*LEks.zz[i][j][k])) ;
 }}}



 free_dVector3D(dLEksxx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksxx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksxx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksxy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksxy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksxy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksxz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksxz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksxz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksyy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksyy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksyy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksyz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksyz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEksyz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEkszz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEkszz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(dLEkszz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

 return;
}

