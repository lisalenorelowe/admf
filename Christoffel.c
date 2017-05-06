#include "DeclareFunctions.h"

void Christoffel(Tensor gup,Vector dgxx,Vector dgyy,Vector dgzz,Vector dgxy,Vector dgxz,Vector dgyz,Params Par,Connection C2){

 int i,j,k,m;
 double ***C1xxx,***C1xxy,***C1xxz,***C1xyy,***C1xyz,***C1xzz,***C1yxx,***C1yxy;
 double ***C1yxz,***C1yyy,***C1yyz,***C1yzz,***C1zxx,***C1zxy,***C1zxz,***C1zyy;
 double ***C1zyz,***C1zzz;
/*----------------------------------------------------------------*/

C1xxx= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1xxy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1xxz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1yxy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1xyy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1xyz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1yxz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1zyy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1xzz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1zxy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1yyy= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1yyz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1yxx= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1zxx= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1yzz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1zxz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1zyz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
C1zzz= dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
       C1xxx[i][j][k] = 0.5*dgxx.x[i][j][k];
       C1xxy[i][j][k] = 0.5*dgxx.y[i][j][k];
       C1xxz[i][j][k] = 0.5*dgxx.z[i][j][k];
       C1xyy[i][j][k] = dgxy.y[i][j][k] - 0.5*dgyy.x[i][j][k];
       C1xyz[i][j][k] = 0.5*dgxy.z[i][j][k] + 0.5*dgxz.y[i][j][k] - 0.5*dgyz.x[i][j][k];
       C1xzz[i][j][k] = dgxz.z[i][j][k] - 0.5*dgzz.x[i][j][k];

       C1yxx[i][j][k] = dgxy.x[i][j][k] - 0.5*dgxx.y[i][j][k];
       C1yxy[i][j][k] = 0.5*dgyy.x[i][j][k];
       C1yxz[i][j][k] = 0.5*dgxy.z[i][j][k] + 0.5*dgyz.x[i][j][k] - 0.5*dgxz.y[i][j][k];
       C1yyy[i][j][k] = 0.5*dgyy.y[i][j][k];
       C1yyz[i][j][k] = 0.5*dgyy.z[i][j][k];
       C1yzz[i][j][k] = dgyz.z[i][j][k] - 0.5*dgzz.y[i][j][k];

       C1zxx[i][j][k] = dgxz.x[i][j][k] - 0.5*dgxx.z[i][j][k];
       C1zxy[i][j][k] = 0.5*dgxz.y[i][j][k] + 0.5*dgyz.x[i][j][k] - 0.5*dgxy.z[i][j][k];
       C1zxz[i][j][k] = 0.5*dgzz.x[i][j][k];
       C1zyy[i][j][k] = dgyz.y[i][j][k] - 0.5*dgyy.z[i][j][k];
       C1zyz[i][j][k] = 0.5*dgzz.y[i][j][k];
       C1zzz[i][j][k] = 0.5*dgzz.z[i][j][k];

       C2.xxx[i][j][k] = gup.xx[i][j][k]*C1xxx[i][j][k]+gup.xy[i][j][k]*C1yxx[i][j][k]+gup.xz[i][j][k]*C1zxx[i][j][k]; 

       C2.xxy[i][j][k] = gup.xx[i][j][k]*C1xxy[i][j][k]+gup.xy[i][j][k]*C1yxy[i][j][k]+gup.xz[i][j][k]*C1zxy[i][j][k];

       C2.xxz[i][j][k] = gup.xx[i][j][k]*C1xxz[i][j][k]+gup.xy[i][j][k]*C1yxz[i][j][k]+gup.xz[i][j][k]*C1zxz[i][j][k];

       C2.xyy[i][j][k] = gup.xx[i][j][k]*C1xyy[i][j][k]+gup.xy[i][j][k]*C1yyy[i][j][k]+gup.xz[i][j][k]*C1zyy[i][j][k];

       C2.xyz[i][j][k] = gup.xx[i][j][k]*C1xyz[i][j][k]+gup.xy[i][j][k]*C1yyz[i][j][k]+gup.xz[i][j][k]*C1zyz[i][j][k]; 

       C2.xzz[i][j][k] = gup.xx[i][j][k]*C1xzz[i][j][k]+gup.xy[i][j][k]*C1yzz[i][j][k]+gup.xz[i][j][k]*C1zzz[i][j][k];

       C2.yxx[i][j][k] = gup.xy[i][j][k]*C1xxx[i][j][k]+gup.yy[i][j][k]*C1yxx[i][j][k]+gup.yz[i][j][k]*C1zxx[i][j][k];

       C2.yxy[i][j][k] = gup.xy[i][j][k]*C1xxy[i][j][k]+gup.yy[i][j][k]*C1yxy[i][j][k]+gup.yz[i][j][k]*C1zxy[i][j][k];

       C2.yxz[i][j][k] = gup.xy[i][j][k]*C1xxz[i][j][k]+gup.yy[i][j][k]*C1yxz[i][j][k]+gup.yz[i][j][k]*C1zxz[i][j][k];

       C2.yyy[i][j][k] = gup.xy[i][j][k]*C1xyy[i][j][k]+gup.yy[i][j][k]*C1yyy[i][j][k]+gup.yz[i][j][k]*C1zyy[i][j][k];

       C2.yyz[i][j][k] = gup.xy[i][j][k]*C1xyz[i][j][k]+gup.yy[i][j][k]*C1yyz[i][j][k]+gup.yz[i][j][k]*C1zyz[i][j][k];

       C2.yzz[i][j][k] = gup.xy[i][j][k]*C1xzz[i][j][k]+gup.yy[i][j][k]*C1yzz[i][j][k]+gup.yz[i][j][k]*C1zzz[i][j][k];

       C2.zxx[i][j][k] = gup.xz[i][j][k]*C1xxx[i][j][k]+gup.yz[i][j][k]*C1yxx[i][j][k]+gup.zz[i][j][k]*C1zxx[i][j][k];

       C2.zxy[i][j][k] = gup.xz[i][j][k]*C1xxy[i][j][k]+gup.yz[i][j][k]*C1yxy[i][j][k]+gup.zz[i][j][k]*C1zxy[i][j][k];

       C2.zxz[i][j][k] = gup.xz[i][j][k]*C1xxz[i][j][k]+gup.yz[i][j][k]*C1yxz[i][j][k]+gup.zz[i][j][k]*C1zxz[i][j][k];

       C2.zyy[i][j][k] = gup.xz[i][j][k]*C1xyy[i][j][k]+gup.yz[i][j][k]*C1yyy[i][j][k]+gup.zz[i][j][k]*C1zyy[i][j][k];

       C2.zyz[i][j][k] = gup.xz[i][j][k]*C1xyz[i][j][k]+gup.yz[i][j][k]*C1yyz[i][j][k]+gup.zz[i][j][k]*C1zyz[i][j][k];

       C2.zzz[i][j][k] = gup.xz[i][j][k]*C1xzz[i][j][k]+gup.yz[i][j][k]*C1yzz[i][j][k]+gup.zz[i][j][k]*C1zzz[i][j][k];

}}}

free_dVector3D(C1xxx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1xxy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1xxz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1yxy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1xyy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1xyz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1yxz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1zyy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1xzz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1zxy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1yyy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1yyz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1yxx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1zxx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1yzz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1zxz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1zyz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C1zzz,0,Par.nxb,0,Par.nyb,0,Par.nzb);


return;
}
