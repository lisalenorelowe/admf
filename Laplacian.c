#include "DeclareFunctions.h"

void Laplacian(Scalar Phi, Scalar LapPhi, Tensor gup, Connection C2, Params Par){

 int i, j, k;
 Vector dPhi;
 Tensor d2Phi;

d2Phi.xx     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
d2Phi.yy     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
d2Phi.zz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
d2Phi.xy     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
d2Phi.xz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
d2Phi.yz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dPhi.x     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dPhi.y     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dPhi.z     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

 FunctionFirstDerivs(Phi.s,&dPhi,Par);
 FunctionDerivs(Phi.s,&d2Phi,Par);

 for(i=0;i<Par.nxb;i++){
   for(j=0;j<Par.nyb;j++){
      for(k=0;k<Par.nzb;k++){
         LapPhi.s[i][j][k] = gup.xx[i][j][k]*(d2Phi.xx[i][j][k] - C2.xxx[i][j][k]*dPhi.x[i][j][k] 
                                                             - C2.yxx[i][j][k]*dPhi.y[i][j][k]
                                                             - C2.zxx[i][j][k]*dPhi.z[i][j][k])
                        + gup.yy[i][j][k]*(d2Phi.yy[i][j][k] - C2.xyy[i][j][k]*dPhi.x[i][j][k] 
                                                             - C2.yyy[i][j][k]*dPhi.y[i][j][k]
                                                             - C2.zyy[i][j][k]*dPhi.z[i][j][k])
                        + gup.zz[i][j][k]*(d2Phi.zz[i][j][k] - C2.xzz[i][j][k]*dPhi.x[i][j][k] 
                                                             - C2.yzz[i][j][k]*dPhi.y[i][j][k]
                                                             - C2.zzz[i][j][k]*dPhi.z[i][j][k])
                        + 2.0*gup.xy[i][j][k]*(d2Phi.xy[i][j][k] - C2.xxy[i][j][k]*dPhi.x[i][j][k]
                                                                 - C2.yxy[i][j][k]*dPhi.y[i][j][k]
                                                                 - C2.zxy[i][j][k]*dPhi.z[i][j][k])
                        + 2.0*gup.xz[i][j][k]*(d2Phi.xz[i][j][k] - C2.xxz[i][j][k]*dPhi.x[i][j][k]
                                                                 - C2.yxz[i][j][k]*dPhi.y[i][j][k]
                                                                 - C2.zxz[i][j][k]*dPhi.z[i][j][k])
                        + 2.0*gup.yz[i][j][k]*(d2Phi.yz[i][j][k] - C2.xyz[i][j][k]*dPhi.x[i][j][k]
                                                                 - C2.yyz[i][j][k]*dPhi.y[i][j][k]
                                                                 - C2.zyz[i][j][k]*dPhi.z[i][j][k]);
 }}}

free_dVector3D(d2Phi.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2Phi.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2Phi.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2Phi.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2Phi.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2Phi.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dPhi.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dPhi.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dPhi.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

 return;
}

