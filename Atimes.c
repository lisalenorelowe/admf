#include "DeclareFunctions.h"

void Atimes(Scalar Phi,Vector Eks,Scalar Atimes0,Vector Atimes,double Coeff0[XMAX][YMAX][ZMAX],double Coeffx[XMAX][YMAX][ZMAX],double Coeffy[XMAX][YMAX][ZMAX],double Coeffz[XMAX][YMAX][ZMAX],Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,Connection C2,Params Par){

 int i, j, k;
 Scalar LapPhi;
 Tensor LEks;
 Vector DelEks;

 LapPhi.s = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 LEks.xx = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 LEks.xy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 LEks.xz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 LEks.yy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 LEks.yz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 LEks.zz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 DelEks.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 DelEks.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 DelEks.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

 Laplacian(Phi,LapPhi,gup,C2,Par);
 Loperator(Eks,LEks,gij,gup,C2,Par);
 Deltaop(LEks,DelEks,gij,gup,C2,Par);

 for(i=0;i<Par.nxb;i++){
   for(j=0;j<Par.nyb;j++){
     for(k=0;k<Par.nzb;k++){
       Atimes0.s[i][j][k] = Coeff0[i][j][k]*Phi.s[i][j][k] - 8.0*LapPhi.s[i][j][k] 
                      - 2.0*(Kup.xx[i][j][k]*LEks.xx[i][j][k] + Kup.yy[i][j][k]*LEks.yy[i][j][k] 
                            +Kup.zz[i][j][k]*LEks.zz[i][j][k] + 2.0*(Kup.xy[i][j][k]*LEks.xy[i][j][k]
                            +Kup.xz[i][j][k]*LEks.xz[i][j][k] + Kup.yz[i][j][k]*LEks.yz[i][j][k]));
       Atimes.x[i][j][k] = Coeffx[i][j][k]*Phi.s[i][j][k] + DelEks.x[i][j][k];
       Atimes.y[i][j][k] = Coeffy[i][j][k]*Phi.s[i][j][k] + DelEks.y[i][j][k];
       Atimes.z[i][j][k] = Coeffz[i][j][k]*Phi.s[i][j][k] + DelEks.z[i][j][k];
 }}}


 free_dVector3D(LapPhi.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(LEks.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(LEks.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(LEks.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(LEks.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(LEks.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(LEks.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(DelEks.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(DelEks.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(DelEks.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

return;
}

