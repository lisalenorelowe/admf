#include "DeclareFunctions.h"

void GetRicci(Tensor gij,Tensor gup,Tensor *Ricci,Params Par,Connection C2){
 int i,j,k,m;
 double ***Riccixx,***Ricciyy,***Riccizz,***Riccixy,***Riccixz,***Ricciyz;
 Tensor d2gxx,d2gyy,d2gzz,d2gxy,d2gxz,d2gyz;
 Vector dgxx,dgyy,dgzz,dgxy,dgxz,dgyz;
/*----------------------------------------------------------------*/

Riccixx=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Ricciyy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Riccizz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Riccixy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Riccixz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
Ricciyz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

d2gxx.xx=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxx.yy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxx.zz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxx.xy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxx.xz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxx.yz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgxx.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgxx.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgxx.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyy.xx=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyy.yy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyy.zz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyy.xy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyy.xz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyy.yz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgyy.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgyy.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgyy.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gzz.xx=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gzz.yy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gzz.zz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gzz.xy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gzz.xz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gzz.yz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgzz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgzz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgzz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxy.xx=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxy.yy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxy.zz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxy.xy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxy.xz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxy.yz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgxy.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgxy.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgxy.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyz.xx=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyz.yy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyz.zz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyz.xy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyz.xz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gyz.yz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgyz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgyz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgyz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxz.xx=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxz.yy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxz.zz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxz.xy=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxz.xz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);d2gxz.yz=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgxz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgxz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);dgxz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Find the Ricci Tensor...
FunctionDerivs(gij.xx,&d2gxx,Par);
FunctionFirstDerivs(gij.xx,&dgxx,Par);

FunctionDerivs(gij.yy,&d2gyy,Par);
FunctionFirstDerivs(gij.yy,&dgyy,Par);

FunctionDerivs(gij.zz,&d2gzz,Par);
FunctionFirstDerivs(gij.zz,&dgzz,Par);

FunctionDerivs(gij.xy,&d2gxy,Par);
FunctionFirstDerivs(gij.xy,&dgxy,Par);

FunctionDerivs(gij.xz,&d2gxz,Par);
FunctionFirstDerivs(gij.xz,&dgxz,Par);

FunctionDerivs(gij.yz,&d2gyz,Par);
FunctionFirstDerivs(gij.yz,&dgyz,Par);

//printf("Dgxx at origin = %2.10f %2.10f %2.10f \n", dgxx.x[Par.nxb/2][Par.nyb/2][Par.nzb/2], 
//   dgxx.y[Par.nxb/2][Par.nyb/2][Par.nzb/2], dgxx.z[Par.nxb/2][Par.nyb/2][Par.nzb/2]);
//printf("DDgxx at origin = %2.10f %2.10f %2.10f \n", d2gxx.xx[Par.nxb/2][Par.nyb/2][Par.nzb/2], 
//   d2gxx.yy[Par.nxb/2][Par.nyb/2][Par.nzb/2], d2gxx.zz[Par.nxb/2][Par.nyb/2][Par.nzb/2]);


Rxx(Riccixx,gij,d2gxx,d2gyy,d2gzz,d2gxy,d2gxz,d2gyz,dgxx,dgyy,dgzz,dgxy,dgxz,dgyz,Par);
Ryy(Ricciyy,gij,d2gxx,d2gyy,d2gzz,d2gxy,d2gxz,d2gyz,dgxx,dgyy,dgzz,dgxy,dgxz,dgyz,Par);
Rzz(Riccizz,gij,d2gxx,d2gyy,d2gzz,d2gxy,d2gxz,d2gyz,dgxx,dgyy,dgzz,dgxy,dgxz,dgyz,Par);
Rxy(Riccixy,gij,d2gxx,d2gyy,d2gzz,d2gxy,d2gxz,d2gyz,dgxx,dgyy,dgzz,dgxy,dgxz,dgyz,Par);
Rxz(Riccixz,gij,d2gxx,d2gyy,d2gzz,d2gxy,d2gxz,d2gyz,dgxx,dgyy,dgzz,dgxy,dgxz,dgyz,Par);
Ryz(Ricciyz,gij,d2gxx,d2gyy,d2gzz,d2gxy,d2gxz,d2gyz,dgxx,dgyy,dgzz,dgxy,dgxz,dgyz,Par);

for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
     Ricci->xx[i][j][k]= Riccixx[i][j][k];
     Ricci->yy[i][j][k]= Ricciyy[i][j][k];
     Ricci->zz[i][j][k]= Riccizz[i][j][k];
     Ricci->xy[i][j][k]= Riccixy[i][j][k];
     Ricci->xz[i][j][k]= Riccixz[i][j][k];
     Ricci->yz[i][j][k]= Ricciyz[i][j][k];
}}}


Christoffel(gup,dgxx,dgyy,dgzz,dgxy,dgxz,dgyz,Par,C2);


free_dVector3D(d2gxx.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxx.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxx.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxx.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxx.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxx.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

free_dVector3D(d2gyy.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gyy.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gyy.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gyy.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gyy.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gyy.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

free_dVector3D(d2gzz.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gzz.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gzz.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gzz.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gzz.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gzz.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgzz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgzz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgzz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

free_dVector3D(d2gxy.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxy.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxy.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxy.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxy.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxy.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

free_dVector3D(d2gxz.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxz.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxz.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxz.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxz.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gxz.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

free_dVector3D(d2gyz.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gyz.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gyz.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gyz.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gyz.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2gyz.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

free_dVector3D(Riccixx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricciyy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Riccizz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Riccixy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Riccixz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricciyz,0,Par.nxb,0,Par.nyb,0,Par.nzb);

return;
}
