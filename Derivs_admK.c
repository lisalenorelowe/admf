#include "DeclareFunctions.h"

void Derivs_admK(double t,Scalar alpha,Vector shift,Tensor gij,Tensor gup,Tensor Kij,Scalar trK, Tensor *dKij,Params Par,Connection C2,Tensor Ricci){
 int i,j,k,m;
 double x,y,z;
 double ***dummy;
 Vector dalpha,dbx,dby,dbz,dKxx,dKyy,dKzz,dKxy,dKxz,dKyz;
 Tensor d2alpha;
/*
Call ginv, HamCon, and MomCon before calling Derivs_admK.
*/

dummy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
//K = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

d2alpha.xx     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
d2alpha.yy     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
d2alpha.zz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
d2alpha.xy     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
d2alpha.xz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
d2alpha.yz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

dalpha.x     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dalpha.y     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dalpha.z     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbx.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbx.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbx.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dby.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dby.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dby.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

dKxx.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKxx.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKxx.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKyy.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKyy.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKyy.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKzz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKzz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKzz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKxy.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKxy.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKxy.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKyz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKyz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKyz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKxz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKxz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dKxz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//alpha derivs
 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
   for(k=0;k<Par.nzb;k++){
    dummy[i][j][k] = alpha.s[i][j][k];
 }}}
 FunctionFirstDerivs(dummy,&dalpha,Par);
 FunctionDerivs(dummy,&d2alpha,Par);


// Get Derivatives of Shift Vector \beta^i and Kij
 if(Par.ishift!=0){    // But not if shift is zero

// (shift derivs)
 FunctionFirstDerivs(shift.x,&dbx,Par);
 FunctionFirstDerivs(shift.y,&dby,Par);
 FunctionFirstDerivs(shift.z,&dbz,Par);
// End Shift Derivs

//Get Derivs of Kij
FunctionFirstDerivs(Kij.xx,&dKxx,Par);
FunctionFirstDerivs(Kij.yy,&dKyy,Par);
FunctionFirstDerivs(Kij.zz,&dKzz,Par);
FunctionFirstDerivs(Kij.xy,&dKxy,Par);
FunctionFirstDerivs(Kij.xz,&dKxz,Par);
FunctionFirstDerivs(Kij.yz,&dKyz,Par);
//End Kij Derivs

} //End if Shift=0 statment


for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
        trK.s[i][j][k] = gup.xx[i][j][k]*Kij.xx[i][j][k] + gup.yy[i][j][k]*Kij.yy[i][j][k] +gup.zz[i][j][k]*Kij.zz[i][j][k] +
        2.*gup.xy[i][j][k]*Kij.xy[i][j][k] +2.*gup.xz[i][j][k]*Kij.xz[i][j][k] +2.*gup.yz[i][j][k]*Kij.yz[i][j][k] ;

	dKij->xx[i][j][k] =  alpha.s[i][j][k]*trK.s[i][j][k]*Kij.xx[i][j][k]; 
        dKij->yy[i][j][k] =  alpha.s[i][j][k]*trK.s[i][j][k]*Kij.yy[i][j][k]; 
        dKij->zz[i][j][k] =  alpha.s[i][j][k]*trK.s[i][j][k]*Kij.zz[i][j][k];  
        dKij->xy[i][j][k] =  alpha.s[i][j][k]*trK.s[i][j][k]*Kij.xy[i][j][k]; 
        dKij->xz[i][j][k] =  alpha.s[i][j][k]*trK.s[i][j][k]*Kij.xz[i][j][k]; 
        dKij->yz[i][j][k] =  alpha.s[i][j][k]*trK.s[i][j][k]*Kij.yz[i][j][k]; 

}}}


for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
        dKij->xx[i][j][k] = dKij->xx[i][j][k] -2.*alpha.s[i][j][k]*(
Kij.xx[i][j][k]*gup.xx[i][j][k]*Kij.xx[i][j][k]+Kij.xx[i][j][k]*gup.xy[i][j][k]*Kij.xy[i][j][k]+Kij.xx[i][j][k]*gup.xz[i][j][k]*Kij.xz[i][j][k]    +
Kij.xy[i][j][k]*gup.xy[i][j][k]*Kij.xx[i][j][k]+Kij.xy[i][j][k]*gup.yy[i][j][k]*Kij.xy[i][j][k]+Kij.xy[i][j][k]*gup.yz[i][j][k]*Kij.xz[i][j][k]    +
Kij.xz[i][j][k]*gup.xz[i][j][k]*Kij.xx[i][j][k]+Kij.xz[i][j][k]*gup.yz[i][j][k]*Kij.xy[i][j][k]+Kij.xz[i][j][k]*gup.zz[i][j][k]*Kij.xz[i][j][k]);


        dKij->yy[i][j][k] = dKij->yy[i][j][k] -2.*alpha.s[i][j][k]*(
Kij.xy[i][j][k]*gup.xx[i][j][k]*Kij.xy[i][j][k]+Kij.xy[i][j][k]*gup.xy[i][j][k]*Kij.yy[i][j][k]+Kij.xy[i][j][k]*gup.xz[i][j][k]*Kij.yz[i][j][k]    +
Kij.yy[i][j][k]*gup.xy[i][j][k]*Kij.xy[i][j][k]+Kij.yy[i][j][k]*gup.yy[i][j][k]*Kij.yy[i][j][k]+Kij.yy[i][j][k]*gup.yz[i][j][k]*Kij.yz[i][j][k]    +
Kij.yz[i][j][k]*gup.xz[i][j][k]*Kij.xy[i][j][k]+Kij.yz[i][j][k]*gup.yz[i][j][k]*Kij.yy[i][j][k]+Kij.yz[i][j][k]*gup.zz[i][j][k]*Kij.yz[i][j][k]);

        dKij->zz[i][j][k] = dKij->zz[i][j][k] -2.*alpha.s[i][j][k]*(
Kij.xz[i][j][k]*gup.xx[i][j][k]*Kij.xz[i][j][k]+Kij.xz[i][j][k]*gup.xy[i][j][k]*Kij.yz[i][j][k]+Kij.xz[i][j][k]*gup.xz[i][j][k]*Kij.zz[i][j][k]    +
Kij.yz[i][j][k]*gup.xy[i][j][k]*Kij.xz[i][j][k]+Kij.yz[i][j][k]*gup.yy[i][j][k]*Kij.yz[i][j][k]+Kij.yz[i][j][k]*gup.yz[i][j][k]*Kij.zz[i][j][k]    +
Kij.zz[i][j][k]*gup.xz[i][j][k]*Kij.xz[i][j][k]+Kij.zz[i][j][k]*gup.yz[i][j][k]*Kij.yz[i][j][k]+Kij.zz[i][j][k]*gup.zz[i][j][k]*Kij.zz[i][j][k]);

        dKij->xy[i][j][k] = dKij->xy[i][j][k] -2.*alpha.s[i][j][k]*(
Kij.xx[i][j][k]*gup.xx[i][j][k]*Kij.xy[i][j][k]+Kij.xx[i][j][k]*gup.xy[i][j][k]*Kij.yy[i][j][k]+Kij.xx[i][j][k]*gup.xz[i][j][k]*Kij.yz[i][j][k]    +
Kij.xy[i][j][k]*gup.xy[i][j][k]*Kij.xy[i][j][k]+Kij.xy[i][j][k]*gup.yy[i][j][k]*Kij.yy[i][j][k]+Kij.xy[i][j][k]*gup.yz[i][j][k]*Kij.yz[i][j][k]    +
Kij.xz[i][j][k]*gup.xz[i][j][k]*Kij.xy[i][j][k]+Kij.xz[i][j][k]*gup.yz[i][j][k]*Kij.yy[i][j][k]+Kij.xz[i][j][k]*gup.zz[i][j][k]*Kij.yz[i][j][k]);

        dKij->xz[i][j][k] = dKij->xz[i][j][k] -2.*alpha.s[i][j][k]*(
Kij.xx[i][j][k]*gup.xx[i][j][k]*Kij.xz[i][j][k]+Kij.xx[i][j][k]*gup.xy[i][j][k]*Kij.yz[i][j][k]+Kij.xx[i][j][k]*gup.xz[i][j][k]*Kij.zz[i][j][k]    +
Kij.xy[i][j][k]*gup.xy[i][j][k]*Kij.xz[i][j][k]+Kij.xy[i][j][k]*gup.yy[i][j][k]*Kij.yz[i][j][k]+Kij.xy[i][j][k]*gup.yz[i][j][k]*Kij.zz[i][j][k]    +
Kij.xz[i][j][k]*gup.xz[i][j][k]*Kij.xz[i][j][k]+Kij.xz[i][j][k]*gup.yz[i][j][k]*Kij.yz[i][j][k]+Kij.xz[i][j][k]*gup.zz[i][j][k]*Kij.zz[i][j][k]);

        dKij->yz[i][j][k] = dKij->yz[i][j][k] -2.*alpha.s[i][j][k]*(
Kij.xy[i][j][k]*gup.xx[i][j][k]*Kij.xz[i][j][k]+Kij.xy[i][j][k]*gup.xy[i][j][k]*Kij.yz[i][j][k]+Kij.xy[i][j][k]*gup.xz[i][j][k]*Kij.zz[i][j][k]    +
Kij.yy[i][j][k]*gup.xy[i][j][k]*Kij.xz[i][j][k]+Kij.yy[i][j][k]*gup.yy[i][j][k]*Kij.yz[i][j][k]+Kij.yy[i][j][k]*gup.yz[i][j][k]*Kij.zz[i][j][k]    +
Kij.yz[i][j][k]*gup.xz[i][j][k]*Kij.xz[i][j][k]+Kij.yz[i][j][k]*gup.yz[i][j][k]*Kij.yz[i][j][k]+Kij.yz[i][j][k]*gup.zz[i][j][k]*Kij.zz[i][j][k]);

}}}



//Now find the Ricci Tensor...
//GetRicci(gij,gup,&Ricci,Par,C2);


for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
        dKij->xx[i][j][k] = dKij->xx[i][j][k] + alpha.s[i][j][k]*Ricci.xx[i][j][k];
        dKij->yy[i][j][k] = dKij->yy[i][j][k] + alpha.s[i][j][k]*Ricci.yy[i][j][k];
        dKij->zz[i][j][k] = dKij->zz[i][j][k] + alpha.s[i][j][k]*Ricci.zz[i][j][k];
        dKij->xy[i][j][k] = dKij->xy[i][j][k] + alpha.s[i][j][k]*Ricci.xy[i][j][k];
        dKij->xz[i][j][k] = dKij->xz[i][j][k] + alpha.s[i][j][k]*Ricci.xz[i][j][k];
        dKij->yz[i][j][k] = dKij->yz[i][j][k] + alpha.s[i][j][k]*Ricci.yz[i][j][k];

}}}

//Add dalpha terms:
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
        dKij->xx[i][j][k] = dKij->xx[i][j][k] -(d2alpha.xx[i][j][k]-dalpha.x[i][j][k]*C2.xxx[i][j][k]-dalpha.y[i][j][k]*C2.yxx[i][j][k]-dalpha.z[i][j][k]*C2.zxx[i][j][k]);
        dKij->yy[i][j][k] = dKij->yy[i][j][k] -(d2alpha.yy[i][j][k]-dalpha.x[i][j][k]*C2.xyy[i][j][k]-dalpha.y[i][j][k]*C2.yyy[i][j][k]-dalpha.z[i][j][k]*C2.zyy[i][j][k]);
        dKij->zz[i][j][k] = dKij->zz[i][j][k] -(d2alpha.zz[i][j][k]-dalpha.x[i][j][k]*C2.xzz[i][j][k]-dalpha.y[i][j][k]*C2.yzz[i][j][k]-dalpha.z[i][j][k]*C2.zzz[i][j][k]);
        dKij->xy[i][j][k] = dKij->xy[i][j][k] -(d2alpha.xy[i][j][k]-dalpha.x[i][j][k]*C2.xxy[i][j][k]-dalpha.y[i][j][k]*C2.yxy[i][j][k]-dalpha.z[i][j][k]*C2.zxy[i][j][k]);
        dKij->xz[i][j][k] = dKij->xz[i][j][k] -(d2alpha.xz[i][j][k]-dalpha.x[i][j][k]*C2.xxz[i][j][k]-dalpha.y[i][j][k]*C2.yxz[i][j][k]-dalpha.z[i][j][k]*C2.zxz[i][j][k]);
        dKij->yz[i][j][k] = dKij->yz[i][j][k] -(d2alpha.yz[i][j][k]-dalpha.x[i][j][k]*C2.xyz[i][j][k]-dalpha.y[i][j][k]*C2.yyz[i][j][k]-dalpha.z[i][j][k]*C2.zyz[i][j][k]);
}}}

if(Par.ishift!=0){
/* Add shift parts:
    $ \partial K_{ij} = -D_iD_j\alpha + \alpha(R_{ij} + KK_{ij} - 2K_{im}K^m_j) 
>>>  + \beta^k\partial_kK_{ij} + K_{ik}\partial_j\beta^k + K_{kj}\partial_i\beta^k $ <<<
*/

// \beta^k \partial_k K_{ij}
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
        dKij->xx[i][j][k] = dKij->xx[i][j][k] + shift.x[i][j][k]*dKxx.x[i][j][k] + shift.y[i][j][k]*dKxx.y[i][j][k] + shift.z[i][j][k]*dKxx.z[i][j][k];
        dKij->yy[i][j][k] = dKij->yy[i][j][k] + shift.x[i][j][k]*dKyy.x[i][j][k] + shift.y[i][j][k]*dKyy.y[i][j][k] + shift.z[i][j][k]*dKyy.z[i][j][k];
        dKij->zz[i][j][k] = dKij->zz[i][j][k] + shift.x[i][j][k]*dKzz.x[i][j][k] + shift.y[i][j][k]*dKzz.y[i][j][k] + shift.z[i][j][k]*dKzz.z[i][j][k];
        dKij->xy[i][j][k] = dKij->xy[i][j][k] + shift.x[i][j][k]*dKxy.x[i][j][k] + shift.y[i][j][k]*dKxy.y[i][j][k] + shift.z[i][j][k]*dKxy.z[i][j][k];
        dKij->xz[i][j][k] = dKij->xz[i][j][k] + shift.x[i][j][k]*dKxz.x[i][j][k] + shift.y[i][j][k]*dKxz.y[i][j][k] + shift.z[i][j][k]*dKxz.z[i][j][k];
        dKij->yz[i][j][k] = dKij->yz[i][j][k] + shift.x[i][j][k]*dKyz.x[i][j][k] + shift.y[i][j][k]*dKyz.y[i][j][k] + shift.z[i][j][k]*dKyz.z[i][j][k];
}}}

// K_{ik} \partial_j \beta^k
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
        dKij->xx[i][j][k] = dKij->xx[i][j][k] + Kij.xx[i][j][k]*dbx.x[i][j][k] + Kij.xy[i][j][k]*dby.x[i][j][k] + Kij.xz[i][j][k]*dbz.x[i][j][k];
        dKij->yy[i][j][k] = dKij->yy[i][j][k] + Kij.xy[i][j][k]*dbx.y[i][j][k] + Kij.yy[i][j][k]*dby.y[i][j][k] + Kij.yz[i][j][k]*dbz.y[i][j][k];
        dKij->zz[i][j][k] = dKij->zz[i][j][k] + Kij.xz[i][j][k]*dbx.z[i][j][k] + Kij.yz[i][j][k]*dby.z[i][j][k] + Kij.zz[i][j][k]*dbz.z[i][j][k];
        dKij->xy[i][j][k] = dKij->xy[i][j][k] + Kij.xx[i][j][k]*dbx.y[i][j][k] + Kij.xy[i][j][k]*dby.y[i][j][k] + Kij.xz[i][j][k]*dbz.y[i][j][k];
        dKij->xz[i][j][k] = dKij->xz[i][j][k] + Kij.xx[i][j][k]*dbx.z[i][j][k] + Kij.xy[i][j][k]*dby.z[i][j][k] + Kij.xz[i][j][k]*dbz.z[i][j][k];
        dKij->yz[i][j][k] = dKij->yz[i][j][k] + Kij.xy[i][j][k]*dbx.z[i][j][k] + Kij.yy[i][j][k]*dby.z[i][j][k] + Kij.yz[i][j][k]*dbz.z[i][j][k];
}}}

//K_{kj}\partial_i\beta^k
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
        dKij->xx[i][j][k] = dKij->xx[i][j][k] + Kij.xx[i][j][k]*dbx.x[i][j][k] + Kij.xy[i][j][k]*dby.x[i][j][k] + Kij.xz[i][j][k]*dbz.x[i][j][k]; 
        dKij->yy[i][j][k] = dKij->yy[i][j][k] + Kij.xy[i][j][k]*dbx.y[i][j][k] + Kij.yy[i][j][k]*dby.y[i][j][k] + Kij.yz[i][j][k]*dbz.y[i][j][k];
        dKij->zz[i][j][k] = dKij->zz[i][j][k] + Kij.xz[i][j][k]*dbx.z[i][j][k] + Kij.yz[i][j][k]*dby.z[i][j][k] + Kij.zz[i][j][k]*dbz.z[i][j][k];
        dKij->xy[i][j][k] = dKij->xy[i][j][k] + Kij.xy[i][j][k]*dbx.x[i][j][k] + Kij.yy[i][j][k]*dby.x[i][j][k] + Kij.yz[i][j][k]*dbz.x[i][j][k];
        dKij->xz[i][j][k] = dKij->xz[i][j][k] + Kij.xz[i][j][k]*dbx.x[i][j][k] + Kij.yz[i][j][k]*dby.x[i][j][k] + Kij.zz[i][j][k]*dbz.x[i][j][k];
        dKij->yz[i][j][k] = dKij->yz[i][j][k] + Kij.xz[i][j][k]*dbx.y[i][j][k] + Kij.yz[i][j][k]*dby.y[i][j][k] + Kij.zz[i][j][k]*dbz.y[i][j][k];
}}}


}//End Shift Terms


free_dVector3D(d2alpha.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2alpha.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2alpha.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2alpha.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2alpha.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(d2alpha.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);

free_dVector3D(dalpha.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dalpha.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dalpha.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dby.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dby.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dby.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

free_dVector3D(dKxx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKxx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKxx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKyy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKyy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKyy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKzz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKzz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKzz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKxy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKxy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKxy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKxz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKxz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKxz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKyz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKyz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dKyz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dummy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
//free_dVector3D(K,0,Par.nxb,0,Par.nyb,0,Par.nzb);

return;
}
