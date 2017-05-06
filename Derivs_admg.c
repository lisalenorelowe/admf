#include "DeclareFunctions.h"

void Derivs_admg(double t,Scalar alpha,Vector shift, Tensor Kij,Tensor gij,Tensor *dgij,Params Par){
 int i,j,k,m;
 Vector dbx,dby,dbz,dgxx,dgyy,dgzz,dgxy,dgxz,dgyz;
/*----------------------------------------------------------------*/


if(Par.ishift!=0){   
dbx.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbx.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbx.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

dby.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dby.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dby.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

dbz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

dgxx.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgxx.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgxx.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

dgyy.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgyy.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgyy.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

dgzz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgzz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgzz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

dgxy.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgxy.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgxy.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

dgyz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgyz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgyz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

dgxz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgxz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dgxz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
}



                     // Get Derivatives of Shift Vector \beta_i
if(Par.ishift!=0){    // But not if shift is zero
 FunctionFirstDerivs(shift.x,&dbx,Par);
 FunctionFirstDerivs(shift.y,&dby,Par);
 FunctionFirstDerivs(shift.z,&dbz,Par);
}

		   // Get Derivatives of Metric
if(Par.ishift!=0){    // But not if shift is zero
FunctionFirstDerivs(gij.xx,&dgxx,Par);
FunctionFirstDerivs(gij.yy,&dgyy,Par);
FunctionFirstDerivs(gij.zz,&dgzz,Par);
FunctionFirstDerivs(gij.xy,&dgxy,Par);
FunctionFirstDerivs(gij.xz,&dgxz,Par);
FunctionFirstDerivs(gij.yz,&dgyz,Par);
}

// RHS Gdot Equations
//$\partial g{ij}= -2\alpha K_{ij} + \beta^k\partial_kg_{ij} + g_{ik}\partial_j\beta^k + g_{kj}\partial_i\beta^k$

for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
//  First Part: $ -2\alpha K_{ij} $
	dgij->xx[i][j][k] = - 2.*alpha.s[i][j][k]*Kij.xx[i][j][k];
        dgij->yy[i][j][k] = - 2.*alpha.s[i][j][k]*Kij.yy[i][j][k]; 
        dgij->zz[i][j][k] = - 2.*alpha.s[i][j][k]*Kij.zz[i][j][k];
        dgij->xy[i][j][k] = - 2.*alpha.s[i][j][k]*Kij.xy[i][j][k];
        dgij->xz[i][j][k] = - 2.*alpha.s[i][j][k]*Kij.xz[i][j][k];
        dgij->yz[i][j][k] = - 2.*alpha.s[i][j][k]*Kij.yz[i][j][k];
}}}

		     //Put in Shift Terms of equation
if(Par.ishift!=0){    // But not if shift is zero
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
// Shift Part: $\beta^k \partial_k g_{ij} + g_{ik} \partial_j \beta^k + g_{kj} \partial_i \beta^k$

     // \beta^k $\partial_k g_{ij}$
        dgij->xx[i][j][k] = dgij->xx[i][j][k] + shift.x[i][j][k]*dgxx.x[i][j][k] + shift.y[i][j][k]*dgxx.y[i][j][k] + shift.z[i][j][k]*dgxx.z[i][j][k];
        dgij->yy[i][j][k] = dgij->yy[i][j][k] + shift.x[i][j][k]*dgyy.x[i][j][k] + shift.y[i][j][k]*dgyy.y[i][j][k] + shift.z[i][j][k]*dgyy.z[i][j][k];
        dgij->zz[i][j][k] = dgij->zz[i][j][k] + shift.x[i][j][k]*dgzz.x[i][j][k] + shift.y[i][j][k]*dgzz.y[i][j][k] + shift.z[i][j][k]*dgzz.z[i][j][k]; 
        dgij->xy[i][j][k] = dgij->xy[i][j][k] + shift.x[i][j][k]*dgxy.x[i][j][k] + shift.y[i][j][k]*dgxy.y[i][j][k] + shift.z[i][j][k]*dgxy.z[i][j][k];
        dgij->xz[i][j][k] = dgij->xz[i][j][k] + shift.x[i][j][k]*dgxz.x[i][j][k] + shift.y[i][j][k]*dgxz.y[i][j][k] + shift.z[i][j][k]*dgxz.z[i][j][k];
        dgij->yz[i][j][k] = dgij->yz[i][j][k] + shift.x[i][j][k]*dgyz.x[i][j][k] + shift.y[i][j][k]*dgyz.y[i][j][k] + shift.z[i][j][k]*dgyz.z[i][j][k];

     // g_{ik} \partial_j \beta^k
        dgij->xx[i][j][k] = dgij->xx[i][j][k] + gij.xx[i][j][k]*dbx.x[i][j][k] + gij.xy[i][j][k]*dby.x[i][j][k] + gij.xz[i][j][k]*dbz.x[i][j][k];
        dgij->yy[i][j][k] = dgij->yy[i][j][k] + gij.xy[i][j][k]*dbx.y[i][j][k] + gij.yy[i][j][k]*dby.y[i][j][k] + gij.yz[i][j][k]*dbz.y[i][j][k];
        dgij->zz[i][j][k] = dgij->zz[i][j][k] + gij.xz[i][j][k]*dbx.z[i][j][k] + gij.yz[i][j][k]*dby.z[i][j][k] + gij.zz[i][j][k]*dbz.z[i][j][k];
        dgij->xy[i][j][k] = dgij->xy[i][j][k] + gij.xx[i][j][k]*dbx.y[i][j][k] + gij.xy[i][j][k]*dby.y[i][j][k] + gij.xz[i][j][k]*dbz.y[i][j][k];
        dgij->xz[i][j][k] = dgij->xz[i][j][k] + gij.xx[i][j][k]*dbx.z[i][j][k] + gij.xy[i][j][k]*dby.z[i][j][k] + gij.xz[i][j][k]*dbz.z[i][j][k];
        dgij->yz[i][j][k] = dgij->yz[i][j][k] + gij.xy[i][j][k]*dbx.z[i][j][k] + gij.yy[i][j][k]*dby.z[i][j][k] + gij.yz[i][j][k]*dbz.z[i][j][k];

     // g_{kj} \partial_i \beta^k
        dgij->xx[i][j][k] = dgij->xx[i][j][k] + gij.xx[i][j][k]*dbx.x[i][j][k] + gij.xy[i][j][k]*dby.x[i][j][k] + gij.xz[i][j][k]*dbz.x[i][j][k];
        dgij->yy[i][j][k] = dgij->yy[i][j][k] + gij.xy[i][j][k]*dbx.y[i][j][k] + gij.yy[i][j][k]*dby.y[i][j][k] + gij.yz[i][j][k]*dbz.y[i][j][k];
        dgij->zz[i][j][k] = dgij->zz[i][j][k] + gij.xz[i][j][k]*dbx.z[i][j][k] + gij.yz[i][j][k]*dby.z[i][j][k] + gij.zz[i][j][k]*dbz.z[i][j][k];
        dgij->xy[i][j][k] = dgij->xy[i][j][k] + gij.xy[i][j][k]*dbx.x[i][j][k] + gij.yy[i][j][k]*dby.x[i][j][k] + gij.yz[i][j][k]*dbz.x[i][j][k];
        dgij->xz[i][j][k] = dgij->xz[i][j][k] + gij.xz[i][j][k]*dbx.x[i][j][k] + gij.yz[i][j][k]*dby.x[i][j][k] + gij.zz[i][j][k]*dbz.x[i][j][k];
        dgij->yz[i][j][k] = dgij->yz[i][j][k] + gij.xz[i][j][k]*dbx.y[i][j][k] + gij.yz[i][j][k]*dby.y[i][j][k] + gij.zz[i][j][k]*dbz.y[i][j][k];
}}}


}


//Boundary Conditions on Metric Derivatives
//bc_on_dT(t,gij,dgij,Par);

if(Par.ishift!=0){
free_dVector3D(dbx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dby.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dby.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dby.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

free_dVector3D(dgxx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgzz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgzz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgzz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgxz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dgyz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
}

return;
/*******************************************************/

}
