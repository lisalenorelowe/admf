#include "DeclareFunctions.h"

void Derivs_admlapse(double t,Scalar alpha,Vector shift,Tensor gup,Tensor Kij,double ***dalpha,Params Par){
 int i,j,k,m;
 double K[XMAX][YMAX][ZMAX],alpha2;
/*----------------------------------------------------------------*/

if(Par.elapse==0){  //Do not evolve lapse
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
	dalpha[i][j][k] = 0.;
}}}
}

if(Par.elapse==1){  //Harmonic lapse evolution  alpha dot = -alpha^2 K
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
        K[i][j][k] = gup.xx[i][j][k]*Kij.xx[i][j][k] + gup.yy[i][j][k]*Kij.yy[i][j][k] +gup.zz[i][j][k]*Kij.zz[i][j][k] +
        2.*gup.xy[i][j][k]*Kij.xy[i][j][k] +2.*gup.xz[i][j][k]*Kij.xz[i][j][k] +2.*gup.yz[i][j][k]*Kij.yz[i][j][k] ;

        dalpha[i][j][k] = - alpha.s[i][j][k]*alpha.s[i][j][k]*K[i][j][k];

//        bc3_on_dT(t,alpha,dalpha,Par);

}}}
}


return;
/*******************************************************/

}
