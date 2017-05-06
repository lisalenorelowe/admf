#include "DeclareFunctions.h"

int ResetLapse(double t,Scalar alpha,Params Par){
int i,j,k;
double x,y,z,rinv,A;

//Reset Lapse 
if(Par.ilapse==30){  //Plane Wave (AWA)
 A = 0.02;
 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
   for(k=0;k<Par.nzb;k++){
    x = Par.x[i];
    alpha.s[i][j][k] = sqrt(1. + A*sin((x-t)));
}}}
}

return 0;
}
