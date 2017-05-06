#include "DeclareFunctions.h"

double dotproduct(Scalar A0,Vector A,Scalar B0,Vector B,Params Par){

 int i, j, k;
 double dp;

   dp = 0.0;
   for(i=0;i<Par.nxb;i++){
     for(j=0;j<Par.nyb;j++){
       for(k=0;k<Par.nzb;k++){
         dp += A0.s[i][j][k]*B0.s[i][j][k] + A.x[i][j][k]*B.x[i][j][k] 
             + A.y[i][j][k]*B.y[i][j][k] + A.z[i][j][k]*B.z[i][j][k];
   }}}

return dp;
}

