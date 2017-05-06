#include "DeclareFunctions.h"
//find the inverse of metric, called gup

int GetDetg(Tensor gij,Scalar detg,Params Par)
//input- metric, output- gup
{
  int i,j,k;

//find the metric determinant g
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
   detg.s[i][j][k] = gij.xx[i][j][k]*gij.yy[i][j][k]*gij.zz[i][j][k] +
    2.*gij.xy[i][j][k]*gij.xz[i][j][k]*gij.yz[i][j][k] -
    gij.xx[i][j][k]*gij.yz[i][j][k]*gij.yz[i][j][k] - 
    gij.yy[i][j][k]*gij.xz[i][j][k]*gij.xz[i][j][k] -
    gij.zz[i][j][k]*gij.xy[i][j][k]*gij.xy[i][j][k];
}}}

return 0;
}
