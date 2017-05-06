#include "DeclareFunctions.h"
//find the inverse of metric, called gup

int ginv(Tensor gij,Tensor *gup,Params Par)
//input- metric, output- gup
{
  double temp[XMAX][YMAX][ZMAX];
  int i,j,k;

//find the metric determinant g, then 1/g
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
   temp[i][j][k] = gij.xx[i][j][k]*gij.yy[i][j][k]*gij.zz[i][j][k] +
    2.*gij.xy[i][j][k]*gij.xz[i][j][k]*gij.yz[i][j][k] -
    gij.xx[i][j][k]*gij.yz[i][j][k]*gij.yz[i][j][k] - 
    gij.yy[i][j][k]*gij.xz[i][j][k]*gij.xz[i][j][k] -
    gij.zz[i][j][k]*gij.xy[i][j][k]*gij.xy[i][j][k];
    if(temp[i][j][k]>0.) temp[i][j][k] = 1./temp[i][j][k];
    if(temp[i][j][k]<=0.) printf("METRIC DETERMINANT IS ZERO or NEGATIVE!  %1.8e\n",temp);
    if(temp[i][j][k]<=0.) exit(1);
}}}

//g^{xx}
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
   gup->xx[i][j][k] = temp[i][j][k]*(
    gij.yy[i][j][k]*gij.zz[i][j][k] -
    gij.yz[i][j][k]*gij.yz[i][j][k]);
}}}

//g^{yy}
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
   gup->yy[i][j][k] = temp[i][j][k]*(
    gij.xx[i][j][k]*gij.zz[i][j][k] -
    gij.xz[i][j][k]*gij.xz[i][j][k]);
}}}

//g^{zz}
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
   gup->zz[i][j][k] = temp[i][j][k]*(
    gij.xx[i][j][k]*gij.yy[i][j][k] -
    gij.xy[i][j][k]*gij.xy[i][j][k]);
}}}

//g^{xy}
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
   gup->xy[i][j][k] = -temp[i][j][k]*(
    gij.zz[i][j][k]*gij.xy[i][j][k] -
    gij.xz[i][j][k]*gij.yz[i][j][k]);
}}}

//g^{xz}
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
   gup->xz[i][j][k] = temp[i][j][k]*(
    gij.xy[i][j][k]*gij.yz[i][j][k] -
    gij.yy[i][j][k]*gij.xz[i][j][k]);
}}}
  
//g^{yz}
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
   gup->yz[i][j][k] = -temp[i][j][k]*(
    gij.xx[i][j][k]*gij.yz[i][j][k] -
    gij.xy[i][j][k]*gij.xz[i][j][k]);
}}}
     

return 0;
}
