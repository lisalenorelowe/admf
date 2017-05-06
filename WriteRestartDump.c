#include "DeclareFunctions.h"

int WriteRestartDump(int timestep,Scalar alpha,Vector shift,Tensor gij,Tensor Kij,Params Par){
 FILE *fp;
 int i,j,k;
 char filename[100];

//Dump gij
sprintf(filename,"%s.ts=%d","Restartgij",timestep);
fp = fopen(filename,"w");
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
 fprintf(fp,"%1.18e %1.18e %1.18e %1.18e %1.18e %1.18e\n",gij.xx[i][j][k],gij.yy[i][j][k],gij.zz[i][j][k],gij.xy[i][j][k],gij.xz[i][j][k],gij.yz[i][j][k]);
}}}
fclose(fp);

//Dump Kij
sprintf(filename,"%s.ts=%d","RestartKij",timestep);
fp = fopen(filename,"w");
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
 fprintf(fp,"%1.18e %1.18e %1.18e %1.18e %1.18e %1.18e\n",Kij.xx[i][j][k],Kij.yy[i][j][k],Kij.zz[i][j][k],Kij.xy[i][j][k],Kij.xz[i][j][k],Kij.yz[i][j][k]);}}}
fclose(fp);

//Dump Lapse 
sprintf(filename,"%s.ts=%d","RestartLapse",timestep);
fp = fopen(filename,"w");
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
 fprintf(fp,"%1.18e\n",alpha.s[i][j][k]);}}}
fclose(fp);

//Dump Shift 
sprintf(filename,"%s.ts=%d","RestartShift",timestep);
fp = fopen(filename,"w");
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
 fprintf(fp,"%1.18e %1.18e %1.18e\n",shift.x[i][j][k],shift.y[i][j][k],shift.z[i][j][k]);}}}
fclose(fp);

 return 0;
}
