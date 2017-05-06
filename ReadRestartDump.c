#include "DeclareFunctions.h"

int ReadRestartDump(Scalar *alpha,Vector *shift,Tensor *gij,Tensor *Kij,Params Par){
 FILE *fp;
 int i,j,k;
 char filename[100];
 double a,b,c,d,e,f;

//Allocate Structure Tensor Kij
gij->xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gij->yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gij->zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gij->xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gij->xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gij->yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Allocate Structure Tensor Kij
 Kij->xx     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Kij->yy     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Kij->zz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Kij->xy     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Kij->xz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Kij->yz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Allocate Structure Scalar alpha
 alpha->s   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Allocate Structure Vector shift 
 shift->x     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 shift->y     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 shift->z     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Read gij
sprintf(filename,"%s.ts=%d","Restartgij",Par.start);
fp = fopen(filename,"r");
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
 fscanf(fp,"%lf %lf %lf %lf %lf %lf\n",&gij->xx[i][j][k],&gij->yy[i][j][k],&gij->zz[i][j][k],&gij->xy[i][j][k],&gij->xz[i][j][k],&gij->yz[i][j][k]);
}}}
fclose(fp);

//Read Kij
sprintf(filename,"%s.ts=%d","RestartKij",Par.start);
fp = fopen(filename,"r");
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
 fscanf(fp,"%lf %lf %lf %lf %lf %lf\n",&Kij->xx[i][j][k],&Kij->yy[i][j][k],&Kij->zz[i][j][k],&Kij->xy[i][j][k],&Kij->xz[i][j][k],&Kij->yz[i][j][k]);}}}
fclose(fp);

//Read Lapse 
sprintf(filename,"%s.ts=%d","RestartLapse",Par.start);
fp = fopen(filename,"r");
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
 fscanf(fp,"%lf\n",&alpha->s[i][j][k]);}}}
fclose(fp);

//Read Shift 
sprintf(filename,"%s.ts=%d","RestartShift",Par.start);
fp = fopen(filename,"r");
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){
 fscanf(fp,"%lf %lf %lf\n",&shift->x[i][j][k],&shift->y[i][j][k],&shift->z[i][j][k]);}}}
fclose(fp);

 return 0;
}
