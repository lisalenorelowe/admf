#include "DeclareFunctions.h"

int OutputStep_gij(char outfile[],double time,Tensor gij,Tensor Kij,Params Par){
 FILE *fp;
 int i,j,k;
 char filename[100],filename2[100];
 double x,y,z;
 double A, gxxexact;

  sprintf(filename,"%s.t=%d",outfile,(int)(100*(time+0.001)));

  fp = fopen(filename,"w");

 j = Par.nyb/2;
 k = Par.nzb/2;

 for(i=0;i<Par.nxb;i++){
//  for(j=0;j<Par.nyb;j++){
//   for(k=0;k<Par.nzb;k++){

 x = Par.x[i];
//   y = Par.y[j];
//   z = Par.z[k];

//For the cosine bump (initial data case 8)
 A = 0.1;
 gxxexact = 1.0 - A*A*cos(Par.z[k])*cos(Par.z[k])*cos(Par.y[j])*cos(Par.y[j])*sin(Par.x[i])*sin(Par.x[i]) ;
 fprintf(fp,"%lf %1.8e %1.8e %1.8e\n",x,gij.xx[i][j][k],Kij.xx[i][j][k],gij.xx[i][j][k] - gxxexact);

//}}}
}

 fclose(fp);

 return 0;
}
