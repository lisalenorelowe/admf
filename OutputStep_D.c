#include "DeclareFunctions.h"
//Output Double

int OutputStep_D(char outfile[],double time,double T[XMAX][YMAX][ZMAX],Params Par){
 FILE *fp;
// FILE *fp2;
 int i,j,k;
 char filename[100],filename2[100];
 double x,y,z;

//  sprintf(filename,"%s.t=%1.3f",outfile,time);
  sprintf(filename, "%s.t=%d",outfile,(int)(100*(time+0.001)));
  fp = fopen(filename,"w");
//  fp2 = fopen("HCfile","a");

//fprintf(fp2,"%f, %2.10f \n",time,T[Par.nxb/2][Par.nyb/2][Par.nzb/2]);
//fclose(fp2);

//  i = 12;
  j = Par.nyb/2;
  k = Par.nzb/2;

 for(i=0;i<Par.nxb;i++){
// for(j=0;j<Par.nyb;j++){
// for(k=0;k<Par.nzb;k++){

   x = Par.x[i];
//   y = Par.y[j];
//   z = Par.z[k];

   fprintf(fp,"%lf %1.8e\n",x,T[i][j][k]);
//}}}
}

 fclose(fp);

 return 0;
}
