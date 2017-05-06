#include "DeclareFunctions.h"
//Output Double

int OutputStep_S(char outfile[],double time,Scalar T,Params Par){
 FILE *fp;
// FILE *fp2;
 int i,j,k;
 char filename[100],filename2[100];
 double x,y,z;

  sprintf(filename, "%s.t=%d",outfile,(int)(100*(time+0.001)));
  fp = fopen(filename,"w");

  j = Par.nyb/2;
  k = Par.nzb/2;

 for(i=0;i<Par.nxb;i++){
   x = Par.x[i];
   fprintf(fp,"%lf %1.8e\n",x,T.s[i][j][k]);
}

 fclose(fp);

 return 0;
}
