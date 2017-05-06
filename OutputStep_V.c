#include "DeclareFunctions.h"

//Output Tensor (V for Vector)
int OutputStep_V(char outfile[],double time,Vector T,Params Par){
 FILE *fp;
 int i,j,k;
 char filename[100],filename2[100];
 double x,y,z;

  sprintf(filename,"%s.t=%d",outfile,(int)(100*(time+0.001)));
  fp = fopen(filename,"w");

   j = Par.nyb/2;
   k = Par.nzb/2; 

 for(i=0;i<Par.nxb;i++){
   x = Par.x[i];
   fprintf(fp,"%lf %1.8e  %1.8e %1.8e\n",Par.x[i],T.x[i][j][k],T.y[i][j][k],T.z[i][j][k]);
}


 fclose(fp);

 return 0;
}
