#include "DeclareFunctions.h"

//Output TwoTensor
int OutputStep_T(char outfile[],double time,Tensor T,Params Par){
 FILE *fp;
 int i,j,k;
 char filename[100],filename2[100];
 double error,gerr,x,y,z;

  sprintf(filename,"%s.t=%1.3f",outfile,time);

  fp = fopen(filename,"w");
  error = 0.;

// for(k=0;k<Par.nzb;k++){
// for(j=0;j<Par.nyb;j++){
 for(i=0;i<Par.nxb;i++){

   x = Par.x[i];
//   y = Par.x[j];
//   z = Par.x[k];
   fprintf(fp,"%lf %1.8e  %1.8e %1.8e %1.8e  %1.8e %1.8e\n",x,T.xx[i][12][12],T.yy[12][i][12],T.zz[12][12][i],T.xy[i][12][12],T.xz[i][12][12],T.yz[i][12][12]);
//}}}
}

 fclose(fp);

 return 0;
}
