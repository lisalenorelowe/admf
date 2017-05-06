#include "DeclareFunctions.h"

int InitializeLapse(Scalar *alpha,Params Par){
int i,j,k;
double x,y,z,rinv,A,r2,grr,sig,sig2,sig4;
 double cx2,cy2,cz2,cx4,cy4,cz4;

alpha->s   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Initialize Lapse 

if(Par.ilapse==0){ //Lapse  \alpha = 1
 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
   for(k=0;k<Par.nzb;k++){
    alpha->s[i][j][k]=1.;
//     alpha->s[i][j][k]=0.0;
 }}}
}

if(Par.ilapse==30){  //Plane Wave (AWA)
                    //Parameter A is also defined in "ResetLapse.c"
 A = 0.02;
 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
   for(k=0;k<Par.nzb;k++){
    x = Par.x[i];
    y = Par.y[j];
    z = Par.z[k];
    alpha->s[i][j][k] = sqrt(1. + A*sin(x));
}}}
}

 if(Par.ilapse == 8){  //Cosine bump
   for(i=0;i<Par.nxb;i++){
     for(j=0;j<Par.nyb;j++){
       for(k=0;k<Par.nzb;k++){
	 x = Par.x[i];
	 y = Par.y[j];
	 z = Par.z[k];
	 A = 0.1;
	 alpha->s[i][j][k] = 1/sqrt(1 - A*A*(cos(x)*cos(x)*cos(y)*cos(y)*sin(z)*sin(z) 
				     + cos(x)*cos(x)*sin(y)*sin(y)*cos(z)*cos(z) 
				     + sin(x)*sin(x)*cos(y)*cos(y)*cos(z)*cos(z)));
       }}}
 }

 if(Par.ilapse == 81){  //Cosine squared bump
   for(i=0;i<Par.nxb;i++){
     for(j=0;j<Par.nyb;j++){
       for(k=0;k<Par.nzb;k++){
	 x = Par.x[i];
	 y = Par.y[j];
	 z = Par.z[k];
	 A = 0.2;
	 cx2 = cos(x/2.0)*cos(x/2.0);
	 cy2 = cos(y/2.0)*cos(y/2.0);
	 cz2 = cos(z/2.0)*cos(z/2.0);
	 cx4 = cx2*cx2;
	 cy4 = cy2*cy2;
	 cz4 = cz2*cz2;

	 alpha->s[i][j][k] = 1/sqrt(1.0 + 3.0*A*A*cx4*cy4*cz4 
				    - A*A*cx4*cy4*cz2 - A*A*cx4*cy2*cz4
				    - A*A*cx2*cy4*cz4);
       }}}
 }


if(Par.ilapse==90){ //Flat trK distortion
 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
   for(k=0;k<Par.nzb;k++){
    alpha->s[i][j][k]=1.;
 }}}
}


return 0;
}
