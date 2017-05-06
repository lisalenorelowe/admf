#include "DeclareFunctions.h"

int InitializeShift(Vector *shift,Params Par){
int i,j,k;
double x,y,z,r2,A,sig,sig2,sig4,betar,grr,alp2;
 double cx2,cy2,cz2,cx4,cy4,cz4;

 shift->x     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 shift->y     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 shift->z     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Initialize Shift 

if(Par.ishift==0){ //Zero shift (\beta^i = 0)
 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
   for(k=0;k<Par.nzb;k++){
    shift->x[i][j][k] = 0. ;
    shift->y[i][j][k] = 0. ;
    shift->z[i][j][k] = 0. ;
 }}}
}

if(Par.ishift==30){ //Plane Wave
 printf("USE ISHIFT = 0 INSTEAD \n");
 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
   for(k=0;k<Par.nzb;k++){
    shift->x[i][j][k] = 0. ;
    shift->y[i][j][k] = 0. ;
    shift->z[i][j][k] = 0. ;
 }}}
}


 if(Par.ishift == 8){  // Cosine bump
   for(i=0;i<Par.nxb;i++){
     for(j=0;j<Par.nyb;j++){
       for(k=0;k<Par.nzb;k++){
	 x = Par.x[i];
	 y = Par.y[j];
	 z = Par.z[k];
	 A = 0.1;
	 alp2 = 1/(1 - A*A*(cos(x)*cos(x)*cos(y)*cos(y)*sin(z)*sin(z) 
				     + cos(x)*cos(x)*sin(y)*sin(y)*cos(z)*cos(z) 
				     + sin(x)*sin(x)*cos(y)*cos(y)*cos(z)*cos(z)));
	 shift->x[i][j][k] = alp2*A*sin(x)*cos(y)*cos(z);
	 shift->y[i][j][k] = alp2*A*cos(x)*sin(y)*cos(z);
	 shift->z[i][j][k] = alp2*A*cos(x)*cos(y)*sin(z);
       }}}
 }


 if(Par.ishift == 81){  // Cosine squared bump
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
	 alp2 = 1/(1.0 + 3.0*A*A*cx4*cy4*cz4 
				    - A*A*cx4*cy4*cz2 - A*A*cx4*cy2*cz4
				    - A*A*cx2*cy4*cz4);
	 shift->x[i][j][k] = alp2*A*cos(x/2.0)*cy2*cz2*sin(x/2.0);
	 shift->y[i][j][k] = alp2*A*cos(y/2.0)*cx2*cz2*sin(y/2.0);
	 shift->z[i][j][k] = alp2*A*cos(z/2.0)*cx2*cy2*sin(z/2.0);
       }}}
 }

if(Par.ishift==90){ //Flat trK distortion
 printf("USE ISHIFT = 0 INSTEAD \n");
 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++){ 
    shift->x[i][j][k] = 0. ;
    shift->y[i][j][k] = 0. ;
    shift->z[i][j][k] = 0. ;
 }}}
}


//End Initialize

return 0;
}
