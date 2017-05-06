#include "DeclareFunctions.h"

int InitializeKij(Tensor *Kij,Tensor *gij,Params Par){

int i,j,k;
double x,y,z,rinv,A,r2,sig,sig2,sig4,grr,Krr,Ktt,eps,alp;
 double cx2,cy2,cz2,cx4,cy4,cz4,f1,f2,f3;

//Allocate Structure Tensor Kij 
 Kij->xx     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Kij->yy     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Kij->zz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Kij->xy     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Kij->xz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Kij->yz     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Initialize Kij 
 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
   for(k=0;k<Par.nzb;k++){

    x = Par.x[i];
    y = Par.y[j];
    z = Par.z[k];

   rinv = sqrt(Par.x[i]*Par.x[i]+Par.y[j]*Par.y[j]+Par.z[k]*Par.z[k]);
   rinv = 1./rinv;

if(Par.iKij==0){ //Kij=0 'moment of time symmetry'
    Kij->xx[i][j][k] = 0.;
    Kij->yy[i][j][k] = 0.;
    Kij->zz[i][j][k] = 0.;
    Kij->xy[i][j][k] = 0.;
    Kij->xz[i][j][k] = 0.;
    Kij->yz[i][j][k] = 0.;
}

if(Par.iKij==1){ //Kij=0 with random perturbations
    Kij->xx[i][j][k] = 1.e-9*rand()/(RAND_MAX+1.);
    Kij->yy[i][j][k] = 1.e-9*rand()/(RAND_MAX+1.);
    Kij->zz[i][j][k] = 1.e-9*rand()/(RAND_MAX+1.);
    Kij->xy[i][j][k] = 1.e-9*rand()/(RAND_MAX+1.);
    Kij->xz[i][j][k] = 1.e-9*rand()/(RAND_MAX+1.);
    Kij->yz[i][j][k] = 1.e-9*rand()/(RAND_MAX+1.);
}

if(Par.iKij==30){  //Plane Wave (AWA)
    A = 0.02;
    eps = 0.01;
    Kij->xx[i][j][k] = 0.5*A*cos(x)/sqrt(1. + A*sin(x)) + eps*gij->xx[i][j][k];
    Kij->yy[i][j][k] = eps*gij->yy[i][j][k];
    Kij->zz[i][j][k] = eps*gij->zz[i][j][k];
    Kij->xy[i][j][k] = eps*gij->xy[i][j][k];
    Kij->xz[i][j][k] = eps*gij->xz[i][j][k];
    Kij->yz[i][j][k] = eps*gij->yz[i][j][k];
}

 if(Par.iKij == 8){ //Cosine bump
   A = 0.1;
   eps = 0.01;
   alp = 1/sqrt(1 - A*A*(cos(x)*cos(x)*cos(y)*cos(y)*sin(z)*sin(z) 
				     + cos(x)*cos(x)*sin(y)*sin(y)*cos(z)*cos(z) 
				     + sin(x)*sin(x)*cos(y)*cos(y)*cos(z)*cos(z)));
   Kij->xx[i][j][k] = alp*A*cos(x)*cos(y)*cos(z) + eps*gij->xx[i][j][k];
   Kij->yy[i][j][k] = alp*A*cos(x)*cos(y)*cos(z) + eps*gij->yy[i][j][k];
   Kij->zz[i][j][k] = alp*A*cos(x)*cos(y)*cos(z) + eps*gij->zz[i][j][k];
   Kij->xy[i][j][k] = -alp*A*sin(x)*sin(y)*cos(z) + eps*gij->xy[i][j][k];
   Kij->xz[i][j][k] = -alp*A*sin(x)*cos(y)*sin(z) + eps*gij->xz[i][j][k];
   Kij->yz[i][j][k] = -alp*A*cos(x)*sin(y)*sin(z) + eps*gij->yz[i][j][k];
 }


 if(Par.iKij == 81){ //Cosine squared bump
   A = 0.2;
   cx2 = cos(x/2.0)*cos(x/2.0);
   cy2 = cos(y/2.0)*cos(y/2.0);
   cz2 = cos(z/2.0)*cos(z/2.0);
   cx4 = cx2*cx2;
   cy4 = cy2*cy2;
   cz4 = cz2*cz2;
   alp = 1/sqrt(1.0 + 3.0*A*A*cx4*cy4*cz4 
				    - A*A*cx4*cy4*cz2 - A*A*cx4*cy2*cz4
				    - A*A*cx2*cy4*cz4);
   Kij->xx[i][j][k] = 0.5*A*cy2*cz2*(-1.0 + 2.0*cx2)*alp;
   Kij->yy[i][j][k] = 0.5*A*cx2*cz2*(-1.0 + 2.0*cy2)*alp;
   Kij->zz[i][j][k] = 0.5*A*cx2*cy2*(-1.0 + 2.0*cz2)*alp;
   Kij->xy[i][j][k] = -A*cos(x/2.0)*cos(y/2.0)*cz2*sin(x/2.0)*sin(y/2.0)*alp;
   Kij->xz[i][j][k] = -A*cos(x/2.0)*cos(z/2.0)*cy2*sin(x/2.0)*sin(z/2.0)*alp;
   Kij->yz[i][j][k] = -A*cos(y/2.0)*cos(z/2.0)*cx2*sin(y/2.0)*sin(z/2.0)*alp;
 }


 if(Par.iKij == 90){ //Flat trK distortion
   A = 0.01;
//   sig = 0.25;
//   f1 = A*exp(-(y*y*y*y + z*z*z*z)/(sig*sig)); //function of y and z
//   f2 = A*exp(-(x*x*x*x + z*z*z*z)/(sig*sig)); //function of x and z
//   f3 = A*exp(-(x*x*x*x + y*y*y*y)/(sig*sig)); //function of x and y
   f1 = A*(cos(y)*cos(y) + cos(z)*cos(z));
   f2 = A*(cos(x)*cos(x) + cos(z)*cos(z));
   f3 = A*(cos(x)*cos(x) + cos(y)*cos(y));
   Kij->xx[i][j][k] = -f1 + f2 + f3;
   Kij->yy[i][j][k] = f1 - f2 + f3;
   Kij->zz[i][j][k] = f1 + f2 - f3;
   Kij->xy[i][j][k] = 0.0;
   Kij->xz[i][j][k] = 0.0;
   Kij->yz[i][j][k] = 0.0;
 }

 if(Par.iKij == 95){ //DiagK
   A = 0.01*cos(x);
   Kij->xx[i][j][k] = (A/3.0)*gij->xx[i][j][k];
   Kij->yy[i][j][k] = (A/3.0)*gij->yy[i][j][k];
   Kij->zz[i][j][k] = (A/3.0)*gij->zz[i][j][k];
   Kij->xy[i][j][k] = (A/3.0)*gij->xy[i][j][k];
   Kij->xz[i][j][k] = (A/3.0)*gij->xz[i][j][k];
   Kij->yz[i][j][k] = (A/3.0)*gij->yz[i][j][k];
 }


}}}//end i,j,k loop

//End Initialize

return 0;
}
