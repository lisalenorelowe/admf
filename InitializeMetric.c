#include "DeclareFunctions.h"

int InitializeMetric(Tensor *gij,Params Par){

int i,j,k;
double x,y,z,rinv,A,r2,sig,sig2,sig4,grr,Krr,Ktt,eps,Aeps2;
 double cx2,cy2,cz2,cx3,cy3,cz3,cx4,cy4,cz4;

gij->xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gij->yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gij->zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gij->xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gij->xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
gij->yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Initialize gij 
 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
   for(k=0;k<Par.nzb;k++){

    x = Par.x[i];
    y = Par.y[j];
    z = Par.z[k];

  rinv = sqrt(Par.x[i]*Par.x[i]+Par.y[j]*Par.y[j]+Par.z[k]*Par.z[k]);
  rinv = 1./rinv;

if(Par.igij==0){  //Flat Metric
gij->xx[i][j][k] = 1.; // + 0.01*cos(x);
gij->yy[i][j][k] = 1.; // + 0.001*cos(y);
gij->zz[i][j][k] = 1.; // + 0.001*cos(z);
gij->xy[i][j][k] = 0.;
gij->xz[i][j][k] = 0.;
gij->yz[i][j][k] = 0.;
}

if(Par.igij==1){  //Flat Metric + random perturbations
gij->xx[i][j][k] = 1.+ 1.e-9*(2.*rand()/(RAND_MAX+1.)-1.);
gij->yy[i][j][k] = 1.+ 1.e-9*(2.*rand()/(RAND_MAX+1.)-1.);
gij->zz[i][j][k] = 1.+ 1.e-9*(2.*rand()/(RAND_MAX+1.)-1.);
gij->xy[i][j][k] = 1.e-9*(2.*rand()/(RAND_MAX+1.)-1.);
gij->xz[i][j][k] = 1.e-9*(2.*rand()/(RAND_MAX+1.)-1.);
gij->yz[i][j][k] = 1.e-9*(2.*rand()/(RAND_MAX+1.)-1.);
}

if(Par.igij==30){  //Plane Wave (AWA)
A = 0.02;
gij->xx[i][j][k] = 1. + A*sin(x);
gij->yy[i][j][k] = 1.;
gij->zz[i][j][k] = 1.;
gij->xy[i][j][k] = 0.;
gij->xz[i][j][k] = 0.;
gij->yz[i][j][k] = 0.;
}

 if(Par.igij==8){  //Cosine bump
   A = 0.1;
//   eps = -0.0001;
   eps = 0.0;
   Aeps2 = (A + eps)*(A + eps);
   gij->xx[i][j][k] = 1.0 - Aeps2*cos(y)*cos(y)*cos(z)*cos(z)*sin(x)*sin(x);
   gij->yy[i][j][k] = 1.0 - Aeps2*cos(x)*cos(x)*cos(z)*cos(z)*sin(y)*sin(y);
   gij->zz[i][j][k] = 1.0 - Aeps2*cos(x)*cos(x)*cos(y)*cos(y)*sin(z)*sin(z);
   gij->xy[i][j][k] = -Aeps2*sin(x)*cos(x)*sin(y)*cos(y)*cos(z)*cos(z);
   gij->xz[i][j][k] = -Aeps2*sin(x)*cos(x)*sin(z)*cos(z)*cos(y)*cos(y);
   gij->yz[i][j][k] = -Aeps2*sin(y)*cos(y)*sin(z)*cos(z)*cos(x)*cos(x);
 }
   
 if(Par.igij==81){  //Cosine squared bump
   A = 0.2;
   //eps = 0.001;
   eps = 0.0;
   Aeps2 = (A + eps)*(A + eps);
   cx2 = cos(x/2.0)*cos(x/2.0);
   cy2 = cos(y/2.0)*cos(y/2.0);
   cz2 = cos(z/2.0)*cos(z/2.0);
   cx3 = cx2*cos(x/2.0);
   cy3 = cy2*cos(y/2.0);
   cz3 = cz2*cos(z/2.0);
   cx4 = cx2*cx2;
   cy4 = cy2*cy2;
   cz4 = cz2*cz2;
   gij->xx[i][j][k] = -Aeps2*cx2*cy4*cz4 + Aeps2*cx4*cy4*cz4 + 1.0;
   gij->yy[i][j][k] = -Aeps2*cx4*cy2*cz4 + Aeps2*cx4*cy4*cz4 + 1.0;
   gij->zz[i][j][k] = -Aeps2*cx4*cy4*cz2 + Aeps2*cx4*cy4*cz4 + 1.0;
   gij->xy[i][j][k] = -Aeps2*cx3*cy3*cz4*sin(x/2.0)*sin(y/2.0);
   gij->xz[i][j][k] = -Aeps2*cx3*cy4*cz3*sin(x/2.0)*sin(z/2.0);
   gij->yz[i][j][k] = -Aeps2*cx4*cy3*cz3*sin(y/2.0)*sin(z/2.0);
 }

 if(Par.igij==90){  //Flat trK distortion
   gij->xx[i][j][k] = 1.0;
   gij->yy[i][j][k] = 1.0;
   gij->zz[i][j][k] = 1.0;
   gij->xy[i][j][k] = 0.0;
   gij->xz[i][j][k] = 0.0;
   gij->yz[i][j][k] = 0.0;
 }


}}} //End ijk loop


//End Initialize

return 0;
}
