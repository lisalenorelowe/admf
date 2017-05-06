#include "DeclareFunctions.h"

void FunctionDerivs(double ***T,Tensor *dT,Params Par){
 int i,j,k,m;
 double x,y,z;
 double phi[XMAX], phibar[XMAX], phitilde[XMAX];
 double psibar[XMAX], psitilde[XMAX], psi2bar[XMAX], psi2tilde[XMAX];
 double dTdx[XMAX][YMAX][ZMAX],dTdy[XMAX][YMAX][ZMAX];
/*----------------------------------------------------------------*/

//The boxes have to be square
//X Derivs
for(j=0;j<Par.nyb;j++){
 for(k=0;k<Par.nzb;k++){
  for(i=0;i<Par.nxb;i++) phi[i]= T[i][j][k] ;
  fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
  fourierder(Par.xleftend,Par.xrightend,phibar,phitilde,psibar,psitilde,Par.nxb);
  fourierder(Par.xleftend,Par.xrightend,psibar,psitilde,psi2bar,psi2tilde,Par.nxb);
  for(i=0;i<Par.nxb;i++){
	dTdx[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,psibar,psitilde,Par.nxb,i,Par);
	dT->xx[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,psi2bar,psi2tilde,Par.nxb,i,Par);
  }
/*
  mychebft(Par.xleftend,Par.xrightend,fa,Par.nxb,f,Par);
  mychder(Par.xleftend,Par.xrightend,fa,fder,Par.nxb);
  mychder(Par.xleftend,Par.xrightend,fder,f2der,Par.nxb);
  for(i=0;i<Par.nxb;i++) {
   dTdx[i][j][k]= mychebev(Par.xleftend,Par.xrightend,fder,Par.nxb,Par.x[i]);
   dT->xx[i][j][k]= mychebev(Par.xleftend,Par.xrightend,f2der,Par.nxb,Par.x[i]);
  }
*/
 }
}

                                                                                
//Y Derivs
for(i=0;i<Par.nxb;i++){
 for(k=0;k<Par.nzb;k++){
  for(j=0;j<Par.nyb;j++) phi[j]= T[i][j][k] ;
  fourierfit(Par.yleftend,Par.yrightend,phibar,phitilde,Par.nyb,phi,Par);
  fourierder(Par.yleftend,Par.yrightend,phibar,phitilde,psibar,psitilde,Par.nyb);
  fourierder(Par.yleftend,Par.yrightend,psibar,psitilde,psi2bar,psi2tilde,Par.nyb);
  for(j=0;j<Par.nyb;j++){
	dTdy[i][j][k] = fouriereval(Par.yleftend,Par.yrightend,psibar,psitilde,Par.nyb,j,Par);
	dT->yy[i][j][k] = fouriereval(Par.yleftend,Par.yrightend,psi2bar,psi2tilde,Par.nyb,j,Par);
  }

/*
  mychebft(Par.yleftend,Par.yrightend,fa,Par.nyb,f,Par);
  mychder(Par.yleftend,Par.yrightend,fa,fder,Par.nyb);
  mychder(Par.yleftend,Par.yrightend,fder,f2der,Par.nyb);
  for(j=0;j<Par.nyb;j++) {
   dTdy[i][j][k]= mychebev(Par.yleftend,Par.yrightend,fder,Par.nyb,Par.y[j]);
   dT->yy[i][j][k]= mychebev(Par.yleftend,Par.yrightend,f2der,Par.nyb,Par.y[j]);
  }
*/
       
//Mixed Derivs, XY
  for(j=0;j<Par.nyb;j++) phi[j]= dTdx[i][j][k] ;
  fourierfit(Par.yleftend,Par.yrightend,phibar,phitilde,Par.nyb,phi,Par);
  fourierder(Par.yleftend,Par.yrightend,phibar,phitilde,psibar,psitilde,Par.nyb);
  for(j=0;j<Par.nyb;j++){
	dT->xy[i][j][k] = fouriereval(Par.yleftend,Par.yrightend,psibar,psitilde,Par.nyb,j,Par);
  }

}
}

/*
  mychebft(Par.yleftend,Par.yrightend,fa,Par.nyb,f,Par);
  mychder(Par.yleftend,Par.yrightend,fa,fder,Par.nyb);
  for(j=0;j<Par.nyb;j++) {
  dT->xy[i][j][k]= mychebev(Par.yleftend,Par.yrightend,fder,Par.nyb,Par.y[j]);
  }
  
 }
}
*/

//Z Derivs
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){

  for(k=0;k<Par.nzb;k++) phi[k]= T[i][j][k] ;
  fourierfit(Par.zleftend,Par.zrightend,phibar,phitilde,Par.nzb,phi,Par);
  fourierder(Par.zleftend,Par.zrightend,phibar,phitilde,psibar,psitilde,Par.nzb);
  fourierder(Par.zleftend,Par.zrightend,psibar,psitilde,psi2bar,psi2tilde,Par.nzb);
  //mychebft(Par.zleftend,Par.zrightend,fa,Par.nzb,f,Par);
  //mychder(Par.zleftend,Par.zrightend,fa,fder,Par.nzb);
  //mychder(Par.zleftend,Par.zrightend,fder,f2der,Par.nzb);
  for(k=0;k<Par.nzb;k++) {
   dT->zz[i][j][k]= fouriereval(Par.zleftend,Par.zrightend,psi2bar,psi2tilde,Par.nzb,k,Par);
   //dT->zz[i][j][k]= mychebev(Par.zleftend,Par.zrightend,f2der,Par.nzb,Par.z[k]);
  }

//Mixed Derivs, YZ
  for(k=0;k<Par.nzb;k++) phi[k]= dTdy[i][j][k] ;
  fourierfit(Par.zleftend,Par.zrightend,phibar,phitilde,Par.nzb,phi,Par);
  fourierder(Par.zleftend,Par.zrightend,phibar,phitilde,psibar,psitilde,Par.nzb);
  for(k=0;k<Par.nzb;k++) {
   dT->yz[i][j][k]= fouriereval(Par.zleftend,Par.zrightend,psibar,psitilde,Par.nzb,k,Par);
  }

//Mixed Derivs, XZ
  for(k=0;k<Par.nzb;k++) phi[k]= dTdx[i][j][k] ;
  fourierfit(Par.zleftend,Par.zrightend,phibar,phitilde,Par.nzb,phi,Par);
  fourierder(Par.zleftend,Par.zrightend,phibar,phitilde,psibar,psitilde,Par.nzb);
  for(k=0;k<Par.nzb;k++) {
   dT->xz[i][j][k]= fouriereval(Par.zleftend,Par.zrightend,psibar,psitilde,Par.nzb,k,Par);
  }


 }
}


return;
}
