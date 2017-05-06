#include "DeclareFunctions.h"

void FunctionFirstDerivs(double ***T,Vector *dT,Params Par){
 int i,j,k,m;
 double x,y,z;
 double phi[XMAX], phibar[XMAX], phitilde[XMAX];
 double psibar[XMAX], psitilde[XMAX];
/*----------------------------------------------------------------*/
//The boxes have to be square
//X Derivs
for(j=0;j<Par.nyb;j++){
 for(k=0;k<Par.nzb;k++){
  for(i=0;i<Par.nxb;i++) phi[i]= T[i][j][k] ;
  fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
  fourierder(Par.xleftend,Par.xrightend,phibar,phitilde,psibar,psitilde,Par.nxb);
  for(i=0;i<Par.nxb;i++) {
   dT->x[i][j][k]= fouriereval(Par.xleftend,Par.xrightend,psibar,psitilde,Par.nxb,i,Par);
  }
 }
}
                                                                                 

//Y Derivs
for(i=0;i<Par.nxb;i++){
 for(k=0;k<Par.nzb;k++){
  for(j=0;j<Par.nyb;j++) phi[j]= T[i][j][k] ;
  fourierfit(Par.yleftend,Par.yrightend,phibar,phitilde,Par.nyb,phi,Par);
  fourierder(Par.yleftend,Par.yrightend,phibar,phitilde,psibar,psitilde,Par.nyb);
  for(j=0;j<Par.nyb;j++) {
   dT->y[i][j][k]= fouriereval(Par.yleftend,Par.yrightend,psibar,psitilde,Par.nyb,j,Par);
  }
 }
}


//Z Derivs
for(i=0;i<Par.nxb;i++){
 for(j=0;j<Par.nyb;j++){
  for(k=0;k<Par.nzb;k++) phi[k]= T[i][j][k] ;
  fourierfit(Par.zleftend,Par.zrightend,phibar,phitilde,Par.nzb,phi,Par);
  fourierder(Par.zleftend,Par.zrightend,phibar,phitilde,psibar,psitilde,Par.nzb);
  for(k=0;k<Par.nzb;k++) {
   dT->z[i][j][k]= fouriereval(Par.zleftend,Par.zrightend,psibar,psitilde,Par.nzb,k,Par);
  }
 }
}

return;
}
