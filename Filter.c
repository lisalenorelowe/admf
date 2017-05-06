#include "DeclareFunctions.h"

void Filter(Tensor gij,Tensor Kij,Params Par,int imod3)
{

double phi[XMAX], chi[XMAX], phibar[XMAX], phitilde[XMAX], chibar[XMAX], chitilde[XMAX];
int i, j, k;

//Par.nfilter = number of modes to remove. 
//modes removed are 0, 1, 2,...,Par.nfilter-1 and N, N-1, N-2,...,N-Par.nfilter+1.
//where N=nxb=nyb=nzb.
//
// Mode counting example with N=10 and Par.nfilter=2 (the "zero mode" is k=N/2=5): 
// k =     0  1  2  3  4  5  6  7  8  9 
// remove: X  X                       X

printf("calling Filter with imod3 = %d %d \n", imod3, Par.nfilter);

if (imod3 == 0){
	for(j=0;j<Par.nyb;j++){
	for(k=0;k<Par.nzb;k++){
//xx
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.xx[i][j][k];
			chi[i] = Kij.xx[i][j][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.xx[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.xx[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//yy		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.yy[i][j][k];
			chi[i] = Kij.yy[i][j][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.yy[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.yy[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//zz		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.zz[i][j][k];
			chi[i] = Kij.zz[i][j][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.zz[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.zz[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//xy		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.xy[i][j][k];
			chi[i] = Kij.xy[i][j][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.xy[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.xy[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//xz		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.xz[i][j][k];
			chi[i] = Kij.xz[i][j][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.xz[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.xz[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//yz		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.yz[i][j][k];
			chi[i] = Kij.yz[i][j][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.yz[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.yz[i][j][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
	
	}}
}


if (imod3 == 1){
	for(j=0;j<Par.nyb;j++){
	for(k=0;k<Par.nzb;k++){
//xx
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.xx[j][i][k];
			chi[i] = Kij.xx[j][i][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.xx[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.xx[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//yy		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.yy[j][i][k];
			chi[i] = Kij.yy[j][i][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.yy[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.yy[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//zz		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.zz[j][i][k];
			chi[i] = Kij.zz[j][i][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.zz[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.zz[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//xy		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.xy[j][i][k];
			chi[i] = Kij.xy[j][i][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.xy[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.xy[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//xz		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.xz[j][i][k];
			chi[i] = Kij.xz[j][i][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.xz[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.xz[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//yz		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.yz[j][i][k];
			chi[i] = Kij.yz[j][i][k];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.yz[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.yz[j][i][k] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
	
	}}
}

if (imod3 == 2){
	for(j=0;j<Par.nyb;j++){
	for(k=0;k<Par.nzb;k++){
//xx
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.xx[j][k][i];
			chi[i] = Kij.xx[j][k][i];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.xx[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.xx[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//yy		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.yy[j][k][i];
			chi[i] = Kij.yy[j][k][i];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.yy[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.yy[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//zz		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.zz[j][k][i];
			chi[i] = Kij.zz[j][k][i];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.zz[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.zz[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//xy		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.xy[j][k][i];
			chi[i] = Kij.xy[j][k][i];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.xy[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.xy[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//xz		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.xz[j][k][i];
			chi[i] = Kij.xz[j][k][i];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.xz[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.xz[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
//yz		
		for(i=0;i<Par.nxb;i++){ 
			phi[i] = gij.yz[j][k][i];
			chi[i] = Kij.yz[j][k][i];
		}
		fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
		fourierfit(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,chi,Par);
		for(i=0;i<Par.nfilter;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=Par.nxb-Par.nfilter+1;i<Par.nxb;i++){
			phibar[i] = 0.0;
			phitilde[i] = 0.0;
			chibar[i] = 0.0;
			chitilde[i] = 0.0;
		}
		for(i=0;i<Par.nxb;i++){
			gij.yz[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,i,Par);
			Kij.yz[j][k][i] = fouriereval(Par.xleftend,Par.xrightend,chibar,chitilde,Par.nxb,i,Par);
		}
	
	}}
}


  return;
}
  
