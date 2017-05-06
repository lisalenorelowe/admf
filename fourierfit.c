#include "DeclareFunctions.h"

void fourierfit(double a, double b, double phibar[],double phitilde[], int n,double phi[], Params Par)
{

int i, k, nrm;
double twopin;

nrm = 8;  //number of modes to remove in addition to the critical mode (k = 0)
twopin = 2.0*M_PI/n;

for(k=0;k<n;k++){
	phibar[k] = 0.0;
	phitilde[k] = 0.0;
	for(i=0;i<n;i++){
		phibar[k] += Par.cosfac[k][i]*phi[i];
		phitilde[k] += Par.sinfac[k][i]*phi[i];
	}
	phibar[k] = twopin*phibar[k];
	phitilde[k] = -twopin*phitilde[k];
//	if(fabs(phibar[k])<=1.0e-12) phibar[k] = 0.0;
//	if(fabs(phitilde[k])<=1.0e-12) phitilde[k] = 0.0;
}

/*
//FILTER...
for(k=0;k<=nrm;k++){
	phibar[k] = 0.0;
	phitilde[k] = 0.0;
}
for(k=n-nrm;k<n;k++){
	phibar[k] = 0.0;
	phitilde[k] = 0.0;
}
*/

}
 
