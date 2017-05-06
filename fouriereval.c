#include "DeclareFunctions.h"

double fouriereval(double a, double b, double phibar[],double phitilde[],int m, int i,Params Par)
{
double phi;
int k;

 phi = 0.0;
 for(k=0;k<m;k++){
	phi += Par.cosfac[k][i]*phibar[k] - Par.sinfac[k][i]*phitilde[k];
 }
 phi = phi/(2.0*M_PI);
	
return phi;

}
