void fourierder(double a,double b,double phibar[],double phitilde[],double psibar[],double psitilde[],int n)
{
 int k;
 double no2;

 no2 = n/2.0;
 
 psibar[0] = 0.0;
 psitilde[0] = 0.0;
 for(k=1;k<n;k++){
	psibar[k] = -(k-no2)*phitilde[k];
 	psitilde[k] = (k-no2)*phibar[k];
 }

}

