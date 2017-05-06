#include "DeclareFunctions.h"

void AllocateStructs(Params *Par)
{
  int i,j,k;
  char string[30];
  FILE *fp;

//This reads params.dat and stores the values in the Params Structure 'Par'
//For options look at README_params.txt

  fp = fopen("params.dat", "r");   

// fscanf reads the number plus a string (to bypass the comment line)
  fscanf(fp,"%d %s\n",&Par->nxb,string);
  fscanf(fp,"%d %s\n",&Par->nyb,string);
  fscanf(fp,"%d %s\n",&Par->nzb,string);
  fscanf(fp,"%d %s\n",&Par->nfilter,string);
  fscanf(fp,"%lf %s\n",&Par->mass,string);
  fscanf(fp,"%d %s\n",&Par->ilapse,string);
  fscanf(fp,"%d %s\n",&Par->elapse,string);
  fscanf(fp,"%d %s\n",&Par->ishift,string);
  fscanf(fp,"%d %s\n",&Par->eshift,string);
  fscanf(fp,"%d %s\n",&Par->igij,string);
  fscanf(fp,"%d %s\n",&Par->iKij,string);
  fscanf(fp,"%d %s\n",&Par->bcflag_lapse,string);
  fscanf(fp,"%d %s\n",&Par->bcflag_shift,string);
  fscanf(fp,"%d %s\n",&Par->bcflag_gij,string);
  fscanf(fp,"%d %s\n",&Par->bcflag_Kij,string);
  fscanf(fp,"%d %s\n",&Par->start,string);
  fscanf(fp,"%d %s\n",&Par->timesteps,string);
  fscanf(fp,"%d %s\n",&Par->outsteps,string);
  fscanf(fp,"%d %s\n",&Par->dumpsteps,string);
  fscanf(fp,"%d %s\n",&Par->iadm,string);
  fscanf(fp,"%d %s\n",&Par->icon,string);
  fscanf(fp,"%lf %s\n",&Par->dt,string);
  fscanf(fp,"%lf %s\n",&Par->xleftend,string);
  fscanf(fp,"%lf %s\n",&Par->xrightend,string);
  fscanf(fp,"%lf %s\n",&Par->yleftend,string);
  fscanf(fp,"%lf %s\n",&Par->yrightend,string);
  fscanf(fp,"%lf %s\n",&Par->zleftend,string);
  fscanf(fp,"%lf %s\n",&Par->zrightend,string);
  fclose(fp);

  printf("%d %d %d %lf\n",Par->nxb,Par->nyb,Par->nxb,Par->mass);
  printf("%d %d %d %d %d\n",Par->ilapse,Par->ishift,Par->igij,Par->iKij,Par->bcflag_gij);
  printf("%d %d %lf \n",Par->timesteps,Par->outsteps,Par->dt);
  printf("%lf %lf %lf %lf %lf %lf\n",Par->xleftend,Par->xrightend,Par->yleftend,Par->yrightend,Par->zleftend,Par->zrightend);
  if(Par->iadm == 0) printf("Smarr-York equations\n");
  if(Par->iadm == 1) printf("True ADM equations\n");
  if(Par->icon == 0) printf("Constraint control: off\n");
  if(Par->icon == 1) printf("Constraint control: on\n");

  Par->scalex = 4./pow(Par->xrightend - Par->xleftend,2);
  Par->scaley = 4./pow(Par->yrightend - Par->yleftend,2);
  Par->scalez = 4./pow(Par->zrightend - Par->zleftend,2);

 Par->cosfac = dVector2D(0,Par->nxb,0,Par->nxb);
 Par->sinfac = dVector2D(0,Par->nxb,0,Par->nxb);
 Par->x = dVector1D(0,Par->nxb);
 Par->y = dVector1D(0,Par->nxb);
 Par->z = dVector1D(0,Par->nxb);

 for (j=0;j<Par->nxb;j++){
  for(k=0;k<Par->nxb;k++){
  Par->cosfac[j][k] = cos(2.0*M_PI*(j-(Par->nxb)/2.0)*(k-(Par->nxb)/2.0)/Par->nxb);
  Par->sinfac[j][k] = sin(2.0*M_PI*(j-(Par->nxb)/2.0)*(k-(Par->nxb)/2.0)/Par->nxb);
 }
}

for (i=0;i<Par->nxb;i++){
 Par->x[i] = 2.0*M_PI*(i-(Par->nxb)/2.0)/(Par->nxb);
 Par->y[i] = 2.0*M_PI*(i-(Par->nxb)/2.0)/(Par->nxb);
 Par->z[i] = 2.0*M_PI*(i-(Par->nxb)/2.0)/(Par->nxb);
}

printf("%lf %lf 0 nxb\n",Par->x[0],Par->x[Par->nxb]);

}
