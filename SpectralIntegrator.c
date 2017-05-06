#include "DeclareFunctions.h"

int main(int argc, char *argv[])
{
 int i,j,start, imod3;
 char restart,outfile[100];
 double x,t,tp1;
 double KupKij[XMAX][YMAX][ZMAX];
 double DxK[XMAX][YMAX][ZMAX], DyK[XMAX][YMAX][ZMAX], DzK[XMAX][YMAX][ZMAX];
 double phi[XMAX], phibar[XMAX], phitilde[XMAX];
 Params Par;
 Scalar alpha, Phi;
 Vector shift, Eks;
 Tensor Kij,gij,Ricci,gup,Kup;
 Scalar hc;
 Vector mc; 
 Connection C2;
 FILE *fpc;
 FILE *fpm;
 FILE *fgxx;
 FILE *fKxx;
 FILE *fhc;
 FILE *fmcx;
 double hcL2,mcL2;


// Error Checking for Input
 if(argc<2){
  printf("Usage: a.out outputfile\n");
 return 0;
}

// Dynamically Allocate Structures
  AllocateStructs(&Par);

if(Par.start==0){
  InitializeLapse(&alpha,Par);
  InitializeShift(&shift,Par);
  InitializeMetric(&gij,Par);
  InitializeKij(&Kij,&gij,Par);
  InitializePhiandEks(&Phi,&Eks,Par);
  hc.s = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  mc.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  mc.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  mc.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.xxx = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.xxy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.xxz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.xyy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.xyz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.xzz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.yxx = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.yxy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.yxz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.yyy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.yyz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.yzz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.zxx = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.zxy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.zxz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.zyy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.zyz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  C2.zzz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Ricci.xx = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Ricci.xy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Ricci.xz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Ricci.yy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Ricci.yz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Ricci.zz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  gup.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  gup.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  gup.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  gup.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  gup.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  gup.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Kup.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Kup.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Kup.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Kup.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Kup.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  Kup.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
}else ReadRestartDump(&alpha,&shift,&gij,&Kij,Par);
 

// Print Initial Waveform
  sprintf(outfile,"%s",argv[1]);
  OutputStep_gij(outfile,Par.start*Par.dt,gij,Kij,Par);
  printf("\n");
  printf("t =       hcL2 =           mcL2 =\n");
  ginv(gij,&gup,Par);
  hcL2 = HamCon(gij,gup,Kij,Kup,hc,Par,KupKij,C2,Ricci);
  mcL2 = MomCon(gij,gup,Kij,Kup,mc,Par,DxK,DyK,DzK,C2,Ricci);
  printf("%lf %1.13e %1.13e \n",Par.start*Par.dt,hcL2,mcL2);
  OutputStep_S("hamcon",Par.start*Par.dt,hc,Par);
  OutputStep_V("momcon",Par.start*Par.dt,mc,Par);

  fpc = fopen("HCfile","a");
  fprintf(fpc,"%f, %2.14f \n", Par.start*Par.dt,hcL2);
  fclose(fpc);
  fpm = fopen("MCfile","a");
  fprintf(fpm,"%f, %2.14f \n", Par.start*Par.dt,mcL2);
  fclose(fpm);

// Loop over time intervals
  for(i=Par.start;i<Par.timesteps;i++){
  t = i*Par.dt;
  printf("\n");
//  printf("lapse at origin = %2.10f \n", alpha.s[Par.nxb/2][Par.nyb/2][Par.nzb/2]);
//  printf("shift at origin = %2.10f %2.10f %2.10f \n", shift.x[Par.nxb/2][Par.nyb/2][Par.nzb/2], 
//	shift.y[Par.nxb/2][Par.nyb/2][Par.nzb/2], shift.z[Par.nxb/2][Par.nyb/2][Par.nzb/2]);
  IntegrateStep_adm(alpha,shift,gij,gup,Kij,Kup,t,Par,Phi,Eks,hc,mc,C2,Ricci);
  imod3 = i%3;
  Filter(gij,Kij,Par,imod3);
  imod3 = (i+1)%3;
  Filter(gij,Kij,Par,imod3);
  imod3 = (i+2)%3;
  Filter(gij,Kij,Par,imod3);

//Print New Values
   if((i+1)%Par.outsteps==0) {
    tp1 = (i+1)*Par.dt;
    OutputStep_gij(outfile,tp1,gij,Kij,Par);
    ginv(gij,&gup,Par);
    hcL2 = HamCon(gij,gup,Kij,Kup,hc,Par,KupKij,C2,Ricci);
    mcL2 = MomCon(gij,gup,Kij,Kup,mc,Par,DxK,DyK,DzK,C2,Ricci);
    printf("%lf %1.13e %1.13e \n",tp1,hcL2,mcL2);
    OutputStep_S("hamcon",tp1,hc,Par);
    OutputStep_V("momcon",tp1,mc,Par);
    fpc = fopen("HCfile","a");
    fprintf(fpc,"%f, %2.14f \n", tp1,hcL2);
    fclose(fpc);
    fpm = fopen("MCfile","a");
    fprintf(fpm,"%f, %2.14f \n", tp1,mcL2);
    fclose(fpm);


    fgxx = fopen("SCgxxfile","a");
    fprintf(fgxx,"%f ",tp1);
    for(j=0;j<Par.nxb;j++) phi[j] = gij.xx[j][Par.nyb/2][Par.nzb/2];
    fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
    for(j=0;j<Par.nxb;j++){
	 phi[j] = sqrt(phibar[j]*phibar[j] + phitilde[j]*phitilde[j]);
         fprintf(fgxx,"%2.10f ",phi[j]);
	 }
    fprintf(fgxx,"\n");
    fclose(fgxx);

    fKxx = fopen("SCKxxfile","a");
    fprintf(fKxx,"%f ",tp1);
    for(j=0;j<Par.nxb;j++) phi[j] = Kij.xx[j][Par.nyb/2][Par.nzb/2];
    fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
    for(j=0;j<Par.nxb;j++){
	 phi[j] = sqrt(phibar[j]*phibar[j] + phitilde[j]*phitilde[j]);
         fprintf(fKxx,"%2.10f ",phi[j]);
	 }
    fprintf(fKxx,"\n");
    fclose(fKxx);

    fhc = fopen("SCHCfile","a");
    fprintf(fhc,"%f ",tp1);
    for(j=0;j<Par.nxb;j++) phi[j] = hc.s[j][Par.nyb/2][Par.nzb/2];
    fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
    for(j=0;j<Par.nxb;j++){
	 phi[j] = sqrt(phibar[j]*phibar[j] + phitilde[j]*phitilde[j]);
         fprintf(fhc,"%2.10f ",phi[j]);
	 }
    fprintf(fhc,"\n");
    fclose(fhc);

    fmcx = fopen("SCMCxfile","a");
    fprintf(fmcx,"%f ",tp1);
    for(j=0;j<Par.nxb;j++) phi[j] = mc.x[j][Par.nyb/2][Par.nzb/2];
    fourierfit(Par.xleftend,Par.xrightend,phibar,phitilde,Par.nxb,phi,Par);
    for(j=0;j<Par.nxb;j++){
	 phi[j] = sqrt(phibar[j]*phibar[j] + phitilde[j]*phitilde[j]);
         fprintf(fmcx,"%2.10f ",phi[j]);
	 }
    fprintf(fmcx,"\n");
    fclose(fmcx);


   }
   if((i+1)%Par.dumpsteps==0) WriteRestartDump(i+1,alpha,shift,gij,Kij,Par);
  }// End Loop over time intervals



free_dVector1D(Par.x,0,Par.nxb);
free_dVector1D(Par.y,0,Par.nyb);
free_dVector1D(Par.z,0,Par.nzb);
free_dVector2D(Par.cosfac,0,Par.nxb,0,Par.nyb);
free_dVector3D(gij.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gij.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gij.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gij.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gij.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gij.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kij.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kij.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kij.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kij.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kij.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Kij.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(alpha.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(shift.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(shift.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(shift.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(hc.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(mc.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(mc.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(mc.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.xxx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.xxy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.xxz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.xyy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.xyz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.xzz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.yxx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.yxy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.yxz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.yyy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.yyz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.yzz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.zxx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.zxy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.zxz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.zyy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.zyz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(C2.zzz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricci.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricci.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricci.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricci.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricci.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(Ricci.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(gup.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);

return 0;

}
