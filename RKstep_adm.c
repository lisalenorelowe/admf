#include "DeclareFunctions.h"

void Rkstep_adm(double t,double dt,Scalar alpha,Vector shift,Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,Params Par,Scalar Phi,Vector Eks,Scalar hc,Vector mc,Connection C2,Tensor Ricci)
{
  Tensor g1,g2,g3,g4,Tgij,Tgup;
  Tensor k1,k2,k3,k4,TKij;
  Scalar Talpha,trK;
  Tensor LEks;
  double ***a1,***a2,***a3,***a4;
  int m,i,j,k, admflag, ccflag;
  double hcL2, mcL2, tmp;
  double KupKij[XMAX][YMAX][ZMAX];
  double DxK[XMAX][YMAX][ZMAX], DyK[XMAX][YMAX][ZMAX], DzK[XMAX][YMAX][ZMAX];

 a1    = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 a2    = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 a3    = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 a4    = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Talpha.s = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 trK.s = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k1.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k1.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k1.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k1.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k1.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k1.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k2.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k2.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k2.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k2.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k2.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k2.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k3.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k3.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k3.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k3.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k3.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k3.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k4.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k4.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k4.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k4.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k4.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 k4.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g1.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g1.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g1.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g1.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g1.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g1.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g2.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g2.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g2.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g2.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g2.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g2.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g3.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g3.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g3.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g3.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g3.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g3.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g4.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g4.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g4.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g4.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g4.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 g4.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

 Tgij.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Tgij.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Tgij.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Tgij.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Tgij.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Tgij.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

 TKij.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 TKij.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 TKij.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 TKij.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 TKij.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 TKij.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

 Tgup.xx   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Tgup.yy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Tgup.zz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Tgup.xy   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Tgup.xz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Tgup.yz   = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

 LEks.xx = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 LEks.xy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 LEks.xz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 LEks.yy = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 LEks.yz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 LEks.zz = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);


     admflag = 0;
     if(Par.iadm == 1) admflag = 1;

     ccflag = 0;
     if(Par.icon==1) ccflag = 1;

     if(ccflag == 0){
        for(i=0;i<Par.nxb;i++){
        for(j=0;j<Par.nyb;j++){
        for(k=0;k<Par.nzb;k++){
           LEks.xx[i][j][k] = 0.0;
           LEks.xy[i][j][k] = 0.0;
           LEks.xz[i][j][k] = 0.0;
           LEks.yy[i][j][k] = 0.0;
           LEks.yz[i][j][k] = 0.0;
           LEks.zz[i][j][k] = 0.0;
        }}}
     }     

     if(Par.ilapse==30) ResetLapse(t,alpha,Par);
     ginv(gij,&gup,Par);
     hcL2 = HamCon(gij,gup,Kij,Kup,hc,Par,KupKij,C2,Ricci);
     mcL2 = MomCon(gij,gup,Kij,Kup,mc,Par,DxK,DyK,DzK,C2,Ricci);
     Derivs_admg(t, alpha ,shift, Kij, gij, &g1, Par);
     Derivs_admK(t, alpha ,shift, gij, gup, Kij, trK, &k1, Par, C2, Ricci);
     Derivs_admlapse(t, alpha , shift, gup, Kij, a1, Par);
     if(Par.icon==1){
          printf("Calling ConCon \n");
          ConCon(gij,gup,Kij,Kup,alpha,shift,hc,mc,KupKij,DxK,DyK,DzK,Phi,Eks,Par,C2,Ricci);
          Loperator(Eks,LEks,gij,gup,C2,Par);
          }
        for(i=0;i<Par.nxb;i++){
        for(j=0;j<Par.nyb;j++){
        for(k=0;k<Par.nzb;k++){
          tmp = admflag*0.25*alpha.s[i][j][k]*hc.s[i][j][k];
	  k1.xx[i][j][k] = k1.xx[i][j][k] - tmp*gij.xx[i][j][k] + LEks.xx[i][j][k] - 2.0*(Kij.xx[i][j][k] - trK.s[i][j][k]*gij.xx[i][j][k])*Phi.s[i][j][k];
	  k1.yy[i][j][k] = k1.yy[i][j][k] - tmp*gij.yy[i][j][k] + LEks.yy[i][j][k] - 2.0*(Kij.yy[i][j][k] - trK.s[i][j][k]*gij.yy[i][j][k])*Phi.s[i][j][k];
	  k1.zz[i][j][k] = k1.zz[i][j][k] - tmp*gij.zz[i][j][k] + LEks.zz[i][j][k] - 2.0*(Kij.zz[i][j][k] - trK.s[i][j][k]*gij.zz[i][j][k])*Phi.s[i][j][k];
	  k1.xy[i][j][k] = k1.xy[i][j][k] - tmp*gij.xy[i][j][k] + LEks.xy[i][j][k] - 2.0*(Kij.xy[i][j][k] - trK.s[i][j][k]*gij.xy[i][j][k])*Phi.s[i][j][k];
	  k1.xz[i][j][k] = k1.xz[i][j][k] - tmp*gij.xz[i][j][k] + LEks.xz[i][j][k] - 2.0*(Kij.xz[i][j][k] - trK.s[i][j][k]*gij.xz[i][j][k])*Phi.s[i][j][k];
	  k1.yz[i][j][k] = k1.yz[i][j][k] - tmp*gij.yz[i][j][k] + LEks.yz[i][j][k] - 2.0*(Kij.yz[i][j][k] - trK.s[i][j][k]*gij.yz[i][j][k])*Phi.s[i][j][k];

          TKij.xx[i][j][k] = Kij.xx[i][j][k] + .5*dt*k1.xx[i][j][k];
          TKij.yy[i][j][k] = Kij.yy[i][j][k] + .5*dt*k1.yy[i][j][k];
          TKij.zz[i][j][k] = Kij.zz[i][j][k] + .5*dt*k1.zz[i][j][k];
          TKij.xy[i][j][k] = Kij.xy[i][j][k] + .5*dt*k1.xy[i][j][k];
          TKij.xz[i][j][k] = Kij.xz[i][j][k] + .5*dt*k1.xz[i][j][k];
          TKij.yz[i][j][k] = Kij.yz[i][j][k] + .5*dt*k1.yz[i][j][k];

          Tgij.xx[i][j][k] = gij.xx[i][j][k] + .5*dt*(g1.xx[i][j][k] + 4.0*Phi.s[i][j][k]*gij.xx[i][j][k]);
          Tgij.yy[i][j][k] = gij.yy[i][j][k] + .5*dt*(g1.yy[i][j][k] + 4.0*Phi.s[i][j][k]*gij.yy[i][j][k]);
          Tgij.zz[i][j][k] = gij.zz[i][j][k] + .5*dt*(g1.zz[i][j][k] + 4.0*Phi.s[i][j][k]*gij.zz[i][j][k]);
          Tgij.xy[i][j][k] = gij.xy[i][j][k] + .5*dt*(g1.xy[i][j][k] + 4.0*Phi.s[i][j][k]*gij.xy[i][j][k]);
          Tgij.xz[i][j][k] = gij.xz[i][j][k] + .5*dt*(g1.xz[i][j][k] + 4.0*Phi.s[i][j][k]*gij.xz[i][j][k]);
          Tgij.yz[i][j][k] = gij.yz[i][j][k] + .5*dt*(g1.yz[i][j][k] + 4.0*Phi.s[i][j][k]*gij.yz[i][j][k]);

          Talpha.s[i][j][k]  = alpha.s[i][j][k] + .5*dt*a1[i][j][k];
	}}}

     if(Par.ilapse==30) ResetLapse(t+.5*dt,Talpha,Par);
     ginv(Tgij,&Tgup,Par);
     hcL2 = HamCon(Tgij,Tgup,TKij,Kup,hc,Par,KupKij,C2,Ricci);
     mcL2 = MomCon(Tgij,Tgup,TKij,Kup,mc,Par,DxK,DyK,DzK,C2,Ricci);
     Derivs_admg(t+.5*dt, Talpha, shift, TKij, Tgij, &g2, Par);
     Derivs_admK(t+.5*dt, Talpha, shift, Tgij, Tgup, TKij, trK, &k2, Par, C2,Ricci);
     Derivs_admlapse(t+.5*dt, Talpha, shift, Tgup, TKij, a2, Par);
     if(Par.icon==1){
          printf("Calling ConCon\n");
          ConCon(Tgij,Tgup,TKij,Kup,Talpha,shift,hc,mc,KupKij,DxK,DyK,DzK,Phi,Eks,Par,C2,Ricci);
          Loperator(Eks,LEks,Tgij,Tgup,C2,Par);
          }
        for(i=0;i<Par.nxb;i++){
        for(j=0;j<Par.nyb;j++){
        for(k=0;k<Par.nzb;k++){
        tmp = admflag*0.25*Talpha.s[i][j][k]*hc.s[i][j][k];
	k2.xx[i][j][k] = k2.xx[i][j][k] - tmp*Tgij.xx[i][j][k] + LEks.xx[i][j][k] - 2.0*(TKij.xx[i][j][k] - trK.s[i][j][k]*Tgij.xx[i][j][k])*Phi.s[i][j][k];
	k2.yy[i][j][k] = k2.yy[i][j][k] - tmp*Tgij.yy[i][j][k] + LEks.yy[i][j][k] - 2.0*(TKij.yy[i][j][k] - trK.s[i][j][k]*Tgij.yy[i][j][k])*Phi.s[i][j][k];
	k2.zz[i][j][k] = k2.zz[i][j][k] - tmp*Tgij.zz[i][j][k] + LEks.zz[i][j][k] - 2.0*(TKij.zz[i][j][k] - trK.s[i][j][k]*Tgij.zz[i][j][k])*Phi.s[i][j][k];
	k2.xy[i][j][k] = k2.xy[i][j][k] - tmp*Tgij.xy[i][j][k] + LEks.xy[i][j][k] - 2.0*(TKij.xy[i][j][k] - trK.s[i][j][k]*Tgij.xy[i][j][k])*Phi.s[i][j][k];
	k2.xz[i][j][k] = k2.xz[i][j][k] - tmp*Tgij.xz[i][j][k] + LEks.xz[i][j][k] - 2.0*(TKij.xz[i][j][k] - trK.s[i][j][k]*Tgij.xz[i][j][k])*Phi.s[i][j][k];
	k2.yz[i][j][k] = k2.yz[i][j][k] - tmp*Tgij.yz[i][j][k] + LEks.yz[i][j][k] - 2.0*(TKij.yz[i][j][k] - trK.s[i][j][k]*Tgij.yz[i][j][k])*Phi.s[i][j][k];

        TKij.xx[i][j][k] = Kij.xx[i][j][k] + .5*dt*k2.xx[i][j][k];
        TKij.yy[i][j][k] = Kij.yy[i][j][k] + .5*dt*k2.yy[i][j][k];
        TKij.zz[i][j][k] = Kij.zz[i][j][k] + .5*dt*k2.zz[i][j][k];
        TKij.xy[i][j][k] = Kij.xy[i][j][k] + .5*dt*k2.xy[i][j][k];
        TKij.xz[i][j][k] = Kij.xz[i][j][k] + .5*dt*k2.xz[i][j][k];
        TKij.yz[i][j][k] = Kij.yz[i][j][k] + .5*dt*k2.yz[i][j][k];

        Tgij.xx[i][j][k] = gij.xx[i][j][k] + .5*dt*(g2.xx[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.xx[i][j][k]);
        Tgij.yy[i][j][k] = gij.yy[i][j][k] + .5*dt*(g2.yy[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.yy[i][j][k]);
        Tgij.zz[i][j][k] = gij.zz[i][j][k] + .5*dt*(g2.zz[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.zz[i][j][k]);
        Tgij.xy[i][j][k] = gij.xy[i][j][k] + .5*dt*(g2.xy[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.xy[i][j][k]);
        Tgij.xz[i][j][k] = gij.xz[i][j][k] + .5*dt*(g2.xz[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.xz[i][j][k]);
        Tgij.yz[i][j][k] = gij.yz[i][j][k] + .5*dt*(g2.yz[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.yz[i][j][k]);

        Talpha.s[i][j][k] = alpha.s[i][j][k] + .5*dt*a2[i][j][k];

	}}}

     if(Par.ilapse==30) ResetLapse(t+.5*dt,Talpha,Par);
     ginv(Tgij,&Tgup,Par);
     hcL2 = HamCon(Tgij,Tgup,TKij,Kup,hc,Par,KupKij,C2,Ricci);
     mcL2 = MomCon(Tgij,Tgup,TKij,Kup,mc,Par,DxK,DyK,DzK,C2,Ricci);
     Derivs_admg(t+.5*dt, Talpha, shift, TKij, Tgij, &g3, Par);
     Derivs_admK(t+.5*dt, Talpha, shift, Tgij, Tgup, TKij, trK, &k3, Par, C2,Ricci);
     Derivs_admlapse(t+.5*dt, Talpha, shift, Tgup, TKij, a3, Par);
     if(Par.icon==1){
          printf("Calling ConCon\n");
          ConCon(Tgij,Tgup,TKij,Kup,Talpha,shift,hc,mc,KupKij,DxK,DyK,DzK,Phi,Eks,Par,C2,Ricci);
          Loperator(Eks,LEks,Tgij,Tgup,C2,Par);
          }
        for(i=0;i<Par.nxb;i++){
        for(j=0;j<Par.nyb;j++){
        for(k=0;k<Par.nzb;k++){
        tmp = admflag*0.25*Talpha.s[i][j][k]*hc.s[i][j][k];
	k3.xx[i][j][k] = k3.xx[i][j][k] - tmp*Tgij.xx[i][j][k] + LEks.xx[i][j][k] - 2.0*(TKij.xx[i][j][k] - trK.s[i][j][k]*Tgij.xx[i][j][k])*Phi.s[i][j][k];
	k3.yy[i][j][k] = k3.yy[i][j][k] - tmp*Tgij.yy[i][j][k] + LEks.yy[i][j][k] - 2.0*(TKij.yy[i][j][k] - trK.s[i][j][k]*Tgij.yy[i][j][k])*Phi.s[i][j][k];
	k3.zz[i][j][k] = k3.zz[i][j][k] - tmp*Tgij.zz[i][j][k] + LEks.zz[i][j][k] - 2.0*(TKij.zz[i][j][k] - trK.s[i][j][k]*Tgij.zz[i][j][k])*Phi.s[i][j][k];
	k3.xy[i][j][k] = k3.xy[i][j][k] - tmp*Tgij.xy[i][j][k] + LEks.xy[i][j][k] - 2.0*(TKij.xy[i][j][k] - trK.s[i][j][k]*Tgij.xy[i][j][k])*Phi.s[i][j][k];
	k3.xz[i][j][k] = k3.xz[i][j][k] - tmp*Tgij.xz[i][j][k] + LEks.xz[i][j][k] - 2.0*(TKij.xz[i][j][k] - trK.s[i][j][k]*Tgij.xz[i][j][k])*Phi.s[i][j][k];
	k3.yz[i][j][k] = k3.yz[i][j][k] - tmp*Tgij.yz[i][j][k] + LEks.yz[i][j][k] - 2.0*(TKij.yz[i][j][k] - trK.s[i][j][k]*Tgij.yz[i][j][k])*Phi.s[i][j][k];

        TKij.xx[i][j][k] = Kij.xx[i][j][k] + dt*k3.xx[i][j][k];
        TKij.yy[i][j][k] = Kij.yy[i][j][k] + dt*k3.yy[i][j][k];
        TKij.zz[i][j][k] = Kij.zz[i][j][k] + dt*k3.zz[i][j][k];
        TKij.xy[i][j][k] = Kij.xy[i][j][k] + dt*k3.xy[i][j][k];
        TKij.xz[i][j][k] = Kij.xz[i][j][k] + dt*k3.xz[i][j][k];
        TKij.yz[i][j][k] = Kij.yz[i][j][k] + dt*k3.yz[i][j][k];

        Tgij.xx[i][j][k] = gij.xx[i][j][k] + dt*(g3.xx[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.xx[i][j][k]);
        Tgij.yy[i][j][k] = gij.yy[i][j][k] + dt*(g3.yy[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.yy[i][j][k]);
        Tgij.zz[i][j][k] = gij.zz[i][j][k] + dt*(g3.zz[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.zz[i][j][k]);
        Tgij.xy[i][j][k] = gij.xy[i][j][k] + dt*(g3.xy[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.xy[i][j][k]);
        Tgij.xz[i][j][k] = gij.xz[i][j][k] + dt*(g3.xz[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.xz[i][j][k]);
        Tgij.yz[i][j][k] = gij.yz[i][j][k] + dt*(g3.yz[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.yz[i][j][k]);

        Talpha.s[i][j][k]  = alpha.s[i][j][k]  + dt*a3[i][j][k];

      }}}

     if(Par.ilapse==30) ResetLapse(t+dt,Talpha,Par);
     ginv(Tgij,&Tgup,Par);
     hcL2 = HamCon(Tgij,Tgup,TKij,Kup,hc,Par,KupKij,C2,Ricci);
     mcL2 = MomCon(Tgij,Tgup,TKij,Kup,mc,Par,DxK,DyK,DzK,C2,Ricci);
     Derivs_admg(t+dt,Talpha, shift, TKij, Tgij, &g4, Par);
     Derivs_admK(t+dt, Talpha, shift, Tgij, Tgup, TKij, trK, &k4, Par, C2,Ricci);
     Derivs_admlapse(t+dt, Talpha, shift, Tgup, TKij, a4, Par);
     if(Par.icon==1){
          printf("Calling ConCon \n");
          ConCon(Tgij,Tgup,TKij,Kup,Talpha,shift,hc,mc,KupKij,DxK,DyK,DzK,Phi,Eks,Par,C2,Ricci);
          Loperator(Eks,LEks,Tgij,Tgup,C2,Par);
          }
        for(i=0;i<Par.nxb;i++){
        for(j=0;j<Par.nyb;j++){
        for(k=0;k<Par.nzb;k++){
	tmp = admflag*0.25*Talpha.s[i][j][k]*hc.s[i][j][k];
	k4.xx[i][j][k] = k4.xx[i][j][k] - tmp*Tgij.xx[i][j][k] + LEks.xx[i][j][k] - 2.0*(TKij.xx[i][j][k] - trK.s[i][j][k]*Tgij.xx[i][j][k])*Phi.s[i][j][k];
	k4.yy[i][j][k] = k4.yy[i][j][k] - tmp*Tgij.yy[i][j][k] + LEks.yy[i][j][k] - 2.0*(TKij.yy[i][j][k] - trK.s[i][j][k]*Tgij.yy[i][j][k])*Phi.s[i][j][k];
	k4.zz[i][j][k] = k4.zz[i][j][k] - tmp*Tgij.zz[i][j][k] + LEks.zz[i][j][k] - 2.0*(TKij.zz[i][j][k] - trK.s[i][j][k]*Tgij.zz[i][j][k])*Phi.s[i][j][k];
	k4.xy[i][j][k] = k4.xy[i][j][k] - tmp*Tgij.xy[i][j][k] + LEks.xy[i][j][k] - 2.0*(TKij.xy[i][j][k] - trK.s[i][j][k]*Tgij.xy[i][j][k])*Phi.s[i][j][k];
	k4.xz[i][j][k] = k4.xz[i][j][k] - tmp*Tgij.xz[i][j][k] + LEks.xz[i][j][k] - 2.0*(TKij.xz[i][j][k] - trK.s[i][j][k]*Tgij.xz[i][j][k])*Phi.s[i][j][k];
	k4.yz[i][j][k] = k4.yz[i][j][k] - tmp*Tgij.yz[i][j][k] + LEks.yz[i][j][k] - 2.0*(TKij.yz[i][j][k] - trK.s[i][j][k]*Tgij.yz[i][j][k])*Phi.s[i][j][k];

        Kij.xx[i][j][k] = Kij.xx[i][j][k]+ dt*(1./6.*k1.xx[i][j][k] + 1./3.*k2.xx[i][j][k] + 1./3.*k3.xx[i][j][k] + 1./6.*k4.xx[i][j][k]);
        Kij.yy[i][j][k] = Kij.yy[i][j][k]+ dt*(1./6.*k1.yy[i][j][k] + 1./3.*k2.yy[i][j][k] + 1./3.*k3.yy[i][j][k] + 1./6.*k4.yy[i][j][k]);
        Kij.zz[i][j][k] = Kij.zz[i][j][k]+ dt*(1./6.*k1.zz[i][j][k] + 1./3.*k2.zz[i][j][k] + 1./3.*k3.zz[i][j][k] + 1./6.*k4.zz[i][j][k]);
        Kij.xy[i][j][k] = Kij.xy[i][j][k]+ dt*(1./6.*k1.xy[i][j][k] + 1./3.*k2.xy[i][j][k] + 1./3.*k3.xy[i][j][k] + 1./6.*k4.xy[i][j][k]);
        Kij.xz[i][j][k] = Kij.xz[i][j][k]+ dt*(1./6.*k1.xz[i][j][k] + 1./3.*k2.xz[i][j][k] + 1./3.*k3.xz[i][j][k] + 1./6.*k4.xz[i][j][k]);
        Kij.yz[i][j][k] = Kij.yz[i][j][k]+ dt*(1./6.*k1.yz[i][j][k] + 1./3.*k2.yz[i][j][k] + 1./3.*k3.yz[i][j][k] + 1./6.*k4.yz[i][j][k]);

        g4.xx[i][j][k] = g4.xx[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.xx[i][j][k];
        g4.yy[i][j][k] = g4.yy[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.yy[i][j][k];
        g4.zz[i][j][k] = g4.zz[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.zz[i][j][k];
        g4.xy[i][j][k] = g4.xy[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.xy[i][j][k];
        g4.xz[i][j][k] = g4.xz[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.xz[i][j][k];
        g4.yz[i][j][k] = g4.yz[i][j][k] + 4.0*Phi.s[i][j][k]*Tgij.yz[i][j][k];

        gij.xx[i][j][k] = gij.xx[i][j][k]+ dt*(1./6.*g1.xx[i][j][k] + 1./3.*g2.xx[i][j][k] + 1./3.*g3.xx[i][j][k] + 1./6.*g4.xx[i][j][k]);
        gij.yy[i][j][k] = gij.yy[i][j][k]+ dt*(1./6.*g1.yy[i][j][k] + 1./3.*g2.yy[i][j][k] + 1./3.*g3.yy[i][j][k] + 1./6.*g4.yy[i][j][k]);
        gij.zz[i][j][k] = gij.zz[i][j][k]+ dt*(1./6.*g1.zz[i][j][k] + 1./3.*g2.zz[i][j][k] + 1./3.*g3.zz[i][j][k] + 1./6.*g4.zz[i][j][k]);
        gij.xy[i][j][k] = gij.xy[i][j][k]+ dt*(1./6.*g1.xy[i][j][k] + 1./3.*g2.xy[i][j][k] + 1./3.*g3.xy[i][j][k] + 1./6.*g4.xy[i][j][k]);
        gij.xz[i][j][k] = gij.xz[i][j][k]+ dt*(1./6.*g1.xz[i][j][k] + 1./3.*g2.xz[i][j][k] + 1./3.*g3.xz[i][j][k] + 1./6.*g4.xz[i][j][k]);
        gij.yz[i][j][k] = gij.yz[i][j][k]+ dt*(1./6.*g1.yz[i][j][k] + 1./3.*g2.yz[i][j][k] + 1./3.*g3.yz[i][j][k] + 1./6.*g4.yz[i][j][k]);

            if(Par.elapse!=0) alpha.s[i][j][k] = alpha.s[i][j][k]+ dt*(1./6.*a1[i][j][k] + 1./3.*a2[i][j][k] + 1./3.*a3[i][j][k] + 1./6.*a4[i][j][k]);

     }}}

     if(Par.ilapse==30) ResetLapse(t+dt,alpha,Par);


  free_dVector3D(k1.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k1.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k1.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k1.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k1.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k1.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k2.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k2.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k2.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k2.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k2.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k2.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k3.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k3.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k3.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k3.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k3.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k3.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k4.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k4.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k4.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k4.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k4.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(k4.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g1.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g1.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g1.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g1.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g1.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g1.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g2.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g2.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g2.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g2.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g2.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g2.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g3.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g3.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g3.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g3.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g3.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g3.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g4.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g4.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g4.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g4.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g4.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(g4.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);

  free_dVector3D(Tgij.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(Tgij.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(Tgij.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(Tgij.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(Tgij.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(Tgij.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);

  free_dVector3D(TKij.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(TKij.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(TKij.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(TKij.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(TKij.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(TKij.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(Tgup.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(Tgup.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(Tgup.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(Tgup.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(Tgup.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(Tgup.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);

 free_dVector3D(a1,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(a2,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(a3,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(a4,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(Talpha.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(trK.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);

 free_dVector3D(LEks.xx,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(LEks.yy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(LEks.zz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(LEks.xy,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(LEks.xz,0,Par.nxb,0,Par.nyb,0,Par.nzb);
 free_dVector3D(LEks.yz,0,Par.nxb,0,Par.nyb,0,Par.nzb);


 return;

}

