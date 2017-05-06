#include "DeclareFunctions.h"

void ConCon(Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,Scalar alpha,Vector shift,Scalar hc,Vector mc,double KupKij[XMAX][YMAX][ZMAX],double DxK[XMAX][YMAX][ZMAX],double DyK[XMAX][YMAX][ZMAX],double DzK[XMAX][YMAX][ZMAX],Scalar Phi,Vector Eks,Params Par, Connection C2,Tensor Ricci)
{
 int i, j, k;
 Scalar detg;
 double Coeff0[XMAX][YMAX][ZMAX], Coeffx[XMAX][YMAX][ZMAX], Coeffy[XMAX][YMAX][ZMAX], Coeffz[XMAX][YMAX][ZMAX];
 Vector dhc, dmcx, dmcy, dmcz, dalpha, dbx, dby, dbz, dMupx, dMupy, dMupz;
 double Lam0[XMAX][YMAX][ZMAX], Lamx[XMAX][YMAX][ZMAX], Lamy[XMAX][YMAX][ZMAX], Lamz[XMAX][YMAX][ZMAX];
 double trK, MDalpha, divM, rtg, Liehc, Liemcx, Liemcy, Liemcz; 

detg.s = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dalpha.x     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dalpha.y     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dalpha.z     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbx.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbx.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbx.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dby.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dby.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dby.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dbz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dhc.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dhc.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dhc.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dmcx.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dmcx.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dmcx.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dmcy.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dmcy.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dmcy.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dmcz.x=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dmcz.y=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
dmcz.z=dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Find the derivatives of alpha, shift, hc, and mc*
  FunctionFirstDerivs(alpha.s,&dalpha,Par);
  FunctionFirstDerivs(shift.x,&dbx,Par);
  FunctionFirstDerivs(shift.y,&dby,Par);
  FunctionFirstDerivs(shift.z,&dbz,Par);
  FunctionFirstDerivs(hc.s,&dhc,Par);
  FunctionFirstDerivs(mc.x,&dmcx,Par);
  FunctionFirstDerivs(mc.y,&dmcy,Par);
  FunctionFirstDerivs(mc.z,&dmcz,Par);

//Construct the right-hand sides. First the terms that come from the old evolution equations 
//for the constraints.
 if(Par.iadm == 0){
   for(i=0;i<Par.nxb;i++){
      for(j=0;j<Par.nyb;j++){
         for(k=0;k<Par.nzb;k++){
            trK = gup.xx[i][j][k]*Kij.xx[i][j][k] + gup.yy[i][j][k]*Kij.yy[i][j][k]
                         + gup.zz[i][j][k]*Kij.zz[i][j][k] + 2.0*gup.xy[i][j][k]*Kij.xy[i][j][k] 
                         + 2.0*gup.xz[i][j][k]*Kij.xz[i][j][k] + 2.0*gup.yz[i][j][k]*Kij.yz[i][j][k];
	    MDalpha = mc.x[i][j][k]*gup.xx[i][j][k]*dalpha.x[i][j][k]
                                + mc.x[i][j][k]*gup.xy[i][j][k]*dalpha.y[i][j][k]
                                + mc.x[i][j][k]*gup.xz[i][j][k]*dalpha.z[i][j][k]
                                + mc.y[i][j][k]*gup.xy[i][j][k]*dalpha.x[i][j][k]
                                + mc.y[i][j][k]*gup.yy[i][j][k]*dalpha.y[i][j][k]
                                + mc.y[i][j][k]*gup.yz[i][j][k]*dalpha.z[i][j][k]
                                + mc.z[i][j][k]*gup.xz[i][j][k]*dalpha.x[i][j][k]
                                + mc.z[i][j][k]*gup.yz[i][j][k]*dalpha.y[i][j][k]
                                + mc.z[i][j][k]*gup.zz[i][j][k]*dalpha.z[i][j][k];
            divM = gup.xx[i][j][k]*(dmcx.x[i][j][k] - C2.xxx[i][j][k]*mc.x[i][j][k] 
                                                    - C2.yxx[i][j][k]*mc.y[i][j][k]
                                                    - C2.zxx[i][j][k]*mc.z[i][j][k])
                 + gup.yy[i][j][k]*(dmcy.y[i][j][k] - C2.xyy[i][j][k]*mc.x[i][j][k] 
                                                    - C2.yyy[i][j][k]*mc.y[i][j][k]
                                                    - C2.zyy[i][j][k]*mc.z[i][j][k])
                 + gup.zz[i][j][k]*(dmcz.z[i][j][k] - C2.xzz[i][j][k]*mc.x[i][j][k] 
                                                    - C2.yzz[i][j][k]*mc.y[i][j][k]
                                                    - C2.zzz[i][j][k]*mc.z[i][j][k])
                 + gup.xy[i][j][k]*(dmcy.x[i][j][k] + dmcx.y[i][j][k] 
                                    - 2.0*C2.xxy[i][j][k]*mc.x[i][j][k]
                                    - 2.0*C2.yxy[i][j][k]*mc.y[i][j][k]
                                    - 2.0*C2.zxy[i][j][k]*mc.z[i][j][k])
                 + gup.xz[i][j][k]*(dmcz.x[i][j][k] + dmcx.z[i][j][k] 
                                    - 2.0*C2.xxz[i][j][k]*mc.x[i][j][k]
                                    - 2.0*C2.yxz[i][j][k]*mc.y[i][j][k]
                                    - 2.0*C2.zxz[i][j][k]*mc.z[i][j][k])
                 + gup.yz[i][j][k]*(dmcz.y[i][j][k] + dmcy.z[i][j][k] 
                                    - 2.0*C2.xyz[i][j][k]*mc.x[i][j][k]
                                    - 2.0*C2.yyz[i][j][k]*mc.y[i][j][k]
                                    - 2.0*C2.zyz[i][j][k]*mc.z[i][j][k]);
            Liehc = shift.x[i][j][k]*dhc.x[i][j][k] + shift.y[i][j][k]*dhc.y[i][j][k] 
                  + shift.z[i][j][k]*dhc.z[i][j][k];
            Liemcx = shift.x[i][j][k]*dmcx.x[i][j][k] + shift.y[i][j][k]*dmcx.y[i][j][k] 
                   + shift.z[i][j][k]*dmcx.z[i][j][k] + mc.x[i][j][k]*dbx.x[i][j][k] 
                   + mc.y[i][j][k]*dby.x[i][j][k] + mc.z[i][j][k]*dbz.x[i][j][k];
            Liemcy = shift.x[i][j][k]*dmcy.x[i][j][k] + shift.y[i][j][k]*dmcy.y[i][j][k] 
                   + shift.z[i][j][k]*dmcy.z[i][j][k] + mc.x[i][j][k]*dbx.y[i][j][k] 
                   + mc.y[i][j][k]*dby.y[i][j][k] + mc.z[i][j][k]*dbz.y[i][j][k];
            Liemcz = shift.x[i][j][k]*dmcz.x[i][j][k] + shift.y[i][j][k]*dmcz.y[i][j][k] 
                   + shift.z[i][j][k]*dmcz.z[i][j][k] + mc.x[i][j][k]*dbx.z[i][j][k] 
                   + mc.y[i][j][k]*dby.z[i][j][k] + mc.z[i][j][k]*dbz.z[i][j][k];
	    Lam0[i][j][k] = 2.0*alpha.s[i][j][k]*trK*hc.s[i][j][k] 
                          - 4.0*MDalpha -2.0*alpha.s[i][j][k]*divM + Liehc;
            Lamx[i][j][k] = -hc.s[i][j][k]*dalpha.x[i][j][k] + alpha.s[i][j][k]*trK*mc.x[i][j][k] 
	                   - 0.5*alpha.s[i][j][k]*dhc.x[i][j][k] + Liemcx;                     
            Lamy[i][j][k] = -hc.s[i][j][k]*dalpha.y[i][j][k] + alpha.s[i][j][k]*trK*mc.y[i][j][k] 
	                   - 0.5*alpha.s[i][j][k]*dhc.y[i][j][k] + Liemcy;                     
            Lamz[i][j][k] = -hc.s[i][j][k]*dalpha.z[i][j][k] + alpha.s[i][j][k]*trK*mc.z[i][j][k] 
	                   - 0.5*alpha.s[i][j][k]*dhc.z[i][j][k] + Liemcz;                     
   }}}
 }
 if(Par.iadm == 1){
   for(i=0;i<Par.nxb;i++){
      for(j=0;j<Par.nyb;j++){
         for(k=0;k<Par.nzb;k++){
            trK = gup.xx[i][j][k]*Kij.xx[i][j][k] + gup.yy[i][j][k]*Kij.yy[i][j][k]
                         + gup.zz[i][j][k]*Kij.zz[i][j][k] + 2.0*gup.xy[i][j][k]*Kij.xy[i][j][k] 
                         + 2.0*gup.xz[i][j][k]*Kij.xz[i][j][k] + 2.0*gup.yz[i][j][k]*Kij.yz[i][j][k];
	    MDalpha = mc.x[i][j][k]*gup.xx[i][j][k]*dalpha.x[i][j][k]
                                + mc.x[i][j][k]*gup.xy[i][j][k]*dalpha.y[i][j][k]
                                + mc.x[i][j][k]*gup.xz[i][j][k]*dalpha.z[i][j][k]
                                + mc.y[i][j][k]*gup.xy[i][j][k]*dalpha.x[i][j][k]
                                + mc.y[i][j][k]*gup.yy[i][j][k]*dalpha.y[i][j][k]
                                + mc.y[i][j][k]*gup.yz[i][j][k]*dalpha.z[i][j][k]
                                + mc.z[i][j][k]*gup.xz[i][j][k]*dalpha.x[i][j][k]
                                + mc.z[i][j][k]*gup.yz[i][j][k]*dalpha.y[i][j][k]
                                + mc.z[i][j][k]*gup.zz[i][j][k]*dalpha.z[i][j][k];
            divM = gup.xx[i][j][k]*(dmcx.x[i][j][k] - C2.xxx[i][j][k]*mc.x[i][j][k] 
                                                    - C2.yxx[i][j][k]*mc.y[i][j][k]
                                                    - C2.zxx[i][j][k]*mc.z[i][j][k])
                 + gup.yy[i][j][k]*(dmcy.y[i][j][k] - C2.xyy[i][j][k]*mc.x[i][j][k] 
                                                    - C2.yyy[i][j][k]*mc.y[i][j][k]
                                                    - C2.zyy[i][j][k]*mc.z[i][j][k])
                 + gup.zz[i][j][k]*(dmcz.z[i][j][k] - C2.xzz[i][j][k]*mc.x[i][j][k] 
                                                    - C2.yzz[i][j][k]*mc.y[i][j][k]
                                                    - C2.zzz[i][j][k]*mc.z[i][j][k])
                 + gup.xy[i][j][k]*(dmcy.x[i][j][k] + dmcx.y[i][j][k] 
                                    - 2.0*C2.xxy[i][j][k]*mc.x[i][j][k]
                                    - 2.0*C2.yxy[i][j][k]*mc.y[i][j][k]
                                    - 2.0*C2.zxy[i][j][k]*mc.z[i][j][k])
                 + gup.xz[i][j][k]*(dmcz.x[i][j][k] + dmcx.z[i][j][k] 
                                    - 2.0*C2.xxz[i][j][k]*mc.x[i][j][k]
                                    - 2.0*C2.yxz[i][j][k]*mc.y[i][j][k]
                                    - 2.0*C2.zxz[i][j][k]*mc.z[i][j][k])
                 + gup.yz[i][j][k]*(dmcz.y[i][j][k] + dmcy.z[i][j][k] 
                                    - 2.0*C2.xyz[i][j][k]*mc.x[i][j][k]
                                    - 2.0*C2.yyz[i][j][k]*mc.y[i][j][k]
                                    - 2.0*C2.zyz[i][j][k]*mc.z[i][j][k]);
            Liehc = shift.x[i][j][k]*dhc.x[i][j][k] + shift.y[i][j][k]*dhc.y[i][j][k] 
                  + shift.z[i][j][k]*dhc.z[i][j][k];
            Liemcx = shift.x[i][j][k]*dmcx.x[i][j][k] + shift.y[i][j][k]*dmcx.y[i][j][k] 
                   + shift.z[i][j][k]*dmcx.z[i][j][k] + mc.x[i][j][k]*dbx.x[i][j][k] 
                   + mc.y[i][j][k]*dby.x[i][j][k] + mc.z[i][j][k]*dbz.x[i][j][k];
            Liemcy = shift.x[i][j][k]*dmcy.x[i][j][k] + shift.y[i][j][k]*dmcy.y[i][j][k] 
                   + shift.z[i][j][k]*dmcy.z[i][j][k] + mc.x[i][j][k]*dbx.y[i][j][k] 
                   + mc.y[i][j][k]*dby.y[i][j][k] + mc.z[i][j][k]*dbz.y[i][j][k];
            Liemcz = shift.x[i][j][k]*dmcz.x[i][j][k] + shift.y[i][j][k]*dmcz.y[i][j][k] 
                   + shift.z[i][j][k]*dmcz.z[i][j][k] + mc.x[i][j][k]*dbx.z[i][j][k] 
                   + mc.y[i][j][k]*dby.z[i][j][k] + mc.z[i][j][k]*dbz.z[i][j][k];
	    Lam0[i][j][k] = alpha.s[i][j][k]*trK*hc.s[i][j][k] 
                          - 4.0*MDalpha -2.0*alpha.s[i][j][k]*divM + Liehc;
            Lamx[i][j][k] = -0.5*hc.s[i][j][k]*dalpha.x[i][j][k] + alpha.s[i][j][k]*trK*mc.x[i][j][k] 
	                    + Liemcx;                     
            Lamy[i][j][k] = -0.5*hc.s[i][j][k]*dalpha.y[i][j][k] + alpha.s[i][j][k]*trK*mc.y[i][j][k] 
	                    + Liemcy;                     
            Lamz[i][j][k] = -0.5*hc.s[i][j][k]*dalpha.z[i][j][k] + alpha.s[i][j][k]*trK*mc.z[i][j][k] 
	                    + Liemcz;                     
   }}}
 }

//Now subtract these from the new terms. This completes construction of the right-hand sides.
   for(i=0;i<Par.nxb;i++){
      for(j=0;j<Par.nyb;j++){
         for(k=0;k<Par.nzb;k++){
	    Lam0[i][j][k] = -0.1*hc.s[i][j][k] - Lam0[i][j][k];
            Lamx[i][j][k] = -0.1*mc.x[i][j][k] - Lamx[i][j][k];
            Lamy[i][j][k] = -0.1*mc.y[i][j][k] - Lamy[i][j][k];
            Lamz[i][j][k] = -0.1*mc.z[i][j][k] - Lamz[i][j][k];
   }}}

//Now solve the elliptic problem. First compute some of the coefficients for the linear operator:
  for(i=0;i<Par.nxb;i++){
    for(j=0;j<Par.nyb;j++){
      for(k=0;k<Par.nzb;k++){
  	Coeff0[i][j][k] = -4.0*hc.s[i][j][k] + 8.0*KupKij[i][j][k];
        Coeffx[i][j][k] = -6.0*mc.x[i][j][k] - 4.0*DxK[i][j][k];
        Coeffy[i][j][k] = -6.0*mc.y[i][j][k] - 4.0*DyK[i][j][k]; 
        Coeffz[i][j][k] = -6.0*mc.z[i][j][k] - 4.0*DzK[i][j][k]; 
  }}}

//Call Esolve
  Esolve(Phi,Eks,Coeff0,Coeffx,Coeffy,Coeffz,Lam0,Lamx,Lamy,Lamz,C2,gij,gup,Kij,Kup,Par);

free_dVector3D(detg.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dalpha.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dalpha.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dalpha.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dby.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dby.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dby.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dbz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dhc.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dhc.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dhc.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dmcx.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dmcx.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dmcx.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dmcy.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dmcy.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dmcy.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dmcz.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dmcz.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
free_dVector3D(dmcz.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

return;
}
