#include "DeclareFunctions.h"

void Esolve(Scalar Phi,Vector Eks,double Coeff0[XMAX][YMAX][ZMAX],double Coeffx[XMAX][YMAX][ZMAX],double Coeffy[XMAX][YMAX][ZMAX],double Coeffz[XMAX][YMAX][ZMAX],double Lam0[XMAX][YMAX][ZMAX],double Lamx[XMAX][YMAX][ZMAX],double Lamy[XMAX][YMAX][ZMAX],double Lamz[XMAX][YMAX][ZMAX],Connection C2,Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,Params Par){

 int i,j,k;
 double dp,dp2,alp,ome,bet;
 Scalar rres0, rhat0, p0, v0, t0;
 Vector rres, rhat, p, v, t;
 double rho[101];
 int kay, kaymax;

  rres0.s = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  rhat0.s = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  p0.s = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  v0.s = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  t0.s = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  rres.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  rhat.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  p.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  v.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  t.x = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  rres.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  rhat.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  p.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  v.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  t.y = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  rres.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  rhat.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  p.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  v.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
  t.z = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

 kaymax = 10;


//See Tim Kelley's book, page 50.
//Step 1.
 Atimes(Phi,Eks,rres0,rres,Coeff0,Coeffx,Coeffy,Coeffz,gij,gup,Kij,Kup,C2,Par);
 for(i=0;i<Par.nxb;i++){
   for(j=0;j<Par.nyb;j++){
     for(k=0;k<Par.nzb;k++){
       rres0.s[i][j][k] = Lam0[i][j][k] - rres0.s[i][j][k];
       rhat0.s[i][j][k] = rres0.s[i][j][k];
       v0.s[i][j][k] = 0.0;
       p0.s[i][j][k] = 0.0;
       rres.x[i][j][k] = Lamx[i][j][k] - rres.x[i][j][k];
       rhat.x[i][j][k] = rres.x[i][j][k];
       v.x[i][j][k] = 0.0;
       p.x[i][j][k] = 0.0;
       rres.y[i][j][k] = Lamy[i][j][k] - rres.y[i][j][k];
       rhat.y[i][j][k] = rres.y[i][j][k];
       v.y[i][j][k] = 0.0;
       p.y[i][j][k] = 0.0;
       rres.z[i][j][k] = Lamz[i][j][k] - rres.z[i][j][k];
       rhat.z[i][j][k] = rres.z[i][j][k];
       v.z[i][j][k] = 0.0;
       p.z[i][j][k] = 0.0;
 }}}
 rho[0] = 1.0;
 alp = 1.0;
 ome = 1.0;
 kay = 0;
 rho[1] = dotproduct(rhat0,rhat,rres0,rres,Par);
//Step 2.
 while(kay < kaymax){
   kay = kay + 1;                              //(a)
   bet = (rho[kay]/rho[kay-1])*(alp/ome);      //(b)
   for(i=0;i<Par.nxb;i++){                     //(c)
     for(j=0;j<Par.nyb;j++){
       for(k=0;k<Par.nzb;k++){
         p0.s[i][j][k] = rres0.s[i][j][k] + bet*(p0.s[i][j][k] - ome*v0.s[i][j][k]);
         p.x[i][j][k] = rres.x[i][j][k] + bet*(p.x[i][j][k] - ome*v.x[i][j][k]);
         p.y[i][j][k] = rres.y[i][j][k] + bet*(p.y[i][j][k] - ome*v.y[i][j][k]);
         p.z[i][j][k] = rres.z[i][j][k] + bet*(p.z[i][j][k] - ome*v.z[i][j][k]);
   }}}
   Atimes(p0,p,v0,v,Coeff0,Coeffx,Coeffy,Coeffz,gij,gup,Kij,Kup,C2,Par);  //(d)
   dp = dotproduct(rhat0,rhat,v0,v,Par);       //(e)
   alp = rho[kay]/dp;
   for(i=0;i<Par.nxb;i++){                     //(f)
     for(j=0;j<Par.nyb;j++){
       for(k=0;k<Par.nzb;k++){
         rres0.s[i][j][k] = rres0.s[i][j][k] - alp*v0.s[i][j][k];
         rres.x[i][j][k] = rres.x[i][j][k] - alp*v.x[i][j][k];
         rres.y[i][j][k] = rres.y[i][j][k] - alp*v.y[i][j][k];
         rres.z[i][j][k] = rres.z[i][j][k] - alp*v.z[i][j][k];
   }}}
   Atimes(rres0,rres,t0,t,Coeff0,Coeffx,Coeffy,Coeffz,gij,gup,Kij,Kup,C2,Par);
   dp = dotproduct(t0,t,rres0,rres,Par);       //(g)
   dp2 = dotproduct(t0,t,t0,t,Par);
   ome = dp/dp2;
   dp = dotproduct(rhat0,rhat,t0,t,Par);
   rho[kay+1] = -ome*dp;
   for(i=0;i<Par.nxb;i++){                     //(h)
     for(j=0;j<Par.nyb;j++){
       for(k=0;k<Par.nzb;k++){
         Phi.s[i][j][k] = Phi.s[i][j][k] + alp*p0.s[i][j][k] + ome*rres0.s[i][j][k];
         Eks.x[i][j][k] = Eks.x[i][j][k] + alp*p.x[i][j][k] + ome*rres.x[i][j][k];
         Eks.y[i][j][k] = Eks.y[i][j][k] + alp*p.y[i][j][k] + ome*rres.y[i][j][k];
         Eks.z[i][j][k] = Eks.z[i][j][k] + alp*p.z[i][j][k] + ome*rres.z[i][j][k];
   }}}
   for(i=0;i<Par.nxb;i++){                     //(i)
     for(j=0;j<Par.nyb;j++){
       for(k=0;k<Par.nzb;k++){
         rres0.s[i][j][k] = rres0.s[i][j][k] - ome*t0.s[i][j][k];
         rres.x[i][j][k] = rres.x[i][j][k] - ome*t.x[i][j][k];
         rres.y[i][j][k] = rres.y[i][j][k] - ome*t.y[i][j][k];
         rres.z[i][j][k] = rres.z[i][j][k] - ome*t.z[i][j][k];
   }}}
   dp = dotproduct(rres0,rres,rres0,rres,Par);
   printf("Residual at iteration %d equals %12.8f \n",kay,dp);
 }



  free_dVector3D(rres0.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(rhat0.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(p0.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(v0.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(t0.s,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(rres.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(rhat.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(p.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(v.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(t.x,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(rres.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(rhat.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(p.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(v.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(t.y,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(rres.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(rhat.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(p.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(v.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);
  free_dVector3D(t.z,0,Par.nxb,0,Par.nyb,0,Par.nzb);

return;
}
