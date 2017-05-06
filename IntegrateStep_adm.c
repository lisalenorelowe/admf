#include "DeclareFunctions.h"

int IntegrateStep_adm(Scalar alpha,Vector shift,Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,double t,Params Par,Scalar Phi, Vector Eks,Scalar hc,Vector mc,Connection C2,Tensor Ricci)
{

  Rkstep_adm(t,Par.dt,alpha,shift,gij,gup,Kij,Kup,Par,Phi,Eks,hc,mc,C2,Ricci);

  return 0;
}  
