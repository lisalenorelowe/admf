#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utility.h"
#include "Parameters.h"
#include "Components.h"
#define XMAX 50 
#define YMAX 50 
#define ZMAX 50 

void AllocateStructs(Params *Par);

void Derivs_admg(double t,Scalar alpha,Vector shift,Tensor Kij,Tensor gij,Tensor *dgij,Params Par);

void Derivs_admK(double t,Scalar alpha,Vector shift,Tensor gij,Tensor gup,Tensor Kij,Scalar trK, Tensor *dKij,Params Par,Connection C2,Tensor Ricci);

void Derivs_admlapse(double t,Scalar alpha,Vector shift,Tensor gup,Tensor Kij,double ***dalpha,Params Par);

void FunctionDerivs(double ***T,Tensor *dT,Params Par);

void FunctionFirstDerivs(double ***T,Vector *dT,Params Par);

void GetRicci(Tensor gij,Tensor gup,Tensor *Ricci,Params Par,Connection C2);

void fourierfit(double a, double b, double phibar[], double phitilde[], int n,double phi[], Params Par);

void fourierder(double a,double b,double phibar[],double phitilde[],double psibar[],double psitilde[],int n);

void Rkstep_adm(double t,double h,Scalar alpha,Vector shift,Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,Params Par,Scalar Phi,Vector Eks,Scalar hc,Vector mc,Connection C2,Tensor Ricci);

void Rxx(double ***Riccixx,Tensor gij,Tensor d2gxx,Tensor d2gyy,Tensor d2gzz,Tensor d2gxy,Tensor d2gxz,Tensor d2gyz,Vector dgxx,Vector dgyy, Vector dgzz, Vector dgxy, Vector dgxz, Vector dgyz,Params Par);

void Ryy(double ***Ricciyy,Tensor gij,Tensor d2gxx,Tensor d2gyy,Tensor d2gzz,Tensor d2gxy,Tensor d2gxz,Tensor d2gyz,Vector dgxx,Vector dgyy, Vector dgzz, Vector dgxy, Vector dgxz, Vector dgyz,Params Par);

void Rzz(double ***Riccizz,Tensor gij,Tensor d2gxx,Tensor d2gyy,Tensor d2gzz,Tensor d2gxy,Tensor d2gxz,Tensor d2gyz,Vector dgxx,Vector dgyy, Vector dgzz, Vector dgxy, Vector dgxz, Vector dgyz,Params Par);

void Rxy(double ***Riccixy,Tensor gij,Tensor d2gxx,Tensor d2gyy,Tensor d2gzz,Tensor d2gxy,Tensor d2gxz,Tensor d2gyz,Vector dgxx,Vector dgyy, Vector dgzz, Vector dgxy, Vector dgxz, Vector dgyz,Params Par);

void Rxz(double ***Riccixz,Tensor gij,Tensor d2gxx,Tensor d2gyy,Tensor d2gzz,Tensor d2gxy,Tensor d2gxz,Tensor d2gyz,Vector dgxx,Vector dgyy, Vector dgzz, Vector dgxy, Vector dgxz, Vector dgyz,Params Par);

void Ryz(double ***Ricciyz,Tensor gij,Tensor d2gxx,Tensor d2gyy,Tensor d2gzz,Tensor d2gxy,Tensor d2gxz,Tensor d2gyz,Vector dgxx,Vector dgyy, Vector dgzz, Vector dgxy, Vector dgxz, Vector dgyz,Params Par);

void Christoffel(Tensor gup,Vector dgxx,Vector dgyy,Vector dgzz,Vector dgxy,Vector dgxz,Vector dgyz,Params Par,Connection C2);

int ginv(Tensor gij,Tensor *gup,Params Par);

int InitializeMetric(Tensor *gij,Params Par);

int InitializeKij(Tensor *Kij,Tensor *gij,Params Par);

int InitializeShift(Vector *shift,Params Par);

int InitializeLapse(Scalar *alpha,Params Par);

int IntegrateStep_adm(Scalar alpha,Vector shift,Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,double t,Params Par,Scalar Phi, Vector Eks,Scalar hc,Vector mc, Connection C2,Tensor Ricci);

int OutputStep_gij(char outfile[],double time,Tensor gij,Tensor Kij,Params Par);

int OutputStep_T(char outfile[],double time,Tensor T,Params Par);

int OutputStep_S(char outfile[],double time,Scalar T,Params Par);

int OutputStep_D(char outfile[],double time, double T[][][], Params Par);

int readinput(int *r1,int *r2,int *r3,double *r4,double *r5,double *r6,double *r7,double *r8,double *r9,double *r10);

int ResetLapse(double t,Scalar alpha,Params Par);

int ReadRestartDump(Scalar *alpha,Vector *shift,Tensor *gij,Tensor *Kij,Params Par);

int WriteRestartDump(int timestep,Scalar alpha,Vector shift,Tensor gij,Tensor Kij,Params Par);

double HamCon(Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,Scalar hc,Params Par,double KupKij[][][], Connection C2, Tensor Ricci);

double MomCon(Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,Vector mc,Params Par,double DxK[][][], double DyK[][][], double DzK[][][],Connection C2, Tensor Ricci) ;

double fouriereval(double a, double b, double phibar[],double phitilde[],int m, int i,Params Par);

void Filter(Tensor gij,Tensor Kij,Params Par,int imod3);

int InitializePhiandEks(Scalar *Phi, Vector *Eks, Params Par);

void ConCon(Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,Scalar alpha,Vector shift,Scalar hc,Vector mc,double KupKij[][][],double DxK[][][],double DyK[][][],double DzK[][][],Scalar Phi,Vector Eks,Params Par, Connection C2,Tensor Ricci);

void Esolve(Scalar Phi,Vector Eks,double Coeff0[][][],double Coeffx[][][],double Coeffy[][][],double Coeffz[][][],double Lam0[][][],double Lamx[][][],double Lamy[][][],double Lamz[][][],Connection C2,Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,Params Par);

void Laplacian(Scalar Phi, Scalar LapPhi, Tensor gup, Connection C2, Params Par);

void Loperator(Vector Eks,Tensor LEks, Tensor gij, Tensor gup, Connection C2, Params Par);

double dotproduct(Scalar A0,Vector A,Scalar B0,Vector B,Params Par);

void Atimes(Scalar Phi,Vector Eks,Scalar Atimes0,Vector Atimes,double Coeff0[][][],double Coeffx[][][],double Coeffy[][][],double Coeffz[][][],Tensor gij,Tensor gup,Tensor Kij,Tensor Kup,Connection C2,Params Par);

void Deltaop(Tensor LEks,Vector DelEks, Tensor gij, Tensor gup, Connection C2, Params Par);

