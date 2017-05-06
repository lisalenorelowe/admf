#include "DeclareFunctions.h"

int InitializePhiandEks(Scalar *Phi, Vector *Eks, Params Par){
int i,j,k;

 Phi->s     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Eks->x     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Eks->y     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);
 Eks->z     = dVector3D(0,Par.nxb,0,Par.nyb,0,Par.nzb);

//Initialize 

 for(i=0;i<Par.nxb;i++){
  for(j=0;j<Par.nyb;j++){
   for(k=0;k<Par.nzb;k++){
    Phi->s[i][j][k] = 0.0;
    Eks->x[i][j][k] = 0.0;
    Eks->y[i][j][k] = 0.0;
    Eks->z[i][j][k] = 0.0;
 }}}


//End Initialize

return 0;
}
