#include <stdio.h>
#include <math.h>

int readinput(int *r1,int *r2,int *r3,double *r4,double *r5,double *r6,double *r7,double *r8,double *r9,double *r10){
 FILE *finput;
 char blah[100];

 finput = fopen("input","r"); //Open input file
//Read Data:
  fscanf(finput,"%d %d %d\n",r1,r2,r3);  
//End Read Data.
  *r4 = 0.05;
  *r5 = -1.;
  *r6 = 1.; 
  *r7 = -1.;
  *r8 = 1.;
  *r9 = -1.;
  *r10 = 1.;

 fclose(finput); //close input file
 return 0;
}
