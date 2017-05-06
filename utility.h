#ifndef UTILITY_H
#define UTILITY_H 
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <malloc.h>
/*#include "vec.h"*/

#define PI 3.141592653589793
 
/* memory allocation for 1d and 2d arrays */

void nerror(char error_text[]);

int *iVector1D(int, int);
int **iVector2D(int, int, int, int);
int ***iVector3D(int, int, int, int, int, int);
double *dVector1D(int, int);
double **dVector2D(int, int, int, int);
double ***dVector3D(int, int, int, int, int, int);

/*Vec *VecArray(int, int);*/

void free_iVector1D(int *, int, int);
void free_iVector2D(int **, int, int, int, int);
void free_iVector3D(int ***, int, int, int, int, int, int);
void free_dVector1D(double *, int, int);
void free_dVector2D(double **, int, int, int, int);
void free_dVector3D(double ***, int, int, int, int, int, int);

/*void free_VecArray(Vec *, int, int);*/

#endif





