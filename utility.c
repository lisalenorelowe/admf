#include "utility.h"

void nerror(char error_text[])
{
  fprintf(stderr, "Run time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
 
}

/* routine to allocate memory for 1d and 2d arrays */

int *iVector1D(int low, int high)
{
  /* allocate an integer vector with subscript [low,high] */ 
  int *p;
  p = (int *)malloc((unsigned) ((high-low+1)*sizeof(int)));
  if(!p) nerror("Allocation failure in iVector1D()");
  return p-low;
}

int **iVector2D(int row_low, int row_high, int column_low, int column_high)
{
  /*** 
     allocate an two-domensional array with subscript 
     [row_low,row_high][column_low, column_high] 
  ***/
  
  int i, nrow=row_high-row_low+1, ncolumn=column_high-column_low+1;
  int **p;

  /* allocate pointers to rows */
  p = (int **) malloc((unsigned)((nrow)*sizeof(int *)));
  if (!p) nerror("allocation failure 1 in iVector2D()");
  p -= row_low;

  /* allocate rows and set pointers to them */
  for(i=row_low;i<=row_high;i++) {
    p[i]=(int *)malloc((unsigned) (ncolumn)*sizeof(int));
    if (!p[i]) nerror("allocation failure 2 in iVector2D()");
    p[i] -= column_low;
  }
  return p;
}

int ***iVector3D(int nrl, int nrh, int ncl, int nch, int nclrl, int nclrh)
{
  int     i, j;
  int  ***m;

  m=(int ***) malloc((unsigned) (nrh-nrl+1)*sizeof(int**));
  if (!m) nerror("allocation failure 1 in iVector3D()");
  m -= nrl;

  for(i=nrl; i<=nrh; i++) {
    m[i]=(int **) malloc((unsigned) (nch-ncl+1)*sizeof(int*));
    if (!m[i]) nerror("allocation failure 2 in iVector3D()");
    m[i] -= ncl;
  }

  for(i=nrl; i<=nrh; i++)
    for(j=ncl; j<=nch; j++) {
      m[i][j]=(int *) malloc((unsigned) (nclrh-nclrl+1)*sizeof(int));
      if (!m[i][j]) nerror("allocation failure 3 in iVector3D()");
      m[i][j] -= nclrl;
    }
  return m;
}

double *dVector1D(int nl, int nh)
{
  /* allocate a one dimensional double array with subscript [nl,nh] */
  double *v;
  
  v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  if (!v) nerror("allocation failure in dVector1D()");
  return v-nl;
}

double **dVector2D(int nrl,int nrh,int ncl,int nch)
{
  /* allocate a two dimensional double array with subscript [nrl,nrh][ncl,nch] */
  int i;
  double **m;
  
  m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
  if (!m) nerror("allocation failure 1 in dVector2D()");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
    if (!m[i]) nerror("allocation failure 2 in dVector2D()");
    m[i] -= ncl;
  }
  return m;
}

double ***dVector3D(int nrl, int nrh, int ncl, int nch, int nclrl, int nclrh)
{
  int     i, j;
  double  ***m;

  m=(double ***) malloc((unsigned) (nrh-nrl+1)*sizeof(double**));
  if (!m) nerror("allocation failure 1 in dVector3D()");
  m -= nrl;

  for(i=nrl; i<=nrh; i++) {
    m[i]=(double **) malloc((unsigned) (nch-ncl+1)*sizeof(double*));
    if (!m[i]) nerror("allocation failure 2 in dVector3D()");
    m[i] -= ncl;
  }

  for(i=nrl; i<=nrh; i++)
    for(j=ncl; j<=nch; j++) {
      m[i][j]=(double *) malloc((unsigned) (nclrh-nclrl+1)*sizeof(double));
      if (!m[i][j]) nerror("allocation failure 3 in dVector3D()");
      m[i][j] -= nclrl;
    }
  return m;
}



/*Vec *VecArray(int low, int high)
{
  ccc allocate an Vec array [low, high] ccc
  Vec *p;
  p = (Vec *)malloc((unsigned)((high-low+1)*sizeof(Vec)));
  if(!p) nerror("Allocation failure in iVector1D()");
  return p-low;
}*/

void free_iVector1D(int *v, int nl, int nh)
{
	free((char*) (v+nl));
}

void free_iVector2D(int **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_iVector3D(int ***m, int nrl, int nrh, int ncl, int nch, int nclrl, int nclrh)
{
  int i, j;

  for (i=nrh; i>=nrl; i--)
    for (j=nch; j>=ncl; j--)
      free((char*) (m[i][j]+nclrl));

  for (i=nrh; i>=nrl; i--) free((char*) (m[i]+ncl));

  free((char*) (m+nrl));
}

void free_dVector1D(double *v, int nl, int nh)
{
	free((char*) (v+nl));
}


void free_dVector2D(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dVector3D(double ***m, int nrl, int nrh, int ncl, int nch, int nclrl, int nclrh)
{
  int i, j;

  for (i=nrh; i>=nrl; i--)
    for (j=nch; j>=ncl; j--)
      free((char*) (m[i][j]+nclrl));

  for (i=nrh; i>=nrl; i--) free((char*) (m[i]+ncl));

  free((char*) (m+nrl));
}



/*void free_VecArray(Vec *p, int low, int high)
{
  int i, ierr;

  for (i=low; i<=high; i++) {
    ierr = VecDestroy(p[i]); CHKERRA(ierr);
  }
}*/

    



