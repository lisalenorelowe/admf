#include <stdio.h>
#include <math.h>
#define NMAX 100

int main()
{
  double fields[NMAX][30][30][30];
  int i,j,k,n,loops;

  for(loops=0;loops<1000;loops++)
    {
      printf("loops = %d \n",loops);
      for(n=0;n<NMAX;n++)
	{
	  for(i=0;i<30;i++)
	    {
	      for(j=0;j<30;j++)
		{
		  for(k=0;k<30;k++)
		    {
		      fields[n][i][j][k] = cos(1.0*i*j*k*n);
		    }
		}
	    }
	}
    }

  printf("result = %f \n", fields[NMAX-1][29][29][29]);
  return 7;

}
