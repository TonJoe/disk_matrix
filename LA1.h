#define NR_END 1
#define LOG2 0.693147180559945
#include <complex>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#define TINY    1.5e-300
#define TINYL    1.5e-1000L

#define FREE_ARG char*
using namespace std;

void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run_time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}
void cpxdbl_ludcmp0(complex<double> *a,int n,double *d)
{
	double veevee[n+1];
	int i,imax,j,k;
	double big,dum,temp;
	complex<double> sum,dum2;

  *d=1.0;
  for(i=0;i<n;i++)
  {
    big=0.0;
    for(j=0;j<n;j++)
      if( (temp=abs(a[i*n+j]))>big)       big=temp;
    if(big==0.0)      	{
				nrerror("Singular matrix in routine cpx_ludcmp\n");
			}
    veevee[i]=1.0/big;      /* Save the scaling */
  }
  for(j=0;j<n;j++){
    for(i=0;i<j;i++){
      sum= a[i*n+j];
      for(k=0;k<i;k++)     sum=sum-a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
    }
    big=0.0;
    for(i=j;i<n;i++){
      sum=a[i*n+j];
      for(k=0;k<j;k++)
        sum=sum-a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;

      /* Is the figure of matrix for the pivot better than the best so far? */
      if( (dum=veevee[i]*abs(sum)) >= big){
        big=dum;
        imax=i;
      }
    }

    /*Do we need to interchange rows? Yes, do so.. */
    if(j!=imax){
      for(k=0;k<n;k++){
        dum2=a[imax*n+k];
        a[imax*n+k]=a[j*n+k];
  a[imax*n+k]=a[j*n+k];
        a[j*n+k]=dum2;
      }
      *d=-(*d);
      veevee[imax]=veevee[j];
    }
    if(abs(a[j*n+j])==0.0)     a[j*n+j]=TINY+polar(0.0,1.0)*TINY;
    if(j!=n){
      dum2=1.0/a[j*n+j];
      for(i=j+1;i<n;i++)      a[i*n+j]=a[i*n+j]*dum2;
    }
  }
}

complex<double> cpxdbl_det0(complex<double> *a,int n)
{
  int i;
  double d;
  complex<double> det;

  cpxdbl_ludcmp0(a,n,&d);
  det=d;
  for(i=0;i<n;i++){
    det=det*a[i*n+i];
  }
return det;
}