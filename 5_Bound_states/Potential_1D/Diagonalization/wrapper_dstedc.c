#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* 
   C wrapper to the LAPACK routine dstedc:
   double
   symmetric
   triangular
   eigenproblem
   divide &
   conquer

   compile with whatever works:
   gcc -O2 program.c -lopenblas -lm
   or
   gcc -O2 program.c -llapack -lopenblas -lm
   or
   gcc -O2 -framework Accelerate program.c -lm (macOS)
   or
   gcc -O2 program.c -llapack -lblas -lm
   
*/
void dstedc_(char *compz, int *n, double *D, double *E,
             double *Z, int *ldz, double *work,
             int *lwork, int *iwork, int *liwork, int *info);

int diagonalize_triangular_lapack(double *d, double *e, int n, double *vec)
{
  /*
    of a double precision symmetric triangular nxn matrix
    On input d is a vector with diagonal elements, and e is a vector with
    the subdiagonal elements.

    On output:
    d contains the eigenvalues
    e is destroyed
    vec contains the eigenvectors as rows

    returns INFO. Diagonalization ok if INFO == 0
  */

  char COMPZ='I';
  int LDZ = n;  
  double *WORK;
  double tmpwork;  /* optimal LWORK */
  int LWORK = -1;  /* query */
  int *IWORK;
  int tmpiwork;    /* optimal LIWORK */
  int LIWORK = -1; /* query */
  int INFO;

  /* query */
  dstedc_(&COMPZ, &n, d, e, vec, &LDZ,          
          &tmpwork, &LWORK, &tmpiwork, &LIWORK, &INFO);

  LWORK = (int)tmpwork;
  WORK = (double *)malloc(LWORK * sizeof(double));
  LIWORK = tmpiwork;
  IWORK = (int *)malloc(LIWORK * sizeof(int));

  /* diagonalize */
  dstedc_(&COMPZ, &n, d, e, vec, &LDZ,
          WORK, &LWORK, IWORK, &LIWORK, &INFO);
  
  
  free(WORK);
  free(IWORK);

  return(INFO);
}

int main(int argc, char *argv[])
{
  int n = 1024;
  double L = 5.0;
  double dx = 2.0*L/(n-1);
  double *x   = (double *)calloc(n, sizeof(double));  
  double *d   = (double *)calloc(n, sizeof(double));
  double *e   = (double *)calloc(n, sizeof(double));
  double *vec = (double *)calloc(n*n, sizeof(double));

  int i,j;

  /* prepare the vectors */
  for(i=0;i<n;i++)
    {
      x[i] = -L + i*dx;
      d[i] = 1.0/(dx*dx) + 0.5*(x[i]*x[i]);
      e[i] = -0.5/(dx*dx);
    }
    
  if(diagonalize_triangular_lapack(d,e,n,vec))
    {
      /* check for errors */
    }

  printf("# autovalori: ");
  for(i=0;i<10;i++) printf("%g ",d[i]);
  printf("\n");
  
  for(i=0;i<5;i++)
    {
      printf("# autovettore %d:\n",i);
      for(j=0;j<n;j++) printf("%g %g\n",x[j], vec[n*i+j]/sqrt(dx));
      printf("&\n");
    }
  
}