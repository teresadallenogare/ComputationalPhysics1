#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
  C wrapper to the LAPACK routine dsyev:
  double
  symmetric
  eigen
  values
  
  compile with whatever works:
  gcc -O2 program.c -lopenblas -lm
  or
  gcc -O2 program.c -llapack -lopenblas -lm
  or
  gcc -O2 -framework Accelerate program.c -lm (macOS)
  or
  gcc -O2 program.c -llapack -lblas -lm
*/

void dsyev_(char *jobz, char *uplo, int *n,
            double *m, int *lda, double *w, double *work,
            int *lwork, int *info);

int diagonalize_symmetric_lapack(double *m, double *w, int n)
{
  /*
    diagonalization of the SYMMETRIC matrix m (double precision).
    m[0..N-1][0..N-1];


    On output the matrix m has the eigenvectors as its ROWS
    Eigenvalues are returned in w.

    Inly the UPPER part of m need to be filled (in C convention, that is
    the first row, the second row apart from the first entry, and so on)

    returns INFO. diagonalization ok if INFO == 0.
  */

  char JOBZ='V'; /* eigenvalues and eigenvectors */
  char UPLO='L'; /* lower triangular FORTRAN = upper triangular C */
  int LDA = n;
  double *WORK;
  int LWORK;
  int INFO;
  double tmp = 0.0;

  /* get the optimal work array size */
  LWORK=-1;
  dsyev_(&JOBZ,&UPLO,&n,m,&LDA, w,&tmp,&LWORK,&INFO);
  
  LWORK = (int)(tmp);
  WORK = (double *)malloc(LWORK * sizeof(double));
  
  /* call to the LAPACK subroutine */
  dsyev_(&JOBZ,&UPLO,&n,m,&LDA, w,WORK,&LWORK,&INFO);
  
  free(WORK);

  return(INFO);
}

int main(int argc, char *argv[])
{
  int n = 3;
  double w[3];
  double m[3*3];
  int i,j;
  
  memset(m, 0, n*n*sizeof(double));

  /*
    Let's diagonalize the matrix
    
    1 2 3
    2 4 5
    3 5 6

    only triangular superior part
  */
  m[0] = 1.0;
  m[1] = 2.0;
  m[2] = 3.0;
  m[4] = 4.0;
  m[5] = 5.0;
  m[8] = 6.0;

  if(diagonalize_symmetric_lapack(m,w,n))
    {
      /* check for errors */
    }

  printf("autovalori: ");
  for(i=0;i<3;i++) printf("%g ",w[i]);
  printf("\n");
  
  for(i=0;i<3;i++)
    {
      printf("autovettore %d: ",i);
      for(j=0;j<3;j++) printf("%g ",m[3*i+j]);
      printf("\n");
    }
  
}