#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#define N 1000

double xi[3] = {0.025, 0.005, 0.0025}; 
double R = 5.0; 
double h; 
double x[N]; 
double H_k[N][N], H_p[N][N], H[N][N]; 
double w[N]; 

void set_zero (unsigned int r, unsigned int c, double v[r][c]) {
  for (int kk=0; kk<r; kk++) {
    for (int jj=0; jj<c; jj++) {
        v[kk][jj] = 0.0; 
    }
  }
}
void dsyev_ (char *jobz, char *uplo, int *n, double *m, int *lda, double *w, double *work, int *lwork, int *info);
int diagonalize_symmetric_lapack (double *m, double *w, int n) {
  /*
    diagonalization of the SYMMETRIC matrix m (double precision).
    m[0..N-1][0..N-1];


    On output the matrix m has the eigenvectors as its ROWS
    Eigenvalues are returned in w.

    Only the UPPER part of m need to be filled (in C convention, that is
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
double potential (double x) {
    return -pow(cosh(x),-4.0) ;
}
double eigenfunction (double x, double n) {
    double arg = (n*M_PI)/(2*R); 
    return (1/sqrt(R)) * sin(arg*(x+R));  
} 
double integrand (double x, int n1, int n2) {
    double p1 = eigenfunction(x,n1); 
    double p2 = eigenfunction(x,n2);
    return p1 * potential(x) * p2;  
}
double integrate (double (*f)(double, int, int), double a, double b, int n1, int n2){
    double ris = f(a,n1,n2) + f(b,n1,n2);
    for (int ii=1; ii<(N-1); ii++) {
        ris = ris + 2*f(x[ii],n1,n2);
    }
    return ris * h/2;
}

int main(int argc, const char * argv[]){

    /* EX 1 PART II: Bound states with diagonalization */

    h = (2*R)/(N-1); 
    for (int kk=0; kk<N; kk++) {
        x[kk] = -R + kk*h; 
    }

    memset(H, 0, N*N*sizeof(double));

    // kinetic contribute to H 
    memset(H_k, 0, N*N*sizeof(double));
    for (int ii=0; ii<N; ii++) {
        H_k[ii][ii] = (M_PI*M_PI*(ii+1)*(ii+1)) / (4*R*R); 
    }

    // potential contribute to H
    memset(H_p, 0, N*N*sizeof(double));
    for (int ii=0; ii<N; ii++) {
        for (int jj=ii; jj<N; jj++) {
            H_p[ii][jj] = integrate(&integrand, -R, R, ii+1, jj+1); 
        }
    }

    // complete H for xi[0]
    memset(H, 0, N*N*sizeof(double));
    double info1; 
    for (int ii=0; ii<N; ii++) {
        for (int jj=ii; jj<N; jj++) {
            H[ii][jj] = xi[0] * H_k[ii][jj] + H_p[ii][jj];  
        }
    }
    info1 = diagonalize_symmetric_lapack(&H[0][0],w,N); 
    printf("Eigenvalues for xi = %lf: \n", xi[0]);
    for(int jj=0; jj<N; jj++) {
        if( w[jj] < 0.0) {
            printf("%g \n",w[jj]);
        }
    }
    printf("\n");

    // complete H for xi[1]
    memset(H, 0, N*N*sizeof(double));
    double info2; 
    for (int ii=0; ii<N; ii++) {
        for (int jj=ii; jj<N; jj++) {
            H[ii][jj] = xi[1] * H_k[ii][jj] + H_p[ii][jj];  
        }
    }
    info2 = diagonalize_symmetric_lapack(&H[0][0],w,N); 
    printf("Eigenvalues for xi = %lf: \n", xi[1]);
    for(int jj=0; jj<N; jj++) {
        if( w[jj] < 0.0) {
            printf("%g \n",w[jj]);
        }
    }
    printf("\n");

    // complete H for xi[2]
    memset(H, 0, N*N*sizeof(double));
    double info3; 
    for (int ii=0; ii<N; ii++) {
        for (int jj=ii; jj<N; jj++) {
            H[ii][jj] = xi[2] * H_k[ii][jj] + H_p[ii][jj];  
        }
    }
    info3 = diagonalize_symmetric_lapack(&H[0][0],w,N); 
    printf("Eigenvalues for xi = %lf: \n", xi[2]);
    for(int jj=0; jj<N; jj++) {
        if( w[jj] < 0.0) {
            printf("%g \n",w[jj]);
        }
    }
    printf("\n");

}