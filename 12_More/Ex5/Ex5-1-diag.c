#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#define N 2048 // THIS MUST BE EQUAL TO N IN Ex2-2.c!
#define xi 0.015

double L = 2.5; // -L is our -infinity
double x_initial; // starting point 
double h; // spatial step
double x[N], phi[N]; // position and wavefunction (real)

void dstedc_(char *compz, int *n, double *D, double *E, double *Z, int *ldz, double *work, int *lwork, int *iwork, int *liwork, int *info); 

int diagonalize_triangular_lapack(double *d, double *e, int n, double *vec) {
    // diagonalizes a trnagular matrix n x n with 
    // - d vector with the diagonal elements
    // - e is a vector with the subdiagonal elements
    // on outpul d contains the eigenvalues and vec contains the eigenvectors as rows

  char COMPZ='I';
  int LDZ = n;  
  double *WORK;
  double tmpwork;  /* optimal LWORK */
  int LWORK = -1;  /* query */
  int *IWORK;
  int tmpiwork;    /* optimal LIWORK */
  int LIWORK = -1; /* query */
  int INFO;

  // query 
  dstedc_(&COMPZ, &n, d, e, vec, &LDZ,          
          &tmpwork, &LWORK, &tmpiwork, &LIWORK, &INFO);

  LWORK = (int)tmpwork;
  WORK = (double *)malloc(LWORK * sizeof(double));
  LIWORK = tmpiwork;
  IWORK = (int *)malloc(LIWORK * sizeof(int));

  // diagonalize 
  dstedc_(&COMPZ, &n, d, e, vec, &LDZ,
          WORK, &LWORK, IWORK, &LIWORK, &INFO);
  
  
  free(WORK);
  free(IWORK);

  return(INFO);
}

double potential_L (double x) {
    return  -x * x * (x + 1.0);
}

double integrate (double f[N]){
    // integral using given values of a function f
    double ris = f[0] + f[N-1];
    for (int ii=1; ii<(N-1); ii++) {
        ris += 2*f[ii];
    }
    return ris * h/2.0;
}

int main(int argc, const char * argv[]){ 
    
    // initial conditions 
    x_initial = -L; 
    h = L/(N-1); 
    double *x = (double *)calloc(N, sizeof(double));  
    double *d = (double *)calloc(N, sizeof(double));
    double *e = (double *)calloc(N, sizeof(double));
    double *v = (double *)calloc(N*N, sizeof(double));
    for (int kk=0; kk<N; kk++) {
        x[kk] = x_initial + kk*h; 
        d[kk] = potential_L(x[kk]) + 2 * (xi/(h*h));
        e[kk] = -xi/(h*h);  
    }
    
    int inf = diagonalize_triangular_lapack(d, e, N, v); 

    // energy of the ground state
    double E0 = d[0]; 
    printf("Ground state energy: %lf \n", E0); 
    
    // normalization
    double f[N]; 
    for (int kk=0; kk<N; kk++) {
        f[kk] = fabs(v[kk]) * fabs(v[kk]); 
    }
    double mod = integrate(f); 
    for (int kk=0; kk<N; kk++) {
        phi[kk] = (1.0/sqrt(mod)) * v[kk]; 
    } 

    // ground normalized eigenfunction 
    FILE *fd;     
    fd = fopen("wf_gs.txt","w"); 
    for (int kk=0; kk<N; kk++) {
        fprintf(fd, "%lf \t %lf \n", x[kk], phi[kk]); 
    }
    fclose(fd);
   
}
