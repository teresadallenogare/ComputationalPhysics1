#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#define N 500
#define M 80
#define xi 1

double h;
double dK; 
double x[N];
double K[M];  
double complex H[N][N]; 
double w[N]; 
double E[M][3]; 

double potential (double x) {
    double esp = pow(x-0.5,2) / (2*0.3*0.3); 
    return -exp(-esp); 
}
void zheev_(char *jobz, char *uplo, int *n, complex double *m, int *lda, double *w, double complex *work, int *lwork, double *r, int *info);
int diagonalize_hermitean_lapack(double complex *m, double *w, int n) {
    char JOBZ='V'; /* eigenvalues and eigenvectors */
    char UPLO='L'; /* lower triangular FORTRAN = upper triangular C */
    int LDA = n; 
    double *RWORK;
    double complex *WORK;  
    int LWORK;
    int INFO; 
    double complex tmp = 0.0;
    double tmp1 = 0.0; 

    /* get the optimal work array size */
    LWORK=-1;
    zheev_(&JOBZ,&UPLO,&n,m,&LDA, w,&tmp,&LWORK,&tmp1,&INFO);
    
    LWORK = (int)(tmp);
    WORK = (complex double *)malloc(LWORK * sizeof(double complex));
    RWORK = (double *)malloc((3*n-2) * sizeof(double)); 
    
    /* call to the LAPACK subroutine */
    zheev_(&JOBZ,&UPLO,&n,m,&LDA, w, WORK,&LWORK,RWORK,&INFO);
    
    free(WORK);
    free(RWORK); 
    
    return(INFO);
}
void fileData (double K[M], double E[M][3]) {
    FILE *fd;     
    fd = fopen("dataE_K.txt","w"); 
    for (int kk=0; kk<M; kk++) {
        fprintf (fd, "%lf \t ", K[kk]); 
        for (int jj=0; jj<3; jj++) {
            fprintf (fd, "%lf \t", E[kk][jj]); 
        }
        fprintf (fd, "\n"); 
    }
  fclose(fd);
} 


int main(int argc, const char * argv[]){ 
    h = 1.0/N; 
    for (int ii=0; ii<N; ii++) {
        x[ii] = ii*h; 
    }
    dK = (2*M_PI)/(M-1);
    for (int jj=0; jj<M; jj++) {
        K[jj] = - M_PI + dK*jj; 
    }

    int info;
    for (int ii=0; ii<M; ii++) {
         memset(H, 0, N*N*sizeof(double complex));
    for (int kk=0; kk<N; kk++) {
        H[kk][kk] = potential(x[kk]) + 2.0*N*N; 
    }
    for (int kk=0; kk<N-1; kk++) {
        H[kk][kk+1] = -1.0*N*N; 
    }
    H[0][N-1] = -N*N*cexp(-I*K[ii]); 
        info = diagonalize_hermitean_lapack(&H[0][0],w,N);
        for (int jj=0; jj<3; jj++) {
            E[ii][jj] = w[jj]; 
        } 
    }
    fileData(K,E); 
}