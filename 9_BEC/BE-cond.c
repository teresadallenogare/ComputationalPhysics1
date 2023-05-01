#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h> 

#define M 1024

double xi[3] = {1.0, 10.0, 100.0}; 
double r[M], dr, x_max = 6.0; 
double u[M], Vsc[M], mu, phi[M][3]; 

double gaussian (double r) {
    return exp(-r*r); 
}

void SCpotential (double u[M], double r[M], double xi, double Vsc[M]) {
    for (int ii=0; ii<M; ii++) {
        Vsc[ii] = 2.0 * xi * fabs(u[ii]/r[ii]) * fabs(u[ii]/r[ii]); 
    }
}

void dstedc_(char *compz, int *n, double *D, double *E, double *Z, int *ldz, double *work, int *lwork, int *iwork, int *liwork, int *info);
int diagonalize_triangular(double *d, double *e, int n, double *vec) {
    // diagonalizes a triangular matrix n x n with 
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
double integrate (int dim, double f[dim]){
    // integral using given values of a function f
    double ris = f[0] + f[dim-1];
    for (int ii=1; ii<(dim-1); ii++) {
        ris += 2.0 * f[ii];
    }
    return ris * dr/2.0;
}

int main(int argc, const char * argv[]){ 

    // initial conditions 
    dr = x_max/(M-1); 
    double m[M];
    for (int ii=0; ii<M; ii++) {
        r[ii] = (ii+1) * dr;  
    }
    double mod; 

    // wave functions 
    double *d = (double *)calloc(M, sizeof(double));
    double *e = (double *)calloc(M, sizeof(double));
    double *v = (double *)calloc(M*M, sizeof(double));  

    double alpha = 0.1;  
    double mu_i, mu_f, mu[3]; 
    for (int jj=0; jj<3; jj++) {
        mu_i = 1.0;
        mu_f = 0.0; 
        for (int ii=0; ii<M; ii++) {
            u[ii] = sqrt(4.0 * M_PI) * r[ii] * gaussian(r[ii]); 
            m[ii] = u[ii] * u[ii]; 
        } 
        mod  = integrate(M,m); 
        for (int ii=0; ii<M; ii++) {
            u[ii] = 1.0/(sqrt(mod)) * u[ii]; 
        }
        while (fabs(mu_f - mu_i) > 1e-5) {
            SCpotential(u,r,xi[jj],Vsc);  
            mu_i = mu_f;
            for (int kk=0; kk<M; kk++) {
                d[kk] = r[kk] * r[kk] + Vsc[kk] + 2.0 /(dr*dr); 
            }
            for (int kk=0; kk<M; kk++) {
                e[kk] = -1.0/(dr*dr); 
            }
            int info = diagonalize_triangular(d,e,M,v); 
            mu_f = d[0]; 
            // printf("Mu = %lf \n ", mu_f);
            for (int kk=0; kk<M; kk++) {
                u[kk] = alpha * fabs(v[kk]) / sqrt(dr) + (1.0 - alpha) * u[kk]; 
                m[kk] = u[kk] * u[kk]; 
            }
            mod = integrate(M,m); 
            for (int kk=0; kk<M; kk++) {
                u[kk] = (1.0/sqrt(mod)) * u[kk]; 
                phi[kk][jj] = u[kk] / (sqrt(4.0 * M_PI) * r[kk]);  
            }
        }
        mu[jj] = mu_f; 
        // printf("Mu%d = %lf \n", jj+1, mu_f); 
    }
    FILE *fd;     
    fd = fopen("wf_gs.txt","w"); 
    for (int kk=0; kk<M; kk++) {
        fprintf(fd, "%lf \t ", r[kk]);
        for (int jj=0; jj<3; jj++) {
            fprintf(fd, "%lf \t", phi[kk][jj]);
        } 
        fprintf(fd, "\n"); 
    }
    fclose(fd);

    // energies 
    // 1) kinetic contribute
    double H_kin; 
    double der2u[M-2], u_kin[M-2]; 
    for (int ii=0; ii<(M-2); ii++) {
        der2u[ii] = (u[ii+2] + u[ii] - 2.0 * u[ii+1]) / (dr * dr); 
        u_kin[ii] = u[ii+1] * der2u[ii]; 
    }
    H_kin = -integrate(M-2,u_kin); 
    printf("Kinetic: %lf \n", H_kin); 

    // 2) harmonic contribute
    double H_har; 
    double u_har[M]; 
    for (int ii=0; ii<M; ii++) {
        u_har[ii] = u[ii] * r[ii] * u[ii]; 
    }
    H_har = integrate(M-2,u_har); 
    printf("Harmonic: %lf \n", H_har);

    // 3) self-consistent contribute
    double H_sc[3]; 
    double u_sc[M]; 
    for(int jj=0; jj<3; jj++) {
        for (int ii=0; ii<M; ii++) {
            u_sc[ii] = (pow(u[ii],4) / r[ii]) * (2.0 * xi[jj]); 
        }
        H_sc[jj] = integrate(M,u_sc); 
        printf("Self-consistent %d: %lf \n", jj, H_sc[jj]);
    }

}