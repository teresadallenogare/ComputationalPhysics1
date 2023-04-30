// to compile: gcc -O2 -o ris derivative.c -lgsl -lgslcblas -lm
// to execute: ./ris

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define N 1024
double complex f[N], tr[N], v[N], der[N], Delta; 
double x[N], k[N];  
double L = 100.0; 

double complex DeltaK (int k) {
    if (k<N/2.0) {
        return I*(2*M_PI*k)/L; 
    } else {
        return I*(2*M_PI*(k-N))/L;
    }
}
double complex func (double x) {
    return cos(x) * cos(x);
}
double complex der_an (double x) {
    return 2.0 * (-cos(x) * cos(x) + sin(x) * sin(x)); 
}

int main(int argc, const char * argv[]){ 

    double dx = L/N;
    for (int jj=0; jj<N; jj++) {
        x[jj] = -L/2.0 + jj * dx; 
        f[jj] = func(x[jj]);
    }
    double dk = (2*M_PI)/L;  
    for (int jj=0; jj<N; jj++) {
        k[jj] = dk * (jj<N/2.0 ?  jj : jj-N); 
    }

    gsl_fft_complex_radix2_forward((double *)f, 1, N);

    for (int jj=0; jj<N; jj++) {
        tr[jj] = f[jj];
    }

    for (int jj=0; jj<N; jj++) {
        Delta = DeltaK(jj) * DeltaK(jj);  
        v[jj] = Delta * tr[jj];  
    }

    gsl_fft_complex_radix2_inverse((double *)v, 1, N);
    
    for (int jj=0; jj<N; jj++) {
        f[jj] = func(x[jj]);
        der[jj] = der_an(x[jj]); 
    }
    
    FILE *fd;     
    fd = fopen("data.txt","w"); 
    for (int kk=0; kk<N; kk++) {
        fprintf(fd, "%lf \t %lf \t %lf \t %lf \n", x[kk], creal(f[kk]), creal(der[kk]), creal(v[kk])); 
    } 
    fclose(fd);


}