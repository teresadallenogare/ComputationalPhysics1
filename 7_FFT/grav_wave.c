// to compile: gcc -O2 -o ris grav_wave.c -lgsl -lgslcblas -lm
// to execute: ./ris

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define N 8192
#define g 9.81
#define M 1000
#define EV_TIME 1000    // [s]


double complex f[N], A[N], Ampl[N]; 
double x[N], k[N], w[N], Ev[N];
double L = 100.0;        // period 
double h = 10.0;  

double dispersion_relation (double k, double h) {
    return sqrt(g * k * tanh(k*h)); 
}
void fileData (double x[N], double Ev[N]) {
    // generates file.txt
    FILE *fd;    
    fd = fopen("data_Ev.txt","w");
    for (int kk=0; kk<N; kk++) {
        fprintf(fd, "%lf \t %lf \n", x[kk], Ev[kk]); 
    }
    fclose(fd);
}
int maximum(double complex A[N]) {
    // finds the position of the maximum 
    int pos = 0; 
    double ris = creal(A[0]); 
    for (int kk = 1; kk < N; kk++) {
        if (creal(A[kk]) != 0.00000 && creal(A[kk]) > ris) {
            ris = creal(A[kk]); 
            pos = kk; 
        }
    }
    return pos; 
}

int main(int argc, const char * argv[]){

    // initial data (wavepacket and dispersion relation)
    double dx = L/N; 
    double x0 = L/2.0; 
    double dk = (2*M_PI)/L; 
    for(int ii=0; ii<N; ii++) {
        x[ii] = ii * dx;        // x in [0; L)
        f[ii] = exp(-pow((x[ii]-x0),2.0) / 8.0) * cexp(2*M_PI*I*(x[ii]-x0));
        k[ii] = dk * (ii <= N/2 ? ii : ii-N); 
        w[ii] = dispersion_relation(k[ii],h); 
        // k[ii] = ii * dk;        // k in [0; 2pi N/L )
    }

    // plot of the initial data (wavepacket in direct space and dispersion relation in reciprocal space) 
    FILE *fd1;     
    fd1 = fopen("data_initial.txt","w"); 
    for (int kk=0; kk<N; kk++) {
        fprintf(fd1, "%lf \t %lf \t %lf \t %lf \n", x[kk], creal(f[kk]), k[kk], w[kk]); 
    } 
    fclose(fd1);

    // function A(k) = spatial FT of f at time 0 
    gsl_fft_complex_radix2_forward ((double *)f, 1, N);
    for (int ii=0; ii<N; ii++) {
        A[ii] = f[ii]; 
    }
    /*
    FILE *fd2;     
    fd2 = fopen("data_A.txt","w"); 
    for (int kk=0; kk<N; kk++) {
        fprintf(fd2, "%lf \t %lf \n", k[kk], creal(A[kk])); 
    }
    fclose(fd2); 
    */

    // group velocity
    int pos_k_max = maximum(A); 
    double k_max = k[pos_k_max]; 
    // printf("K_max: %lf \n", k_max); 
    double v_gA = 0.5 * pow(dispersion_relation(k_max,h),-1.0) * g * ( tanh(k_max*h) + k_max*h*pow(cosh(k_max*h),-2.0) );       // analytic derivate
    double v_gN = ( dispersion_relation(k_max+dk,h) - dispersion_relation(k_max-dk,h) ) / (2*dk);       // numerical derivate 
    printf("Group velocity (analytic): %lf \n", v_gA);
    printf("Group velocity (numerical): %lf \n", v_gN);

    // evolution 
    double final_time = EV_TIME/4.0;    // [s] (request: final_time = EV_TIME = 1000 s)
    double d_g = v_gA * final_time;      // group displacement 
    
    for (int ii=0; ii<N; ii++) {    // space cycle  
        Ampl[ii] = A[ii] * cexp(-I*w[ii]*final_time) * cexp(I*d_g*k[ii]);
    }
    gsl_fft_complex_radix2_inverse ((double *)Ampl, 1, N);
    for (int ii=0; ii<N; ii++) {
        Ev[ii] = creal(Ampl[ii]); 
    }
    fileData(x,Ev); 

} 