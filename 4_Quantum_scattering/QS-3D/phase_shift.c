#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_sf.h>

#define N 50000
#define M 200
#define xi 1

double h; 
double r0 = 1.0;  
double r_max = 10.0;
double dE = 0.01; 
double E[M], delta[M], delta_T[M]; 
double r[N]; 
double complex u[N]; 
int l;


double potentialAM (double r, int l) {
    return (l*(l+1)) / (r*r); 
}
double F (double r, double E, int l) {
    return potentialAM(r,l) - E/xi; 
}
void Numerov (double r[N], double complex u[N], double E) { 
    double complex num1, num2, den; 
    for (int kk = 2; kk<N; kk++) {
        num1 = ( 2.0 + (5.0/6.0)*h*h*F(r[kk-1],E,l) ) * u[kk-1]; 
        num2 = (1 - (h*h/12.0)*F(r[kk-2],E,l) )  * u[kk-2]; 
        den = (1 - (h*h/12.0)*F(r[kk],E,l));  
        u[kk] = (num1 - num2)/den; 
    }
} 
void fileDataPS (double E[M], double delta[M], double delta_T[M]) {
  FILE *fd;     
  fd = fopen("dataPS.txt","w"); 
  for (int kk=0; kk<M; kk++) {
      fprintf(fd, "%g \t %g \t %g \n", E[kk], delta[kk], delta_T[kk]); 
  }
  fclose(fd);
}

int main(int argc, const char * argv[]){ 
    h = (r_max-r0) / (N-1); 
    for (int jj=0; jj<N; jj++) {
        r[jj] = r0 + jj*h; 
    }
    for (int ii=0; ii<M; ii++) {
        E[ii] = (ii+1)*dE; 
    }
    l = 1; 
    double k; 
    double r1 = r[N-10]; 
    double r2 = r[N-5]; 
    double complex phi1, phi2, K, num, den; 
    for (int ii=0; ii<M; ii++) {
        k = sqrt(E[ii]/xi); 
        u[0] = 0.0; 
        u[1] = cexp(I*k*r[1]);  
        Numerov(r,u,E[ii]); 
        phi1 = u[N-10]/r1; 
        phi2 = u[N-5]/r2; 
        K = phi2/phi1; 
        num = gsl_sf_bessel_jl(l,k*r2) - K * gsl_sf_bessel_jl(l,k*r1); 
        den = gsl_sf_bessel_yl(l,k*r2) - K * gsl_sf_bessel_yl(l,k*r1); 
        delta[ii] = num/den; 
        delta[ii] = atan(delta[ii]); 
        num = pow(gsl_sf_bessel_jl(l,k),2); 
        den = pow(gsl_sf_bessel_jl(l,k),2) + pow(gsl_sf_bessel_yl(l,k),2);
        delta_T[ii] = sqrt(num/den); 
        delta_T[ii] = asin(-delta_T[ii]);  
    }
    fileDataPS(E,delta,delta_T); 
}
