#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_sf.h>

#define N 20000 // number of positions
#define M 800 // number of energies
#define L 25 // l = 0, 1, ..., 24
#define xi 0.024595 
#define eps 4.688 // [meV]
#define sig 3.075e-10 // [m]

// starting parameters (in r. u.) and vectors
double r_min = 0.5;
double r_max = 8.0;
double r[N];
double E_min = 0.10666; 
double E_max = 1.0666;
double h, dE, delta; 
double E[M], Sigma[M]; // energies and cross section
double complex u[N]; // radial function

double potential (double r, int l) {
    // LJ potential in r. u. 
    return (xi*l*(l+1)) / (r*r) + 4*( pow(r,-12) - pow(r,-6) ); 
}
double F (double r, double E, int l) {
    // F function for Numerov's method
    return (1/xi) * (potential(r,l) - E); 
}
void Numerov (double r[N], double complex u[N], double E, int l) { 
    // performs one step of Numerov's algorithm for diven values of E and l
    double complex num1, num2, den; 
    for (int kk = 2; kk<N; kk++) {
        num1 = ( 2.0 + (5.0/6.0)*h*h*F(r[kk-1],E,l) ) * u[kk-1]; 
        num2 = (1 - (h*h/12.0)*F(r[kk-2],E,l) )  * u[kk-2]; 
        den = (1 - (h*h/12.0)*F(r[kk],E,l));  
        u[kk] = (num1 - num2)/den; 
    }
} 
void fileData (double E[M], double Sigma[M]) {
    // file .txt for plots 
    FILE *fd;     
    fd = fopen("dataCS.txt","w"); 
    for (int kk=0; kk<M; kk++) {
        fprintf(fd, "%e \t %e \n", E[kk] * eps, Sigma[kk] * sig * sig * 1e20); 
    }
    fclose(fd);
}

int main(int argc, const char * argv[]){ 

    // positions
    h = (r_max-r_min) / (N-1); 
    for (int jj=0; jj<N; jj++) {
        r[jj] = r_min + jj*h; 
    }

    // energies
    dE = (E_max-E_min) / (M-1); 
    for (int ii=0; ii<M; ii++) {
        E[ii] = E_min + ii*dE; 
    } 

    double r1 = r[N-10]; // point near +inf
    double r2 = r[N-5]; // point near +inf
    double complex phi1, phi2, K, num, den; 
    double k; 
    for (int ii=0; ii<M; ii++) {
        k = sqrt(E[ii]/xi); 
        u[0] = cexp(-sqrt(4.0/(25.0*xi)) * pow(r[0],-5)); 
        u[1] = cexp(-sqrt(4.0/(25.0*xi)) * pow(r[1],-5));
        Sigma[ii] = 0.0; 
        for (int ll=0; ll<L; ll++) {
            Numerov(r,u,E[ii],ll); 
            phi1 = u[N-10]/r1; 
            phi2 = u[N-5]/r2; 
            K = phi2/phi1; 
            num = gsl_sf_bessel_jl(ll,k*r2) - K * gsl_sf_bessel_jl(ll,k*r1); 
            den = gsl_sf_bessel_yl(ll,k*r2) - K * gsl_sf_bessel_yl(ll,k*r1); 
            delta = atan(num/den); 
            Sigma[ii] += (2*ll+1) * sin(delta) * sin(delta); 
        }
        Sigma[ii] = (4*M_PI*Sigma[ii]) / (k*k); 
    }
    fileData(E,Sigma); 

}