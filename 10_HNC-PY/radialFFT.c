#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define N 1024 /* number of points of divisision */
#define L 4.0

double h[N]; // cn

double r[N], k[N];
double dr, dk;
double fbar_k[N], tbar_n[N];
double S[N];
double rho = 0.5;
double T = 1.3;

double LJpotential(double r){
return 4.0 * (pow(r, -12.0) - pow(r, -6.0));
}

double gfunction(double r){
    if (r < dr )
    {
       return 1.0;
    } else {
    return exp(- 1.0/T * LJpotential(r));
    }
    
}


void radial_FFT_forward( int dim, double fun[dim], double c_k[dim] ){
    /* performs the radia fft forward of the function given in input f */
 // ck bar
    double Fbar_k[2*dim];
    double complex f[dim]; 

    double sum = 0.0;
    
for (int ii = 0; ii < N; ii++) {
        f[ii] = fun[ii] * ii * L/N;
    }
 for (int ii = 0; ii < 2* dim; ii++)
 {
     if (ii < dim)
    {
        f[ii] = f[ii];
    } else {
        if (ii == dim)
        {
            f[ii] = 0.0;
        }else{
        f[ii] = - f[2*dim-ii];
        }   
    }
 }

 for (int n = 0; n < dim; n++)
        {
             sum += (n * L/dim) *(n * L/dim) * f[n]; 
        }
    c_k[0] = 2.0 * sum;
    
 gsl_fft_complex_radix2_forward((double *)f, 1, 2*dim);
 for (int ii = 0; ii < 2*dim; ii++) 
 {
     Fbar_k[ii] = - cimag(f[ii]); 
 }
 for (int k = 1; k < dim; k++) 
 {
    c_k[k] = L/(M_PI * k) * Fbar_k[k];
 }

}


void radial_FFT_inverse(int dim, double fun[dim], double c_n[dim]){
 double Tbar_n[2*dim];
 double complex t[dim];
 double sum = 0.0;
    for (int ii = 0; ii<N; ii++) {
        t[ii] = fun[ii] * M_PI * ii/L;
    }
// cn bar
for (int ii = 0; ii < 2* dim; ii++)
 {
     if (ii < dim)
    {
        t[ii] = t[ii];
    } else {
        if (ii == dim)
        {
            t[ii] = 0.0;
        }else{
        t[ii] = - t[2*dim-ii];
        }   
    }
 }
for (int k = 0; k < N; k++)
 {
     sum += (M_PI*k/L) * (M_PI*k/L) * fun[k];
 }
c_n[0] = 1.0/N * sum;

  gsl_fft_complex_radix2_forward((double *)t, 1, 2*dim);
for (int ii = 0; ii < 2*dim; ii++)
 {
     Tbar_n[ii] = -1.0/(2.0 * N) * cimag(t[ii]); 
 }

for (int n = 1; n < dim; n++)
 {
    c_n[n] = N/(n*L) * Tbar_n[n];

 } 
} 

//------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {

FILE *pf1, *pf2, *pf3;
pf1 = fopen("data_forward.txt", "w");
pf2 = fopen("data_inverse.txt", "w");
pf3 = fopen("data_S.txt", "w");
dr = L/N;
dk = M_PI/L;
for (int ii = 0; ii < N; ii++) 
 {
    r[ii] = ii * dr;
    k[ii] = dk * (ii<= N/2 ? ii: ii-N); 
    h[ii] = gfunction(r[ii]) - 1.0;
 }
 
radial_FFT_forward(N, h, fbar_k);

fprintf(pf1, "# k\t\t fbar_k\n");
for (int ii = 0; ii < N; ii++)
{
    fprintf(pf1, "%lf\t %lf\n", k[ii], fbar_k[ii]);
}

for (int ii = 0; ii < N; ii++)
{
    // S[ii] = 1.0 + rho * fbar_k[ii];
    S[ii] = fbar_k[ii];

}
fprintf(pf3, "# k\t\t S(k)\n");
for (int ii = 0; ii < N; ii++)
{
    fprintf(pf3,"%lf\t %lf\n", k[ii], S[ii] );
}

for (int ii = 0; ii < N; ii++)
{
    k[ii] = dk * (ii<= N/2 ? ii: ii-N);    

}

radial_FFT_inverse(N, fbar_k, tbar_n);
fprintf(pf2, "# r\t\t tbar_n\n");
for (int ii = 0; ii < N; ii++)
{
    fprintf(pf2, "%lf\t %lf\n", r[ii], tbar_n[ii]);
}

fclose(pf1);
fclose(pf2);
fclose(pf3);
}
