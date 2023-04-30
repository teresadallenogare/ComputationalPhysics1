#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_fft_complex.h>  /* da fare: /opt/homebrew/Cellar/gsl/2.6/include/ in includePath */

#define N 2048
#define DIM_rho 20
#define POS_rho 0
#define POS_T 0

double rho[DIM_rho];
double drho;
double T = 1.3; // reduced temperature
double L = 5.0;
double h0[N], c0[N];
double h[N], c[N];
double F[N];
double r[N], k[N];
double dr, dk;
double P, meanU;

double fbar_k[N];
double alpha = 0.05;

double diff_h[N];
double hh0;

double g[N];
double LJ(double r){
  double x = (r<0.3 ? 0.3: r);
    return 4.0 * (pow (x, -12.0) - pow(x, -6.0));
}

double LJ_repulsive(double x){
    return 4.0 * pow(x, -12.0);
}
double derLJ(double r){
  double x = (r<0.3? 0.3: r);
  return 24.0 * (-2.0 * pow(x, -10.0) + pow(x, -4.0) );
}
double derLJ_repulsive(double x){
    return - 48.0 * pow(x, -10.0);
}
double integrand_meanU_LJ( double x){
 return LJ(x) * 4.0 * M_PI * x * x;
}

double integrand_meanU_LJ_repulsive(double x){
 return LJ_repulsive(x) * 4.0 * M_PI * x * x;
}

void setToZeroV(int dim, double v[dim]){
for (int ii = 0; ii < dim; ii++)
{
    v[ii] = 0.0;
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
for (int k = 0; k < dim; k++)
 {
     sum += (M_PI*k/L) * (M_PI*k/L) * fun[k];
 }
c_n[0] = 1.0/dim * sum;

gsl_fft_complex_radix2_forward((double *)t, 1, 2*dim);
for (int ii = 0; ii < 2*dim; ii++)
 {
     Tbar_n[ii] = -1.0/(2.0 * dim) * cimag(t[ii]); 
 }

for (int n = 1; n < dim; n++)
 {
    c_n[n] = dim/(n*L) * Tbar_n[n];

 } 
} 

double integral( double f[N], double a, double cutoff){
 
    double h,f_a, f_b;
    double trapez_sum= 0;
    h = cutoff/(N-1) ; // lunghezza intervallino
    f_a = f[0];
    f_b = f[N-1];

 for (int jj = 1; jj < N; jj++)
 {
       trapez_sum += f[jj];
   
 }
  trapez_sum = trapez_sum + 1.0/2.0 * (f_a+f_b);
return h* trapez_sum;
}

double gIntegral(double g[N], double (*fun) (double),  double a, double cutoff ){
    
    double h,  f_a, f_b;
    double trapez_sum= 0;
    h = cutoff/N ; // lunghezza intervallino

    f_a = fun(r[0])   * g[0];
    f_b = fun(r[N-1]) * g[N-1];

    for (int jj = 1; jj < N; jj++)
      {
        trapez_sum += fun(r[jj]) * g[jj];
      }
    trapez_sum = h * (trapez_sum + 0.5* (f_a+f_b));
  
    return trapez_sum;
}

//------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
    double cbar_k[N], hbar_k[N];
    FILE *pf_h, *pf_g;
    pf_h = fopen("data_h.txt", "w");
    pf_g = fopen("data_g.txt", "w");
    dr = L/N;
    dk = 2*M_PI/L;
    rho[0] = 0.05; // initial density
    rho[DIM_rho-1] = 0.3;
    drho = (rho[DIM_rho-1] - rho[0])/(DIM_rho-1.0); // step nel passare gradualmente da una densità bassa a una desiserata

    for (int ii = 1; ii < DIM_rho; ii++)
    {
        rho[ii] = rho[ii-1] + drho;
    }

    /*  --------------------------------------- LJ potential ------------------------------------------ */
    // stime iniziali che valgono per bassa densità
    fprintf(pf_h, "# r\t\t h\n");
    for (int ii = 0; ii < N; ii++)
    {
        r[ii] = (ii+1)* dr;
        k[ii] = dk * (ii<= N/2.0 ? ii : ii-N);
        h0[ii] = exp( - LJ(r[ii])/T ) - 1.0; 
        c0[ii] = h0[ii]; 
        fprintf(pf_h, "%lf\t%lf\n", r[ii], h0[ii]);
        c[ii] = c0[ii];
        h[ii] = h0[ii];
    }

for (int dd = 0; dd < DIM_rho ; dd++) // ripeto procedura per i diversi valori di rho che crescono gradualmente da rho0 alla rho desisderata
{
    do
    {
        radial_FFT_forward(N, c, cbar_k);
        radial_FFT_forward(N, h, hbar_k);

        for (int ii = 0; ii < N; ii++)
        {
            fbar_k[ii] = cbar_k[ii] * hbar_k[ii];
        }
        radial_FFT_inverse(N, fbar_k, F);
        for (int ii = 0; ii < N; ii++)
        {
            F[ii] = rho[dd] * 2.0 * M_PI*L/N * F[ii];
        }

        for (int ii = 0; ii < N; ii++)
        {
            c0[ii] = c[ii]; // c_i
            h0[ii] = h[ii]; // h_i
            c[ii] = (F[ii] + 1.0) * (exp( - LJ(r[ii])/T ) -1.0); // c_i+1
            h[ii] = c[ii] + F[ii]; // h_i+1
        }

        for (int ii = 0; ii < N; ii++)
        {
            c[ii] = alpha * c[ii] + (1.0 - alpha) * c0[ii];
            h[ii] = alpha * h[ii] + (1.0 - alpha) * h0[ii];
            // g[ii] = exp( - LJ(r[ii])/T ) * (1.0 + h[ii] - c[ii]);
            g[ii] = h[ii] + 1.0;
            diff_h[ii] = fabs(h[ii] - h0[ii]);
        }
        hh0 = integral(diff_h, 10e-6, L);
        // printf(" rho: %lf\t hh0: %lf\n", rho[dd], hh0);
    } while(hh0> 1e-6);

}
fprintf(pf_g,"#r\t\t\t c\t\t\t h\t\t\t g\n");
for (int ii = 0; ii < N; ii++)
{
    fprintf(pf_g, "%lf\t %lf\t %lf\t %lf\n", r[ii], c[ii], h[ii], g[ii]);
} 
// Evaluation of P and mean(U) 
P = gIntegral(g,derLJ, 1e-6, L );
P = -2.0* M_PI / 3.0* rho[DIM_rho-1] * rho[DIM_rho-1] * P;

meanU = gIntegral(g,integrand_meanU_LJ, 1e-6, L);
meanU = 1.0/2.0 * rho[DIM_rho-1] * meanU;

printf("\n\n LJ potential\n T: %lf\t rho: %lf\t pressure: %lf \t meanU: %lf\n", T, rho[DIM_rho-1], P, meanU);

fclose(pf_h);
fclose(pf_g);


/* ---------------------------------------- LJ repulsive potential ---------------------------------------- */
FILE *pf_gr, *pf_hr;

pf_hr = fopen("data_h_repulsive.txt", "w");
pf_gr = fopen("data_g_repulsive.txt", "w");

for (int ii = 0; ii < N; ii++)
    {
        r[ii] = (ii+1)* dr;
        k[ii] = dk * (ii<= N/2.0 ? ii : ii-N);
        h0[ii] = exp( - LJ_repulsive(r[ii])/T ) - 1.0; 
        c0[ii] = h0[ii]; 
        fprintf(pf_hr, "%lf\t%lf\n", r[ii], h0[ii]);
        c[ii] = c0[ii];
        h[ii] = h0[ii];
    }
int kk;
for (int dd = 0; dd < DIM_rho; dd++) // ripeto procedura per i diversi valoori di rho che crescono gradualmente da rho0 alla rho desisderata
{
    kk= 0;
   do
    {
        radial_FFT_forward(N, c, cbar_k);
        radial_FFT_forward(N, h, hbar_k);

        for (int ii = 0; ii < N; ii++)
        {
            fbar_k[ii] = cbar_k[ii] * hbar_k[ii];
        }
        radial_FFT_inverse(N, fbar_k, F);
        for (int ii = 0; ii < N; ii++)
        {
            F[ii] = rho[dd] * 2.0 * M_PI*L/N * F[ii];
        }

        for (int ii = 0; ii < N; ii++)
        {
            c0[ii] = c[ii]; // c_i
            h0[ii] = h[ii]; // h_i
            c[ii] = (F[ii] + 1.0) * (exp( - LJ_repulsive(r[ii])/T ) -1.0); // c_i+1
            h[ii] = c[ii] + F[ii]; // h_i+1
        }

        for (int ii = 0; ii < N; ii++)
        {
            c[ii] = alpha * c[ii] + (1.0 - alpha) * c0[ii];
            h[ii] = alpha * h[ii] + (1.0 - alpha) * h0[ii];
            // g[ii] = exp( - LJ_repulsive(r[ii])/T ) * (1.0 + h[ii] - c[ii]);
            g[ii] = h[ii] + 1.0;
            diff_h[ii] = fabs(h[ii] - h0[ii]);
        }
        kk++;
    hh0 = integral(diff_h, 10e-5, L);
    // printf(" rho: %lf\t hh0: %lf\n", rho[dd], hh0);
    }  while(hh0 > 1e-6  && kk<1000);
}

fprintf(pf_gr,"#r\t\t\t c\t\t\t h\t\t\t g\n");
for (int ii = 0; ii < N; ii++)
{
    fprintf(pf_gr, "%lf\t %lf\t %lf\t %lf\n", r[ii], c[ii], h[ii], g[ii]);
}
/* Evaluation of P and mean(U) */
P = gIntegral(g,derLJ, 1e-6, L );
P = - 2.0* M_PI / 3.0* rho[DIM_rho-1] * rho[DIM_rho-1] * P;

meanU = gIntegral(g,integrand_meanU_LJ, 1e-6, L);
meanU = 1.0/2.0 * rho[DIM_rho-1] * meanU;

printf("\n\n LJ potential\n T: %lf\t rho: %lf\t pressure: %lf \t meanU: %lf\n", T, rho[DIM_rho-1], P, meanU);

fclose(pf_hr);
fclose(pf_gr);
}
// -----------------------------------------------------------------------------------
