#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_fft_complex.h>  /* da fare: /opt/homebrew/Cellar/gsl/2.6/include/ in includePath */

#define N 4096
#define DIM_rho 10  /* number of densities from rho0 to rho desired */ 

#define POS_rhoDes 0

double rho[DIM_rho], drho;
double rho_desid[4] = {0.1, 0.3, 0.5, 0.7};
double T = 1.35; /* fixed temperature */

double L = 5.0; /* cutoff */
double r[N], dr;

double alpha = 20.0;
double a = 0.1; /* gentle factor */

double pot[N], der_pot[N];

double g[N]; /* radial distribution function */
double h0[N], h[N]; /* pair correlation function */
double c0[N], c[N]; /* direct correlation function*/
double F[N]; /* integral in OZ equation */

double cbar_k[N], hbar_k[N], fbar_k[N]; /* fourier radial transforms */
double diff_h[N], hh0;
double P, meanU; /* pressure and u mean */

double f_function(double x){
    /* expresses the reduced potential */
    double pt1, pt2;
    if (x <= 0.5)
    {
        pt1 = 1.0 / (1.0 - 6.0 * pow(alpha, -1.0));
        pt2 = ( 6.0 / alpha * exp(alpha * 0.5) - pow(0.5, -6.0) );
        return pt1 * pt2;
    } else {
        pt1 = 1.0 / (1.0 - 6.0 * pow(alpha, -1.0));
        pt2 = ( 6.0 / alpha * exp(alpha - alpha * x) - pow(x, -6.0) );
        return pt1 * pt2;
    }
}

double g_function(double x){
    /* limit at low density of g */
    return exp(- 1.0/T * f_function(x) );
}

double der_f_function(double x){
    /* performs numerical derivative of potential */
    return (f_function(x + dr) - f_function(x-dr))/(2 * dr);
}

double integrand_Pex(double x){
    return der_f_function(x) * x * x * x;
}

double integrand_meanU(double x){
    return f_function(x) * 4.0 * M_PI * x * x;
}


void radial_FFT_forward( int dim, double fun[dim], double c_k[dim] ){
    /* performs the radial fft forward of the function given in input fun */
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

void radial_FFT_inverse( int dim, double fun[dim], double c_n[dim]){
    /* performs the radial fft inverse of the function given in input fun */
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
    /* perfroms integration with trapezoids method */
    double h, f_a, f_b;
    double trapez_sum= 0;
    h = fabs(cutoff - a)/(N-1.0); 
    f_a = f[0];
    f_b = f[N-1];

 for (int jj = 1; jj < N; jj++)
 {
       trapez_sum += f[jj];
 }
  trapez_sum = trapez_sum + 0.5 * (f_a + f_b);
return h * trapez_sum;
}

double gIntegral(double g[N], double (*fun) (double),  double a, double cutoff ){
    /* performs integration with trapezoids method with g function */
    double h, f_a, f_b;
    double trapez_sum= 0;
    h = fabs(cutoff - a) /(N - 1.0);

    f_a = fun(r[0]) * g[0];
    f_b = fun(r[N-1]) * g[N-1];

    for (int jj = 1; jj < N; jj++)
      {
        trapez_sum += fun(r[jj]) * g[jj];
      }
    trapez_sum = trapez_sum + 0.5 * (f_a+f_b);

    return h * trapez_sum;
}

//------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {

FILE *pf_h0, *pf_g;
pf_h0 = fopen("data_h0.txt", "w");
pf_g = fopen("data_g.txt", "w");

rho[0] = 0.05; /* initial value of density */
rho[DIM_rho - 1] = rho_desid[POS_rhoDes]; /* desidered value of density */
drho = (rho[DIM_rho - 1] - rho[0])/ (DIM_rho - 1.0);
for (int ii = 1; ii < DIM_rho; ii++)
{
   rho[ii] = rho[ii-1] + drho;
}

/* initial conditions for low density */
dr = L/N;
fprintf(pf_h0, "#r \t\t v(r) \t\t der_v(r) \t\t h0(r)\n");
for (int ii = 0; ii < N; ii++)
{
    r[ii] = ii * dr;
    pot[ii] = f_function(r[ii]); 
    der_pot[ii] = der_f_function(r[ii]);
    h0[ii] = g_function(r[ii]) - 1.0; /* initial esteem of h */
    c0[ii] = h0[ii]; /* initial esteem of c */
    fprintf(pf_h0, "%lf\t %lf\t %lf\t %lf\n", r[ii], pot[ii], der_pot[ii], h0[ii]);
    h[ii] = h0[ii];
    c[ii] = c0[ii];
}
fclose(pf_h0);

for (int dd = 0; dd < DIM_rho; dd++) /* cycle on densities from rho0 to the desidered one */
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
            c[ii] = (g_function(r[ii]) - 1.0) * exp(F[ii]) + (exp(F[ii]) - 1.0 - F[ii]); // c_i+1
            h[ii] = c[ii] + F[ii]; // h_i+1
        }
        for (int ii = 0; ii < N; ii++)
        {
            c[ii] = a * c[ii] + (1.0 - a) * c0[ii];
            h[ii] = a * h[ii] + (1.0 - a) * h0[ii];
            g[ii] = g_function(r[ii]) * exp(h[ii] - c[ii]); // g[ii] = h[ii] + 1.0; (same result with both expressions )
            diff_h[ii] = fabs(h[ii] - h0[ii]);
        }
    hh0 = integral(diff_h, 1e-6, L);
    } while (hh0 > 1e-6);  /* cycle until I reach the desidered resliution */
}

fprintf(pf_g,"#r\t\t\t c\t\t\t h\t\t\t g\n");
for (int ii = 0; ii < N; ii++)
{
    fprintf(pf_g, "%lf\t %lf\t %lf\t %lf\n", r[ii], c[ii], h[ii], g[ii]);
} 
fclose(pf_g);

/* Evaluation of excess pressure and mean potential energy */
P = gIntegral(g, integrand_Pex, 1e-6, L );
P = - 2.0* M_PI / 3.0* rho[DIM_rho-1] * rho[DIM_rho-1] * P;

meanU = gIntegral(g,integrand_meanU, 1e-6, L);
meanU = 1.0/2.0 * rho[DIM_rho-1] * meanU;

printf("\n\n LJ potential\n T: %lf\t rho: %lf\t pressure: %lf \t meanU: %lf\n", T, rho[DIM_rho-1], P, meanU);

}