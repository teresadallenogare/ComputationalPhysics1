#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_fft_complex.h>


#define N 1024 /* N, potenza di 2 */
#define DIM_t 200

#define g 9.81 //[m/s^2]

double h = 10.0; //[m] profondità
double lambda = 1.0; //[m]
double x[N];
double dx;
double L = 100.0; // lunghezza cella 

double k[N]; // vettori d'onda
double dk;

double w[N]; // omega

double t[DIM_t]; //[s] /(Ev_time = 1000.0 )
double dt;
double T = 1000.0;

double data[N][DIM_t];

double complex f[N], A[N], phi[N], Ampl[N];

double dispersionRelation(double k){
return sqrt(g * k * tanh(k*h));
}

double derDispersionRelation(double k){
    double pt1, pt2;
    pt1 = 0.5 * pow( g*k*tanh(h*k), -1.0/2.0);
    pt2 = g*tanh(h*k)+ g*k*h* 1.0/(cosh(h*k)*cosh(h*k));
    return pt1*pt2;
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


//------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
    FILE *pf_ev;
    pf_ev = fopen( "data_ev.txt", "w");

double pt1, x0;
double complex pt2;


dt = T/DIM_t;
for (int ii = 0; ii < DIM_t; ii++)
{
    t[ii] = ii*dt;
}

x0 = 0.0;
dx = L/N;
dk = 2.0*M_PI/L;
for (int ii = 0; ii < N; ii++)
{
    k[ii] =  dk * (ii<= N/2 ?ii : ii-N); 
    w[ii] = dispersionRelation(k[ii]); 
}

// servono quando intro velocità di gruppo
int pos_k_max = maximum(A); 
double k_max = k[pos_k_max];     
double v_gA = derDispersionRelation(k_max);    // analytic derivate
double d_g = v_gA * T;      // group displacement 


for (int kk = 0; kk < DIM_t; kk++) // ad ogni istante ho un pacchetto 
{
    for (int jj = 0; jj < N; jj++)
    {
        x[jj] = jj*dx - L/2.0;
        pt1 = exp( -(x[jj]-x0)*(x[jj]-x0)/(8.0*lambda*lambda) );
        pt2 = cexp( I*2*M_PI*(x[jj]-x0)/lambda );
        f[jj] = pt1 *pt2;
    }
    gsl_fft_complex_radix2_forward((double *)f, 1, N);
    for (int jj = 0; jj < N; jj++)
    {
        A[jj] = f[jj];
        Ampl[jj] = A[jj] * cexp(-I*w[jj]*t[kk]); //* cexp(I*d_g*k[jj]); // senza velocità di gruppo
    }
    gsl_fft_complex_radix2_inverse((double *)Ampl, 1, N);
    for (int jj = 0; jj < N; jj++)
    {
        data[jj][kk] = creal(Ampl[jj]);
    }
}

for (int kk = 0; kk < N; kk++)
{
    fprintf(pf_ev, "%lf\t", x[kk]);
    for (int jj = 0; jj < DIM_t; jj++)
    {
     fprintf(pf_ev, "%lf \t",data[kk][jj]);   
    }
    fprintf(pf_ev, "\n");
}




fclose(pf_ev);
}
//------------------------------------------------------------------------------------