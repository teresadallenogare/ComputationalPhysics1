#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_fft_complex.h>  /* da fare: /opt/homebrew/Cellar/gsl/2.6/include/ in includePath */

#define N 2048
#define DIM_E 2048
#define POS_x0 1407

double xi = 0.015;
double x[N], dx;
double x0; /* convergence point */

double phiF[N], phiB[N];
double phiL[N], phiR[N];
double Delta[DIM_E];
double E[DIM_E], dE;

double phi0[N];

double V_L(double x){
    if (x < 0.0)
    {
        return - x * x * (x + 1.0);
    } else {
        return 0.0; 
    }
}

double V_R(double x){
    if (x < 0.0)
    {
        return 0.0;
    } else {
        return x * x * (x - 1.0);
    }
    
}

double init_potential(double x){
    /* shape of the potential */
if (x <= 0.0)
    {
        return V_L(x);
    } else {
        return 0.0;
    }
}

double F(double x, double E){
 return - 1.0/xi * (E - init_potential(x) );
}

void NumerovF( double E){
    /* performs the Numerov evolution of phi_F */
    double num1, num2, den;
    int ii = 2;
    while (x[ii]<= x0)
    {
        num1 = (2.0 + (5.0/6.0) * dx * dx * F(x[ii-1], E)) * phiF[ii-1]; 
        num2 = (1.0 - dx * dx /12.0 * F(x[ii-2], E))  * phiF[ii-2]; 
        den = (1.0 - (dx * dx /12.0) * F(x[ii],E));  
        phiF[ii] = (num1 - num2)/den; 
        ii++;
    }
    
}

void NumerovB( double E){
    /* performs the Numerov evolution of phi_B*/
    double num1, num2, den;
    int ii = N-3;
while (x[ii]>= x0)
    {
        num1 = (2.0 + (5.0/6.0)* dx * dx *F(x[ii+1],E)) * phiB[ii+1]; 
        num2 = (1.0 - dx * dx /12.0 * F(x[ii+2],E))  * phiB[ii+2]; 
        den = (1.0 - (dx * dx/12.0) * F(x[ii],E));  
        phiB[ii] = (num1 - num2)/den; 
        ii--;
    }
}

double deltaFunction( double E){
    /* implementation of the Delta function */
    double num1, num2;
    num1 = phiL[POS_x0-1]+ phiR[POS_x0+1];
    num2 = (2.0 + dx * dx * F(x[POS_x0], E)) * phiL[POS_x0];
    return (num1 - num2)/dx;
}

double bisectionMethod (double delta[DIM_E], double E[DIM_E]) {
    /* function which finds the initial value and which performs bisection method (it returns the energy) */
    double Ea, Eb, Ec, delta_Ea, delta_Eb, delta_Ec; 
    int jj = 1;  
        while ( delta[jj-1] * delta[jj] > 0.0 && jj < DIM_E ){
            jj++; 
        }
        if (jj == DIM_E ) {
            return 0.0; /* Energy with bisection method */ 
        } else {
            Ea = E[jj-1]; /* primo estremo */
            Eb = E[jj];  /* secondo estremo */
            while ( fabs(Eb - Ea) > 1e-5 ) {
                Ec = (Ea + Eb)/2.0; 
                delta_Ea = deltaFunction (Ea);
                delta_Eb = deltaFunction (Eb);
                delta_Ec = deltaFunction (Ec);
                if ( delta_Ea * delta_Ec < 0.0 ) {
                    Eb = Ec; 
                } else {
                    Ea = Ec; 
                }
            }
            return Ec;   /* energy with bisection method */
        }
} 

double integral( double f[N], double a, double cutoff){
    /* perfroms integration with trapezoids method */
    double h, f_a, f_b;
    double trapez_sum= 0.0;
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

//------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
double En = 0.0;
FILE  *pf_phi0;

pf_phi0 = fopen("data_phi0.txt", "w");

/* initialization of x */
x[0] = - 2.5; // -infty
x[N-1] = 0.0;
dx = (x[N-1] - x[0])/(N-1.0);
for (int ii = 1; ii < N; ii++)
{
    x[ii] = x[ii-1] + dx;
}
x0 = x[POS_x0];
printf("x0: %lf\n", x0);

/* initialization of E */
E[0] = - 4.0/27.0; // minimum value of the potential 
E[DIM_E - 1] = 0.0;
dE = (E[DIM_E -1] - E[0])/ (DIM_E -1.0);
for (int ii = 1; ii < DIM_E; ii++)
{
     E[ii] = E[ii-1] + dE;
}

/* Evolution of phiF and phiB */
for (int ii = 0 ; ii < DIM_E; ii++)
{ /* initial conditions of phiF*/
   phiF[0] = 1e-7;
   phiF[1] = 1e-4;
   NumerovF(E[ii]); 
   /* definition of phi_L */
   for (int jj = 0; jj < N && x[jj] <= x0; jj++)
    {
        phiL[jj] = phiF[jj];
    }
    /* initial conditions of phi_B */
   phiB[N -1] = 0.0;
   phiB[N-2] = 1e-7;
   NumerovB(E[ii]);
   /* definition of phi_R*/
   for (int jj = N-1; jj>= POS_x0 ; jj--)
    {
        phiR[jj] = phiB[jj] * phiF[POS_x0]/ phiB[POS_x0];
    }
    Delta[ii] = deltaFunction(E[ii]);
}

En = bisectionMethod(Delta, E); /* eigenvalue of the fundamental state */ 
printf("Eigenvalue En: %lf  \n", En);

/* build eigenfunction phi0 */
   phiF[0] = 0.0;
   phiF[1] = 1e-10;
   NumerovF(En);
   for (int jj = 0; jj < N && x[jj] <= x0; jj++)
   { 
        phiL[jj] = phiF[jj];
   }
   phiB[N -1] = 0.0000;
   phiB[N-2] = 1e-10;
   NumerovB(En);
   for (int jj = N-1; jj>= POS_x0 ; jj--)
    {
        phiR[jj] = phiB[jj] * phiF[POS_x0]/ phiB[POS_x0];
    }

for (int ii = 0; ii < N; ii++)
{
    if (x[ii]<= x0)
    {
        phi0[ii] = phiL[ii];
    } else {
        phi0[ii] = phiR[ii];
    }
}
/* normalization of phi0 */
double modPhi0, mm[N];
for (int ii = 0; ii < N; ii++)
{
    mm[ii] = fabs(phi0[ii]) * fabs(phi0[ii]); 
}
modPhi0 = integral(mm, x[0], 0.0);
printf("modPhi: %lf\n", modPhi0);

fprintf(pf_phi0, "#r\t\t phi0\n");
for (int ii = 0; ii < N; ii++)
{
    phi0[ii] = phi0[ii] * 1.0/ sqrt(modPhi0);
fprintf(pf_phi0, "%lf\t%g\n", x[ii], phi0[ii]);
}

fclose(pf_phi0);
}
