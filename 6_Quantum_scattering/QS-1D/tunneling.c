#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>

#define DIM_E 3000
#define DIM_xi 3
#define N 50000/* number of Numerov steps is equal to the number of steps */

double xi2[DIM_xi] = {1.0, 3.0, 7.0 }; /* xi2 = 2ma^2V0/hbar^2 */
double x[N]; /* potsitions */
double h = 1e-2; /* interval length */
double x0 = -5.0; /* -inf */
double E[DIM_E]; /* energies (>0) */
double dE = 1e-3; /* energy gap */
double k[DIM_E];
double complex phi[N];

void setToZeroV(int dim, double complex V[dim]);

double rectPotential(double x);
double gaussianPotential(double x);
double asymmetricPotential1(double x);
double asymmetricPotential2(double x);
double F(double (*potential)(double), double x, double xi, double E);
void Numerov(double (*potential)(double), double xi, double E);
double complex transmissionCoefficient(double k);

void fileDataT(double E[DIM_E], double T[DIM_E][DIM_xi], double R[DIM_E][DIM_xi]);

// ----------------------------------------------------------------------------------------
int main(int argc, const char * argv[]){

double xi[DIM_xi];
double T[DIM_E][DIM_xi], R[DIM_E][DIM_xi];
double complex t;
for (int kk = 0; kk < DIM_xi; kk++)
{
    xi[kk] = xi2[kk]* xi2[kk]; /* valori richiesti di xi */
}

x[0] = x0;
for (int kk = 1; kk < N; kk++)
{
    x[kk] = x[kk-1] + h;
}

E[0] = 0.0;
for (int kk = 1; kk < DIM_E; kk++)
{
    E[kk] = E[kk-1] + dE;
}

/* Rectangular barrier */
for (int kk = 0; kk < DIM_xi; kk++)
{
    for (int jj = 0; jj < DIM_E; jj++)
    {
        k[jj] = sqrt( xi[kk] * E[jj]);
        phi[0] = cexp(I * k[jj] * x0);
        phi[1] = cexp(I * k[jj] * (x0+h));
        Numerov(rectPotential, xi[kk], E[jj]);
        t = transmissionCoefficient(k[jj]);
        T[jj][kk] = cabs(t) * cabs(t);
        R[jj][kk] = 1 - T[jj][kk];
    }
    // fileDataT(E, T, R);
}

/* Gaussian barrier */
for (int kk = 0; kk < DIM_xi; kk++)
{
    for (int jj = 0; jj < DIM_E; jj++)
    {
        k[jj] = sqrt( xi[kk] * E[jj]);
        phi[0] = cexp(I * k[jj] * x0);
        phi[1] = cexp(I * k[jj] * (x0+h));
        Numerov(gaussianPotential, xi[kk], E[jj]);
        t = transmissionCoefficient(k[jj]);
        T[jj][kk] = cabs(t) * cabs(t);
        R[jj][kk] = 1 - T[jj][kk];
    }
    //fileDataT(E, T, R);
}

/* Asymmetric barrier 1 */
for (int kk = 0; kk < DIM_xi; kk++)
{
    for (int jj = 0; jj < DIM_E; jj++)
    {
        k[jj] = sqrt( xi[kk] * E[jj]);
        phi[0] = cexp(I * k[jj] * x0);
        phi[1] = cexp(I * k[jj] * (x0+h));
        Numerov(asymmetricPotential1, xi[kk], E[jj]);
        t = transmissionCoefficient(k[jj]);
        T[jj][kk] = cabs(t) * cabs(t);
        R[jj][kk] = 1 - T[jj][kk];
    }
    //fileDataT(E, T, R);
}

/* Asymmetric barrier 2 */
for (int kk = 0; kk < DIM_xi; kk++)
{
    for (int jj = 0; jj < DIM_E; jj++)
    {
        k[jj] = sqrt( xi[kk] * E[jj]);
        phi[0] = cexp(I * k[jj] * x0);
        phi[1] = cexp(I * k[jj] * (x0+h));
        Numerov(asymmetricPotential2, xi[kk], E[jj]);
        t = transmissionCoefficient(k[jj]);
        T[jj][kk] = cabs(t) * cabs(t);
        R[jj][kk] = 1 - T[jj][kk];
    }
     fileDataT(E, T, R);
}


}
// ----------------------------------------------------------------------------------------

void setToZeroV(int dim, double complex V[dim]){
  for (int kk=0; kk < dim; kk++) {
            V[kk]=0;
        }
}

double rectPotential(double x){
    if ( x < 0.5 && x > -0.5)
    {
       return 1.0;
    } else {
       return 0.0;
    }
}
double gaussianPotential(double x){
   return exp(-x*x/2); 
}

double asymmetricPotential1(double x){
  if (x < 0.0 && x > -0.5)
    {
        return 1.0/2.0;
    } 
  if(x < 0.5 && x > 0 ){
        return 1.0;
  }
  else {
        return 0.0;
    }
}

double asymmetricPotential2(double x){
  if (x < 0.0 && x > -0.5)
    {
        return 1.0;
    }
  if(x < 0.5 && x > 0 ){
        return 1.0/2.0;
  }
  else {
        return 0.0;
    }
}

double F(double (*potential)(double), double x, double xi, double E){
    return xi * ( potential(x) - E);
}

void Numerov(double (*potential)(double), double xi, double E){
    double complex num1 = 0, num2 = 0, den = 0;
    
    for (int kk = 2; kk < N; kk++)
    {
        num1 = ( 2.0 + 5.0/6.0 * h*h * F(potential, x[kk-1], xi, E) ) * phi[kk-1];
        num2 = (1.0 - h*h/12.0 * F(potential, x[kk-2], xi, E)) * phi[kk-2];
        den = 1.0 - h*h/12.0 * F(potential, x[kk], xi, E);
        phi[kk] = (num1 - num2)/ den;
    }

}

double complex transmissionCoefficient(double k){
    double complex t;
    double x1, x2; 
    double complex phiD1, phiD2;
    double complex num, den;
    x1 = x[N-11]; 
    x2 = x[N-1]; 
    phiD1 = phi[N-11];
    phiD2 = phi[N-1];
    num = cexp(I*k*(x1-x2)) - cexp(-I*k*(x1-x2)) ;
    den = phiD1* cexp(-I*k*x2) - phiD2* cexp(-I*k*x1);
    t = num/den;
    return t;
}

void fileDataT(double E[DIM_E], double T[DIM_E][DIM_xi], double R[DIM_E][DIM_xi]){

FILE *pf;
pf = fopen("DataT.txt", "w");
fprintf(pf, "# Data transmission coefficient and reflection coefficient depending on E \n");
fprintf(pf, "# E[kk] \t T[csi0] \t R[csi0] \t T[csi1] \t R[csi1] \t T[csi2] \t R[csi2]\n");
for (int jj = 0; jj < DIM_E; jj++)
{   fprintf(pf, "%lf \t", E[jj]);
   for (int kk = 0; kk < DIM_xi; kk++)
   {
       fprintf(pf, "%lf\t %lf\t ", T[jj][kk], R[jj][kk]);
   }
   fprintf(pf, "\n");
}

fclose(pf);

}