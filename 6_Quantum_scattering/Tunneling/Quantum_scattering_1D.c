#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#define N 50000 /* number of Numerov algorithm evolutions */
#define NUM_EN 500
#define SEL_X1 900

double x0 = -5.0; /* -inf */
double h = 0.01; /* step between two positions */
double X[N]; /* positions in which I evaluate phi */
double x1, x2;

double E[NUM_EN]; /* energy values */

double complex phi[N]; /* wave function values */
double complex phiD1, phiD2; /* wave function values for x1 and x2 */
double V; /* potential */
double V0; 
double d = 1.0;
double k; /* vettore d'onda che dipende dall'energia */
double lambda; /* wavelenght */


double rectangularPotential(double x);
double gaussianPotential(double x);
double asymmetricPotential1(double x);
double asymmetricPotential2(double x);
double F(double x, double En, double (*f) (double x));
void setToZeroV(int dim, double complex V[dim]);
void Numerov(double En, double (*f) (double x));
double complex transCoeff();
void fileDataTrasmission(double E[NUM_EN], double coeffTr1[NUM_EN], double coeffTr2[NUM_EN], double coeffTr3[NUM_EN], double coeffTr4[NUM_EN]);
void fileDataPhi( double X[N], double complex phi[N]);

// ----------------------------------------------------------------------------------------
int main(int argc, const char * argv[]){

double coeffTrR[NUM_EN], coeffTrG[NUM_EN], coeffTrA1[NUM_EN], coeffTrA2[NUM_EN];
double complex TR, TG, TA1, TA2;
/* energy values */
V0 =  d*d/2.0; 
// V0 = 49.2;
E[0] = 0.0;
    for (int kk=0; kk<NUM_EN; kk++) {
        E[kk] = V0*kk*0.01;
    }
/* position values */
X[0] = x0;
for (int kk = 1; kk < N; kk++)
{
    X[kk] =  X[kk-1] + h;
}
/* evolution of phi for each value of energy and calculation of T^2 */
for (int kk = 0; kk < NUM_EN; kk++)
{
  k = sqrt(2.0*E[kk]);
  phi[0] = cexp(I*k*x0);
  phi[1] = cexp(I*k*(x0+h)); 
//phi[1] = I*h*k* cexp(I*k*x0) + cexp(I*k*x0);
  
/* Rectangular barrier */
  Numerov(E[kk], &rectangularPotential);
  TR = transCoeff();
  coeffTrR[kk] = cabs(TR) * cabs(TR);

  /* Gaussian barrier */
  Numerov(E[kk], &gaussianPotential);
  TG = transCoeff();
  coeffTrG[kk] = cabs(TG) * cabs(TG);

  /* Asymmetric barrier 1 */

  Numerov(E[kk], &asymmetricPotential1);
  TA1 = transCoeff();
  coeffTrA1[kk] = cabs(TA1) * cabs(TA1);

    /* Asymmetric barrier 2 */
  Numerov(E[kk], &asymmetricPotential2);
  TA2 = transCoeff();
  coeffTrA2[kk] = cabs(TA2) * cabs(TA2);
}
fileDataPhi(X,phi);
fileDataTrasmission(E, coeffTrR, coeffTrG, coeffTrA1, coeffTrA2);


}
// --------------------------------------------------------------------------------------

double rectangularPotential(double x){
    if (x < 0.5 && x > -0.5)
    {
        return V0;
    }else {
        return 0.0;
    }
}
double gaussianPotential(double x){
   return V0 * exp(-x*x/2); 
}

double  F(double x, double En, double (*f) (double x)){
    return 2.0 * ( f(x) - En);
}

double asymmetricPotential1(double x){
  if (x < 0.0 && x > -0.5)
    {
        return V0/2.0;
    } 
  if(x < 0.5 && x > 0 ){
        return V0;
  }
  else {
        return 0.0;
    }

  
}

double asymmetricPotential2(double x){
  if (x < 0.0 && x > -0.5)
    {
        return V0;
    }
  if(x < 0.5 && x > 0 ){
        return V0/2.0;
  }
  else {
        return 0.0;
    }
}

void Numerov(double En, double (*f) (double x)){
double complex den;
double complex num1, num2;

for (int kk = 2; kk < N; kk++){
    num1 = (2.0 + 5.0/6.0 * h * h * F(X[kk-1], En, f)) * phi[kk-1];
    num2 = (1.0 - h * h /12.0 * F(X[kk-2], En, f)) * phi[kk-2];
    den = 1.0 - h * h /12.0 * F(X[kk], En, f);
    phi[kk] = (num1 - num2) / den;

  }
}

void setToZeroV(int dim, double complex V[dim]){
  for (int kk=0; kk < dim; kk++) {
            V[kk]=0;
        }
}

double complex transCoeff(){
double complex num, den;
double complex T;
  /*lambda = (2.0 * M_PI)/k;
  x1 = X[SEL_X1]; // fissato io (per ora a caso) 
  x2 = x1 - lambda / 2.0;
  //printf("x1: %lf  x2: %lf", x1, x2);
  int sel_x2 = 0;
  int trovato = 0;
  for (int jj = 0; jj < N && trovato == 0 ; jj++)
  {
    if(X[jj] < x2 ){
      sel_x2++;
      
    } else{
      trovato = 1;
    }
  }
  phiD1 = phi[SEL_X1];
  phiD2 = phi[sel_x2]; 
num = cexp(I*k*(x1-x2)) - cexp(-I*k*(x1-x2)); 
den = phiD2 * cexp(-I*k*x2) - phiD1 * cexp(-I*k*x1);
return num/den;*/
    x1 = X[N-11]; 
    x2 = X[N-1];  
    num = cexp(I*k*(x1-x2)) - cexp(-I*k*(x1-x2)); 
    den = cexp(-I*k*x2)*phi[N-11] - cexp(-I*k*x1)*phi[N-1];
    return num/den; 
}

void fileDataTrasmission(double E[NUM_EN], double coeffTr1[NUM_EN], double coeffTr2[NUM_EN], double coeffTr3[NUM_EN], double coeffTr4[NUM_EN]){
  FILE *pf;
pf = fopen("dataTrans.txt", "w");

fprintf(pf, "# E[kk] \t abs(TR)^2[kk]\t abs(TG)^2[kk]\tabs(TA1)^2[kk] \t abs(TA2)^2[kk] \n");
  for (int kk = 0; kk < NUM_EN; kk++) 
  {
    fprintf(pf, " %lf \t %lf \t %lf \t %lf \t %lf \n", E[kk]/V0, coeffTr1[kk], coeffTr2[kk], coeffTr3[kk], coeffTr4[kk])  ; 
  }
fclose(pf);
}

void fileDataPhi( double X[N], double complex phi[N]){
FILE *pf;
pf = fopen("dataPhi.txt", "w");

fprintf(pf, "# X[kk] \t\t phi[kk] \n");
  for (int kk = 0; kk < N; kk++) 
  {
    fprintf(pf, " %lf \t %lf %lf \n", X[kk], creal(phi[kk]), cimag(phi[kk]) )  ; 
  }
fclose(pf);
}
