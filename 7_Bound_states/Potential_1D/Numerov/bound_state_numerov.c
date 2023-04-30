#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#define DIM_x 1002
#define DIM_E 2000
#define DIM_xi 3
#define POS_xi 0
#define POS_x0 530
#define DIM_n 15 /* maximum number of possible zeroes */

#define N 1000


double x[DIM_x];
double x0; /* convergence point */
double xin_F = -5.0;/* -inf */
double xin_B = 5.0;
double h; /* interval length*/
double E[DIM_E];
double dE = 0.0005;
double xi[DIM_xi]= {0.05, 0.01, 0.005}; /* xi = hbar^2 / ma^2|vo|*/
double k[DIM_E]; 
double phiF[DIM_x];
double phiB[DIM_x];
double phiL[DIM_x];
double phiR[DIM_x];
double Delta[DIM_E];
double eps[3]= {2.0, 0.12, 0.008};

double R = 5.0; 


double potential(double x);
double F(double x, double xi, double E);

void NumerovF(double xi, double E);
void NumerovB(double xi, double E);
double deltaFunction(double E);
void bisectionMethod(double delta[DIM_E], double E[DIM_E], double xi, double E_bm[DIM_n]);


//------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {


/* NUMEROV */
double num1, num2;
int n[DIM_n];
double En[DIM_n];

FILE  *pf_E;
pf_E = fopen("dataE.txt", "w");

x[0] = xin_F;
h = (xin_B - xin_F)/ (DIM_x-1) ; 
for (int kk = 1; kk < DIM_x; kk++)
{
    x[kk] = x[kk-1] + h;
}
x0 = x[POS_x0];

E[0] = -1.0;
fprintf(pf_E, "#E[kk] \t Delta[kk] \n");
for (int kk = 1; kk < DIM_E; kk++)
{
    E[kk] = E[kk-1] + dE;
}


for (int jj = 0; jj < DIM_E; jj++)
{
    k[jj] = sqrt(2.0 * fabs(E[jj])/ xi[POS_xi]);
    phiF[0] = exp(k[jj] * xin_F);
    phiF[1] = exp(k[jj] * (xin_F + h));
    NumerovF(xi[POS_xi], E[jj]);
    for (int ii = 0; ii < DIM_x && x[ii]<= x0; ii++)
    {
         phiL[ii] = phiF[ii];
    }


    phiB[DIM_x-1] = exp(- k[jj]* xin_B);
    phiB[DIM_x-2] = exp(-k[jj]* (xin_B-h));
    NumerovB(xi[POS_xi], E[jj]);
    for (int ii = DIM_x-1; ii>= POS_x0 ; ii--)
    {
        phiR[ii] = phiB[ii]*phiF[POS_x0]/ phiB[POS_x0];
    }
      
     Delta[jj] = deltaFunction(E[jj]);
     fprintf(pf_E, "%lf\t %lf \n", E[jj], Delta[jj]);
}

fclose(pf_E);



}

//------------------------------------------------------------------------------------

double potential(double x){
    return - 1/ pow(cosh(x), 4.0);
}

double F(double x, double xi, double E){
    return - 2/xi * ( -potential(x) + E );
}

void NumerovF(double xi, double E){
    double num1, num2, den;
    int ii = 2;
    while (x[ii]<= x0)
    {
        num1 = (2.0 + (5.0/6.0)*h*h*F(x[ii-1],xi,E)) * phiF[ii-1]; 
        num2 = (1 - h*h/12.0*F(x[ii-2],xi,E))  * phiF[ii-2]; 
        den = (1 - (h*h/12.0)*F(x[ii],xi,E));  
        phiF[ii] = (num1 - num2)/den; 
        ii++;
    }
    
}

void NumerovB(double xi, double E){
    double num1, num2, den;
    int ii = DIM_x-3;
while (x[ii]>= x0)
    {
        num1 = (2.0 + (5.0/6.0)*h*h*F(x[ii+1],xi,E)) * phiB[ii+1]; 
        num2 = (1 - h*h/12.0*F(x[ii+2],xi,E))  * phiB[ii+2]; 
        den = (1 - (h*h/12.0)*F(x[ii],xi,E));  
        phiB[ii] = (num1 - num2)/den; 
        ii--;
    }
}

double deltaFunction( double E){
    double num1, num2;
     num1 = phiL[POS_x0-1]+ phiR[POS_x0+1];
     num2 = (2.0 + h*h*F(x[POS_x0], xi[POS_xi], E)) * phiL[POS_x0];
      return (num1 - num2)/h;
}

void bisectionMethod (double delta[DIM_E], double E[DIM_E], double xi, double E_bm[DIM_n]) {
    /* function which finds the initial value and which performs bisection method (it returns the energy) */
    double Ea, Eb, Ec, delta_Ea, delta_Eb, delta_Ec; 
    int jj; 
    for (int kk = 0; kk < DIM_n && jj<DIM_E; kk++) {  
         jj++;
        while ( delta[jj-1] * delta[jj] > 0.0 && jj < DIM_E ){
            jj++; 
        }
        if (jj == DIM_E || fabs(delta[jj-1]- delta[jj]) > eps[POS_xi]) {
            E_bm[kk] = 0.0; /* Energy with bisection method */ 
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
            E_bm[kk] = Ec;   /* energy with bisection method */
        }
            
    }
} 

