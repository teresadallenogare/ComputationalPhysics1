#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#define DIM_l 10
#define DIM_r 1000
#define DIM_E 2000
#define POS_l 0
#define POS_r0 90 // a occhio
#define DIM_n 20 /* number of possible energy levels */

double r[DIM_r];
double r0;
double rin_F = 0.6;
double rin_B = 10.0;
double h;

double xi = 8.65e-3; 

int l[DIM_l];
double E[DIM_E];
double dE = 5e-4;
double k[DIM_E];

double uF[DIM_r];
double uB[DIM_r];
double uL[DIM_r];
double uR[DIM_r];

double Delta[DIM_E];

double asySol(double r);
void initialConditions(double u[DIM_r]);
double potential(double r, int l);
double derPotential(double r, int l);
double minimumPotential( int l);
double F(double r, int l, double E);
void NumerovF( int l, double E);
void NumerovB( int l, double E);
double deltaFunction( int l, double E);
void bisectionMethod (double delta[DIM_E], int l, double E[DIM_E], double E_bm[DIM_n]);

//------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
double r_min;
FILE *pf, *pf1;

pf = fopen("dataE.txt", "w");
pf1 = fopen("dataX.txt", "w");

r[0] = rin_F;

h = (rin_B - rin_F)/ (DIM_r-1);
for (int kk = 1; kk < DIM_r; kk++)
{
    r[kk] = r[kk-1] + h;
    // printf("%d\t %lf\n", kk, r[kk]);
}
r0 = r[POS_r0];
printf("r0: %lf\n \n", r0);
for (int kk = 0; kk < DIM_l; kk++)
{
    l[kk] = kk;
}
E[0] = minimumPotential(l[POS_l]);
for (int kk = 1; kk < DIM_E; kk++)
{
    E[kk] = E[kk-1] + dE;
}

for (int jj = 0; jj < DIM_E ; jj++)
{
    k[jj] = sqrt( fabs(E[jj])/ xi);
    uF[0] = asySol(rin_F);
    uF[1] = asySol(rin_F+h);
    NumerovF( l[POS_l], E[jj]);
    for (int ii = 0; ii < DIM_r && r[ii]<= r0; ii++)
    {
         uL[ii] = uF[ii];
    }
    uB[DIM_r-1] = exp(- k[jj]* rin_B);
    uB[DIM_r-2] = exp(- k[jj]* (rin_B-h));
    NumerovB(l[POS_l], E[jj]);
    for (int ii = DIM_r-1; ii>= POS_r0 ; ii--)
    {
        uR[ii] = uB[ii]*uF[POS_r0]/ uB[POS_r0];

    }

    for (int ii = 0; ii < DIM_r; ii++)
    {
       fprintf(pf1, "%lf \t %lg \t %lg\n", r[ii], uL[ii], uB[ii]);
    }
    
      
    Delta[jj] = deltaFunction(l[POS_l], E[jj]);
    fprintf(pf, "%lf \t %lf \n", E[jj], Delta[jj] );

}

double En[DIM_n];
bisectionMethod(Delta, l[POS_l], E, En);
printf("Energies of the bound states for l = %d \n", l[POS_l]); 
for (int kk = 0; kk < DIM_n; kk++)
{
    printf("%lf  \n", En[kk]);
}


fclose(pf1);
fclose(pf);
}
//------------------------------------------------------------------------------------

double asySol(double r){
    return exp( -sqrt(4.0/(25.0*xi)) * pow(r, -5.0)  );
}


double potential(double r, int l){
double pot1, pot2;
pot1 = xi * (l*(l+1))/(r*r);
pot2 = 4.0 * (pow(r, -12.0) - pow(r, -6.0));
return pot1+ pot2;
}

double derPotential(double r, int l){
    return -xi*l*(l+1) + 12.0*pow(r, -4.0) - 24.0*pow(r,-10.0);
}

double minimumPotential( int l){
    /* find the minimum value of the potential with respect to l */
double a, b, c, fa, fb, fc;
a = 1.0;    b = 2.0;

while (b - a > 1e-5)
{
    c = (a + b)/2.0;
    fa = derPotential(a, l);
    fb = derPotential(b, l);
    fc = derPotential(c, l);

    if (fa * fc > 0)
    {
        a = c;
    } else{
        b = c;
    }
}
    return potential(c, l);
}

double F(double r, int l, double E){
    return -E/xi + potential(r, l)/xi;
}

void NumerovF( int l, double E){
    double num1, num2, den;
    int ii = 2;
    while (r[ii]<= r0)
    {
        num1 = (2.0 + (5.0/6.0)*h*h*F(r[ii-1],l,E)) * uF[ii-1]; 
        num2 = (1 - h*h/12.0*F(r[ii-2],l,E))  * uF[ii-2]; 
        den = (1 - (h*h/12.0)*F(r[ii],l,E));  
        uF[ii] = (num1 - num2)/den; 
        ii++;
    }
    
}

void NumerovB(int l, double E){
    double num1, num2, den;
    int ii = DIM_r-3;
while (r[ii]>= r0)
    {
        num1 = (2.0 + (5.0/6.0)*h*h*F(r[ii+1], l, E)) * uB[ii+1]; 
        num2 = (1 - h*h/12.0*F(r[ii+2], l, E))  * uB[ii+2]; 
        den = (1 - (h*h/12.0)*F(r[ii], l, E));  
        uB[ii] = (num1 - num2)/den; 
        ii--;
    }
}

double deltaFunction( int l, double E){
    double num1, num2;
     num1 = uL[POS_r0-1]+ uR[POS_r0+1];
     num2 = (2.0 + h*h*F(r[POS_r0], l,  E)) * uL[POS_r0];
      return (num1 - num2)/h;
}


void bisectionMethod (double delta[DIM_E], int l, double E[DIM_E], double E_bm[DIM_n]) {
    /* function which finds the initial value and which performs bisection method (it returns the energy) */
    double Ea, Eb, Ec, delta_Ea, delta_Eb, delta_Ec; 
    int jj; 
    for (int kk = 0; kk < DIM_n && jj<DIM_E; kk++) {  
         jj++;
        while ( delta[jj-1] * delta[jj] > 0.0 && jj < DIM_E ){
            jj++; 
        }
        if (jj == DIM_E || E[jj]>0/*|| fabs(delta[jj-1]- delta[jj]) > eps */) {
            E_bm[kk] = 0.0; /* Energy with bisection method */ 
        } else {
            Ea = E[jj-1]; /* primo estremo */
            Eb = E[jj];  /* secondo estremo */
            while ( fabs(Eb - Ea) > 1e-5 ) {
                Ec = (Ea + Eb)/2.0; 
                delta_Ea = deltaFunction (l, Ea);
                delta_Eb = deltaFunction (l, Eb);
                delta_Ec = deltaFunction (l, Ec);
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

