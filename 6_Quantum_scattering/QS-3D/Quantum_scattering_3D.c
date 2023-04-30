#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_sf.h>


#define DIM_r 1000
#define DIM_E 2000
#define DIM_l 15

double r[DIM_r];
double h;
double r0 = 1.0;
double r_inf = 6.0;
double u[DIM_r];
double E[DIM_E];
double dE;
double Emax = 3.0;
double k;
double phi1, phi2;
double r1, r2;
double xi = 1.0;

double sigma;
double Sigma[DIM_E];



double potential(double r, int l){
  return xi* l*(l+1)/(r*r);
}

double F(double r, int l, double E){
  return - 1/xi * (E - potential(r,l) );
}

double asySol(double r){
  return exp(-sqrt(2.0 * 4.0/25.0) * pow(r, -5.0) );
}

void Numerov( int l, double E){
    double num1, num2, den;
    int ii = 2;
    while (ii< DIM_r)
    {
        num1 = (2.0 + (5.0/6.0)*h*h*F(r[ii-1],l,E)) * u[ii-1]; 
        num2 = (1 - h*h/12.0*F(r[ii-2],l,E))  * u[ii-2]; 
        den = (1 - (h*h/12.0)*F(r[ii],l,E));  
        u[ii] = (num1 - num2)/den; 
        ii++;
    }
}

double verifyDelta(int l, double k){
  double num, den;
  double temp;
  num = gsl_sf_bessel_jl(l, k) * gsl_sf_bessel_jl(l, k);
  den = gsl_sf_bessel_jl(l, k) * gsl_sf_bessel_jl(l, k) + gsl_sf_bessel_yl(l, k) * gsl_sf_bessel_yl(l, k);
  temp = sqrt(num/den);
  return asin(-temp);
}

// ----------------------------------------------------------------------------------------
int main(int argc, const char * argv[]){

  FILE *pf;

  pf = fopen("dataSigma.txt", "w");

r[0] = r0;
h = (r_inf-r0)/(DIM_r-1);
for (int kk = 1; kk < DIM_r; kk++)
{
  r[kk] = r[kk-1] + h;
}
r1 = r[DIM_r - 10];
r2 = r[DIM_r - 5];

E[0] = 0.0001;
dE = (Emax-E[0])/(DIM_E-1);
for (int kk = 1; kk < DIM_E; kk++)
{
  E[kk] = E[kk-1] + dE;
}

double K;
double jl_r1, jl_r2, nl_r1, nl_r2;
double num, den;

double delta_l;
double delta_lVerify[DIM_E];

for (int jj = 0; jj < DIM_E; jj++)
{
  k = sqrt( E[jj]/ xi);
  sigma = 0;
  for (int ll = 0; ll < DIM_l; ll++)
  {
    u[0] = 0.0;
    u[1] = gsl_sf_bessel_jl(ll, k);
    Numerov(ll, E[jj]);
    phi1 = u[DIM_r - 10]/ r1;
    phi2 = u[DIM_r - 5]/ r2;
    K = phi2/phi1;
    jl_r1 = gsl_sf_bessel_jl(ll, k * r1);
    jl_r2 = gsl_sf_bessel_jl(ll, k * r2);
    nl_r1 = gsl_sf_bessel_yl(ll, k * r1);
    nl_r2 = gsl_sf_bessel_yl(ll, k * r2);
    num = jl_r2 - K * jl_r1;
    den = nl_r2 - K * nl_r1;
    delta_l = atan(num/den);
    // delta_lVerify[jj] = verifyDelta(ll,k);
    sigma += (2*ll + 1) * sin(delta_l) * sin(delta_l);
  }
  Sigma[jj] = 4*M_PI*sigma/(k*k);
  // fprintf(pf, "%lf\t %lf \t %lf\n", E[jj], delta_l[jj], delta_lVerify[jj]); 
fprintf(pf, "%lf\t %lf\n", E[jj], Sigma[jj]); 

// Limit to low energies

}

fclose(pf);
}
// ----------------------------------------------------------------------------------------

