#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define R_sun 7.0e8 // [m]
#define M_sun 2.0e30 // [kg]
#define rho0 1.622e5 // [kg m^-3]
#define G 6.67428e-11
#define DIM_xi 500
#define DIM_n 1000

double xi[DIM_xi];
double xi_inf = 10.0;
double h;

double xi_0;
int pos_xi_0;

double theta[DIM_xi];
double eta[DIM_xi];

double n[DIM_n];
double n_max = 3.7;
double dn;

double c[DIM_n];

double F1(double eta, double xi){
    return -eta/(xi * xi);
}

double F2(double theta, double n, double xi){
    return pow(theta, n) * xi*xi;
}

void RK4(double n ){
double K1_t, K2_t, K3_t, K4_t;
double K1_e, K2_e, K3_e, K4_e;
for (int jj = 1; jj < DIM_xi; jj++)
{
    K1_t = F1(eta[jj-1], xi[jj-1]);
    K1_e = F2(theta[jj-1], n, xi[jj-1]);

    K2_t = F1(eta[jj-1] + h/2.0 * K1_t, xi[jj-1] + h/2.0); 
    K2_e = F2(theta[jj-1] + h/2.0 * K1_e, n, xi[jj-1] + h/2.0); 

    K3_t = F1(eta[jj-1] + h/2.0 * K2_t, xi[jj-1] + h/2.0);
    K3_e = F2(theta[jj-1] + h/2.0 * K2_e, n, xi[jj-1] + h/2.0);

    K4_t = F1(eta[jj-1] + h * K3_t, xi[jj-1] + h);
    K4_e = F2(theta[jj-1] + h * K3_e, n, xi[jj-1] + h);
    theta[jj] = theta[jj-1] + h * (1.0/6.0 * K1_t + 1.0/3.0 * K2_t + 1.0/3.0 * K3_t + 1.0/6.0 * K4_t);
    eta[jj] = eta[jj-1] + h * (1.0/6.0 * K1_e + 1.0/3.0 * K2_e + 1.0/3.0 * K3_e + 1.0/6.0 * K4_e);

}

}

double find_xi_0(){
    double a, b, f_a, f_b;
    int jj = 1;
while (theta[jj-1] * theta[jj] > 0 && jj< DIM_xi)
{
    jj++;
}
if (jj == DIM_xi)
{
    return 0.0; /* valore di xi_0 ad un dato n */
} else {
a = xi[jj-1];
b = xi[jj];
f_a = theta[jj-1];
f_b = theta[jj];

 return a - f_a * (b - a)/(f_b - f_a);
}
}

int position_xi_0(){ /* circa  non è proprio giusta */
    int jj = 1;
    while (theta[jj-1] * theta[jj] > 0 && jj< DIM_xi)
    {
        jj++;
    } if(jj == DIM_xi){
        return 0;
    } else
    {
        return jj-1;
    }
 
}

double integral(double xi_0, double n){
double a, b, f_a, f_b;
double  xi_jj;
double trapez_sum = 0;
int dim_xi;

a = xi[0];
b = xi_0;
f_a = 4 * M_PI* xi[0] * xi[0] * pow(theta[0], n); 
f_b = 0.0; /* theta(xi_0) = 0 */
dim_xi = position_xi_0(); /* definisce il limite di integrazione per theta */
//printf("dimXi: % d\n", dim_xi);
h = (xi_inf - xi[0])/ (DIM_xi -1);
//printf("h:%lf\n", h);
for (int jj = 1; jj < dim_xi-1; jj++)
{
    trapez_sum += 4 * M_PI * xi[jj] * xi[jj] * pow(theta[jj], n);
}
return  h * ((f_a + f_b)/2.0 + trapez_sum);
}

int minimum(double c[DIM_n]) {
    int pos = 0; 
    double ris = c[0]; 
    for (int kk = 1; kk < DIM_n; kk++) {
        if (c[kk] != 0.00000 && c[kk] < ris) {
            ris = c[kk]; 
            pos = kk; 
        }
    }
    return pos; 
}

// ----------------------------------------------------------------------------------
int main(int argc, const char * argv[]){ 

    FILE *pf;

    pf = fopen("dataTh_Eta.txt", "w");

xi[0] = 1e-6;
h = (xi_inf - xi[0])/ (DIM_xi -1);
for (int kk = 1; kk < DIM_xi; kk++)
{
    xi[kk] = xi[kk-1] + h;
}

n[0] = 2.8;
dn = (n_max - n[0])/(DIM_n);
for (int kk = 1; kk < DIM_n; kk++)
{
    n[kk] = n[kk-1] + dn;
    //printf("%lf\n", n[kk]);
}

theta[0] = 1.0;
eta[0] = 0.0;
double eps = 5e-2;
double rappSun;
double n_Sun;
double integr;

double pressure = 0.0;
double K;
double num,den;
for (int kk = 0; kk < DIM_n; kk++) /* evoluzione con RK per tutti i valori di n */
{
    RK4(n[kk]);
    xi_0 = find_xi_0();
    //pos_xi_0 = position_xi_0();
    if (xi_0 != 0.0)
    {
    integr = integral(xi_0, n[kk]);
    rappSun = M_sun * pow(xi_0, 3.0)/(rho0 * pow(R_sun, 3.0));
    c[kk] = fabs(rappSun-integr);
    printf("n: %lf\t xi0: %lf\t I: %lf\t Rapp Sun: %lf \t c: %lf\n", n[kk], xi_0, integr, rappSun, c[kk]);

        if (integr < rappSun + eps && integr > rappSun - eps)
        {
        n_Sun = n[kk];
        num = 4*M_PI*R_sun*R_sun* G;
        den =  (n_Sun+1) * pow(rho0, 1.0/(n_Sun-1.0)) * xi_0*xi_0;
        K = num/den;
        pressure = K * pow(rho0, (n_Sun+1.0)/n_Sun) * pow(theta[0],n_Sun+1.0 ); 
        printf("press: %lf\n", pressure);
        }
    }
}
    int pos_n = minimum(c); 
    printf("Polytropic index: %lf \n", n[pos_n]); 

for (int kk = 0; kk < DIM_xi; kk++)
{
      fprintf(pf, "%lf\t %lf \t %lf\n", xi[kk], theta[kk], eta[kk] );
}





fclose(pf);

}
// ----------------------------------------------------------------------------------
