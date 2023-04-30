#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define R_sun 7.0e8 
#define M_sun 2.0e30 
#define rho0 1.622e5 
#define G 6.6726e-11
#define N 5000
#define DIM_n 1000
#define eps 5e-3

double xi[N], theta[N], eta[N];
double xi_max = 10.0;
double h;
double n[DIM_n];
double n_min = 2.8; 
double n_max = 3.5;
double dn;
int pos_xi0; 

double F_theta(double eta, double xi) {
    return -eta / (xi * xi);
}
double F_eta(double theta, double n, double xi) {
    return pow(theta, n) * xi * xi;
}
void RK4(double n) {
double K1_theta, K2_theta, K3_theta, K4_theta, K1_eta, K2_eta, K3_eta, K4_eta;
for (int jj = 1; jj < N; jj++)
{
    K1_theta = F_theta(eta[jj-1], xi[jj-1]);
    K1_eta = F_eta(theta[jj-1], n, xi[jj-1]);

    K2_theta = F_theta(eta[jj-1] + h/2.0 * K1_theta, xi[jj-1] + h/2.0); 
    K2_eta = F_eta(theta[jj-1] + h/2.0 * K1_eta, n, xi[jj-1] + h/2.0); 

    K3_theta = F_theta(eta[jj-1] + h/2.0 * K2_theta, xi[jj-1] + h/2.0);
    K3_eta = F_eta(theta[jj-1] + h/2.0 * K2_eta, n, xi[jj-1] + h/2.0);

    K4_theta = F_theta(eta[jj-1] + h * K3_theta, xi[jj-1] + h);
    K4_eta = F_eta(theta[jj-1] + h * K3_eta, n, xi[jj-1] + h);

    theta[jj] = theta[jj-1] + h * (1.0/6.0 * K1_theta + 1.0/3.0 * K2_theta + 1.0/3.0 * K3_theta + 1.0/6.0 * K4_theta);
    eta[jj] = eta[jj-1] + h * (1.0/6.0 * K1_eta + 1.0/3.0 * K2_eta + 1.0/3.0 * K3_eta + 1.0/6.0 * K4_eta);

}

}
int find_xi0(double theta[N]){
    double a, b, f_a, f_b, x;
    int jj = 1;
    while (theta[jj-1] * theta[jj] > 0 && jj< N) {
        jj++;
    }
    if (jj == N) {
        return 0; 
    } else {
        a = xi[jj-1];
        b = xi[jj];
        f_a = theta[jj-1];
        f_b = theta[jj];
    }
    x = a - f_a * (b - a)/(f_b - f_a);
    if (x <= (a+h/2.0)) { 
        return jj-1; 
    } else {
        return jj; 
    }
}
double integral(double n){
    double a, f_a, f_b, tr;
    a = xi[0];
    f_a = xi[0] * xi[0];  
    f_b = 0.0; //xi[pos_xi0] * xi[pos_xi0] * pow(theta[pos_xi0],n);
    tr = (f_a + f_b) / 2.0; 
    for (int jj = 1; jj < pos_xi0; jj++) {
        tr +=  xi[jj] * xi[jj] * pow(theta[jj],n);
    }
    return 4 * M_PI * h * tr; 
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

int main(int argc, const char * argv[]){ 

    // initial conditions
    h = (xi_max - xi[0])/(N-1);
    xi[0] = eps; 
    for (int kk = 1; kk < N; kk++){
        xi[kk] = xi[0] + kk * h;
    } 
    theta[0] = 1.0; 
    eta[0] = 0.0; 
    n[0] = n_min;
    dn = (n_max - n_min)/(DIM_n);
    for (int kk = 1; kk < DIM_n; kk++) {
        n[kk] = n_min + kk * dn;
    }
    
    // polytropic index
    double integr[DIM_n], rappSun[DIM_n], c[DIM_n], xi0[DIM_n];
    for (int kk = 0; kk < DIM_n; kk++) {
        RK4(n[kk]);
        pos_xi0 = find_xi0(theta);
        xi0[kk] = xi[pos_xi0]; 
        if (pos_xi0 != 0) {
            integr[kk] = integral(n[kk]);
            rappSun[kk] = M_sun * pow(xi0[kk], 3.0) / (rho0 * pow(R_sun, 3.0));
            c[kk] =  fabs(rappSun[kk] -  integr[kk]); 
             printf("n: %lf \t xi0: %lf \t I: %lf \t Rapp Sun: %lf \t c: %lf \t %d \n", n[kk], xi0[kk], integr[kk], rappSun[kk], c[kk], kk);
        }
    }
    int pos_n = minimum(c); 
    printf("Polytropic index: %lf \n", n[pos_n]); 
    // printf("Xi0: %lf \n", xi0[pos_n]); 

    // plot of theta for some n (integers)
    double theta_int[N][6]; 
    for (int jj=0; jj<6; jj++) {
        RK4(jj); 
        for (int kk=0; kk<N; kk++) {
            theta_int[kk][jj] = theta[kk];
        }   
    }
    FILE *pf;
    pf = fopen("data_Theta.txt", "w");
    for (int kk=0; kk<N; kk++) {
        fprintf(pf, "%lf \t ", xi[kk]); 
        for (int jj=0; jj<6; jj++) {
            fprintf(pf, "%lf \t ", theta_int[kk][jj]);
        }
        fprintf(pf, "\n"); 
    }
    
    // pressure at the center
    double alpha2 = (R_sun*R_sun) / (xi0[pos_n] * xi0[pos_n]); 
    double espo1 = 1.0/(n[pos_n] - 1.0); 
    double K = ( 4*M_PI*alpha2*G ) / ( (n[pos_n]+1) * pow(rho0,espo1) );
    double espo2 = 1.0/n[pos_n] + 1.0; 
    double pressure = K*pow(rho0,espo2); 
    printf("Pressure: %g \n", pressure); 

}