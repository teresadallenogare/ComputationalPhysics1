#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define N 5000 
#define R0 2 // adimensional parameter

double t[N], s[N], i[N], r[N]; // functions and time
double t_max = 30.0; // chosen from the plot
double h; // time step 

double F_s(double i, double s) {
    // function for RK4 for s(t)
    return -R0 * i * s;
}
double F_i(double i, double s) {
    // function for RK4 for i(t)
    return R0 * i * s - i;
}
double F_r(double i) {
    // function for RK4 for r(t)
    return i;
}
void RK4() {
    // performs Runge-Kutta method (4th order) of s(t), i(t) and r()
    double K1_s, K2_s, K3_s, K4_s, K1_i, K2_i, K3_i, K4_i, K1_r, K2_r, K3_r, K4_r;
    for (int jj = 1; jj < N; jj++) {
        
        K1_s = F_s(i[jj-1], s[jj-1]);
        K1_i = F_i(i[jj-1], s[jj-1]);
        K1_r = F_r(i[jj-1]);

        K2_s = F_s(i[jj-1] + (h/2.0) * K1_i, s[jj-1] + (h/2.0) * K1_s);
        K2_i = F_i(i[jj-1] + (h/2.0) * K1_i, s[jj-1] + (h/2.0) * K1_s);
        K2_r = F_r(i[jj-1] + (h/2.0) * K1_i);

        K3_s = F_s(i[jj-1] + (h/2.0) * K2_i, s[jj-1] + (h/2.0) * K2_s);
        K3_i = F_i(i[jj-1] + (h/2.0) * K2_i, s[jj-1] + (h/2.0) * K2_s);
        K3_r = F_r(i[jj-1] + (h/2.0) * K2_i);

        K4_s = F_s(i[jj-1] + h * K3_i, s[jj-1] + h * K3_s);
        K4_i = F_i(i[jj-1] + h * K3_i, s[jj-1] + h * K3_s);
        K4_r = F_r(i[jj-1] + h * K3_i);

        s[jj] = s[jj-1] + (h/6.0) * (K1_s + 2.0 * K2_s + 2.0 * K3_s + K4_s);
        i[jj] = i[jj-1] + (h/6.0) * (K1_i + 2.0 * K2_i + 2.0 * K3_i + K4_i);
        r[jj] = r[jj-1] + (h/6.0) * (K1_r + 2.0 * K2_r + 2.0 * K3_r + K4_r);
    }
}

int main(int argc, const char * argv[]){ 

    // initial conditions
    t[0] = 0.0; 
    h = (t_max - t[0])/(N-1);
    for (int kk = 1; kk < N; kk++){
        t[kk] =  kk * h;
    } 
    s[0] = 1.0 - 1e-5; 
    i[0] = 1e-5; 
    r[0] = 0.0; 

    // evolution 
    RK4();
    
    FILE *pf;
    pf = fopen("data.txt", "w"); //plot 
    for (int kk=0; kk<N; kk++) {
        fprintf(pf, "%lf \t %lf \t %lf \t %lf \n ", t[kk], s[kk], i[kk], r[kk]); 
    }
}