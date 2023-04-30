#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>


#define DIM_t 1000000
// #define N 200 /* number of RK4 evolutions */
#define H0 7.2082e-11 //[ years-1]
#define Omega_0 0.27
#define Omega_lamb 0.73
 
double t[DIM_t];
double t_inf = -13.8e9; // [years]
double a[DIM_t]; /* universal scale factor */

double h;

double F(double a, double t){
    return H0 * sqrt(Omega_0/a + Omega_lamb * a*a );
}
void RK4( ){
double K1, K2, K3, K4;

for (int jj = 1; jj < DIM_t; jj++)
{
    K1 = F(a[jj-1], t[jj-1]);
    K2 = F(a[jj-1] + h/2.0 * K1, t[jj-1] + h/2.0); 
    K3 = F(a[jj-1] + h/2.0 * K2, t[jj-1] + h/2.0);
    K4 = F(a[jj-1] + h * K3, t[jj-1] + h);
    a[jj] = a[jj-1] + h * (1.0/6.0 * K1 + 1.0/3.0 * K2 + 1.0/3.0 * K3 + 1.0/6.0 * K4);
}
}


// ----------------------------------------------------------------------------------
int main(int argc, const char * argv[]){ 
    FILE *pf, *pf1;
    pf = fopen("data_a.txt", "w");
    pf1 = fopen("data_aFuture.txt", "w");
    t[0] = 0.0;
    h = (t_inf - t[0] )/DIM_t;
    for (int kk = 1; kk < DIM_t; kk++)
    {
        t[kk] = t[kk-1] + h; /* vado indietro nel tempo: C.I. presente */
    }
    a[0] = 1.0;
    RK4();

    int trovato = 0;
    fprintf(pf, "#t[years]\t a(t)\n");
    
    for (int kk = 0; kk < DIM_t; kk++)
    {
            fprintf(pf, "%lf\t %lf\n", t[kk], a[kk]);
        
    }

    for (int kk = 0; kk < DIM_t && trovato == 0; kk++)
    {
        while (a[kk]>=0)
        {
            kk++;
        }
        trovato = 1;
        printf("t_min: %lg years\n", t[kk-1]);
        
    }
    trovato = 0;
    for (int kk = 0; kk < DIM_t && trovato == 0; kk++)
    {
        while(a[kk] > 0.5){
            kk++;
        }
        trovato = 1;
        printf("t_1/2: %lg years\n", t[kk-1]);
        
    }

    double t_future[DIM_t];
    t_future[0] = 0.0;
    h = (t_future[0] - t_inf )/DIM_t;
    for (int kk = 1; kk < DIM_t; kk++)
    {
        t_future[kk] = t_future[kk-1] + h; /* vado indietro nel tempo: C.I. presente */
    }
    a[0] = 1.0;
    RK4();

    for (int kk = 0; kk < DIM_t; kk++)
    {
        fprintf(pf1, "%lf\t %lf\n", t_future[kk], a[kk]);
    }

     trovato = 0;
    for (int kk = 0; kk < DIM_t && trovato == 0; kk++)
    {
        while(a[kk] <2.0){
            kk++;
        }
        trovato = 1;
        printf("t_2: %lg years\n", t_future[kk]);
        
    }
    
    


fclose(pf);
fclose(pf1);
}
// ----------------------------------------------------------------------------------
