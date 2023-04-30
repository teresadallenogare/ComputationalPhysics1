#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#define N 1000 /* number of discretizations of position */
#define DIM_K 100
#define DIM_E 1000
#define DIM_n 3

double K[DIM_K];
double dK;

double x[N];
double h;

double complex phi[N];
double  phi1[N];
double  phi2[N];

double xi = 1.0;
double n[DIM_n];

double E[DIM_E];
double dE;

double complex beta;

 double potential(double x){
    return - exp(-(x - 0.5)* (x - 0.5)/ (2.0 * 0.3 * 0.3));
}

double F(double x, double E){
    return - 1.0/xi * (E - potential(x));
}

void Numerov(double phi_N[N], double E){
    double  num1, num2, den;
    int ii = 2;
    while (ii < N)
    {
     num1 = (2.0 + 5.0/6.0 * h*h * F(x[ii-1], E)) * phi_N[ii-1];
     num2 = (1.0 - h*h/12.0 * F(x[ii-2], E)) * phi_N[ii-2];
     den = 1.0 - h*h/12.0 * F(x[ii], E);
     phi_N[ii] = (num1 - num2)/den;
     ii++;
    }
}

double deltaFunction(double K, double E){ 
return cabs( (cexp(I*K) * phi[1] + phi[N-2] - 2.0 * phi[N-1]) / h );
}


void findMinimum (double delta[DIM_E], double E[DIM_E], double Eminimum[DIM_n], double K){ /* ricerca del minimo di Delta*/
/* function which finds the initial value and which performs bisection method (it returns the energy) */
    double derDelta1, derDelta2;
    double Emin;
    int trovato;
    int pos = 1;
for (int kk = 0; kk < DIM_n; kk++)
{ 
    Emin = 0;
    trovato = 0;
        for (int ii = pos; ii < DIM_E && trovato == 0; ii++)
    {
       derDelta1 = (delta[ii+1] - delta[ii])/dE;
       derDelta2 = (delta[ii+2] - delta[ii+1])/dE;

       if (derDelta1 < 0 && derDelta2 > 0)
       {
           Emin = E[ii+1];
           pos = ii+1;
           trovato = 1;
       }
    } 
 Eminimum[kk] = Emin; 
}   


}

//------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
    FILE *pf, *pf1, *pf2;
    pf = fopen("dataE_K_numerov.txt", "w");
    pf1 = fopen("dataDelta.txt", "w");
    pf2 = fopen("dataPhi.txt", "w");
/* initialization of K */
    K[0] = -M_PI; /* considero K tra 0 e Pi perchè relazione di dispersione è simmetrica */
    K[DIM_K-1] = M_PI;
    dK = (K[DIM_K-1] - K[0])/(DIM_K - 1);
    for (int kk = 1; kk < DIM_K; kk++)
    {
        K[kk] = K[kk-1] + dK;
    }
    
/* initialization of x */
    x[0] = 0.0;
    x[N-1] = 1.0;
    h = (x[N-1] - x[0])/(N-1);
    printf("%lf\n", h);
    for (int kk = 1; kk < N; kk++)
    {
        x[kk] = x[kk-1] + h;
    }
    
/* initialization of E */
    E[0] = - 1.0;
    E[DIM_E-1] = 150.0;
    dE = (E[DIM_E-1] - E[0])/(DIM_E -1);
    for (int kk = 1; kk < DIM_E; kk++)
    {
        E[kk] = E[kk-1] + dE;
    }

    double Delta[DIM_E]; /* Per ogni K fissato provo DIM_E valori di energia e calcolo il delta */
    double E_bm[DIM_n];
    double En[DIM_K][DIM_n];

    for (int kk = 0; kk < DIM_K/* DIM_K */; kk++) /* Fisso un valore di K. Per ogni K devo ciclare su tutte le energie */
    { 
        for (int jj = 0; jj < DIM_E/* DIM_E*/ ; jj++) /* Scorro su tutti i possibili valori di Energia per un dato K fissato */
        { /* C.I: */
          phi1[0] = 0.0;
          phi1[1] = h;
          phi2[0] = 1.0;
          phi2[1] = 1.0;
         Numerov(phi1, E[jj]);
         Numerov(phi2, E[jj]);
        beta = phi1[N-1]/ (cexp(I*K[kk]) - phi2[N-1]);
         // printf("%lf\t %lf \n \n", creal(beta), cimag (beta) );
        for (int ii = 0; ii < N; ii++)
        {
            phi[ii] = phi1[ii] + beta * phi2[ii];
          fprintf(pf2, "%lf\t %lf\n", x[ii], pow(cabs(phi[ii]), 2.0));
        }

        Delta[jj] = deltaFunction(K[kk], E[jj]);
        }
   findMinimum(Delta, E, E_bm, K[kk]);
    En[kk][0] = E_bm[0];    En[kk][1] = E_bm[1];    En[kk][2] = E_bm[2];

    }



    for(int ii = 0; ii < DIM_E; ii++)
    {
        fprintf(pf1, "%lf\t %lf \n", E[ii] , Delta[ii]);
    }

    fprintf(pf, "Energies of the bound states for csi = 1.0 \n"); 
    fprintf(pf, "K\t \t \t En0\t \t \t En1\t \t \t En2\n");
    for (int kk = 0; kk < DIM_K; kk++)
    {
     fprintf(pf, "%lf \t %lf \t %lf\t %lf \n", K[kk], En[kk][0], En[kk][1], En[kk][2]);
    }

fclose(pf);
fclose(pf1);
fclose(pf2);
}
//------------------------------------------------------------------------------------
