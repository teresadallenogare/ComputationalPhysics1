#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_fft_complex.h>  /* da fare: /opt/homebrew/Cellar/gsl/2.6/include/ in includePath */
#define N 1024
#define DIM_E 100
#define DIM_t 300

double x[N], dx; 
double x1, x2; // estremi del potenziale
double k[N], dk;
double sigma, x0;
double xi = 1.0;
double complex V[N]; // potential

double L = 500.0; // lenght of the periodic box

double T = 1.0; // evolution of the system for T
double dt; // dt = t/T [s]
// double T = 5.0; //interval of one evolution [s]

double complex rho_n[N]; // density in real space
double complex xi_k[N]; // density in complex space
double complex psi_n[N];
double complex Psi_k[N];

double complex phi[N]; // wavepacket

double data[N][DIM_t];

double E[DIM_E], dE; // energies

double Tr_coeff[DIM_E], Refl_coeff[DIM_E];
double modPhi;

double complex Delta_k(double k){
    if (k<= N/2) // attenta (ci sto dando non il k ma la posizione )
    {
        return I * 2.0 * M_PI * k / L;
    } else {
        return I * 2.0 * M_PI * (k- N) /L;
    }
    
}

// barriera lunga  centrata in L/2 normalizzata con V0 = 1
double potential_barrier(double x){
    if (x <= x2 && x >= x1) // cambiando larghezza potenziale cambia la probabilità
    {  
         return 1.0;
    } else {
        return 0.0;
    }
}

double complex potential(double x){
    double complex pot1, pot2;
    double V0 = 100.0;
    double a = 2.0;
    double x1 = 2.0;
    double x2 = 490.0;
    pot1 = -I * V0 / pow(cosh((x - x1) /a ), 2.0);
    pot2 = - I * V0 / pow(cosh((x - x2) /a), 2.0);
    return potential_barrier(x) + pot1 + pot2;

}

double integral(int index_a, int index_b ){
double f_a, f_b, trapez_sum = 0.0;
f_a = cabs( phi[index_a] * phi[index_a] );
f_b = cabs( phi[index_b] * phi[index_b] );
for (int ii = index_a+1; ii < index_b ; ii++)
{
    trapez_sum += cabs( phi[ii] * phi[ii] );
}
return dx * (f_a/2.0 + f_b/2.0 + trapez_sum);  
}

int maximum_WP1() {
    // finds the position of the maximum 
    int pos = 0; 
    double ris = creal(phi[0]); 
    for (int kk = 1; kk < N/2 - 1; kk++) {
        if (creal(phi[kk]) != 0.00000 && creal(phi[kk]) > ris) {
            ris = creal(phi[kk]); 
            pos = kk; 
        }
    }
    return pos; 
}

int maximum_WP2() {
    // finds the position of the maximum 
    int pos = 0; 
    double ris = creal(phi[0]); 
    for (int kk = N/2 + 1; kk < N ; kk++) {
        if (creal(phi[kk]) != 0.00000 && creal(phi[kk]) > ris) {
            ris = creal(phi[kk]); 
            pos = kk; 
        }
    }
    return pos; 
}

int ctrlDoublePacket(){
    int trovato = 0;
    double eps = 10.0;
    double eps1 = 70.0;
    int pos_wp1_max, pos_wp2_max;
    pos_wp1_max = maximum_WP1();
    pos_wp2_max = maximum_WP2();
        if ( ( x[pos_wp1_max] < L/2.0 - eps1 && x[pos_wp1_max] > eps ) && creal(phi[pos_wp1_max]) > 0.000239  ) // controllo se massimo del I wp è in finestra [eps, L/2 - eps ]
        {
            if ( (x[pos_wp2_max] < L - eps && x[pos_wp2_max] > L/2.0 + eps1) && creal(phi[pos_wp2_max]) > 0.001037  ) // controllo se massimo del II wp è in finestra [L/2 + eps, L -eps ]
            {
                if ( fabs(creal(phi[15]))< 1e-2 && fabs(creal(phi[N/2 - 15])) < 1e-2) // estremi I wp compatibili con 0
                {
                    if (fabs(creal(phi[N/2 + 15])) < 1e-2 && fabs(creal(phi[N-15])) < 1e-2 ) // estremi I wp compatibili con 0
                    {
                        trovato = 1;
                    }
                    else{
                        trovato = 0;
                    }
                    
                }
                
            }
            
        }
   return trovato;     
}

//------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {

FILE *pf, *pf1;
pf = fopen("data_ev.txt", "w");
pf1 = fopen("data_potential.txt", "w");

/* initialization */
x1 = L/2.0 - 1.0/2.0;
x2 = L/2.0 + 1.0/2.0;
dx = L/N;
dk = 2*M_PI/L;
E[0] = 2.0;
dt = 1e-2; // !!!!!!!
x0 = 20.0;
sigma = 10.0;
for (int ii = 0; ii < N; ii++)
{
    x[ii] = dx * ii; // [0, L]
    V[ii] = potential(x[ii]);
}

for (int ii = 0; ii < N; ii++)
{
    rho_n[ii] = cexp (- I * potential(x[ii]) * dt);
    xi_k[ii] = cexp(I * xi * dt * Delta_k(ii) * Delta_k(ii) );
    phi[ii] = exp( -pow( (x[ii] - x0 ), 2.0)/(4.0 * sigma * sigma) ) * cexp( I * sqrt( 4.0* E[0]) * (x[ii] - x0)) ; // gaussian wavepacket (initial condition )
}
   /* modPhi = integral(0, N) ;
    for (int ii = 0; ii < N; ii++)
    {
         phi[ii] = phi[ii] /modPhi;
    }*/

/* evolution of the wp */
double tt;
for (int tau = 0; tau < DIM_t; tau++) 
{   tt = dt;
    while(tt < T ){ // evolution of the system for T
        for (int ii = 0; ii < N; ii++)
        {   
            psi_n[ii] = rho_n[ii] * phi[ii];
        }
        gsl_fft_complex_radix2_forward((double * )psi_n, 1, N);
        for (int ii = 0; ii < N; ii++)
        {
            Psi_k[ii] = psi_n[ii] * xi_k[ii]; // Csi_k [tau + 1]
        }
        gsl_fft_complex_radix2_inverse((double * )Psi_k, 1, N);
        for (int ii = 0; ii < N; ii++)
        {
            psi_n[ii] = Psi_k[ii];
            phi[ii] = psi_n[ii] * rho_n[ii];
        }
        tt += dt;
    }
    for (int ii = 0; ii < N; ii++)
    {
        data[ii][tau] = creal(phi[ii]);
    }

}
for (int ii = 0; ii < N; ii++)
{
    fprintf(pf,"%lf\t ", x[ii]);
    for (int jj = 0; jj < DIM_t; jj++)
    {
       fprintf(pf, "%lf\t", data[ii][jj]);
    }
    fprintf(pf,"\n");
    
}
fprintf(pf, "#x\t\t V(x)\n");
for (int ii = 0; ii < N; ii++)
{
    fprintf(pf1, "%lf \t %lf \n", x[ii], creal(V[ii]));
}

fclose(pf1);
fclose(pf);


 /* transmission and reflection coefficients */
 
FILE *pf_phi, *pf_TR;
pf_phi = fopen("data_phi.txt", "w");
pf_TR = fopen("dataTR.txt", "w");

E[0] = 0.1;
E[DIM_E-1] = 5.0; //1.0 E_max = V0 per stati di tunneling 
dE = (E[DIM_E-1] - E[0])/(DIM_E -1);
for (int ii = 1; ii < DIM_E-1; ii++)
{
    E[ii] = E[ii-1] + dE;
}


double Evolution[2] = {500.0, 150.0}; // massima evoluzione
double Ev_time;
int trovato = 0;
int pos_max1, pos_max2;

for (int kk = 0; kk < DIM_E ; kk++) // cycle on energies
{
    for (int ii = 0; ii < N; ii++)
    {    
        rho_n[ii] = cexp (- I * potential(x[ii]) * dt);
        xi_k[ii] = cexp(I/2.0 * dt * Delta_k(ii) * Delta_k(ii) );
        phi[ii] = exp( -pow( (x[ii] - x0 ), 2.0)/(2.0 * sigma * sigma) ) * cexp( I * sqrt(4.0* E[kk]) * (x[ii] - x0)) ; // ?? OK ma perchè??

    }

    tt = dt;
    trovato = 0;
    if (E[kk] < 2.2)
    {
        Ev_time = Evolution[0];
    }else {
        Ev_time = Evolution[1];
    }
    // begin evoulution over time 
    while(tt < Ev_time && trovato == 0)
    { // evolution of the system for T
        for (int ii = 0; ii < N; ii++)
        {   
           psi_n[ii] = rho_n[ii] * phi[ii];
        }
        gsl_fft_complex_radix2_forward((double * )psi_n, 1, N);
        for (int ii = 0; ii < N; ii++)
        {
           Psi_k[ii] = psi_n[ii] * xi_k[ii]; // Csi_k [tau + 1]
        }
        gsl_fft_complex_radix2_inverse((double * )Psi_k, 1, N);
        for (int ii = 0; ii < N; ii++)
        {
          psi_n[ii] = Psi_k[ii];
          phi[ii] = psi_n[ii] * rho_n[ii];
        }
        trovato = ctrlDoublePacket();
        pos_max1 = maximum_WP1();
        pos_max2 = maximum_WP2();
        tt += dt;
    }
    //end evolution over time 
  printf(" E:%lf\t phi[pos_M1]: %lf \t phi[pos_M2]: %lf \t tt: %lf \t kk: %d \t trv: %d\n",E[kk], creal(phi[pos_max1]), creal( phi[pos_max2]), tt, kk, trovato);
if (trovato ==1)
{
modPhi = integral(0, N) ;
Refl_coeff[kk] = integral(0, N/2-1)/ modPhi;
Tr_coeff[kk] = integral(N/2+1, N)/ modPhi;
} else{
    Refl_coeff[kk] = 0.0;
    Tr_coeff[kk] = 0.0;
}
}
// plot oh phi at Ev_time 
for (int ii = 0; ii < N; ii++)
{
    fprintf(pf_phi, "%lf\t %lf\n", x[ii], creal(phi[ii]) );
}

fprintf(pf_TR, "#E\t\t |R|^2\t\t |T|^2\t\t\t |T|^2 + |R|^2\n");
for (int ii = 0; ii < DIM_E; ii++)
{
    fprintf(pf_TR, "%lf\t%lf\t%lf\t\t %lf \n", E[ii], Refl_coeff[ii], Tr_coeff[ii], Refl_coeff[ii]+Tr_coeff[ii]);
}




fclose(pf_phi);
fclose(pf_TR); 


}
//------------------------------------------------------------------------------------