#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>

#define dim_n 15 
#define dim_E 200 
#define N 50000 /* number of Numerov integrations */
#define eps 1e-8 
#define dim_alpha 5000 

double csi[3] = {0.05, 0.01, 0.005};
int n[dim_n]; 
double E[dim_E]; /* energies for bound states */
double E_s[dim_E]; /* enerigies for scattering */
double dE = 0.005; 
double alpha[dim_alpha]; /* variational parameter */
double dalpha = 0.0001; 
double f_alpha[dim_alpha][3]; 
double csi_s[3] = {10, 40, 200}; 
double X[N];  /* trajectory */
double x0 = -2.0; /* initial position (-inf) */
double k[dim_E]; 
double f_E0[dim_E][dim_n]; 
double f_E1[dim_E][dim_n];
double f_E2[dim_E][dim_n];
double h = 0.001; /* Numerov spatial step */
double complex Phi[N]; 

double potential (double x) {
    return -pow(cosh(x),-4.0) ;
}
double integrand_BS (double x, double E) {
    return sqrt( (E - potential(x)) );
}
double integrand_VP (double x, double alpha, double csi ) {
    double p1 = exp(-2*(x*x/(alpha*alpha))); 
    double p2 = csi * ( -2/(alpha*alpha) + 4*x*x/pow(alpha, 4) ); 
    return p1* (p2 - potential(x)); 
}
double integrate_BS (double (*f)(const double, const double), const double a, const double b, const int n, const double E){
    /* function which integrates f from a to b with n passes for a given value of energy E */
    double h = (b-a)/(n-1);
    double x = f(a,E) + f(b,E);
    for (int ii=1; ii<n-1; ii++) {
        x = x + 2*f(a+ii*h,E);
    }
    return x * h/2;
}
double integrate_VP (double (*f)(const double, const double, const double), const double a, const double b, const int n, const double alpha, const double csi){
    /* function which integrates f from a to b with n passes for given values of alpha and csi */
    double h = (b-a)/(n-1);
    double x = f(a,alpha,csi) + f(b,alpha,csi);
    for (int ii=1; ii<n-1; ii++) {
        x = x + 2*f(a+ii*h,alpha,csi);
    }
    return x * h/2;
}
void fileData_BS (double f_E[dim_E][dim_n]) {
    /* creates a file of data */
    FILE *pf;
    pf = fopen("f_E.txt", "w");
    fprintf(pf, "# sulle righe ci sono gli n, sulle colonne i valori di E \n # \t E[jj] \t n=0 \t \t n=1 \t \t n=2 \t \t  ... \n"); 
    
    for (int jj=0; jj<dim_E; jj++) {
        fprintf(pf, "%lf \t", E[jj]); 
        for (int kk=0; kk<dim_n; kk++) {
            fprintf(pf, "%lf \t", f_E[jj][kk]); 
        }
        fprintf(pf, "\n"); 
    }     
    fclose(pf); 
} 
double f_BS (double E, int n, double csi) {
    /* function which represents BS quantization rule */
    double a, b, integral; 
    a = - acosh ( pow( 1.0 / fabs(E), 1.0/4.0) ) + eps; 
    b = acosh ( pow( 1.0 / fabs(E), 1.0/4.0) ) - eps; 
    integral = integrate_BS (&integrand_BS, a, b, 10000, E); 
    return integral - M_PI * sqrt(csi/2.0) * (n + 0.5); 
}
void bisectionMethod (double f_E[dim_E][dim_n], double E[dim_E], int n[dim_n], double csi, double E_bm[dim_n]) {
    /* function which finds the initial value and which performs bisection method (it returns the energy) */
    double Ea, Eb, Ec, f_Ea, f_Eb,f_Ec; 
    int jj; 
    for (int kk=0; kk<dim_n; kk++) {
        jj=1;
        while ( f_E[jj-1][kk]*f_E[jj][kk] > 0.0 && jj < dim_E ){
            jj++; 
        }
        if (jj == dim_E) {
            E_bm[kk] = 0.0; 
        } else {
            Ea = E[jj-1];
            Eb = E[jj];  
            while ( fabs(Eb-Ea) > 1e-5 ) {
                Ec = (Ea+Eb)/2.0; 
                f_Ea = f_BS (Ea, n[kk], csi);
                f_Eb = f_BS (Eb, n[kk], csi);
                f_Ec = f_BS (Ec, n[kk], csi);
                if ( f_Ea * f_Ec < 0.0 ) {
                    Eb = Ec; 
                } else {
                    Ea = Ec; 
                }
            }
            E_bm[kk] = Ec;   
        }
            
    }
} 
void fileData_VP (double alpha[dim_alpha], double f_alpha[dim_alpha][3], int dim) {
    /* creates a file of data */
    FILE *pf_a;
    pf_a = fopen("f_alpha.txt", "w");
    for (int jj=0; jj<dim_alpha; jj++) {
            fprintf(pf_a, "%lf \t %lf \n", alpha[jj], f_alpha[jj][dim]); 
    } 
    fclose(pf_a);   
}
void minimum (int pos[3]) {
    /* function which finds the minimun of energy and it returns the corresponding position */
    double min; 
    pos[0] = pos[1] = pos[2] = 0.0; 
    for (int kk=0; kk<3; kk++) {
        min = f_alpha[0][kk];
        for (int jj=1; jj<dim_alpha; jj++) {
            if (f_alpha[jj][kk] < min) {
                min = f_alpha[jj][kk]; 
                pos[kk] = jj; 
            }
        }
    } 
}
double F (const double x, const double E, double csi) {
    /* function for Numerov method */
    return csi * (-potential(x) - E); 
}
void Numerov (double X[N], double complex Phi[N], double E, double csi) { 
    /* performs Numerov method to calculated phi at different points */
    int kk=2; 
    while (kk<N) {
        double complex num1, num2, den; 
        num1 = (2.0 + (5.0/6.0)*h*h*F(X[kk-1],E,csi)) * Phi[kk-1]; 
        num2 = (1 - h*h/12.0*F(X[kk-2],E,csi))  * Phi[kk-2]; 
        den = (1 - (h*h/12.0)*F(X[kk],E,csi));  
        Phi[kk] = (num1 - num2)/den; 
        kk++; 
    }
} 
double complex transmissionCoefficient (double const X[N], double complex const Phi[N], double k) {
    /* calculates the transmission coefficient for a given value of energy (k)*/
    double complex num, den;
    double x1, x2;
    x1 = X[N-15]; 
    x2 = X[N-8];  
    num = cexp(I*k*(x1-x2)) - cexp(-I*k*(x1-x2)); 
    den = cexp(-I*k*x2)*Phi[N-15] - cexp(-I*k*x1)*Phi[N-8];
    return num/den; 
}
void fileDataT (double E[dim_E], double T[dim_E][3]) {
    /* creates a file of data */
    FILE *fd;     
    fd = fopen("dataT.txt","w"); 
    for (int kk=0; kk<dim_E; kk++) {
      fprintf (fd, "%lf \t ", E[kk]); 
        for (int jj=0; jj<3; jj++) {
            fprintf(fd, "%lf \t ", T[kk][jj]); 
        }    
        fprintf(fd, "\n"); 
    }
  fclose(fd);
}

int main(int argc, const char * argv[]) {

    // FIRST PART: bound states - semiclassical approach
    
    /* initialization */
    for (int ii=0; ii<dim_n; ii++) {
        n[ii] = ii; 
    }
    for (int ii=0; ii<dim_E; ii++) {
        E[ii] =  -1.0 + ii*dE; 
    }
 
    /* evaluation of the integrals */
    double a, b, integral; 
    for (int kk=0; kk<dim_n; kk++) { 
        for (int jj=0; jj<dim_E; jj++) {
            a = - acosh ( pow( 1.0 / fabs(E[jj]), 1.0/4.0) ) + eps; 
            b = acosh ( pow( 1.0 / fabs(E[jj]), 1.0/4.0) ) - eps; 
            integral = integrate_BS (&integrand_BS, a, b, 10000, E[jj]); 
            f_E0[jj][kk] = integral - M_PI * sqrt(csi[0]/2.0) * (n[kk] + 0.5); 
            f_E1[jj][kk] = integral - M_PI * sqrt(csi[1]/2.0) * (n[kk] + 0.5); 
            f_E2[jj][kk] = integral - M_PI * sqrt(csi[2]/2.0) * (n[kk] + 0.5); 
        }
    }
    // fileData_BS (f_E0); 

    /* bisection method */
    double E_bm0[dim_n], E_bm1[dim_n], E_bm2[dim_n]; 
    bisectionMethod (f_E0, E, n, csi[0]/2.0, E_bm0); 
    bisectionMethod (f_E1, E, n, csi[1]/2.0, E_bm1); 
    bisectionMethod (f_E2, E, n, csi[2]/2.0, E_bm2); 
    
    printf("SEMICLASSICAL APPROACH \n"); 
    printf("Energies of the bound states for csi = %lf \n", csi[0]/2.0); 
    for (int ii=0; ii<4; ii++) {
        printf ("%lf \t n = %d \n", E_bm0[ii], ii); 
    } 
    printf("Energies of the bound states for csi = %lf \n", csi[1]/2.0); 
    for (int ii=0; ii<8; ii++) {
        printf ("%lf \t n = %d \n", E_bm1[ii], ii); 
    } 
    printf("Energies of the bound states for csi = %lf \n", csi[2]/2.0); 
    for (int ii=0; ii<12; ii++) {
        printf ("%lf \t n = %d \n", E_bm2[ii], ii); 
    } 

    // SECOND PART: ground state - variational principle

    double ex1 = -5.0; /* first extreme */
    double ex2 = 5.0; /* second extreme */
    double A[dim_alpha]; /* constant factor */
    for (int ii=0; ii<dim_alpha; ii++) {
        alpha[ii] = 0.3 + ii*dalpha;  
        A[ii] = sqrt( 2/(M_PI* alpha[ii]*alpha[ii]) ); 
    }
    for (int ii=0; ii<dim_alpha; ii++) {
        f_alpha[ii][0] = -A[ii]*integrate_VP (&integrand_VP, ex1, ex2, 50000, alpha[ii], csi[0]/2.0); 
        f_alpha[ii][1] = -A[ii]*integrate_VP (&integrand_VP, ex1, ex2, 50000, alpha[ii], csi[1]/2.0); 
        f_alpha[ii][2] = -A[ii]*integrate_VP (&integrand_VP, ex1, ex2, 50000, alpha[ii], csi[2]/2.0); 
    }

    int pos[3]; 
    minimum(pos); 
    printf("VARIATIONAL PRINCIPLE \n"); 
    printf("Ground state energy for csi = %lf : \n %lf \t (for alpha = %lf) \n", csi[0]/2.0, f_alpha[pos[0]+1][0], alpha[pos[0]+1]); 
    printf("Ground state energy for csi = %lf : \n %lf \t (for alpha = %lf) \n", csi[1]/2.0, f_alpha[pos[1]+1][1], alpha[pos[1]+1]); 
    printf("Ground state energy for csi = %lf : \n %lf \t (for alpha = %lf) \n", csi[2]/2.0, f_alpha[pos[1]+1][2], alpha[pos[2]+1]); 

    // fileData_VP (alpha, f_alpha, 0); 

    // THIRD PART: 1D scattering

    for (int ii=0; ii<dim_E; ii++) {
        E_s[ii] =  ii*dE*3.0; 
    }
    for (int ii=0; ii<N; ii++) {
        X[ii] = x0 + ii*h; 
    }
    double complex t;
    double T[dim_E][3]; 
    int ii; 
    for (int jj=0; jj<3; jj++) {
        ii = 0;
        while (ii<dim_E) {
            k[ii] = sqrt(csi_s[jj] * E_s[ii]);
            Phi[0] = cexp(I*k[ii]*x0); 
            Phi[1] = cexp(I*k[ii]*(x0+h)); 
            Numerov (X,Phi,E_s[ii],csi_s[jj]); 
            t = transmissionCoefficient(X,Phi,k[ii]); 
            T[ii][jj] = cabs(t)*cabs(t); 
            ii++; 
        } 
    }
    fileDataT(E_s,T); 
}