#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define N 50000
#define dimE 500

double h = 0.01; 
double k;  
double d = 1.0; 
double V0; 
double  E[dimE]; 
double x0 = -5.0; 
double X[N]; 
double complex Phi[N]; 

double rectangularBarrier (const double x) {
    if (x>-0.5 && x<0.5) {
        return V0; 
    } else {
        return 0.0; 
    }
}
double gaussianBarrier (const double x) {
    return V0*exp(-x*x/2); 
}
double F (const double x, const double E, double (*f)(const double)) {
    return 2*(f(x) - E); 
}
void setToZero(unsigned int r, unsigned int c, double v[r][c]) {
  for (int kk=0; kk<r; kk++) {
    for (int jj=0; jj<c; jj++) {
        v[kk][jj] = 0.0; 
    }
  }
}
void NumerovRectangular (double X[N], double complex Phi[N], double E) { 
    int kk=2; 
    while (kk<N) {
        double complex num1, num2, den; 
        num1 = (2.0 + (5.0/6.0)*h*h*F(X[kk-1],E,&rectangularBarrier)) * Phi[kk-1]; 
        num2 = (1 - h*h/12.0*F(X[kk-2],E,&rectangularBarrier))  * Phi[kk-2]; 
        den = (1 - (h*h/12.0)*F(X[kk],E,&rectangularBarrier));  
        Phi[kk] = (num1 - num2)/den; 
        kk++; 
    }
} 
void NumerovGaussian (double X[N], double complex Phi[N], double E) { 
    int kk=2; 
    while (kk<N) {
        double complex num1, num2, den; 
        num1 = (2.0 + (5.0/6.0)*h*h*F(X[kk-1],E,&gaussianBarrier)) * Phi[kk-1]; 
        num2 = (1 - h*h/12.0*F(X[kk-2],E,&gaussianBarrier))  * Phi[kk-2]; 
        den = (1 - (h*h/12.0)*F(X[kk],E,&gaussianBarrier));  
        Phi[kk] = (num1 - num2)/den; 
        kk++; 
    }
} 
void fileDataPhi (double X[N], double complex Phi[N]) {
  FILE *fd;     
  fd = fopen("dataPhi_ilaria.txt","w"); 
  for (int kk=0; kk<N; kk++) {
      fprintf(fd, "%lf \t %f \t %lf \t %lf \n", X[kk], creal(Phi[kk]), cimag(Phi[kk]), cabs(Phi[kk])); 
  }
  fclose(fd);
}
void fileDataT (double E[dimE], double T[dimE]) {
  FILE *fd;     
  fd = fopen("dataT_ilaria.txt","w"); 
  for (int kk=0; kk<dimE; kk++) {
      fprintf(fd, "%lf \t %f \n", kk*0.01, T[kk]); 
  }
  fclose(fd);
}
double complex transmissionCoefficient (double const X[N], double complex const Phi[N]) {
    double complex num, den;
    double x1, x2;
    x1 = X[N-11]; 
    x2 = X[N-1];  
    num = cexp(I*k*(x1-x2)) - cexp(-I*k*(x1-x2)); 
    den = cexp(-I*k*x2)*Phi[N-11] - cexp(-I*k*x1)*Phi[N-1];
    return num/den; 
}

int main(int argc, const char * argv[]){

    V0 = d*d/2.0; 
    for (int kk=0; kk<dimE; kk++) {
        E[kk] = V0*kk*0.01; 
    }
    for (int ii=0; ii<N; ii++) {
        X[ii] = x0 + ii*h; 
    }

    /* Rectangular barrier */
    double complex t;
    double T_Rec[dimE]; 
    int ii=0; 
    while (ii<dimE) {
        k = sqrt(2*E[ii]);
        Phi[0] = cexp(I*k*x0); 
        Phi[1] = cexp(I*k*(x0+h)); 
        NumerovRectangular (X,Phi,E[ii]); 
        t = transmissionCoefficient(X,Phi); 
        T_Rec[ii] = cabs(t)*cabs(t); 
        printf("Trec: %lf", T_Rec[ii]);
        ii++; 

    } 


    /* Gaussian barrier */
    // double complex t;
    double T_Gau[dimE]; 
     ii=0; 
    while (ii<dimE) {
        k = sqrt(2*E[ii]);
        Phi[0] = cexp(I*k*x0); 
        Phi[1] = cexp(I*k*(x0+h)); 
        NumerovGaussian (X,Phi,E[ii]); 
        t = transmissionCoefficient(X,Phi); 
        T_Gau[ii] = cabs(t)*cabs(t); 
        ii++; 
    } 

}