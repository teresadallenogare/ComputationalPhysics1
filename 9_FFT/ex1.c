#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define N 256
double f_1[2*N], f_2[2*N];
double x[N]; 
double L = 12.0; 
double R = 8.0; 

void fileData (double x[N], double f_1[2*N], double f_2[2*N]) {
  FILE *fd;     
  fd = fopen("data_f.txt","w"); 
  for (int kk=0; kk<N; kk++) {
      fprintf(fd, "%lf \t %lf \t %lf \t %lf \t %lf \n", x[kk], f_1[2*kk], f_1[2*kk+1], f_2[2*kk], f_2[2*kk+1]); 
  }
  fclose(fd);
}

int main(int argc, const char * argv[]){
double dx = L/N; 
double sigma; 
for(int n=0; n<N; n++) {
  x[n] = n * dx;
  sigma = 1.0;
  f_1[2*n] = exp(-pow((x[n] - R/2.0),2.0)/(2.0*sigma*sigma));
  sigma = 3.0;
  f_1[2*n] += 0.25*exp(-pow((x[n] - 2.0*R/3.0),2.0)/(2.0*sigma*sigma));
  f_1[2*n+1] = 0.0;
}
gsl_fft_complex_radix2_forward (f_1, 1, N);

double c, s, zre, zim;
for(int k=0;k<N;k++) {
  c = cos((2.0*M_PI/N)*k*64.0);
  s = sin((2.0*M_PI/N)*k*64.0);
  zre = f_1[2*k];
  zim = f_1[2*k+1];
  f_2[2*k]   = zre * c - zim * s;
  f_2[2*k+1] = zre * s + zim * c;
}
gsl_fft_complex_radix2_inverse (f_2, 1, N);

for(int n=0; n<N; n++) {
  x[n] = n * dx;
  sigma = 1.0;
  f_1[2*n] = exp(-pow((x[n] - R/2.0),2.0)/(2.0*sigma*sigma));
  sigma = 3.0;
  f_1[2*n] += 0.25*exp(-pow((x[n] - 2.0*R/3.0),2.0)/(2.0*sigma*sigma));
  f_1[2*n+1] = 0.0;
}
fileData(x,f_1,f_2); 

}