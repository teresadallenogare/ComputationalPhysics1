#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h> 

#define N 2048 // THIS MUST BE EQUAL TO N IN Ex2-1.c!
#define M 4096 // number of total points 
#define xi 0.015 // adimesional parameter
#define DIM_t 75000 // number of evolution steps
 
double L = 5.0; 
double x[M], dx; // positions 
double complex phi[M], rho_n[M], xi_k[M]; 
double dt = 0.1; // time step
double complex data[M][DIM_t]; 
double p[DIM_t];

double potential (double x) {
    // returns the value of the potential 
    if (x<0.0) {
        return - x * x * (x + 1.0); 
    } else {
        return x * x * (x - 1.0); 
    }
}

double complex Delta_k (int k){
    if (k<= M/2) {
        return I * 2.0 * M_PI * k / L;
    } else {
        return I * 2.0 * M_PI * (k-M) / L;
    }
}

double integrate_data ( int pos_c ) {
    // perfroms integration for the column in pos_c of the matrix 'data'
    double f_a, f_b, sum;
    f_a = cabs(data[N-1][pos_c]) * cabs(data[N-1][pos_c]); 
    f_b = cabs(data[M-1][pos_c]) * cabs(data[M-1][pos_c]);
    sum = f_a + f_b; 
    for (int jj = N; jj < M-1 ; jj++) {
        sum += 2.0 * cabs(data[jj][pos_c]) * cabs(data[jj][pos_c]);
    }
    return (dx/2.0) * sum;
}

double integrate (double f[M]){
    // integral using given values of a function f
    double ris = f[0] + f[M-1];
    for (int ii=1; ii<(M-1); ii++) {
        ris += 2*f[ii];
    }
    return ris * dx/2.0;
}

void pos_MaxMin (int dim, double f[dim], int pos_max[dim], int pos_min[dim], int n[2]) {
    // returns two arrays with the positions of maxima and minima with dimension n[0] and n[1] respectively 
    for (int ii=0; ii<dim; ii++) {
        pos_max[ii] = 0; 
        pos_min[ii] = 0; 
    }
    int jj=0; 
    int kk=0; 
    for (int ii=2; ii<dim-1; ii++) {
        if (f[ii] > f[ii-1] && f[ii] > f[ii+1] && f[ii] > 0.89) { // 0.89 is chosen from the plot
            pos_max[jj] = ii; 
            jj++; 
        }
        if (f[ii] < f[ii-1] && f[ii] < f[ii+1] && f[ii] < 0.11) { // 0.11 is chosen from the plot 
            pos_min[kk] = ii; 
            kk++; 
        }
    }
    n[0] = jj--; 
    n[1] = kk--; 
}

int minimum (int dim, double f[dim], int pos_min[dim]) {
    // finds local minima of f
    for (int ii=0; ii<dim; ii++) {
        pos_min[ii] = 0; 
    } 
    int kk=0; 
    for (int ii=2; ii<dim-1; ii++) {
        if (f[ii] < f[ii-1] && f[ii] < f[ii+1] && f[ii] < 0.93) { // 0.93 is chosen from the plot 
            pos_min[kk] = ii; 
            kk++; 
        }
    }
    return kk--;  
}

int maximum (int dim, double f[dim], int pos_max[dim]) {
    // finds local maxima of f
    for (int ii=0; ii<dim; ii++) {
        pos_max[ii] = 0; 
    }
    int jj=0; 
    int kk=0; 
    for (int ii=2; ii<dim-1; ii++) {
        if (f[ii] > f[ii-1] && f[ii] > f[ii+1] && f[ii] > 0.08 && f[ii] < 0.109) { // 0.109 is chosen from the plot 
            pos_max[kk] = ii; 
            kk++; 
        }
    }
    return kk--;  
}


int main(int argc, const char * argv[]){ 

    // initial vectors

    dx = L/(2*(N-1));
    for (int ii = 0; ii < M; ii++) {
        x[ii] = -L/2.0 + ii * dx;
    }

    double m[N][2]; // matrix to store data from 'wf_gs.txt'
    FILE *fd;
    fd = fopen("wf_gs.txt", "r");
    for (int ii=0; ii<N; ii++) {
        for (int jj=0; jj<2; jj++) {
            fscanf(fd, "%lf", &m[ii][jj]); 
        }
    }
    fclose(fd);
    for (int ii=0; ii<M; ii++) {
        if (ii<N) {
            phi[ii] = m[ii][1]; 
        } else {
            phi[ii] = 0.0; 
        }
    }

    FILE *fd1;     
    fd1 = fopen("data_initial.txt","w");  // to show initial situation 
    for (int kk=0; kk<M; kk++) {
        fprintf(fd1, "%lf \t %lf \t %lf \n", x[kk], creal(phi[kk]), potential(x[kk])); 
    }
    fclose(fd1);
    
    // evolution 

    double time[DIM_t]; 
    for (int ii = 0; ii < M; ii++) {
        rho_n[ii] = cexp(- I * potential(x[ii]) * dt);
        xi_k[ii] = cexp(I * dt * xi * Delta_k(ii) * Delta_k(ii));
    }

    for (int tau = 0; tau < DIM_t; tau++) {
        for (int ii = 0; ii < M; ii++) {
            phi[ii] *= rho_n[ii];
        }
        
        gsl_fft_complex_radix2_forward((double * )phi, 1, M);
        
        for (int ii = 0; ii < M; ii++) {
            phi[ii] *= xi_k[ii];
        }
        
        gsl_fft_complex_radix2_inverse((double * )phi, 1, M);
        
        for (int ii = 0; ii < M; ii++) {
            data[ii][tau] = phi[ii];
        }
        time[tau] = (tau+1) * dt; 
        p[tau] = integrate_data(tau); 
        
    }
    FILE *fd2; 
    fd2 = fopen("evol.txt", "w"); // to show p(t)
    for (int ii = 0; ii < DIM_t; ii++) {
        fprintf(fd2, "%lf \t  %lf \n ", time[ii], p[ii]);
    }
    fclose(fd2);

    // "small" period

    int pos_max[DIM_t], pos_min[DIM_t], n[2];
    pos_MaxMin(DIM_t, p, pos_max, pos_min, n); 
    double t_max[n[0]], f_max[n[0]], t_min[n[1]], f_min[n[1]];
    for (int ii=0; ii<n[0]; ii++) {
        t_max[ii] = time[pos_max[ii]]; // time corresponding to local maxima
        f_max[ii] = p[pos_max[ii]]; // local maxima
    }
    for (int ii=0; ii<n[1]; ii++) {
        t_min[ii] = time[pos_min[ii]]; // time corresponding to local minima
        f_min[ii] = p[pos_min[ii]];  // local minima
    }
    FILE *fd3, *fd4; 
    fd3 = fopen("max.txt", "w");  
    fd4 = fopen("min.txt", "w"); 
    for (int ii = 0; ii < n[0]; ii++) {
        fprintf(fd3, "%lf \t  %lf \n ", t_max[ii], f_max[ii]); 
    }
    for (int ii = 0; ii < n[1]; ii++) {
        fprintf(fd4, "%lf \t  %lf \n ", t_min[ii], f_min[ii]);
    }
    fclose(fd3);
    fclose(fd4);

    int d = n[0] + n[1] -2; 
    double mean_t = 0.0; // mean value of the small period 
    double dist_t[d];  // distance between local maxima/minima
    for (int ii=0; ii<d; ii++) {
        if (ii<n[0]-1) {
            dist_t[ii] = t_max[ii+1] - t_max[ii];
        }
        if (ii>=n[0]-1) {
            dist_t[ii] = t_min[ii-n[0]+2] - t_min[ii-n[0]+1]; 
        }
    }
    int jj = 0; 
    for (int ii=0; ii<d; ii++) {
        if (dist_t[ii] > 150.0) { // some maxima/minima are too close (we exlude them)
            mean_t += dist_t[ii]; 
            jj++; 
        }
    }
    mean_t = mean_t/jj; 
    printf("\n Small period: %lf \n", mean_t); 

    // "large" period
    
    int n_max, p_min[n[0]]; 
    n_max = minimum(n[0], f_max, p_min); 
    double time_max[n_max]; // time correspondig to minima of local maxima 
    for (int ii=0; ii<n_max; ii++) {
        time_max[ii] = t_max[p_min[ii]]; 
    }
    int n_min, p_max[n[1]];  
    n_min = maximum(n[1], f_min, p_max); 
    double time_min[n_min];  // time correspondig to maxima of local minima
    for (int ii=0; ii<n_min; ii++) {
        time_min[ii] = t_min[p_max[ii]]; 
    }
    FILE *fd5, *fd6; 
    fd5 = fopen("max_loc.txt", "w"); // just to check if the points are correct 
    fd6 = fopen("min_loc.txt", "w"); // just to check if the points are correct 
    for (int ii = 0; ii < n_max; ii++) {
        fprintf(fd5, "%lf \t  %lf \n ", time_max[ii], f_max[p_min[ii]]);
    }
    for (int ii = 0; ii < n_min; ii++) {
        fprintf(fd6, "%lf \t  %lf \n ", time_min[ii], f_min[p_max[ii]]);
    }
    fclose(fd5);
    fclose(fd6);

    double dist_T[n_max + n_min -2]; // similar to the previous part 
    for (int ii=0; ii<(n_max + n_min -2); ii++) {
        if (ii<n_max-1) {
            dist_T[ii] = time_max[ii+1] - time_max[ii]; 
        } 
        if (ii>=n_max-1) {
            dist_T[ii] = time_min[ii-n_max+2] - time_min[ii-n_max+1]; 
        } 
    }
    double mean_T = 0.0; 
    for (int ii=0; ii<(n_max + n_min -2); ii++) {
        mean_T += dist_T[ii]; 
    }
    mean_T = mean_T/(n_max + n_min -2); 
    printf(" Large period: %lf \n \n", mean_T);
    
}