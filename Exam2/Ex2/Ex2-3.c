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
#define DIM_t 65000 // number of evolution steps
 
double L = 5.0; 
double x[M], dx; // positions
double complex phi[M], rho_n[M], xi_k[M], psi[M];
double dt = 0.1; // time step
double complex data[M][DIM_t];
double p[DIM_t]; 

double potential (double x) {
    if (x<0.0) {
        return - x * x * (x + 1.0); 
    } else {
        return  2.0 * x * x * (x - 1.0); 
    }
}

double complex Delta_k(int k){
    if (k<= M/2) {
        return I * 2.0 * M_PI * k / L;
    } else {
        return I * 2.0 * M_PI * (k-M) / L;
    }
}

double integrate_data(int pos_c) {
    // perfroms integration with trapezoids method 
    double f_a, f_b, sum;
    f_a = cabs(data[N-1][pos_c]) * cabs(data[N-1][pos_c]); 
    f_b = cabs(data[M-1][pos_c]) * cabs(data[M-1][pos_c]);
    sum = f_a + f_b; 
    for (int jj = N; jj < M-1 ; jj++) {
        sum += 2.0 * cabs(data[jj][pos_c]) * cabs(data[jj][pos_c]);
    }
    return (dx/2.0) * sum;
}

int pos_Max (int dim, double f[dim], int pos_max[dim], double ex1, double ex2) {
    // returns the number of maxima between ex1 and ex2 
    for (int ii=0; ii<dim; ii++) {
        pos_max[ii] = 0; 
    }
    int jj=0;  
    for (int ii=2; ii<dim-1; ii++) {
        if (f[ii] > f[ii-1] && f[ii] > f[ii+1] && f[ii] > ex1 && f[ii] < ex2) {
            pos_max[jj] = ii; 
            jj++; 
        }
    }
    return jj--; 
}

int minimum (int dim, double f[dim], int pos_min[dim], double ex) {
    // finds local minima of f (just for f < ex)
    for (int ii=0; ii<dim; ii++) {
        pos_min[ii] = 0; 
    } 
    int kk=0; 
    for (int ii=2; ii<dim-1; ii++) {
        if (f[ii] < f[ii-1] && f[ii] < f[ii+1] && f[ii] < ex) { 
            pos_min[kk] = ii; 
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
    fd1 = fopen("data_initial.txt","w"); // to show initial situation 
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

    FILE *pf2; 
    pf2 = fopen("evol.txt", "w"); // to show p(t)
    for (int ii = 0; ii < DIM_t; ii++) {
        fprintf(pf2, "%lf \t  %lf \n ", time[ii], p[ii]);
    }
    fclose(pf2);


    // "small" period

    int pos_max[DIM_t];
    int ns = pos_Max (DIM_t, p, pos_max, 0.1, 0.5);
    double t_maxs[ns], f_maxs[ns];
    for (int ii=0; ii<ns; ii++) {
        t_maxs[ii] = time[pos_max[ii]]; // time corresponding to local maxima
        f_maxs[ii] = p[pos_max[ii]]; // local maxima
    }
    FILE *fdMaxs; 
    fdMaxs = fopen("max_small.txt", "w");  
    for (int ii = 0; ii < ns; ii++) {
        fprintf(fdMaxs, "%lf \t  %lf \n ", t_maxs[ii], f_maxs[ii]); 
    }
    fclose(fdMaxs);

    double mean_t = 0.0; // mean value of the small period 
    double diff_t[ns], dist_t[ns];  // distance between local maxima
    int jj=0; 
    for (int ii=0; ii<(ns-1); ii++) {
        diff_t[ii] = t_maxs[ii+1] - t_maxs[ii];
        if (diff_t[ii] > 40.0) {
            dist_t[jj] = diff_t[ii];
            jj++; 
        }
    }
    for (int kk=0; kk<jj--; kk++) {
        mean_t += dist_t[kk]; 
    }
    mean_t = mean_t/jj--; 
    printf("\n Small period: %lf \n", mean_t);

    // "large" period

    double dist_T[5]; 
    int n1 = pos_Max (DIM_t, p, pos_max, 0.183, 0.218);
    double t_max1[n1], f_max1[n1]; 
    for (int ii=0; ii<n1; ii++) {
        t_max1[ii] = time[pos_max[ii]]; // time corresponding to local maxima
        f_max1[ii] = p[pos_max[ii]]; // local maxima
    }
    FILE *fdMax1; 
    fdMax1 = fopen("max1.txt", "w");   
    for (int ii = 0; ii < n1; ii++) {
        fprintf(fdMax1, "%lf \t  %lf \n ", t_max1[ii], f_max1[ii]); 
    }
    fclose(fdMax1);
    int pos_min1[n1]; 
    int dim1 = minimum (n1, f_max1, pos_min1, 0.185); 
    FILE *fdMin1; 
    fdMin1 = fopen("min1.txt", "w");   
    for (int ii = 0; ii < dim1; ii++) {
        fprintf(fdMin1, "%lf \t  %lf \n ", t_max1[pos_min1[ii]], f_max1[pos_min1[ii]]); 
    }
    fclose(fdMin1);
    dist_T[0] = t_max1[pos_min1[1]] - t_max1[pos_min1[0]];  

    int n2 = pos_Max (DIM_t, p, pos_max, 0.210, 0.247); 
    double t_max2[n2], f_max2[n2]; 
    for (int ii=0; ii<n2; ii++) {
        t_max2[ii] = time[pos_max[ii]]; // time corresponding to local maxima
        f_max2[ii] = p[pos_max[ii]]; // local maxima
    }
    FILE *fdMax2; 
    fdMax2 = fopen("max2.txt", "w");  
    for (int ii = 0; ii < n2; ii++) {
        fprintf(fdMax2, "%lf \t  %lf \n ", t_max2[ii], f_max2[ii]); 
    }
    fclose(fdMax2);
    int pos_min2[n1]; 
    int dim2 = minimum (n2, f_max2, pos_min2, 0.212); 
    FILE *fdMin2; 
    fdMin2 = fopen("min2.txt", "w");   
    for (int ii = 0; ii < dim2; ii++) {
        fprintf(fdMin2, "%lf \t  %lf \n ", t_max2[pos_min2[ii]], f_max2[pos_min2[ii]]); 
    }
    fclose(fdMin2);
    dist_T[1] = t_max2[pos_min2[2]] - t_max1[pos_min2[1]];
    dist_T[2] = t_max2[pos_min2[1]] - t_max2[pos_min2[0]]; 
    
    int n3 = pos_Max (DIM_t, p, pos_max, 0.249, 0.5);
    double t_max3[n3], f_max3[n3];
    for (int ii=0; ii<n3; ii++) {
        t_max3[ii] = time[pos_max[ii]]; // time corresponding to local maxima
        f_max3[ii] = p[pos_max[ii]]; // local maxima
    }
    FILE *fdMax3; 
    fdMax3 = fopen("max3.txt", "w");  
    for (int ii = 0; ii < n3; ii++) {
        fprintf(fdMax3, "%lf \t  %lf \n ", t_max3[ii], f_max3[ii]); 
    }
    fclose(fdMax3);

    int pos_min3[n3]; 
    int dim3 = minimum (n3, f_max3, pos_min3, 0.252); 
    FILE *fdMin3; 
    fdMin3 = fopen("min3.txt", "w");   
    for (int ii = 0; ii < dim3; ii++) {
        fprintf(fdMin3, "%lf \t  %lf \n ", t_max3[pos_min3[ii]], f_max3[pos_min3[ii]]); 
    }
    fclose(fdMin3);
    dist_T[3] = t_max3[pos_min3[1]] - t_max3[pos_min3[0]]; 
    dist_T[4] = t_max3[pos_min3[4]] - t_max3[pos_min3[1]]; 

    double mean_T = 0.0; 
    for (int ii=0; ii<5; ii++) {
        mean_T += dist_T[ii]; 
    }
    mean_T = mean_T/5.0; 
    printf(" Large period: %lf \n \n", mean_T);

}

