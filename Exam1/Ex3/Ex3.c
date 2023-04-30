/* DALLE NOGARE TERESA - 201558, DELBONO ILARIA - 201849  */
/* EXERCISE 3: A FLUID OF REPULSIVE PARTICLES */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#define N 108 /* number of particles */
#define NUM_EV 2000 /* number of evolutions of for the equilibration phase */
#define NUM_SAMP 100 /* number of samples for production phase */
#define NUM_RHO 3 
#define NUM_T 2
#define POS_T 1 /* change in order to have different T */
#define POS_L 2 /* change in order to have different rho */
#define NUM_SLOT 1000 /* number of slots of the g function */


double rho[NUM_RHO] = {0.5, 0.6, 0.7}; /* density - reduced units */
double T[NUM_T] = {0.5, 0.8}; /* temperature - reduced units */
double L[NUM_RHO]; /* length of the box side (depends on rho) */
double cutoff[NUM_RHO]; /* cutoff for the pair interaction calculation */
double X[N][3]; /* positions */
double V[N][3]; /* velocities */
double F[N][3]; /* forces */
double dt= 0.005; /* time step */
double Dt = 1.0; /* evolution time of the system among samples */ 
double g[NUM_SLOT]; /* radial distribution function*/
double dr; /* lenght of a slot interval of g */
double Pex = 0.0; /* excess pressure */

void generate_FCC( unsigned int numPart, double L, double X[N][3]);
void setToZeroV(unsigned int numPart, double V[N]);
void setToZeroM(unsigned int numPart, double M[N][3]);
double repulsiveLJ(double r);
double derRepulsiveLJ(double x, double y, double z);
void calculateForces(double X[N][3], double F[N][3], double L, double cutoff);
void velocity_Verlet(double X[N][3], double V[N][3], double F[N][3], double L, double cutoff);
double kineticEnergySystem(double V[N][3]);
double potentialEnergySystem(double X[N][3], double L, double cutoff);
void rescale_velocities( double V[N][3], double T);
void evolution(double X[N][3], double V[N][3], double F[N][3], double time, double L, double cutoff);
void gFunction(double X[N][3], double g[NUM_SLOT], double L, double cutoff);
double pressureViral(double X[N][3], double L);
double pressureIntegral( double g[NUM_SLOT], double a, double cutoff, double rho );

void fileDataPotential (double time[NUM_EV], double U[NUM_EV],double K_sys[NUM_EV], double E_sys[NUM_EV]);
void fileDatag(double r[NUM_SLOT], double g[NUM_SLOT]);
void fileDataPotential_ev(double time_ev[NUM_SAMP], double U_sys_ev[NUM_SAMP], double K_sys_ev[NUM_SAMP],double E_sys_ev[NUM_SAMP]);
//----------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, const char * argv[]){

  double U_sys[NUM_EV], K_sys[NUM_EV], E_sys[NUM_EV]; /* Energies of the system */
  double time[NUM_EV]; 

/* Definition of box lenght and cutoff */
    for (int kk = 0; kk < NUM_RHO; kk++)
    {
        L[kk] = pow( N/rho[kk], 1.0/3.0);

        if(L[kk] /2.0 < 3.0)
      {
        cutoff[kk] = L[kk]/2.0;
      }
       else
      {
       cutoff[kk] = 3.0;
      }
      printf("rho: %lf\t L: %lf\t cutoff: %lf\n", rho[kk], L[kk], cutoff[kk]);
    }

/* Initial positions, velocities and forces */
generate_FCC(N, L[POS_L], X);
setToZeroM(N, V);
calculateForces(X, F, L[POS_L], cutoff[POS_L]);
  time[0]=0;
  U_sys[0] = potentialEnergySystem(X, L[POS_L], cutoff[POS_L]);
  K_sys[0] = kineticEnergySystem(V);
  E_sys[0] = U_sys[0] + K_sys[0];

/* Equilibration of the system */
for (int kk = 1; kk < NUM_EV; kk++)
{
  velocity_Verlet(X,V,F, L[POS_L], cutoff[POS_L]); 
  U_sys[kk] = potentialEnergySystem(X, L[POS_L], cutoff[POS_L] );
  K_sys[kk] = kineticEnergySystem(V);
  E_sys[kk] = U_sys[kk] + K_sys[kk];
  time[kk] = kk * dt;
  if(kk % 100 == 0){
   rescale_velocities(V, T[POS_T]);
  }
}
fileDataPotential(time, U_sys, K_sys, E_sys);


/* Production: calculation of the mean value */
/* I start to sample after NUM_EV evolutions of the system */
double meanU_sys, meanK_sys;
double Pex_viral = 0.0; /* Excess pressure calucated via Viral theorem */
double U_sys_ev[NUM_SAMP], K_sys_ev[NUM_SAMP], E_sys_ev[NUM_SAMP];
double time_ev[NUM_SAMP];
  for (int kk= 0; kk < NUM_SAMP; kk++)
  {
    evolution(X,V,F,Dt, L[POS_L], cutoff[POS_L]); /* evolution of the system for Dt */
    meanU_sys = meanU_sys + potentialEnergySystem(X, L[POS_L], cutoff[POS_L]);
    meanK_sys = meanK_sys + kineticEnergySystem(V);
    U_sys_ev[kk] = potentialEnergySystem(X, L[POS_L], cutoff[POS_L] );
    K_sys_ev[kk] = kineticEnergySystem(V);
    E_sys_ev[kk] = U_sys_ev[kk] + K_sys_ev[kk];
    time_ev[kk] = NUM_EV* dt+ kk * Dt;
    gFunction(X,g, L[POS_L], cutoff[POS_L]); 
    Pex_viral = pressureViral(X, L[POS_L]);  
  }

  fileDataPotential_ev(time_ev, U_sys_ev, K_sys_ev, E_sys_ev);

  meanU_sys = meanU_sys/ NUM_SAMP;
  meanK_sys = meanK_sys/NUM_SAMP;
  printf("meanU_sys: %lf\n", meanU_sys);
  printf("meanK_sys: %lf\n", meanK_sys);

/* Calculation of the radial distribution function g(s) */
  double r[NUM_SLOT]; /* vector of sphere radii to evaluate g(s) */
  setToZeroV(NUM_SLOT,r);
  dr = cutoff[POS_L]/NUM_SLOT; /* length of the interval */
  for (int s = 0; s < NUM_SLOT; s++)
  {
    r[s] = (s+1)*dr;
  }
  for (int s = 0; s < NUM_SLOT; s++)
  {
    g[s] = g[s]/(4* M_PI* N * NUM_SAMP * rho[POS_L]* r[s]*r[s]); /* final normalization of the pair distribution function
    ( evaluated for different NUM_SAMP configurations )*/
  }
  fileDatag(r,g);

  /* Pressure via viral theorem */
  double V = L[POS_L] * L[POS_L] * L[POS_L];
  Pex_viral = Pex_viral / (6.0 * V * NUM_SAMP);
  printf("Pex_viral: %lf\n", Pex_viral);

  /* Presure via integral */
  double Pex_integral; /* Excess pressure calucated via integral */
  // Pex_integral = pressureIntegral(r, g, 1e-8, cutoff[POS_L], rho[POS_L] );
  Pex_integral = pressureIntegral(g, 1e-8, cutoff[POS_L], rho[POS_L] );
 
  printf("Pex_integral: %lf\n", Pex_integral);

  printf("Accuracy: %lf\n", Pex_viral- Pex_integral);

}
//--------------------------------------------------------------------------------------------------------------------------------------------

void generate_FCC(unsigned int numPart, double L, double X[N][3])

{
  /*
    This function fills the matrix X[N][3] with the coordinates of
    particles in a FCC lattice within a box of side L whose center is the
    center of the coordinate system.
    The box in filled uniformly if N = 4 n^3 with n integer hence for
    N = 4     32    108    256    500    864   1372   2048   2916
    4000   5324 6912   8788  10976  13500  16384 ...
  */

  int i,j,k,m,n,p,c;
  
  /* position within the primary cell */
  const double rFCC[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
                             {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
  double b, rCell[3];
    
  for (c = 1; ; c++)
    if (4*c*c*c >= N)
      break;

  b = L / (double)c;            /* side length of the primary cell */
  p = 0;                        /* particles placed so far */
  for (i = 0; i < c; i++)
    {
      rCell[0] = i;
      for (j = 0; j < c; j++)
        {
          rCell[1] = j;
          for (k = 0; k < c; k++)
            {
              rCell[2] = k;
              for (m = 0; m < 4; m++) /* 4 particles in cell */
                if (p < N)
                  {
                    
                    /* add the com to each bead, and project to the real cell */
                    for(n=0;n<3;n++)
                      {
                        X[p][n] = b * (rCell[n] + rFCC[m][n]);
                        X[p][n] -= L * rint(X[p][n]/L);
                      }
                ++p;
                  }
            }
        }
    }
}

void setToZeroV(unsigned int numPart, double V[N]){
  for (int kk=0; kk<numPart; kk++) {
            V[kk]=0;
        }
}

void setToZeroM(unsigned int numPart, double M[N][3]){
    for (int jj=0; jj<numPart; jj++) {
        for (int kk=0; kk<3; kk++) {
            M[jj][kk]=0;
        }
    
    }
}

double repulsiveLJ(double r){
    return 4*pow(r,-12.0);
}

double derRepulsiveLJ(double x, double y, double z){
    return -48* x* pow( x*x + y*y + z*z ,-7.0 );
}

void calculateForces(double X[N][3], double F[N][3], double L, double cutoff){
    /* zeroes out a[N][3], loop on all the pairs, calculates the distance,
   checks the cutoff, update forces */
   double dx = 0, dy = 0, dz = 0, r2 = 0;
   setToZeroM(N, F);
   for (int kk = 1; kk < N; kk++) /* cicle on all particles */
   {
       for (int jj = 0; jj < kk; jj++)
       {   /* calculates the smallest distance between two particles in every direction */
               dx = X[kk][0] - X[jj][0];    dy = X[kk][1] - X[jj][1];    dz = X[kk][2] - X[jj][2];
               dx = dx - L*rint(dx/L);  dy = dy - L*rint(dy/L);  dz = dz - L*rint(dz/L);
               r2 = dx*dx + dy*dy + dz*dz;

               if (r2<cutoff* cutoff)
               { /* controls if the distance is into the sphere with radius: cutoff */
                 double f;
                 f = derRepulsiveLJ(dx, dy, dz);
                   F[kk][0] = F[kk][0] - f;
                   F[jj][0] = F[jj][0] + f;
                 f = derRepulsiveLJ(dy, dz, dx);
                   F[kk][1] = F[kk][1] - f;
                   F[jj] [1] = F[jj][1] + f;
                 f = derRepulsiveLJ(dz, dx, dy);
                   F[kk][2] = F[kk][2] - f;
                   F[jj][2] = F[jj] [2] +f;
               }
           
       }
       
   }

}

void velocity_Verlet(double X[N][3], double V[N][3], double F[N][3], double L, double cutoff){
    /* performs one step of the velocity Verlet algorithm, calling the previous
function to calculate the forces */
  for (int kk = 0; kk<N; kk++) {
    for (int jj = 0; jj < 3; jj++)
    {
      V[kk][jj] += F[kk][jj] * dt/2.0;
      X[kk][jj] += V[kk][jj] * dt; 
    }
  }
  calculateForces(X, F, L, cutoff);

  for (int kk= 0 ; kk<N; kk++) {
    for(int jj= 0; jj < 3; jj++){
      V[kk][jj] += F[kk][jj] * dt/2.0;
     }
  }
}

double kineticEnergySystem(double V[N][3]){
  double K_tot = 0, v[N];
  /* initialization of the velocity vector */
   setToZeroV(N, v);
  for (int kk = 0; kk < N; kk++)
  {
    v[kk] = sqrt(pow(V[kk][0], 2) + pow(V[kk][1], 2) + pow(V[kk][2], 2));
    K_tot = K_tot + 0.5* pow(v[kk],2);
  }
  return K_tot;
}

double potentialEnergySystem(double X[N][3], double L, double cutoff){
  double U = 0;
  double dx = 0, dy = 0, dz = 0, r = 0;
  for (int kk = 0; kk < N; kk++)
  {
    for (int jj = 0; jj < kk; jj++)
    {
      /* calculates the smallest distance between two particles in every direction */
         dx = X[kk][0] - X[jj][0];    dy = X[kk][1] - X[jj][1];    dz = X[kk][2] - X[jj][2];
         dx = dx - L*rint(dx/L);  dy = dy - L*rint(dy/L);  dz = dz - L*rint(dz/L);
         r = sqrt(dx*dx + dy*dy + dz*dz);
        if (r < cutoff ) {
          U = U + repulsiveLJ(r);
      }
    }  
  }
  return U;
}

void rescale_velocities( double V[N][3], double T){
  /* rescale velocities to match the assigned temperatures */
  double K_ist = 0, T_ist = 0, alpha = 0;
  K_ist = kineticEnergySystem(V); /* total kinetic energy of the sy stem at a specific t*/
  T_ist = 2.0*K_ist/(3.0*N);
  alpha = sqrt(T/T_ist); 
  for (int kk = 0; kk < N; kk++){
    for (int jj = 0; jj < 3.0; jj++)
    {
     V[kk][jj] = alpha * V[kk][jj];
    }
  }
}

void evolution(double X[N][3], double V[N][3], double F[N][3], double time, double L, double cutoff){
  /* cicles velocityVerlet until kk*dt= time */ 
   for (int kk=0; kk*dt < time; kk++)
  {   
     velocity_Verlet(X,V,F,L,cutoff);
  }
  
}

void gFunction(double X[N][3], double g[NUM_SLOT], double L, double cutoff){
  double dx, dy, dz, dist;
  int s; 
  int trovato = 0; 
  //setToZeroV(S,g);
  dr = cutoff/NUM_SLOT ; /* fisso dr come la lunghezza dell'intervallino entro cui voglio calcolare il numero di particelle:
               lo definisco come la lunghezza totale cutoff fratto il numero di slot possibili */
  for (int kk = 0; kk < N; kk++){
    for (int jj = 0; jj < N-1; jj++){ 
      s = 0; /* slot in cui sono */
      trovato = 0;
      if (jj != kk){
        dx = X[kk][0] - X[jj][0];    dy = X[kk][1] - X[jj][1];    dz = X[kk][2] - X[jj][2];
        dx = dx - L*rint(dx/L);  dy = dy - L*rint(dy/L);  dz = dz - L*rint(dz/L);
        dist = sqrt(dx*dx + dy*dy + dz*dz); /* distance between particle kk and jj */
        if (dist < cutoff){
          while ( s < NUM_SLOT && trovato == 0){
             if (dist < (s+1)*dr && dist > s*dr){
              g[s] = g[s] + 1/dr ;
              trovato = 1;
              }
             s++;
            }
          }
        }
      }
    }

 }

double pressureViral(double X[N][3], double L){
double dx, dy, dz;
double r;
double f;
for (int kk = 0; kk < N; kk++)
{
  for (int jj = 0; jj < N; jj++)
  {
    if( jj != kk){
        /* calculates the smallest distance between two particles in every direction */
         dx = X[kk][0] - X[jj][0];    dy = X[kk][1] - X[jj][1];    dz = X[kk][2] - X[jj][2];
         dx = dx - L*rint(dx/L);  dy = dy - L*rint(dy/L);  dz = dz - L*rint(dz/L);
         r = sqrt(dx*dx + dy*dy + dz*dz);
            f = 48.0 * pow(r , -13.0); // dalla derivata esce un meno, poi uso f = -grad(v)
          Pex += r* f;
    }
  }
}
  return Pex;
}


double pressureIntegral(double g[NUM_SLOT], double a, double cutoff, double rho ){
 
    double h, h_small, f_a, f_b;
    double trapez_sum= 0;
    double Pex;
    double r_kk = 0;
    h = cutoff/NUM_SLOT ; // lunghezza intervallino
    h_small = h/100.0;
    f_a = - 48.0 * pow( a, -10.0) * g[0];
    f_b = - 48.0 * pow( cutoff, -10.0) * g[NUM_SLOT-1];

 for (int jj = 1; jj < NUM_SLOT; jj++)
 {
    for (int kk = 1; kk < 100 ; kk++) {
       r_kk = h*(jj-1) + kk * h_small;
       trapez_sum += - 48.0 * pow( r_kk, -10.0) * g[jj];
    }
   
 }
  trapez_sum = h_small * (trapez_sum + 1/2 * (f_a+f_b));

  Pex = -2.0* M_PI / 3.0* rho * rho *trapez_sum;


    return Pex;
}


void fileDataPotential (double time[NUM_EV], double U_sys[NUM_EV], double K_sys[NUM_EV], double E_sys[NUM_EV]){

FILE *pf;
pf = fopen("dataU.txt", "w");
fprintf(pf, "# File with the system potential over time \n");
fprintf(pf, "#t[kk]\t\t U[kk] \t\t K[kk]\t\t E[kk]\n");
  for (int kk=0; kk < NUM_EV; kk++) 
  {
    fprintf(pf, "%lf \t %lf \t %lf\t %lf\n", time[kk], U_sys[kk], K_sys[kk], E_sys[kk]) ; 
  }
fclose(pf);
}

void fileDatag(double r[NUM_SLOT], double g[NUM_SLOT]){
  FILE *pf;
pf = fopen("dataG.txt", "w");
fprintf(pf, "# File with the g function of the system \n");
fprintf(pf, "#r[kk]\t\t g[kk] \n");
  for (int s = 0; s < NUM_SLOT; s++) 
  {
    fprintf(pf, "%lf \t %lf\n", r[s], g[s]) ; 
  }
fclose(pf);

}

void fileDataPotential_ev (double time_ev[NUM_SAMP], double U_sys_ev[NUM_SAMP], double K_sys_ev[NUM_SAMP], double E_sys_ev[NUM_SAMP]){

FILE *pf;
pf = fopen("dataU_ev.txt", "w");
fprintf(pf, "# File with the system potential over time \n");
fprintf(pf, "#t[kk]\t\t U[kk] \t\t K[kk]\t\t E[kk]\n");
  for (int kk=0; kk < NUM_SAMP; kk++) 
  {
    fprintf(pf, "%lf \t %lf \t %lf\t %lf\n", time_ev[kk], U_sys_ev[kk], K_sys_ev[kk], E_sys_ev[kk]) ; 
  }
fclose(pf);
}
