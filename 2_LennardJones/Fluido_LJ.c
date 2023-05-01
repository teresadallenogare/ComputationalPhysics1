
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#define N 108   /* number of particles */
#define NUM_EV 2000 /* number of evolutions of for the equilibration phase */
#define NUM_SAMP 500 /* number of samples for production phase */
    double rho= 0.5;     /* density - reduced units */
    double T = 1.3;       /* temperature - reduced units (1 means I have a liquid) */
    double X[N][3]; /* positions */
    double V[N][3]; /* velocities */
    double F[N][3]; /* forces */
    double a[N][3]; /* accelerations = forces */
    double L;      /* box length */
    double cutoff;  /* cutoff for the pair interaction calculation: con particelle lontane forza zero: considero solo quelle relativamente più vicine: assumo che forza sia zero circa 3sigma */
    double dt= 0.005;      /* time step */
    double Dt = 1; /* evolution time of the system among samples */    
#define S 100
    double g[S];
    double dr;
    double T_ist; /* istantaneal temperature*/

void setToZeroM(unsigned int numPart, double M[N][3]);
void setToZeroV(unsigned int numPart, double V[N]);

void generate_FCC( unsigned int numPart, double L, double X[N][3]);
void calculateForces(double X[N][3], double F[N][3]);
void velocity_Verlet(double X[N][3], double V[N][3], double F[N][3]);
double kineticEnergy(double V[N][3]);
void rescale_velocities( double V[N][3]);
void evolution(double X[N][3], double V[N][3], double F[N][3], double time);

double LJpotential(double r);
double derLennardJones ( double x, double y, double z);

double potentialEnergySystem(double X[N][3]);
void gFunction(double X[N][3], double g[S]);

void fileData ( double X[N][3], double V[N][3], double F[N][3]);
void fileDataPotential (double time[NUM_EV], double U[NUM_EV],double K_sys[NUM_EV], double E_sys[NUM_EV]);
void fileDatag(double r[S], double g[S]);
//----------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, const char * argv[]){
  
  double U_sys[NUM_EV], K_sys[NUM_EV], time[NUM_EV];
  double E_sys[NUM_EV]; /* total energy of the system */
  double meanU_sys = 0, meanK_sys = 0;
  double X0[N][3];
 
  L = pow(N/rho, 1.0/3.0);
  printf("L: %lf\n", L);
  //cutoff = (L/2.0 < 3.0 ? L/2.0 : 3.0);
  if(L /2.0 < 3.0)
   {
     cutoff = L/2.0;
   }
 else
   {
     cutoff = 3.0;
   }

// initial positions, velocities and forces
generate_FCC(N, L, X);

FILE *ppf;
ppf = fopen("dataCI.txt", "w");
for (int i = 0; i < N; i++)
{
    fprintf(ppf,"%lf \t %lf \t %lf\n", X[i][0], X[i][1], X[i][2]);
}
fclose(ppf);

setToZeroM(N, V);
calculateForces(X,F);
  time[0]=0;
  U_sys[0] = potentialEnergySystem(X);
  K_sys[0] = kineticEnergy(V);
  E_sys[0] = U_sys[0] + K_sys[0];

// Equilibration of the system
for (int kk = 1; kk < NUM_EV; kk++)
{
  velocity_Verlet(X,V,F); 
  U_sys[kk] = potentialEnergySystem(X);
  K_sys[kk] = kineticEnergy(V);
  E_sys[kk] = U_sys[kk] + K_sys[kk];
  time[kk] = kk * 0.005;
  if(kk % 100 == 0){
   rescale_velocities(V);
  }
}
fileDataPotential(time, U_sys, K_sys, E_sys);

// Production: calculation of the mean value
// I start to sample after NUM_EV evolutions of the system
  for (int kk= 0; kk < NUM_SAMP; kk++)
  {
    evolution(X,V,F,Dt); /* faccio evolvere sistema fino per un intervallo Dt e poi lo blocco*/
    meanU_sys = meanU_sys + potentialEnergySystem(X);
    meanK_sys = meanK_sys + kineticEnergy(V);
    gFunction(X,g);
  }
  meanU_sys = meanU_sys/ NUM_SAMP;
  meanK_sys = meanK_sys/NUM_SAMP;
  printf("meanU_sys: %lf\n", meanU_sys);
  printf("meanK_sys: %lf\n", meanK_sys);
 
 double r[S];
  setToZeroV(S,r);
  dr = cutoff/S;
  for (int s = 0; s < S; s++)
  {
    r[s] = s*dr;
  }
  for (int s = 0; s < S; s++)
  {
    g[s] = g[s]/(4* M_PI* N * NUM_SAMP * rho* r[s]*r[s]);
  }
  fileDatag(r,g);

}
//----------------------------------------------------------------------------------------------------------------------------------------------------------

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

void calculateForces(double X[N][3], double F[N][3]){
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
                 f = derLennardJones(dx, dy, dz);
                   F[kk][0] = F[kk][0] - f;
                   F[jj][0] = F[jj][0] + f;
                 f = derLennardJones(dy, dz, dx);
                   F[kk][1] = F[kk][1] - f;
                   F[jj] [1] = F[jj][1] + f;
                 f = derLennardJones(dz, dx, dy);
                   F[kk][2] = F[kk][2] - f;
                   F[jj][2] = F[jj] [2] +f;
               }
           
       }
       
   }

}

double derLennardJones ( double x, double y, double z) {
    // correct sign
    return -24*x*( 2*pow( (pow(x,2)+pow(y,2)+ pow(z, 2)),-7) -  pow( (pow(x,2)+pow(y,2)+pow(z, 2)),-4) );
}

void velocity_Verlet(double X[N][3], double V[N][3], double F[N][3]){
    /* performs one step of the velocity Verlet algorithm, calling the previous
function to calculate the forces */
  for (int kk = 0; kk<N; kk++) {
    for (int jj = 0; jj < 3; jj++)
    {
      V[kk][jj] += F[kk][jj] * dt/2.0;
      X[kk][jj] += V[kk][jj] * dt; 
    }
  }
  calculateForces(X, F);

  for (int kk= 0 ; kk<N; kk++) {
    for(int jj= 0; jj < 3; jj++){
      V[kk][jj] += F[kk][jj] * dt/2.0;
     }
  }
}

double kineticEnergy(double V[N][3]){
  double K_tot = 0, v[N];
  // initialization of the velocity vector
   setToZeroV(N, v);
  for (int kk = 0; kk < N; kk++)
  {
    v[kk] = sqrt(pow(V[kk][0], 2) + pow(V[kk][1], 2) + pow(V[kk][2], 2));
    K_tot = K_tot + 0.5* pow(v[kk],2);
  }
  return K_tot;
}

void rescale_velocities( double V[N][3]){
  /* rescale velocities to match the assigned temperatures */
  double K_ist = 0, T_ist = 0, alpha = 0;
  K_ist = kineticEnergy(V); /* total kinetic energy of the sy stem at a specific t*/
  T_ist = 2.0*K_ist/(3.0*N);
  alpha = sqrt(T/T_ist); 
  for (int kk = 0; kk < N; kk++){
    for (int jj = 0; jj < 3.0; jj++)
    {
     V[kk][jj] = alpha * V[kk][jj];
    }
  }
}

void evolution(double X[N][3], double V[N][3], double F[N][3], double time){
  /* cycles velocityVerlet until kk*dt= time */ 
   for (int kk=0; kk*dt < time; kk++)
  {   
     velocity_Verlet(X,V,F);
  }
  
}

double LJpotential(double r){
  return 4*(pow(r,-12) - pow(r,-6));
}

double potentialEnergySystem(double X[N][3]){
  double U = 0;
  double dx = 0, dy = 0, dz = 0, r = 0;
  //double r[N-1];
  //double dx[N-1], dy[N-1], dz[N-1];
  //setToZeroV(N-1, dx); setToZeroV(N-1, dy);  setToZeroV(N-1, dz);
  //setToZeroV(N-1, r);
  for (int kk = 0; kk < N; kk++)
  {
    for (int jj = 0; jj < kk; jj++)
    {
      /* calculates the smallest distance between two particles in every direction */
         dx = X[kk][0] - X[jj][0];    dy = X[kk][1] - X[jj][1];    dz = X[kk][2] - X[jj][2];
         dx = dx - L*rint(dx/L);  dy = dy - L*rint(dy/L);  dz = dz - L*rint(dz/L);
         r = sqrt(dx*dx + dy*dy + dz*dz);
        if (r < cutoff ) {
          U = U + LJpotential(r);
      }
    }
    
  }
  return U;
}

void gFunction(double X[N][3], double g[S]){
  double dx, dy, dz, dist;
  int s; 
  int trovato = 0; 
  //setToZeroV(S,g);
  dr = cutoff/S ; /* fisso dr come la lunghezza dell'intervallino entro cui voglio calcolare il numero di particelle:
               lo definisco come la lunghezza totale cutoff fratto il numero di slot possibili S */
  for (int kk = 0; kk < N; kk++){
    for (int jj = 0; jj < N-1; jj++){ 
      s = 0; /* slot in cui sono */
      trovato = 0;
      if (jj != kk){
        dx = X[kk][0] - X[jj][0];    dy = X[kk][1] - X[jj][1];    dz = X[kk][2] - X[jj][2];
        dx = dx - L*rint(dx/L);  dy = dy - L*rint(dy/L);  dz = dz - L*rint(dz/L);
        dist = sqrt(dx*dx + dy*dy + dz*dz); /* distance between particle kk and jj */
        if (dist < cutoff){
          while ( s < S && trovato == 0){
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


void fileData ( double X[N][3], double V[N][3], double F[N][3]) {

  FILE *fdX;     
  fdX = fopen("dataX.txt","w"); 
  fprintf(fdX, "# File with positions of N particles. \n");
  fprintf(fdX, "#x\t\t y\t\t z \n");
  for (int k=0; k<N; k++) 
  {
    fprintf(fdX, "%lf \t %lf \t %lf \n", X[k][0], X[k][1], X[k][2]) ; 
  }
  fclose(fdX);
  FILE *fdV;    
  fdV = fopen("dataV.txt","w");
  fprintf(fdV, "# File with velocities of N particles. \n");
  fprintf(fdV, "#vx\t\t vy\t\t vz \n");
  for (int k=0; k<N; k++) 
  {
    fprintf(fdV, "%lf \t %lf \t %lf \n", V[k][0], V[k][1], V[k][2]) ; 
  }
  fclose(fdV); 
  FILE *fdF;     
  fdF = fopen("dataF.txt","w"); 
  fprintf(fdF, "# File with forces of N particles. \n");
  fprintf(fdF, "#x\t\t y\t\t z \n");
  for (int k=0; k<N; k++) 
  {
    fprintf(fdF, "%lf \t %lf \t %lf \n", F[k][0], F[k][1], F[k][2]) ; 
  }
  fclose(fdF);
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
void fileDatag(double r[S], double g[S]){
  FILE *pf;
pf = fopen("dataG.txt", "w");
fprintf(pf, "# File with the g function of the system \n");
fprintf(pf, "#r[kk]\t\t g[kk] \n");
  for (int s = 0; s < S; s++) 
  {
    fprintf(pf, "%lf \t %lf\n", r[s], g[s]) ; 
  }
fclose(pf);

}
