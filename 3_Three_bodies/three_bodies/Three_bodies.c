/* Consider three bodies of equal mass m = 1 iteracting gravitationally and use a system of units where G = 1. 
Plot trajectories in xy plane. */
/* Reduced units: G = 1, m = 1*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#define NUM_EV 80000
#define NUM_BODIES 3

double dt = 1e-4;
double alpha1 = 0.306893; /* 0.347111; 0.306893; 0.464445; */
double alpha2 =  0.125507; /* 0.532728; 0.125507; 0.396060; */


void initialConditions(double X_sys[NUM_BODIES][2], double V_sys[NUM_BODIES][2]);
void setToZeroV(int dim, double V[dim]);
void setToZeroM(int dim, double M[dim][2]);
void velocityVerlet(double X_sys[NUM_BODIES][2], double V_sys[NUM_BODIES][2], double F_sys[NUM_BODIES][2]);
void calculateForces(double X_sys[NUM_BODIES][2], double F_sys[NUM_BODIES][2]);
double derPot_X(double Xa[2], double Xb[2], double Xc[2]);
double derPot_Y(double Xa[2], double Xb[2], double Xc[2]);

void fileDataX( double X1[NUM_BODIES][2], double X2[NUM_BODIES][2], double X3[NUM_BODIES][2], double time[NUM_EV]);
//----------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, const char * argv[]){
double X_sys[NUM_BODIES][2]; /* Positions of each body: row 0: (x1, y1); 1: (x2, y2); 2: (x3, y3)*/
double V_sys[NUM_BODIES][2]; /* Velocities of each body: row 0: (vx1, vy1); 1: (vx2, vy2); 2: (vx3, vy3)*/
double F_sys[NUM_BODIES][2]; /* Forces of each body: row 0: (Fx1, Fy1); 1: (Fx2, Fy2); 2: (Fx3, Fy3)*/
double X1[NUM_EV][2], X2[NUM_EV][2], X3[NUM_EV][2]; /* Positions of each body: each raw is a new evolution */
double V1[NUM_EV][2], V2[NUM_EV][2], V3[NUM_EV][2]; /* Velocities of each body: each raw is a new evolution */
double time[NUM_EV];

/* Initial conditions */
initialConditions (X_sys, V_sys); 
for (int jj = 0; jj < 2; jj++)
{
    X1[0][jj] = X_sys[0][jj];  X2[0][jj] = X_sys[1][jj];  X3[0][jj] = X_sys[2][jj];
    V1[0][jj] = V_sys[0][jj];  V2[0][jj] = V_sys[1][jj];  V3[0][jj] = V_sys[2][jj];
}
calculateForces(X_sys, F_sys); /* Initialize forces to give to velocityverlet */
time[0] = 0;

/* Evolution of each body */
for (int ev = 1; ev < NUM_EV; ev++)
{
 velocityVerlet(X_sys, V_sys, F_sys);
  for (int jj = 0; jj < 2; jj++)
  {
      X1[ev][jj] = X_sys[0][jj];
      X2[ev][jj] = X_sys[1][jj];
      X3[ev][jj] = X_sys[2][jj];
  }
  time[ev] = ev * dt;
}
  fileDataX(X1, X2, X3, time);


}
//----------------------------------------------------------------------------------------------------------------------------------------------------------

void initialConditions(double X_sys[NUM_BODIES][2], double V_sys[NUM_BODIES][2]){
setToZeroM(NUM_BODIES, X_sys);   setToZeroM(NUM_BODIES, V_sys);
// body 1
X_sys[0][0] = -1;  X_sys[0][1] = 0;
V_sys[0][0] = alpha1; V_sys[0][1] = alpha2;
// body 2
X_sys[1][0] = 1; X_sys[1][1] = 0;
V_sys[1][0] = alpha1; V_sys[1][1] = alpha2;
// body 3
X_sys[2][0] = 0;  X_sys[2][1] = 0;
V_sys[2][0] = -2*alpha1;  V_sys[2][1] = -2*alpha2;

}

void setToZeroV(int dim, double V[dim]){
  for (int kk=0; kk < dim; kk++) {
            V[kk]=0;
        }
}
void setToZeroM(int dim, double M[dim][2]){
     for (int kk=0; kk < dim; kk++) {
           for (int jj = 0; jj < 2; jj++)
           {
               M[kk][jj] = 0;
           }
           
        }
}

double derPot_X(double Xa[2], double Xb[2], double Xc[2]){
    /* calculates the derivative of the gravitational potential with respect to x */
double pr1, pr2, Dpot;
 pr1 = pow( (Xa[0]-Xb[0])*(Xa[0]-Xb[0]) + (Xa[1]-Xb[1])*(Xa[1]-Xb[1]) , -3.0/2.0 ) * (Xa[0]-Xb[0]);
 pr2 = pow( (Xc[0]-Xa[0])*(Xc[0]-Xa[0]) + (Xc[1]-Xa[1])*(Xc[1]-Xa[1]) , -3.0/2.0 ) * (Xc[0]-Xa[0]);
 Dpot = pr1 - pr2;
return Dpot;
}
double derPot_Y(double Xa[2], double Xb[2], double Xc[2]){
    /* calculates the derivative of the gravitational potential with respect to y */
double pr1, pr2, Dpot;
 pr1 = pow( (Xa[0]-Xb[0])*(Xa[0]-Xb[0]) + (Xa[1]-Xb[1])*(Xa[1]-Xb[1]) , -3.0/2.0 ) * (Xa[1]-Xb[1]);
 pr2 = pow( (Xc[0]-Xa[0])*(Xc[0]-Xa[0]) + (Xc[1]-Xa[1])*(Xc[1]-Xa[1]) , -3.0/2.0 ) * (Xc[1]-Xa[1]);
 Dpot = pr1 - pr2;
return Dpot;
}

void calculateForces(double X_sys[NUM_BODIES][2], double F_sys[NUM_BODIES][2]){
setToZeroM(NUM_BODIES, F_sys);

double Xa[2], Xb[2], Xc[2];
for (int jj = 0; jj < 2; jj++)
{
    Xa[jj] = X_sys[0][jj];  Xb[jj] = X_sys[1][jj];  Xc[jj] = X_sys[2][jj];
}
// body 1
     F_sys[0][0] = - derPot_X(Xa, Xb, Xc); 
     F_sys[0][1] = - derPot_Y(Xa, Xb, Xc);
 // body 2
     F_sys[1][0] = - derPot_X(Xb, Xc, Xa); 
     F_sys[1][1] = - derPot_Y(Xb, Xc, Xa);
 // body 3
     F_sys[2][0] = - derPot_X(Xc, Xa, Xb); 
     F_sys[2][1] = - derPot_Y(Xc, Xa, Xb);
}

void velocityVerlet(double X[NUM_BODIES][2], double V[NUM_BODIES][2], double F[NUM_BODIES][2]){
    /* performs one step of the velocity Verlet algorithm */
 for (int kk = 0; kk < NUM_BODIES; kk++)
 {  
     for (int jj = 0; jj < 2; jj++)
    {
      V[kk][jj] += F[kk][jj] * dt/2.0;
      X[kk][jj] += V[kk][jj] * dt; 
    }
  }
  calculateForces(X, F);

  for (int kk = 0 ; kk < NUM_BODIES; kk++) {
    for(int jj = 0; jj < 2; jj++){
      V[kk][jj] += F[kk][jj] * dt/2.0;
     }
    
 }
}

void fileDataX( double X1[NUM_BODIES][2], double X2[NUM_BODIES][2], double X3[NUM_BODIES][2], double time[NUM_EV]){
  FILE *pf;
pf = fopen("dataX.txt", "w");
fprintf(pf, "# File with the trajectories of the 3 bodies \n");
fprintf(pf, "# x1[kk] \t\t y1[kk] \t\t x2[kk] \t\t y2[kk] \t\t x3[kk] \t\t y3[kk] \t\t time[kk] \n");
  for (int kk = 0; kk < NUM_EV; kk++) 
  {
    fprintf(pf, " %lf \t %lf \t\t\t %lf \t %lf \t\t\t %lf \t %lf \t %lf \n", X1[kk][0], X1[kk][1], X2[kk][0], X2[kk][1], X3[kk][0], X3[kk][1], time[kk]) ; 
  }
fclose(pf);
}
