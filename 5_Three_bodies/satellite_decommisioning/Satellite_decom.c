#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#define STEPS 10600 /* number of Verlet_pro steps*/
#define M_Earth 5.972e24 //[kg]
#define R_Earth 6378e3 //[m]
#define G 6.6743e-11 // [N*m^2/ kg^2] 6.6743e-11 
#define Cd 8e-4 //[m^2/kg]
#define k1 1.2e4 //[m]
#define k2 2.2e4 //[m]


double X[2]; /* position of the satellite in cartesian coordinates */
double V[2]; /* velocity of the satellite in cartesian coordinates*/
double V_prov[2];
double Fg[2], Fdrag[2]; /* acceleration of the satellite in cartesian coordinates */
double F[2];
double dt = 1; /* Verlet_pro integration interval */

void initialConditions();
double module(double a, double b);
void setToZeroV(int dim, double V[dim]);
double rho( double X[2] );
void gravityForce( );
void dragForce();
void calculateForces();

void velocityVerlet_pro();

void fileData_r(double time[STEPS], double x[STEPS], double y[STEPS], double r[STEPS] );
void fileData_V(double time[STEPS], double vx[STEPS], double vy[STEPS], double v[STEPS], double a[STEPS]);
// --------------------------------------------------------------------------------------
int main(int argc, const char * argv[]){

double r[STEPS], v[STEPS], a[STEPS], time[STEPS];
double x[STEPS], y[STEPS], vx[STEPS], vy[STEPS];
double timeBreak;
int Break;
/* initial conditions */
initialConditions();
x[0] = X[0]; y[0] = X[1];
vx[0] = V[0]; vy[0] =V[1];
calculateForces();

r[0] = module(X[0],X[1]);
v[0] = module(V[0], V[1]);
a[0] = module(F[0], F[1]);
time[0] = 0;

/* evolution of the satellite   */
 for (int kk = 1; kk < STEPS && r[kk-1] > R_Earth ; kk++)
 {
   velocityVerlet_pro();


    x[kk] = X[0]; y[kk] = X[1];
    vx[kk] = V[0]; vy[kk] =V[1];
  
    r[kk] = module(X[0], X[1]);
    v[kk] = module(V[0], V[1]);
    a[kk] = module(F[0], F[1]);

    time[kk] = 30 * kk * dt; // dt = 1s  
    Break = kk;
}
timeBreak = time[Break];
printf("time of break: %lf s \n", timeBreak);


fileData_r(time, x, y, r);
fileData_V(time, vx, vy, v, a);
}
// --------------------------------------------------------------------------------------

double module(double a, double b){
    return sqrt(a*a + b*b);
}

void setToZeroV(int dim, double V[dim]){
    for (int kk = 0; kk < dim; kk++)
    {
        V[kk] = 0;
    }
    
}

void initialConditions(){
    /* initial positions and velocities of the satellite in cartesian coordinates */
   setToZeroV(2, X); setToZeroV(2, V);
   X[0] = R_Earth + 120e3; // [m] ;
   X[1] = 0;
   V[0] = 0; V[1] = sqrt( G * M_Earth/ X[0]);

}

 double rho( double X[2]){
     double h;
     h = module(X[0], X[1]) - R_Earth;
     return 1.225 * exp( - h/k1 - pow(h/k2, 3.0/2.0) ); //[kg/m^3]
}

void gravityForce( ){
    setToZeroV(2, Fg);
  Fg[0] = - G * M_Earth * X[0] / pow( module(X[0], X[1]), 3.0);
  Fg[1] = - G * M_Earth * X[1] / pow( module(X[0], X[1]), 3.0);
// printf("Fg[0]: %lf, Fg[1]: %lf\n ", Fg[0], Fg[1]);

}
void dragForce( ){
   double Rho;
   Rho = rho(X);
   Fdrag[0] = - Cd * Rho * V[0] * module(V[0], V[1]);
   Fdrag[1] = - Cd * Rho * V[1] * module(V[0], V[1]);
//printf("Fdrag[0]: %lf, Fdrag[1]: %lf\n\n ", Fdrag[0], Fdrag[1]);
}

void calculateForces(){
    /* calculates the total force acting on the body */
    gravityForce(); dragForce();
    F[0] = Fg[0] + Fdrag[0];
    F[1] = Fg[1] + Fdrag[1];
    
}


void velocityVerlet_pro(){
    /* performs 30 step of the algorithm VelocityVerlet with dependence on velocity */
for(int c = 0; c < 30; c++ ){
X[0] += V[0] * dt + F[0] * dt*dt/2; // x(n+1)
X[1] += V[1] * dt + F[1] * dt*dt/2;
V_prov[0] = V[0] + F[0] * dt/2; // v(n+1/2)
V_prov[1] = V[1] + F[1] * dt/2;

V[0] = V_prov[0];
V[1] = V_prov[1];

calculateForces();

V[0] += F[0] * dt/2; // v(n+1)
V[1] += F[1] * dt/2;

calculateForces();

V[0] = V_prov[0] + F[0] * dt/2;
V[1] = V_prov[1] + F[1] * dt/2;
}

}



void fileData_r(double time[STEPS], double x[STEPS], double y[STEPS], double r[STEPS] ){
/* insert data in FILE .txt*/
    FILE *pf;

    pf = fopen("dataR.txt", "w");
    fprintf(pf, "#time[kk]\t\t\t x[kk] \t\t y[kk] \t\t r[kk]\n");
    for (int kk = 0; kk < STEPS; kk++)
    {
        fprintf(pf, "%lf \t %lf\t %lf\t %lf \n", time[kk], x[kk], y[kk] , r[kk]);
    }
    fclose(pf);
}

void fileData_V(double time[STEPS], double vx[STEPS], double vy[STEPS], double v[STEPS], double a[STEPS] ){
/* insert data in FILE .txt*/
    FILE *pf;

    pf = fopen("dataV.txt", "w");
    fprintf(pf, "#time[kk]\t\t vx[kk] \t vy[kk] \t\t v[kk]\t\t a[kk]\n");
    for (int kk = 0; kk < STEPS; kk++)
    {
        fprintf(pf, "%lf \t %lf\t %lf\t %lf\t %lf \n", time[kk], vx[kk], vy[kk] , v[kk], a[kk]);
    }
    fclose(pf);
}

