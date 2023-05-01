#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define dimB 10000
#define dimZ 1000
#define dimE 6
#define N 50000 /* number of Verlet evolutions */
#define POS_E 4

double b[dimB]; /* impact parameters */
double z[dimZ]; /* positions on the screen*/
double dz = 0.02;
double db = 0.001;  
double theta[2] = {0, M_PI/4.0}; 
double E[dimE] = {0.1, 0.3, 0.5, 1.0, 2.0, 5.0}; 

double X[N][2]; /* positions of the particle */
double dt = 0.01; /* step interval Verlet */ 
double ist[dimZ]; /* probability distribution for a specific value of E */ 

double derPot_X (double x, double z, double theta) {
    double r2 = x*x + z*z; 
    double a = x*sin(theta) - z*cos(theta); 
    return 1/(4*M_PI) * ( 3*x*pow(r2, -5.0/2.0)* a - sin(theta)* pow(r2, -3.0/2.0) ); 
}
double derPot_Z (double x, double z, double theta) {
    double r2 = x*x + z*z; 
    double a = x*sin(theta) - z*cos(theta); 
    return 1/(4*M_PI) * ( 3*z*pow(r2, -5.0/2.0)*a + cos(theta)* pow(r2, -3.0/2.0) ); 
}
void setToZeroV(unsigned int dim, double V[N]){
  for (int kk=0; kk<dim; kk++) {
            V[kk]=0;
        }
}
void setToZeroM(unsigned int dim, double M[N][2]){
    for (int jj=0; jj<dim; jj++) {
        for (int kk=0; kk<2; kk++) {
            M[jj][kk]=0;
        }
    
    }
} 
double evolVerlet (double E, double b, double theta) { 
    /* performs Verlet algorithm until x of the particle trajectory is 5 and returns the corrisponding z */
    setToZeroM (N,X); 
    X[0][0] = -5.0; 
    X[0][1] = b; 
    double V[2]; 
    V[0] = sqrt(2*E); 
    V[1] = 0.0; 
    double A[2]; 
    A[0] = -derPot_X(X[0][0], X[0][1], theta); 
    A[1] = -derPot_Z(X[0][0], X[0][1], theta); 
    int ii = 1; 
    while (ii<N && X[ii-1][0] < 5.0) {
            // v(t+dt/2) = v(t) + a(t)*(dt/2)
            V[0] += A[0]* (dt/2.0);     
            V[1] += A[1]* (dt/2.0);  
            // r(t+dt) = r(t) + v(t+dt/2)*dt
            X[ii][0] = X[ii-1][0] + V[0]*dt;  
            X[ii][1] = X[ii-1][1] + V[1]*dt;
            // a(t+dt)
            A[0] = -derPot_X(X[ii][0], X[ii][1], theta); 
            A[1] = -derPot_Z(X[ii][0], X[ii][1], theta);
            // v(t+dt) = v(t+dt/2) + a(t+dt)*(dt/2)
            V[0] += A[0]*(dt/2.0);
            V[1] += A[1]*(dt/2.0);  
            ii++;       
    }
    if (ii != N) {
        return X[ii-1][1]; 
    } else {
        return -100.0; 
    }
}
void counter (double z_pos, double ist[dimZ]) {
    /* counts the number of particles that go into a certain bin  */
    int trovato = 0;
    int s = 0;
while (s< dimZ && trovato == 0)
{
    if ( z_pos < -10.0+(s+1)*dz && z_pos >-10.0+ s*dz ){
         ist[s] = ist[s] + 1;
         trovato = 1;
    }
    s++;
}
   /*if (z_pos > -20.0) {
        int kk=0; 
        while ( z_pos > z[kk]) {
            kk++;   
        }
        ist[kk] = ist[kk] + 1;
    } */
}

void fileData0(double zz[dimZ], double ist[dimZ]) {
    FILE *fd;     
    fd = fopen( "dataIst_theta0.txt","w"); 
    fprintf(fd, "#z \t ist\n");
    for (int kk = 0; kk < dimZ; kk++) {
        fprintf (fd, "%lf \t %lf\n", zz[kk], ist[kk]/dimB); 

    }
}
void fileData1(double zz[dimZ], double ist[dimZ]) {
    FILE *fd;     
    fd = fopen( "dataIst_theta1.txt","w"); 
    fprintf(fd, "#z \t ist\n");
    for (int kk = 0; kk < dimZ; kk++) {
        fprintf (fd, "%lf \t %lf\n", zz[kk], ist[kk]/dimB); 

    }
}

int main(int argc, const char * argv[]) { 
   double z_pos; /* value of z in which the particle impacts */
    
    /* initialization of the z vector (from -10.0 to 10.0) */
    for (int ii=0; ii<dimZ; ii++) {
        z[ii] = -10.0 + ii*dz; 
    }
    /* initialization of the b vector ( from -5.0 to 5.0 )*/
    for (int ii=0; ii<dimB; ii++) {
        b[ii] = -5.0 + ii*db; 
    }

    /* histogram for theta = 0 */
    setToZeroV (dimZ, ist); 
    for (int jj = 0; jj<dimB; jj++) {
         z_pos = evolVerlet (E[POS_E], b[jj], theta[0]); 
        counter (z_pos, ist); 
    }
    
    fileData0(z, ist); 

    /* histogram for theta = pi/4 */
    setToZeroV (dimZ, ist); 
        for (int jj=0; jj<dimB; jj++) {
            z_pos = evolVerlet (E[POS_E], b[jj], theta[1]); 
            counter (z_pos, ist); 
        }

     fileData1(z, ist); 

}