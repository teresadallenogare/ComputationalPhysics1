#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 200
#define M 30
#define P 150
#define T 2000
#define NUM_EN 20

double u_eff(double x, double b, double E);
double integrand( double x, double b, double E);
double radicand(double x, double b, double E);
double g_function( double x, double b, double E);
double derivative(double (*f)(double, double, double), double x0, double h, double b, double E);
double min_absc_g(double b, double E);
double min_absc_radicand(double b, double E);
double integral(double (*fun)( double, double, double), double R, int n, double b, double E );
//dx in integrale: E0_1b1, E05b1, E1b1: dx=1e-5;
//    E2b1: dx=0.5e-1;
//    E3b1: dx=1e-1;
//    E5b1: dx=1.53e-1
double theta_function(double b, double E);
double theta_4pi(double b, double E);
double find_b_zero(double E);
double find_b_orbit(double E);


double derLJ_potentialX(double x, double y);
double derLJ_potentialY(double x, double y);
double theta_Verlet( double b, double E);


double conta_part( double dth, double b, double E);

void fileData_b_E( double b[P], double theta_b[P][6]);
void fileData_bzero(double b[6], double E[6]);
void fileData2(double b[P], double theta[P][6]);
void fileDataCrossSection (double E[NUM_EN], double cs[NUM_EN]);

int main(int argc, const char * argv[]) {
   
    // Dati: energia ridotta e parametri d'impatto
    double E_star[]= {0.1, 0.5, 1, 2, 3, 5};
    double b_try[]= {0.2, 1, 3};
    int ind_b= 1;
    int ind_E= 2;
    
    double R= 50;
    
    double x[N], dx; // indica le posizioni
    double integranda[N];
    double u_efficace [N];
    double theta;
    double g[N];
    
// Considero dei valori di x e valuto il potenziale efficace e l'integranda per uno specifico valore di E
    x[0] = 0.80;   dx = 0.01; // scelti
//    Inizializzazione vettori
        u_efficace[0]= u_eff(x[0], b_try[ind_b], E_star[ind_E]);
        integranda[0]= integrand( x[0], b_try[ind_b], E_star[ind_E]);
        g[0]= g_function(x[0],b_try[ind_b],E_star[ind_E]);
        for (int kk=1 ; kk<N; kk++) {
                x[kk]= x[kk-1]+dx;
                u_efficace[kk]= u_eff(x[kk], b_try[ind_b], E_star[ind_E]);
                integranda[kk]= integrand( x[kk], b_try[ind_b], E_star[ind_E]);
               g[kk]= g_function(x[kk],b_try[ind_b],E_star[ind_E]);
    }
    
//    FILE *pf;
//
//    pf= fopen("dat_Eb.txt", "w");
//    fprintf(pf, "#Theta in funzione di b= %lf e di E= %lf: %lf\n",b_try[ind_b], E_star[ind_E], theta);
//    fprintf(pf,"# x[kk] \t\t u_eff[kk]  \t g[kk] \t integranda[kk]\n");
//    for (int kk=2; kk<N; kk++) {
//        fprintf (pf, "%lf\t   %lf\t  %lf\t %lf \n",x[kk], u_efficace[kk], g[kk], integranda[kk]);
//    }
//    fclose(pf);

//---------    METODO 1: QUADRATURA Studio theta in funzione di b per i valori di E usati sopra   -----------------------------
    double b[P], db;
    double theta_b[P][6];
  
    b[0] = 1e-6;   db = 0.025;
    for (int kk = 1; kk < P; kk++) {
        b[kk] = b[kk-1]+db;
    }
    for (int kk=0; kk < P; kk++) { // ciclo su b
        for (int jj = 0; jj < 6; jj++) { // ciclo su energia
            theta_b[kk][jj] = M_PI - 2*integral(&integrand, R, 1000000, b[kk], E_star[jj]);
            //theta_b[kk][jj] *= 180/M_PI;
            }
    }
    fileData_b_E(b, theta_b);

 // Find bzero
    double b_zero[6];
    for (int kk=0; kk<6; kk++) {
        b_zero[kk] = find_b_zero(E_star[kk]);
    }
    fileData_bzero(E_star, b_zero);
    

// -------------   METODO 2: RISOLUZIONE EQ. NEWTON CON ALGORITMO DI VERLET ------------------------------------------------
    double theta_2[P][6];

    for (int kk = 0; kk < P; kk++) { // ciclo su b (righe)
            for (int jj = 0; jj < 6; jj++) { // ciclo su energia (colonne)
                theta_2[kk][jj] = theta_Verlet(b[kk], E_star[jj]);    
                theta_2[kk][jj] *= 180/M_PI;    
        }
    }
    fileData2(b, theta_2);
    
    


    
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    printf("\nOrbiting: \n");
    double b_4pi;
    b_4pi= find_b_orbit(E_star[0]);
    printf("b for 4pi orbiting: %lf\n", b_4pi);
    
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    printf("\n2. CROSS SECTION:\n");
    
    double cross_section[NUM_EN];
    double E[NUM_EN], dE= 0.2;
    double thetabE[P][NUM_EN];
    double dtheta= M_PI/180;
    b[0] = 1e-6;   db = 1e-3;
    for (int kk=1; kk<P; kk++) {
        b[kk]=b[kk-1]+db;
    }
    E[0] = 0.1;
    for (int kk = 1; kk < NUM_EN; kk++) {
        E[kk] = E[kk-1]+dE;
    }


    for (int jj = 0; jj < NUM_EN; jj++) {// ciclo su energie
        for (int kk = 0; kk<P; kk++) {// ciclo su b
        thetabE[kk][jj] = M_PI - 2*integral(&integrand, R, 1000000, b[kk], E[jj]);
            cross_section[jj]+= 1/sin(thetabE[kk][jj]) *b[kk]*db/dtheta;
        }
    }
    
    printf("Cross section:\n");
    printf("E[kk]\t cross_section[kk]\n");
    for (int kk=0; kk<NUM_EN; kk++) {
        printf("%lf\t %lf\n", E[kk], cross_section[kk]);
    }
   fileDataCrossSection(E, cross_section);
    
   return 0;
}

double u_eff(double x, double b, double E){
    double v, u;
    v= 4*(pow(x, -12)-pow(x, -6)); // potenziale di L-J in unità ridotte
    u= v/E + pow(b,2)/pow(x, 2); // potenziale efficace
    return u;
}

double integrand(double x, double b, double E){
    double funz;
    double v, u;
    v= 4*(pow(x, -12)-pow(x, -6));
    u= v/E + pow(b,2)/pow(x, 2);
    funz= b/pow(x, 2) * 1/sqrt(1-u);
    return funz;
}

double radicand(double x, double b, double E){
    double rad;
    double v, u;
    v= 4*(pow(x, -12)-pow(x, -6));
    u= v/E + pow(b,2)/pow(x, 2);
    rad= 1-u;
    return rad;
}

double g_function( double x, double b, double E){
    double g;
    double v, u;
    v= 4*(pow(x, -12)-pow(x, -6));
    u= v/E + pow(b,2)/pow(x, 2);
    g= u-E;
    return g;
}

double derivative(double (*f)(double, double, double), double x0, double h, double b, double E){
    double df;
    df= (f(x0+h, b, E )-f(x0-h, b, E))/(2*h);
    return df;
}

// ricerca di x_min con metodo di Newton-Raphson
double min_absc_g( double b,  double E){
    double g[M], dg[M], x_min;
    double x[M];
    double h= 1e-5; // delta nella derivata (scelto io)
    x[0]= 0.5;// x0: stima iniziale dello zero della funzione;
    g[0]= g_function(x[0], b, E);
    dg[0]= derivative(&g_function, x[0], h, b, E);
//    dg[0]= (g_function(x[0]+h, b, E)- g_function(x[0]-h, b, E))/(2*h);
    for (int kk=1; kk<M; kk++) {
        x[kk]= x[kk-1]- g[kk-1]/dg[kk-1]; // punto successivo su cui valutare g
        g[kk]= g_function(x[kk],b,E);
//        dg[kk]= (g_function( x[kk]+h, b, E)- g_function( x[kk]-h, b, E))/(2*h);
        dg[kk]= derivative(&g_function, x[kk], h, b, E);
    }
    x_min= x[M-1];

    return x_min;
}

double min_absc_radicand_tdg(double b, double E){
    double rad[M], drad[M], x_min;
    double x[M];
    double h= 1e-5; // delta nella derivata (scelto io)
    x[0]= 0.5;// x0: stima iniziale dello zero della funzione;
    rad[0]= radicand(x[0], b, E);
    drad[0]= derivative(&radicand, x[0], h, b, E);
    for (int kk=1; kk<M; kk++) {
        x[kk]= x[kk-1]- rad[kk-1]/drad[kk-1]; // punto successivo su cui valutare g
        rad[kk]= radicand(x[kk],b,E);
//        dg[kk]= (g_function( x[kk]+h, b, E)- g_function( x[kk]-h, b, E))/(2*h);
        drad[kk]= derivative(&radicand, x[kk], h, b, E);
    }
    x_min= x[M-1];

    return x_min;
    
}

double min_absc_radicand(double b, double E){
  // start from a large value and go backwards (GG)
  double x_min = 10.0;
  double h = 1e-4;

  do {
    x_min -= h;
    
  } while(radicand(x_min, b, E) > 0.0);
  return(x_min + h); // so that the radicand is positive
    
}


//double min_absc_radicand(double b, double E){
//    double rad[M], drad[M], x_min;
//    double x[M];
//    double h= 1e-5; // delta nella derivata (scelto io)
//    x[0]= 0.5;// x0: stima iniziale dello zero della funzione;
//    rad[0]= radicand(x[0], b, E);
//    drad[0]= derivative(&radicand, x[0], h, b, E);
//    for (int kk=1; kk<M; kk++) {
//        x[kk]= x[kk-1]- rad[kk-1]/drad[kk-1]; // punto successivo su cui valutare g
//        rad[kk]= radicand(x[kk],b,E);
////        dg[kk]= (g_function( x[kk]+h, b, E)- g_function( x[kk]-h, b, E))/(2*h);
//        drad[kk]= derivative(&radicand, x[kk], h, b, E);
//    }
//    x_min= x[M-1];
//
//    return x_min;
//
//}

double integral(double (*fun)( double, double, double), double R,  int n, double b, double E){
    double f_min, f_R, x_kk, h;
    double trapez_sum= 0;
    double result;
    double  x_min;
    //  ricerca di x_min usando funzione min_absc (annullamento del radicando)
    //    altra opzione: min_absc_g (quando u_eff=E)
        x_min= min_absc_radicand( b, E);
        //x_min2= x_min+dx; // comincio a integrare dopo un poco altrimenti divergenza dell'integrando
        //  integrazione con metodo dei trapezi
        // x_min2 non serve visto che ho cambiato min_absc_radicand() - GG
        h= (R-x_min)/(n-1);
        f_min= (*fun)(x_min,b,E); // comincio a integrare da x_min+dx
        f_R= (*fun)(R,b,E);
        for (int kk=1; kk<(n-1); kk++) {
          x_kk= x_min+(kk-0.5)*h; // midpoint (GG)
            trapez_sum+= (*fun)(x_kk,b,E);
        }
        
        //result= h*(trapez_sum+1/2*(f_min+f_R));
        // il metodo del punto medio è forse meglio a causa della divergenza -GG
        // nota che più sopra (1/2) = 0 - GG
        result = h * trapez_sum;
        return result;
}


double theta_function(double b, double E){
    double th;
    double R= 50;
    th= M_PI- 2*integral(&integrand, R, 1000, b, E);
    return th;
}

double find_b_zero(double E){
    double th[M], dth[M], b_zero;
    double b[M];
    double h= 1e-5; // delta nella derivata (scelto io)
    b[0]= 1;// b0: stima iniziale dello zero della funzione;
    th[0]= theta_function(b[0], E);
    dth[0]=(theta_function(b[0]+h, E)- theta_function(b[0]-h, E))/(2*h);
    for (int kk=1; kk<M; kk++) {
        b[kk]= b[kk-1]- th[kk-1]/dth[kk-1]; // punto successivo su cui valutare g
        th[kk]= theta_function(b[kk], E);
        dth[kk]=(theta_function(b[kk]+h, E)- theta_function(b[kk]-h, E))/(2*h);
    }
    b_zero= b[M-1];
    return b_zero;
}

double derLJ_potentialX(double x, double y){
    return -24*x*(2*pow((pow(x,2)+pow(y,2)), -7) - pow((pow(x,2)+pow(y,2)), -4));

}
double derLJ_potentialY(double x, double y){
    return -24*y*(2*pow((pow(x,2)+pow(y,2)), -7) - pow((pow(x,2)+pow(y,2)), -4));

}

double theta_Verlet( double b, double E){
    double dt = 0.005;
    double x[T], y[T];
    double vx[T], vy[T];
    double ax[T], ay[T];
    double vx_provv = 0, vy_provv = 0;
    double x_provv, y_provv; 

    double dx = 0, dy = 0, dr = 0;
    double dphi = 0;
    double ipotenusa = 0;
    double theta = 0;
    double phi = 0;

    x[0] = -5;   y[0] = b;
    vx[0] = sqrt(2*E); vy[0] = 0;
    ax[0] = -derLJ_potentialX(x[0], y[0]);
    ay[0] = -derLJ_potentialY(x[0], y[0]);

    x_provv = x[0];  y_provv = b;
    
    for (int kk = 1; kk < T; kk++) {
        
        vx_provv = vx[kk-1] + ax[kk-1]*dt/2;
        vy_provv = vy[kk-1] + ay[kk-1]*dt/2;

        x[kk] = x[kk-1] + vx_provv*dt;
        y[kk] = y[kk-1] + vy_provv*dt;

        ax[kk]= -derLJ_potentialX(x[kk], y[kk]);
        ay[kk]= -derLJ_potentialY(x[kk], y[kk]);
        
        vx[kk]= vx_provv+ ax[kk]*dt/2;
        vy[kk]= vy_provv+ ay[kk]*dt/2;

        dx = x[kk] - x[kk-1]; 
        dy = y[kk] - y[kk-1]; 
        dr = sqrt( pow(dx,2) + pow(dy,2) );

        x_provv += dx;   y_provv += dy ;
        ipotenusa = sqrt(pow(x_provv,2) + pow(y_provv,2));
        dphi = dr /ipotenusa;
        
        phi += dphi;
    }
    // atan2
    theta = M_PI - phi;
    return theta;
}

double theta_4pi(double b, double E){
    double th;
    th= theta_function(b,E)+4*M_PI;
    return th;
}
double find_b_orbit( double E){
    double th[M], dth[M], b_4pi;
    double b[M];
    double h= 1e-5; // delta nella derivata (scelto io)
    b[0]= 1;// b0: stima iniziale dello zero della funzione;
    th[0]= theta_4pi(b[0], E);
    dth[0]= (theta_4pi(b[0]+h, E)- theta_4pi(b[0]-h, E))/(2*h);
    for (int kk=1; kk<M; kk++) {
        b[kk]= b[kk-1]- th[kk-1]/dth[kk-1]; // punto successivo su cui valutare g
        th[kk]= theta_4pi(b[kk], E);
        dth[kk]=(theta_4pi(b[kk]+h, E)- theta_4pi(b[kk]-h, E))/(2*h);
    }
    b_4pi= b[M-1];
    return b_4pi;
}

/* Insert data in files*/
void fileData_b_E ( double b[P], double theta_b[P][6]) {

  FILE *pf_b;
    pf_b= fopen("dati_theta_b_E.txt", "w");
    fprintf(pf_b, "\n#Deflection angle: 1. resolution in quadrature.\n");
    fprintf(pf_b, "#Anadamento di Theta in funzione di b. ");
    fprintf(pf_b, "# theta[b][E]: righe: b, colonne: E\n");
    fprintf(pf_b, "#b\t\t E=0.1\t\t E=0.5\t\t E=1\t\t E=2\t\t E=3\t\t E=5\n");
    for (int kk=0; kk < P; kk++) { // ciclo su b (righe)
      fprintf(pf_b, "%lf\t", b[kk]);
        for (int jj=0; jj < 6; jj++) { // ciclo su energia (colonne)
         fprintf(pf_b, " \t%lf",theta_b[kk][jj]);
        }
      fprintf(pf_b, "\n");
    }
    fclose(pf_b);
}
void fileData_bzero(double b[6], double E[6]){
  FILE *pf_b;
    pf_b= fopen("dati_theta_bzero.txt", "w");
    fprintf(pf_b, "\n#Deflection angle: 2. no deflection\n" );
    fprintf(pf_b,"\n# E[kk]\t b_zero[kk]:\n");
    for (int kk=0; kk<6; kk++) {
        fprintf(pf_b, "%lf\t %lf\n", E[kk], b[kk]);
    }
}
void fileData2(double b[P], double theta[P][6]){
FILE *pf_2;
    pf_2= fopen("dati_theta2.txt", "w");
    fprintf(pf_2, "\n#Deflection angle: 2. resolution with Verlet's algorithm.\n");
    fprintf(pf_2, "#Anadamento di Theta in funzione di b. ");
    fprintf(pf_2, "# theta2[b][E]\n");
    fprintf(pf_2, "#b\t\t E=0.1\t\t E=0.5\t\t E=1\t\t E=2\t\t E=3\t\t E=5\n");
    for (int kk=0; kk<P; kk++) { 
        fprintf(pf_2, "%lf\t", b[kk]);
            for (int jj=0; jj<6; jj++) { 
            fprintf(pf_2, " \t%lf",theta[kk][jj]);
        }
        fprintf(pf_2, "\n");
    }

    fclose(pf_2);
}

void fileDataCrossSection (double E[NUM_EN], double cs[NUM_EN]){
    FILE *pf_2;
    pf_2= fopen("datiCrossSection.txt", "w");
    fprintf(pf_2, "# E[kk]\t\t cs[kk]\n");
    for (int kk=0; kk<NUM_EN; kk++) { 
            fprintf(pf_2, " %lf \t%lf\n",E[kk], cs[kk]);
        }

    fclose(pf_2);

}
