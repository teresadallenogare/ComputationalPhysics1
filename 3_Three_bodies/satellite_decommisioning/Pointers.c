#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


int main(int argc, const char * argv[]){
    
    double r_ist;
    double v_ist;
    double *pr_ist = &r_ist;
    double *pv_ist = &v_ist;
    //void initialConditions(double *r, double *v);
    // initialConditions(pr_ist, pv_ist);
    //  initialConditions(&r_ist, &v_ist);  printf("%lf, %lf \n", r_ist, v_ist);

   
    
    //void initialConditions(double *r, double *v){
        /* assegno il valore 120000 nella cella di indirizzo di una variabile pr */
      // *r = 120000.0; //[m]
      // *v = sqrt(M_Earth/ R);

   // }
int Array[10] = {3, 6, 9, 12, 15, 18, 21, 24, 27, 30};

//int *pLocation6 = &Array[6];
//int *pLocation0 = &Array[0];
//printf("pLocation0: %p\n", pLocation0);
//printf("pLocation6: %p\n", pLocation6);

int *pLocation0 = &Array[0];

for (int i = 0; i < 10; i++)
{
    printf("%p = %d\n", Array + i, *(Array + i));
}
// Array is the variable that holds the memory address
// *(Array) value stored at the memory location



}
