/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>

/* vector operations */
void quicksort(float *arr, int dummy, int elements);
int find_min_positive_num(int n1,int n2,float *a);
float global_percentiles(float *fullvec, int n1, int n2, float p1);
void V_write(char modfile[256],float *x_vector,int nn);
void V_read(char modfile[256],float *x_vector,int nn);
float V_max_abs(int n1, int n2, float * x);
float V_normL2(int n1, int n2, float * x);
/* operate vectors element-wise.  */
void V_init(int n1, int n2, float *x);
void V_copy(int n1, int n2, float *x, float *r);
void V_add(int n1, int n2, float *x, float *y,float *r);
void V_subtract(int n1, int n2, float *x, float *y,float *r);
void V_multiply(int n1, int n2, float *x, float *y,float *r);
void V_divide(int n1, int n2, float *x, float *y,float *r);
void V_square_sum(int n1, int n2, float *x, float *y,float *r);
void V_square_sum_a(int n1, int n2, float *x, float *y,float *r);
void V_AXPY(int n1, int n2, float a, float * x, float * y,float *r);
void V_AXPb(int n1, int n2, float a, float * x,  float b, float *r);
void V_scaling(int n1, int n2, float *x, float a);
void V_sqrt(int n1, int n2, float *x, float *r);
void V_multiply_sqrt(int n1, int n2, float *x, float *y, float *r);
void V_multiply_a(int n1, int n2,float *a, float *b,float *cor);
void V_multiply_a2(int n1, int n2, float *a, float *b, float *c, 
                   float *d, float *cor1, float *cor2);