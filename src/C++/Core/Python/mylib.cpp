#include <stdio.h>
#include "mylib.h"

// !clang++ -std=c++11 -Wall -Wextra -pedantic -c -fPIC mylib.cpp -o mylib.o
// !clang++ -shared mylib.o -o mylib.dylib
        
 void test_empty(void) {
      puts("Hello from C");
 }
  
float test_add(float x, float y) {
     return x + y;
}

double test_add_double(double x, double y) {
     return x + y;
}

void test_passing_array(int *data, int len) {
     printf("Data as received from Python\n");
     for(int i = 0; i < len; ++i) {
         printf("%d ", data[i]);
     }
     puts("");
 
     // Modifying the array
     for(int i = 0; i < len; ++i) {
         data[i] = -i;
     }
}

void test_add_array(double *x, double *y, double *z, int n) {
    for (int i=0; i<n; i++)
        x[i] = y[i] + z[i];
}


template <typename T> void test_axpby(T *z, T *x, T *y, T a, T b, int n) {    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        z[i] = a*x[i] + b*y[i];        
}
template void test_axpby(int*, int*, int*, int, int, int);
template void test_axpby(float*, float*, float*, float, float, int);
template void test_axpby(double*, double*, double*, double, double, int);
void test_axpby_int(int* z, int* x, int* y, int a, int b, int n) {
    test_axpby(z, x, y, a, b, n);
}
void test_axpby_float(float* z, float* x, float* y, float a, float b, int n) {
    test_axpby(z, x, y, a, b, n);
}
void test_axpby_double(double* z, double* x, double* y, double a, double b, int n) {
    test_axpby(z, x, y, a, b, n);
}





