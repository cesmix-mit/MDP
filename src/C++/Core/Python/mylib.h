#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void test_empty(void);
float test_add(float x, float y);
double test_add_double(double x, double y);
void test_passing_array(int *data, int len);
void test_add_array(double *x, double *y, double *z, int n);

void test_axpby_int(int*, int*, int*, int, int, int);
void test_axpby_float(float*, float*, float*, float, float, int);
void test_axpby_double(double*, double*, double*, double, double, int);

#ifdef __cplusplus
}
#endif


