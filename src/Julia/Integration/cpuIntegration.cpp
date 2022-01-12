/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

// !clang++ -std=c++11 -Wall -Wextra -pedantic -c -fPIC cpuIntegration.cpp -o cpuIntegration.o
// !clang++ -shared cpuIntegration.o -o cpuIntegration.dylib

#ifndef CPUINTEGRATION
#define CPUINTEGRATION

#include "cpuIntegration.h"

#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <cmath>

#define BIG 1.0e30
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
using std::min;
using std::max;

double cpuArraySum(double *a, int n)
{
    double b = a[0];
    for (int i=1; i<n; i++)
            b = b + a[i];    
    return b;
}

void cpuArraySetValue(double *y, double a, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = a;        
}

void cpuArrayMultiplyScalar(double *y, double a, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = a*y[i];        
}

void sub3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[0] - v2[0];
  ans[1] = v1[1] - v2[1];
  ans[2] = v1[2] - v2[2];
}

/* ----------------------------------------------------------------------
   length of vector v
------------------------------------------------------------------------- */

double len3(const double *v)
{
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/* ----------------------------------------------------------------------
   squared length of vector v, or dot product of v with itself
------------------------------------------------------------------------- */

double lensq3(const double *v)
{
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

/* ----------------------------------------------------------------------
   ans = distance squared between pts v1 and v2
------------------------------------------------------------------------- */

double distsq3(const double *v1, const double *v2)
{
  double dx = v1[0] - v2[0];
  double dy = v1[1] - v2[1];
  double dz = v1[2] - v2[2];
  return dx * dx + dy * dy + dz * dz;
}

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

double dot3(const double *v1, const double *v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

void cross3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
  ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
  ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void cpuArraySumEveryColumn(double *b, double *a, int m, int n)
{
    // m: number of columns
    for (int j=0; j<m; j++) {
        int k = n*j;
        b[j] = a[0 + k];
        for (int i=1; i<n; i++)
            b[j] += a[i + k];    
    }
}

void cpuArrayPlusAtColumnIndex(double *A, double *B, int *colind, int m, int n)
{  
    for (int ii=0; ii<n; ii++) { // loop over each column
        int i = colind[ii];   
        for (int j=0; j<m; j++) // loop over each row
            A[ii*m+j] += B[i*m+j]; 
    }
}

   //cpuArrayDistSquareSum
void cpuArrayDistSquareSum(double*y, double*x1, double*x2, int m, int n)
{    
    //#pragma omp parallel 
    for (int j=0; j<n; j++)
    {        
        y[j] = 0;        
        for (int i=0; i<m; i++)
            y[j] += (x2[i + m*j]-x1[i + m*j])*(x2[i + m*j]-x1[i + m*j]);
    }        
}

//#include "cpuCoordinateTransformations.cpp"
#include "cpuRandom.cpp"
#include "cpuDomain.cpp"
#include "cpuLattice.cpp"
#include "cpuAtom.cpp"
#include "cpuRegion.cpp"
#include "cpuComputes.cpp"
#include "cpuVelocityIntegration.cpp"
#include "cpuFixes.cpp"

#endif

