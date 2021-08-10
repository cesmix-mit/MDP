/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __CPUCORE
#define __CPUCORE

#include <math.h>
#include <algorithm>

// #define min(a,b) ((a) < (b) ? (a) : (b))
// #define max(a,b) ((a) > (b) ? (a) : (b))
using std::min;
using std::max;
#define BIG 1.0e30

template <typename T> void print2d(T* a, int m, int n)
{        
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++)
            printf("%g   ", a[j*m+i]);                
        printf("\n");
    }
    printf("\n");
}

template <typename T> void sub3(const T *v1, const T *v2, T *ans)
{
  ans[0] = v1[0] - v2[0];
  ans[1] = v1[1] - v2[1];
  ans[2] = v1[2] - v2[2];
}

/* ----------------------------------------------------------------------
   length of vector v
------------------------------------------------------------------------- */

template <typename T> T len3(const T *v)
{
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/* ----------------------------------------------------------------------
   squared length of vector v, or dot product of v with itself
------------------------------------------------------------------------- */

template <typename T> T lensq3(const T *v)
{
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

/* ----------------------------------------------------------------------
   ans = distance squared between pts v1 and v2
------------------------------------------------------------------------- */

template <typename T> T distsq3(const T *v1, const T *v2)
{
  T dx = v1[0] - v2[0];
  T dy = v1[1] - v2[1];
  T dz = v1[2] - v2[2];
  return dx * dx + dy * dy + dz * dz;
}

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

template <typename T> T dot3(const T *v1, const T *v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

template <typename T> void cross3(const T *v1, const T *v2, T *ans)
{
  ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
  ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
  ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

#include "cpuArrayOperations.cpp"
#include "cpuArrayPermute.cpp"
#include "cpuSort.cpp"
#include "cpuCoordinateTransformations.cpp"
#include "cpuNeighborList2D.cpp"
#include "cpuNeighborList3D.cpp"
#include "cpuNeighCount.cpp"
#include "cpuEnergyForceTally.cpp"
#include "cpuVirialTally.cpp"
#include "cpuClebschGordan.cpp"
#include "cpuSphericalHarmonicsBessel.cpp"
#include "cpuRadialSphericalHarmonics.cpp"
#include "cpuSnap.cpp"  
#include "cpuRandom.cpp"
#include "cpuDomain.cpp"
#include "cpuLattice.cpp"
#include "cpuAtom.cpp"
#include "cpuRegion.cpp"
#include "cpuComputes.cpp"
#include "cpuVelocityIntegration.cpp"
#include "cpuFixes.cpp"

#endif


