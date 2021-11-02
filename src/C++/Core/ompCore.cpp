/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_OMPCORE
#define MDP_OMPCORE

//#include <omp.h>
#include <math.h>
#include <algorithm>

// #define min(a,b) ((a) < (b) ? (a) : (b))
// #define max(a,b) ((a) > (b) ? (a) : (b))
using std::min;
using std::max;

#include "ompArrayOperations.cpp"
#include "ompArrayPermute.cpp"
#include "ompSort.cpp"
#include "ompCoordinateTransformations.cpp"
#include "ompNeighborList2D.cpp"
#include "ompNeighborList3D.cpp"
#include "ompNeighCount.cpp"
#include "ompForceDecomposition.cpp"
#include "ompClebschGordan.cpp"
#include "ompSphericalHarmonicsBessel.cpp"
#include "ompRadialSphericalHarmonics.cpp"

#endif


