#ifndef __CPUCORE
#define __CPUCORE

#include <math.h>
#include <algorithm>

// #define min(a,b) ((a) < (b) ? (a) : (b))
// #define max(a,b) ((a) > (b) ? (a) : (b))
using std::min;
using std::max;

#include "cpuArrayOperations.cpp"
#include "cpuArrayPermute.cpp"
#include "cpuSort.cpp"
#include "cpuCoordinateTransformations.cpp"
#include "cpuNeighborList.cpp"
#include "cpuNeighborList2D.cpp"
#include "cpuNeighborList3D.cpp"
#include "cpuNeighCount.cpp"
#include "cpuForceCalculation.cpp"
#include "cpuSphericalHarmonicFunctions.cpp"
#include "cpuClebschGordan.cpp"
#include "cpuSphericalBessel.cpp"
#include "cpuSphericalHarmonics.cpp"
#include "cpuSphericalHarmonicsBessel.cpp"
#include "cpuRadialSphericalHarmonics.cpp"

#endif


