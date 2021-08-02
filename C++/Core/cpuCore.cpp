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
#include "cpuComputes.cpp"
#include "cpuVelocityIntegration.cpp"
#include "cpuFixes.cpp"

#endif


