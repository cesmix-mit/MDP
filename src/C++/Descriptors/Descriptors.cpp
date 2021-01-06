#ifndef __DESCRIPTORS
#define __DESCRIPTORS

#include "SphericalHarmonics.cpp"
#include "Descriptors.h"

// destructor 
CDescriptors::~CDescriptors()
{    
     nb.freememory(backend);
     sys.freememory(backend);
}

void CDescriptors::NeighborList(CConfiguration& conf, Int confignum)
{
    
}

void CDescriptors::Energy(dstype* coeff)
{
}

void CDescriptors::Forces(dstype* coeff)
{
}

void CDescriptors::Stresses(dstype* coeff)
{
}
        
void CDescriptors::BasisFunctions()
{
}

void CDescriptors::BasisFunctionsDeriv()
{
}

#endif        

