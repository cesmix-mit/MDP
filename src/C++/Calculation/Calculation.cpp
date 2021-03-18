#ifndef __CALCULATION
#define __CALCULATION

#include "Calculation.h"

// constructor 
CCalculation::CCalculation(string filein, string fileout, Int mpiprocs, Int mpirank, Int backend) 
{
    implReadInputFiles(app, config, common, filein, fileout, mpiprocs, mpirank, backend);    
}

// destructor        
CCalculation::~CCalculation()
{
    tmp.freememory(common.backend);
    nb.freememory(common.backend);
    sys.freememory(common.backend);     
    app.freememory(common.backend);       
    sh.freememory(common.backend);       
    config.freememory(); // always in cpu memory    
    common.freememory(); // always in cpu memory
}

void CCalculation::SetConfiguration(Int ci)
{
    implSetConfiguration(nb, tmp, sys, app, config, common, ci);               
}

void CCalculation::GetPositions(dstype *x, Int ci)
{    
    implGetPositions(x, common, config, ci);
}

void CCalculation::GetAtomtypes(Int *atomtype, Int ci)
{    
    implGetAtomtypes(atomtype, common, config, ci);
}

void CCalculation::GetVelocities(dstype *v, Int ci)
{    
    implGetVelocities(v, common, config, ci);
}

void CCalculation::GetForces(dstype *f, Int ci)
{    
    implGetForces(f, common, config, ci);
}

void CCalculation::GetEnergy(dstype *e, Int ci)
{    
    implGetEnergy(e, common, config, ci);
}

void CCalculation::NeighborList(dstype* x)
{
    implNeighborList(nb, common, app, tmp, x, common.inum);    
}

void CCalculation::SetSphericalHarmonics()
{
    InitSphericalHarmonics(sh, common);            
}

void CCalculation::NonbondedSingleEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    implNonbondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

void CCalculation::BondedSingleEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    implBondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

void CCalculation::NonbondedPairEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    implNonbondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

void CCalculation::BondedPairEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    implBondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

void CCalculation::BondOrderPairEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    implBoPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

void CCalculation::NonbondedTripletEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    implNonbondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

void CCalculation::BondedTripletEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    implBondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

void CCalculation::BondOrderTripletEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    implBoTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

void CCalculation::NonbondedQuadrupletEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    implNonbondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

void CCalculation::BondedQuadrupletEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    implBondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

void CCalculation::EmpiricalPotentialEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    implEmpiricalPotentialEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

void CCalculation::RadialSphericalHarmonicDescriptors(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    if (common.descriptor==0)
        implSphericalHarmonicBesselDescriptors(e, f, nb, common, app, tmp, sh, x, q, param, nparam);              
}

void CCalculation::RadialSphericalHarmonicEnergyForce(dstype *e, dstype *f, dstype* x, dstype *coeff, dstype *q, dstype *param, Int nparam) 
{    
    if (common.descriptor==0)
        implSphericalHarmonicBesselEnergyForce(e, f, nb, common, app, tmp, sh, x, coeff, q, param, nparam);              
}


#endif        
