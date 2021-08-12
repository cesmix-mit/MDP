/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __CALCULATION
#define __CALCULATION

#include "Calculation.h"

// constructor 
CCalculation::CCalculation(string filein, string fileout, Int mpiprocs, Int mpirank, Int backend) 
{
    if (backend == 2) {
        appstruct happ; 
        implReadInputFiles(happ, config, common, filein, fileout, mpiprocs, mpirank, backend);    
        implSetAppStruct(app, happ, backend);  
#ifdef HAVE_CUDA            
        // create cuda event handle
        CUDA_CHECK(cudaEventCreate(&common.eventHandle));

        // create cublas handle
        CHECK_CUBLAS(cublasCreate(&common.cublasHandle));
        CHECK_CUBLAS(cublasSetPointerMode(common.cublasHandle, CUBLAS_POINTER_MODE_HOST));                     
        //     CHECK_CUBLAS(cublasSetPointerMode(common.cublasHandle, CUBLAS_POINTER_MODE_DEVICE)); 
#endif        
    }
    else
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
    if (common.K > 0)
        InitSphericalHarmonics(sh, common);            
    implSetTempStruct(tmp, common);  
        
    int M = common.Ncoeff + common.Nempot;
    //printf("%i %i %i\n",M, common.Ncoeff, common.Nempot);
    TemplateMalloc(&sys.c, M, common.backend);  
    TemplateMalloc(&sys.d, M, common.backend);   
    TemplateMalloc(&sys.ee, common.Nempot*common.inummax, common.backend);  
    if (common.dftdata > 1)
        TemplateMalloc(&sys.dd, common.dim*common.inummax*M, common.backend);   
    else
        TemplateMalloc(&sys.dd, common.dim*common.inummax*common.Nempot, common.backend);       
    
    if (common.potential==0)
        common.M = common.Nempot;
    else if (common.potential==1)
        common.M = common.Ncoeff;
    else if (common.potential==2)
        common.M = M;
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
    // build neighbor list
    implNeighborList(nb, common, app, tmp, x, common.inum);    
    
    // copy x to xhold to check neighbor list
    ArrayCopy(sys.xhold, sys.x, common.dim*common.inum, common.backend);
    //Unmap(sys.xhold, sys.x, app.dom.h, sys.image, common.triclinic, common.dim, common.inum, common.backend);
}

void CCalculation::EmpiricalPotentialDescriptors(dstype *e, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    if (common.backend == 1)
        cpuEmpiricalPotentialDescriptors(sys.ee, nb, common, app, tmp, x, q, param, nparam);    
#ifdef USE_OMP
    if (common.backend == 4)
        ompEmpiricalPotentialDescriptors(sys.ee, nb, common, app, tmp, x, q, param, nparam);    
#endif                          
#ifdef USE_HIP
    if (common.backend == 3)
        hipEmpiricalPotentialDescriptors(sys.ee, nb, common, app, tmp, x, q, param, nparam);    
#endif                                  
#ifdef USE_CUDA            
    if (common.backend == 2)
        gpuEmpiricalPotentialDescriptors(sys.ee, nb, common, app, tmp, x, q, param, nparam);    
#endif                          
    
    dstype *onevec =  &tmp.tmpmem[0];  
    ArraySetValue(onevec, 1.0, common.inum, common.backend);
    PGEMTV(common.cublasHandle, common.inum, common.Nempot, &one, sys.ee, common.inum, onevec, inc1, &one, e, inc1, common.backend);                                    
}

void CCalculation::EmpiricalPotentialDescriptors(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    ArraySetValue(sys.ee, 0.0, common.inum*common.Nempot, common.backend);  
    //implEmpiricalPotentialDescriptors(sys.ee, f, nb, common, app, tmp, x, q, param, nparam);              
    if (common.backend == 1)
        cpuEmpiricalPotentialDescriptors(sys.ee, f, nb, common, app, tmp, x, q, param, nparam);              
#ifdef USE_OMP
    if (common.backend == 4)
        ompEmpiricalPotentialDescriptors(sys.ee, f, nb, common, app, tmp, x, q, param, nparam);              
#endif                          
#ifdef USE_HIP
    if (common.backend == 3)
        hipEmpiricalPotentialDescriptors(sys.ee, f, nb, common, app, tmp, x, q, param, nparam);              
#endif                                  
#ifdef USE_CUDA            
    if (common.backend == 2)
        gpuEmpiricalPotentialDescriptors(sys.ee, f, nb, common, app, tmp, x, q, param, nparam);              
#endif                          
        
    dstype *onevec =  &tmp.tmpmem[0];  
    ArraySetValue(onevec, 1.0, common.inum, common.backend);
    PGEMTV(common.cublasHandle, common.inum, common.Nempot, &one, sys.ee, common.inum, onevec, inc1, &one, e, inc1, common.backend);                                        
}

void CCalculation::EmpiricalPotentialEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int *nparam) 
{        
    //implEmpiricalPotentialEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);       
    if (common.backend == 1)
        cpuEmpiricalPotentialEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);       
#ifdef USE_OMP
    if (common.backend == 4)
        ompEmpiricalPotentialEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);       
#endif                          
#ifdef USE_HIP
    if (common.backend == 3)
        hipEmpiricalPotentialEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);       
#endif                                  
#ifdef USE_CUDA            
    if (common.backend == 2)
        gpuEmpiricalPotentialEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);       
#endif                                  
}

void CCalculation::EmpiricalPotentialEnergyForce(dstype *e, dstype *f, dstype* x, dstype *coeff, dstype *q, dstype *param, Int *nparam) 
{        
    Int dim = common.dim;
    Int inum = common.inum;
    Int Nempot = common.Nempot;    
    
    ArraySetValue(sys.ee, 0.0, inum*Nempot, common.backend);  
    ArraySetValue(sys.dd, 0.0, dim*inum*Nempot, common.backend);  
    
    //implEmpiricalPotentialDescriptors(sys.ee, sys.dd, nb, common, app, tmp, x, q, param, nparam);       
    if (common.backend == 1)
        cpuEmpiricalPotentialDescriptors(sys.ee, sys.dd, nb, common, app, tmp, x, q, param, nparam);       
#ifdef USE_OMP
    if (common.backend == 4)
        ompEmpiricalPotentialDescriptors(sys.ee, sys.dd, nb, common, app, tmp, x, q, param, nparam);       
#endif                          
#ifdef USE_HIP
    if (common.backend == 3)
        hipEmpiricalPotentialDescriptors(sys.ee, sys.dd, nb, common, app, tmp, x, q, param, nparam);       
#endif                                  
#ifdef USE_CUDA            
    if (common.backend == 2)
        gpuEmpiricalPotentialDescriptors(sys.ee, sys.dd, nb, common, app, tmp, x, q, param, nparam);       
#endif                                     
    PGEMNV(common.cublasHandle, inum, Nempot, &one, sys.ee, inum, coeff, inc1, &zero, e, inc1, common.backend);    
    PGEMNV(common.cublasHandle, dim*inum, Nempot, &minusone, sys.dd, dim*inum, coeff, inc1, &zero, f, inc1, common.backend);    
}

void CCalculation::RadialSphericalHarmonicDescriptors(dstype *e, dstype* x, dstype *q, dstype *param, Int nparam) 
{        
    if ((common.descriptor==0) && (common.K > 0))
       implSphericalHarmonicBesselDescriptors(e, nb, common, app, tmp, sh, x, q, param, nparam);              
}

void CCalculation::RadialSphericalHarmonicDescriptors(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{            
    if ((common.descriptor==0) && (common.K > 0))
       implSphericalHarmonicBesselDescriptors(e, f, nb, common, app, tmp, sh, x, q, param, nparam);                         
}

void CCalculation::RadialSphericalHarmonicEnergyForce(dstype *e, dstype *f, dstype* x, dstype *coeff, dstype *q, dstype *param, Int nparam) 
{    
    if ((common.descriptor==0) && (common.K > 0))
        implSphericalHarmonicBesselEnergyForce(e, f, nb, common, app, tmp, sh, x, coeff, q, param, nparam);              
}

void CCalculation::PotentialDescriptors(dstype *e, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    if (common.potential == 0)
        this->EmpiricalPotentialDescriptors(e, x, q, param, nparam);              
    else if (common.potential == 1)
        this->RadialSphericalHarmonicDescriptors(e, x, q, param, 0);            
    else if (common.potential == 2) {
        this->EmpiricalPotentialDescriptors(e, x, q, param, nparam);              
        this->RadialSphericalHarmonicDescriptors(&e[common.Nempot], x, q, param, 0);                    
    }    
    else
        error("Potential is not implemented");    
}

void CCalculation::PotentialDescriptors(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    if (common.potential == 0)
        this->EmpiricalPotentialDescriptors(e, f, x, q, param, nparam);        
    else if (common.potential == 1)
        this->RadialSphericalHarmonicDescriptors(e, f, x, q, param, 0);            
    else if (common.potential == 2) {
        this->EmpiricalPotentialDescriptors(e, f, x, q, param, nparam);              
        this->RadialSphericalHarmonicDescriptors(&e[common.Nempot], &f[common.Nempot*common.dim*common.inum], x, q, param, 0);                    
    }    
    else
        error("Potential is not implemented");            
}

void CCalculation::PotentialEnergyForce(dstype *e, dstype *f, dstype *x, dstype *coeff, dstype *q, dstype *param, Int *nparam) 
{     
    ArraySetValue(e, 0.0, common.inum, common.backend);  
    ArraySetValue(f, 0.0, common.inum*common.dim, common.backend);  

    if (common.potential == 0)
        this->EmpiricalPotentialEnergyForce(e, f, x, coeff, q, param, nparam);              
    else if (common.potential == 1)
        this->RadialSphericalHarmonicEnergyForce(e, f, x, coeff, q, param, 0);            
    else if (common.potential == 2) {
        this->EmpiricalPotentialEnergyForce(e, f, x, coeff, q, param, nparam);   
        this->RadialSphericalHarmonicEnergyForce(e, f, x, &coeff[common.Nempot], q, param, 0);  
    }    
    else
        error("Potential is not implemented");        
}


// void CCalculation::NonbondedSingleEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
// {    
//     implNonbondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
// }
// 
// void CCalculation::BondedSingleEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
// {    
//     implBondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
// }
// 
// void CCalculation::NonbondedPairEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
// {    
//     implNonbondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
// }
// 
// void CCalculation::BondedPairEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
// {    
//     implBondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
// }
// 
// void CCalculation::BondOrderPairEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
// {    
//     implBoPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
// }
// 
// void CCalculation::NonbondedTripletEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
// {    
//     implNonbondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
// }
// 
// void CCalculation::BondedTripletEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
// {    
//     implBondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
// }
// 
// void CCalculation::BondOrderTripletEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
// {    
//     implBoTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
// }
// 
// void CCalculation::NonbondedQuadrupletEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
// {    
//     implNonbondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
// }
// 
// void CCalculation::BondedQuadrupletEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
// {    
//     implBondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
// }

#endif        
