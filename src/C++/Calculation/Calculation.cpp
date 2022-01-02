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
    sna.freememory(common.backend);       
    //config.freememory(); // always in cpu memory    
    //common.freememory(); // always in cpu memory
}

void CCalculation::SetConfiguration(Int ci)
{
    implSetConfiguration(nb, tmp, sys, app, config, common, ci);       
    
    if (common.descriptor == 0) {// spherical harmonic
        InitSphericalHarmonics(sh, common);      
    }
    else if (common.descriptor == 1) {// snap
        InitSnap(sna, common);     
        common.Ncoeff = sna.ntypes*sna.ncoeff;
    }
    else 
        common.Ncoeff = 0;
    
    implSetTempStruct(tmp, common, sna, sh);              
    
    int M = common.Ncoeff + common.Nempot;
    common.M = M;    
    
    //printf("%i %i %i\n",M, common.Ncoeff, common.Nempot);
    if (common.training>0) {
        TemplateMalloc(&sys.c, M, common.backend);  
        TemplateMalloc(&sys.d, M, common.backend);   
        TemplateMalloc(&sys.ee, common.Nempot*common.inummax, common.backend);  
        TemplateMalloc(&sys.vv, 6*common.Nempot*common.inummax, common.backend);      
//     if (common.dftdata > 1)
//         TemplateMalloc(&sys.dd, common.dim*common.inummax*M, common.backend);   
//     else
//         TemplateMalloc(&sys.dd, common.dim*common.inummax*common.Nempot, common.backend);               
    }
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

void CCalculation::EmpiricalPotentialDescriptors(dstype *e, dstype *f, dstype *v, dstype* x, dstype *q, dstype *param, Int *nparam) 
{     
    ArraySetValue(sys.ee, 0.0, common.inum*common.Nempot, common.backend);   
    ArraySetValue(f, 0.0, common.dim*common.inum*common.Nempot, common.backend);   
    ArraySetValue(sys.vv, 0.0, 6*common.inum*common.Nempot, common.backend);  
    if (common.backend == 1)
        cpuEmpiricalPotentialDescriptors(sys.ee, f, sys.vv, nb, common, app, tmp, x, q, param, nparam);              
#ifdef USE_OMP
    if (common.backend == 4)
        ompEmpiricalPotentialDescriptors(sys.ee, f, sys.vv, nb, common, app, tmp, x, q, param, nparam);              
#endif                          
#ifdef USE_HIP
    if (common.backend == 3)
        hipEmpiricalPotentialDescriptors(sys.ee, f, sys.vv, nb, common, app, tmp, x, q, param, nparam);              
#endif                                  
#ifdef USE_CUDA            
    if (common.backend == 2)
        gpuEmpiricalPotentialDescriptors(sys.ee, f, sys.vv, nb, common, app, tmp, x, q, param, nparam);              
#endif                          
            
    dstype *onevec =  &tmp.tmpmem[0];  
    ArraySetValue(onevec, 1.0, 6*common.inum, common.backend);
    PGEMTV(common.cublasHandle, common.inum, common.Nempot, &one, sys.ee, common.inum, onevec, inc1, &one, e, inc1, common.backend);                                        
    PGEMTV(common.cublasHandle, common.inum, 6*common.Nempot, &one, sys.vv, common.inum, onevec, inc1, &one, v, inc1, common.backend);                                        
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

void CCalculation::RadialSphericalHarmonicEnergyForceVirial(dstype *e, dstype *f, dstype *v, dstype* x, dstype *coeff, dstype *q, dstype *param, Int nparam) 
{    
    if ((common.descriptor==0) && (common.K > 0))
        implSphericalHarmonicBesselEnergyForceVirial(e, f, v, nb, common, app, tmp, sh, x, coeff, q, param, nparam);              
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

void CCalculation::PotentialEnergyForce(dstype *e, dstype *f, dstype *x, dstype *param, Int *nparam) 
{     
    ArraySetValue(e, 0.0, common.inum, common.backend);  
    ArraySetValue(f, 0.0, common.inum*common.dim, common.backend);  

    if (common.potential == 0)
        this->EmpiricalPotentialEnergyForce(e, f, x, sys.q, param, nparam);              
//     else if (common.potential == 1)
//         this->RadialSphericalHarmonicEnergyForce(e, f, x, coeff, sys.q, param, 0);            
//     else if (common.potential == 2) {
//         this->EmpiricalPotentialEnergyForce(e, f, x, coeff, sys.q, param, nparam);   
//         this->RadialSphericalHarmonicEnergyForce(e, f, x, &coeff[common.Nempot], sys.q, param, 0);  
//     }    
//     else
//         error("Potential is not implemented");        
}

void CCalculation::EmpiricalPotentialEnergyForceVirial(dstype *e, dstype *f, dstype *v, dstype* x, dstype *q, dstype *param, Int *nparam) 
{        
    //implEmpiricalPotentialEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);       
    if (common.backend == 1)
        cpuEmpiricalPotentialEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, nparam);       
#ifdef USE_OMP
    if (common.backend == 4)
        ompEmpiricalPotentialEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, nparam);       
#endif                          
#ifdef USE_HIP
    if (common.backend == 3)
        hipEmpiricalPotentialEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, nparam);       
#endif                                  
#ifdef USE_CUDA            
    if (common.backend == 2)
        gpuEmpiricalPotentialEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, nparam);       
#endif                
}

void CCalculation::EmpiricalPotentialEnergyForceVirial(dstype *e, dstype *f, dstype *v, dstype* x, 
        dstype *coeff, dstype *q, dstype *param, Int *nparam) 
{        
    Int dim = common.dim;
    Int inum = common.inum;
    Int Nempot = common.Nempot;    
    
    ArraySetValue(sys.ee, 0.0, inum*Nempot, common.backend);  
    ArraySetValue(sys.dd, 0.0, dim*inum*Nempot, common.backend);  
    ArraySetValue(sys.vv, 0.0, 6*inum*Nempot, common.backend);  
    
    //implEmpiricalPotentialDescriptors(sys.ee, sys.dd, nb, common, app, tmp, x, q, param, nparam);       
    if (common.backend == 1)
        cpuEmpiricalPotentialDescriptors(sys.ee, sys.dd, sys.vv, nb, common, app, tmp, x, q, param, nparam);       
#ifdef USE_OMP
    if (common.backend == 4)
        ompEmpiricalPotentialDescriptors(sys.ee, sys.dd, sys.vv, nb, common, app, tmp, x, q, param, nparam);       
#endif                          
#ifdef USE_HIP
    if (common.backend == 3)
        hipEmpiricalPotentialDescriptors(sys.ee, sys.dd, sys.vv, nb, common, app, tmp, x, q, param, nparam);       
#endif                                  
#ifdef USE_CUDA            
    if (common.backend == 2)
        gpuEmpiricalPotentialDescriptors(sys.ee, sys.dd, sys.vv, nb, common, app, tmp, x, q, param, nparam);       
#endif                                     
    PGEMNV(common.cublasHandle, inum, Nempot, &one, sys.ee, inum, coeff, inc1, &zero, e, inc1, common.backend);    
    PGEMNV(common.cublasHandle, dim*inum, Nempot, &minusone, sys.dd, dim*inum, coeff, inc1, &zero, f, inc1, common.backend);    
    PGEMNV(common.cublasHandle, 6*inum, Nempot, &minusone, sys.vv, 6*inum, coeff, inc1, &zero, f, inc1, common.backend);    
}

void CCalculation::PotentialEnergyForceVirial(dstype *e, dstype *f, dstype *v, dstype *x, dstype *coeff, dstype *q, dstype *param, Int *nparam) 
{     
    ArraySetValue(e, 0.0, common.inum, common.backend);  
    ArraySetValue(f, 0.0, common.inum*common.dim, common.backend);  
    ArraySetValue(v, 0.0, common.inum*6, common.backend);  
    
    this->EmpiricalPotentialEnergyForceVirial(e, f, v, x, coeff, sys.q, param, nparam);              
    
    if (common.descriptor == 1) // snap
        ComputePairSnap(e, f, v, x, sna, common, sys, nb, tmp);  
    else if (common.descriptor == 0) // spherical harmonic
        this->RadialSphericalHarmonicEnergyForceVirial(e, f, v, x,  &coeff[common.Nempot], sys.q, param, 0);            

//     if (common.potential == 0) {
//         this->EmpiricalPotentialEnergyForceVirial(e, f, v, x, coeff, q, param, nparam);              
// //     else if (common.potential == 1)
// //         this->RadialSphericalHarmonicEnergyForce(e, f, v, x, coeff, q, param, 0);            
// //     else if (common.potential == 2) {
// //         this->EmpiricalPotentialEnergyForce(e, f, v, x, coeff, q, param, nparam);   
// //         this->RadialSphericalHarmonicEnergyForce(e, f, v, x, &coeff[common.Nempot], q, param, 0);  
//     }    
//     else
//         error("Potential is not implemented");        
}

void CCalculation::PotentialEnergyForceVirial(dstype *e, dstype *f, dstype *v, dstype *x, dstype *param, Int *nparam) 
{     
    ArraySetValue(e, 0.0, common.inum, common.backend);  
    ArraySetValue(f, 0.0, common.inum*common.dim, common.backend);  
    ArraySetValue(v, 0.0, common.inum*6, common.backend);  
     
//     if (common.potential == 0) 
//         this->EmpiricalPotentialEnergyForceVirial(e, f, v, x, sys.q, param, nparam);              
//     else if (common.potential == 1)
//         this->RadialSphericalHarmonicEnergyForce(e, f, v, x, coeff, q, param, 0);            
//     else if (common.potential == 2) {
//         this->EmpiricalPotentialEnergyForce(e, f, v, x, coeff, q, param, nparam);   
//         this->RadialSphericalHarmonicEnergyForce(e, f, v, x, &coeff[common.Nempot], q, param, 0);  
//     }   
    
    INIT_TIMING;
    
    START_TIMING;
    this->EmpiricalPotentialEnergyForceVirial(e, f, v, x, sys.q, param, nparam);              
    END_TIMING(10);
    
    int dim = common.dim;
    int inum = common.inum;    
    int backend = common.backend;
    dstype *tmpmem = tmp.tmpmem;       
    ArraySumEveryColumn(tmpmem, e, 1, inum, backend);  
    ArraySumEveryColumn(&tmpmem[1], sys.vatom, 6, inum, backend);    
    printArray2D(tmpmem, 1, 7, backend);
    error("here")
    
    START_TIMING;
    if (common.descriptor == 1) { // snap
        if (common.backend<=1)
            ComputePairSnap3(e, f, v, x, sna, common, sys, nb, tmp);  
        else
            ComputePairSnap3(e, f, v, x, sna, common, sys, nb, tmp);                      
    }
    else if (common.descriptor == 0) // spherical harmonic
        this->RadialSphericalHarmonicEnergyForceVirial(e, f, v, x, sys.c, sys.q, param, 0);                               
    END_TIMING(11); 
}

void CCalculation::ThermoOutput(int flag)
{
    int dim = common.dim;
    int inum = common.inum;    
    int backend = common.backend;
    int *atomtype = nb.atomtype;            
    int *ilist = nb.alist;    
    
    dstype *tmpmem = tmp.tmpmem;    
    dstype *mass = app.atommass;    
    dstype *scalars = common.scalars;
    dstype mvv2e = common.mvv2e;            
    dstype boltz = common.boltz;            
    dstype nktv2p = common.nktv2p;    
    dstype tdof = (inum-1)*dim;    
    dstype tfactor = mvv2e/(tdof * boltz);
    dstype inv_volume = 1.0 / (common.dom.h[0] * common.dom.h[1] * common.dom.h[2]);  
    
    int nout;
    if (flag == 0) {
        nout=8;
        ArraySumEveryColumn(tmpmem, sys.e, 1, inum, backend);    
        ArraySumEveryColumn(&tmpmem[1], sys.vatom, 6, inum, backend);    
        ComputeKEAtom(&tmpmem[inum+nout], mass, sys.v, 2.0, atomtype, ilist, dim, inum, backend);    
        ArraySumEveryColumn(&tmpmem[7], &tmpmem[inum+nout], 1, inum, backend);
        TemplateCopytoHost(scalars, tmpmem, nout, backend);        
        common.pe = scalars[0]/inum;                
        common.temp = tfactor*scalars[7];    
        common.ke = 0.5 * tdof * boltz * common.temp/inum;                          
        common.pres = (tdof * boltz * common.temp + scalars[1] + scalars[2] + scalars[3]) / 3.0 * inv_volume * nktv2p;
    }
}

void CCalculation::ReferencePotentialDescriptors(dstype *ev, dstype *e, dstype *f, dstype *v, 
        dstype *bi, dstype *bd, dstype *bv, dstype *x, dstype *param, Int *nparam) 
{         
    ArraySetValue(e, 0.0, common.inum, common.backend);  
    ArraySetValue(f, 0.0, common.inum*common.dim, common.backend);  
    ArraySetValue(v, 0.0, common.inum*6, common.backend);  
             
    this->EmpiricalPotentialEnergyForceVirial(e, f, v, x, sys.q, param, nparam);                         
    ArraySumEveryColumn(ev, e, 1, common.inum, common.backend);    
    ArraySumEveryColumn(&ev[1], v, 6, common.inum, common.backend);              
    
    if (common.descriptor == 1) // snap
        ComputeSnap(bi, bd, bv, x, sna, common, sys, nb, tmp);  
    //else if (common.descriptor == 0) // spherical harmonic
        //this->RadialSphericalHarmonicEnergyForceVirial(e, f, v, x, sys.c, sys.q, param, 0);                    
}

void CCalculation::PotentialDescriptors(dstype *e, dstype *f, dstype *v, 
        dstype *bi, dstype *bd, dstype *bv, dstype *x, dstype *param, Int *nparam) 
{         
    if (common.Nempot>0) {
        ArraySetValue(e, 0.0, common.Nempot, common.backend);  
        ArraySetValue(f, 0.0, common.inum*common.dim*common.Nempot, common.backend);  
        ArraySetValue(v, 0.0, 6*common.Nempot, common.backend);               
        this->EmpiricalPotentialDescriptors(e, f, v, x, sys.q, param, nparam);                         
    }
    
    //void CCalculation::EmpiricalPotentialDescriptors(dstype *e, dstype *f, dstype *v, dstype* x, dstype *q, dstype *param, Int *nparam) 
    
    if (common.descriptor == 1) // snap
        ComputeSnap(bi, bd, bv, x, sna, common, sys, nb, tmp);  
    //else if (common.descriptor == 0) // spherical harmonic
        //this->RadialSphericalHarmonicEnergyForceVirial(e, f, v, x, sys.c, sys.q, param, 0);                    
}

#endif        
