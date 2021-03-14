#ifndef __CONFIGURATION
#define __CONFIGURATION

#include "Configuration.h"
#include "errormsg.cpp"
#include "ioutilities.cpp"
#include "Configuration_impl.cpp"
#include "Calculation_impl.cpp"

#ifdef HAVE_CUDA
#include "gpuDeviceInfo.cpp"
#endif

// constructor 
CConfiguration::CConfiguration(string filein, string fileout, Int mpiprocs, Int mpirank, Int backend) 
{
    implReadInputFiles(app, config, common, filein, fileout, mpiprocs, mpirank, backend);    
}

// destructor        
CConfiguration::~CConfiguration()
{
    tmp.freememory(common.backend);
    nb.freememory(common.backend);
    sys.freememory(common.backend);     
    app.freememory(common.backend);       
    config.freememory(); // always in cpu memory    
    common.freememory(); // always in cpu memory
}

void CConfiguration::SetConfiguration(Int ci)
{
    implSetConfiguration(nb, tmp, sys, app, config, common, ci);               
}

void CConfiguration::GetPositions(dstype *x, Int ci)
{    
    implGetPositions(x, common, config, ci);
}

void CConfiguration::GetAtomtypes(Int *atomtype, Int ci)
{    
    implGetAtomtypes(atomtype, common, config, ci);
}

void CConfiguration::GetVelocities(dstype *v, Int ci)
{    
    implGetVelocities(v, common, config, ci);
}

void CConfiguration::GetForces(dstype *f, Int ci)
{    
    implGetForces(f, common, config, ci);
}

void CConfiguration::GetEnergy(dstype *e, Int ci)
{    
    implGetEnergy(e, common, config, ci);
}

void CConfiguration::NeighborList(dstype* x)
{
    implNeighborList(nb, common, app, tmp, x, common.inum);    
}

void CConfiguration::NonbondedPairEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam) 
{    
    implNonbondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

// void implNonbondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
//         dstype* x, dstype *q, dstype* param, Int nparam)



// void CConfiguration::Init(string filein, string fileout, Int mpiprocs, Int mpirank, Int backend)
// {    
//     if (backend==1) {
//         implReadAppStruct(app, filein, mpiprocs, mpirank, backend);
//         implReadConfigStruct(config, filein, mpiprocs, mpirank, backend);
//         implSetCommonStruct(common, app, config, filein, fileout, mpiprocs,  mpirank, backend);                
//         this->SetNeighborStruct(0);
//         this->SetTempStruct();
//         this->SetSysStruct();            
//     }
//     else if (backend==2) {
//         appstruct happ;
//         implReadAppStruct(happ, filein, mpiprocs, mpirank, backend);
//         implReadConfigStruct(config, filein, mpiprocs, mpirank, backend);
//         implSetCommonStruct(common, happ, config, filein, fileout, mpiprocs,  mpirank, backend);
//         implSetAppStruct(app, happ, backend);
//         this->SetNeighborStruct(0);
//         this->SetTempStruct();
//         this->SetSysStruct();            
//     }    
// }

#endif

