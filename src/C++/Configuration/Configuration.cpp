#ifndef __CONFIGURATION
#define __CONFIGURATION

#include "Configuration.h"
#include "errormsg.cpp"
#include "ioutilities.cpp"
#include "Configuration_impl.cpp"

#ifdef HAVE_CUDA
#include "gpuDeviceInfo.cpp"
#endif

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

void CConfiguration::SetNeighborStruct(Int ci)
{           
    nb.freememory(common.backend);
    implSetNeighborStruct(nb, common, config, ci);    
}

void CConfiguration::SetTempStruct()
{
    
}

void CConfiguration::SetSysStruct()
{
    
}

void CConfiguration::NeighborList(dstype* x, Int *atomtype, Int inum)
{
    implNeighborList(nb, common, app, tmp, x, atomtype, inum);    
}

void CConfiguration::Init(string filein, string fileout, Int mpiprocs, Int mpirank, Int backend)
{    
    if (backend==1) {
        implReadAppStruct(app, filein, mpiprocs, mpirank, backend);
        implReadConfigStruct(config, filein, mpiprocs, mpirank, backend);
        implSetCommonStruct(common, app, config, filein, fileout, mpiprocs,  mpirank, backend);                
        this->SetNeighborStruct(0);
        this->SetTempStruct();
        this->SetSysStruct();            
    }
    else if (backend==2) {
        appstruct happ;
        implReadAppStruct(happ, filein, mpiprocs, mpirank, backend);
        implReadConfigStruct(config, filein, mpiprocs, mpirank, backend);
        implSetCommonStruct(common, happ, config, filein, fileout, mpiprocs,  mpirank, backend);
        implSetAppStruct(app, happ, backend);
        this->SetNeighborStruct(0);
        this->SetTempStruct();
        this->SetSysStruct();            
    }    
}

// constructor 
CConfiguration::CConfiguration(string filein, string fileout, Int mpiprocs, Int mpirank, Int backend) 
{
    this->Init(filein, fileout, mpiprocs, mpirank, backend);
}


#endif

