#ifndef __CONFIGURATION
#define __CONFIGURATION

#include "Configuration.h"
#include "errormsg.cpp"
#include "ioutilities.cpp"
#include "readbinaryfiles.cpp"

#ifdef HAVE_CUDA
#include "gpuDeviceInfo.cpp"
#endif

// destructor        
CConfiguration::~CConfiguration()
{
    app.freememory(common.backend);
    config.freememory();    
    common.freememory();
}

void CConfiguration::ReadConfigStruct(string filename, Int mpiprocs, Int mpirank, Int backend)
{
    
}

void CConfiguration::ReadAppStruct(string filename, Int mpiprocs, Int mpirank, Int backend)
{
    
}

void CConfiguration::SetCommonStruct(Int mpiprocs, Int mpirank, Int backend)
{
    
}

void CConfiguration::Init(string configfile, string appfile, Int mpiprocs, Int mpirank, Int backend)
{
    this->ReadConfigStruct(configfile, mpiprocs, mpirank, backend);
    this->ReadAppStruct(appfile, mpiprocs, mpirank, backend);
    this->SetCommonStruct(mpiprocs,  mpirank, backend);
}

// constructor 
CConfiguration::CConfiguration(string configfile, string appfile, Int mpiprocs, Int mpirank, Int backend) 
{
    this->Init(configfile, appfile, mpiprocs, mpirank, backend);
}

#endif

