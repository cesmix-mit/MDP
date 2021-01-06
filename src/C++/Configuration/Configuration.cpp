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
    tmp.freememory(common.backend);
    config.freememory();
    common.freememory();
}

void CConfiguration::ReadConfigStruct(string filename, Int backend)
{
    
}

void CConfiguration::ReadAppStruct(string filename, Int backend)
{
    
}

void CConfiguration::SetCommonStruct(Int backend)
{
    
}

void CConfiguration::SetTmpStruct(Int backend)
{
    
}

void CConfiguration::SetStructs(string configfile, string appfile, Int mpiprocs, Int mpirank, Int backend)
{
    this->ReadConfigStruct(configfile, backend);
    this->ReadAppStruct(appfile, backend);
    this->SetCommonStruct(backend);
    this->SetTmpStruct(backend);
}

// constructor 
CConfiguration::CConfiguration(string configfile, string appfile, Int mpiprocs, Int mpirank, Int backend) 
{
    this->SetStructs(configfile, appfile, mpiprocs, mpirank, backend);
}


#endif

