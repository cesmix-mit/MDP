#ifndef __CONFIGURATION_H__
#define __CONFIGURATION_H__

class CConfiguration {
private:
public:          
    configstruct config;
    appstruct app;
    commonstruct common;
        
    // default constructor 
    CConfiguration(){}; 
    
    // constructor 
    CConfiguration(string configfile, string appfile, Int mpiprocs, Int mpirank, Int backend); 
    
    // destructor        
    ~CConfiguration(); 
    
    void ReadConfigStruct(string filename, Int mpiprocs, Int mpirank, Int backend);
    void ReadAppStruct(string filename, Int mpiprocs, Int mpirank, Int backend);    
    void SetCommonStruct(Int mpiprocs, Int mpirank, Int backend);
    void Init(string configfile, string appfile, Int mpiprocs, Int mpirank, Int backend);     
};

#endif        

