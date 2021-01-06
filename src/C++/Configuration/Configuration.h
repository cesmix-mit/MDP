#ifndef __CONFIGURATION_H__
#define __CONFIGURATION_H__

class CConfiguration {
private:
public:          
    commonstruct common;   
    configstruct config;
    appstruct app;
    tempstruct tmp;
    
    // default constructor 
    CConfiguration(){}; 
    
    // constructor 
    CConfiguration(string configfile, string appfile, Int mpiprocs, Int mpirank, Int backend); 
    
    // destructor        
    ~CConfiguration(); 
    
    void ReadConfigStruct(string filename, Int backend);
    void ReadAppStruct(string filename, Int backend);    
    void SetCommonStruct(Int backend);
    void SetTmpStruct(Int backend);         
    void SetStructs(string configfile, string appfile, Int mpiprocs, Int mpirank, Int backend);     
};

#endif        

