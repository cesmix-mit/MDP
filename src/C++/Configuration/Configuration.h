#ifndef __CONFIGURATION_H__
#define __CONFIGURATION_H__

class CConfiguration {
private:
public:          
    appstruct app;
    configstruct config;    
    commonstruct common;
    neighborstruct nb; 
    sysstruct sys;     
    tempstruct tmp;    
    
    // default constructor 
    CConfiguration(){}; 
    
    // constructor 
    CConfiguration(string filein, string fileout, Int mpiprocs, Int mpirank, Int backend); 
    
    // destructor        
    ~CConfiguration(); 
    
    void Init(string filein, string fileout, Int mpiprocs, Int mpirank, Int backend) ;     
    void SetNeighborStruct(Int ci);
    void SetSysStruct();
    void SetTempStruct();     
    void NeighborList(dstype* x, Int *atomtype, Int inum);        
};

#endif        

