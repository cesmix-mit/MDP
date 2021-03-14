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
    
    //void Init(string filein, string fileout, Int mpiprocs, Int mpirank, Int backend) ;     
    void SetConfiguration(Int ci);    
    void GetPositions(dstype* x, Int ci);   
    void GetAtomtypes(Int* atomtype, Int ci);   
    void GetVelocities(dstype* v, Int ci);   
    void GetForces(dstype* f, Int ci);   
    void GetEnergy(dstype* e, Int ci);   
    void NeighborList(dstype* x);        
    void NonbondedPairEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, dstype *param, Int nparam); 

};


#endif        

