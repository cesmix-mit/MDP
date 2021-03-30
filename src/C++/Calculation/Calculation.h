#ifndef __CALCULATION_H__
#define __CALCULATION_H__

class CCalculation {
private:
public:          
    appstruct app;
    configstruct config;    
    commonstruct common;
    neighborstruct nb; 
    sysstruct sys;     
    tempstruct tmp;    
    shstruct sh;
    
    // default constructor 
    CCalculation(){}; 
    
    // constructor 
    CCalculation(string filein, string fileout, Int mpiprocs, Int mpirank, Int backend); 
    
    // destructor        
    ~CCalculation(); 
    
    // Read input files
    void SetConfiguration(Int ci);    
    void GetPositions(dstype* x, Int ci);   
    void GetAtomtypes(Int* atomtype, Int ci);   
    void GetVelocities(dstype* v, Int ci);   
    void GetForces(dstype* f, Int ci);   
    void GetEnergy(dstype* e, Int ci);   
    
    // Build neighbor list
    void NeighborList(dstype* x);        
    
    // Empirical potentials
//     void NonbondedSingleEnergyForce(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int nparam); 
//     void BondedSingleEnergyForce(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int nparam); 
     void NonbondedPairEnergyForce(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int nparam); 
//     void BondedPairEnergyForce(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int nparam); 
//     void BondOrderPairEnergyForce(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int nparam); 
//     void NonbondedTripletEnergyForce(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int nparam); 
//     void BondedTripletEnergyForce(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int nparam); 
//     void BondOrderTripletEnergyForce(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int nparam); 
//     void NonbondedQuadrupletEnergyForce(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int nparam); 
//     void BondedQuadrupletEnergyForce(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int nparam); 
    
    // Empirical potentials
    void EmpiricalPotentialEnergyForce(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int *nparam); 
    
    // Machine learning potentials
    void RadialSphericalHarmonicDescriptors(dstype *e, dstype *x, dstype *q, dstype *param, Int nparam);     
    void RadialSphericalHarmonicDescriptors(dstype *e, dstype *f, dstype *x, dstype *q, dstype *param, Int nparam);     
    void RadialSphericalHarmonicEnergyForce(dstype *e, dstype *f, dstype *x, dstype *coeff, dstype *q, dstype *param, Int nparam);     
};


#endif        

