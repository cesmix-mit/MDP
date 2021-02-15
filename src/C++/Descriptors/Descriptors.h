#ifndef __DESCRIPTORS_H__
#define __DESCRIPTORS_H__

class CDescriptors {
private:
public:
    CSphericalHarmonics csh;
    
    CDescriptors(CConfiguration& cconf)
        : csh(cconf.common.K, cconf.common.L, cconf.common.backend) {  };    
        
    // destructor        
    ~CDescriptors(){};     
    
    
//     Int backend;      // computing platform: 1 -> CPU, 2-> (CUDA) GPU
//     Int descriptor;   // descriptor flag: 0 -> Spherical Harmonics Bessel
//     Int spectrum;     // spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
//     Int natomtypes;   // number of atom types;    
//     Int niatoms;      // number of i atoms
//     Int nijatoms;     // number of (i,j) atom pairs    
//     Int nimax;        // upper bound of niatoms
//     Int njmax=100;    // maximum number of neighbors
//     Int nbasis;       // number of basis functions
//     dstype rcut;      // cutoff radius
    
//     cudaEvent_t eventHandle;
//     cublasHandle_t cublasHandle;
//             
//     neighborstruct nb;
//     sysstruct sys;     
//     tempstruct tmp;    
//    CSphericalHarmonics csh;
    
//     void SetNeighborStruct(CConfiguration& cconf);
//     void SetSysStruct(CConfiguration& cconf);
//     void SetTempStruct(CConfiguration& cconf);
//     void Init(CConfiguration& cconf);
//     
//     void BuildVerletList(dstype* x, Int* atomtype, Int numatoms);        
//     void NeighborList(dstype* x, Int* atomtype, Int numatoms);    
//     int  NeighborPairs(dstype* x, Int* atomtype, Int *ilist, Int ilistsize, Int typei, Int typej);    
//     
//     void BasisFunctions(dstype* d, dstype* x, Int* atomtype, Int numatoms);
//     void BasisFunctionsDeriv(dstype* dd, dstype* x, Int* atomtype, Int numatoms);
//     void BasisFunctionsWithDeriv(dstype* d, dstype* dd, dstype* x, Int* atomtype, Int numatoms);
//     
//     void Energy(dstype e, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff);
//     void Forces(dstype* f, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff);
//     void Stresses(dstype* s, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff);
//     void EnergyForces(dstype e, dstype* f, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff);
//     void EnergyForcesStresses(dstype e, dstype* f, dstype* s, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff);
        
    // constructor 
//     CDescriptors(CConfiguration& cconf)
//         : csh(cconf.common.K, cconf.common.L, cconf.common.backend)
//         {
//             this->Init(cconf);                        
//         };    
//         
//     // destructor        
//     ~CDescriptors(){};     
};


#endif        

