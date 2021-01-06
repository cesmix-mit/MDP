#ifndef __DESCRIPTORS_H__
#define __DESCRIPTORS_H__

class CDescriptors {
private:
public:
    Int backend;         // computing platform
    neighborstruct nb;
    sysstruct sys;         
    CSphericalHarmonics shd;
                
    // constructor 
    CDescriptors(CConfiguration& conf, Int confignum)
        : shd(conf.common.K, conf.common.L, conf.common.backend)
        {backend = conf.common.backend;};    
        
    // destructor        
    ~CDescriptors(); 
    
    void NeighborList(CConfiguration& conf, Int confignum);
    void Energy(dstype* coeff);
    void Forces(dstype* coeff);
    void Stresses(dstype* coeff);    
    void BasisFunctions();
    void BasisFunctionsDeriv();
};


#endif        

