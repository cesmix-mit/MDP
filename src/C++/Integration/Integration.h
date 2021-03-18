#ifndef __INTEGRATION_H__
#define __INTEGRATION_H__

class CIntegration {
private:
public:
    // constructor 
    CIntegration(CCalculation &CCal);          
    
    // destructor        
    ~CIntegration(){  }; 
    
    void InitialCondition();
    void VelocityVerlet();
    void OutputSampling();
};



#endif        

