#ifndef __INTEGRATION_H__
#define __INTEGRATION_H__

class CIntegration {
private:
public:
    CConfiguration cconf;   // configuration class
    CDescriptors cdesc;     // descriptors class
        
    // constructor 
    CIntegration(string configfile, string appfile, Int mpiprocs, Int mpirank, Int backend)   
       : cconf(configfile, appfile, mpiprocs, mpirank, backend),
         cdesc(cconf){ };        
    
    // destructor        
    ~CIntegration(){  }; 
    
    void InitialCondition();
    void VelocityVerlet();
    void OutputSampling();
};



#endif        

