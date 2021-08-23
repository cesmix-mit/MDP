/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __INTEGRATION_H__
#define __INTEGRATION_H__

class CIntegration {
private:
public:
    // constructor 
    CIntegration(CCalculation &CCal);          
    
    // destructor        
    ~CIntegration(){  }; 
    
    void VelocityVerlet(CCalculation &CCal);
    void IntegrationSetup(CCalculation &CCal);
    void InitialIntegration(CCalculation &CCal);
    void PostInitialIntegration(CCalculation &CCal);
    void DynamicLoadBalancing(CCalculation &CCal);
    void ProcessorsCommunication(CCalculation &CCal);
    void NeighborListRebuild(CCalculation &CCal, dstype *xstart, int istep);
    //void ForceComputation(CCalculation &CCal);
    void PostForceComputation(CCalculation &CCal);
    void FinalIntegration(CCalculation &CCal);
    void PostFinalIntegration(CCalculation &CCal);
    void TimestepCompletion(CCalculation &CCal);
};


#endif        

