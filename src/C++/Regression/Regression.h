#ifndef __REGRESSION_H__
#define __REGRESSION_H__

class CRegression {
private:
public:
    CConfiguration cconf;   // configuration class
    CDescriptors cdesc;     // descriptors class
        
    // constructor 
    CRegression(string configfile, string appfile, Int mpiprocs, Int mpirank, Int backend)   
       : cconf(configfile, appfile, mpiprocs, mpirank, backend),
         cdesc(cconf){ };        
    
    // destructor        
    ~CRegression(){  }; 
    
    void LinearRegression();
    void GaussianRegression();
    void NeuralNetRegression();
};


#endif        

