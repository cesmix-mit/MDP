#ifndef __REGRESSION_H__
#define __REGRESSION_H__

class CRegression {
private:
public:
    CConfiguration conf;   // configuration class
    CDescriptors desc;     // descriptors class
        
    // constructor 
    CRegression(string configfile, string appfile, Int mpiprocs, Int mpirank, Int backend)   
       : conf(configfile, appfile, mpiprocs, mpirank, backend),
         desc(conf, backend){ };        
    
    // destructor        
    ~CRegression(){  }; 
    
    void LinearRegression();
    void GaussianRegression();
    void NeuralNetRegression();
};


#endif        

