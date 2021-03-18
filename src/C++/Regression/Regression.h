#ifndef __REGRESSION_H__
#define __REGRESSION_H__

class CRegression {
private:
public:        
    // constructor 
    CRegression(CCalculation &CCal);       
    
    // destructor        
    ~CRegression(){  }; 
    
    void LinearRegression();
    void GaussianRegression();
    void NeuralNetRegression();
};


#endif        

