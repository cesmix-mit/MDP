#ifndef __REGRESSION_H__
#define __REGRESSION_H__

class CRegression {
private:
public:        
    // constructor 
    CRegression(CCalculation &CCal);       
    
    // destructor        
    ~CRegression(){  }; 
    
    void LinearRegression(CCalculation &CCal);
    void GaussianRegression(CCalculation &CCal);
    void NeuralNetRegression(CCalculation &CCal);
};


#endif        

