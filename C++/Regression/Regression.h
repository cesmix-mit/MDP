/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

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
    void ValidateLinearRegression(CCalculation &CCal);    
    void ReferencePotentialDescriptors(CCalculation &CCal);    
};



#endif        

