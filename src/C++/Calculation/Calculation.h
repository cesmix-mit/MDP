#ifndef __CALCULATION_H__
#define __CALCULATION_H__

class CCalculation {
private:
public:                  
    CConfiguration cconf;   // configuration class
    CDescriptors cdesc;     // descriptors class
                
    // constructor 
    CCalculation(string filein, string fileout, Int mpiprocs, Int mpirank, Int backend)   
       : cconf(filein, fileout, mpiprocs, mpirank, backend),
         cdesc(cconf){ };        
        
    // destructor        
    ~CCalculation();    
    
    void NeighPairs(dstype* x, dstype *q, Int *atomtype, Int istart, Int iend);   
    void NeighPairs(dstype* x, dstype *q, Int *atomtype, Int typei, Int istart, Int iend);   
    void NeighPairs(dstype* x, dstype *q, Int *atomtype, Int typei, Int typej, Int istart, Int iend);   
        
    void NeighTriplets(dstype* x, dstype *q, Int *atomtype, Int istart, Int iend);      
    void NeighTriplets(dstype* x, dstype *q, Int *atomtype, Int typei, Int typej, Int typek, Int istart, Int iend);              
    
    
};


#endif        

