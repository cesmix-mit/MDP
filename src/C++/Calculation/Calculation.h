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
    void NeighFullTriplets(dstype* x, dstype *q, Int *atomtype, Int istart, Int iend);     
    void NeighTriplets(dstype* x, dstype *q, Int *atomtype, Int typei, Int typej, Int typek, Int istart, Int iend);              
    void NeighFullTriplets(dstype* x, dstype *q, Int *atomtype, Int typei, Int typej, Int typek, Int istart, Int iend);              
  
//     void HalfForceDecomposition(T *fi, T *fij, int *ai, int *aj, int ijnum);
//     void FullForceDecomposition(T *fi, T *fij, int *ai, int *aj, int ijnum);
//     void FullAtomDecomposition(T *fi, T *fij, int *ilist, int *anumsum, int inum);
//     void HalfAtomDecomposition(T *fi, T *fij, int *ilist, int *anumsum, int *bnumsum, int *index, int inum);
//     
//     void ForceDecompositionTriplet(T *fi, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum);
//     void AtomDecompositionTriplet(T *fi, T *fij, T *fik, int *ilist, int *aj, int *ak, int *anumsum, int inum);    
};


#endif        

