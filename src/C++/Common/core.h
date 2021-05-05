#ifndef __CORE_H__
#define __CORE_H__

#include "cpucore.h"

#ifdef HAVE_CUDA      
#include "gpucore.h"
#endif                

inline void Kron(dstype *C, dstype *A, dstype *B, Int M1, Int M2, Int backend)
{
    if (backend == 1)  
        cpuKron(C, A, B, M1, M2);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuKron(C, A, B, M1, M2);
#endif                  
}

inline void GetArrayAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
    if (backend == 1)  
        cpuGetArrayAtIndex(y, x, ind, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuGetArrayAtIndex(y, x, ind, n);
#endif                  
}

inline void PutArrayAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
    if (backend == 1)  
        cpuPutArrayAtIndex(y, x, ind, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuPutArrayAtIndex(y, x, ind, n);
#endif                  
}

inline void ArrayPlusXAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayPlusXAtIndex(y, x, ind, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayPlusXAtIndex(y, x, ind, n);
#endif                  
}

inline void ArrayMinusXAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayMinusXAtIndex(y, x, ind, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayMinusXAtIndex(y, x, ind, n);
#endif                  
}

inline void ArraySquare(dstype *y, dstype *x, Int n, Int backend)
{
    if (backend == 1)  
        cpuArraySquare(y, x, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArraySquare(y, x, n);
#endif                  
}

inline void ArraySetValue(dstype *y, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArraySetValue(y, a, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArraySetValue(y, a, n);
#endif                  
}

inline void ArraySetValueAtIndex(dstype *y, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        y[n] = a;
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArraySetValueAtIndex(y, a, n);
#endif                  
}

inline dstype ArrayGetValueAtIndex(dstype *y, Int n, Int backend)
{
    dstype val; 
    if (backend == 1)  
        val = y[n];
#ifdef HAVE_CUDA            
    if (backend == 2)  
        val = gpuArrayGetValueAtIndex(y, n);
#endif                      
    return val;
}

inline dstype ArraySum(dstype *a, dstype *b, Int n, Int backend)
{
    dstype val; 
    if (backend == 1)  
        val = cpuArraySum(a, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        val = gpuArraySum(a, b, n);
#endif                      
    return val;
}

inline Int IntArrayGetValueAtIndex(Int *y, Int n, Int backend)
{
    Int val; 
    if (backend == 1)  
        val = y[n];
#ifdef HAVE_CUDA            
     if (backend == 2)  
         val = gpuArrayGetValueAtIndex(y, n);
#endif                      
    return val;
}

inline void ArrayFill(int *a, int start, int n, Int backend)
{
    if (backend == 1)  
        cpuArrayFill(a, start, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayFill(a, start, n);
#endif                  
}

inline void ArrayTripletnum(int* output, int* input, int length, int backend) 
{
    if (backend == 1)  
        cpuTripletnum(output, input, length); 
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuTripletnum(output, input, length); 
#endif                  
}
inline void ArrayQuadrupletnum(int* output, int* input, int length, int backend) 
{
    if (backend == 1)  
        cpuQuadrupletnum(output, input, length); 
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuQuadrupletnum(output, input, length); 
#endif                  
}


inline void ArrayAddScalar(dstype *y, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAddScalar(y, a, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAddScalar(y, a, n);
#endif                  
}

inline void ArrayMultiplyScalar(dstype *y, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayMultiplyScalar(y, a, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayMultiplyScalar(y, a, n);
#endif                  
}

inline void ArrayCopy(dstype *y, dstype *x, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayCopy(y, x, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayCopy(y, x, n);
#endif                  
}

inline void ArrayCopy(Int *y, Int *x, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayCopy(y, x, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayCopy(y, x, n);
#endif                  
}

inline void ArrayMinus(dstype *y, dstype *x, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayMinus(y, x, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayMinus(y, x, n);
#endif                  
}

inline void ArrayAbs(dstype *y, dstype *x, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAbs(y, x, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAbs(y, x, n);
#endif                  
}

inline void ArraySqrt(dstype *y, dstype *x, Int n, Int backend)
{
    if (backend == 1)  
        cpuArraySqrt(y, x, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArraySqrt(y, x, n);
#endif                  
}

inline void ArrayMultiplyScalarDiagonal(dstype *C, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayMultiplyScalarDiagonal(C, a, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayMultiplyScalarDiagonal(C, a, n);
#endif                  
}

inline void ArrayAddVectorToDiagonal(dstype *C, dstype *x, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAddVectorToDiagonal(C, x, a, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAddVectorToDiagonal(C, x, a, n);
#endif                  
}

inline void ArrayAXPB(dstype *y, dstype *x, dstype a, dstype b, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAXPB(y, x, a, b, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAXPB(y, x, a, b, n);
#endif                  
}

inline void ArrayAXPBY(dstype *z, dstype *x, dstype *y, dstype a, dstype b, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAXPBY(z, x, y, a, b, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAXPBY(z, x, y, a, b, n);
#endif                  
}

inline void ArrayAXY(dstype *z, dstype *x, dstype *y, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAXY(z, x, y, a, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAXY(z, x, y, a, n);
#endif                  
}

inline void ArrayAXYZ(dstype *s, dstype *x, dstype *y, dstype *z, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAXYZ(s, x, y, z, a, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAXYZ(s, x, y, z, a, n);
#endif                  
}

inline void ArrayAXYPBZ(dstype *s, dstype *x, dstype *y, dstype *z, dstype a, dstype b, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAXYPBZ(s, x, y, z, a, b, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAXYPBZ(s, x, y, z, a, b, n);
#endif                  
}

inline void ArrayAdd3Vectors(dstype *s, dstype *x, dstype *y, dstype *z, dstype a, dstype b, dstype c, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAdd3Vectors(s, x, y, z, a, b, c, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAdd3Vectors(s, x, y, z, a, b, c, n);
#endif                  
}

inline void ArrayExtract(dstype *un, dstype *u, Int I, Int J, Int K, 
        Int i1, Int i2, Int j1, Int j2, Int k1, Int k2, Int backend)
{
    if (backend == 1)  
        cpuArrayExtract(un, u, I, J, K, i1, i2, j1, j2, k1, k2);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayExtract(un, u, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                  
}

inline void ArrayInsert(dstype *u, dstype *un, Int I, Int J, Int K, 
        Int i1, Int i2, Int j1, Int j2, Int k1, Int k2, Int backend)
{
    if (backend == 1)  
        cpuArrayInsert(u, un, I, J, K, i1, i2, j1, j2, k1, k2);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayInsert(u, un, I, J, K, i1, i2, j1, j2, k1, k2);    
#endif                  
}

inline void ArrayGemmBatch(dstype *C, dstype *A, dstype *B, Int I, Int J, Int K, Int S, Int backend)
{
    if (backend == 1)  
        cpuArrayGemmBatch(C, A, B, I, J, K, S);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayGemmBatch(C, A, B, I, J, K, S);    
#endif                  
}

inline void ArrayGemmBatch1(dstype *C, dstype *A, dstype *B, Int I, Int J, Int K, Int S, Int backend)
{
    if (backend == 1)  
        cpuArrayGemmBatch1(C, A, B, I, J, K, S);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayGemmBatch1(C, A, B, I, J, K, S);    
#endif                  
}

inline void ArrayDG2CG(dstype *ucg, dstype *udg, Int *cgent2dgent, Int *rowent2elem, Int nent, Int backend)
{
    if (backend == 1)  
        cpuArrayDG2CG(ucg, udg, cgent2dgent, rowent2elem, nent);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayDG2CG(ucg, udg, cgent2dgent, rowent2elem, nent);    
#endif                  
}

inline void ArrayDG2CG2(dstype *ucg, dstype *udg, Int *colent2elem, Int *rowent2elem, Int nent, Int npe, Int backend)
{
    if (backend == 1)  
        cpuArrayDG2CG2(ucg, udg, colent2elem, rowent2elem, nent, npe);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayDG2CG2(ucg, udg, colent2elem, rowent2elem, nent, npe);    
#endif                  
}

inline void Cart2Sphere(dstype *the, dstype *phi, dstype *r, dstype *x, dstype *y, dstype *z, Int N, Int backend)
{
    if (backend == 1)  
        cpuCart2Sphere(the, phi, r, x, y, z, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuCart2Sphere(the, phi, r, x, y, z, N);    
#endif                  
}

inline void Cart2SphereDeriv(dstype *the, dstype *phi, dstype *r, dstype *thex, dstype *they, dstype *thez, dstype *phix, 
        dstype *phiy, dstype *phiz, dstype *rx, dstype *ry, dstype *rz, dstype *x, dstype *y, dstype *z, Int N, Int backend)
{
    if (backend == 1)  
        cpuCart2SphereDeriv(the, phi, r, thex, they, thez, phix, phiy, phiz,
            rx, ry, rz, x, y, z, N);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuCart2SphereDeriv(the, phi, r, thex, they, thez, phix, phiy, phiz,
            rx, ry, rz, x, y, z, N);
#endif                  
}

inline void Sphere2Cart(dstype *x, dstype *y, dstype *z, dstype *the, dstype *phi, dstype *r, Int N, Int backend)
{
    if (backend == 1)  
        cpuSphere2Cart(x, y, z, the, phi, r, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuSphere2Cart(x, y, z, the, phi, r, N);    
#endif                  
}

inline void Euler2Rotm(dstype *R11, dstype *R12, dstype *R13, dstype *R21, 
                dstype *R22, dstype *R23, dstype *R31, dstype *R32, dstype *R33, dstype *alpha, dstype *beta, dstype *gamma, Int N, Int backend)
{
    if (backend == 1)  
        cpuEuler2Rotm(R11, R12, R13, R21, R22, R23, R31, R32, R33, alpha, beta, gamma, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuEuler2Rotm(R11, R12, R13, R21, R22, R23, R31, R32, R33, alpha, beta, gamma, N);    
#endif                  
}

inline void Rotc(dstype *X, dstype *Y, dstype *Z, dstype *R, dstype *x, dstype *y, dstype *z, Int N, Int backend)
{
    if (backend == 1)  
        cpuRotc(X, Y, Z, R, x, y, z, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRotc(X, Y, Z, R, x, y, z, N);    
#endif                  
}

inline void Cumsum(int *d_out, int *d_in, int *d_sums, int *d_incr, int length, Int backend)
{
    if (backend == 1)  
        cpuCumsum(d_out, d_in, length);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuCumsum(d_out, d_in, d_sums, d_incr, length);
#endif                      
}

inline void AtomList2D(int *alist,  int *inside, int *glistnumsum, int *glistnum, int *d_sums, int *d_incr, 
        dstype *x, dstype *pimages, dstype *wc, dstype *s2rmap, int inum, int pnum, int dim, Int backend)
{
    if (backend == 1) { 
        cpuAtomList2D(alist, inside, glistnumsum, glistnum, x, pimages, wc, 
                s2rmap, inum, pnum, dim);      
    }
#ifdef HAVE_CUDA            
    if (backend == 2) { 
        gpuAtomList2D(alist, inside, glistnumsum, glistnum, d_sums, d_incr, x, pimages, 
                wc, s2rmap, inum, pnum, dim);
    }
#endif                          
}        

inline void AtomList3D(int *alist,  int *inside, int *glistnumsum, int *glistnum, int *d_sums, int *d_incr, 
        dstype *x, dstype *pimages, dstype *wc, dstype *s2rmap, int inum, int pnum, int dim, int backend)
{
    if (backend == 1) { 
        cpuAtomList3D(alist, inside, glistnumsum, glistnum, x, pimages, wc, 
                s2rmap, inum, pnum, dim);      
    }
#ifdef HAVE_CUDA            
    if (backend == 2) { 
        gpuAtomList3D(alist, inside, glistnumsum, glistnum, d_sums, d_incr, x, pimages, 
                wc, s2rmap, inum, pnum, dim);
    }
#endif                          
}        
        
inline void FullNeighborList2D(int *neighlist, int *neighnum, dstype *x, dstype *rcutsq, 
        int anum, int inum, int jnum, int dim, int backend)
{
    if (backend == 1) { 
        cpuFullNeighborList2D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);
    }
#ifdef HAVE_CUDA            
    if (backend == 2) { 
        gpuFullNeighborList2D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);        
    }
#endif                          
}                

inline void FullNeighborList3D(int *neighlist, int *neighnum, dstype *x, dstype *rcutsq, 
        int anum, int inum, int jnum, int dim, int backend)
{
    if (backend == 1) { 
        cpuFullNeighborList3D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);
    }
#ifdef HAVE_CUDA            
    if (backend == 2) { 
        gpuFullNeighborList3D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);        
    }
#endif                          
}                

inline void CellList2D(int *clist, dstype *x, dstype *eta1, dstype *eta2, dstype *eta3, dstype *s2rmap, 
        int *nc, int inum, int natom, int dim, int backend)
{
    if (backend == 1) { 
        cpuCellList2D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);  
    }
#ifdef HAVE_CUDA            
    if (backend == 2) { 
        gpuCellList2D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);  
    }
#endif                          
}                

inline void CellList3D(int *clist, dstype *x, dstype *eta1, dstype *eta2, dstype *eta3, dstype *s2rmap, 
        int *nc, int inum, int natom, int dim, int backend)
{
    if (backend == 1) { 
        cpuCellList3D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);  
    }
#ifdef HAVE_CUDA            
    if (backend == 2) { 
        gpuCellList3D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);  
    }
#endif                          
}                

inline void Cell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int *dtemp, 
        int natom, int ncell, int backend)
{
    if (backend == 1) { 
        cpuCell2AtomList(c2alist, c2anumsum, c2anum, clist, natom, ncell);    
    }
#ifdef HAVE_CUDA            
    if (backend == 2) { 
        gpuCell2AtomList(c2alist, c2anumsum, c2anum, clist, dtemp, natom, ncell);    
    }
#endif                          
}                

inline void FullNeighborList2D(int *neighlist, int *neighnum, dstype *x, dstype *rcutsq,
        int *alist, int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim, int backend)
{
    if (backend == 1) { 
        cpuFullNeighborList2D(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);     
    }
#ifdef HAVE_CUDA            
    if (backend == 2) { 
        gpuFullNeighborList2D(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);     
    }
#endif                          
}                

inline void FullNeighborList3D(int *neighlist, int *neighnum, dstype *x, dstype *rcutsq,
        int *alist, int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim, int backend)
{
    if (backend == 1) { 
        cpuFullNeighborList3D(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);     
    }
#ifdef HAVE_CUDA            
    if (backend == 2) { 
        gpuFullNeighborList3D(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);     
    }
#endif                          
}                

inline int FindAtomType(int *ilist, int* olist, int *atomtype, int *t0, int *t1, int typei, int na, int backend)
{
    if (backend == 1) { 
        return cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
    }
#ifdef HAVE_CUDA            
    if (backend == 2) { 
        return gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
    }
#endif                          
    return 0;
}                

inline void FullNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim, Int backend)
{
    if (backend == 1)  
        cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, neighnum, inum, jnum, dim);        
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, neighnum, inum, jnum, dim);        
#endif                  
}

inline void FullNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *atomtype, int *ilist, 
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim, Int backend)
{
    if (backend == 1)  
        cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#endif                  
}

inline void FullNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *atomtype, int *ilist, 
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int typek, int dim, Int backend)
{
    if (backend == 1)  
        cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, dim);        
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, dim);        
#endif                  
}

inline void FullNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *atomtype, int *ilist, 
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int typek, int typel, int dim, Int backend)
{
    if (backend == 1)  
        cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, typel, dim);        
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, typel, dim);        
#endif                  
}

inline void HalfNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *ilist, int *alist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim, Int backend)
{
    if (backend == 1)  
        cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, alist, neighlist, neighnum, inum, jnum, dim);        
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, alist, neighlist, neighnum, inum, jnum, dim);        
#endif                  
}

inline void HalfNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *atomtype, int *ilist, 
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim, Int backend)
{
    if (backend == 1)  
        cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#endif                  
}

inline void NeighTripletList(int *tripletlist, int *tripletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum, Int backend)
{
    if (backend == 1)  
        cpuNeighTripletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuNeighTripletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#endif                  
}

inline void NeighQuadrupletList(int *tripletlist, int *tripletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum, Int backend)
{
    if (backend == 1)  
        cpuNeighQuadrupletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuNeighQuadrupletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#endif                  
}

// template <typename T> void cpuHalfNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *alist, 
//         int *neighlist, int *neighnum, int inum, int jnum, int dim);
// template <typename T> void cpuHalfNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *atomtype, int *ilist, 
//         int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim);

inline void NeighSingles(dstype *xi, dstype *qi, dstype *x, dstype *q, int *ai, 
      int *ti, int *ilist, int *atomtype, int inum, int ncq, int dim, Int backend)
{
    if (backend == 1)  
        cpuNeighSingles(xi, qi, x, q, ai, ti, ilist, atomtype, inum, ncq, dim);
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuNeighSingles(xi, qi, x, q, ai, ti, ilist, atomtype, inum, ncq, dim);
    }
#endif             
}

inline void NeighPairs(dstype *xij, dstype *qi, dstype *qj, dstype *x, dstype *q, int *ai, int *aj,  
      int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, int *alist,  
      int *atomtype, int inum, int jnum, int ncq, int dim, int backend)
{
    if (backend == 1)  
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, alist, 
                atomtype, inum, jnum, ncq, dim);       
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, alist, 
                atomtype, inum, jnum, ncq, dim);       
#endif                  
}

inline void NeighTriplets(dstype *xij, dstype *xik, dstype *qi, dstype *qj, dstype *qk, dstype *x, dstype *q, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *ilist, int *alist,  
      int *atomtype, int inum, int ncq, int dim, int backend)
{
    if (backend == 1)  
        cpuNeighTriplets(xij, xik, qi, qj, qk, x, q, ai, aj, ak, ti, tj, tk, tripletnum, tripletlist, tripletnumsum, 
                ilist, alist, atomtype, inum, ncq, dim);       
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuNeighTriplets(xij, xik, qi, qj, qk, x, q, ai, aj, ak, ti, tj, tk, tripletnum, tripletlist, tripletnumsum, 
                ilist, alist, atomtype, inum, ncq, dim);       
#endif                  
}

inline void NeighQuadruplets(dstype *xij, dstype *xik, dstype *xil, dstype *qi, dstype *qj, dstype *qk, dstype *ql, dstype *x, dstype *q, 
      int *ai, int *aj, int *ak, int *al, int *ti, int *tj, int *tk, int *tl, int *quadrupletnum, int *quadrupletlist, 
      int *quadrupletnumsum, int *ilist, int *alist,  int *atomtype, int inum, int ncq, int dim, int backend)
{
    if (backend == 1)  
        cpuNeighQuadruplets(xij, xik, xil, qi, qj, qk, ql, x, q, ai, aj, ak, al, ti, tj, tk, tl, 
                quadrupletnum, quadrupletlist, quadrupletnumsum, ilist, alist, atomtype, inum, ncq, dim);       
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuNeighQuadruplets(xij, xik, xil, qi, qj, qk, ql, x, q, ai, aj, ak, al, ti, tj, tk, tl, 
                quadrupletnum, quadrupletlist, quadrupletnumsum, ilist, alist, atomtype, inum, ncq, dim);       
#endif                  
}

inline void SphericalHarmonicsBessel(dstype *Sr, dstype *Si, dstype *xij, dstype *x0, 
        dstype *P, dstype *tmp, dstype *f, dstype *fac, dstype pi, Int L, Int K, Int N, Int backend)
{
    if (backend == 1)  
        cpuSphericalHarmonicsBessel(Sr, Si, xij, x0, P, tmp, f, fac, pi, L, K, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuSphericalHarmonicsBessel(Sr, Si, xij, x0, P, tmp, f, fac, pi, L, K, N);    
    }
#endif                  
}

inline void RadialSphericalHarmonicsSpectrum(dstype *c, dstype *ar, dstype *ai, dstype *sr, dstype *si, dstype *cg, int *indk, 
        int *indl, int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsSpectrum(c, ar, ai, sr, si, cg, indk, indl, indm, 
                rowm, Nnb, Nub, Ncg, Na, L, K, spectrum); 
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuRadialSphericalHarmonicsSpectrum(c, ar, ai, sr, si, cg, indk, indl, indm, 
                rowm, Nnb, Nub, Ncg, Na, L, K, spectrum); 
    }
#endif                  
}

inline void SphericalHarmonicsBesselWithDeriv(dstype *Sr, dstype *Si, dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, dstype *xij, 
                dstype *x0, dstype *P, dstype *tmp, dstype *f, dstype *dP, dstype *dtmp, dstype *df, dstype *fac, dstype pi, int L, int K, int N, Int backend)
{
    if (backend == 1)  
        cpuSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, xij,
                 x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, xij,
                 x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);
    }
#endif                  
}

inline void RadialSphericalHarmonicsSpectrumDeriv(dstype *cd, dstype *ar, dstype *ai, 
        dstype *srx, dstype *six, dstype *sry, dstype *siy, dstype *srz, dstype *siz, dstype *cg, int *indk, int *indl, int *indm, 
        int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsSpectrumDeriv(cd, ar, ai, srx, six, sry, siy, srz, siz, cg, indk, 
                indl, indm, rowm, Nnb, Nub, Ncg, Na, L, K, spectrum);                        
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuRadialSphericalHarmonicsSpectrumDeriv(cd, ar, ai, srx, six, sry, siy, srz, siz, cg, indk, 
                indl, indm, rowm, Nnb, Nub, Ncg, Na, L, K, spectrum);                        
    }
#endif                  
}

inline void ForceDecomposition(dstype *f, dstype *fij, int *ai, int *aj, int inum, int ijnum, int Nbf, Int backend)
{
    if (backend == 1)  
        cpuForceDecomposition(f, fij, ai, aj, inum, ijnum, Nbf);
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuForceDecomposition(f, fij, ai, aj, inum, ijnum, Nbf);
    }
#endif                  
}

inline void ForceDecomposition(dstype *f, dstype *fij, int *ai, int *aj, int ijnum, Int backend)
{
    if (backend == 1)  
        cpuForceDecomposition(f, fij, ai, aj, ijnum);
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuForceDecomposition(f, fij, ai, aj, ijnum);
    }
#endif                  
}

inline void CenterAtomDecomposition(dstype *f, dstype *fij, int *ilist, int *anumsum, 
        int inum, int ijnum, int Na, int Nbf, Int backend)
{
    if (backend == 1)  
        cpuCenterAtomDecomposition(f, fij, ilist, anumsum, inum, ijnum, Na, Nbf);
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuCenterAtomDecomposition(f, fij, ilist, anumsum, inum, ijnum, Na, Nbf);
    }
#endif                  
}

inline void CenterAtomDecomposition(dstype *e, dstype *ei, int *ilist, int Na, Int backend)
{
    if (backend == 1)  
        cpuCenterAtomDecomposition(e, ei, ilist, Na);
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuCenterAtomDecomposition(e, ei, ilist, Na);
    }
#endif                  
}

inline void CenterAtomDecomposition(dstype *f, dstype *fij, int *ilist, int *anumsum, 
        int inum, Int backend)
{
    if (backend == 1)  
        cpuCenterAtomDecomposition(f, fij, ilist, anumsum, inum);
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuCenterAtomDecomposition(f, fij, ilist, anumsum, inum);
    }
#endif                  
}

inline void NeighborAtomDecomposition(dstype *f, dstype *fij, int *ilist, int *anumsum, int *index, 
        int inum, int ijnum, int jnum, int Nbf, Int backend)
{
    if (backend == 1)  
        cpuNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum, ijnum, jnum, Nbf);
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum, ijnum, jnum, Nbf);
    }
#endif                  
}

inline void NeighborAtomDecomposition(dstype *f, dstype *fij, int *ilist, int *anumsum, int *index, 
        int inum, Int backend)
{
    if (backend == 1)  
        cpuNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum);
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        gpuNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum);
    }
#endif                  
}

inline int UniqueSort(int *b, int *c, int *d, int *e, int *a, int *p, int *t, int *q, int n, Int backend)
{
    if (backend == 1)  
        return cpuUniqueSort(b, c, d, e, a, p, t, q, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  {
        return gpuUniqueSort(b, c, d, e, a, p, t, q, n);
    }
#endif             
    return 0;
}

// inline void FullNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *ilist, int *neighlist, 
//         int *neighnum, int inum, int jnum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, neighnum, inum, jnum, dim);                
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
//         gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, neighnum, inum, jnum, dim);                
//     }
// #endif             
// }
// 
// inline void FullNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *atomtype, int *ilist, 
//         int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);                
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
//         gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);                
//     }
// #endif             
// }

// ---------- //

// inline void GhostAtoms(int *glistnum, dstype *x, dstype *pimages, dstype *wc, dstype *B2C, int n, int m, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuGhostAtoms(glistnum, x, pimages, wc, B2C, n, m, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuGhostAtoms(glistnum, x, pimages, wc, B2C, n, m, dim);    
// // #endif                  
// }
// 
// inline void CreateAtomList(int *ilist, int *glistnumsum, int *glistnum, int *atomtype, dstype *x, dstype 
//         *pimages, dstype *wc, dstype *B2C, int n, int m, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuCreateAtomList(ilist, glistnumsum, glistnum, atomtype, x, pimages, wc, B2C, n, m, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuCreateAtomList(ilist, glistnumsum, glistnum, atomtype, x, pimages, wc, B2C, n, m, dim);    
// // #endif                  
// }
// 
// inline void CellList(int *clist, int *c2inum, dstype *x, dstype *eta1, dstype *eta2, dstype *eta3, dstype *B2C, int *nc, 
//         int inum, int pnum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuCellList(clist, c2inum, x, eta1, eta2, eta3, B2C, nc, inum, pnum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuCellList(clist, c2inum, x, eta1, eta2, eta3, B2C, nc, inum, pnum, dim);    
// // #endif                  
// }
// 
// inline void Cell2AtomList(int *c2ilist, int *c2inumsum, int *c2inum, int *clist, int *nc, 
//         int inum, int pnum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuCell2AtomList(c2ilist, c2inumsum, c2inum, clist, nc, inum, pnum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuCell2AtomList(c2ilist, c2inumsum, c2inum, clist, nc, inum, pnum, dim);    
// // #endif                  
// }
// 
// inline void VerletAtoms(int *verletnum, dstype *x, dstype *ellipsoid, int *ilist, int *clist, 
//         int *c2ilist, int *c2inum, int *c2inumsum, int *nc, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuVerletAtoms(verletnum, x, ellipsoid, ilist, clist, c2ilist, c2inum, c2inumsum, 
//                 nc, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuVerletAtoms(verletnum, x, ellipsoid, ilist, clist, c2ilist, c2inum, c2inumsum, 
// //                 nc, inum, dim);    
// // #endif                  
// }
// 
// inline void CreateVerletList(int *verletlist, dstype *x, dstype *ellipsoid, int *verletnum,  int *verletnumsum, 
//     int *ilist, int *clist, int *c2ilist, int *c2inum, int *c2inumsum, int *nc, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuCreateVerletList(verletlist, x, ellipsoid, verletnum, verletnumsum, ilist, clist, c2ilist, 
//                 c2inum, c2inumsum, nc, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuCreateVerletList(verletlist, x, ellipsoid, verletnum, verletnumsum, ilist, clist, c2ilist, 
// //                 c2inum, c2inumsum, nc, inum, dim);     
// // #endif                  
// }
// 
// inline void FullNeighNum(int *neighnum, dstype *x, dstype* ellipsoid, 
//          int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuFullNeighNum(neighnum, x, ellipsoid, ilist, verletlist, verletnum, verletnumsum, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuFullNeighNum(neighnum, x, ellipsoid, ilist, verletlist, verletnum, verletnumsum, inum, dim);    
// // #endif                  
// }
// 
// inline void FullNeighList(int *neighlist, dstype *x, dstype* ellipsoid, int *neighnum, int *neighnumsum,
//          int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuFullNeighList(neighlist, x, ellipsoid, neighnum, neighnumsum, ilist, verletlist, verletnum, verletnumsum, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuFullNeighList(neighlist, x, ellipsoid, neighnum, neighnumsum, ilist, verletlist, verletnum, verletnumsum, inum, dim);    
// // #endif                  
// }
// 
// inline void HalfNeighNum(int *neighnum, dstype *x, dstype* ellipsoid, 
//          int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuHalfNeighNum(neighnum, x, ellipsoid, ilist, verletlist, verletnum, verletnumsum, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuHalfNeighNum(neighnum, x, ellipsoid, ilist, verletlist, verletnum, verletnumsum, inum, dim);    
// // #endif                  
// }
// 
// inline void HalfNeighList(int *neighlist, dstype *x, dstype* ellipsoid, int *neighnum, int *neighnumsum,
//          int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuHalfNeighList(neighlist, x, ellipsoid, neighnum, neighnumsum, ilist, verletlist, verletnum, verletnumsum, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuHalfNeighList(neighlist, x, ellipsoid, neighnum, neighnumsum, ilist, verletlist, verletnum, verletnumsum, inum, dim);    
// // #endif                  
// }
// 
// inline void VerletAtoms(int *verletnum, dstype *x, dstype *ellipsoid, int *atomtype, int *ilist, int *clist, 
//         int *c2ilist, int *c2inum, int *c2inumsum, int *nc, int ntype, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuVerletAtoms(verletnum, x, ellipsoid, atomtype, ilist, clist, c2ilist, c2inum, c2inumsum, 
//                 nc, ntype, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuVerletAtoms(verletnum, x, ellipsoid, atomtype, ilist, clist, c2ilist, c2inum, c2inumsum, 
// //                 nc, ntype, inum, dim);    
// // #endif                  
// }
// 
// inline void CreateVerletList(int *verletlist, dstype *x, dstype *ellipsoid, int *verletnum,  int *verletnumsum, int *atomtype, 
//     int *ilist, int *clist, int *c2ilist, int *c2inum, int *c2inumsum, int *nc, int ntype, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuCreateVerletList(verletlist, x, ellipsoid, verletnum, verletnumsum, atomtype, ilist, clist, c2ilist, 
//                 c2inum, c2inumsum, nc, ntype, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuCreateVerletList(verletlist, x, ellipsoid, verletnum, verletnumsum, atomtype, ilist, clist, c2ilist, 
// //                 c2inum, c2inumsum, nc, ntype, inum, dim);     
// // #endif                  
// }
// 
// inline void FullNeighNum(int *neighnum, dstype *x, dstype* ellipsoid, int *atomtype, 
//          int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuFullNeighNum(neighnum, x, ellipsoid, atomtype, ilist, verletlist, verletnum, verletnumsum, ntype, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuFullNeighNum(neighnum, x, ellipsoid, atomtype, ilist, verletlist, verletnum, verletnumsum, ntype, inum, dim);    
// // #endif                  
// }
// 
// inline void FullNeighList(int *neighlist, dstype *x, dstype* ellipsoid, int *neighnum, int *neighnumsum,
//          int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuFullNeighList(neighlist, x, ellipsoid, neighnum, neighnumsum, atomtype, ilist, verletlist, verletnum, verletnumsum, ntype, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuFullNeighList(neighlist, x, ellipsoid, neighnum, neighnumsum, atomtype, ilist, verletlist, verletnum, verletnumsum, ntype, inum, dim);    
// // #endif                  
// }
// 
// inline void HalfNeighNum(int *neighnum, dstype *x, dstype* ellipsoid, int *atomtype, 
//          int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuHalfNeighNum(neighnum, x, ellipsoid, atomtype, ilist, verletlist, verletnum, verletnumsum, ntype, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuHalfNeighNum(neighnum, x, ellipsoid, atomtype, ilist, verletlist, verletnum, verletnumsum, ntype, inum, dim);    
// // #endif                  
// }
// 
// inline void HalfNeighList(int *neighlist, dstype *x, dstype* ellipsoid, int *neighnum, int *neighnumsum,
//          int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuHalfNeighList(neighlist, x, ellipsoid, neighnum, neighnumsum, atomtype, ilist, verletlist, verletnum, verletnumsum, ntype, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuHalfNeighList(neighlist, x, ellipsoid, neighnum, neighnumsum, atomtype, ilist, verletlist, verletnum, verletnumsum, ntype, inum, dim);    
// // #endif                  
// }
// 
// inline void GetNeighPairs(dstype *xij, dstype *x, int *ti, int *tj, int *ilist, int *neighlist,  
//          int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuGetNeighPairs(xij, x, ti, tj, ilist, neighlist, neighnum, neighnumsum, atomtype, ntype, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  
// //         gpuGetNeighPairs(xij, x, ti, tj, ilist, neighlist, neighnum, neighnumsum, atomtype, ntype, inum, dim);    
// // #endif                  
// }
// 
// inline void GetNeighPairs(dstype *xij, dstype *x, int *anum, int *anumsum, int *ai, int *aj, int *ti, 
//       int *tj, int *ilist, int *neighlist, int *neighnum, int *neighnumsum, int *atomtype, 
//       int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuGetNeighPairs(xij, x, anum, anumsum, ai, aj, ti, tj, ilist, neighlist, neighnum, 
//                 neighnumsum, atomtype, inum, dim);    
// // #ifdef HAVE_CUDA            
// //     if (backend == 2)  {
// //         //gpuGetNeighPairs(xij, x, anum, anumsum, ai, aj, ti, tj, ilist, neighlist, neighnum, 
// //        //         neighnumsum, atomtype, inum, dim);    
// //     }
// // #endif                  
// }
// 
// inline int GetNeighPairs(dstype *xij, dstype *x, int *anum, int *anumsum, int *ai, int *aj, 
//       int *ti, int *tj, int *ilist, int *tlist, int *neighlist, int *neighnum, int *neighnumsum, 
//       int *atomtype, int typei, int inum, int dim, Int backend)
// {
//     int tnum=0;
//     if (backend == 1)  
//         tnum = cpuGetNeighPairs(xij, x, anum, anumsum, ai, aj, ti, tj, ilist, tlist, neighlist, neighnum, 
//                 neighnumsum, atomtype, typei, inum, dim);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //         tnum = gpuGetNeighPairs(xij, x, anum, anumsum, ai, aj, ti, tj, ilist, tlist, neighlist, neighnum, 
// //                 neighnumsum, atomtype, typei, inum, dim);    
//     }
// #endif        
//     
//     return tnum;
// }
// 
// inline int GetNeighPairs(dstype *xij, dstype *x, int *anum, int *anumsum, int *ai, int *aj, 
//       int *ti, int *tj, int *ilist, int *tlist, int *neighlist, int *neighnum, int *neighnumsum, 
//       int *atomtype, int typei, int typej, int inum, int dim, Int backend)
// {
//     int tnum=0;
//     if (backend == 1)  
//         tnum = cpuGetNeighPairs(xij, x, anum, anumsum, ai, aj, ti, tj, ilist, tlist, neighlist, neighnum, 
//                 neighnumsum, atomtype, typei, typej, inum, dim);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //         tnum = gpuGetNeighPairs(xij, x, anum, anumsum, ai, aj, ti, tj, ilist, tlist, neighlist, neighnum, 
// //                 neighnumsum, atomtype, typei, typej, inum, dim);    
//     }
// #endif         
//     
//     return tnum;
// }
// 
// inline void GetNeighPairs(dstype *xi, dstype *xj, dstype *x, int *ti, int *tj, int *ilist, int *neighlist,  
//          int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim, Int backend)
// {
//     if (backend == 1)  
//         cpuGetNeighPairs(xi, xj, x, ti, tj, ilist, neighlist, neighnum, neighnumsum, atomtype, ntype, inum, dim);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
//         //gpuGetNeighPairs(xi, xj, x, ti, tj, ilist, neighlist, neighnum, neighnumsum, atomtype, ntype, inum, dim);    
//     }
// #endif                  
// }
// 
// inline void coreSphericalHarmonicsBessel(dstype *Sr, dstype *Si, dstype *x, dstype *y, dstype *z, 
//                 dstype *x0, dstype *P, dstype *tmp, dstype *f, dstype *fac, dstype pi, Int L, Int K, Int N, Int backend)
// {
//     if (backend == 1)  
//         cpuSphericalHarmonicsBessel(Sr, Si, x, y, z, x0, P, tmp, f, fac, pi, L, K, N);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
//         //gpuSphericalHarmonicsBessel(Sr, Si, x, y, z, x0, P, tmp, f, fac, pi, L, K, N);    
//     }
// #endif                  
// }
// 
// 
// inline void coreSphericalHarmonicsBesselDeriv(dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, 
//       dstype *x, dstype *y, dstype *z, dstype *x0, dstype *P, dstype *tmp, dstype *f, dstype *dP, dstype *dtmp, dstype *df, dstype *fac, dstype pi, Int L, Int K, Int N, Int backend)
// {
//     if (backend == 1)  
//         cpuSphericalHarmonicsBesselDeriv(Srx, Six, Sry, Siy, Srz, Siz, x, y, z, 
//             x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //         gpuSphericalHarmonicsBesselDeriv(Srx, Six, Sry, Siy, Srz, Siz, x, y, z, 
// //             x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);    
//     }
// #endif                  
// }
// 
// inline void coreSphericalHarmonicsBesselWithDeriv(dstype *Sr, dstype *Si, dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, 
//       dstype *x, dstype *y, dstype *z, dstype *x0, dstype *P, dstype *tmp, dstype *f, dstype *dP, dstype *dtmp, dstype *df, dstype *fac, dstype pi, Int L, Int K, Int N, Int backend)
// {
//     if (backend == 1)  
//         cpuSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, x, y, z, 
//             x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //         gpuSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, x, y, z, 
// //             x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);    
//     }
// #endif                  
// }
// 
// inline void coreSphericalHarmonicsBesselWithDeriv(dstype *Sr, dstype *Si, dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, 
//       dstype *xij, dstype *x0, dstype *P, dstype *tmp, dstype *f, dstype *dP, dstype *dtmp, dstype *df, dstype *fac, dstype pi, Int L, Int K, Int N, Int backend)
// {
//     if (backend == 1)  
//         cpuSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, xij,
//             x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //         gpuSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, xij,
// //             x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);    
//     }
// #endif                  
// }
// 
// inline void coreRadialSphericalHarmonicsSum(dstype *ar, dstype *ai, dstype *Sr, dstype *Si, Int *Nnb, Int Na, Int L, Int K, Int backend)
// {
//     if (backend == 1)  
//         cpuRadialSphericalHarmonicsSum(ar, ai, Sr, Si, Nnb, Na, L, K);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
//         //gpuRadialSphericalHarmonicsSum(ar, ai, Sr, Si, Nnb, Na, L, K);    
//     }
// #endif                  
// }
// 
// inline void coreRadialSphericalHarmonicsPower(dstype *p, dstype *ar, dstype *ai, Int *indk, Int Na, Int L, Int K, Int backend)
// {
//     if (backend == 1)  
//         cpuRadialSphericalHarmonicsPower(p, ar, ai, indk, Na, L, K);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //        gpuRadialSphericalHarmonicsPower(p, ar, ai, indk, Na, L, K);    
//     }
// #endif                  
// }
// 
// inline void coreRadialSphericalHarmonicsPowerDeriv(dstype *px, dstype *py, dstype *pz, dstype *ar, dstype *ai, 
//         dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *indk, Int *Nnb, Int Na, Int L, Int K, Int backend)
// {
//     if (backend == 1)  
//         cpuRadialSphericalHarmonicsPowerDeriv(px, py, pz, ar, ai, arx, aix, 
//             ary, aiy, arz, aiz, indk, Nnb, Na, L, K);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //         gpuRadialSphericalHarmonicsPowerDeriv(px, py, pz, ar, ai, arx, aix, 
// //             ary, aiy, arz, aiz, indk, Nnb, Na, L, K);
//     }
// #endif                  
// }
// 
// inline void coreRadialSphericalHarmonicsPowerDeriv2(dstype *pd, dstype *ar, dstype *ai, dstype *arx, dstype *aix,
//         dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *indk, Int *Nnb, Int Na, Int L, Int K, Int backend)
// {
//     if (backend == 1)  
//         cpuRadialSphericalHarmonicsPowerDeriv2(pd, ar, ai, arx, aix, 
//             ary, aiy, arz, aiz, indk, Nnb, Na, L, K);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //         gpuRadialSphericalHarmonicsPowerDeriv2(pd, ar, ai, arx, aix, 
// //             ary, aiy, arz, aiz, indk, Nnb, Na, L, K);
//     }
// #endif                  
// }
// 
// inline void coreRadialSphericalHarmonicsBispectrum(dstype *b, dstype *ar, dstype *ai, dstype *cg, Int *indk, 
//         Int *indl, Int *indm, Int *rowm, Int Nub, Int Ncg, Int Na, Int L, Int K, Int backend)
// {
//     if (backend == 1)  
//         cpuRadialSphericalHarmonicsBispectrum(b, ar, ai, cg, indk, indl, indm, rowm,
//             Nub, Ncg, Na, L, K);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //         gpuRadialSphericalHarmonicsBispectrum(b, ar, ai, cg, indk, indl, indm, rowm,
// //             Nub, Ncg, Na, L, K);
//     }
// #endif                  
// }
// 
// inline void coreRadialSphericalHarmonicsBispectrumDeriv(dstype *bx, dstype *by, dstype *bz, 
//         dstype *ar, dstype *ai, dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, dstype*cg, Int *indk, Int *indl,
//         Int *indm, Int *rowm, Int *Nnb, Int Na, Int Nub, Int Ncg, Int K, Int backend)
// {
//     if (backend == 1)  
//         cpuRadialSphericalHarmonicsBispectrumDeriv(bx, by, bz, ar, ai, arx, aix, 
//             ary, aiy, arz, aiz, cg, indk, indl, indm, rowm, Nnb, Na, Nub, Ncg, K);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //         gpuRadialSphericalHarmonicsBispectrumDeriv(bx, by, bz, ar, ai, arx, aix, 
// //             ary, aiy, arz, aiz, cg, indk, indl, indm, rowm, Nnb, Na, Nub, Ncg, K);
//     }
// #endif                  
// }
// 
// inline void coreRadialSphericalHarmonicsBispectrumDeriv2(dstype *bd, dstype *ar, dstype *ai, 
//         dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, dstype*cg, 
//         Int *indk, Int *indl, Int *indm, Int *rowm, Int *Nnb, Int Na, Int Nub, Int Ncg, Int K, Int backend)
// {
//     if (backend == 1)  
//         cpuRadialSphericalHarmonicsBispectrumDeriv2(bd, ar, ai, arx, aix, 
//             ary, aiy, arz, aiz, cg, indk, indl, indm, rowm, Nnb, Na, Nub, Ncg, K);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //         gpuRadialSphericalHarmonicsBispectrumDeriv2(bd, ar, ai, arx, aix, 
// //             ary, aiy, arz, aiz, cg, indk, indl, indm, rowm, Nnb, Na, Nub, Ncg, K);
//     }
// #endif                  
// }
// 
// inline void coreRadialSphericalHarmonicsBasis(dstype *d, dstype *c, Int *atomtype, 
//         Int Ntype, Int Na, Int Nbf, Int backend)
// {
//     if (backend == 1)  
//         cpuRadialSphericalHarmonicsBasis(d, c, atomtype, Ntype, Na, Nbf);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
//         //gpuRadialSphericalHarmonicsBasis(d, c, atomtype, Ntype, Na, Nbf);
//     }
// #endif                  
// }
// 
// inline void coreRadialSphericalHarmonicsBasisDeriv(dstype *dx, dstype *dy, dstype *dz, 
//         dstype *cx, dstype *cy, dstype *cz, Int *atomtype, Int *neighlist, Int *Nnb, Int Ntype, Int Na, Int Nbf, Int backend)
// {
//     if (backend == 1)  
//         cpuRadialSphericalHarmonicsBasisDeriv(dx, dy, dz, cx, cy, cz, atomtype, 
//             neighlist, Nnb, Ntype, Na, Nbf);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
// //         gpuRadialSphericalHarmonicsBasisDeriv(dx, dy, dz, cx, cy, cz, atomtype, 
// //             neighlist, Nnb, Ntype, Na, Nbf);
//     }
// #endif                  
// }
// 
// inline void coreRadialSphericalHarmonicsBasisDeriv2(dstype *dd, dstype *cd, 
//         Int *atomtype, Int *neighlist, Int *Nnb, Int Ntype, Int Na, Int Nbf, Int backend)
// {
//     if (backend == 1)  
//         cpuRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, neighlist, Nnb, Ntype, Na, Nbf);    
// #ifdef HAVE_CUDA            
//     if (backend == 2)  {
//        // gpuRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, neighlist, Nnb, Ntype, Na, Nbf);
//     }
// #endif                  
// }

#endif  

