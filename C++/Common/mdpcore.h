/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_CORE_H
#define MDP_CORE_H

#include "mdpcpucore.h"

#ifdef USE_CUDA      
#include "mdpgpucore.h"
#endif                

#ifdef USE_HIP      
#include "mdphipcore.h"
#endif                

#ifdef USE_OMP      
#include "mdpompcore.h"
#endif                

inline void Kron(dstype *C, dstype *A, dstype *B, int M1, int M2, int backend)
{
    if (backend == 1)
        cpuKron(C, A, B, M1, M2);
#ifdef USE_OMP
    if (backend == 4)
        ompKron(C, A, B, M1, M2);
#endif                          
#ifdef USE_HIP
    if (backend == 3)
        hipKron(C, A, B, M1, M2);
#endif                                  
#ifdef USE_CUDA            
    if (backend == 2)
        gpuKron(C, A, B, M1, M2);
#endif                          
}

inline void GetArrayAtIndex(dstype *y, dstype *x, int *ind, int n, int backend)
{
    if (backend == 1)
        cpuGetArrayAtIndex(y, x, ind, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompGetArrayAtIndex(y, x, ind, n);
#endif                            
#ifdef USE_HIP  
    if (backend == 3)
        hipGetArrayAtIndex(y, x, ind, n);
#endif                                    
#ifdef USE_CUDA   
    if (backend == 2)
        gpuGetArrayAtIndex(y, x, ind, n);
#endif                  
}

inline void PutArrayAtIndex(dstype *y, dstype *x, int *ind, int n, int backend)
{
    if (backend == 1)
        cpuPutArrayAtIndex(y, x, ind, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompPutArrayAtIndex(y, x, ind, n);
#endif        
#ifdef USE_HIP  
    if (backend == 3)
        hipPutArrayAtIndex(y, x, ind, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuPutArrayAtIndex(y, x, ind, n);
#endif                  
}

inline void ArrayPlusXAtIndex(dstype *y, dstype *x, int *ind, int n, int backend)
{
    if (backend == 1)
        cpuArrayPlusXAtIndex(y, x, ind, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayPlusXAtIndex(y, x, ind, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayPlusXAtIndex(y, x, ind, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayPlusXAtIndex(y, x, ind, n);
#endif                  
}

inline void ArrayMinusXAtIndex(dstype *y, dstype *x, int *ind, int n, int backend)
{
    if (backend == 1)
        cpuArrayMinusXAtIndex(y, x, ind, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayMinusXAtIndex(y, x, ind, n);
#endif        
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayMinusXAtIndex(y, x, ind, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayMinusXAtIndex(y, x, ind, n);
#endif                  
}

inline void ArraySquare(dstype *y, dstype *x, int n, int backend)
{
    if (backend == 1)
        cpuArraySquare(y, x, n);   
#ifdef USE_OMP  
    if (backend == 4)
        ompArraySquare(y, x, n);   
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipArraySquare(y, x, n);   
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArraySquare(y, x, n);
#endif                  
}

inline void ArraySetValue(dstype *y, dstype a, int n, int backend)
{
    if (backend == 1)
        cpuArraySetValue(y, a, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArraySetValue(y, a, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipArraySetValue(y, a, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArraySetValue(y, a, n);
#endif                  
}

inline void ArraySetValueAtIndex(dstype *y, dstype a, int n, int backend)
{
    if (backend == 1)
        y[n] = a;
#ifdef USE_OMP     
    if (backend == 4)
        y[n] = a;
#endif                 
#ifdef USE_HIP                
    if (backend == 3)
        hipArraySetValueAtIndex(y, a, n);
#endif                         
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArraySetValueAtIndex(y, a, n);
#endif                  
}

inline dstype ArrayGetValueAtIndex(dstype *y, int n, int backend)
{
    dstype val; 
    if (backend == 1)
        val = y[n];
#ifdef USE_OMP  
    if (backend == 4)
        val = y[n];
#endif                             
#ifdef USE_HIP    
    if (backend == 3)
        val = hipArrayGetValueAtIndex(y, n);
#endif                              
#ifdef USE_CUDA   
    if (backend == 2)
        val = gpuArrayGetValueAtIndex(y, n);
#endif                      
    return val;
}

inline void ArrayTranspose(dstype *A, dstype *B, int m, int n, int backend)
{
	if (backend == 1)
		cpuArrayTranspose(A, B, m, n);
#ifdef USE_OMP
	if (backend == 4)
		ompArrayTranspose(A, B, m, n);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipArrayTranspose(A, B, m, n);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuArrayTranspose(A, B, m, n);
#endif
}

inline void GetArrayAtColumnIndex(dstype *A, dstype *B, int *colind, int m, int n, int backend)
{
	if (backend == 1)
		cpuGetArrayAtColumnIndex(A, B, colind, m, n);
#ifdef USE_OMP
	if (backend == 4)
		ompGetArrayAtColumnIndex(A, B, colind, m, n);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipGetArrayAtColumnIndex(A, B, colind, m, n);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuGetArrayAtColumnIndex(A, B, colind, m, n);
#endif
}

inline void ArrayTransposeAtColumnIndex(dstype *A, dstype *B, int *colind, int m, int n, int backend)
{
	if (backend == 1)
		cpuArrayTransposeAtColumnIndex(A, B, colind, m, n);
#ifdef USE_OMP
	if (backend == 4)
		ompArrayTransposeAtColumnIndex(A, B, colind, m, n);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipArrayTransposeAtColumnIndex(A, B, colind, m, n);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuArrayTransposeAtColumnIndex(A, B, colind, m, n);
#endif
}

inline void GetArrayAtRowIndex(dstype *A, dstype *B, int *rowind, int m, int n, int backend)
{
	if (backend == 1)
		cpuGetArrayAtRowIndex(A, B, rowind, m, n);
#ifdef USE_OMP
	if (backend == 4)
		ompGetArrayAtRowIndex(A, B, rowind, m, n);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipGetArrayAtRowIndex(A, B, rowind, m, n);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuGetArrayAtRowIndex(A, B, rowind, m, n);
#endif
}

inline void ArrayTransposeAtRowIndex(dstype *A, dstype *B, int *rowind, int m, int n, int backend)
{
	if (backend == 1)
		cpuArrayTransposeAtRowIndex(A, B, rowind, m, n);
#ifdef USE_OMP
	if (backend == 4)
		ompArrayTransposeAtRowIndex(A, B, rowind, m, n);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipArrayTransposeAtRowIndex(A, B, rowind, m, n);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuArrayTransposeAtRowIndex(A, B, rowind, m, n);
#endif
}

inline void ArrayRowSum(dstype *y, dstype *x, int m, int n, int backend)
{
	if (backend == 1)
		cpuArrayRowSum(y, x, m, n);
#ifdef USE_OMP
	if (backend == 4)
		ompArrayRowSum(y, x, m, n);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipArrayRowSum(y, x, m, n);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuArrayRowSum(y, x, m, n);
#endif
}

inline void ArrayRowSquareSum(dstype *y, dstype *x, int m, int n, int backend)
{
	if (backend == 1)
		cpuArrayRowSquareSum(y, x, m, n);
#ifdef USE_OMP
	if (backend == 4)
		ompArrayRowSquareSum(y, x, m, n);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipArrayRowSquareSum(y, x, m, n);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuArrayRowSquareSum(y, x, m, n);
#endif
}

inline void ArrayDistSquareSum(dstype *y, dstype *x1, dstype *x2, int m, int n, int backend)
{
	if (backend == 1)
		cpuArrayDistSquareSum(y, x1, x2, m, n);
#ifdef USE_OMP
	if (backend == 4)
		ompArrayDistSquareSum(y, x1, x2, m, n);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipArrayDistSquareSum(y, x1, x2, m, n);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuArrayDistSquareSum(y, x1, x2, m, n);
#endif
}

inline dstype ArraySum(dstype *a, dstype *b, int n, int backend)
{
    dstype val; 
    if (backend == 1)
        val = cpuArraySum(a, n);
#ifdef USE_OMP  
    if (backend == 4)
        val = ompArraySum(a, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        val = hipArraySum(a, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        val = gpuArraySum(a, b, n);
#endif                      
    return val;
}

inline dstype ArrayMax(dstype *a, dstype *b, int n, int backend)
{
    dstype val; 
    if (backend == 1)
        val = cpuArrayMax(a, n);
#ifdef USE_OMP  
    if (backend == 4)
        val = ompArrayMax(a, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        val = hipArrayMax(a, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        val = gpuArrayMax(a, b, n);
#endif                      
    return val;
}


inline dstype ArrayMin(dstype *a, dstype *b, int n, int backend)
{
    dstype val; 
    if (backend == 1)
        val = cpuArrayMin(a, n);
#ifdef USE_OMP  
    if (backend == 4)
        val = ompArrayMin(a, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        val = hipArrayMin(a, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        val = gpuArrayMin(a, b, n);
#endif                      
    return val;
}

inline void ArraySumEveryColumn(dstype *b, dstype *a, int ncol, int n, int backend)
{
    if (backend == 1)
        cpuArraySumEveryColumn(b, a, ncol, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArraySumEveryColumn(b, a, ncol, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipArraySumEveryColumn(b, a, ncol, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArraySumEveryColumn(b, a, ncol, n);
#endif                          
}

inline void ArrayMaxEveryColumn(dstype *b, dstype *a, int ncol, int n, int backend)
{
    if (backend == 1)
        cpuArrayMaxEveryColumn(b, a, ncol, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayMaxEveryColumn(b, a, ncol, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayMaxEveryColumn(b, a, ncol, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayMaxEveryColumn(b, a, ncol, n);
#endif                          
}

inline void ArrayMinEveryColumn(dstype *b, dstype *a, int ncol, int n, int backend)
{
    if (backend == 1)
        cpuArrayMinEveryColumn(b, a, ncol, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayMinEveryColumn(b, a, ncol, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayMinEveryColumn(b, a, ncol, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayMinEveryColumn(b, a, ncol, n);
#endif                          
}

inline void ArraySumEveryRow(dstype *b, dstype *a, dstype *c, int nrow, int n, int backend)
{
    if (backend == 1)
        cpuArraySumEveryRow(b, a, nrow, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArraySumEveryRow(b, a, nrow, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipArraySumEveryRow(b, a, c, nrow, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArraySumEveryRow(b, a, c, nrow, n);
#endif                          
}

inline void ArrayMaxEveryRow(dstype *b, dstype *a, dstype *c, int nrow, int n, int backend)
{
    if (backend == 1)
        cpuArrayMaxEveryRow(b, a, nrow, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayMaxEveryRow(b, a, nrow, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayMaxEveryRow(b, a, c, nrow, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayMaxEveryRow(b, a, c, nrow, n);
#endif                          
}

inline void ArrayMinEveryRow(dstype *b, dstype *a, dstype *c, int nrow, int n, int backend)
{
    if (backend == 1)
        cpuArrayMinEveryRow(b, a, nrow, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayMinEveryRow(b, a, nrow, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayMinEveryRow(b, a, c, nrow, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayMinEveryRow(b, a, c, nrow, n);
#endif                          
}

inline int IntArrayGetValueAtIndex(int *y, int n, int backend)
{
    int val; 
    if (backend == 1)
        val = y[n];
#ifdef USE_OMP  
    if (backend == 4)
        val = y[n];
#endif                                     
#ifdef USE_HIP    
    if (backend == 3)
        val = hipArrayGetValueAtIndex(y, n);
#endif                              
#ifdef USE_CUDA   
    if (backend == 2)
        val = gpuArrayGetValueAtIndex(y, n);
#endif                      
    return val;
}

inline void ArrayFill(int *a, int start, int n, int backend)
{
    if (backend == 1)
        cpuArrayFill(a, start, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayFill(a, start, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayFill(a, start, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayFill(a, start, n);
#endif                  
}

inline void ArrayTripletnum(int* output, int* input, int length, int backend) 
{
    if (backend == 1)
        cpuTripletnum(output, input, length); 
#ifdef USE_OMP  
    if (backend == 4)
        ompTripletnum(output, input, length); 
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipTripletnum(output, input, length); 
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuTripletnum(output, input, length); 
#endif                  
}
inline void ArrayQuadrupletnum(int* output, int* input, int length, int backend) 
{
    if (backend == 1)
        cpuQuadrupletnum(output, input, length); 
#ifdef USE_OMP  
    if (backend == 4)
        ompQuadrupletnum(output, input, length); 
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipQuadrupletnum(output, input, length); 
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuQuadrupletnum(output, input, length); 
#endif                  
}

inline void ArrayAddScalar(dstype *y, dstype a, int n, int backend)
{
    if (backend == 1)
        cpuArrayAddScalar(y, a, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayAddScalar(y, a, n);
#endif                             
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayAddScalar(y, a, n);
#endif                                     
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayAddScalar(y, a, n);
#endif                  
}

inline void ArrayMultiplyScalar(dstype *y, dstype a, int n, int backend)
{
    if (backend == 1)
        cpuArrayMultiplyScalar(y, a, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayMultiplyScalar(y, a, n);
#endif                                  
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayMultiplyScalar(y, a, n);
#endif                                          
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayMultiplyScalar(y, a, n);
#endif                  
}

inline void ArrayCopy(dstype *y, dstype *x, int n, int backend)
{
    if (backend == 1)
        cpuArrayCopy(y, x, n);   
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayCopy(y, x, n);   
#endif                        
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayCopy(y, x, n);   
#endif                                
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayCopy(y, x, n);
#endif                  
}

inline void ArrayCopy(int *y, int *x, int n, int backend)
{
    if (backend == 1)
        cpuArrayCopy(y, x, n);   
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayCopy(y, x, n);   
#endif                        
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayCopy(y, x, n);   
#endif                                
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayCopy(y, x, n);
#endif                  
}

inline void ArrayMinus(dstype *y, dstype *x, int n, int backend)
{
    if (backend == 1)
        cpuArrayMinus(y, x, n);  
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayMinus(y, x, n);  
#endif                        
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayMinus(y, x, n);  
#endif                                
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayMinus(y, x, n);
#endif                  
}

inline void ArrayAbs(dstype *y, dstype *x, int n, int backend)
{
    if (backend == 1)
        cpuArrayAbs(y, x, n);   
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayAbs(y, x, n);   
#endif                        
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayAbs(y, x, n);   
#endif                                
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayAbs(y, x, n);
#endif                  
}

inline void ArraySqrt(dstype *y, dstype *x, int n, int backend)
{
    if (backend == 1)
        cpuArraySqrt(y, x, n);   
#ifdef USE_OMP  
    if (backend == 4)
        ompArraySqrt(y, x, n);   
#endif                        
#ifdef USE_HIP  
    if (backend == 3)
        hipArraySqrt(y, x, n);   
#endif                                
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArraySqrt(y, x, n);
#endif                  
}

inline void ArrayMultiplyScalarDiagonal(dstype *C, dstype a, int n, int backend)
{
    if (backend == 1)
        cpuArrayMultiplyScalarDiagonal(C, a, n);   
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayMultiplyScalarDiagonal(C, a, n);   
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayMultiplyScalarDiagonal(C, a, n);   
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayMultiplyScalarDiagonal(C, a, n);
#endif                  
}

inline void ArrayAddVectorToDiagonal(dstype *C, dstype *x, dstype a, int n, int backend)
{
    if (backend == 1)
        cpuArrayAddVectorToDiagonal(C, x, a, n);   
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayAddVectorToDiagonal(C, x, a, n);   
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayAddVectorToDiagonal(C, x, a, n);   
#endif                                        
#ifdef USE_CUDA                
    if (backend == 2)
        gpuArrayAddVectorToDiagonal(C, x, a, n);
#endif                  
}

inline void ArrayAXPB(dstype *y, dstype *x, dstype a, dstype b, int n, int backend)
{
    if (backend == 1)
        cpuArrayAXPB(y, x, a, b, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayAXPB(y, x, a, b, n);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayAXPB(y, x, a, b, n);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayAXPB(y, x, a, b, n);
#endif                  
}

inline void ArrayAXPBY(dstype *z, dstype *x, dstype *y, dstype a, dstype b, int n, int backend)
{
    if (backend == 1)
        cpuArrayAXPBY(z, x, y, a, b, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayAXPBY(z, x, y, a, b, n);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayAXPBY(z, x, y, a, b, n);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayAXPBY(z, x, y, a, b, n);
#endif                  
}

inline void ArrayAXY(dstype *z, dstype *x, dstype *y, dstype a, int n, int backend)
{
    if (backend == 1)
        cpuArrayAXY(z, x, y, a, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayAXY(z, x, y, a, n);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayAXY(z, x, y, a, n);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayAXY(z, x, y, a, n);
#endif                  
}

inline void ArrayAXYZ(dstype *s, dstype *x, dstype *y, 
        dstype *z, dstype a, int n, int backend)
{
    if (backend == 1)
        cpuArrayAXYZ(s, x, y, z, a, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayAXYZ(s, x, y, z, a, n);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayAXYZ(s, x, y, z, a, n);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayAXYZ(s, x, y, z, a, n);
#endif                  
}

inline void ArrayAXYPBZ(dstype *s, dstype *x, dstype *y, dstype *z, dstype a,
        dstype b, int n, int backend)
{
    if (backend == 1)
        cpuArrayAXYPBZ(s, x, y, z, a, b, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayAXYPBZ(s, x, y, z, a, b, n);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayAXYPBZ(s, x, y, z, a, b, n);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayAXYPBZ(s, x, y, z, a, b, n);
#endif                  
}

inline void ArrayAdd3Vectors(dstype *s, dstype *x, dstype *y, dstype *z, 
        dstype a, dstype b, dstype c, int n, int backend)
{
    if (backend == 1)
        cpuArrayAdd3Vectors(s, x, y, z, a, b, c, n);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayAdd3Vectors(s, x, y, z, a, b, c, n);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayAdd3Vectors(s, x, y, z, a, b, c, n);
#endif            
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayAdd3Vectors(s, x, y, z, a, b, c, n);
#endif                  
}

inline void ArrayExtract(dstype *un, dstype *u, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2, int backend)
{
    if (backend == 1)
        cpuArrayExtract(un, u, I, J, K, i1, i2, j1, j2, k1, k2);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayExtract(un, u, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayExtract(un, u, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayExtract(un, u, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                  
}

inline void ArrayInsert(dstype *u, dstype *un, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2, int backend)
{
    if (backend == 1)
        cpuArrayInsert(u, un, I, J, K, i1, i2, j1, j2, k1, k2);
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayInsert(u, un, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayInsert(u, un, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayInsert(u, un, I, J, K, i1, i2, j1, j2, k1, k2);    
#endif                  
}

inline void ArrayGemmBatch(dstype *C, dstype *A, dstype *B, 
        int I, int J, int K, int S, int backend)
{
    if (backend == 1)
        cpuArrayGemmBatch(C, A, B, I, J, K, S);    
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayGemmBatch(C, A, B, I, J, K, S);    
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayGemmBatch(C, A, B, I, J, K, S);    
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayGemmBatch(C, A, B, I, J, K, S);    
#endif                  
}

inline void ArrayGemmBatch1(dstype *C, dstype *A, dstype *B, 
        int I, int J, int K, int S, int backend)
{
    if (backend == 1)
        cpuArrayGemmBatch1(C, A, B, I, J, K, S);    
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayGemmBatch1(C, A, B, I, J, K, S);    
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayGemmBatch1(C, A, B, I, J, K, S);    
#endif                                       
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayGemmBatch1(C, A, B, I, J, K, S);    
#endif                  
}

inline void ArrayDG2CG(dstype *ucg, dstype *udg, int *cgent2dgent, int *rowent2elem, int nent, int backend)
{
    if (backend == 1)
        cpuArrayDG2CG(ucg, udg, cgent2dgent, rowent2elem, nent);    
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayDG2CG(ucg, udg, cgent2dgent, rowent2elem, nent);    
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayDG2CG(ucg, udg, cgent2dgent, rowent2elem, nent);    
#endif                                       
#ifdef USE_CUDA   
    if (backend == 2)
        gpuArrayDG2CG(ucg, udg, cgent2dgent, rowent2elem, nent);    
#endif                      
}

inline void ArrayDG2CG2(dstype *ucg, dstype *udg, int *colent2elem, int *rowent2elem, int nent, int npe, int backend)
{
    if (backend == 1)
        cpuArrayDG2CG2(ucg, udg, colent2elem, rowent2elem, nent, npe);    
#ifdef USE_OMP  
    if (backend == 4)
        ompArrayDG2CG2(ucg, udg, colent2elem, rowent2elem, nent, npe);    
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipArrayDG2CG2(ucg, udg, colent2elem, rowent2elem, nent, npe);    
#endif                                       
#ifdef USE_CUDA                
    if (backend == 2)
        gpuArrayDG2CG2(ucg, udg, colent2elem, rowent2elem, nent, npe);    
#endif                  
}

inline void Cart2Sphere(dstype *the, dstype *phi, dstype *r, dstype *x, 
        dstype *y, dstype *z, int N, int backend)
{
    if (backend == 1)
        cpuCart2Sphere(the, phi, r, x, y, z, N);    
#ifdef USE_OMP  
    if (backend == 4)
        ompCart2Sphere(the, phi, r, x, y, z, N);    
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipCart2Sphere(the, phi, r, x, y, z, N);    
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuCart2Sphere(the, phi, r, x, y, z, N);    
#endif                  
}

inline void Cart2SphereDeriv(dstype *the, dstype *phi, dstype *r, dstype *thex, dstype *they, dstype *thez, dstype *phix, 
        dstype *phiy, dstype *phiz, dstype *rx, dstype *ry, dstype *rz, dstype *x, dstype *y, dstype *z, int N, int backend)
{
    if (backend == 1)
        cpuCart2SphereDeriv(the, phi, r, thex, they, thez, phix, phiy, phiz,
            rx, ry, rz, x, y, z, N);
#ifdef USE_OMP  
    if (backend == 4)
        ompCart2SphereDeriv(the, phi, r, thex, they, thez, phix, phiy, phiz,
            rx, ry, rz, x, y, z, N);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipCart2SphereDeriv(the, phi, r, thex, they, thez, phix, phiy, phiz,
            rx, ry, rz, x, y, z, N);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuCart2SphereDeriv(the, phi, r, thex, they, thez, phix, phiy, phiz,
            rx, ry, rz, x, y, z, N);
#endif                  
}

inline void Sphere2Cart(dstype *x, dstype *y, dstype *z, dstype *the, 
        dstype *phi, dstype *r, int N, int backend)
{
    if (backend == 1)
        cpuSphere2Cart(x, y, z, the, phi, r, N);    
#ifdef USE_OMP  
    if (backend == 4)
        ompSphere2Cart(x, y, z, the, phi, r, N);    
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipSphere2Cart(x, y, z, the, phi, r, N);    
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuSphere2Cart(x, y, z, the, phi, r, N);    
#endif                  
}

inline void Euler2Rotm(dstype *R11, dstype *R12, dstype *R13, dstype *R21, 
                dstype *R22, dstype *R23, dstype *R31, dstype *R32, dstype *R33,
        dstype *alpha, dstype *beta, dstype *gamma, int N, int backend)
{
    if (backend == 1)
        cpuEuler2Rotm(R11, R12, R13, R21, R22, R23, R31, R32, R33, alpha, beta, gamma, N);    
#ifdef USE_OMP  
    if (backend == 4)
        ompEuler2Rotm(R11, R12, R13, R21, R22, R23, R31, R32, R33, alpha, beta, gamma, N);    
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipEuler2Rotm(R11, R12, R13, R21, R22, R23, R31, R32, R33, alpha, beta, gamma, N);    
#endif                                
#ifdef USE_CUDA   
    if (backend == 2)
        gpuEuler2Rotm(R11, R12, R13, R21, R22, R23, R31, R32, R33, alpha, beta, gamma, N);    
#endif                  
}

inline void Rotc(dstype *X, dstype *Y, dstype *Z, dstype *R, dstype *x, 
        dstype *y, dstype *z, int N, int backend)
{
    if (backend == 1)
        cpuRotc(X, Y, Z, R, x, y, z, N);    
#ifdef USE_OMP  
    if (backend == 4)
        ompRotc(X, Y, Z, R, x, y, z, N);    
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipRotc(X, Y, Z, R, x, y, z, N);    
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuRotc(X, Y, Z, R, x, y, z, N);    
#endif                  
}

inline void Cumsum(int *d_out, int *d_in, int *d_sums, int *d_incr, int length, int backend)
{
    if (backend == 1)
        cpuCumsum(d_out, d_in, length);
#ifdef USE_OMP  
    if (backend == 4)
        ompCumsum(d_out, d_in, length);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipCumsum(d_out, d_in, d_sums, d_incr, length);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuCumsum(d_out, d_in, d_sums, d_incr, length);
#endif                      
}


inline int FindAtomType(int *ilist, int* olist, int *atomtype, int *t0, int *t1, int typei, int na, int backend)
{
    if (backend == 1) { 
        return cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
    }
#ifdef USE_OMP            
    if (backend == 4) { 
        return ompFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
    }
#endif                                  
#ifdef USE_HIP            
    if (backend == 3) { 
        return hipFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
    }
#endif                              
#ifdef USE_CUDA            
    if (backend == 2) { 
        return gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
    }
#endif                          
    return 0;
}                

inline int UniqueSort(int *b, int *c, int *d, int *e, int *a, int *p, int *t, int *q, int n, int backend)
{
    if (backend == 1)  
        return cpuUniqueSort(b, c, d, e, a, p, t, q, n);
#ifdef USE_OMP            
    if (backend == 4) { 
        return ompUniqueSort(b, c, d, e, a, p, t, q, n);
    }
#endif                                  
#ifdef USE_HIP            
    if (backend == 3) { 
        return hipUniqueSort(b, c, d, e, a, p, t, q, n);
    }
#endif                                 
#ifdef USE_CUDA            
    if (backend == 2)  {
        return gpuUniqueSort(b, c, d, e, a, p, t, q, n);
    }
#endif             
    return 0;
}

inline void AtomList2D(int *alist,  int *inside, int *glistnumsum, int *glistnum, int *d_sums, int *d_incr, 
        dstype *x, dstype *pimages, dstype *wc, dstype *s2rmap, int inum, int pnum, int dim, Int backend)
{
    if (backend == 1) 
        cpuAtomList2D(alist, inside, glistnumsum, glistnum, x, pimages, wc, 
                s2rmap, inum, pnum, dim);   
#ifdef USE_OMP             
    if (backend == 4) 
        ompAtomList2D(alist, inside, glistnumsum, glistnum, x, pimages, wc, 
                s2rmap, inum, pnum, dim);              
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipAtomList2D(alist, inside, glistnumsum, glistnum, d_sums, d_incr, x, pimages, 
                wc, s2rmap, inum, pnum, dim);    
#endif                              
#ifdef USE_CUDA            
    if (backend == 2) 
        gpuAtomList2D(alist, inside, glistnumsum, glistnum, d_sums, d_incr, x, pimages, 
                wc, s2rmap, inum, pnum, dim);    
#endif                          
}        

inline void AtomList3D(int *alist,  int *inside, int *glistnumsum, int *glistnum, int *d_sums, int *d_incr, 
        dstype *x, dstype *pimages, dstype *wc, dstype *s2rmap, int inum, int pnum, int dim, int backend)
{
    if (backend == 1)  
        cpuAtomList3D(alist, inside, glistnumsum, glistnum, x, pimages, wc, 
                s2rmap, inum, pnum, dim);     
#ifdef USE_OMP             
    if (backend == 4) 
        ompAtomList3D(alist, inside, glistnumsum, glistnum, x, pimages, wc, 
                s2rmap, inum, pnum, dim);     
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipAtomList3D(alist, inside, glistnumsum, glistnum, d_sums, d_incr, x, pimages, 
                wc, s2rmap, inum, pnum, dim);    
#endif                                  
#ifdef USE_CUDA            
    if (backend == 2) 
        gpuAtomList3D(alist, inside, glistnumsum, glistnum, d_sums, d_incr, x, pimages, 
                wc, s2rmap, inum, pnum, dim);    
#endif                          
}        

inline void CellList2D(int *clist, dstype *x, dstype *eta1, dstype *eta2, dstype *eta3, dstype *s2rmap, 
        int *nc, int inum, int natom, int dim, int backend)
{
    if (backend == 1) 
        cpuCellList2D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);  
#ifdef USE_OMP             
    if (backend == 4) 
        ompCellList2D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);  
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipCellList2D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);  
#endif                                      
#ifdef USE_CUDA            
    if (backend == 2) 
        gpuCellList2D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);      
#endif                          
}                

inline void CellList3D(int *clist, dstype *x, dstype *eta1, dstype *eta2, dstype *eta3, dstype *s2rmap, 
        int *nc, int inum, int natom, int dim, int backend)
{
    if (backend == 1) 
        cpuCellList3D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);  
#ifdef USE_OMP             
    if (backend == 4) 
        ompCellList3D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);  
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipCellList3D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);  
#endif                                          
#ifdef USE_CUDA            
    if (backend == 2) 
        gpuCellList3D(clist, x, eta1, eta2, eta3, s2rmap, nc, inum, natom, dim);      
#endif                          
}                

inline void Cell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int *dtemp, 
        int natom, int ncell, int backend)
{
    if (backend == 1) 
        cpuCell2AtomList(c2alist, c2anumsum, c2anum, clist, natom, ncell);        
#ifdef USE_OMP             
    if (backend == 4) 
        ompCell2AtomList(c2alist, c2anumsum, c2anum, clist, natom, ncell);        
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipCell2AtomList(c2alist, c2anumsum, c2anum, clist, dtemp, natom, ncell);        
#endif                                              
#ifdef USE_CUDA            
    if (backend == 2) 
        gpuCell2AtomList(c2alist, c2anumsum, c2anum, clist, dtemp, natom, ncell);        
#endif                          
}                

inline void FullNeighborList2D(int *neighlist, int *neighnum, dstype *x, dstype *rcutsq, 
        int anum, int inum, int jnum, int dim, int backend)
{
    if (backend == 1) 
        cpuFullNeighborList2D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);    
#ifdef USE_OMP             
    if (backend == 4) 
        ompFullNeighborList2D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);    
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipFullNeighborList2D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);        
#endif                                                  
#ifdef USE_CUDA            
    if (backend == 2) { 
        gpuFullNeighborList2D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);        
    }
#endif                          
}                

inline void FullNeighborList3D(int *neighlist, int *neighnum, dstype *x, dstype *rcutsq, 
        int anum, int inum, int jnum, int dim, int backend)
{
    if (backend == 1) 
        cpuFullNeighborList3D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);    
#ifdef USE_OMP             
    if (backend == 4) 
        ompFullNeighborList3D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);    
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipFullNeighborList3D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);        
#endif                                                      
#ifdef USE_CUDA            
    if (backend == 2) 
        gpuFullNeighborList3D(neighlist, neighnum, x, rcutsq, anum, inum, jnum, dim);            
#endif                          
}                

inline void FullNeighborList2D(int *neighlist, int *neighnum, dstype *x, dstype *rcutsq,
        int *alist, int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim, int backend)
{
    if (backend == 1) 
        cpuFullNeighborList2D(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);     
#ifdef USE_OMP             
    if (backend == 4) 
        ompFullNeighborList2D(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);     
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipFullNeighborList2D(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);         
#endif                                                          
#ifdef USE_CUDA            
    if (backend == 2) 
        gpuFullNeighborList2D(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);         
#endif                          
}                

inline void FullNeighborList3D(int *neighlist, int *neighnum, dstype *x, dstype *rcutsq,
        int *alist, int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim, int backend)
{
    if (backend == 1) 
        cpuFullNeighborList3D(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);         
#ifdef USE_OMP             
    if (backend == 4) 
        ompFullNeighborList3D(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);     
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipFullNeighborList23(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);         
#endif                                                              
#ifdef USE_CUDA            
    if (backend == 2) 
        gpuFullNeighborList3D(neighlist, neighnum, x, rcutsq, alist, clist,   
                        c2alist, c2anumsum, nc, inum, jnum, dim);         
#endif                          
}                

inline void FullNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim, Int backend)
{
    if (backend == 1)  
        cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, neighnum, inum, jnum, dim);        
#ifdef USE_OMP             
    if (backend == 4) 
        ompFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, neighnum, inum, jnum, dim);        
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, neighnum, inum, jnum, dim);        
#endif                                                                  
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, neighnum, inum, jnum, dim);        
#endif                  
}

inline void FullNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *atomtype, int *ilist, 
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim, Int backend)
{
    if (backend == 1)  
        cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#ifdef USE_OMP             
    if (backend == 4) 
        ompFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#endif                                                                      
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#endif                  
}

inline void FullNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *atomtype, int *ilist, 
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int typek, int dim, Int backend)
{
    if (backend == 1)  
        cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, dim);        
#ifdef USE_OMP             
    if (backend == 4) 
        ompFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, dim);        
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, dim);        
#endif                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, dim);        
#endif                  
}

inline void FullNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *atomtype, int *ilist, 
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int typek, int typel, int dim, Int backend)
{
    if (backend == 1)  
        cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, typel, dim);        
#ifdef USE_OMP             
    if (backend == 4) 
        ompFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, typel, dim);        
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, typel, dim);        
#endif                                                                              
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, typek, typel, dim);        
#endif                  
}

inline void HalfNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *ilist, int *alist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim, Int backend)
{
    if (backend == 1)  
        cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, alist, neighlist, neighnum, inum, jnum, dim);        
#ifdef USE_OMP             
    if (backend == 4) 
        ompHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, alist, neighlist, neighnum, inum, jnum, dim);        
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, alist, neighlist, neighnum, inum, jnum, dim);        
#endif                                                                                  
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, alist, neighlist, neighnum, inum, jnum, dim);        
#endif                  
}

inline void HalfNeighPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *atomtype, int *ilist, 
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim, Int backend)
{
    if (backend == 1)  
        cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#ifdef USE_OMP             
    if (backend == 4) 
        ompHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#endif                                                                                      
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, alist, neighlist, neighnum, inum, jnum, typej, dim);        
#endif                  
}

inline void NeighTripletList(int *tripletlist, int *tripletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum, Int backend)
{
    if (backend == 1)  
        cpuNeighTripletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#ifdef USE_OMP             
    if (backend == 4) 
        ompNeighTripletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipNeighTripletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuNeighTripletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#endif                  
}

inline void NeighQuadrupletList(int *tripletlist, int *tripletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum, Int backend)
{
    if (backend == 1)  
        cpuNeighQuadrupletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#ifdef USE_OMP             
    if (backend == 4) 
        ompNeighQuadrupletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipNeighQuadrupletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuNeighQuadrupletList(tripletlist, tripletnumsum, pairnum,  pairlist, ilist, alist, inum, jnum);
#endif                  
}

inline void NeighSingles(dstype *xi, dstype *qi, dstype *x, dstype *q, int *ai, 
      int *ti, int *ilist, int *atomtype, int inum, int ncq, int dim, Int backend)
{
    if (backend == 1)  
        cpuNeighSingles(xi, qi, x, q, ai, ti, ilist, atomtype, inum, ncq, dim);
#ifdef USE_OMP             
    if (backend == 4) 
        ompNeighSingles(xi, qi, x, q, ai, ti, ilist, atomtype, inum, ncq, dim);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipNeighSingles(xi, qi, x, q, ai, ti, ilist, atomtype, inum, ncq, dim);
#endif                                                                                          
#ifdef USE_CUDA            
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
#ifdef USE_OMP             
    if (backend == 4) 
        ompNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, alist, 
                atomtype, inum, jnum, ncq, dim);       
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, alist, 
                atomtype, inum, jnum, ncq, dim);       
#endif                                                                                          
#ifdef USE_CUDA            
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
#ifdef USE_OMP             
    if (backend == 4) 
        ompNeighTriplets(xij, xik, qi, qj, qk, x, q, ai, aj, ak, ti, tj, tk, tripletnum, tripletlist, tripletnumsum, 
                ilist, alist, atomtype, inum, ncq, dim);       
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipNeighTriplets(xij, xik, qi, qj, qk, x, q, ai, aj, ak, ti, tj, tk, tripletnum, tripletlist, tripletnumsum, 
                ilist, alist, atomtype, inum, ncq, dim);       
#endif                                                                                          
#ifdef USE_CUDA            
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
#ifdef USE_OMP             
    if (backend == 4) 
        ompNeighQuadruplets(xij, xik, xil, qi, qj, qk, ql, x, q, ai, aj, ak, al, ti, tj, tk, tl, 
                quadrupletnum, quadrupletlist, quadrupletnumsum, ilist, alist, atomtype, inum, ncq, dim);       
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipNeighQuadruplets(xij, xik, xil, qi, qj, qk, ql, x, q, ai, aj, ak, al, ti, tj, tk, tl, 
                quadrupletnum, quadrupletlist, quadrupletnumsum, ilist, alist, atomtype, inum, ncq, dim);       
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuNeighQuadruplets(xij, xik, xil, qi, qj, qk, ql, x, q, ai, aj, ak, al, ti, tj, tk, tl, 
                quadrupletnum, quadrupletlist, quadrupletnumsum, ilist, alist, atomtype, inum, ncq, dim);       
#endif                  
}

inline void PackIntProperty(dstype *buf, int *prop, int *ilist, int m, int mvalues, 
        int n, int nvalues, int inum, int backend)
{
    if (backend == 1)  
        cpuPackIntProperty(buf, prop, ilist, m, mvalues, n, nvalues, inum);
#ifdef USE_OMP             
    if (backend == 4) 
        ompPackIntProperty(buf, prop, ilist, m, mvalues, n, nvalues, inum);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipPackIntProperty(buf, prop, ilist, m, mvalues, n, nvalues, inum);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuPackIntProperty(buf, prop, ilist, m, mvalues, n, nvalues, inum);
#endif                      
}

inline void PackIntProperty(dstype *buf, int *prop, int *type, int *ilist, 
         int n, int nvalues, int inum, int backend)
{
    if (backend == 1)  
        cpuPackIntProperty(buf, prop, type, ilist, n, nvalues, inum);
#ifdef USE_OMP             
    if (backend == 4) 
        ompPackIntProperty(buf, prop, type, ilist, n, nvalues, inum);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipPackIntProperty(buf, prop, type, ilist, n, nvalues, inum);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuPackIntProperty(buf, prop, type, ilist, n, nvalues, inum);
#endif                      
}

inline void PackFloatProperty(dstype *buf, dstype *prop, int *ilist, int m, int mvalues, 
        int n, int nvalues, int inum, int backend)
{
    if (backend == 1)  
        cpuPackFloatProperty(buf, prop, ilist, m, mvalues, n, nvalues, inum);
#ifdef USE_OMP             
    if (backend == 4) 
        ompPackFloatProperty(buf, prop, ilist, m, mvalues, n, nvalues, inum);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipPackFloatProperty(buf, prop, ilist, m, mvalues, n, nvalues, inum);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuPackFloatProperty(buf, prop, ilist, m, mvalues, n, nvalues, inum);
#endif                      
}

inline void PackFloatProperty(dstype *buf, dstype *prop, dstype a, dstype b, int *ilist, 
        int m, int mvalues, int n, int nvalues, int inum, int backend)
{
    if (backend == 1)  
        cpuPackFloatProperty(buf, prop, a, b, ilist, m, mvalues, n, nvalues, inum);
#ifdef USE_OMP             
    if (backend == 4) 
        ompPackFloatProperty(buf, prop, a, b, ilist, m, mvalues, n, nvalues, inum);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipPackFloatProperty(buf, prop, a, b, ilist, m, mvalues, n, nvalues, inum);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuPackFloatProperty(buf, prop, a, b, ilist, m, mvalues, n, nvalues, inum);
#endif                      
}

inline void PackFloatProperty(dstype *buf, dstype *prop, int *type, int *ilist, 
         int n, int nvalues, int inum, int backend)
{
    if (backend == 1)  
        cpuPackFloatProperty(buf, prop, type, ilist, n, nvalues, inum);
#ifdef USE_OMP             
    if (backend == 4) 
        ompPackFloatProperty(buf, prop, type, ilist, n, nvalues, inum);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipPackFloatProperty(buf, prop, type, ilist, n, nvalues, inum);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuPackFloatProperty(buf, prop, type, ilist, n, nvalues, inum);
#endif                      
}

inline dstype ComputeMass(dstype *amass, dstype *mass, dstype *tmp, int *type, 
        int *ilist, int inum, int backend)
{
    if (backend == 1)  
        return cpuComputeMass(amass, mass, type, ilist, inum);
#ifdef USE_OMP             
    if (backend == 4) 
        return ompComputeMass(amass, mass, tmp, type, ilist, inum);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        return hipComputeMass(amass, mass, tmp, type, ilist, inum);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        return gpuComputeMass(amass, mass, tmp, type, ilist, inum);
#endif                      
    return 0;
}

inline void ComputeXCM(dstype *xcm, dstype *axcm, dstype *x, dstype *tmp, dstype *mass, dstype *box, dstype masstotal, int *ilist, int *type, int *image, int triclinic, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeXCM(xcm, x, mass, box, masstotal, ilist, type, image, triclinic, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeXCM(xcm, axcm, x, tmp, mass, box, masstotal, ilist, type, image, triclinic, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeXCM(xcm, axcm, x, tmp, mass, box, masstotal, ilist, type, image, triclinic, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeXCM(xcm, axcm, x, tmp, mass, box, masstotal, ilist, type, image, triclinic, dim, inum);
#endif
}

inline void ComputeVCM(dstype *vcm, dstype *avcm, dstype *v, dstype *tmp, dstype *mass, dstype masstotal, int *ilist, int *type, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeVCM(vcm, v, mass, masstotal, ilist, type, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeVCM(vcm, avcm, v, tmp, mass, masstotal, ilist, type, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeVCM(vcm, avcm, v, tmp, mass, masstotal, ilist, type, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeVCM(vcm, avcm, v, tmp, mass, masstotal, ilist, type, dim, inum);
#endif
}

inline dstype ComputeGyration(dstype * ag, dstype *xcm, dstype *x, dstype *tmp, dstype *mass, dstype *box, dstype masstotal, int *ilist, int *type, int *image, int triclinic, int dim, int inum, int backend)
{
	if (backend == 1)
		return cpuComputeGyration(xcm, x, mass, box, masstotal, ilist, type, image, triclinic, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		return ompComputeGyration(ag, xcm, x, tmp, mass, box, masstotal, ilist, type, image, triclinic, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		return hipComputeGyration(ag, xcm, x, tmp, mass, box, masstotal, ilist, type, image, triclinic, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		return gpuComputeGyration(ag, xcm, x, tmp, mass, box, masstotal, ilist, type, image, triclinic, dim, inum);
#endif
    return 0;
}

inline void ComputeAngmom(dstype *lmom, dstype *p, dstype *xcm, dstype *x, dstype *v, dstype *tmp, dstype *mass, dstype *box, int *ilist, int *type, int *image, int triclinic, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeAngmom(lmom, xcm, x, v, mass, box, ilist, type, image, triclinic, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeAngmom(lmom, p, xcm, x, v, tmp, mass, box, ilist, type, image, triclinic, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeAngmom(lmom, p, xcm, x, v, tmp, mass, box, ilist, type, image, triclinic, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeAngmom(lmom, p, xcm, x, v, tmp, mass, box, ilist, type, image, triclinic, dim, inum);
#endif
}

inline void ComputeTorque(dstype *tq, dstype *q, dstype *xcm, dstype *x, dstype *f, dstype *tmp, dstype *box, int *ilist, int *image, int triclinic, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeTorque(tq, xcm, x, f, box, ilist, image, triclinic, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeTorque(tq, q, xcm, x, f, tmp, box, ilist, image, triclinic, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeTorque(tq, q, xcm, x, f, tmp, box, ilist, image, triclinic, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeTorque(tq, q, xcm, x, f, tmp, box, ilist, image, triclinic, dim, inum);
#endif
}

inline void ComputeInertia(dstype *inertia, dstype *ione, dstype *xcm, dstype *x, dstype *tmp, dstype *mass, dstype *box, int *ilist, int *type, int *image, int triclinic, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeInertia(inertia, xcm, x, mass, box, ilist, type, image, triclinic, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeInertia(inertia, ione, xcm, x, tmp, mass, box, ilist, type, image, triclinic, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeInertia(inertia, ione, xcm, x, tmp, mass, box, ilist, type, image, triclinic, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeInertia(inertia, ione, xcm, x, tmp, mass, box, ilist, type, image, triclinic, dim, inum);
#endif
}

inline void ComputeKEAtom(dstype *ke, dstype *mass, dstype *v, dstype mvv2e, int *type, int *ilist, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeKEAtom(ke, mass, v, mvv2e, type, ilist, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeKEAtom(ke, mass, v, mvv2e, type, ilist, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeKEAtom(ke, mass, v, mvv2e, type, ilist, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeKEAtom(ke, mass, v, mvv2e, type, ilist, dim, inum);
#endif
}

inline void ComputeStressAtom(dstype *stress, dstype *mass, dstype *vatom, dstype *v, dstype mvv2e, dstype nktv2p, int *type, int *ilist,  int vflag, int keflag, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeStressAtom(stress, mass, vatom, v, mvv2e, nktv2p, type, ilist,  vflag, keflag, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeStressAtom(stress, mass, vatom, v, mvv2e, nktv2p, type, ilist,  vflag, keflag, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeStressAtom(stress, mass, vatom, v, mvv2e, nktv2p, type, ilist,  vflag, keflag, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeStressAtom(stress, mass, vatom, v, mvv2e, nktv2p, type, ilist,  vflag, keflag, dim, inum);
#endif
}

inline void ComputeCentroidStressAtom(dstype *stress, dstype *mass, dstype *cvatom, dstype *v, dstype mvv2e, dstype nktv2p, int *type, int *ilist,  int vflag, int keflag, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeCentroidStressAtom(stress, mass, cvatom, v, mvv2e, nktv2p, type, ilist,  vflag, keflag, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeCentroidStressAtom(stress, mass, cvatom, v, mvv2e, nktv2p, type, ilist,  vflag, keflag, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeCentroidStressAtom(stress, mass, cvatom, v, mvv2e, nktv2p, type, ilist,  vflag, keflag, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeCentroidStressAtom(stress, mass, cvatom, v, mvv2e, nktv2p, type, ilist,  vflag, keflag, dim, inum);
#endif
}

inline void ComputeDisplaceAtom(dstype *displace, dstype *x, dstype *xoriginal, dstype *box, int *pbc, int *ilist,  int triclinic, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeDisplaceAtom(displace, x, xoriginal, box, pbc, ilist,  triclinic, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeDisplaceAtom(displace, x, xoriginal, box, pbc, ilist,  triclinic, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeDisplaceAtom(displace, x, xoriginal, box, pbc, ilist,  triclinic, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeDisplaceAtom(displace, x, xoriginal, box, pbc, ilist,  triclinic, dim, inum);
#endif
}

inline void ComputeTempSymTensor(dstype *ke_tensor, dstype *stress, dstype *v, dstype *tmp, dstype *mass, dstype tfactor, int *type, int *ilist, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeTempSymTensor(ke_tensor, v, mass, tfactor, type, ilist, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeTempSymTensor(ke_tensor, stress, v, tmp, mass, tfactor, type, ilist, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeTempSymTensor(ke_tensor, stress, v, tmp, mass, tfactor, type, ilist, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeTempSymTensor(ke_tensor, stress, v, tmp, mass, tfactor, type, ilist, dim, inum);
#endif
}

inline dstype ComputeTempScalar(dstype *ke, dstype *v, dstype *tmp, dstype *mass, dstype tfactor, int *type, int *ilist, int dim, int inum, int backend)
{
	if (backend == 1)
		return cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		return ompComputeTempScalar(ke, v, tmp, mass, tfactor, type, ilist, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		return hipComputeTempScalar(ke, v, tmp, mass, tfactor, type, ilist, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		return gpuComputeTempScalar(ke, v, tmp, mass, tfactor, type, ilist, dim, inum);
#endif
    return 0;
}


inline dstype ComputePressureScalar(dstype *virial, dstype volume, dstype temp, dstype tempdof, dstype boltz, dstype nktv2p, int dim, int backend)
{
	if (backend == 1)
		return cpuComputePressureScalar(virial, volume, temp, tempdof, boltz, nktv2p, dim);
#ifdef USE_OMP
	if (backend == 4)
		return ompComputePressureScalar(virial, volume, temp, tempdof, boltz, nktv2p, dim);
#endif
#ifdef USE_HIP
	if (backend == 3)
		return hipComputePressureScalar(virial, volume, temp, tempdof, boltz, nktv2p, dim);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		return gpuComputePressureScalar(virial, volume, temp, tempdof, boltz, nktv2p, dim);
#endif
    return 0;
}


inline void ComputePressureSymTensor(dstype *vector, dstype *virial, dstype *ke_tensor, dstype volume, dstype nktv2p, int dim, int backend)
{
	if (backend == 1)
		cpuComputePressureSymTensor(vector, virial, ke_tensor, volume, nktv2p, dim);
#ifdef USE_OMP
	if (backend == 4)
		ompComputePressureSymTensor(vector, virial, ke_tensor, volume, nktv2p, dim);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputePressureSymTensor(vector, virial, ke_tensor, volume, nktv2p, dim);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputePressureSymTensor(vector, virial, ke_tensor, volume, nktv2p, dim);
#endif
}

inline void ComputeHeatFlux(dstype *vector, dstype *jc, dstype *ke, dstype *pe, dstype *stress, dstype *v, dstype *tmp, dstype nktv2p, int *ilist,  int pressatomflag, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeHeatFlux(vector, ke, pe, stress, v, nktv2p, ilist,  pressatomflag, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeHeatFlux(vector, jc, ke, pe, stress, v, tmp, nktv2p, ilist,  pressatomflag, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeHeatFlux(vector, jc, ke, pe, stress, v, tmp, nktv2p, ilist,  pressatomflag, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeHeatFlux(vector, jc, ke, pe, stress, v, tmp, nktv2p, ilist,  pressatomflag, dim, inum);
#endif
}

inline void ComputeOrientOrderAtom(dstype *qnarray, dstype *x, dstype *rlist, dstype *cglist, dstype *fac, dstype *qnm_r, dstype *qnm_i, dstype *distsq, dstype cutsq, dstype MY_EPSILON, dstype QEPSILON, dstype MY_4PI, int *neighlist, int *neighnum,  int *ilist, int *qlist, int *nearest,  int nqlist, int qmax, int wlflag, int wlhatflag,  int qlcompflag, int iqlcomp, int qlcomp, int nnn, int jnum, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuComputeOrientOrderAtom(qnarray, x, rlist, cglist, fac, qnm_r, qnm_i, distsq, cutsq, MY_EPSILON, QEPSILON, MY_4PI, neighlist, neighnum,  ilist, qlist, nearest,  nqlist, qmax, wlflag, wlhatflag,  qlcompflag, iqlcomp, qlcomp, nnn, jnum, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeOrientOrderAtom(qnarray, x, rlist, cglist, fac, qnm_r, qnm_i, distsq, cutsq, MY_EPSILON, QEPSILON, MY_4PI, neighlist, neighnum,  ilist, qlist, nearest,  nqlist, qmax, wlflag, wlhatflag,  qlcompflag, iqlcomp, qlcomp, nnn, jnum, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeOrientOrderAtom(qnarray, x, rlist, cglist, fac, qnm_r, qnm_i, distsq, cutsq, MY_EPSILON, QEPSILON, MY_4PI, neighlist, neighnum,  ilist, qlist, nearest,  nqlist, qmax, wlflag, wlhatflag,  qlcompflag, iqlcomp, qlcomp, nnn, jnum, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeOrientOrderAtom(qnarray, x, rlist, cglist, fac, qnm_r, qnm_i, distsq, cutsq, MY_EPSILON, QEPSILON, MY_4PI, neighlist, neighnum,  ilist, qlist, nearest,  nqlist, qmax, wlflag, wlhatflag,  qlcompflag, iqlcomp, qlcomp, nnn, jnum, dim, inum);
#endif
}

inline void ComputeCoordAtomCutoff(int *cvec, dstype *x, dstype *rcutsq, int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, int *jgroupbit, int dim, int ntypes, int jnum, int inum, int backend)
{
	if (backend == 1)
		cpuComputeCoordAtomCutoff(cvec, x, rcutsq, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, dim, ntypes, jnum, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeCoordAtomCutoff(cvec, x, rcutsq, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, dim, ntypes, jnum, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeCoordAtomCutoff(cvec, x, rcutsq, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, dim, ntypes, jnum, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeCoordAtomCutoff(cvec, x, rcutsq, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, dim, ntypes, jnum, inum);
#endif
}

inline void ComputeCoordAtomCutoff(int *carray, dstype *x, dstype *rcutsq, int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, int *jgroupbit, int ncol, int dim, int ntypes, int jnum, int inum, int backend)
{
	if (backend == 1)
		cpuComputeCoordAtomCutoff(carray, x, rcutsq, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, ncol, dim, ntypes, jnum, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeCoordAtomCutoff(carray, x, rcutsq, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, ncol, dim, ntypes, jnum, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeCoordAtomCutoff(carray, x, rcutsq, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, ncol, dim, ntypes, jnum, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeCoordAtomCutoff(carray, x, rcutsq, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, ncol, dim, ntypes, jnum, inum);
#endif
}

inline void ComputeCoordAtomOrient(int *cvec, dstype *x, dstype *rcutsq, dstype *normv, dstype threshold, int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, int *jgroupbit, int dim, int ntypes, int nqlist, int ncol, int l, int jnum, int inum, int backend)
{
	if (backend == 1)
		cpuComputeCoordAtomOrient(cvec, x, rcutsq, normv, threshold, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, dim, ntypes, nqlist, ncol, l, jnum, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeCoordAtomOrient(cvec, x, rcutsq, normv, threshold, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, dim, ntypes, nqlist, ncol, l, jnum, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeCoordAtomOrient(cvec, x, rcutsq, normv, threshold, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, dim, ntypes, nqlist, ncol, l, jnum, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeCoordAtomOrient(cvec, x, rcutsq, normv, threshold, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, dim, ntypes, nqlist, ncol, l, jnum, inum);
#endif
}

inline void ComputeMSD(dstype *msd, dstype *vec, dstype *x, dstype *xoriginal, dstype *box, dstype *xcm, dstype *tmp, int *ilist, int *image, int naverage, int avflag, int triclinic, int nmsd, int dim,  int inum, int backend)
{
	if (backend == 1)
		cpuComputeMSD(msd, x, xoriginal, box, xcm, ilist, image, naverage, avflag, triclinic, nmsd, dim,  inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeMSD(msd, vec, x, xoriginal, box, xcm, tmp, ilist, image, naverage, avflag, triclinic, nmsd, dim,  inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeMSD(msd, vec, x, xoriginal, box, xcm, tmp, ilist, image, naverage, avflag, triclinic, nmsd, dim,  inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeMSD(msd, vec, x, xoriginal, box, xcm, tmp, ilist, image, naverage, avflag, triclinic, nmsd, dim,  inum);
#endif
}

inline void ComputeVACF(dstype *vacf, dstype *vec, dstype *v, dstype *voriginal, dstype *tmp, int *ilist, int nvacf, int dim,  int inum, int backend)
{
	if (backend == 1)
		cpuComputeVACF(vacf, v, voriginal, ilist, nvacf, dim,  inum);
#ifdef USE_OMP
	if (backend == 4)
		ompComputeVACF(vacf, vec, v, voriginal, tmp, ilist, nvacf, dim,  inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipComputeVACF(vacf, vec, v, voriginal, tmp, ilist, nvacf, dim,  inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuComputeVACF(vacf, vec, v, voriginal, tmp, ilist, nvacf, dim,  inum);
#endif
}

inline void Velocity(dstype *x, dstype *v, dstype *f, dstype *box, dstype *xcm, dstype *vcm, dstype *mass, dstype *second, dstype *omega, dstype *vext, dstype *v_lo, dstype *v_hi, dstype *coord_lo, dstype *coord_hi, dstype t_desired, dstype t_current, int *seed, int *save, int *map, int *image, int *type, int *coord_dim, int *vdim, int sum_flag, int dist_flag, int loop_flag, int rotation_flag, int momentum_flag, int triclinic, int dim, int mpiRank, int vmode, int nlocal, int natoms, int backend)
{
	if (backend == 1)
		cpuVelocity(x, v, f, box, xcm, vcm, mass, second, omega, vext, v_lo, v_hi, coord_lo, coord_hi, t_desired, t_current, seed, save, map, image, type, coord_dim, vdim, sum_flag, dist_flag, loop_flag, rotation_flag, momentum_flag, triclinic, dim, mpiRank, vmode, nlocal, natoms);
#ifdef USE_OMP
	if (backend == 4)
		ompVelocity(x, v, f, box, xcm, vcm, mass, second, omega, vext, v_lo, v_hi, coord_lo, coord_hi, t_desired, t_current, seed, save, map, image, type, coord_dim, vdim, sum_flag, dist_flag, loop_flag, rotation_flag, momentum_flag, triclinic, dim, mpiRank, vmode, nlocal, natoms);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipVelocity(x, v, f, box, xcm, vcm, mass, second, omega, vext, v_lo, v_hi, coord_lo, coord_hi, t_desired, t_current, seed, save, map, image, type, coord_dim, vdim, sum_flag, dist_flag, loop_flag, rotation_flag, momentum_flag, triclinic, dim, mpiRank, vmode, nlocal, natoms);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuVelocity(x, v, f, box, xcm, vcm, mass, second, omega, vext, v_lo, v_hi, coord_lo, coord_hi, t_desired, t_current, seed, save, map, image, type, coord_dim, vdim, sum_flag, dist_flag, loop_flag, rotation_flag, momentum_flag, triclinic, dim, mpiRank, vmode, nlocal, natoms);
#endif
}

inline void SetVelocityInitialIntegrate(dstype *x, dstype *v, dstype *f, dstype *mass, dstype *fparam, dstype dtf, dstype dtv, int *type, int *ilist, int *iparam, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuSetVelocityInitialIntegrate(x, v, f, mass, fparam, dtf, dtv, type, ilist, iparam, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompSetVelocityInitialIntegrate(x, v, f, mass, fparam, dtf, dtv, type, ilist, iparam, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipSetVelocityInitialIntegrate(x, v, f, mass, fparam, dtf, dtv, type, ilist, iparam, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuSetVelocityInitialIntegrate(x, v, f, mass, fparam, dtf, dtv, type, ilist, iparam, dim, inum);
#endif
}

inline void SetVelocityFinalIntegrate(dstype *x, dstype *v, dstype *f, dstype *mass, dstype dtf, int *type, int *ilist, int *iparam, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuSetVelocityFinalIntegrate(x, v, f, mass, dtf, type, ilist, iparam, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompSetVelocityFinalIntegrate(x, v, f, mass, dtf, type, ilist, iparam, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipSetVelocityFinalIntegrate(x, v, f, mass, dtf, type, ilist, iparam, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuSetVelocityFinalIntegrate(x, v, f, mass, dtf, type, ilist, iparam, dim, inum);
#endif
}

inline void InitialIntegrate(dstype *x, dstype *v, dstype *f, dstype *mass, dstype *dtarray, 
        dstype *tarray, dstype *eta_mass, dstype *eta, dstype *eta_dot, dstype *eta_dotdot, 
        dstype *ke, dstype *tmp, dstype vlimitsq, int *type, int *ilist, int eta_mass_flag, 
        int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuInitialIntegrate(x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, vlimitsq, type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, mode, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompInitialIntegrate(x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, ke, tmp, vlimitsq, type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, mode, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipInitialIntegrate(x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, ke, tmp, vlimitsq, type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, mode, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuInitialIntegrate(x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, ke, tmp, vlimitsq, type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, mode, dim, inum);
#endif
}

inline void FinalIntegrate(dstype *x, dstype *v, dstype *f, dstype *mass, dstype *dtarray, dstype *tarray, 
        dstype *eta_mass, dstype *eta, dstype *eta_dot, dstype *eta_dotdot, dstype *ke, dstype *tmp, 
        dstype vlimitsq, int *type, int *ilist, int eta_mass_flag, int biasflag, int mtchain, 
        int nc_tchain, int mode, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFinalIntegrate(x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, vlimitsq, type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, mode, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFinalIntegrate(x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, 
                ke, tmp, vlimitsq, type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, mode, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFinalIntegrate(x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, 
                ke, tmp, vlimitsq, type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, mode, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFinalIntegrate(x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, 
                ke, tmp, vlimitsq, type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, mode, dim, inum);
#endif
}

inline dstype VelocityRescalingThermostat(dstype *v, dstype *mass, dstype *dtarray, dstype *tarray, 
        dstype *ke, dstype *tmp, dstype *second, dstype energy, int *type, int *ilist, int *seed, 
        int *save, int biasflag, int mode, int dim, int inum, int backend)
{
	if (backend == 1)
		return cpuVelocityRescalingThermostat(v, mass, dtarray, tarray, 
                second, energy, type, ilist, seed, save, biasflag, mode, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		return ompVelocityRescalingThermostat(v, mass, dtarray, tarray, ke, tmp,
                second, energy, type, ilist, seed, save, biasflag, mode, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		return hipVelocityRescalingThermostat(v, mass, dtarray, tarray, ke, tmp,
                second, energy, type, ilist, seed, save, biasflag, mode, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		return gpuVelocityRescalingThermostat(v, mass, dtarray, tarray, ke, tmp,
                second, energy, type, ilist, seed, save, biasflag, mode, dim, inum);
#endif
    return 0;
}

inline void PBC(dstype *x, dstype *v, int *image, dstype *boxhi, dstype *boxlo, dstype *hi_lambda, dstype *lo_lambda, dstype *h, dstype *h_inv, dstype *h_rate, int *pbc, int vdeform, int triclinic, int dim, int nlocal, int backend)
{
	if (backend == 1)
		cpuPBC(x, v, image, boxhi, boxlo, hi_lambda, lo_lambda, h, h_inv, h_rate, pbc, vdeform, triclinic, dim, nlocal);
#ifdef USE_OMP
	if (backend == 4)
		ompPBC(x, v, image, boxhi, boxlo, hi_lambda, lo_lambda, h, h_inv, h_rate, pbc, vdeform, triclinic, dim, nlocal);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipPBC(x, v, image, boxhi, boxlo, hi_lambda, lo_lambda, h, h_inv, h_rate, pbc, vdeform, triclinic, dim, nlocal);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuPBC(x, v, image, boxhi, boxlo, hi_lambda, lo_lambda, h, h_inv, h_rate, pbc, vdeform, triclinic, dim, nlocal);
#endif
}

inline void FixSetForce(dstype *x, dstype *v, dstype *f, dstype *eatom, dstype *vatom, dstype *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFixSetForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFixSetForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFixSetForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFixSetForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
}

inline void FixLineForce(dstype *x, dstype *v, dstype *f, dstype *eatom, dstype *vatom, dstype *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFixLineForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFixLineForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFixLineForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFixLineForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
}

inline void FixPlaneForce(dstype *x, dstype *v, dstype *f, dstype *eatom, dstype *vatom, dstype *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFixPlaneForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFixPlaneForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFixPlaneForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFixPlaneForce(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
}

inline void FixAddForce(dstype *x, dstype *v, dstype *f, dstype *eatom, dstype *vatom, dstype *fparam, dstype *box, int *iparam, int *ilist, int *image, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFixAddForce(x, v, f, eatom, vatom, fparam, box, iparam, ilist, image, triclinic, eflag_atom, vflag_atom, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFixAddForce(x, v, f, eatom, vatom, fparam, box, iparam, ilist, image, triclinic, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFixAddForce(x, v, f, eatom, vatom, fparam, box, iparam, ilist, image, triclinic, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFixAddForce(x, v, f, eatom, vatom, fparam, box, iparam, ilist, image, triclinic, eflag_atom, vflag_atom, dim, inum);
#endif
}

inline void FixDragForce(dstype *x, dstype *v, dstype *f, dstype *eatom, dstype *vatom, dstype *fparam, dstype *box, int *iparam, int *ilist, int *pbc, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFixDragForce(x, v, f, eatom, vatom, fparam, box, iparam, ilist, pbc, triclinic, eflag_atom, vflag_atom, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFixDragForce(x, v, f, eatom, vatom, fparam, box, iparam, ilist, pbc, triclinic, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFixDragForce(x, v, f, eatom, vatom, fparam, box, iparam, ilist, pbc, triclinic, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFixDragForce(x, v, f, eatom, vatom, fparam, box, iparam, ilist, pbc, triclinic, eflag_atom, vflag_atom, dim, inum);
#endif
}

inline void FixWallReflect(dstype *x, dstype *v, dstype *f, dstype *eatom, dstype *vatom, dstype *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFixWallReflect(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFixWallReflect(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFixWallReflect(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFixWallReflect(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
}

inline void FixWallHarmonic(dstype *x, dstype *v, dstype *f, dstype *eatom, dstype *vatom, dstype *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFixWallHarmonic(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFixWallHarmonic(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFixWallHarmonic(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFixWallHarmonic(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
}

inline void FixWallLJ93(dstype *x, dstype *v, dstype *f, dstype *eatom, dstype *vatom, dstype *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFixWallLJ93(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFixWallLJ93(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFixWallLJ93(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFixWallLJ93(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
}

inline void FixWallLJ126(dstype *x, dstype *v, dstype *f, dstype *eatom, dstype *vatom, dstype *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFixWallLJ126(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFixWallLJ126(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFixWallLJ126(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFixWallLJ126(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
}

inline void FixWallLJ1043(dstype *x, dstype *v, dstype *f, dstype *eatom, dstype *vatom, dstype *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFixWallLJ1043(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFixWallLJ1043(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFixWallLJ1043(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFixWallLJ1043(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
}

inline void FixWallMorse(dstype *x, dstype *v, dstype *f, dstype *eatom, dstype *vatom, dstype *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum, int backend)
{
	if (backend == 1)
		cpuFixWallMorse(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#ifdef USE_OMP
	if (backend == 4)
		ompFixWallMorse(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipFixWallMorse(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuFixWallMorse(x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum);
#endif
}

inline void SphericalHarmonicsBessel(dstype *Sr, dstype *Si, dstype *xij, dstype *x0, 
        dstype *P, dstype *tmp, dstype *f, dstype *fac, dstype pi, Int L, Int K, Int N, Int backend)
{
    if (backend == 1)  
        cpuSphericalHarmonicsBessel(Sr, Si, xij, x0, P, tmp, f, fac, pi, L, K, N);    
#ifdef USE_OMP             
    if (backend == 4) 
        ompSphericalHarmonicsBessel(Sr, Si, xij, x0, P, tmp, f, fac, pi, L, K, N);    
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipSphericalHarmonicsBessel(Sr, Si, xij, x0, P, tmp, f, fac, pi, L, K, N);    
#endif                                                                                              
#ifdef USE_CUDA            
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
#ifdef USE_OMP             
    if (backend == 4) 
        ompRadialSphericalHarmonicsSpectrum(c, ar, ai, sr, si, cg, indk, indl, indm, 
                rowm, Nnb, Nub, Ncg, Na, L, K, spectrum); 
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipRadialSphericalHarmonicsSpectrum(c, ar, ai, sr, si, cg, indk, indl, indm, 
                rowm, Nnb, Nub, Ncg, Na, L, K, spectrum); 
#endif                                                                                              
#ifdef USE_CUDA            
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
#ifdef USE_OMP             
    if (backend == 4) 
        ompSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, xij,
                 x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipNSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, xij,
                 x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);
#endif                                                                                              
#ifdef USE_CUDA            
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
#ifdef USE_OMP             
    if (backend == 4) 
        ompRadialSphericalHarmonicsSpectrumDeriv(cd, ar, ai, srx, six, sry, siy, srz, siz, cg, indk, 
                indl, indm, rowm, Nnb, Nub, Ncg, Na, L, K, spectrum);                        
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipRadialSphericalHarmonicsSpectrumDeriv(cd, ar, ai, srx, six, sry, siy, srz, siz, cg, indk, 
                indl, indm, rowm, Nnb, Nub, Ncg, Na, L, K, spectrum);                        
#endif                                                                                          
#ifdef USE_CUDA            
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
#ifdef USE_OMP             
    if (backend == 4) 
        ompForceDecomposition(f, fij, ai, aj, inum, ijnum, Nbf);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipForceDecomposition(f, fij, ai, aj, inum, ijnum, Nbf);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  {
        gpuForceDecomposition(f, fij, ai, aj, inum, ijnum, Nbf);
    }
#endif                  
}

inline void ForceDecomposition(dstype *f, dstype *fij, int *ai, int *aj, int ijnum, Int backend)
{
    if (backend == 1)  
        cpuForceDecomposition(f, fij, ai, aj, ijnum);
#ifdef USE_OMP             
    if (backend == 4) 
        ompForceDecomposition(f, fij, ai, aj, ijnum);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipForceDecomposition(f, fij, ai, aj, ijnum);
#endif                                                                                          
#ifdef USE_CUDA            
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
#ifdef USE_OMP             
    if (backend == 4) 
        ompCenterAtomDecomposition(f, fij, ilist, anumsum, inum, ijnum, Na, Nbf);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipCenterAtomDecomposition(f, fij, ilist, anumsum, inum, ijnum, Na, Nbf);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuCenterAtomDecomposition(f, fij, ilist, anumsum, inum, ijnum, Na, Nbf);    
#endif                  
}

inline void CenterAtomDecomposition(dstype *e, dstype *ei, int *ilist, int Na, Int backend)
{
    if (backend == 1)  
        cpuCenterAtomDecomposition(e, ei, ilist, Na);
#ifdef USE_OMP             
    if (backend == 4) 
        ompCenterAtomDecomposition(e, ei, ilist, Na);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipCenterAtomDecomposition(e, ei, ilist, Na);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuCenterAtomDecomposition(e, ei, ilist, Na);    
#endif                  
}

inline void CenterAtomDecomposition(dstype *f, dstype *fij, int *ilist, int *anumsum, 
        int inum, Int backend)
{
    if (backend == 1)  
        cpuCenterAtomDecomposition(f, fij, ilist, anumsum, inum);
#ifdef USE_OMP             
    if (backend == 4) 
        ompCenterAtomDecomposition(f, fij, ilist, anumsum, inum);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipCenterAtomDecomposition(f, fij, ilist, anumsum, inum);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuCenterAtomDecomposition(f, fij, ilist, anumsum, inum);    
#endif                  
}

inline void NeighborAtomDecomposition(dstype *f, dstype *fij, int *ilist, int *anumsum, int *index, 
        int inum, int ijnum, int jnum, int Nbf, Int backend)
{
    if (backend == 1)  
        cpuNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum, ijnum, jnum, Nbf);
#ifdef USE_OMP             
    if (backend == 4) 
        ompNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum, ijnum, jnum, Nbf);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum, ijnum, jnum, Nbf);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum, ijnum, jnum, Nbf);
#endif                  
}

inline void NeighborAtomDecomposition(dstype *f, dstype *fij, int *ilist, int *anumsum, int *index, 
        int inum, Int backend)
{
    if (backend == 1)  
        cpuNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum);
#ifdef USE_OMP             
    if (backend == 4) 
        ompNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum);
#endif               
#ifdef USE_HIP            
    if (backend == 3) 
        hipNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum);
#endif                                                                                          
#ifdef USE_CUDA            
    if (backend == 2)  
        gpuNeighborAtomDecomposition(f, fij, ilist, anumsum, index, inum);    
#endif                  
}

inline void InitSna(dstype *rootpqarray, dstype *cglist, dstype *factorial, 
        int *idx_max, int *idxz, int *idxz_block, int *idxb, int *idxb_block, 
        int *idxu_block, int *idxcg_block, int twojmax, int backend)
{
    if (backend == 1)
        cpuInitSna(rootpqarray, cglist, factorial, idx_max, idxz, idxz_block, idxb, idxb_block, 
                idxu_block, idxcg_block, twojmax);
#ifdef USE_OMP  
    if (backend == 4)
        ompInitSna(rootpqarray, cglist, factorial, idx_max, idxz, idxz_block, idxb, idxb_block, 
                idxu_block, idxcg_block, twojmax);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipInitSna(rootpqarray, cglist, factorial, idx_max, idxz, idxz_block, idxb, idxb_block, 
                idxu_block, idxcg_block, twojmax);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuInitSna(rootpqarray, cglist, factorial, idx_max, idxz, idxz_block, idxb, idxb_block, 
                idxu_block, idxcg_block, twojmax);
#endif                          
}

inline void ComputeUij(dstype *ulist_r, dstype *ulist_i, dstype *rootpqarray, dstype *rij, 
        dstype *radelem, dstype rmin0, dstype rfac0, dstype rcutfac, int *idxu_block, 
        int *type, int *ai, int *aj, int twojmax, int idxu_max, int ijnum, int backend)
{
    if (backend == 1)
        cpuComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, idxu_block, 
                        type, ai, aj, twojmax, idxu_max, ijnum);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, idxu_block, 
                        type, ai, aj, twojmax, idxu_max, ijnum);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, idxu_block, 
                        type, ai, aj, twojmax, idxu_max, ijnum);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, idxu_block, 
                        type, ai, aj, twojmax, idxu_max, ijnum);
#endif                          
}

inline void ZeroUarraytot(dstype *ulisttot_r, dstype *ulisttot_i, dstype wself, int *idxu_block, int *type,
        int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
        int twojmax, int inum, int backend)
{
    if (backend == 1)
        cpuZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, type,
                    map, ai, wselfall_flag, chemflag, idxu_max, nelements, twojmax, inum);
#ifdef USE_OMP  
    if (backend == 4)
        ompZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, type,
                    map, ai, wselfall_flag, chemflag, idxu_max, nelements, twojmax, inum);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, type,
                    map, ai, wselfall_flag, chemflag, idxu_max, nelements, twojmax, inum);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, type,
                    map, ai, wselfall_flag, chemflag, idxu_max, nelements, twojmax, inum);
#endif                          
}

inline void AddUarraytot(dstype *ulisttot_r, dstype *ulisttot_i, dstype *ulist_r, dstype *ulist_i, dstype *rij, 
        dstype *wjelem, dstype *radelem, dstype rmin0, dstype rcutfac, int *idxu_block, int *ilist, 
        int *type, int *pairnum, int *pairnumsum, int *map, int *tj, int twojmax, int idxu_max, 
        int nelements, int inum, int switch_flag, int chemflag, int backend)
{
    if (backend == 1)
        cpuAddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, 
                rcutfac, idxu_block, ilist, type, pairnum, pairnumsum, map, tj, twojmax, 
                idxu_max, nelements, inum, switch_flag, chemflag);
#ifdef USE_OMP  
    if (backend == 4)
        ompAddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, 
                rcutfac, idxu_block, ilist, type, pairnum, pairnumsum, map, tj, twojmax, 
                idxu_max, nelements, inum, switch_flag, chemflag);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipAddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, 
                rcutfac, idxu_block, ilist, type, pairnum, pairnumsum, map, tj, twojmax, 
                idxu_max, nelements, inum, switch_flag, chemflag);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuAddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, 
                rcutfac, idxu_block, ilist, type, pairnum, pairnumsum, map, tj, twojmax, 
                idxu_max, nelements, inum, switch_flag, chemflag);
#endif                          
}

inline void ComputeBeta2(dstype *beta, dstype *bispectrum, dstype *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag, int backend)
{
    if (backend == 1)
        cpuComputeBeta2(beta, bispectrum, coeffelem, ilist, map, type, inum, ncoeff, ncoeffall, 
                quadraticflag);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeBeta2(beta, bispectrum, coeffelem, ilist, map, type, inum, ncoeff, ncoeffall, 
                quadraticflag);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeBeta2(beta, bispectrum, coeffelem, ilist, map, type, inum, ncoeff, ncoeffall, 
                quadraticflag);
#endif            
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeBeta2(beta, bispectrum, coeffelem, ilist, map, type, inum, ncoeff, ncoeffall, 
                quadraticflag);
#endif                          
}

inline void ComputeZi(dstype *zlist_r, dstype *zlist_i, dstype *ulisttot_r, dstype *ulisttot_i, dstype *cglist,
        int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, int idxz_max, int nelements, 
        int bnorm_flag, int inum, int backend)
{
    if (backend == 1)
        cpuComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
                        idxcg_block, twojmax, idxu_max, idxz_max, nelements, bnorm_flag, inum);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
                        idxcg_block, twojmax, idxu_max, idxz_max, nelements, bnorm_flag, inum);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
                        idxcg_block, twojmax, idxu_max, idxz_max, nelements, bnorm_flag, inum);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
                        idxcg_block, twojmax, idxu_max, idxz_max, nelements, bnorm_flag, inum);
#endif                          
}

inline void ComputeYi(dstype *ylist_r, dstype *ylist_i, dstype *ulisttot_r, dstype *ulisttot_i, 
        dstype *cglist, dstype *betalist, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag,
        int ncoeff, int inum, int backend)
{
    if (backend == 1)
        cpuComputeYi(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, betalist, idxz, idxb_block, idxu_block, 
                    idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, nelements, bnorm_flag, ncoeff, inum);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeYi(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, betalist, idxz, idxb_block, idxu_block, 
                    idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, nelements, bnorm_flag, ncoeff, inum);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeYi(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, betalist, idxz, idxb_block, idxu_block, 
                    idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, nelements, bnorm_flag, ncoeff, inum);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeYi(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, betalist, idxz, idxb_block, idxu_block, 
                    idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, nelements, bnorm_flag, ncoeff, inum);
#endif                          
}

inline void ComputeBi(dstype *blist, dstype *zlist_r, dstype *zlist_i, dstype *ulisttot_r, dstype *ulisttot_i, 
        dstype *bzero,int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum, int backend)
{
    if (backend == 1)
        cpuComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, ilist, type, map, 
                idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, nelements,
                bzero_flag, wselfall_flag, chemflag, inum);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, ilist, type, map, 
                idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, nelements,
                bzero_flag, wselfall_flag, chemflag, inum);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, ilist, type, map, 
                idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, nelements,
                bzero_flag, wselfall_flag, chemflag, inum);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, ilist, type, map, 
                idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, nelements,
                bzero_flag, wselfall_flag, chemflag, inum);
#endif                          
}

inline void ComputeDbidrj(dstype *dblist, dstype *zlist_r, dstype *zlist_i, dstype *dulist_r,
        dstype *dulist_i, int *idxb, int *idxu_block, int *idxz_block, int *type, int *map, 
        int *ai, int *aj, int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements,
        int bnorm_flag, int chemflag, int ijnum, int backend)
{
    if (backend == 1)
        cpuComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, idxb, idxu_block, 
                idxz_block, type, map, ai, aj, twojmax, idxb_max, idxu_max, idxz_max, 
                nelements, bnorm_flag, chemflag, ijnum);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, idxb, idxu_block, 
                idxz_block, type, map, ai, aj, twojmax, idxb_max, idxu_max, idxz_max, 
                nelements, bnorm_flag, chemflag, ijnum);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, idxb, idxu_block, 
                idxz_block, type, map, ai, aj, twojmax, idxb_max, idxu_max, idxz_max, 
                nelements, bnorm_flag, chemflag, ijnum);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, idxb, idxu_block, 
                idxz_block, type, map, ai, aj, twojmax, idxb_max, idxu_max, idxz_max, 
                nelements, bnorm_flag, chemflag, ijnum);
#endif                          
}

inline void ComputeDuijdrj(dstype *dulist_r, dstype *dulist_i, dstype *ulist_r, dstype *ulist_i, dstype *rootpqarray, 
        dstype* rij, dstype *wjelem, dstype *radelem, dstype rmin0, dstype rfac0, dstype rcutfac, int *idxu_block, 
        int *type, int *ai, int *aj, int twojmax, int idxu_max, int ijnum, int switch_flag, int backend)
{
    if (backend == 1)
        cpuComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
          rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, type, 
          ai, aj, twojmax, idxu_max, ijnum, switch_flag); 
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
          rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, type, 
          ai, aj, twojmax, idxu_max, ijnum, switch_flag); 
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
          rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, type, 
          ai, aj, twojmax, idxu_max, ijnum, switch_flag); 
#endif                                        
#ifdef USE_CUDA                
    if (backend == 2)
        gpuComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
          rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, type, 
          ai, aj, twojmax, idxu_max, ijnum, switch_flag); 
#endif                          
}

inline void ComputeDeidrj(dstype *dedr, dstype *ylist_r, dstype *ylist_i, dstype *dulist_r, dstype *dulist_i,         
        int *idxu_block, int *type, int *map, int *ai, int *aj, int nelements, int twojmax, int idxu_max, 
        int chemflag, int ijnum, int backend)
{
    if (backend == 1)
        cpuComputeDeidrj(dedr, ylist_r, ylist_i, dulist_r, dulist_i, idxu_block, type, map, 
                ai, aj, nelements, twojmax, idxu_max, chemflag, ijnum); 
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeDeidrj(dedr, ylist_r, ylist_i, dulist_r, dulist_i, idxu_block, type, map, 
                ai, aj, nelements, twojmax, idxu_max, chemflag, ijnum); 
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeDeidrj(dedr, ylist_r, ylist_i, dulist_r, dulist_i, idxu_block, type, map, 
                ai, aj, nelements, twojmax, idxu_max, chemflag, ijnum); 
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeDeidrj(dedr, ylist_r, ylist_i, dulist_r, dulist_i, idxu_block, type, map, 
                ai, aj, nelements, twojmax, idxu_max, chemflag, ijnum); 
#endif                          
}

inline void ComputeSna(dstype *sna, dstype *blist, int *ilist, int *mask, 
        int ncoeff, int nperdim, int inum, int quadraticflag, int backend)
{
    if (backend == 1)
        cpuComputeSna(sna, blist, ilist, mask, ncoeff, nperdim, inum, quadraticflag);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeSna(sna, blist, ilist, mask, ncoeff, nperdim, inum, quadraticflag);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeSna(sna, blist, ilist, mask, ncoeff, nperdim, inum, quadraticflag);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeSna(sna, blist, ilist, mask, ncoeff, nperdim, inum, quadraticflag);
#endif                          
}

inline void ComputeSna(dstype *sna, dstype *blist, int *ilist, int *mask, int *type,
        int ncoeff, int ntype, int nperdim, int inum, int quadraticflag, int backend)
{
    if (backend == 1)
        cpuComputeSna(sna, blist, ilist, mask, type, ncoeff, ntype, nperdim, inum, quadraticflag);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeSna(sna, blist, ilist, mask, type, ncoeff, ntype, nperdim, inum, quadraticflag);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeSna(sna, blist, ilist, mask, type, ncoeff, ntype, nperdim, inum, quadraticflag);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeSna(sna, blist, ilist, mask, type, ncoeff, ntype, nperdim, inum, quadraticflag);
#endif                          
}

inline void ComputeSnad(dstype *snad, dstype *dblist, dstype *blist, int *aii, int *ai, int *aj, int *ti, 
        int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag, int backend)
{
    if (backend == 1)
        cpuComputeSnad(snad, dblist, blist, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeSnad(snad, dblist, blist, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeSnad(snad, dblist, blist, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeSnad(snad, dblist, blist, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif                          
}

inline void ComputeSnad(dstype *snad, dstype *dblist, dstype *blist, int *aii, int *ai, int *aj, int *ti, 
        int *mask, int *tag, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag, int backend)
{
    if (backend == 1)
        cpuComputeSnad(snad, dblist, blist, aii, ai, aj, 
          ti, mask, tag, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeSnad(snad, dblist, blist, aii, ai, aj, 
          ti, mask, tag, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeSnad(snad, dblist, blist, aii, ai, aj, 
          ti, mask, tag, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeSnad(snad, dblist, blist, aii, ai, aj, 
          ti, mask, tag, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif                          
}

inline void ComputeSnav(dstype *snav, dstype *dblist, dstype *blist, dstype *x, int *aii, int *ai, int *aj,
        int *ti, int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag, int backend)
{
    if (backend == 1)
        cpuComputeSnav(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeSnav(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeSnav(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeSnav(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif                          
}

inline void ComputeSnav2(dstype *snav, dstype *dblist, dstype *blist, dstype *x, int *aii, int *ai, int *aj,
        int *ti, int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag, int backend)
{
    if (backend == 1)
        cpuComputeSnav2(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#ifdef USE_OMP  
    if (backend == 4)
        ompComputeSnav2(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipComputeSnav2(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif            
#ifdef USE_CUDA   
    if (backend == 2)
        gpuComputeSnav2(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ijnum, quadraticflag);
#endif                          
}

inline void SnapTallyEnergyFull(dstype *eatom, dstype *bispectrum, dstype *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag, int backend)
{
    if (backend == 1)
        cpuSnapTallyEnergyFull(eatom, bispectrum, coeffelem, ilist, map, type, 
                inum, ncoeff, ncoeffall, quadraticflag);
#ifdef USE_OMP  
    if (backend == 4)
        ompSnapTallyEnergyFull(eatom, bispectrum, coeffelem, ilist, map, type, 
                inum, ncoeff, ncoeffall, quadraticflag);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipSnapTallyEnergyFull(eatom, bispectrum, coeffelem, ilist, map, type, 
                inum, ncoeff, ncoeffall, quadraticflag);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuSnapTallyEnergyFull(eatom, bispectrum, coeffelem, ilist, map, type, 
                inum, ncoeff, ncoeffall, quadraticflag);
#endif                          
}

inline void SnapTallyForceFull(dstype *fatom, dstype *fij, int *ai, int *aj, int ijnum, int backend)
{
    if (backend == 1)
        cpuSnapTallyForceFull(fatom, fij, ai, aj, ijnum);
#ifdef USE_OMP  
    if (backend == 4)
        ompSnapTallyForceFull(fatom, fij, ai, aj, ijnum);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipSnapTallyForceFull(fatom, fij, ai, aj, ijnum);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuSnapTallyForceFull(fatom, fij, ai, aj, ijnum);
#endif                          
}

inline void SnapTallyVirialFull(dstype *vatom, dstype *fij, dstype *rij, 
        int *ai, int *aj, int ijnum, int backend)
{
    if (backend == 1)
        cpuSnapTallyVirialFull(vatom, fij, rij, ai, aj, ijnum);
#ifdef USE_OMP  
    if (backend == 4)
        ompSnapTallyVirialFull(vatom, fij, rij, ai, aj, ijnum);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipSnapTallyVirialFull(vatom, fij, rij, ai, aj, ijnum);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuSnapTallyVirialFull(vatom, fij, rij, ai, aj, ijnum);
#endif                          
}

inline void NeighborPairList(int *pairnum, int *pairlist, dstype *x, dstype rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim, int backend)
{
    if (backend == 1)
        cpuNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, 
                            neighnum, inum, jnum, dim);
#ifdef USE_OMP  
    if (backend == 4)
        ompuNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, 
                            neighnum, inum, jnum, dim);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, 
                            neighnum, inum, jnum, dim);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, 
                            neighnum, inum, jnum, dim);
#endif                          
}

inline void NeighborPairList(int *pairnum, int *pairlist, dstype *x, dstype *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int *atomtype, int inum, int jnum, int dim, int ntypes, int backend)
{
    if (backend == 1)
        cpuNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, 
                            neighnum, atomtype, inum, jnum, dim, ntypes);
#ifdef USE_OMP  
    if (backend == 4)
        ompNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, 
                            neighnum, atomtype, inum, jnum, dim, ntypes);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, 
                            neighnum, atomtype, inum, jnum, dim, ntypes);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuNeighPairList(pairnum, pairlist, x, rcutsq, ilist, neighlist, 
                            neighnum, atomtype, inum, jnum, dim, ntypes);
#endif                          
}
        
inline void NeighborPairs(dstype *xij, dstype *x, int *aii, int *ai, int *aj,  
      int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, 
      int *atomtype, int inum, int jnum, int dim, int backend)
{
    if (backend == 1)
        cpuNeighPairs(xij, x, aii, ai, aj,  ti, tj, pairnum, pairlist, pairnumsum, ilist, 
                        atomtype, inum, jnum, dim);
#ifdef USE_OMP  
    if (backend == 4)
        ompNeighPairs(xij, x, aii, ai, aj,  ti, tj, pairnum, pairlist, pairnumsum, ilist, 
                        atomtype, inum, jnum, dim);
#endif                                
#ifdef USE_HIP  
    if (backend == 3)
        hipNeighPairs(xij, x, aii, ai, aj,  ti, tj, pairnum, pairlist, pairnumsum, ilist, 
                        atomtype, inum, jnum, dim);
#endif                                        
#ifdef USE_CUDA   
    if (backend == 2)
        gpuNeighPairs(xij, x, aii, ai, aj,  ti, tj, pairnum, pairlist, pairnumsum, ilist, 
                        atomtype, inum, jnum, dim);
#endif                          
}


#endif  

