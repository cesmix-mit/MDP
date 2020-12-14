#ifndef __CORE_H__
#define __CORE_H__

#include "commoncore.h"

#ifdef HAVE_ONETHREAD
#include "opucore.h"
#endif                

#ifdef HAVE_OPENMP        
#include "cpucore.h"
#endif                

#ifdef HAVE_CUDA      
#include "gpucore.h"
#endif                

// static void ArraySetValueAtIndex(dstype *y, dstype a, Int i, Int backend)
// {
//     if (backend < 2)  
//         cpuArraySetValueAtIndex(y, a, i);
// #ifdef HAVE_CUDA            
//     if (backend == 2)  // CUDA C                
//         gpuArraySetValueAtIndex(y, a, i);
// #endif                  
// }
// 
static void ApplyGivensRotation(dstype *H, dstype *s, dstype *cs, dstype *sn, Int i, Int backend)
{
    if (backend < 2)     
        cpuApplyGivensRotation(H, s, cs, sn, i);
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuApplyGivensRotation(H, s, cs, sn, i);
#endif                  
}

static void BackSolve(dstype *y, dstype *H, dstype *s, Int i, Int n, Int backend)
{
    if (backend < 2)  
        cpuBackSolve(y, H, s, i, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuBackSolve(y, H, s, i, n);
#endif                  
}

static void GetArrayAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD           
    if (backend == 0)        // One thread CPU
        opuGetArrayAtIndex(y, x, ind, n);
#endif            
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuGetArrayAtIndex(y, x, ind, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuGetArrayAtIndex(y, x, ind, n);
#endif                  
}

static void PutArrayAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD           
    if (backend == 0)        // One thread CPU
        opuPutArrayAtIndex(y, x, ind, n);
#endif            
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuPutArrayAtIndex(y, x, ind, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuPutArrayAtIndex(y, x, ind, n);
#endif                  
}

static void ArrayPlusXAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD           
    if (backend == 0)        // One thread CPU
        opuArrayPlusXAtIndex(y, x, ind, n);
#endif            
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayPlusXAtIndex(y, x, ind, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayPlusXAtIndex(y, x, ind, n);
#endif                  
}

static void ArrayMinusXAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD           
    if (backend == 0)        // One thread CPU
        opuArrayMinusXAtIndex(y, x, ind, n);
#endif            
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayMinusXAtIndex(y, x, ind, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayMinusXAtIndex(y, x, ind, n);
#endif                  
}

// static void ArrayAXPYAtIndex(dstype *y, dstype *x, dstype a, Int *ind, Int n, Int backend)
// {
// #ifdef HAVE_ONETHREAD           
//     if (backend == 0)        // One thread CPU
//         opuArrayAXPYAtIndex(y, x, a, ind, n);
// #endif            
// #ifdef HAVE_OPENMP        
//     if (backend == 1)  // Open MP
//         cpuArrayAXPYAtIndex(y, x, a, ind, n);
// #endif                 
// #ifdef HAVE_CUDA            
//     if (backend == 2)  // CUDA C                
//         gpuArrayAXPYAtIndex(y, x, a, ind, n);
// #endif                  
// }

void PutFaceNodes(dstype *udg, dstype *uh, Int *ind1, Int *ind2, Int n, Int opts, Int backend)
{        
    if (opts==0) {
        ArrayMinusXAtIndex(udg, uh, ind1, n, backend);
        ArrayPlusXAtIndex(udg, uh, ind2, n, backend);
    }
    else {
        ArrayMinusXAtIndex(udg, uh, ind1, n, backend);
    }
}

// static void ArrayAverageAtIndex(dstype *y, dstype *x, Int *ind1, Int *ind2, Int n, Int backend)
// {
// #ifdef HAVE_ONETHREAD           
//     if (backend == 0)        // One thread CPU
//         opuArrayAverageAtIndex(y, x, ind1, ind2, n);
// #endif            
// #ifdef HAVE_OPENMP        
// //     if (backend == 1)  // Open MP
// //         cpuArrayAverageAtIndex(y, x, ind1, ind2, n);
// #endif                 
// #ifdef HAVE_CUDA            
// //     if (backend == 2)  // CUDA C                
// //         gpuArrayAverageAtIndex(y, x, ind1, ind2, n);
// #endif                  
// }
// 
static void ArraySquare(dstype *y, dstype *x, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD           
    if (backend == 0)        // One thread CPU
        opuArraySquare(y, x, n);
#endif            
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArraySquare(y, x, n);   
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArraySquare(y, x, n);
#endif                  
}

static void ArraySetValue(dstype *y, dstype a, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD       
    if (backend == 0)        // One thread CPU
        opuArraySetValue(y, a, n);
#endif         
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArraySetValue(y, a, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArraySetValue(y, a, n);
#endif                  
}

static void ArrayAddScalar(dstype *y, dstype a, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD         
    if (backend == 0)        // One thread CPU
        opuArrayAddScalar(y, a, n);
#endif             
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayAddScalar(y, a, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayAddScalar(y, a, n);
#endif                  
}

static void ArrayMultiplyScalar(dstype *y, dstype a, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD         
    if (backend == 0)       // One thread CPU
        opuArrayMultiplyScalar(y, a, n);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayMultiplyScalar(y, a, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayMultiplyScalar(y, a, n);
#endif                  
}

static void ArrayCopy(dstype *y, dstype *x, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD           
    if (backend == 0)        // One thread CPU
        opuArrayCopy(y, x, n);
#endif            
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayCopy(y, x, n);   
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayCopy(y, x, n);
#endif                  
}

static void ArrayMinus(dstype *y, dstype *x, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD           
    if (backend == 0)        // One thread CPU
        opuArrayMinus(y, x, n);
#endif            
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayMinus(y, x, n);   
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayMinus(y, x, n);
#endif                  
}

static void ArrayAbs(dstype *y, dstype *x, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD           
    if (backend == 0)        // One thread CPU
        opuArrayAbs(y, x, n);
#endif            
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayAbs(y, x, n);   
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayAbs(y, x, n);
#endif                  
}

static void ArraySqrt(dstype *y, dstype *x, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD           
    if (backend == 0)        // One thread CPU
        opuArraySqrt(y, x, n);
#endif            
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArraySqrt(y, x, n);   
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArraySqrt(y, x, n);
#endif                  
}

static void ArrayMultiplyScalarDiagonal(dstype *C, dstype a, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD           
    if (backend == 0)        // One thread CPU
        opuArrayMultiplyScalarDiagonal(C, a, n);
#endif            
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayMultiplyScalarDiagonal(C, a, n);   
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayMultiplyScalarDiagonal(C, a, n);
#endif                  
}

static void ArrayAddVectorToDiagonal(dstype *C, dstype *x, dstype a, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD           
    if (backend == 0)        // One thread CPU
        opuArrayAddVectorToDiagonal(C, x, a, n);
#endif            
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayAddVectorToDiagonal(C, x, a, n);   
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayAddVectorToDiagonal(C, x, a, n);
#endif                  
}

// static void ArrayRowAverage(dstype *y, dstype *x, Int m, Int n, Int backend)
// {
// #ifdef HAVE_ONETHREAD           
//     if (backend == 0)        // One thread CPU
//         opuArrayRowAverage(y, x, m, n);
// #endif            
// #ifdef HAVE_OPENMP        
//     if (backend == 1)  // Open MP
//         cpuArrayRowAverage(y, x, m, n);   
// #endif                 
// #ifdef HAVE_CUDA            
//     if (backend == 2)  // CUDA C                
//         gpuArrayRowAverage(y, x, m, n);
// #endif                  
// }

static void ArrayAXPB(dstype *y, dstype *x, dstype a, dstype b, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)      // One thread CPU
        opuArrayAXPB(y, x, a, b, n);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayAXPB(y, x, a, b, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayAXPB(y, x, a, b, n);
#endif                  
}

static void ArrayAXPBY(dstype *z, dstype *x, dstype *y, dstype a, dstype b, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuArrayAXPBY(z, x, y, a, b, n);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayAXPBY(z, x, y, a, b, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayAXPBY(z, x, y, a, b, n);
#endif                  
}

static void ArrayAXY(dstype *z, dstype *x, dstype *y, dstype a, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuArrayAXY(z, x, y, a, n);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayAXY(z, x, y, a, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayAXY(z, x, y, a, n);
#endif                  
}

static void ArrayAXYZ(dstype *s, dstype *x, dstype *y, dstype *z, dstype a, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuArrayAXYZ(s, x, y, z, a, n);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayAXYZ(s, x, y, z, a, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayAXYZ(s, x, y, z, a, n);
#endif                  
}

static void ArrayAXYPBZ(dstype *s, dstype *x, dstype *y, dstype *z, dstype a, dstype b, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuArrayAXYPBZ(s, x, y, z, a, b, n);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayAXYPBZ(s, x, y, z, a, b, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayAXYPBZ(s, x, y, z, a, b, n);
#endif                  
}

static void ArrayAdd3Vectors(dstype *s, dstype *x, dstype *y, dstype *z, dstype a, dstype b, dstype c, Int n, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuArrayAdd3Vectors(s, x, y, z, a, b, c, n);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayAdd3Vectors(s, x, y, z, a, b, c, n);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayAdd3Vectors(s, x, y, z, a, b, c, n);
#endif                  
}

static void ArrayExtract(dstype *un, dstype *u, Int I, Int J, Int K, 
        Int i1, Int i2, Int j1, Int j2, Int k1, Int k2, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuArrayExtract(un, u, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayExtract(un, u, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayExtract(un, u, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                  
}

static void ArrayInsert(dstype *u, dstype *un, Int I, Int J, Int K, 
        Int i1, Int i2, Int j1, Int j2, Int k1, Int k2, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)     // One thread CPU
        opuArrayInsert(u, un, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayInsert(u, un, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayInsert(u, un, I, J, K, i1, i2, j1, j2, k1, k2);    
#endif                  
}

static void ArrayGemmBatch(dstype *C, dstype *A, dstype *B, Int I, Int J, Int K, Int S, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuArrayGemmBatch(C, A, B, I, J, K, S);    
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayGemmBatch(C, A, B, I, J, K, S);    
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayGemmBatch(C, A, B, I, J, K, S);    
#endif                  
}

static void ArrayGemmBatch1(dstype *C, dstype *A, dstype *B, Int I, Int J, Int K, Int S, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuArrayGemmBatch1(C, A, B, I, J, K, S);    
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayGemmBatch1(C, A, B, I, J, K, S);    
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayGemmBatch1(C, A, B, I, J, K, S);    
#endif                  
}

static void ArrayDG2CG(dstype *ucg, dstype *udg, Int *cgent2dgent, Int *rowent2elem, Int nent, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuArrayDG2CG(ucg, udg, cgent2dgent, rowent2elem, nent);    
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayDG2CG(ucg, udg, cgent2dgent, rowent2elem, nent);    
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayDG2CG(ucg, udg, cgent2dgent, rowent2elem, nent);    
#endif                  
}

static void ArrayDG2CG2(dstype *ucg, dstype *udg, Int *colent2elem, Int *rowent2elem, Int nent, Int npe, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuArrayDG2CG2(ucg, udg, colent2elem, rowent2elem, nent, npe);    
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuArrayDG2CG2(ucg, udg, colent2elem, rowent2elem, nent, npe);    
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        gpuArrayDG2CG2(ucg, udg, colent2elem, rowent2elem, nent, npe);    
#endif                  
}

static void GetElemNodes(dstype *un, dstype *u, Int np, Int nc, Int nc1, Int nc2, Int e1, Int e2, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuGetElemNodes(un, u, np, nc, nc1, nc2, e1, e2);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuGetElemNodes(un, u, np, nc, nc1, nc2, e1, e2);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                                    
        gpuGetElemNodes(un, u, np*(e2-e1)*(nc2-nc1), np*(e2-e1), np, e1, nc1, nc);    
#endif                      
}

static void PutElemNodes(dstype *u, dstype *un, Int np, Int nc, Int nc1, Int nc2, Int e1, Int e2, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuPutElemNodes(u, un, np, nc, nc1, nc2, e1, e2);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuPutElemNodes(u, un, np, nc, nc1, nc2, e1, e2);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                                    
        gpuPutElemNodes(u, un, np*(e2-e1)*(nc2-nc1), np*(e2-e1), np, e1, nc1, nc);    
#endif                      
}

static void GetFaceNodes(dstype *uh, dstype *udg, Int *facecon, Int npf, Int ncu, Int npe, Int nc, Int f1, Int f2, Int opts, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuGetFaceNodes(uh, udg, facecon, npf, ncu, npe, nc, f1, f2, opts);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuGetFaceNodes(uh, udg, facecon, npf, ncu, npe, nc, f1, f2, opts);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2) { // CUDA C                                    
        if (opts==0)
            gpuGetFaceNodes0(uh, udg, facecon, npf*(f2-f1)*ncu, npf*(f2-f1), npf, npe, nc, f1);
        else if (opts==1)
            gpuGetFaceNodes1(uh, udg, facecon, npf*(f2-f1)*ncu, npf*(f2-f1), npf, npe, nc, f1);
        else if (opts==2)
            gpuGetFaceNodes2(uh, udg, facecon, npf*(f2-f1)*ncu, npf*(f2-f1), npf, npe, nc, f1);
    }
#endif                      
}

static void PutFaceNodes(dstype *udg, dstype *uh, Int *facecon, Int npf, Int ncu, Int npe, Int nc, Int f1, Int f2, Int opts, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuPutFaceNodes(udg, uh, facecon, npf, ncu, npe, nc, f1, f2, opts);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuPutFaceNodes(udg, uh, facecon, npf, ncu, npe, nc, f1, f2, opts);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2) { // CUDA C                                    
        if (opts==0)
            gpuPutFaceNodes0(udg, uh, facecon, npf*(f2-f1)*ncu, npf*(f2-f1), npf, npe, nc, f1);
        else 
            gpuPutFaceNodes1(udg, uh, facecon, npf*(f2-f1)*ncu, npf*(f2-f1), npf, npe, nc, f1);
    }
#endif                      
}

static void PutFaceNodes(dstype *udg, dstype *uh, int *rowe2f1, int *cole2f1, int *ent2ind1,
        int *rowe2f2, int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int opts, Int backend)
{
#ifdef HAVE_ONETHREAD      
    if (backend == 0)        // One thread CPU
        opuPutFaceNodes(udg, uh, rowe2f1, cole2f1, ent2ind1, rowe2f2, cole2f2, ent2ind2, npf, npe, nc, e1, e2, opts);
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        cpuPutFaceNodes(udg, uh, rowe2f1, cole2f1, ent2ind1, rowe2f2, cole2f2, ent2ind2, npf, npe, nc, e1, e2, opts);
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2) { // CUDA C                                    
        gpuPutFaceNodes(udg, uh, rowe2f1, cole2f1, ent2ind1, rowe2f2, cole2f2, ent2ind2, npf, npe, nc, e1, e2, opts);
    }
#endif                      
}

#endif  

